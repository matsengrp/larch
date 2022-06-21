#include "usher_optimize.hpp"

namespace MAT = Mutation_Annotated_Tree;

extern std::atomic_bool interrupted;
extern bool use_bound;
extern int process_count;
extern int this_rank;
extern uint32_t num_threads;

#define DRIFT_MASK 0x80000000
#define ALL_DIR_MASK 0x40000000

int count_back_mutation(const MAT::Tree &tree);
void get_pos_samples_old_tree(MAT::Tree &tree, std::vector<mutated_t> &output);
void min_back_reassign_state_local(MAT::Tree &tree,
                                   const std::vector<mutated_t> &mutations);
void MPI_min_back_reassign_states(MAT::Tree &tree,
                                  const std::vector<mutated_t> &mutations,
                                  int start_position);

static void make_output_path(std::string &path_template) {
  auto fd = mkstemps(const_cast<char *>(path_template.c_str()), 3);
  close(fd);
}

void InitUsherMPI(int argc, char **argv) {
  int ignored;
  auto init_result = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &ignored);
  if (init_result != MPI_SUCCESS) {
    fprintf(stderr, "MPI init failed\n");
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &process_count);
  fprintf(stderr, "Running with %d processes\n", process_count);
}

void UsherOptimize(Mutation_Annotated_Tree::Tree &t) {
  t.uncondense_leaves();
  t.populate_ignored_range();

  int drift_iterations = 0;
  int radius = -1;
  float search_proportion = 2;
  int rand_sel_seed = 0;
  unsigned int minutes_between_save = 0;
  bool no_write_intermediate = true;

  auto search_end_time = std::chrono::steady_clock::time_point::max();
  auto start_time = std::chrono::steady_clock::now();
  auto save_period = std::chrono::minutes(minutes_between_save);
  FILE *movalbe_src_log = fopen("/dev/null", "w");
  bool log_moves = false;

  auto last_save_time = std::chrono::steady_clock::now();
  size_t new_score;
  size_t score_before;
  int stalled = -1;
  score_before = t.get_parsimony_score();
  new_score = score_before;
  fprintf(stderr, "after state reassignment:%zu\n", score_before);
  fprintf(stderr, "Height:%zu\n", t.get_max_level());

  std::vector<MAT::Node *> nodes_to_search;
  std::vector<MAT::Node *> bfs_ordered_nodes;
  bfs_ordered_nodes = t.breadth_first_expansion();

  bool isfirst = true;
  bool allow_drift = false;
  int iteration = 1;
  int max_round = 1000;
  float min_improvement = 0.0005;
  bool reduce_back_mutations = true;
  num_threads = tbb::task_scheduler_init::default_num_threads();
  tbb::task_scheduler_init init(num_threads);
  std::string intermediate_writing;
  std::string intermediate_pb_base_name = "";
  std::string intermediate_template =
      intermediate_pb_base_name + "temp_eriting_XXXXXX.pb";
  std::string intermediate_nwk_out = "";
  while (stalled < drift_iterations) {
    bfs_ordered_nodes = t.breadth_first_expansion();
    fputs("Start Finding nodes to move \n", stderr);
    bool search_all_nodes = false;
    bool search_all_dir = false;
    if (isfirst || allow_drift) {
      search_all_nodes = true;
      search_all_dir = true;
    } else if (radius < 0 && radius >= -2 * (int)t.max_level) {
      radius *= 2;
      search_all_nodes = true;
      search_all_dir = true;
    }
    find_nodes_to_move(bfs_ordered_nodes, nodes_to_search, search_all_nodes,
                       search_all_dir, radius, t);
    if (search_proportion < 1) {
      std::vector<MAT::Node *> nodes_to_search_temp;
      nodes_to_search_temp.reserve(nodes_to_search.size() * search_proportion);
      std::mt19937_64 rng(rand_sel_seed);
      std::sample(nodes_to_search.begin(), nodes_to_search.end(),
                  std::back_inserter(nodes_to_search_temp),
                  size_t(std::round(nodes_to_search.size() * search_proportion)), rng);
      nodes_to_search.swap(nodes_to_search_temp);
    }
    isfirst = false;
    fprintf(stderr, "%zu nodes to search\n", nodes_to_search.size());
    if (nodes_to_search.empty()) {
      break;
    }
    bool isfirst_this_iter = true;
    // Actual optimization loop
    while (!nodes_to_search.empty()) {
      auto dfs_ordered_nodes = t.depth_first_expansion();
      std::mt19937_64 rng;
      std::shuffle(nodes_to_search.begin(), nodes_to_search.end(), rng);
      bool distribute = (process_count > 1) && (nodes_to_search.size() > 1000);
      if (distribute) {
        MPI_Request req;
        int radius_to_boardcast = abs(radius);
        if (allow_drift) {
          radius_to_boardcast |= DRIFT_MASK;
        }
        if (search_all_dir) {
          radius_to_boardcast |= ALL_DIR_MASK;
          fprintf(stderr, "Search all directions\n");
        }
        MPI_Ibcast(&radius_to_boardcast, 1, MPI_INT, 0, MPI_COMM_WORLD, &req);
        fprintf(stderr, "Sent radius\n");
        MPI_Wait(&req, MPI_STATUS_IGNORE);
        fprintf(stderr, "Start Send tree\n");
        t.MPI_send_tree();
      }
      adjust_all(t);
      use_bound = true;
      std::vector<size_t> nodes_to_search_idx;
      nodes_to_search_idx.reserve(nodes_to_search.size());
      for (const auto node : nodes_to_search) {
        nodes_to_search_idx.push_back(node->dfs_index);
      }
      std::vector<size_t> defered_nodes;
      auto next_save_time = minutes_between_save
                                ? last_save_time + save_period
                                : std::chrono::steady_clock::time_point::max();
      bool do_continue = true;
      auto search_stop_time = next_save_time;
      if (no_write_intermediate || search_end_time < next_save_time) {
        search_stop_time = search_end_time;
      }
      optimize_tree_main_thread(
          nodes_to_search_idx, t, std::abs(radius), movalbe_src_log, allow_drift,
          log_moves ? iteration : -1, defered_nodes, distribute, search_stop_time,
          do_continue, search_all_dir, isfirst_this_iter
#ifdef CHECK_STATE_REASSIGN
          ,
          origin_states
#endif
      );
      isfirst_this_iter = false;
      fprintf(stderr, "Defered %zu nodes\n", defered_nodes.size());
      nodes_to_search.reserve(defered_nodes.size());
      nodes_to_search.clear();
      for (auto idx : defered_nodes) {
        if (t.get_node(idx)) {
          nodes_to_search.push_back(t.get_node(idx));
        }
      }
      auto curr_score = t.get_parsimony_score();
      if (curr_score >= new_score) {
        nodes_to_search.clear();
      }
      new_score = curr_score;
      fprintf(stderr,
              "parsimony score after optimizing: %zu,with radius %d, second from start "
              "%ld \n\n",
              new_score, std::abs(radius),
              std::chrono::duration_cast<std::chrono::seconds>(
                  std::chrono::steady_clock::now() - start_time)
                  .count());
      if (!no_write_intermediate) {
        intermediate_writing = intermediate_template;
        make_output_path(intermediate_writing);
        auto save_start = std::chrono::steady_clock::now();
        t.save_detailed_mutations(intermediate_writing);
        rename(intermediate_writing.c_str(), intermediate_pb_base_name.c_str());
        last_save_time = std::chrono::steady_clock::now();
        fprintf(stderr, "Took %ldsecond to save intermediate protobuf\n",
                std::chrono::duration_cast<std::chrono::seconds>(last_save_time -
                                                                 save_start)
                    .count());
      }
      if (allow_drift) {
        MAT::save_mutation_annotated_tree(
            t, intermediate_nwk_out + std::to_string(iteration) + ".pb.gz");
      }
      if (std::chrono::steady_clock::now() >= search_end_time) {
        break;
      }
      if (interrupted) {
        break;
      }
      if (allow_drift) {
        nodes_to_search.clear();
      }
      search_all_dir = true;
    }
    if (interrupted) {
      break;
    }
    float improvement = 1 - ((float)new_score / (float)score_before);
    fprintf(stderr, "Last round improvement %f\n", improvement);
    if (improvement < min_improvement) {
      fprintf(stderr, "Less than minimium improvement,stalled for %d iterations\n",
              stalled);
      fprintf(stderr, "Will drift for %d iterations \n", drift_iterations);
      stalled++;
      allow_drift = true;
    } else {
      score_before = new_score;
      stalled = -1;
    }
    iteration++;
    if (iteration >= max_round) {
      fprintf(stderr, "Reached %d interations\n", iteration);
      break;
    }
    if (std::chrono::steady_clock::now() >= search_end_time) {
      fprintf(stderr, "Exceeded search time \n");
      break;
    }
  }
  int temp = 0;
  MPI_Request req;
  MPI_Ibcast(&temp, 1, MPI_INT, 0, MPI_COMM_WORLD, &req);
  if (reduce_back_mutations) {
    fprintf(stderr, "Parsimony score before %zu\n", t.get_parsimony_score());
    fprintf(stderr, "Back mutation count before %d\n", count_back_mutation(t));
    std::vector<mutated_t> output(MAT::Mutation::refs.size());
    get_pos_samples_old_tree(t, output);
    if (process_count == 1) {
      min_back_reassign_state_local(t, output);
    } else {
      MPI_min_back_reassign_states(t, output, 0);
    }
    fprintf(stderr, "Back mutation count after %d\n", count_back_mutation(t));
  }
  fprintf(stderr, "Final Parsimony score %zu\n", t.get_parsimony_score());
  fclose(movalbe_src_log);
  MPI_Wait(&req, MPI_STATUS_IGNORE);
}
