#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <chrono>

#include <unistd.h>
#include <sys/wait.h>

#include "arguments.hpp"
#include "larch/dag_loader.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/weight_accumulator.hpp"
#include "larch/subtree/tree_count.hpp"
#include "larch/subtree/parsimony_score_binary.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/spr/spr_view.hpp"
#include "larch/merge/merge.hpp"
#include "larch/merge/leaf_set.hpp"
#include "larch/spr/batching_callback.hpp"
#include "benchmark.hpp"
#include <mpi.h>

#include "larch/usher_glue.hpp"

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "larch-usher -i,--input file -o,--output file [-m,--matopt file] "
               "[-c,--count number]\n";
  std::cout << "  -i,--input   Path to input DAG\n";
  std::cout << "  -r,--MAT-refseq-file   Provide a path to a file containing a "
               "reference sequence\nif input points to MAT protobuf\n";
  std::cout << "  -o,--output  Path to output DAG\n";
  std::cout << "  -m,--matopt  Path to matOptimize executable. Default: matOptimize\n";
  std::cout << "  -l,--logpath Path for logging\n";
  std::cout << "  -c,--count   Number of iterations. Default: 1\n";
  std::cout << "  -s,--switch-subtree           Switch to optimizing subtrees after "
               "the specified "
               "number of iterations (default never)\n";
  std::cout << "  --min-subtree-clade-size      The minimum number of leaves in a "
               "subtree sampled for optimization (default 100, ignored without option "
               "`-s`)\n";
  std::cout << "  --max-subtree-clade-size      The maximum number of leaves in a "
               "subtree sampled for optimization (default 1000, ignored without option "
               "`-s`)\n";
  std::cout << "  --uniform-subtree-root        Choose subtree root node uniformly"
               "from allowed options. Default choice is weighted by (1 + m^2) where m "
               "is minimum "
               "mutations on a node's parent edge\n";
  std::cout
      << "  --move-coeff-nodes   New node coefficient for scoring moves. Default: 1\n";
  std::cout << "  --move-coeff-pscore  Parsimony score coefficient for scoring moves. "
               "Default: 1\n";
  std::cout << "  --sample-any-tree    Sample any tree for optimization, rather than "
               "requiring "
               "the sampled tree to maximize parsimony.\n";
  std::cout << "  --callback-option    Callback configuration choice(default merge all "
               "profitable moves)\n";

  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
}

static int ParseNumber(std::string_view str) {
  int result{};
  std::istringstream stream{std::string{str}};
  stream >> result;
  if (stream.fail()) {
    throw std::runtime_error("Invalid number");
  }
  return result;
}

void check_edge_mutations(MADAG madag);

std::vector<std::vector<const CompactGenome*>> clades_union(
    const std::vector<std::vector<const CompactGenome*>>& lhs,
    const std::vector<std::vector<const CompactGenome*>>& rhs) {
  std::vector<std::vector<const CompactGenome*>> result;

  for (auto [lhs_clade, rhs_clade] : ranges::views::zip(lhs, rhs)) {
    std::vector<const CompactGenome*> clade{lhs_clade};
    clade.insert(clade.end(), rhs_clade.begin(), rhs_clade.end());
    ranges::sort(clade);
    ranges::unique(clade);
    result.push_back(std::move(clade));
  }

  ranges::sort(result);
  return result;
}

std::vector<std::vector<const CompactGenome*>> clades_difference(
    const std::vector<std::vector<const CompactGenome*>>& lhs,
    const std::vector<std::vector<const CompactGenome*>>& rhs) {
  std::vector<std::vector<const CompactGenome*>> result;

  for (auto [lhs_clade, rhs_clade] : ranges::views::zip(lhs, rhs)) {
    std::vector<const CompactGenome*> clade;
    std::set_difference(lhs_clade.begin(), lhs_clade.end(), rhs_clade.begin(),
                        rhs_clade.end(), std::inserter(clade, clade.begin()));
    ranges::sort(clade);
    ranges::unique(clade);
    result.push_back(std::move(clade));
  }

  ranges::sort(result);
  return result;
}

template <typename SampleDAG>
struct Treebased_Move_Found_Callback
    : public BatchingCallback<Treebased_Move_Found_Callback<SampleDAG>, SampleDAG> {
  MOVE_ONLY_VIRT_DTOR(Treebased_Move_Found_Callback);

  Treebased_Move_Found_Callback(Merge& merge, SampleDAG sample_dag,
                                std::pair<int, int> move_score_coeffs)
      : BatchingCallback<Treebased_Move_Found_Callback<SampleDAG>,
                         SampleDAG>{merge, sample_dag},
        move_score_coeffs_{std::move(move_score_coeffs)} {};

  Treebased_Move_Found_Callback(Merge& merge, SampleDAG sample_dag)
      : BatchingCallback<Treebased_Move_Found_Callback<SampleDAG>,
                         SampleDAG>{merge, sample_dag},
        move_score_coeffs_{1, 1} {};

  template <typename SPRView, typename FragmentType>
  bool OnMove(SPRView spr, const FragmentType& fragment, Profitable_Moves& move,
              int /*best_score_change*/, std::vector<Node_With_Major_Allele_Set_Change>&
              /*nodes_with_major_allele_set_change*/) {
    int node_id_map_count = 0;
    if (move_score_coeffs_.first != 0) {
      auto src_leaf_set = this->GetMerge()
                              .GetResultNodeLabels()
                              .at(this->GetMappedStorage().Get(spr.GetMoveSource().GetId()).GetOriginalId())
                              .GetLeafSet()
                              ->GetClades();
      auto dst_leaf_set = this->GetMerge()
                              .GetResultNodeLabels()
                              .at(this->GetMappedStorage().Get(spr.GetMoveTarget().GetId()).GetOriginalId())
                              .GetLeafSet()
                              ->GetClades();
      for (auto hypothetical_node : fragment.GetNodes()) {
        if (hypothetical_node.IsMoveNew()) {
          if (not(this->GetMerge().ContainsLeafset(clades_union(src_leaf_set, dst_leaf_set)))) {
            ++node_id_map_count;
          }
        } else if (hypothetical_node.HasChangedTopology() and 
                   not (hypothetical_node.IsMoveSource() or hypothetical_node.IsMoveTarget())) {
          const auto& current_leaf_sets = this->GetMerge()
                                              .GetResultNodeLabels()
                                              .at(this->GetMappedStorage().Get(hypothetical_node.GetId()).GetOriginalId())
                                              .GetLeafSet()
                                              ->GetClades();
          if (not(this->GetMerge().ContainsLeafset(
                      clades_difference(current_leaf_sets, src_leaf_set)) and
                  this->GetMerge().ContainsLeafset(
                      clades_difference(current_leaf_sets, dst_leaf_set)))) {
            ++node_id_map_count;
          }
        }
      }
    }
    move.score_change = move_score_coeffs_.second * move.score_change -
                        move_score_coeffs_.first * node_id_map_count;
    return move.score_change <= 0;
  }

  void OnRadius(){};

  std::pair<int, int> move_score_coeffs_;
};

template <typename SampleDAG>
struct Merge_All_Moves_Found_Callback
    : public BatchingCallback<Merge_All_Moves_Found_Callback<SampleDAG>, SampleDAG> {
  MOVE_ONLY_VIRT_DTOR(Merge_All_Moves_Found_Callback);

  Merge_All_Moves_Found_Callback(Merge& merge, SampleDAG sample_dag)
      : BatchingCallback<Merge_All_Moves_Found_Callback<SampleDAG>, SampleDAG>{
            merge, sample_dag} {};

  template <typename SPRView, typename FragmentType>
  bool OnMove(SPRView /*spr*/, const FragmentType& /*fragment*/, Profitable_Moves& move,
              int best_score_change, std::vector<Node_With_Major_Allele_Set_Change>&
              /*nodes_with_major_allele_set_change*/) {
    return move.score_change < best_score_change;
  }

  void OnRadius(){};
};

template <typename SampleDAG>
struct Merge_All_Profitable_Moves_Found_Callback
    : public BatchingCallback<Merge_All_Profitable_Moves_Found_Callback<SampleDAG>, SampleDAG> {
  MOVE_ONLY_VIRT_DTOR(Merge_All_Profitable_Moves_Found_Callback);

  Merge_All_Profitable_Moves_Found_Callback(Merge& merge, SampleDAG sample_dag,
                                            std::pair<int, int> move_score_coeffs)
      : BatchingCallback<Merge_All_Profitable_Moves_Found_Callback<SampleDAG>,
                         SampleDAG>{merge, sample_dag},
        move_score_coeffs_{std::move(move_score_coeffs)} {};

  Merge_All_Profitable_Moves_Found_Callback(Merge& merge, SampleDAG sample_dag)
      : BatchingCallback<Merge_All_Profitable_Moves_Found_Callback<SampleDAG>,
                         SampleDAG>{merge, sample_dag},
        move_score_coeffs_{1, 1} {};

  template <typename SPRView, typename FragmentType>
  bool OnMove(SPRView spr, const FragmentType& fragment, Profitable_Moves& move,
              int /*best_score_change*/, std::vector<Node_With_Major_Allele_Set_Change>&
              /*nodes_with_major_allele_set_change*/) {
    int node_id_map_count = 0;
    if (move_score_coeffs_.first != 0) {
      auto src_leaf_set = this->GetMerge()
                              .GetResultNodeLabels()
                              .at(this->GetMappedStorage().Get(spr.GetMoveSource().GetId()).GetOriginalId())
                              .GetLeafSet()
                              ->GetClades();
      auto dst_leaf_set = this->GetMerge()
                              .GetResultNodeLabels()
                              .at(this->GetMappedStorage().Get(spr.GetMoveTarget().GetId()).GetOriginalId())
                              .GetLeafSet()
                              ->GetClades();
      for (auto hypothetical_node : fragment.GetNodes()) {
        if (hypothetical_node.IsMoveNew()) {
          if (not(this->GetMerge().ContainsLeafset(clades_union(src_leaf_set, dst_leaf_set)))) {
            ++node_id_map_count;
          }
        } else if (hypothetical_node.HasChangedTopology() and 
                   not (hypothetical_node.IsMoveSource() or hypothetical_node.IsMoveTarget())) {
          const auto& current_leaf_sets = this->GetMerge()
                                              .GetResultNodeLabels()
                                              .at(this->GetMappedStorage().Get(hypothetical_node.GetId()).GetOriginalId())
                                              .GetLeafSet()
                                              ->GetClades();
          if (not(this->GetMerge().ContainsLeafset(
                      clades_difference(current_leaf_sets, src_leaf_set)) and
                  this->GetMerge().ContainsLeafset(
                      clades_difference(current_leaf_sets, dst_leaf_set)))) {
            ++node_id_map_count;
          }
        }
      }
    }
    move.score_change = move_score_coeffs_.second * move.score_change -
                        move_score_coeffs_.first * node_id_map_count;
    return move.score_change <= 0;
  }

  void OnRadius(){};

  std::pair<int, int> move_score_coeffs_;
};

template <typename SampleDAG>
struct Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback
    : public BatchingCallback<
          Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback<SampleDAG>, SampleDAG> {
  MOVE_ONLY_VIRT_DTOR(Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback);

  Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback(
      Merge& merge, SampleDAG sample_dag, std::pair<int, int> move_score_coeffs)
      : BatchingCallback<
            Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback<SampleDAG>,
            SampleDAG>{merge, sample_dag},
        move_score_coeffs_{std::move(move_score_coeffs)} {};

  Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback(Merge& merge,
                                                       SampleDAG sample_dag)
      : BatchingCallback<
            Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback<SampleDAG>,
            SampleDAG>{merge, sample_dag},
        move_score_coeffs_{1, 1} {};

  template <typename SPRView, typename FragmentType>
  bool OnMove(SPRView spr, const FragmentType& fragment, Profitable_Moves& move,
              int /*best_score_change*/, std::vector<Node_With_Major_Allele_Set_Change>&
              /*nodes_with_major_allele_set_change*/) {
    int node_id_map_count = 0;
    if (move_score_coeffs_.first != 0) {
      auto src_leaf_set = this->GetMerge()
                              .GetResultNodeLabels()
                              .at(this->GetMappedStorage().Get(spr.GetMoveSource().GetId()).GetOriginalId())
                              .GetLeafSet()
                              ->GetClades();
      auto dst_leaf_set = this->GetMerge()
                              .GetResultNodeLabels()
                              .at(this->GetMappedStorage().Get(spr.GetMoveTarget().GetId()).GetOriginalId())
                              .GetLeafSet()
                              ->GetClades();
      for (auto hypothetical_node : fragment.GetNodes()) {
        if (hypothetical_node.IsMoveNew()) {
          if (not(this->GetMerge().ContainsLeafset(clades_union(src_leaf_set, dst_leaf_set)))) {
            ++node_id_map_count;
          }
        } else if (hypothetical_node.HasChangedTopology() and 
                   not (hypothetical_node.IsMoveSource() or hypothetical_node.IsMoveTarget())) {
          const auto& current_leaf_sets = this->GetMerge()
                                              .GetResultNodeLabels()
                                              .at(this->GetMappedStorage().Get(hypothetical_node.GetId()).GetOriginalId())
                                              .GetLeafSet()
                                              ->GetClades();
          if (not(this->GetMerge().ContainsLeafset(
                      clades_difference(current_leaf_sets, src_leaf_set)) and
                  this->GetMerge().ContainsLeafset(
                      clades_difference(current_leaf_sets, dst_leaf_set)))) {
            ++node_id_map_count;
          }
        }
      }
    }
    move.score_change = move_score_coeffs_.second * move.score_change -
                        move_score_coeffs_.first * node_id_map_count;
    return move.score_change <= 0;
  }

  void OnRadius(){};

  std::pair<int, int> move_score_coeffs_;
};

int main(int argc, char** argv) {  // NOLINT(bugprone-exception-escape)
  Arguments args = GetArguments(argc, argv);
  int ignored{};
  std::string input_dag_path;
  std::string output_dag_path;
  std::string matoptimize_path = "matOptimize";
  std::string logfile_path = "optimization_log";
  std::string refseq_path;
  std::string callback_config = "best-moves";
  bool sample_best_tree = true;
  size_t count = 1;
  int move_coeff_nodes = 1;
  int move_coeff_pscore = 1;
  size_t switch_subtrees = std::numeric_limits<size_t>::max();
  size_t min_subtree_clade_size = 100;   // NOLINT
  size_t max_subtree_clade_size = 1000;  // NOLINT
  bool uniform_subtree_root = false;

  for (auto [name, params] : args) {
    if (name == "-h" or name == "--help") {
      Usage();
    } else if (name == "-i" or name == "--input") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      input_dag_path = *params.begin();
    } else if (name == "-o" or name == "--output") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      output_dag_path = *params.begin();
    } else if (name == "-m" or name == "--matopt") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      matoptimize_path = *params.begin();
    } else if (name == "-c" or name == "--count") {
      if (params.empty()) {
        std::cerr << "Count not specified.\n";
        Fail();
      }
      count = static_cast<size_t>(ParseNumber(*params.begin()));
    } else if (name == "-s" or name == "--switch-subtrees") {
      if (params.empty()) {
        std::cerr << "Count not specified.\n";
        Fail();
      }
      switch_subtrees = static_cast<size_t>(ParseNumber(*params.begin()));
    } else if (name == "-l" or name == "--logpath") {
      if (params.empty()) {
        std::cerr << "log path name not specified.\n";
        Fail();
      }
      logfile_path = *params.begin();
    } else if (name == "--move-coeff-pscore") {
      if (params.empty()) {
        std::cerr << "parsimony score move coefficient not specified\n";
        Fail();
      }
      move_coeff_pscore = ParseNumber(*params.begin());
    } else if (name == "--min-subtree-clade-size") {
      if (params.empty()) {
        std::cerr << "minimum subtree clade size not specified\n";
        Fail();
      }
      min_subtree_clade_size = static_cast<size_t>(ParseNumber(*params.begin()));
    } else if (name == "--max-subtree-clade-size") {
      if (params.empty()) {
        std::cerr << "maximum subtree clade size not specified\n";
        Fail();
      }
      max_subtree_clade_size = static_cast<size_t>(ParseNumber(*params.begin()));
    } else if (name == "--move-coeff-nodes") {
      if (params.empty()) {
        std::cerr << "parsimony score move coefficient not specified\n";
        Fail();
      }
      move_coeff_nodes = ParseNumber(*params.begin());
    } else if (name == "--sample-any-tree") {
      sample_best_tree = false;
    } else if (name == "-r" or name == "--MAT-refseq-file") {
      if (params.empty()) {
        std::cerr << "Mutation annotated tree refsequence fasta path not specified.\n";
        Fail();
      }
      refseq_path = *params.begin();
    } else if (name == "--callback-option") {
      if (params.empty()) {
        std::cerr << "Callback configuration not specified.\n";
        Fail();
      }
      callback_config = *params.begin();
    } else {
      std::cerr << "Unknown argument.\n";
      Fail();
    }
  }

  if (input_dag_path.empty()) {
    std::cerr << "Path to input DAG not specified.\n";
    Fail();
  }

  if (output_dag_path.empty()) {
    std::cerr << "Path to output DAG not specified.\n";
    Fail();
  }
  std::filesystem::create_directory(logfile_path);
  std::string logfile_name = logfile_path + "/logfile.csv";

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &ignored);

  std::ofstream logfile;
  logfile.open(logfile_name);
  logfile << "Iteration\tNTrees\tNNodes\tNEdges\tMaxParsimony\tNTreesMaxParsimony\tWors"
             "tParsimony\tSecondsElapsed";

  MADAGStorage input_dag =
      refseq_path.empty()
          ? LoadDAGFromProtobuf(input_dag_path)
          : LoadTreeFromProtobuf(input_dag_path, LoadReferenceSequence(refseq_path));
  auto input_dag_view = input_dag.View();
  input_dag_view.RecomputeCompactGenomes();
  Merge merge{input_dag_view.GetReferenceSequence()};
  merge.AddDAG(input_dag_view);
  std::vector<std::pair<decltype(AddMATConversion(MADAGStorage{{}})), MAT::Tree>>
      optimized_dags;

  auto start_time = std::chrono::high_resolution_clock::now();
  auto time_elapsed = [&start_time]() {
    auto now = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
  };

  auto logger = [&merge, &logfile, &time_elapsed](size_t iteration) {
    SubtreeWeight<BinaryParsimonyScore, MergeDAG> parsimonyscorer{merge.GetResult()};
    SubtreeWeight<MaxBinaryParsimonyScore, MergeDAG> maxparsimonyscorer{
        merge.GetResult()};
    merge.ComputeResultEdgeMutations();
    auto minparsimony =
        parsimonyscorer.ComputeWeightBelow(merge.GetResult().GetRoot(), {});
    auto minparsimonytrees =
        parsimonyscorer.MinWeightCount(merge.GetResult().GetRoot(), {});
    auto maxparsimony =
        maxparsimonyscorer.ComputeWeightBelow(merge.GetResult().GetRoot(), {});
    SubtreeWeight<TreeCount, MergeDAG> treecount{merge.GetResult()};
    auto ntrees = treecount.ComputeWeightBelow(merge.GetResult().GetRoot(), {});
    std::cout << "Best parsimony score in DAG: " << minparsimony << "\n";
    std::cout << "Worst parsimony score in DAG: " << maxparsimony << "\n";
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Total trees in DAG: " << ntrees
              << "\n";
    std::cout << "Optimal trees in DAG: " << minparsimonytrees << "\n";
    logfile << '\n'
            << iteration << '\t' << ntrees << '\t' << merge.GetResult().GetNodesCount()
            << '\t' << merge.GetResult().GetEdgesCount() << '\t' << minparsimony << '\t'
            << minparsimonytrees << '\t' << maxparsimony << '\t' << time_elapsed()
            << std::flush;
  };
  logger(0);

  bool subtrees = false;
  for (size_t i = 0; i < count; ++i) {
    std::cout << "############ Beginning optimize loop " << std::to_string(i)
              << " #######\n";

    subtrees = (i >= switch_subtrees);
    merge.ComputeResultEdgeMutations();
    SubtreeWeight<BinaryParsimonyScore, MergeDAG> weight{merge.GetResult()};
    // choose root node for subtree (if sampling a subtree)
    auto subtree_node = [&]() -> std::optional<NodeId> {
      if (subtrees) {
        std::vector<NodeId> options;
        std::vector<size_t> weights;
        std::set<NodeId> visited;
        auto node_filter = [&](NodeId start_node, auto& node_filter_ref) -> void {
          auto node_instance = weight.GetDAG().Get(start_node);
          if (not visited.count(start_node)) {
            visited.insert(start_node);
            bool is_root = node_instance.IsUA();
            size_t clade_size = merge.GetResultNodeLabels()
                                    .at(node_instance.GetId())
                                    .GetLeafSet()
                                    ->ParentCladeSize();

            if (not is_root and clade_size < min_subtree_clade_size) {
              // terminate recursion because clade size only
              // decreases in children
              return;
            } else {
              // Add any other required conditions here
              if (not is_root and clade_size <= max_subtree_clade_size and
                  not node_instance.IsLeaf() and
                  not node_instance.GetCompactGenome().empty()) {
                if (uniform_subtree_root) {
                  weights.push_back(1);
                } else {
                  size_t min_parent_mutations = std::numeric_limits<size_t>::max();
                  for (auto parent_edge : node_instance.GetParents()) {
                    auto n_edge_muts = parent_edge.GetEdgeMutations().size();
                    if (min_parent_mutations > n_edge_muts) {
                      min_parent_mutations = n_edge_muts;
                    }
                  }
                  // increase preference for many parent mutations,
                  // and ensure weights cannot not be zero:
                  weights.push_back(1 + min_parent_mutations * min_parent_mutations);
                }
                options.push_back(start_node);
              }
              for (auto child_edge : node_instance.GetChildren()) {
                node_filter_ref(child_edge.GetChild().GetId(), node_filter_ref);
              }
            }
          }
        };
        node_filter(weight.GetDAG().GetRoot().GetId(), node_filter);
        if (not options.empty()) {
          std::random_device random_device;
          std::mt19937 random_generator(random_device());
          size_t option_idx = {std::discrete_distribution<size_t>{
              weights.begin(), weights.end()}(random_generator)};
          NodeId chosen_node = options.at(option_idx);
          std::cout << "Chose node " << chosen_node.value
                    << " as subtree root, with score " << weights.at(option_idx) << "\n"
                    << std::flush;
          return weight.GetDAG().Get(chosen_node);
        } else {
          std::cout << "Warning: No suitable subtree root nodes found. Optimizing an "
                       "entire tree.\n"
                    << std::flush;
        }
      }
      return std::nullopt;
    }();

    auto sample = sample_best_tree
                      ? AddMATConversion(weight.MinWeightSampleTree({}, subtree_node))
                      : AddMATConversion(weight.SampleTree({}, subtree_node));
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>> Nodes in sampled (sub)tree: "
              << sample.GetNodesCount() << "\n";
    MAT::Tree mat;
    sample.View().BuildMAT(mat);
    sample.View().GetRoot().Validate(true);
    check_edge_mutations(sample.View().Const());

    if (callback_config == "all-moves") {
      Merge_All_Moves_Found_Callback callback{merge, sample.View()};
      optimized_dags.push_back(
          optimize_dag_direct(sample.View(), callback, callback, callback));
    } else if (callback_config == "best-moves-fixed-tree") {
      Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback callback{
          merge, sample.View(), {move_coeff_nodes, move_coeff_pscore}};
      optimized_dags.push_back(
          optimize_dag_direct(sample.View(), callback, callback, callback));
    } else if (callback_config == "best-moves-treebased") {
      Treebased_Move_Found_Callback callback{
          merge, sample.View(), {move_coeff_nodes, move_coeff_pscore}};
      optimized_dags.push_back(
          optimize_dag_direct(sample.View(), callback, callback, callback));
    } else {
      Merge_All_Profitable_Moves_Found_Callback callback{
          merge, sample.View(), {move_coeff_nodes, move_coeff_pscore}};
      optimized_dags.push_back(
          optimize_dag_direct(sample.View(), callback, callback, callback));
    }

    auto optimized_view = optimized_dags.back().first.View();
    optimized_view.RecomputeCompactGenomes();
    merge.AddDAG(optimized_view);
    logger(i + 1);
  }

  std::cout << "new node coefficient: " << move_coeff_nodes << "\n";
  std::cout << "parsimony score coefficient: " << move_coeff_pscore << "\n";
  logfile.close();
  StoreDAGToProtobuf(merge.GetResult(), output_dag_path);

  return EXIT_SUCCESS;
}
