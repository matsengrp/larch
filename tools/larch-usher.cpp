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
#include "larch/merge/merge.hpp"
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

// NOLINTNEXTLINE(cppcoreguidelines-virtual-class-destructor)
template <typename SampleDAG>
struct Larch_Move_Found_Callback : public Move_Found_Callback {
  Larch_Move_Found_Callback(const Merge<MADAG>& merge, SampleDAG sample,
                            const std::vector<NodeId>& sample_dag_ids)
      : merge_{merge},
        sample_{sample},
        sample_dag_ids_{sample_dag_ids},
        move_score_coeffs_{1, 1} {}
  Larch_Move_Found_Callback(
      const Merge<MADAG>& merge, SampleDAG sample,
      const std::vector<NodeId>& sample_dag_ids,
      std::pair<int, int> move_score_coeffs)  // NOLINT(modernize-pass-by-value)
      : merge_{merge},
        sample_{sample},
        sample_dag_ids_{sample_dag_ids},
        move_score_coeffs_{move_score_coeffs} {}
  bool operator()(Profitable_Moves& move, int /* best_score_change */,
                  std::vector<Node_With_Major_Allele_Set_Change>&
                  /* node_with_major_allele_set_change */) override {
    int node_id_map_count = 0;
    if (move_score_coeffs_.first != 0) {
      NodeId src_id = ToMergedNodeId(move.src->node_id);
      NodeId dst_id = ToMergedNodeId(move.dst->node_id);
      NodeId lca_id = ToMergedNodeId(move.LCA->node_id);

      const auto& src_clades =
          merge_.GetResultNodeLabels().at(src_id.value).GetLeafSet()->GetClades();
      const auto& dst_clades =
          merge_.GetResultNodeLabels().at(dst_id.value).GetLeafSet()->GetClades();

      MAT::Node* curr_node = move.src;
      while (not(curr_node->node_id == lca_id.value)) {
        MergeDAG::NodeView node = merge_.GetResult().Get(NodeId{curr_node->node_id});
        const auto& clades = merge_.GetResultNodeLabels()
                                 .at(node.GetId().value)
                                 .GetLeafSet()
                                 ->GetClades();
        if (not merge_.ContainsLeafset(clades_difference(clades, src_clades))) {
          ++node_id_map_count;
        }
        curr_node = curr_node->parent;
        if (curr_node == nullptr) {
          break;
        }
      }

      curr_node = move.dst;
      while (not(curr_node->node_id == lca_id.value)) {
        MergeDAG::NodeView node = merge_.GetResult().Get(NodeId{curr_node->node_id});
        const auto& clades = merge_.GetResultNodeLabels()
                                 .at(node.GetId().value)
                                 .GetLeafSet()
                                 ->GetClades();
        if (not merge_.ContainsLeafset(clades_union(clades, dst_clades))) {
          ++node_id_map_count;
        }
        curr_node = curr_node->parent;
        if (curr_node == nullptr) {
          break;
        }
      }
    }

    move.score_change = move_score_coeffs_.second * move.score_change -
                        move_score_coeffs_.first * node_id_map_count;
    return move.score_change <= 0;
  }

  void MergeNodeIDs(std::map<NodeId, NodeId>&& node_id_map) {
    node_id_map_.merge(std::forward<decltype(node_id_map)>(node_id_map));
  }

 private:
  NodeId ToMergedNodeId(size_t id) {
    auto it = node_id_map_.find(NodeId{id});
    if (it != node_id_map_.end()) {
      return it->second;
    }
    if (id < sample_dag_ids_.size()) {
      return sample_dag_ids_.at(id);
    }
    return {};
  }

  const Merge<MADAG>& merge_;
  SampleDAG sample_;
  const std::vector<NodeId>& sample_dag_ids_;
  const std::pair<int, int> move_score_coeffs_;
  std::map<NodeId, NodeId> node_id_map_;
};

int main(int argc, char** argv) {  // NOLINT(bugprone-exception-escape)
  Arguments args = GetArguments(argc, argv);
  int ignored{};
  std::string input_dag_path;
  std::string output_dag_path;
  std::string matoptimize_path = "matOptimize";
  std::string logfile_path = "optimization_log";
  std::string refseq_path;
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

  input_dag.View().RecomputeCompactGenomes();
  Merge<MADAG> merge{input_dag.View().GetReferenceSequence()};
  merge.AddDAGs({input_dag.View()});
  std::vector<MADAGStorage> optimized_dags;

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
    std::optional<NodeId> subtree_node;
    auto [sample, dag_ids] = [&] {
      subtree_node = [&]() -> std::optional<NodeId> {
        if (subtrees) {
          std::vector<NodeId> options;
          std::vector<size_t> weights;
          std::set<NodeId> visited;
          auto node_filter = [&](NodeId start_node, auto& node_filter_ref) -> void {
            auto node_instance = weight.GetDAG().Get(start_node);
            if (not visited.count(start_node)) {
              visited.insert(start_node);
              bool is_root = node_instance.IsRoot();
              size_t clade_size = merge.GetResultNodeLabels()
                                      .at(node_instance.GetId().value)
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
                      << " as subtree root, with score " << weights.at(option_idx)
                      << "\n"
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

      if (sample_best_tree) {
        return weight.MinWeightSampleTree({}, subtree_node);
      } else {
        return weight.SampleTree({}, subtree_node);
      }
    }();
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>> Nodes in sampled (sub)tree: "
              << sample.GetNodesCount() << "\n";
    check_edge_mutations(sample.View());
    Larch_Move_Found_Callback callback{
        merge, sample.View(), dag_ids, {move_coeff_nodes, move_coeff_pscore}};
    /* StoreTreeToProtobuf(sample.View(), "before_optimize_dag.pb"); */
    auto radius_callback = [&](MAT::Tree& tree) -> void {
      auto [result, mat_node_map] =
          build_madag_from_mat(tree, merge.GetResult().GetReferenceSequence());
      result.View().RecomputeCompactGenomes();
      optimized_dags.push_back(std::move(result));
      std::map<NodeId, NodeId> full_map = [&, &mat_node_map = mat_node_map] {
        std::map<NodeId, NodeId> merge_node_map;
        if (subtrees) {
          merge_node_map = merge.AddDAG(optimized_dags.back().View(),
                                        merge.GetResult().Get(subtree_node.value()));
        } else {
          merge_node_map = merge.AddDAG(optimized_dags.back().View());
        }
        // mat_node_map is not the identity, so all pairs in mat_node_map must be used
        // to build remaped
        std::map<NodeId, NodeId> remaped;
        for (auto [from, to] : mat_node_map) {
          remaped.insert({to, merge_node_map.at(from)});
        }
        return remaped;
      }();
      callback.MergeNodeIDs(std::move(full_map));
    };
    optimize_dag_direct(sample.View(), callback, radius_callback);

    if (i % 10 == 0) {  // NOLINT
      StoreDAGToProtobuf(merge.GetResult(), logfile_path + "/intermediate_dag.pb");
    }
    logger(i + 1);
  }

  std::cout << "new node coefficient: " << move_coeff_nodes << "\n";
  std::cout << "parsimony score coefficient: " << move_coeff_pscore << "\n";
  logfile.close();
  StoreDAGToProtobuf(merge.GetResult(), output_dag_path);

  return EXIT_SUCCESS;
}
