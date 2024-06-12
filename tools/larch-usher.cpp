#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

#include <unistd.h>
#include <sys/wait.h>

#include "arguments.hpp"
#include "larch/dag_loader.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/weight_accumulator.hpp"
#include "larch/subtree/tree_count.hpp"
#include "larch/subtree/parsimony_score_binary.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/rf_distance.hpp"
#include "larch/spr/spr_view.hpp"
#include "larch/merge/merge.hpp"
#include "larch/merge/leaf_set.hpp"
#include "larch/spr/batching_callback.hpp"
#include "larch/benchmark.hpp"
#include <mpi.h>

#include "larch/usher_glue.hpp"

#include <tbb/global_control.h>

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout
      << "larch-usher -i,--input infile -o,--output outfile [-c,--count number]\n";
  std::cout << "  -i,--input   Path to input Tree/DAG\n";
  std::cout << "  -r,--MAT-refseq-file   Provide a path to a file containing a "
               "reference sequence\nif input points to MAT protobuf\n";
  std::cout << "  -v,--VCF-input-file   Provide a path to a vcf file containing "
               "ambiguous leaf sequence data\n";
  std::cout << "  -o,--output  Path to output DAG\n";
  std::cout << "  -l,--logpath Path for logging\n";
  std::cout << "  -c,--count   Number of iterations. (default: 1)\n";
  std::cout << "  -s,--switch-subtrees          Switch to optimizing subtrees after "
               "the specified number of iterations (default: never)\n";
  std::cout << "  --min-subtree-clade-size      The minimum number of leaves in a "
               "subtree sampled for optimization (default: 100, ignored without option "
               "`-s`)\n";
  std::cout
      << "  --max-subtree-clade-size      The maximum number of leaves in a subtree "
         "sampled for optimization (default: 1000, ignored without option "
         "`-s`)\n";
  std::cout << "  --move-coeff-nodes   New node coefficient for scoring moves. "
               "(default: 1)\n";
  std::cout << "  --move-coeff-pscore  Parsimony score coefficient for scoring moves. "
               "(default: 1)\n";
  std::cout << "  --sample-any-tree    Sample any tree for optimization, rather than "
               "requiring the sampled tree to maximize parsimony.\n";
  std::cout << "  --sample-method      Select method for optimization "
               "[parsimony, random, rf-minsum, rf-maxsum]\n";
  std::cout
      << "  --callback-option    Callback configuration choice (default: merge all "
         "profitable moves) [best-move, best-move-fixed-tree, best-move-treebased, "
         "all-moves]\n";
  std::cout << "  --trim   Trim optimized dag before writing to protobuf\n";
  std::cout
      << "  --keep-fragment-uncollapsed   Optional argument to keep empty fragment "
         "edges\n";
  std::cout << "  --input-format  Specify input file format. (default: inferred) "
               "[dagbin, dag-pb, tree-pb, dag-json]\n";
  std::cout << "  --output-format  Specify output file format. (default: inferred) "
               "[dagbin, dag-pb]\n";
  std::cout << "  --seed     Set seed for random number generation.\n";
  std::cout << "  --thread   Set number of cpu threads.\n";

  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";
  std::exit(EXIT_FAILURE);
}

void check_edge_mutations(MADAG madag);

std::vector<std::vector<const SampleId*>> clades_union(
    const std::vector<std::vector<const SampleId*>>& lhs,
    const std::vector<std::vector<const SampleId*>>& rhs) {
  std::vector<std::vector<const SampleId*>> result;

  for (auto [lhs_clade, rhs_clade] : ranges::views::zip(lhs, rhs)) {
    std::vector<const SampleId*> clade{lhs_clade};
    clade.insert(clade.end(), rhs_clade.begin(), rhs_clade.end());
    ranges::sort(clade);
    ranges::unique(clade);
    result.push_back(std::move(clade));
  }

  ranges::sort(result);
  return result;
}

std::vector<std::vector<const SampleId*>> clades_difference(
    const std::vector<std::vector<const SampleId*>>& lhs,
    const std::vector<std::vector<const SampleId*>>& rhs) {
  std::vector<std::vector<const SampleId*>> result;

  for (auto [lhs_clade, rhs_clade] : ranges::views::zip(lhs, rhs)) {
    std::vector<const SampleId*> clade;
    std::set_difference(lhs_clade.begin(), lhs_clade.end(), rhs_clade.begin(),
                        rhs_clade.end(), std::inserter(clade, clade.begin()));
    ranges::sort(clade);
    ranges::unique(clade);
    result.push_back(std::move(clade));
  }

  ranges::sort(result);
  return result;
}

using Storage = MergeDAGStorage<>;
using ReassignedStatesStorage =
    decltype(AddMappedNodes(AddMATConversion(Storage::EmptyDefault())));

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

  Treebased_Move_Found_Callback(Merge& merge, SampleDAG sample_dag,
                                std::pair<int, int> move_score_coeffs,
                                bool collapse_empty_fragment_edges)
      : BatchingCallback<Treebased_Move_Found_Callback<SampleDAG>,
                         SampleDAG>{merge, sample_dag, collapse_empty_fragment_edges},
        move_score_coeffs_{std::move(move_score_coeffs)} {};

  Treebased_Move_Found_Callback(Merge& merge, SampleDAG sample_dag,
                                bool collapse_empty_fragment_edges)
      : BatchingCallback<Treebased_Move_Found_Callback<SampleDAG>,
                         SampleDAG>{merge, sample_dag, collapse_empty_fragment_edges},
        move_score_coeffs_{1, 1} {};

  template <typename SPRView, typename FragmentType>
  std::pair<bool, bool> OnMove(SPRView spr, const FragmentType& fragment,
                               Profitable_Moves& move, int /*best_score_change*/,
                               std::vector<Node_With_Major_Allele_Set_Change>&
                               /*nodes_with_major_allele_set_change*/) {
    int node_id_map_count = 0;
    if (move_score_coeffs_.first != 0) {
      auto make_leaf_set = [&](std::vector<NodeId> leaf_node_ids) {
        std::vector<const SampleId*> ls;
        for (auto leaf_node : leaf_node_ids) {
          auto sid = SampleId{spr.Const().Get(leaf_node).GetSampleId().value()};
          ls.push_back(&sid);
        }
        ranges::sort(ls);
        ranges::unique(ls);
        std::vector<std::vector<const SampleId*>> to_ret;
        to_ret.push_back(std::move(ls));
        return to_ret;
      };
      auto src_leaf_set =
          spr.GetMoveSource().GetOld().IsCondensedInMAT()
              ? make_leaf_set(spr.GetMoveSources())
              : this->GetMerge()
                    .GetResultNodeLabels()
                    .at(this->GetMappedStorage()
                            .GetNodeFromMAT(spr.GetMoveSource().GetOld().GetMATNode())
                            .GetOriginalId())
                    .GetLeafSet()
                    ->GetClades();
      auto dst_leaf_set =
          spr.GetMoveTarget().GetOld().IsCondensedInMAT()
              ? make_leaf_set(spr.GetMoveTargets())
              : this->GetMerge()
                    .GetResultNodeLabels()
                    .at(this->GetMappedStorage()
                            .GetNodeFromMAT(spr.GetMoveTarget().GetOld().GetMATNode())
                            .GetOriginalId())
                    .GetLeafSet()
                    ->GetClades();
      for (auto hypothetical_node : fragment.GetNodes()) {
        if (hypothetical_node.IsMoveNew()) {
          if (not(this->GetMerge().ContainsLeafset(
                  clades_union(src_leaf_set, dst_leaf_set)))) {
            ++node_id_map_count;
          }
        } else {
          if (not hypothetical_node.GetOld().IsLeaf()) {
            if (hypothetical_node.HasChangedTopology() and
                hypothetical_node.GetOld().HaveMATNode()) {
              if (not(hypothetical_node.GetOld().GetMATNode()->node_id ==
                          move.src->node_id or
                      hypothetical_node.GetOld().GetMATNode()->node_id ==
                          move.dst->node_id or
                      hypothetical_node.GetOld().GetMATNode()->is_root())) {
                const auto& current_leaf_sets =
                    this->GetMerge()
                        .GetResultNodeLabels()
                        .at(this->GetMappedStorage()
                                .GetNodeFromMAT(hypothetical_node.GetOld().GetMATNode())
                                .GetOriginalId())
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
        }
      }
    }
    move.score_change = move_score_coeffs_.second * move.score_change -
                        move_score_coeffs_.first * node_id_map_count;
    return {false, move.score_change <= 0};
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

  Merge_All_Moves_Found_Callback(Merge& merge, SampleDAG sample_dag,
                                 bool collapse_empty_fragment_edges)
      : BatchingCallback<Merge_All_Moves_Found_Callback<SampleDAG>, SampleDAG>{
            merge, sample_dag, collapse_empty_fragment_edges} {};

  template <typename SPRView, typename FragmentType>
  std::pair<bool, bool> OnMove(SPRView /*spr*/, const FragmentType& /*fragment*/,
                               Profitable_Moves& move, int best_score_change,
                               std::vector<Node_With_Major_Allele_Set_Change>&
                               /*nodes_with_major_allele_set_change*/) {
    std::ignore = move;
    std::ignore = best_score_change;
    return {true, true};
  }

  void OnRadius(){};
};

template <typename SampleDAG>
struct Merge_All_Profitable_Moves_Found_Callback
    : public BatchingCallback<Merge_All_Profitable_Moves_Found_Callback<SampleDAG>,
                              SampleDAG> {
  MOVE_ONLY_VIRT_DTOR(Merge_All_Profitable_Moves_Found_Callback);

  Merge_All_Profitable_Moves_Found_Callback(Merge& merge, SampleDAG sample_dag,
                                            std::pair<int, int> move_score_coeffs)
      : BatchingCallback<Merge_All_Profitable_Moves_Found_Callback<SampleDAG>,
                         SampleDAG>{merge, sample_dag},
        sample_dag_{sample_dag},
        merge_{merge},
        move_score_coeffs_{std::move(move_score_coeffs)} {};

  Merge_All_Profitable_Moves_Found_Callback(Merge& merge, SampleDAG sample_dag)
      : BatchingCallback<Merge_All_Profitable_Moves_Found_Callback<SampleDAG>,
                         SampleDAG>{merge, sample_dag},
        sample_dag_{sample_dag},
        merge_{merge},
        move_score_coeffs_{1, 1} {};

  Merge_All_Profitable_Moves_Found_Callback(Merge& merge, SampleDAG sample_dag,
                                            std::pair<int, int> move_score_coeffs,
                                            bool collapse_empty_fragment_edges)
      : BatchingCallback<Merge_All_Profitable_Moves_Found_Callback<SampleDAG>,
                         SampleDAG>{merge, sample_dag, collapse_empty_fragment_edges},
        sample_dag_{sample_dag},
        merge_{merge},
        move_score_coeffs_{std::move(move_score_coeffs)} {};

  Merge_All_Profitable_Moves_Found_Callback(Merge& merge, SampleDAG sample_dag,
                                            bool collapse_empty_fragment_edges)
      : BatchingCallback<Merge_All_Profitable_Moves_Found_Callback<SampleDAG>,
                         SampleDAG>{merge, sample_dag, collapse_empty_fragment_edges},
        sample_dag_{sample_dag},
        merge_{merge},
        move_score_coeffs_{1, 1} {};

  template <typename SPRView, typename FragmentType>
  std::pair<bool, bool> OnMove(SPRView spr, const FragmentType& fragment,
                               Profitable_Moves& move, int /*best_score_change*/,
                               std::vector<Node_With_Major_Allele_Set_Change>&
                               /*nodes_with_major_allele_set_change*/) {
    int node_id_map_count = 0;
    if (move_score_coeffs_.first != 0) {
      auto make_leaf_set = [&](std::vector<NodeId> leaf_node_ids) {
        std::vector<const SampleId*> ls;
        for (auto leaf_node : leaf_node_ids) {
          auto sid = SampleId{spr.Const().Get(leaf_node).GetSampleId().value()};
          ls.push_back(&sid);
        }
        ranges::sort(ls);
        ranges::unique(ls);
        std::vector<std::vector<const SampleId*>> to_ret;
        to_ret.push_back(std::move(ls));
        return to_ret;
      };
      auto src_leaf_set =
          spr.GetMoveSource().GetOld().IsCondensedInMAT()
              ? make_leaf_set(spr.GetMoveSources())
              : this->GetMerge()
                    .GetResultNodeLabels()
                    .at(this->GetMappedStorage()
                            .GetNodeFromMAT(spr.GetMoveSource().GetOld().GetMATNode())
                            .GetOriginalId())
                    .GetLeafSet()
                    ->GetClades();
      auto dst_leaf_set =
          spr.GetMoveTarget().GetOld().IsCondensedInMAT()
              ? make_leaf_set(spr.GetMoveTargets())
              : this->GetMerge()
                    .GetResultNodeLabels()
                    .at(this->GetMappedStorage()
                            .GetNodeFromMAT(spr.GetMoveTarget().GetOld().GetMATNode())
                            .GetOriginalId())
                    .GetLeafSet()
                    ->GetClades();
      for (auto hypothetical_node : fragment.GetNodes()) {
        if (hypothetical_node.IsMoveNew()) {
          if (not(this->GetMerge().ContainsLeafset(
                  clades_union(src_leaf_set, dst_leaf_set)))) {
            ++node_id_map_count;
          }
        } else {
          if (not hypothetical_node.GetOld().IsLeaf()) {
            if (hypothetical_node.HasChangedTopology() and
                hypothetical_node.GetOld().HaveMATNode()) {
              if (not(hypothetical_node.GetOld().GetMATNode()->node_id ==
                          move.src->node_id or
                      hypothetical_node.GetOld().GetMATNode()->node_id ==
                          move.dst->node_id or
                      hypothetical_node.GetOld().GetMATNode()->is_root())) {
                const auto& current_leaf_sets =
                    this->GetMerge()
                        .GetResultNodeLabels()
                        .at(this->GetMappedStorage()
                                .GetNodeFromMAT(hypothetical_node.GetOld().GetMATNode())
                                .GetOriginalId())
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
        }
      }
    }
    move.score_change = move_score_coeffs_.second * move.score_change -
                        move_score_coeffs_.first * node_id_map_count;
    return {move.score_change <= 0, move.score_change <= 0};
  }

  void OnRadius(){};

  SampleDAG sample_dag_;
  Merge& merge_;
  ReassignedStatesStorage reassigned_states_storage_ =
      AddMappedNodes(AddMATConversion(Storage::EmptyDefault()));
  std::atomic<MAT::Tree*> sample_mat_ = nullptr;
  std::mutex merge_mtx_;
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
        sample_dag_{sample_dag},
        merge_{merge},
        move_score_coeffs_{std::move(move_score_coeffs)} {};

  Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback(Merge& merge,
                                                       SampleDAG sample_dag)
      : BatchingCallback<
            Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback<SampleDAG>,
            SampleDAG>{merge, sample_dag},
        sample_dag_{sample_dag},
        merge_{merge},
        move_score_coeffs_{1, 1} {};

  Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback(
      Merge& merge, SampleDAG sample_dag, std::pair<int, int> move_score_coeffs,
      bool collapse_empty_fragment_edges)
      : BatchingCallback<
            Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback<SampleDAG>,
            SampleDAG>{merge, sample_dag, collapse_empty_fragment_edges},
        sample_dag_{sample_dag},
        merge_{merge},
        move_score_coeffs_{std::move(move_score_coeffs)} {};

  Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback(
      Merge& merge, SampleDAG sample_dag, bool collapse_empty_fragment_edges)
      : BatchingCallback<
            Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback<SampleDAG>,
            SampleDAG>{merge, sample_dag, collapse_empty_fragment_edges},
        sample_dag_{sample_dag},
        merge_{merge},
        move_score_coeffs_{1, 1} {};

  template <typename SPRView, typename FragmentType>
  std::pair<bool, bool> OnMove(SPRView spr, const FragmentType& fragment,
                               Profitable_Moves& move, int /*best_score_change*/,
                               std::vector<Node_With_Major_Allele_Set_Change>&
                               /*nodes_with_major_allele_set_change*/) {
    int node_id_map_count = 0;
    if (move_score_coeffs_.first != 0) {
      auto make_leaf_set = [&](std::vector<NodeId> leaf_node_ids) {
        std::vector<const SampleId*> ls;
        for (auto leaf_node : leaf_node_ids) {
          auto sid = SampleId{spr.Const().Get(leaf_node).GetSampleId().value()};
          ls.push_back(&sid);
        }
        ranges::sort(ls);
        ranges::unique(ls);
        std::vector<std::vector<const SampleId*>> to_ret;
        to_ret.push_back(std::move(ls));
        return to_ret;
      };
      auto src_leaf_set =
          spr.GetMoveSource().GetOld().IsCondensedInMAT()
              ? make_leaf_set(spr.GetMoveSources())
              : this->GetMerge()
                    .GetResultNodeLabels()
                    .at(this->GetMappedStorage()
                            .GetNodeFromMAT(spr.GetMoveSource().GetOld().GetMATNode())
                            .GetOriginalId())
                    .GetLeafSet()
                    ->GetClades();
      auto dst_leaf_set =
          spr.GetMoveTarget().GetOld().IsCondensedInMAT()
              ? make_leaf_set(spr.GetMoveTargets())
              : this->GetMerge()
                    .GetResultNodeLabels()
                    .at(this->GetMappedStorage()
                            .GetNodeFromMAT(spr.GetMoveTarget().GetOld().GetMATNode())
                            .GetOriginalId())
                    .GetLeafSet()
                    ->GetClades();
      for (auto hypothetical_node : fragment.GetNodes()) {
        if (hypothetical_node.IsMoveNew()) {
          if (not(this->GetMerge().ContainsLeafset(
                  clades_union(src_leaf_set, dst_leaf_set)))) {
            ++node_id_map_count;
          }
        } else {
          if (not hypothetical_node.GetOld().IsLeaf()) {
            if (hypothetical_node.HasChangedTopology() and
                hypothetical_node.GetOld().HaveMATNode()) {
              if (not(hypothetical_node.GetOld().GetMATNode()->node_id ==
                          move.src->node_id or
                      hypothetical_node.GetOld().GetMATNode()->node_id ==
                          move.dst->node_id or
                      hypothetical_node.GetOld().GetMATNode()->is_root())) {
                const auto& current_leaf_sets =
                    this->GetMerge()
                        .GetResultNodeLabels()
                        .at(this->GetMappedStorage()
                                .GetNodeFromMAT(hypothetical_node.GetOld().GetMATNode())
                                .GetOriginalId())
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
        }
      }
    }
    move.score_change = move_score_coeffs_.second * move.score_change -
                        move_score_coeffs_.first * node_id_map_count;
    return {move.score_change <= 0, false};
  }

  void OnRadius(){};

  SampleDAG sample_dag_;
  Merge& merge_;
  ReassignedStatesStorage reassigned_states_storage_ =
      AddMappedNodes(AddMATConversion(Storage::EmptyDefault()));
  std::atomic<MAT::Tree*> sample_mat_ = nullptr;
  std::mutex merge_mtx_;
  std::pair<int, int> move_score_coeffs_;
};

enum class SampleMethod {
  Random,
  UniformRandom,
  Parsimony,
  UniformParsimony,
  MinSumRFDistance,
  MaxSumRFDistance
};

enum class CallbackMethod {
  BestMoves,
  BestMovesFixedTree,
  BestMovesTreeBased,
  AllMoves
};

int main(int argc, char** argv) {  // NOLINT(bugprone-exception-escape)
  Arguments args = GetArguments(argc, argv);
  int ignored{};
  std::string input_dag_path;
  std::string output_dag_path;
  std::string intermediate_dag_path;
  std::string logfile_path = "optimization_log";
  std::string refseq_path;
  std::string vcf_path;
  CallbackMethod callback_config = CallbackMethod::BestMoves;
  bool write_intermediate_dag = true;
  FileFormat input_format = FileFormat::Infer;
  FileFormat output_format = FileFormat::Infer;
  SampleMethod sample_method = SampleMethod::Parsimony;
  bool sample_uniformly = false;
  size_t count = 1;
  unsigned int thread_count = 0;
  int move_coeff_nodes = 1;
  int move_coeff_pscore = 1;
  size_t switch_subtrees = std::numeric_limits<size_t>::max();
  size_t min_subtree_clade_size = 100;   // NOLINT
  size_t max_subtree_clade_size = 1000;  // NOLINT
  bool uniform_subtree_root = false;
  bool collapse_empty_fragment_edges = true;
  bool final_trim = false;
  std::optional<uint32_t> user_seed = std::nullopt;

  Benchmark total_timer;

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
    } else if (name == "-c" or name == "--count") {
      if (params.empty()) {
        std::cerr << "Count not specified.\n";
        Fail();
      }
      count = static_cast<size_t>(ParseNumber(*params.begin()));
    } else if (name == "-s" or name == "--switch-subtrees") {
      if (params.empty()) {
        std::cerr << "Subtree count not specified.\n";
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
    } else if (name == "--sample-method") {
      if (params.empty()) {
        std::cerr << "sample-method not specified\n";
        Fail();
      }
      auto arg = *params.begin();
      if (arg == "random") {
        sample_method = SampleMethod::Random;
      } else if (arg == "parsimony") {
        sample_method = SampleMethod::Parsimony;
      } else if (arg == "rf-minsum") {
        sample_method = SampleMethod::MinSumRFDistance;
      } else if (arg == "rf-maxsum") {
        sample_method = SampleMethod::MaxSumRFDistance;
      } else {
        std::cerr << "ERROR: Unknown --sample-method argument '" << arg << "'"
                  << std::endl;
        Fail();
      }
    } else if (name == "--sample-any-tree") {
      sample_method = SampleMethod::Random;
    } else if (name == "--sample-uniformly") {
      sample_uniformly = true;
    } else if (name == "-r" or name == "--MAT-refseq-file") {
      if (params.empty()) {
        std::cerr << "Mutation annotated tree refsequence fasta path not specified.\n";
        Fail();
      }
      refseq_path = *params.begin();
    } else if (name == "-v" or name == "--VCF-input-file") {
      if (params.empty()) {
        std::cerr << "VCF file path not specified.\n";
        Fail();
      }
      vcf_path = *params.begin();
    } else if (name == "--callback-option") {
      if (params.empty()) {
        std::cerr << "callback-option not specified\n";
        Fail();
      }
      auto arg = *params.begin();
      if (arg == "best-moves") {
        callback_config = CallbackMethod::BestMoves;
      } else if (arg == "best-moves-fixed-tree") {
        callback_config = CallbackMethod::BestMovesFixedTree;
      } else if (arg == "best-moves-treebased") {
        callback_config = CallbackMethod::BestMovesTreeBased;
      } else if (arg == "all-moves") {
        callback_config = CallbackMethod::AllMoves;
      } else {
        std::cerr << "ERROR: Unknown --callback-option argument '" << arg << "'"
                  << std::endl;
        Fail();
      }
    } else if (name == "--keep-fragment-uncollapsed") {
      collapse_empty_fragment_edges = false;
    } else if (name == "--quiet") {
      write_intermediate_dag = false;
    } else if (name == "--trim") {
      final_trim = true;
    } else if (name == "--input-format") {
      if (params.empty()) {
        std::cerr << "File format not specified.\n";
        Fail();
      }
      std::string ext{*params.begin()};
      input_format = InferFileFormat(ext);
    } else if (name == "--output-format") {
      if (params.empty()) {
        std::cerr << "File format not specified.\n";
        Fail();
      }
      std::string ext{*params.begin()};
      output_format = InferFileFormat(ext);
    } else if (name == "--seed") {
      if (params.empty()) {
        std::cerr << "User seed not specified.\n";
        Fail();
      }
      user_seed = static_cast<uint32_t>(ParseNumber(*params.begin()));
    } else if (name == "--thread") {
      if (params.empty()) {
        std::cerr << "Thread count not specified.\n";
        Fail();
      }
      thread_count = static_cast<unsigned int>(ParseNumber(*params.begin()));
    } else {
      std::cerr << "Unknown argument '" << name << "'.\n";
      Fail();
    }
  }

  if (sample_uniformly) {
    if (sample_method == SampleMethod::Random) {
      sample_method = SampleMethod::UniformRandom;
    } else if (sample_method == SampleMethod::Parsimony) {
      sample_method = SampleMethod::UniformParsimony;
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

  bool is_input_dag = refseq_path.empty();
  if (input_format == FileFormat::Infer) {
    input_format = InferFileFormat(input_dag_path);
    if (input_format == FileFormat::Protobuf) {
      input_format = is_input_dag ? FileFormat::ProtobufDAG : FileFormat::ProtobufTree;
    }
  }
  if (output_format == FileFormat::Infer) {
    output_format = InferFileFormat(output_dag_path);
    if (output_format == FileFormat::Protobuf) {
      output_format = FileFormat::ProtobufDAG;
    }
  }
  if (intermediate_dag_path.empty()) {
    intermediate_dag_path = output_dag_path + ".intermediate";
  }

  RandomNumberGenerator main_rng{user_seed};
  std::cout << "Random seed: " << main_rng.seed_ << "\n";

  std::filesystem::create_directory(logfile_path);
  std::string logfile_name = logfile_path + "/logfile.csv";

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &ignored);
  tbb::global_control c(
      tbb::global_control::max_allowed_parallelism,
      ((thread_count > 0) ? std::min(thread_count, std::thread::hardware_concurrency())
                          : std::thread::hardware_concurrency()));
  std::cout << "Thread count: "
            << c.active_value(tbb::global_control::max_allowed_parallelism) << "\n";

  std::ofstream logfile;
  logfile.open(logfile_name);
  logfile << "Iteration\tNTrees\tNNodes\tNEdges\tMaxParsimony\tNTreesMaxParsimony\tWors"
             "tParsimony\tMinSumRFDistance\tMaxSumRFDistance\tMinSumRFCount\tMaxSumRFCo"
             "unt\tSecondsElapsed";

  Benchmark load_timer;
  std::cout << "Loading input DAG..." << std::flush;
  MADAGStorage<> input_dag = LoadDAG(input_dag_path, input_format, refseq_path);
  load_timer.stop();
  std::cout << "...loaded: " << load_timer.durationFormatMs() << std::endl;

  input_dag.View().RecomputeCompactGenomes(true);
  LoadVCFData(input_dag, vcf_path);
  auto input_dag_view = input_dag.View();
  // if the DAG is from a DAG protobuf file, then it needs to be equipped with
  // SampleIds
  if (vcf_path.empty()) {
    input_dag_view.SampleIdsFromCG();
  }
  Merge merge{input_dag_view.GetReferenceSequence()};
  merge.AddDAG(input_dag_view);
  std::vector<
      std::pair<decltype(AddMATConversion(MADAGStorage<>::EmptyDefault())), MAT::Tree>>
      optimized_dags;
  merge.ComputeResultEdgeMutations();

  Benchmark log_timer;
  auto logger = [&merge, &logfile, &log_timer, &intermediate_dag_path,
                 &write_intermediate_dag, &output_format, &main_rng](size_t iteration) {
    SubtreeWeight<BinaryParsimonyScore, MergeDAG> parsimonyscorer{
        merge.GetResult(), main_rng.GenerateSeed()};
    SubtreeWeight<MaxBinaryParsimonyScore, MergeDAG> maxparsimonyscorer{
        merge.GetResult(), main_rng.GenerateSeed()};
    merge.ComputeResultEdgeMutations();
    auto minparsimony =
        parsimonyscorer.ComputeWeightBelow(merge.GetResult().GetRoot(), {});
    auto minparsimonytrees =
        parsimonyscorer.MinWeightCount(merge.GetResult().GetRoot(), {});
    auto maxparsimony =
        maxparsimonyscorer.ComputeWeightBelow(merge.GetResult().GetRoot(), {});
    SubtreeWeight<TreeCount, MergeDAG> treecount{merge.GetResult(),
                                                 main_rng.GenerateSeed()};
    auto ntrees = treecount.ComputeWeightBelow(merge.GetResult().GetRoot(), {});

    SubtreeWeight<SumRFDistance, MergeDAG> this_min_sum_rf_dist{
        merge.GetResult(), main_rng.GenerateSeed()};
    SumRFDistance this_min_rf_weight_ops{merge, merge, main_rng.GenerateSeed()};
    SubtreeWeight<MaxSumRFDistance, MergeDAG> this_max_sum_rf_dist{
        merge.GetResult(), main_rng.GenerateSeed()};
    MaxSumRFDistance this_max_rf_weight_ops{merge, merge, main_rng.GenerateSeed()};
    auto shiftsum = this_min_rf_weight_ops.GetOps().GetShiftSum();
    auto min_rf_distance = this_min_sum_rf_dist.ComputeWeightBelow(
                               merge.GetResult().GetRoot(), this_min_rf_weight_ops) +
                           shiftsum;
    auto min_rf_count = this_min_sum_rf_dist.MinWeightCount(merge.GetResult().GetRoot(),
                                                            this_min_rf_weight_ops);
    auto max_rf_distance = this_max_sum_rf_dist.ComputeWeightBelow(
                               merge.GetResult().GetRoot(), this_max_rf_weight_ops) +
                           shiftsum;
    auto max_rf_count = this_max_sum_rf_dist.MinWeightCount(merge.GetResult().GetRoot(),
                                                            this_max_rf_weight_ops);

    log_timer.stop();
    std::cout << "Best parsimony score in DAG: " << minparsimony << "\n";
    std::cout << "Worst parsimony score in DAG: " << maxparsimony << "\n";
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Total trees in DAG: " << ntrees
              << "\n";
    std::cout << "Optimal trees in DAG: " << minparsimonytrees << "\n";
    std::cout << "min summed RF distance over trees: " << min_rf_distance << "\n";
    logfile << '\n'
            << iteration << '\t' << ntrees << '\t' << merge.GetResult().GetNodesCount()
            << '\t' << merge.GetResult().GetEdgesCount() << '\t' << minparsimony << '\t'
            << minparsimonytrees << '\t' << maxparsimony << '\t' << min_rf_distance
            << '\t' << max_rf_distance << '\t' << min_rf_count << '\t' << max_rf_count
            << '\t' << log_timer.durationS() << std::flush;
    if (write_intermediate_dag) {
      bool append_changes = (iteration > 0);
      StoreDAG(merge.GetResult(), intermediate_dag_path, output_format, append_changes);
    }
  };
  logger(0);

  bool subtrees = false;
  for (size_t i = 0; i < count; ++i) {
    std::cout << "############ Beginning optimize loop " << std::to_string(i)
              << " #######\n";

    subtrees = (i >= switch_subtrees);
    merge.ComputeResultEdgeMutations();

    SubtreeWeight<BinaryParsimonyScore, MergeDAG> weight{merge.GetResult(),
                                                         main_rng.GenerateSeed()};
    SubtreeWeight<TreeCount, MergeDAG> uniform_sampling_weight{merge.GetResult(),
                                                               main_rng.GenerateSeed()};
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
          std::mt19937 random_generator(main_rng.GenerateSeed());
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

    merge.ComputeResultEdgeMutations();

    SubtreeWeight<SumRFDistance, MergeDAG> min_sum_rf_dist{merge.GetResult(),
                                                           main_rng.GenerateSeed()};
    SumRFDistance min_rf_weight_ops{merge, merge, main_rng.GenerateSeed()};
    SubtreeWeight<MaxSumRFDistance, MergeDAG> max_sum_rf_dist{merge.GetResult(),
                                                              main_rng.GenerateSeed()};
    MaxSumRFDistance max_rf_weight_ops{merge, merge, main_rng.GenerateSeed()};

    auto sample = [&sample_method, &weight, &uniform_sampling_weight, &subtree_node,
                   &min_sum_rf_dist, &min_rf_weight_ops, &max_sum_rf_dist,
                   &max_rf_weight_ops]() {
      switch (sample_method) {
        case SampleMethod::Random:
          return AddMATConversion(weight.SampleTree({}, subtree_node));
        case SampleMethod::UniformRandom:
          return AddMATConversion(
              uniform_sampling_weight.UniformSampleTree({}, subtree_node));
        case SampleMethod::Parsimony:
          return AddMATConversion(weight.MinWeightSampleTree({}, subtree_node));
        case SampleMethod::UniformParsimony:
          return AddMATConversion(weight.MinWeightUniformSampleTree({}, subtree_node));
        case SampleMethod::MinSumRFDistance:
          return AddMATConversion(
              min_sum_rf_dist.MinWeightSampleTree(min_rf_weight_ops, subtree_node));
        case SampleMethod::MaxSumRFDistance:
          return AddMATConversion(
              max_sum_rf_dist.MinWeightSampleTree(max_rf_weight_ops, subtree_node));
        default:
          std::cerr << "ERROR: Invalid SampleMethod" << std::endl;
          std::exit(EXIT_FAILURE);
      }
    }();

    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>> Nodes in sampled (sub)tree: "
              << sample.GetNodesCount() << "\n";
    MAT::Tree mat;
    sample.View().GetRoot().Validate(true);
    sample.View().BuildMAT(mat);
    check_edge_mutations(sample.View().Const());

    if (callback_config == CallbackMethod::AllMoves) {
      Merge_All_Moves_Found_Callback callback{merge, sample.View(),
                                              collapse_empty_fragment_edges};
      optimized_dags.push_back(
          optimize_dag_direct(sample.View(), callback, callback, callback));
    } else if (callback_config == CallbackMethod::BestMovesFixedTree) {
      Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback callback{
          merge,
          sample.View(),
          {move_coeff_nodes, move_coeff_pscore},
          collapse_empty_fragment_edges};
      optimized_dags.push_back(
          optimize_dag_direct(sample.View(), callback, callback, callback));
    } else if (callback_config == CallbackMethod::BestMovesTreeBased) {
      Treebased_Move_Found_Callback callback{merge,
                                             sample.View(),
                                             {move_coeff_nodes, move_coeff_pscore},
                                             collapse_empty_fragment_edges};
      optimized_dags.push_back(
          optimize_dag_direct(sample.View(), callback, callback, callback));
    } else if (callback_config == CallbackMethod::BestMoves) {
      Merge_All_Profitable_Moves_Found_Callback callback{
          merge,
          sample.View(),
          {move_coeff_nodes, move_coeff_pscore},
          collapse_empty_fragment_edges};
      optimized_dags.push_back(
          optimize_dag_direct(sample.View(), callback, callback, callback));
    } else {
      std::cerr << "ERROR: Invalid CallbackMethod" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    auto optimized_view = optimized_dags.back().first.View();
    optimized_view.RecomputeCompactGenomes(false);
    merge.AddDAG(optimized_view);
    logger(i + 1);
  }

  std::cout << "new node coefficient: " << move_coeff_nodes << "\n";
  std::cout << "parsimony score coefficient: " << move_coeff_pscore << "\n";

  logfile.close();

  Benchmark save_timer;
  std::cout << "Saving final DAG..." << std::flush;
  if (final_trim) {
    SubtreeWeight<BinaryParsimonyScore, MergeDAG> parsimonyscorer{
        merge.GetResult(), main_rng.GenerateSeed()};
    merge.ComputeResultEdgeMutations();
    StoreDAG(parsimonyscorer.TrimToMinWeight({}).View(), output_dag_path,
             output_format);
  } else {
    StoreDAG(merge.GetResult(), output_dag_path, output_format);
  }
  save_timer.stop();
  std::cout << "...saved: " << save_timer.durationFormatMs() << std::endl;

  total_timer.stop();
  std::cout << "Total runtime: " << total_timer.durationFormatMs() << std::endl;

  return EXIT_SUCCESS;
}
