#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

#include <unistd.h>
#include <sys/wait.h>

#include "tools_common.hpp"
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
  const std::string program_desc =
      "larch-usher: tool for exploring tree space of DAG/tree through SPR moves";

  const std::vector<std::string> usage_examples = {
      {"larch-usher -i,--input FILE -o,--output FILE [-c,--count INT]"}};

  const std::vector<std::pair<std::string, std::string>> flag_desc_pairs = {
      {"-i,--input FILE", "Path to input DAG/Tree file (REQUIRED)"},
      {"-o,--output FILE", "Path to output DAG file (REQUIRED)"},
      {"-r,--MAT-refseq-file FILE",
       "Path to json reference sequence file \n"
       "(REQUIRED if input file is a MAT protobuf)"},
      {"-v,--VCF-input-file FILE",
       "Path to VCF file, containing ambiguous leaf sequence data"},
      {"-c,--count INT", "Number of iterations (default: 1)"},
      {"--inter-save INT",
       "Saves a new intermediate DAG file once every given number of iterations \n"
       "(default: no intermediate DAG files saved)"},
      {"-s,--switch-subtrees INT",
       "Switch to optimizing subtrees after the specified number of iterations \n"
       "(default: never)"},
      {"--min-subtree-clade-size INT",
       "The minimum number of leaves in a subtree sampled for optimization \n"
       "(default: 100, ignored without option `-s`)"},
      {"--max-subtree-clade-size INT",
       "The maximum number of leaves in a subtree sampled for optimization \n"
       "(default: 1000, ignored without option `-s`)"},
      {"--move-coeff-nodes INT", "New node coefficient for scoring moves (default: 1)"},
      {"--move-coeff-pscore INT",
       "Parsimony score coefficient for scoring moves (default: 1)"},
      {"--sample-any-tree",
       "Sample any tree for optimization, rather than requiring the sampled tree \n"
       "to maximize parsimony"},
      {"--sample-uniformly",
       "Use a uniform distribution to sample trees for optimization, rather than \n"
       "a natural distribution"},
      {"--sample-method OPTION",
       "Select tree sampling method for optimization (default: max parsimony)\n"
       "[parsimony, random, rf-minsum, rf-maxsum]"},
      {"--callback-option OPTION",
       "Callback configuration choice (default: merge all profitable moves) \n"
       "[best-move, best-move-fixed-tree, best-move-treebased, all-moves]"},
      {"--trim", "Trim optimized DAG after final iteration"},
      {"--keep-fragment-uncollapsed",
       "Keep empty fragment edges, rather than collapsing them"},
      {"--input-format OPTION",
       "Specify input file format (default: inferred) \n"
       "[dagbin, dag-pb, tree-pb, dag-json]"},
      {"--output-format OPTION",
       "Specify output file format (default: inferred) \n"
       "[dagbin, dag-pb]"},
      {"--seed INT", "Set seed for random number generation (default: random)"},
      {"--thread INT", "Set number of cpu threads (default: max allowed by system)"}};

  std::cout << FormatUsage(program_desc, usage_examples, flag_desc_pairs);

  std::exit(EXIT_SUCCESS);
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
  std::optional<uint> write_intermediate_every_x_iters = std::nullopt;
  FileFormat input_format = FileFormat::Infer;
  FileFormat output_format = FileFormat::Infer;
  SampleMethod sample_method = SampleMethod::Parsimony;
  bool sample_uniformly = false;
  size_t iter_count = 1;
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
      ParseOption(name, params, input_dag_path, 1);
    } else if (name == "-o" or name == "--output") {
      ParseOption(name, params, output_dag_path, 1);
    } else if (name == "-c" or name == "--count") {
      ParseOption(name, params, iter_count, 1);
    } else if (name == "-s" or name == "--switch-subtrees") {
      ParseOption(name, params, switch_subtrees, 1);
    } else if (name == "-l" or name == "--logpath") {
      ParseOption(name, params, logfile_path, 1);
    } else if (name == "--move-coeff-pscore") {
      ParseOption(name, params, move_coeff_pscore, 1);
    } else if (name == "--min-subtree-clade-size") {
      ParseOption(name, params, min_subtree_clade_size, 1);
    } else if (name == "--max-subtree-clade-size") {
      ParseOption(name, params, max_subtree_clade_size, 1);
    } else if (name == "--move-coeff-nodes") {
      ParseOption(name, params, move_coeff_nodes, 1);
    } else if (name == "--sample-method") {
      std::string temp;
      ParseOption(name, params, temp, 1);
      if (temp == "random") {
        sample_method = SampleMethod::Random;
      } else if (temp == "parsimony") {
        sample_method = SampleMethod::Parsimony;
      } else if (temp == "rf-minsum") {
        sample_method = SampleMethod::MinSumRFDistance;
      } else if (temp == "rf-maxsum") {
        sample_method = SampleMethod::MaxSumRFDistance;
      } else {
        std::cerr << "ERROR: Unknown --sample-method argument '" << temp << "'"
                  << std::endl;
        Fail();
      }
    } else if (name == "--sample-any-tree") {
      ParseOption<false>(name, params, sample_method, 0);
      sample_method = SampleMethod::Random;
    } else if (name == "--sample-uniformly") {
      ParseOption<false>(name, params, sample_uniformly, 0);
      sample_uniformly = true;
    } else if (name == "-r" or name == "--MAT-refseq-file") {
      ParseOption(name, params, refseq_path, 1);
    } else if (name == "-v" or name == "--VCF-input-file") {
      ParseOption(name, params, vcf_path, 1);
    } else if (name == "--callback-option") {
      std::string temp;
      ParseOption(name, params, temp, 1);
      if (temp == "best-moves") {
        callback_config = CallbackMethod::BestMoves;
      } else if (temp == "best-moves-fixed-tree") {
        callback_config = CallbackMethod::BestMovesFixedTree;
      } else if (temp == "best-moves-treebased") {
        callback_config = CallbackMethod::BestMovesTreeBased;
      } else if (temp == "all-moves") {
        callback_config = CallbackMethod::AllMoves;
      } else {
        std::cerr << "ERROR: Unknown --callback-option argument '" << temp << "'"
                  << std::endl;
        Fail();
      }
    } else if (name == "--keep-fragment-uncollapsed") {
      ParseOption<false>(name, params, collapse_empty_fragment_edges, 0);
      collapse_empty_fragment_edges = false;
    } else if (name == "--quiet") {
      ParseOption<false>(name, params, write_intermediate_dag, 0);
      write_intermediate_dag = false;
    } else if (name == "--inter-save") {
      uint temp;
      ParseOption(name, params, temp, 1);
      write_intermediate_every_x_iters = temp;
    } else if (name == "--trim") {
      ParseOption<false>(name, params, final_trim, 0);
      final_trim = true;
    } else if (name == "--input-format") {
      std::string temp;
      ParseOption(name, params, temp, 1);
      input_format = InferFileFormat(temp);
    } else if (name == "--output-format") {
      std::string temp;
      ParseOption(name, params, temp, 1);
      output_format = InferFileFormat(temp);
    } else if (name == "--seed") {
      uint32_t temp;
      ParseOption(name, params, temp, 1);
      user_seed = temp;
    } else if (name == "--thread") {
      ParseOption(name, params, thread_count, 1);
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
  MADAGStorage<> input_dag_storage = LoadDAG(input_dag_path, input_format, refseq_path);
  auto input_dag = input_dag_storage.View();
  load_timer.stop();
  std::cout << "...loaded: " << load_timer.durationFormatMs() << std::endl;

  input_dag.RecomputeCompactGenomes(true);
  LoadVCFData(input_dag_storage, vcf_path);
  // if the DAG is from a DAG protobuf file, then it needs to be equipped with
  // SampleIds
  if (vcf_path.empty()) {
    input_dag.SampleIdsFromCG();
  }
  Merge merge_dag{input_dag.GetReferenceSequence()};
  merge_dag.AddDAG(input_dag);
  std::vector<
      std::pair<decltype(AddMATConversion(MADAGStorage<>::EmptyDefault())), MAT::Tree>>
      optimized_dags;
  merge_dag.ComputeResultEdgeMutations();

  Benchmark log_timer;
  auto logger = [&input_dag, &merge_dag, &logfile, &log_timer, &intermediate_dag_path,
                 &write_intermediate_dag, &write_intermediate_every_x_iters,
                 &output_format, &main_rng](size_t iteration) {
    std::cout << "############ Logging for iteration " << iteration << " #######\n";
    merge_dag.ComputeResultEdgeMutations();

    // Tree count
    SubtreeWeight<TreeCount, MergeDAG> tree_counter{merge_dag.GetResult(),
                                                    main_rng.GenerateSeed()};
    auto tree_count =
        tree_counter.ComputeWeightBelow(merge_dag.GetResult().GetRoot(), {});
    // Min Parsimony score
    SubtreeWeight<BinaryParsimonyScore, MergeDAG> min_parsimony_scorer{
        merge_dag.GetResult()};
    auto min_parsimony_score =
        min_parsimony_scorer.ComputeWeightBelow(merge_dag.GetResult().GetRoot(), {});
    auto min_parsimony_count =
        min_parsimony_scorer.MinWeightCount(merge_dag.GetResult().GetRoot(), {});
    // Max Parsimony score
    SubtreeWeight<MaxBinaryParsimonyScore, MergeDAG> max_parsimony_scorer{
        merge_dag.GetResult()};
    auto max_parsimony_score = max_parsimony_scorer.ComputeWeightBelow(
        merge_dag.GetResult().GetRoot(), MaxBinaryParsimonyScore{});
    // Min Sum RF Distance
    SumRFDistance min_sum_rf_dist_weight_ops{merge_dag, merge_dag};
    auto min_shift_sum = min_sum_rf_dist_weight_ops.GetOps().GetShiftSum();
    SubtreeWeight<SumRFDistance, MergeDAG> min_sum_rf_dist_scorer{
        merge_dag.GetResult()};
    auto min_sum_rf_distance = min_sum_rf_dist_scorer.ComputeWeightBelow(
        merge_dag.GetResult().GetRoot(), min_sum_rf_dist_weight_ops);
    min_sum_rf_distance += min_shift_sum;
    auto min_sum_rf_count = min_sum_rf_dist_scorer.MinWeightCount(
        merge_dag.GetResult().GetRoot(), min_sum_rf_dist_weight_ops);
    // Max Sum RF Distance
    MaxSumRFDistance max_sum_rf_dist_weight_ops{merge_dag, merge_dag};
    auto max_shift_sum = max_sum_rf_dist_weight_ops.GetOps().GetShiftSum();
    SubtreeWeight<MaxSumRFDistance, MergeDAG> max_sum_rf_dist_scorer{
        merge_dag.GetResult()};
    auto max_sum_rf_distance = max_sum_rf_dist_scorer.ComputeWeightBelow(
        merge_dag.GetResult().GetRoot(), max_sum_rf_dist_weight_ops);
    max_sum_rf_distance += max_shift_sum;
    auto max_sum_rf_count = max_sum_rf_dist_scorer.MinWeightCount(
        merge_dag.GetResult().GetRoot(), max_sum_rf_dist_weight_ops);

    log_timer.stop();

    std::cout << "Min parsimony score in DAG: " << min_parsimony_score << "\n";
    std::cout << "Max parsimony score in DAG: " << max_parsimony_score << "\n";
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Total trees in DAG: " << tree_count
              << "\n";
    std::cout << "Optimal trees in DAG: " << min_parsimony_count << "\n";
    std::cout << "Min summed RF distance over trees: " << min_sum_rf_distance << "\n";
    std::cout << "Max summed RF distance over trees: " << max_sum_rf_distance << "\n";

    logfile << '\n'
            << iteration << '\t' << tree_count << '\t'
            << merge_dag.GetResult().GetNodesCount() << '\t'
            << merge_dag.GetResult().GetEdgesCount() << '\t' << min_parsimony_score
            << '\t' << min_parsimony_count << '\t' << max_parsimony_score << '\t'
            << min_sum_rf_distance << '\t' << max_sum_rf_distance << '\t'
            << min_sum_rf_count << '\t' << max_sum_rf_count << '\t'
            << log_timer.durationS() << std::flush;

    if (write_intermediate_dag) {
      bool append_changes = (iteration > 0);
      StoreDAG(merge_dag.GetResult(), intermediate_dag_path, output_format,
               append_changes);
      if (write_intermediate_every_x_iters.has_value() and
          (iteration % write_intermediate_every_x_iters.value() == 0)) {
        std::string intermediate_dag_path_final =
            intermediate_dag_path + "." + std::to_string(iteration);
        std::cout << "############ Saving intermediate DAG file to: "
                  << intermediate_dag_path_final << std::endl;
        std::filesystem::copy_file(intermediate_dag_path, intermediate_dag_path_final,
                                   std::filesystem::copy_options::overwrite_existing);
      }
    }
  };
  logger(0);

  bool subtrees = false;
  for (size_t i = 0; i < iter_count; ++i) {
    std::cout << "############ Beginning optimize loop " << i << " #######\n";

    subtrees = (i >= switch_subtrees);
    merge_dag.ComputeResultEdgeMutations();

    SubtreeWeight<BinaryParsimonyScore, MergeDAG> weight{merge_dag.GetResult(),
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

            size_t clade_size = merge_dag.GetResultNodeLabels()
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

    merge_dag.ComputeResultEdgeMutations();

    auto sample_tree = [&merge_dag, &sample_method, &weight, &subtree_node,
                        &main_rng]() {
      SubtreeWeight<TreeCount, MergeDAG> uniform_sampling_weight{
          merge_dag.GetResult(), main_rng.GenerateSeed()};
      SubtreeWeight<SumRFDistance, MergeDAG> min_sum_rf_dist{merge_dag.GetResult(),
                                                             main_rng.GenerateSeed()};
      SumRFDistance min_rf_weight_ops{merge_dag, merge_dag, main_rng.GenerateSeed()};
      SubtreeWeight<MaxSumRFDistance, MergeDAG> max_sum_rf_dist{
          merge_dag.GetResult(), main_rng.GenerateSeed()};
      MaxSumRFDistance max_rf_weight_ops{merge_dag, merge_dag, main_rng.GenerateSeed()};

      if (sample_method == SampleMethod::Random) {
        return AddMATConversion(weight.SampleTree({}, subtree_node));
      } else if (sample_method == SampleMethod::UniformRandom) {
        return AddMATConversion(
            uniform_sampling_weight.UniformSampleTree({}, subtree_node));
      } else if (sample_method == SampleMethod::Parsimony) {
        return AddMATConversion(weight.MinWeightSampleTree({}, subtree_node));
      } else if (sample_method == SampleMethod::UniformParsimony) {
        return AddMATConversion(weight.MinWeightUniformSampleTree({}, subtree_node));
      } else if (sample_method == SampleMethod::MinSumRFDistance) {
        return AddMATConversion(
            min_sum_rf_dist.MinWeightSampleTree(min_rf_weight_ops, subtree_node));
      } else if (sample_method == SampleMethod::MaxSumRFDistance) {
        return AddMATConversion(
            max_sum_rf_dist.MinWeightSampleTree(max_rf_weight_ops, subtree_node));
      } else {
        std::cerr << "ERROR: Invalid SampleMethod" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }();

    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>> Nodes in sampled (sub)tree: "
              << sample_tree.GetNodesCount() << "\n";
    MAT::Tree mat;
    sample_tree.View().GetRoot().Validate(true);
    sample_tree.View().BuildMAT(mat);
    check_edge_mutations(sample_tree.View().Const());

    if (callback_config == CallbackMethod::AllMoves) {
      Merge_All_Moves_Found_Callback callback{merge_dag, sample_tree.View(),
                                              collapse_empty_fragment_edges};
      optimized_dags.push_back(
          optimize_dag_direct(sample_tree.View(), callback, callback, callback));
    } else if (callback_config == CallbackMethod::BestMovesFixedTree) {
      Merge_All_Profitable_Moves_Found_Fixed_Tree_Callback callback{
          merge_dag,
          sample_tree.View(),
          {move_coeff_nodes, move_coeff_pscore},
          collapse_empty_fragment_edges};
      optimized_dags.push_back(
          optimize_dag_direct(sample_tree.View(), callback, callback, callback));
    } else if (callback_config == CallbackMethod::BestMovesTreeBased) {
      Treebased_Move_Found_Callback callback{merge_dag,
                                             sample_tree.View(),
                                             {move_coeff_nodes, move_coeff_pscore},
                                             collapse_empty_fragment_edges};
      optimized_dags.push_back(
          optimize_dag_direct(sample_tree.View(), callback, callback, callback));
    } else if (callback_config == CallbackMethod::BestMoves) {
      Merge_All_Profitable_Moves_Found_Callback callback{
          merge_dag,
          sample_tree.View(),
          {move_coeff_nodes, move_coeff_pscore},
          collapse_empty_fragment_edges};
      optimized_dags.push_back(
          optimize_dag_direct(sample_tree.View(), callback, callback, callback));
    } else {
      std::cerr << "ERROR: Invalid CallbackMethod" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    auto optimized_view = optimized_dags.back().first.View();
    optimized_view.RecomputeCompactGenomes(false);
    merge_dag.AddDAG(optimized_view);
    logger(i + 1);
  }

  std::cout << "new node coefficient: " << move_coeff_nodes << "\n";
  std::cout << "parsimony score coefficient: " << move_coeff_pscore << "\n";

  logfile.close();

  Benchmark save_timer;
  std::cout << "Saving final DAG..." << std::flush;
  if (final_trim) {
    SubtreeWeight<BinaryParsimonyScore, MergeDAG> parsimonyscorer{
        merge_dag.GetResult(), main_rng.GenerateSeed()};
    merge_dag.ComputeResultEdgeMutations();
    StoreDAG(parsimonyscorer.TrimToMinWeight({}).View(), output_dag_path,
             output_format);
  } else {
    StoreDAG(merge_dag.GetResult(), output_dag_path, output_format);
  }
  save_timer.stop();
  std::cout << "...saved: " << save_timer.durationFormatMs() << std::endl;

  total_timer.stop();
  std::cout << "Total runtime: " << total_timer.durationFormatMs() << std::endl;

  return EXIT_SUCCESS;
}
