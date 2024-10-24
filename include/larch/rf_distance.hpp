/*
computing the summed RF distance over all trees in the DAG.
This routine requires 4 traversals of the nodes in DAG to calculate the rf distance
sums:
- one postorder traversal to assign subtree counts below each node, and
- one in a preorder traversal that uses the subtree counts and computes above-tree
counts
- a loop over the nodes (in any order) accumulates the two types of counts and saves
them in an accumulation based on unique leafsets
- finally, a postorder traversal to calculate the sums of RF distances between trees in
a dag and all the trees in a reference, using the methods provided by a
SubtreeWeight<SumRFDistance, DAGTYPE> object.
*/
#pragma once

#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/tree_count.hpp"
#include "larch/subtree/simple_weight_ops.hpp"
#include "larch/merge/merge.hpp"

template <typename SubtreeWeight, typename Node>
ArbitraryInt ComputeAboveTreeCount(Node node,
                                   std::vector<ArbitraryInt>& above_tree_counts,
                                   SubtreeWeight& cached_below_tree_counts) {
  /* This routine takes an argument that is a cached subtreeWeight
   * `cached_below_tree_counts'. If this subtreeWeight is the TreeCount weight (which
   * returns the number of subtrees below each node), then this routine computes the
   * number of trees "above" a given node. A tree is "above" the node if it contains the
   * node, and taking the graph union with a subtree below that node yields a tree on
   * the full leaf set that belongs in the DAG.
   */

  if (above_tree_counts.at(node.GetId().value) > 0) {
    return above_tree_counts.at(node.GetId().value);
  } else if (node.IsUA()) {
    above_tree_counts.at(node.GetId().value) = 1;
    return 1;
  }

  ArbitraryInt above = 0;
  for (auto parent_edge : node.GetParents()) {
    auto parent_node = parent_edge.GetParent();
    auto current_clade = parent_node.GetClade(parent_edge.GetClade());
    ArbitraryInt below_parent = 1;

    for (auto clade : parent_node.GetClades()) {
      if (not ranges::equal(clade, current_clade)) {
        // sum up the list of values {cached_below_tree_counts[edge] : edge in clade}
        auto below_clade = ranges::accumulate(
            clade | ranges::views::transform([&cached_below_tree_counts](auto edge) {
              return cached_below_tree_counts.ComputeWeightBelow(edge.GetChild(), {});
            }),
            ArbitraryInt{0});
        // multiply these sums across all the (other) clades below the current parent
        // node
        below_parent *= below_clade;
      }
    }
    // sum all values over parents for current node
    above += ComputeAboveTreeCount(parent_node, above_tree_counts,
                                   cached_below_tree_counts) *
             below_parent;
  }
  above_tree_counts.at(node.GetId().value) = above;
  return above;
}

struct LeafSetKey {
  size_t operator()(const LeafSet* key) const {
    size_t hash = 0;
    for (const auto& clade : key->GetClades()) {
      for (const auto* leaf : clade) {
        hash = HashCombine(hash, leaf->Hash());
      }
    }
    return hash;
  }

  bool operator()(const LeafSet* lhs, const LeafSet* rhs) const {
    auto& lhs_clades = lhs->GetClades();
    auto& rhs_clades = rhs->GetClades();
    if (lhs_clades.size() != rhs_clades.size()) {
      return false;
    }
    for (size_t i = 0; i < lhs_clades.size(); ++i) {
      auto& lhs_clade = lhs_clades[i];
      auto& rhs_clade = rhs_clades[i];
      if (lhs_clade.size() != rhs_clade.size()) {
        return false;
      }
      for (size_t j = 0; j < lhs_clade.size(); ++j) {
        if (lhs_clade[j]->sample_id_ != rhs_clade[j]->sample_id_) {
          return false;
        }
      }
    }
    return true;
  }
};

// Create a BinaryOperatorWeightOps for computing sum RF distances to the provided
// reference DAG (using SimpleWeightOps):
struct SumRFDistance_ {
  using Weight = ArbitraryInt;
  static inline Weight Identity = 0;
  ArbitraryInt num_trees_in_dag;
  std::map<std::set<std::string>, ArbitraryInt> leafset_to_full_treecount;
  ArbitraryInt shift_sum_;

 private:
  const Merge* reference_dag_;
  const Merge* compute_dag_;

 public:
  explicit SumRFDistance_(const Merge& reference_dag, const Merge& compute_dag)
      : reference_dag_{&reference_dag}, compute_dag_{&compute_dag} {
    SubtreeWeight<TreeCount, MergeDAG> below_tree_counts{reference_dag.GetResult()};
    std::vector<ArbitraryInt> above_tree_counts;
    above_tree_counts.resize(GetReferenceDAG().GetResult().GetNodesCount());
    num_trees_in_dag = below_tree_counts.ComputeWeightBelow(
        GetReferenceDAG().GetResult().GetRoot(), {});

    // create a list of unique (topologically) nodes in the DAG, and accumulate
    // above_tree_counts[n]*below_tree_counts[n] by adding over all n with identical
    // clade sets
    for (auto node : GetReferenceDAG().GetResult().GetNodes()) {
      ComputeAboveTreeCount(node, above_tree_counts, below_tree_counts);
      if (not node.IsUA()) {
        auto clade = [this, node] {
          auto& label = GetReferenceDAG().GetResultNodeLabels().at(node);
          std::set<std::string> leafs;
          for (auto l : label.GetLeafSet()->ToParentClade(label.GetSampleId())) {
            leafs.insert(l->ToString());
          }
          return leafs;
        }();
        leafset_to_full_treecount[clade] +=
            above_tree_counts[node.GetId().value] *
            below_tree_counts.ComputeWeightBelow(node, {});
      }
    }
    // sum all of the values in leafset_to_full_treecount
    shift_sum_ = ranges::accumulate(leafset_to_full_treecount | ranges::views::values,
                                    ArbitraryInt{0});
  }

  template <typename DAG>
  static Weight ComputeLeaf(DAG, NodeId) {
    return 0;
  }

  template <typename DAG>
  Weight ComputeEdge(DAG dag, EdgeId edge_id) const {
    // clade should be a LeafSet, the key type in the mapping
    // leafset_to_full_treecount:
    auto clade = [this, dag, edge_id] {
      auto edge = dag.Get(edge_id);
      auto& label = GetComputeDAG().GetResultNodeLabels().at(edge.GetChild());
      std::set<std::string> leafs;
      for (auto l : label.GetLeafSet()->ToParentClade(label.GetSampleId())) {
        leafs.insert(l->ToString());
      }
      return leafs;
    }();
    auto record = leafset_to_full_treecount.find(clade);
    if (record == leafset_to_full_treecount.end()) {
      return num_trees_in_dag;
    } else {
      return num_trees_in_dag - (2 * record->second);
    }
  }

  bool Compare(Weight lhs, Weight rhs) const { return lhs < rhs; }

  bool CompareEqual(Weight lhs, Weight rhs) const { return lhs == rhs; }

  Weight Combine(Weight lhs, Weight rhs) const { return lhs + rhs; }

  const ArbitraryInt& GetShiftSum() const { return shift_sum_; }

  const Merge& GetReferenceDAG() const { return *reference_dag_; }

  const Merge& GetComputeDAG() const { return *compute_dag_; }
};

// Not sure if you can substruct a templated struct like this...?
struct SumRFDistance : SimpleWeightOps<SumRFDistance_> {
  explicit SumRFDistance(const Merge& reference_dag, const Merge& compute_dag)
      : SimpleWeightOps<SumRFDistance_>{SumRFDistance_{reference_dag, compute_dag}} {}
};

// Create a WeightOps for computing RF distances to the provided reference tree:
struct RFDistance : SumRFDistance {
  explicit RFDistance(const Merge& reference_dag, const Merge& compute_dag)
      : SumRFDistance{reference_dag, compute_dag} {
    Assert(reference_dag.GetResult().IsTree());
    // now behave exactly like SumRFDistance
  }
};

struct MaxSumRFDistance_ : SumRFDistance_ {
  using Weight = typename SumRFDistance_::Weight;
  explicit MaxSumRFDistance_(const Merge& reference_dag, const Merge& compute_dag)
      : SumRFDistance_{reference_dag, compute_dag} {}
  bool Compare(Weight lhs, Weight rhs) const { return lhs > rhs; }
};

struct MaxSumRFDistance : SimpleWeightOps<MaxSumRFDistance_> {
  explicit MaxSumRFDistance(const Merge& reference_dag, const Merge& compute_dag)
      : SimpleWeightOps<MaxSumRFDistance_>{
            MaxSumRFDistance_{reference_dag, compute_dag}} {}
};

// Create a WeightOps for computing RF distances to the provided reference tree:
struct MaxRFDistance : MaxSumRFDistance {
  explicit MaxRFDistance(const Merge& reference_dag, const Merge& compute_dag)
      : MaxSumRFDistance{reference_dag, compute_dag} {
    Assert(reference_dag.GetResult().IsTree());
    // now behave exactly like MaxSumRFDistance
  }
};
