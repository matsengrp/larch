#pragma once

#include "larch/usher_glue.hpp"
#include "larch/spr/lca.hpp"
#include "larch/contiguous_set.hpp"

struct FitchSet {
  FitchSet(char base) : value_{base} {}
  FitchSet(int base) : value_{static_cast<char>(base)} {}  // TODO
  FitchSet(nuc_one_hot one_hot)
      : value_{[](nuc_one_hot base) {
          static const std::array<char, 4> decode = {'A', 'C', 'G', 'T'};
          return decode.at(one_hot_to_two_bit(base));
        }(one_hot)} {}

  bool find(char base) const {
    switch (base) {
      case 'A':
        return value_ & 1;
      case 'C':
        return value_ & 2;
      case 'G':
        return value_ & 4;
      case 'T':
        return value_ & 8;
    }
    Fail("Unreachable");
  }
  char at(size_t pos) const {
    Assert(pos == 0);
    if (value_ & 1) {
      return 'A';
    } else if (value_ & 2) {
      return 'C';
    } else if (value_ & 4) {
      return 'G';
    } else if (value_ & 8) {
      return 'T';
    }
    Fail("Unreachable");
  }

 private:
  char value_ = 0;
};

struct HypotheticalNode {
  HypotheticalNode Copy() const {
    return {changed_base_sites_.Copy(), has_changed_topology_};
  }

  ContiguousSet<MutationPosition> changed_base_sites_;
  bool has_changed_topology_ = false;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<HypotheticalNode, CRTP, Tag> {
  bool IsMATRoot() const;
  bool IsMoveSource() const;
  bool IsMoveTarget() const;
  bool IsMoveNew() const;

  // returns true if this node is an ancestor of the
  // source-target-LCA (the node returned by HypotheticalTree::GetLCA()),
  // which disqualifies a node from being a nonroot anchor
  // node (this method is used in IsNonrootAnchorNode)
  bool IsLCAAncestor() const;

  bool HasChangedTopology() const;

  auto GetOld() const;

  const ContiguousSet<MutationPosition>& GetChangedBaseSites() const;

  // return a set of site indices at which there are fitch set changes
  [[nodiscard]] ContiguousSet<MutationPosition> GetSitesWithChangedFitchSets() const;

  [[nodiscard]] std::pair<
      MAT::Mutations_Collection,
      std::optional<ContiguousMap<MutationPosition, Mutation_Count_Change>>>
  GetFitchSetParts() const;

  // get the (possibly modified) fitch set at this node at the provided site.
  [[nodiscard]] FitchSet GetFitchSet(MutationPosition site) const;

  // Most of the time this can just return the parent node's
  // changed_base_sites. However, it's different if the node in question is the
  // source node!
  [[nodiscard]] ContiguousSet<MutationPosition> GetParentChangedBaseSites() const;

  [[nodiscard]] CompactGenome ComputeNewCompactGenome() const;

  // A tree fragment has two kinds of anchor nodes: One root anchor node, which
  // is always the parent of HypotheticalTree.GetOldestChangedNode(), and the
  // rest are what I'm calling here "nonroot anchor nodes", whose descendants
  // in the hypothetical tree are guaranteed to also be unchanged from before
  // the SPR move. This method identifies the second kind.
  bool IsNonrootAnchorNode() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<HypotheticalNode, CRTP, Tag> {
  void SetHasChangedTopology() const;
  void PreorderComputeCompactGenome(std::vector<NodeId>& result,
                                    std::vector<EdgeId>& result_edges) const;
};

template <typename DAG>
struct HypotheticalTree {
  struct Data {
    Data(const Profitable_Moves& move, NodeId new_node, bool collapse,
         const std::vector<Node_With_Major_Allele_Set_Change>&
             nodes_with_major_allele_set_change);
    Profitable_Moves move_;
    NodeId new_node_;
    bool collapse_;
    ContiguousMap<MATNodePtr, ContiguousMap<MutationPosition, Mutation_Count_Change>>
        changed_fitch_set_map_;
    ContiguousSet<NodeId> lca_ancestors_;
  };
  std::unique_ptr<Data> data_;  // TODO fixme
};

template <typename DAG, typename CRTP, typename Tag>
struct FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag> {
  // Get the LCA of source and target nodes
  auto GetMoveLCA() const;
  // These return the HypotheticalTreeNodes corresponding to source and target
  // nodes (they're siblings in the hypothetical tree)
  auto GetMoveSource() const;
  auto GetMoveTarget() const;
  auto GetMoveNew() const;
  bool HaveCollapse() const;

  // Returns the HypotheticalTreeNode that used to be the parent of source
  // before the SPR move. TODO: This node may (but need not be) unifurcating
  // after the SPR move, in which case the issue description says it should not
  // appear as a neighbor of any node in the hypothetical tree (it's omitted
  // entirely). However, we will still need to access its data, such as its old
  // compact genome, so this method should return it. the todo is to think
  // through all the consequences of this, and make sure it doesn't cause any
  // problems that we skip it sometimes in the tree traversal.
  auto GetOldSourceParent() const;

  // returns LCA of source and target, or the earliest node with fitch set
  // changes, whichever is higher in the tree.
  [[nodiscard]] auto GetOldestChangedNode() const;

  [[nodiscard]] std::pair<std::vector<NodeId>, std::vector<EdgeId>> GetFragment() const;

  const ContiguousMap<MATNodePtr,
                      ContiguousMap<MutationPosition, Mutation_Count_Change>>&
  GetChangedFitchSetMap() const;

  const ContiguousSet<NodeId>& GetLCAAncestors() const;
};

template <typename DAG, typename CRTP, typename Tag>
struct FeatureMutableView<HypotheticalTree<DAG>, CRTP, Tag> {
  std::pair<NodeId, bool> ApplyMove(NodeId src, NodeId dst) const;
  bool InitHypotheticalTree(const Profitable_Moves& move,
                            const std::vector<Node_With_Major_Allele_Set_Change>&
                                nodes_with_major_allele_set_change);
};

template <typename DAG>
auto SPRStorage(DAG&& dag) {
  auto extend = ExtendStorage(std::forward<DAG>(dag), Extend::Nodes<HypotheticalNode>{},
                              Extend::DAG<HypotheticalTree<std::decay_t<DAG>>>{});
  return OverlayDAGStorage{std::move(extend)};
}

#include "larch/impl/spr/spr_view_impl.hpp"
