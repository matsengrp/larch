#pragma once

#include "larch/usher_glue.hpp"
#include "larch/contiguous_set.hpp"

struct FitchSet {
  FitchSet(char base) : value_{base} { Assert(value_ > 0); }
  FitchSet(int base) : value_{static_cast<char>(base)} { Assert(value_ > 0); }
  bool find(char base) const {
    switch (base) {
      case 'A':
        return (value_ & 1) != 0;
      case 'C':
        return (value_ & 2) != 0;
      case 'G':
        return (value_ & 4) != 0;
      case 'T':
        return (value_ & 8) != 0;  // NOLINT
    }
    Fail("Unreachable");
  }
  char at([[maybe_unused]] size_t pos) const {
    Assert(pos == 0);
    if ((value_ & 1) != 0) {
      return 'A';
    } else if ((value_ & 2) != 0) {
      return 'C';
    } else if ((value_ & 4) != 0) {
      return 'G';
    } else if ((value_ & 8) != 0) {  // NOLINT
      return 'T';
    }
    Fail("Unreachable");
  }

 private:
  char value_ = 0;
};

inline nuc_one_hot base_to_singleton(MutationBase base) {
  switch (base.ToChar()) {
    case 'A':
      return 1;
    case 'C':
      return 2;
    case 'G':
      return 4;
    case 'T':
      return 8;  // NOLINT
  }
  Fail("unrecognized base");
}

struct HypotheticalNode {
  template <typename CRTP>
  HypotheticalNode Copy(const CRTP*) const {
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
  [[nodiscard]] FitchSet GetFitchSetAtSite(MutationPosition site) const;

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
  void PreorderTraverseCollectFragmentEdges(std::vector<EdgeId>& fragment_edges) const;
};

template <typename DAG>
struct HypotheticalTree {
  struct Data {
    Data(Profitable_Moves move, NodeId new_node, bool has_unifurcation_after_move,
         const std::vector<Node_With_Major_Allele_Set_Change>&
             nodes_with_major_allele_set_change);
    Profitable_Moves move_;
    NodeId new_node_;
    bool has_unifurcation_after_move_;
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

  // TODO_DR: Used by larch_usher.cpp
  auto GetMoveSources() const;
  auto GetMoveTargets() const;

  auto GetMoveNew() const;
  bool HasUnifurcationAfterMove() const;

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

  [[nodiscard]] auto MakeFragment() const;
  [[nodiscard]] auto MakeUncollapsedFragment() const;

  [[nodiscard]] std::pair<std::vector<NodeId>, std::vector<EdgeId>>
  CollapseEmptyFragmentEdges(const std::vector<NodeId>& fragment_nodes,
                             const std::vector<EdgeId>& fragment_edges) const;

  const ContiguousMap<MATNodePtr,
                      ContiguousMap<MutationPosition, Mutation_Count_Change>>&
  GetChangedFitchSetMap() const;

  const ContiguousSet<NodeId>& GetLCAAncestors() const;
};

template <typename DAG, typename CRTP, typename Tag>
struct FeatureMutableView<HypotheticalTree<DAG>, CRTP, Tag> {
  std::pair<NodeId, bool> ApplyMove(NodeId lca, NodeId src, NodeId dst) const;
  // TODO_DR: Remove this!
  // std::pair<NodeId, bool> ApplyMove(NodeId lca, std::vector<NodeId> src,
  //                                   std::vector<NodeId> dst) const;
  bool InitHypotheticalTree(const Profitable_Moves& move,
                            const std::vector<Node_With_Major_Allele_Set_Change>&
                                nodes_with_major_allele_set_change) const;
};

template <typename Target>
/**
 * @brief Storage container for SPR (Subtree Prune and Regraft) operations on phylogenetic trees.
 * 
 * SPRStorage extends a target DAG storage with additional features needed for SPR operations,
 * including hypothetical nodes and hypothetical tree structures. It provides the infrastructure
 * for efficiently exploring tree rearrangements by maintaining both the original tree structure
 * and potential modifications that would result from SPR moves.
 */
struct SPRStorage;

template <typename Target>
struct LongNameOf<SPRStorage<Target>> {
  using type_helper = ExtendStorageType<void, Target, Extend::Nodes<HypotheticalNode>,
                                        Extend::DAG<HypotheticalTree<Target>>>;
  using type = OverlayStorageType<SPRStorage<Target>, type_helper>;
};

template <typename Target>
struct SPRStorage : LongNameOf<SPRStorage<Target>>::type {
  MOVE_ONLY(SPRStorage);

  using LongNameType = typename LongNameOf<SPRStorage>::type;
  using LongNameType::LongNameType;

  static SPRStorage EmptyDefault() {
    static_assert(Target::role == Role::Storage);
    using helper = typename LongNameOf<SPRStorage<Target>>::type_helper;
    return SPRStorage{helper::EmptyDefault()};
  }

  static SPRStorage Consume(Target&& target) {
    static_assert(Target::role == Role::Storage);
    using helper = typename LongNameOf<SPRStorage<Target>>::type_helper;
    return SPRStorage{helper::Consume(std::move(target))};
  }

  static SPRStorage FromView(const Target& target) {
    static_assert(Target::role == Role::View);
    using helper = typename LongNameOf<SPRStorage<Target>>::type_helper;
    return SPRStorage{helper::FromView(target)};
  }
};

template <typename Target, typename = std::enable_if_t<Target::role == Role::Storage>>
SPRStorage<Target> AddSPRStorage(Target&& dag) {
  return SPRStorage<Target>::Consume(std::move(dag));
}

template <typename Target, typename = std::enable_if_t<Target::role == Role::View>>
SPRStorage<Target> AddSPRStorage(const Target& dag) {
  return SPRStorage<Target>::FromView(dag);
}

#include "larch/impl/spr/spr_view_impl.hpp"
