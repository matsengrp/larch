#pragma once

#include "larch/spr/scoring_backend.hpp"
#include "larch/contiguous_set.hpp"

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

  // Get Fitch set parts for this node (delegates to backend)
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

/**
 * @brief HypotheticalTree represents a tree after an SPR move has been applied.
 *
 * @tparam DAG The underlying DAG type
 * @tparam Backend The scoring backend type (default: MatOptimizeScoringBackend)
 *
 * The Backend template parameter allows different scoring algorithms to be used:
 * - MatOptimizeScoringBackend: Uses Fitch algorithm for parsimony scoring
 * - MLScoringBackend: Uses ML-based likelihood scoring (stub for now)
 */
template <typename DAG, typename Backend = MatOptimizeScoringBackend<DAG>>
struct HypotheticalTree {
  struct Data {
    // Constructor for matOptimize backend
    template <typename DAGView>
    Data(const DAGView& dag, const Profitable_Moves& move, NodeId new_node,
         bool has_unifurcation_after_move,
         const std::vector<Node_With_Major_Allele_Set_Change>&
             nodes_with_major_allele_set_change);

    // Constructor for ML backend (NodeId-based)
    template <typename DAGView>
    Data(const DAGView& dag, NodeId src, NodeId dst, NodeId lca, NodeId new_node,
         bool has_unifurcation_after_move);

    Backend backend_;
    NodeId new_node_;
    bool has_unifurcation_after_move_;
    ContiguousSet<NodeId> lca_ancestors_;

    // Original move pointers kept for compatibility (matOptimize path only)
    Profitable_Moves move_;
  };
  std::unique_ptr<Data> data_;
};

template <typename DAG, typename Backend, typename CRTP, typename Tag>
struct FeatureConstView<HypotheticalTree<DAG, Backend>, CRTP, Tag> {
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

  // Access to the scoring backend
  const Backend& GetBackend() const;

  // Get score change for this move
  Score GetScoreChange() const;

  const ContiguousSet<NodeId>& GetLCAAncestors() const;

  // Legacy: access to changed fitch set map (matOptimize path)
  // Returns empty map for non-matOptimize backends
  const ContiguousMap<MATNodePtr,
                      ContiguousMap<MutationPosition, Mutation_Count_Change>>&
  GetChangedFitchSetMap() const;
};

template <typename DAG, typename Backend, typename CRTP, typename Tag>
struct FeatureMutableView<HypotheticalTree<DAG, Backend>, CRTP, Tag> {
  std::pair<NodeId, bool> ApplyMove(NodeId lca, NodeId src, NodeId dst) const;

  // Initialize for matOptimize backend
  bool InitHypotheticalTree(const Profitable_Moves& move,
                            const std::vector<Node_With_Major_Allele_Set_Change>&
                                nodes_with_major_allele_set_change) const;

  // Initialize for ML backend (NodeId-based)
  bool InitHypotheticalTree(NodeId src, NodeId dst, NodeId lca) const;
};

// Forward declaration for LongNameOfSPRStorage
template <typename Target, typename Backend>
struct SPRStorage;

template <typename Target, typename Backend>
struct LongNameOfSPRStorage {
  using type_helper = ExtendStorageType<void, Target, Extend::Nodes<HypotheticalNode>,
                                        Extend::DAG<HypotheticalTree<Target, Backend>>>;
  using type = OverlayStorageType<SPRStorage<Target, Backend>, type_helper>;
};

// Specialization of LongNameOf for SPRStorage
template <typename Target, typename Backend>
struct LongNameOf<SPRStorage<Target, Backend>> {
  using type = typename LongNameOfSPRStorage<Target, Backend>::type;
};

/**
 * @brief Storage container for SPR (Subtree Prune and Regraft) operations on
 * phylogenetic trees.
 *
 * SPRStorage extends a target DAG storage with additional features needed for SPR
 * operations, including hypothetical nodes and hypothetical tree structures. It
 * provides the infrastructure for efficiently exploring tree rearrangements by
 * maintaining both the original tree structure and potential modifications that would
 * result from SPR moves.
 *
 * @tparam Target The underlying DAG storage type
 * @tparam Backend The scoring backend type (default: MatOptimizeScoringBackend)
 */
template <typename Target, typename Backend = MatOptimizeScoringBackend<Target>>
struct SPRStorage : LongNameOfSPRStorage<Target, Backend>::type {
  MOVE_ONLY(SPRStorage);

  using LongNameType = typename LongNameOfSPRStorage<Target, Backend>::type;
  using LongNameType::LongNameType;
  using BackendType = Backend;

  static SPRStorage EmptyDefault() {
    static_assert(Target::role == Role::Storage);
    using helper = typename LongNameOfSPRStorage<Target, Backend>::type_helper;
    return SPRStorage{helper::EmptyDefault()};
  }

  static SPRStorage Consume(Target&& target) {
    static_assert(Target::role == Role::Storage);
    using helper = typename LongNameOfSPRStorage<Target, Backend>::type_helper;
    return SPRStorage{helper::Consume(std::move(target))};
  }

  static SPRStorage FromView(const Target& target) {
    static_assert(Target::role == Role::View);
    using helper = typename LongNameOfSPRStorage<Target, Backend>::type_helper;
    return SPRStorage{helper::FromView(target)};
  }
};

// Helper function for matOptimize backend (default)
template <typename Target, typename = std::enable_if_t<Target::role == Role::Storage>>
SPRStorage<Target, MatOptimizeScoringBackend<Target>> AddSPRStorage(Target&& dag) {
  return SPRStorage<Target, MatOptimizeScoringBackend<Target>>::Consume(std::move(dag));
}

template <typename Target, typename = std::enable_if_t<Target::role == Role::View>>
SPRStorage<Target, MatOptimizeScoringBackend<Target>> AddSPRStorage(const Target& dag) {
  return SPRStorage<Target, MatOptimizeScoringBackend<Target>>::FromView(dag);
}

// Helper function with explicit backend
template <typename Backend, typename Target,
          typename = std::enable_if_t<Target::role == Role::Storage>>
SPRStorage<Target, Backend> AddSPRStorageWithBackend(Target&& dag) {
  return SPRStorage<Target, Backend>::Consume(std::move(dag));
}

template <typename Backend, typename Target,
          typename = std::enable_if_t<Target::role == Role::View>>
SPRStorage<Target, Backend> AddSPRStorageWithBackend(const Target& dag) {
  return SPRStorage<Target, Backend>::FromView(dag);
}

#include "larch/impl/spr/spr_view_impl.hpp"
