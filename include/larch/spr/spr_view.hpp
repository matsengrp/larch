#pragma once

#include "larch/usher_glue.hpp"

struct FitchSet {
  FitchSet(char base) : value_{base} {}
  FitchSet(int base) : value_{static_cast<char>(base)} {}  // TODO
  FitchSet(nuc_one_hot one_hot)
      : value_{[](nuc_one_hot base) {
          static const std::array<char, 4> decode = {'A', 'C', 'G', 'T'};
          return decode.at(one_hot_to_two_bit(base));
        }(one_hot)} {}

  bool find(char base) const { return value_ == base; }
  char at(size_t) const { return value_; }

 private:
  char value_ = 0;
};

struct HypotheticalNode {
  std::set<MutationPosition> changed_base_sites_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<HypotheticalNode, CRTP, Tag> {
  bool IsMATRoot() const;
  bool IsSource() const;
  bool IsTarget() const;
  bool IsNew() const;
  auto GetOld() const;

  const std::set<MutationPosition>& GetChangedBaseSites() const;

  // return a set of site indices at which there are fitch set changes
  [[nodiscard]] std::set<MutationPosition> GetSitesWithChangedFitchSets() const;

  [[nodiscard]] std::pair<
      MAT::Mutations_Collection,
      std::optional<std::map<MutationPosition, Mutation_Count_Change>>>
  GetFitchSetParts() const;

  // get the (possibly modified) fitch set at this node at the provided site.
  [[nodiscard]] FitchSet GetFitchSet(MutationPosition site) const;

  // Most of the time this can just return the parent node's
  // changed_base_sites. However, it's different if the node in question is the
  // source node!
  [[nodiscard]] std::set<MutationPosition> GetParentChangedBaseSites() const;

  [[nodiscard]] CompactGenome ComputeNewCompactGenome() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<HypotheticalNode, CRTP, Tag> {};

template <typename DAG>
struct HypotheticalTree {
  struct Data {
    Data(const Profitable_Moves& move,
         const std::vector<Node_With_Major_Allele_Set_Change>&
             nodes_with_major_allele_set_change);
    Profitable_Moves move_;
    std::map<const MAT::Node*, std::map<MutationPosition, Mutation_Count_Change>>
        changed_fitch_set_map_;
  };
  std::unique_ptr<Data> data_;  // TODO fixme
};

template <typename DAG, typename CRTP, typename Tag>
struct FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag> {
  // Get the LCA of source and target nodes
  auto GetLCA() const;
  // These return the HypotheticalTreeNodes corresponding to source and target
  // nodes (they're siblings in the hypothetical tree)
  auto GetSource() const;
  auto GetTarget() const;

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
  auto GetOldestChangedNode() const;

  const std::map<const MAT::Node*, std::map<MutationPosition, Mutation_Count_Change>>&
  GetChangedFitchSetMap() const;
};

template <typename DAG, typename CRTP, typename Tag>
struct FeatureMutableView<HypotheticalTree<DAG>, CRTP, Tag> {
  void InitHypotheticalTree(const Profitable_Moves& move,
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