#pragma once

#include "larch/usher_glue.hpp"

struct FitchSet {
  template <typename X>
  FitchSet(X&&);
};

struct HypotheticalNode {
  std::set<MutationPosition> changed_base_sites_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<HypotheticalNode, CRTP, Tag> {
  const MAT::Node& GetMATNode() const;
  bool IsSource() const;
  bool IsTarget() const;
  bool IsNew() const;

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
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<HypotheticalNode, CRTP, Tag> {};

template <typename DAG>
struct HypotheticalTree {
  struct Data {
    Data(DAG sample_dag, const MAT::Tree& sample_mat, const Profitable_Moves& move,
         const std::vector<Node_With_Major_Allele_Set_Change>&
             nodes_with_major_allele_set_change);
    DAG sample_dag_;
    const MAT::Tree& sample_mat_;
    Profitable_Moves move_;
    std::map<const MAT::Node*, std::map<MutationPosition, Mutation_Count_Change>>
        changed_fitch_set_map_;
  };
  std::unique_ptr<Data> data_;
};

template <typename DAG, typename CRTP, typename Tag>
struct FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag> {
  const MAT::Tree& GetMAT() const;
  auto GetSource() const;
  auto GetTarget() const;
  const std::map<const MAT::Node*, std::map<MutationPosition, Mutation_Count_Change>>&
  GetChangedFitchSetMap() const;
};

template <typename DAG, typename CRTP, typename Tag>
struct FeatureMutableView<HypotheticalTree<DAG>, CRTP, Tag> {
  void InitHypotheticalTree(DAG sample_dag, const MAT::Tree& sample_mat,
                            const Profitable_Moves& move,
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