#pragma once

#include "larch/usher_glue.hpp"

struct FitchSet {
  FitchSet() = default;
  MOVE_ONLY(FitchSet);

  inline bool operator==(const FitchSet& rhs) const noexcept;
  inline bool operator<(const FitchSet& rhs) const noexcept;
  [[nodiscard]] inline size_t Hash() const noexcept;

  [[nodiscard]] inline FitchSet Copy() const;
};

template <>
struct std::hash<FitchSet> {
  inline std::size_t operator()(const FitchSet& fs) const noexcept;
};

template <>
struct std::equal_to<FitchSet> {
  inline bool operator()(const FitchSet& lhs, const FitchSet& rhs) const noexcept;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<FitchSet, CRTP, Tag> {
  const FitchSet& GetFitchSet() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<FitchSet, CRTP, Tag> {
  auto& operator=(FitchSet&& fitch_set);
};

template <typename DAG>
struct HypotheticalTree {
  struct Data {
    Data(DAG sample_dag, const MAT::Tree& sample_mat, const Profitable_Moves& move,
         const std::vector<Node_With_Major_Allele_Set_Change>&
             nodes_with_major_allele_set_change);
    DAG sample_dag_;
    const MAT::Tree& sample_mat_;
    Profitable_Moves move_;
    std::map<MAT::Node*, std::map<size_t, Mutation_Count_Change>>
        changed_fitch_set_map_;
  };
  std::unique_ptr<Data> data_;
};

template <typename DAG, typename CRTP, typename Tag>
struct FeatureConstView<HypotheticalTree<DAG>, CRTP, Tag> {};

template <typename DAG, typename CRTP, typename Tag>
struct FeatureMutableView<HypotheticalTree<DAG>, CRTP, Tag> {
  void InitHypotheticalTree(DAG sample_dag, const MAT::Tree& sample_mat,
                            const Profitable_Moves& move,
                            const std::vector<Node_With_Major_Allele_Set_Change>&
                                nodes_with_major_allele_set_change);
};

template <typename DAG>
auto SPRStorage(DAG&& dag) {
  auto extend =
      ExtendStorage(std::forward<DAG>(dag), Extend::Nodes<Deduplicate<FitchSet>>{},
                    Extend::DAG<HypotheticalTree<std::decay_t<DAG>>>{});
  return OverlayDAGStorage{std::move(extend)};
}

#include "larch/impl/spr/spr_view_impl.hpp"