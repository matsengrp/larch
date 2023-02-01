#pragma once

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
using SPRStorage =
    OverlayDAGStorage<ExtendDAGStorage<DAG, Extend::Nodes<Deduplicate<FitchSet>>>>;

template <typename DAG>
using SPRDAG = DAGView<const SPRStorage<DAG>>;

template <typename DAG>
using MutableSPRDAG = DAGView<SPRStorage<DAG>>;

#include "larch/impl/spr/spr_view_impl.hpp"