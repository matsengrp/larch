#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

#include "larch/common.hpp"

struct NodeId {
  size_t value = NoId;
};

inline bool operator==(NodeId lhs, NodeId rhs);
inline bool operator<(NodeId lhs, NodeId rhs);

template <>
struct std::hash<NodeId> {
  inline size_t operator()(NodeId id) const noexcept;
};

struct EdgeId {
  size_t value = NoId;
};

inline bool operator==(EdgeId lhs, EdgeId rhs);
inline bool operator<(EdgeId lhs, EdgeId rhs);

struct CladeIdx {
  size_t value = NoId;
};

inline bool operator==(CladeIdx lhs, CladeIdx rhs);
inline bool operator<(CladeIdx lhs, CladeIdx rhs);

template <typename Feature, typename View>
class FeatureReader {};
template <typename Feature, typename View>
class FeatureWriter {};

template <typename Feature, typename View>
struct FeatureTraits {
  using Reader = FeatureReader<Feature, View>;
  using Writer = FeatureWriter<Feature, View>;
  template <typename Id>
  using GlobalData = void;
};

#define DAG_FEATURE_FRIENDS     \
  template <typename, typename> \
  friend class FeatureReader;   \
  template <typename, typename> \
  friend class FeatureWriter

#define DAG_VIEW_FRIENDS \
  template <typename>    \
  friend class NodeView; \
  template <typename>    \
  friend class EdgeView; \
  template <typename>    \
  friend class DAGView

template <typename Feature, typename View>
inline const auto& GetFeature(FeatureReader<Feature, View>* reader);

template <typename Feature, typename View>
inline void SetFeature(FeatureWriter<Feature, View>* writer, Feature&& feature);

template <typename Feature, typename Tuple, size_t I = 0>
constexpr auto GlobalDataFor() {
  if constexpr (I >= std::tuple_size_v<Tuple>) {
    return nullptr;
  } else {
    using Type = typename std::tuple_element_t<I, Tuple>::GlobalData::FeatureType;
    if constexpr (std::is_same_v<Type, Feature>) {
      return typename std::tuple_element_t<I, Tuple>::GlobalData{};
    } else {
      return GlobalDataFor<Feature, Tuple, I + 1>();
    }
  }
}
