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
inline auto& GetFeatureStorage(const FeatureReader<Feature, View>* reader);

#define DAG_FEATURE_FRIENDS                                                \
  template <typename _Feature_, typename _View_>                           \
  friend auto& GetFeatureStorage(const FeatureReader<_Feature_, _View_>*); \
  template <typename, typename>                                            \
  friend class FeatureReader;                                              \
  template <typename, typename>                                            \
  friend class FeatureWriter

#define DAG_VIEW_FRIENDS \
  template <typename>    \
  friend class NodeView; \
  template <typename>    \
  friend class EdgeView; \
  template <typename>    \
  friend class DAGView

template <typename Feature>
class NoGlobalData {
};

template <typename, typename = void>
constexpr bool has_global_data{};

template <typename Feature>
constexpr bool has_global_data<Feature, std::void_t<typename Feature::GlobalData> > =
    true;