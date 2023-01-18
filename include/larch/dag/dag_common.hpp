#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * Used by specialization on the Feature parameter to add functions to
 * an attachable feature. Functions declared in FeatureConstView are
 * for read-only access to the Feature's storage, and in FeatureMutableView
 * are for modifying access.
 * @{
 */
template <typename Feature, typename CRTP, typename Tag = Feature>
struct FeatureConstView;
template <typename Feature, typename CRTP, typename Tag = Feature>
struct FeatureMutableView;
/** @} */

template <typename Feature, typename CRTP>
struct ExtraFeatureConstView {};
template <typename Feature, typename CRTP>
struct ExtraFeatureMutableView {};

/**
 * Used by specialization on the Feature parameter to add extra storage for
 * features that are attached to node or edge (not the DAG itself). Extra storage
 * is not per-element, but global for the entire DAG.
 */
template <typename Feature>
struct ExtraFeatureStorage {};

template <typename Id, typename DAGViewType>
struct ElementView;

template <typename DAGStorageType, typename DAGViewType>
struct DefaultViewBase {
  using DAGViewBase = std::conditional_t<
      std::is_const_v<DAGStorageType>,
      typename DAGStorageType::template ConstDAGViewBase<DAGViewType>,
      typename DAGStorageType::template MutableDAGViewBase<DAGViewType>>;

  template <typename Id>
  using ElementViewBase =
      std::conditional_t<std::is_const_v<DAGStorageType>,
                         typename DAGStorageType::template ConstElementViewBase<
                             Id, ElementView<Id, DAGViewType>>,
                         typename DAGStorageType::template MutableElementViewBase<
                             Id, ElementView<Id, DAGViewType>>>;
};

/**
 * Called by functions in FeatureConstView and FeatureMutableView specializations
 * to access the actual storage of the feature. Accepting `this` as sole parameter.
 * @{
 */
template <typename CRTP, typename Feature, typename Tag>
auto& GetFeatureStorage(const FeatureMutableView<Feature, CRTP, Tag>* feature);

template <typename CRTP, typename Feature, typename Tag>
const auto& GetFeatureStorage(const FeatureConstView<Feature, CRTP, Tag>* feature);
/** @} */

struct NodeId {
  size_t value = NoId;
};

inline bool operator==(NodeId lhs, NodeId rhs);
inline bool operator!=(NodeId lhs, NodeId rhs);
inline bool operator<(NodeId lhs, NodeId rhs);

template <>
struct std::hash<NodeId> {
  inline size_t operator()(NodeId id) const noexcept;
};

struct EdgeId {
  size_t value = NoId;
};

inline bool operator==(EdgeId lhs, EdgeId rhs);
inline bool operator!=(EdgeId lhs, EdgeId rhs);
inline bool operator<(EdgeId lhs, EdgeId rhs);

struct CladeIdx {
  size_t value = NoId;
};

inline bool operator==(CladeIdx lhs, CladeIdx rhs);
inline bool operator!=(CladeIdx lhs, CladeIdx rhs);
inline bool operator<(CladeIdx lhs, CladeIdx rhs);

namespace Transform {

inline auto GetParent();
inline auto GetChild();
inline auto GetId();
template <typename DAG>
auto ToNodes(DAG dag);
template <typename DAG>
auto ToEdges(DAG dag);

}  // namespace Transform

template <template <typename...> typename Template, size_t I, typename... Ts>
static constexpr auto select_argument();
template <template <typename...> typename Template, typename... Ts>
using select_argument_t = decltype(select_argument<Template, 0, Ts...>());