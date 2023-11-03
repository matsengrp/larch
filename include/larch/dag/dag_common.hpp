#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/*
+--------------------------------------------------+
| +--------------------+   +---------------------+ |
| | +---------------+  |   |  +---------------+  | |
| | | Node: Element |  |   |  | Edge: Element |  | |
| | +---------------+  |   |  +---------------+  | |
| |                    |   |                     | |
| |  Nodes: Container  |   |  Edges: Container   | |
| +--------------------+   +---------------------+ |
|                                                  |
|                       DAG                        |
+--------------------------------------------------+
*/

enum class Component { Node, Edge, DAG };
enum class Role { Storage, View };

struct NodeId;
struct EdgeId;

template <Component C>
struct IdType;

template <>
struct IdType<Component::Node> {
  using type = NodeId;
};

template <>
struct IdType<Component::Edge> {
  using type = EdgeId;
};

template <Component C>
using Id = typename IdType<C>::type;

template <typename Id>
struct ComponentType;

template <>
struct ComponentType<NodeId> {
  constexpr static const auto value = Component::Node;
};

template <>
struct ComponentType<EdgeId> {
  constexpr static const auto value = Component::Edge;
};

template <typename Id>
inline constexpr const auto ComponentOf = ComponentType<Id>::value;

/**
 * Used by specialization on the Feature parameter to add functions to
 * an attachable feature. Functions declared in FeatureConstView are
 * for read-only access to the Feature's storage, and in FeatureMutableView
 * are for modifying access.
 *
 * The Tag template parameter equals Feature by default, and can be customized
 * for behavior modifiers like Deduplicate.
 * @{
 */
template <typename Feature, typename CRTP, typename Tag = Feature>
struct FeatureConstView;
template <typename Feature, typename CRTP, typename Tag = Feature>
struct FeatureMutableView;
/** @} */

/**
 * A per-element feature can have some additional functions, that are exposed
 * in the DAG view, rather than the element view.
 * @{
 */
template <typename Feature, typename CRTP>
struct ExtraFeatureConstView {};
template <typename Feature, typename CRTP>
struct ExtraFeatureMutableView {};
/** @} */

/**
 * Used by specialization on the Feature parameter to add extra storage for
 * features that are attached to node or edge (not the DAG itself). Extra storage
 * is not per-element, but global for the elements container.
 */
template <typename Feature>
struct ExtraFeatureStorage {};

template <Component C, typename DAGViewType>
struct ElementView;

template <typename DAGStorageType, typename DAGViewType>
struct DefaultViewBase {
 private:
  struct MutableDAGViewBase : DAGStorageType::template ConstDAGViewBase<DAGViewType>,
                              DAGStorageType::template MutableDAGViewBase<DAGViewType> {
  };

  template <Component C>
  struct MutableElementViewBase
      : DAGStorageType::template ConstElementViewBase<C, ElementView<C, DAGViewType>>,
        DAGStorageType::template MutableElementViewBase<C,
                                                        ElementView<C, DAGViewType>> {
    using DAGStorageType::template MutableElementViewBase<
        C, ElementView<C, DAGViewType>>::operator=;
  };

 public:
  using DAGViewBase = std::conditional_t<
      std::is_const_v<DAGStorageType>,
      typename DAGStorageType::template ConstDAGViewBase<DAGViewType>,
      MutableDAGViewBase>;

  template <Component C>
  using ElementViewBase =
      std::conditional_t<std::is_const_v<DAGStorageType>,
                         typename DAGStorageType::template ConstElementViewBase<
                             C, ElementView<C, DAGViewType>>,
                         MutableElementViewBase<C>>;
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

inline std::ostream& operator<<(std::ostream& os, NodeId node_id);
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

inline std::ostream& operator<<(std::ostream& os, EdgeId edge_id);
inline bool operator==(EdgeId lhs, EdgeId rhs);
inline bool operator!=(EdgeId lhs, EdgeId rhs);
inline bool operator<(EdgeId lhs, EdgeId rhs);

template <>
struct std::hash<EdgeId> {
  inline size_t operator()(EdgeId id) const noexcept;
};
struct CladeIdx {
  size_t value = NoId;
};

inline std::ostream& operator<<(std::ostream& os, CladeIdx clade_id);
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
inline auto ToConst();
inline auto ToView();

}  // namespace Transform

template <template <typename...> typename Template, size_t I, typename... Ts>
static constexpr auto select_argument();
template <template <typename...> typename Template, typename... Ts>
using select_argument_t = decltype(select_argument<Template, 0, Ts...>());

template <typename... Ts>
struct CombineBases : Ts... {
  static_assert(std::conjunction_v<std::is_empty<Ts>...>);
};

template <typename... Ts>
struct CombineBases<std::tuple<Ts...>> : Ts... {
  static_assert(std::conjunction_v<std::is_empty<Ts>...>);
};

template <typename T>
auto ViewOf(T&& dag);

template <typename, template <typename, typename> typename = DefaultViewBase>
struct DAGView;

template <typename...>
struct ExtraStorage;

template <typename, typename, typename, typename>
struct DAGStorage;

template <Component, typename, typename...>
struct ElementsContainer;

template <typename...>
struct ElementStorage;

struct Neighbors;

struct Endpoints;

struct Connections;

struct DefaultDAGStorage;

template <typename T, typename = void>
struct ViewTypeOf;

template <typename T>
struct ViewTypeOf<T, std::enable_if_t<T::role == Role::View>> {
  static_assert(T::component == Component::DAG);
  using type = T;
};

template <typename T>
struct ViewTypeOf<T, std::enable_if_t<T::role == Role::Storage>> {
  static_assert(T::component == Component::DAG);
  using type = decltype(std::declval<T>().View());
};

template <typename>
struct LongNameOf;

#define SHORT_NAME(x)                                \
  MOVE_ONLY(x);                                      \
  static x EmptyDefault() {                          \
    static_assert(Target::role == Role::Storage);    \
    return x{Target::EmptyDefault()};                \
  }                                                  \
                                                     \
  static x Consume(Target&& target) {                \
    static_assert(Target::role == Role::Storage);    \
    return x{std::move(target)};                     \
  }                                                  \
                                                     \
  static x FromView(const Target& target) {          \
    static_assert(Target::role == Role::View);       \
    return x{Target{target}};                        \
  }                                                  \
  using LongNameType = typename LongNameOf<x>::type; \
  using LongNameType::LongNameType

template <typename Short, typename Long>
struct IsNameCorrect {
  static constexpr bool value = std::is_same_v<typename LongNameOf<Short>::type, Long>;
};

template <typename Long>
struct IsNameCorrect<void, Long> {
  static constexpr bool value = true;
};
