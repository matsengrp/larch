#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

#include "larch/debug.hpp"

/**
  DAG is the main structure that holds node and edge data, and provides
  various queries.

  Populating with data should be performed by first adding all nodes by the
  AddNode() function, then adding all the edges with AddEdge() and finally
  calling BuildConnections().

  NodeId and EdgeId are strongly typed wrappers around size_t, and data is
  stored internally by the order of its IDs.

  Additional node and edge data may be stored in classes which extend DAG.
  See for example the class MADAG in `mutation_annotated_dag.hpp`.

  GetNodes() and GetEdges() returns a view into the corresponding elements,
  ordered by id.

*/
template <typename Storage, template <typename, typename> typename Base>
struct DAGView : Base<Storage, DAGView<Storage, Base>>::DAGViewBase {
  static_assert(Storage::role == Role::Storage);
  static_assert(Storage::component == Component::DAG);

 public:
  constexpr static const Component component = Component::DAG;
  constexpr static const Role role = Role::View;

  using NodeView = ElementView<Component::Node, DAGView>;
  using EdgeView = ElementView<Component::Edge, DAGView>;
  using BaseType = Base<Storage, DAGView>;
  using StorageType = Storage;
  using MutableType = DAGView<std::remove_const_t<Storage>, Base>;

  template <Component C>
  static constexpr IdContinuity id_continuity =
      StorageType::template Container<C>::id_continuity;

  static const bool is_mutable;
  template <Component C, typename Feature>
  static const bool contains_element_feature;

  explicit DAGView(Storage& dag_storage);

  DAGView(DAGView&& other) = default;
  DAGView(const DAGView& other) = default;

  DAGView& operator=(DAGView&& other) {
    dag_storage_ = std::addressof(other.GetStorage());
    return *this;
  }

  DAGView& operator=(const DAGView& other) {
    dag_storage_ = std::addressof(other.GetStorage());
    return *this;
  }

  bool operator==(const DAGView& other) const {
    LARCH_DEBUG_USE;
    return std::addressof(GetStorage()) == std::addressof(other.GetStorage());
  }

  operator DAGView<const Storage, Base>() const;

  DAGView<const Storage, Base> Const() const;

  /**
   * Get a Node or Edge view by its id
   * @{
   */
  ElementView<Component::Node, DAGView> Get(NodeId id) const;
  ElementView<Component::Edge, DAGView> Get(EdgeId id) const;
  /** @} */

  ElementView<Component::Node, DAGView> AppendNode() const;
  ElementView<Component::Edge, DAGView> AppendEdge() const;

  ElementView<Component::Node, DAGView> AddNode(NodeId id);
  ElementView<Component::Edge, DAGView> AddEdge(EdgeId id);

  ElementView<Component::Edge, DAGView> AddEdge(EdgeId id, NodeId parent, NodeId child,
                                                CladeIdx clade);  // TODO
  ElementView<Component::Edge, DAGView> AppendEdge(NodeId parent, NodeId child,
                                                   CladeIdx clade) const;  // TODO

  size_t GetNodesCount() const;
  size_t GetEdgesCount() const;
  template <Component C>
  size_t GetElementsCount() const {
    LARCH_DEBUG_USE;
    if constexpr (C == Component::Node) {
      return GetNodesCount();
    } else {
      return GetEdgesCount();
    }
  }
  bool empty() const;

  template <Component C, typename VT>
  Id<C> GetNextAvailableId() const {
    LARCH_DEBUG_USE;
    return GetStorage().template GetNextAvailableId<C, VT>();
  }

  template <typename VT>
  NodeId GetNextAvailableNodeId() const {
    LARCH_DEBUG_USE;
    return GetNextAvailableId<Component::Node, VT>();
  }

  template <typename VT>
  EdgeId GetNextAvailableEdgeId() const {
    LARCH_DEBUG_USE;
    return GetNextAvailableId<Component::Edge, VT>();
  }

  template <typename VT>
  bool ContainsId(NodeId id) const {
    LARCH_DEBUG_USE;
    return GetStorage().template ContainsId<VT>(id);
  }

  template <typename VT>
  bool ContainsId(EdgeId id) const {
    LARCH_DEBUG_USE;
    return GetStorage().template ContainsId<VT>(id);
  }

  /**
   * Return a range containing Node views for each node in the DAG
   */
  auto GetNodes() const;

  /**
   * Return a range containing Edge views for each edge in the DAG
   */
  auto GetEdges() const;

  void InitializeNodes(size_t size) const;
  void InitializeEdges(size_t size) const;

  template <typename Feature>
  auto GetFeatureStorage() const;
  template <typename Feature>
  auto GetFeatureStorage(NodeId id) const;
  template <typename Feature>
  auto GetFeatureStorage(EdgeId id) const;
  template <Component C, typename Feature>
  auto GetFeatureExtraStorage() const;

  Storage& GetStorage() const;

 private:
  Storage* dag_storage_ = nullptr;
  LARCH_DEBUG_THIS;
};
