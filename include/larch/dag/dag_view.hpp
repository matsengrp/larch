#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

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
template <typename Storage,
          template <typename, typename> typename Base = DefaultViewBase>
struct DAGView : Base<Storage, DAGView<Storage, Base>>::DAGViewBase {
 public:
  constexpr static const Component component = Component::DAG;
  constexpr static const Role role = Role::View;

  using NodeView = ElementView<Component::Node, DAGView<Storage, Base>>;
  using EdgeView = ElementView<Component::Edge, DAGView<Storage, Base>>;
  using BaseType = Base<Storage, DAGView<Storage, Base>>;
  using StorageType = Storage;
  using MutableType = DAGView<std::remove_const_t<Storage>>;

  static const bool is_mutable;
  template <Component C, typename Feature>
  static const bool contains_element_feature;

  explicit DAGView(Storage& dag_storage);

  operator DAGView<const Storage, Base>() const;

  DAGView<const Storage, Base> Const() const;

  /**
   * Get a Node or Edge view by its id
   * @{
   */
  ElementView<Component::Node, DAGView<Storage, Base>> Get(NodeId id) const;
  ElementView<Component::Edge, DAGView<Storage, Base>> Get(EdgeId id) const;
  /** @} */

  std::optional<ElementView<Component::Edge, DAGView<Storage, Base>>> FindEdge(
      NodeId parent_id, NodeId child_id) const;

  ElementView<Component::Node, DAGView<Storage, Base>> AppendNode() const;
  ElementView<Component::Edge, DAGView<Storage, Base>> AppendEdge() const;

  ElementView<Component::Node, DAGView<Storage, Base>> AddNode(NodeId id);
  ElementView<Component::Edge, DAGView<Storage, Base>> AddEdge(EdgeId id);

  ElementView<Component::Edge, DAGView<Storage, Base>> AddEdge(EdgeId id, NodeId parent,
                                                               NodeId child,
                                                               CladeIdx clade);  // TODO
  ElementView<Component::Edge, DAGView<Storage, Base>> AppendEdge(
      NodeId parent, NodeId child,
      CladeIdx clade) const;  // TODO

  size_t GetNodesCount() const;
  size_t GetEdgesCount() const;
  bool IsEmpty() const;

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

  void ClearNodes() const;
  void ClearEdges() const;

  template <typename Feature>
  auto& GetFeatureStorage() const;
  template <typename Feature>
  auto& GetFeatureStorage(NodeId id) const;
  template <typename Feature>
  auto& GetFeatureStorage(EdgeId id) const;
  template <Component C, typename Feature>
  auto& GetFeatureExtraStorage() const;

  const Storage& GetStorage() const;
  Storage& GetStorage();

 private:
  Storage& dag_storage_;
};
