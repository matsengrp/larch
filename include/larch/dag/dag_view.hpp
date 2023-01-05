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

template <typename DagStorageT>
struct DAGView
    : std::conditional_t<
          std::is_const_v<DagStorageT>,
          typename DagStorageT::template ConstDAGViewBase<DAGView<DagStorageT>>,
          typename DagStorageT::template MutableDAGViewBase<DAGView<DagStorageT>>> {
 public:
  using NodeView = ElementView<NodeId, DAGView<DagStorageT>>;
  using EdgeView = ElementView<EdgeId, DAGView<DagStorageT>>;
  using StorageType = DagStorageT;
  using MutableType = DAGView<std::remove_const_t<DagStorageT>>;

  template <typename Id, typename CRTP>
  using ConstElementViewBase =
      typename DagStorageT::template ConstElementViewBase<Id, CRTP>;
  template <typename Id, typename CRTP>
  using MutableElementViewBase =
      typename DagStorageT::template MutableElementViewBase<Id, CRTP>;

  static const bool is_mutable;
  template <typename Id, typename Feature>
  static const bool contains_element_feature;

  explicit DAGView(DagStorageT& dag_storage);

  operator DAGView<const DagStorageT>() const;

  /**
   * Get a Node or Edge view by its id
   * @{
   */
  ElementView<NodeId, DAGView<DagStorageT>> Get(NodeId id) const;
  ElementView<EdgeId, DAGView<DagStorageT>> Get(EdgeId id) const;
  /** @} */

  ElementView<NodeId, DAGView<DagStorageT>> AppendNode() const;
  ElementView<EdgeId, DAGView<DagStorageT>> AppendEdge() const;

  ElementView<NodeId, DAGView<DagStorageT>> AddNode(NodeId id);
  ElementView<EdgeId, DAGView<DagStorageT>> AddEdge(EdgeId id, NodeId parent,
                                                    NodeId child,
                                                    CladeIdx clade);  // TODO
  ElementView<EdgeId, DAGView<DagStorageT>> AppendEdge(NodeId parent, NodeId child,
                                                       CladeIdx clade) const;  // TODO

  size_t GetNodesCount() const;
  size_t GetEdgesCount() const;

  /**
   * Return a range containing Node views for each node in the DAG
   */
  auto GetNodes() const;

  /**
   * Return a range containing Edge views for each edge in the DAG
   */
  auto GetEdges() const;

  void InitializeNodes(size_t size) const;

  template <typename Feature>
  auto& GetFeatureStorage() const;
  template <typename Feature>
  auto& GetFeatureStorage(NodeId id) const;
  template <typename Feature>
  auto& GetFeatureStorage(EdgeId id) const;
  template <typename Id, typename Feature>
  auto& GetFeatureExtraStorage() const;

 private:
  DagStorageT& dag_storage_;
};