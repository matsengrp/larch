#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename DS>
struct DAGView
    : std::conditional_t<std::is_const_v<DS>,
                         typename DS::template ConstDAGViewBase<DAGView<DS>>,
                         typename DS::template MutableDAGViewBase<DAGView<DS>>> {
 public:
  using NodeView = ElementView<NodeId, DAGView<DS>>;
  using EdgeView = ElementView<EdgeId, DAGView<DS>>;
  using StorageType = DS;

  template <typename Id, typename CRTP>
  using ConstElementViewBase = typename DS::template ConstElementViewBase<Id, CRTP>;
  template <typename Id, typename CRTP>
  using MutableElementViewBase = typename DS::template MutableElementViewBase<Id, CRTP>;

  static const bool is_mutable;
  template <typename Id, typename Feature>
  static const bool contains_element_feature;

  explicit DAGView(DS& dag_storage);

  operator DAGView<const DS>() const;

  ElementView<NodeId, DAGView<DS>> Get(NodeId id) const;
  ElementView<EdgeId, DAGView<DS>> Get(EdgeId id) const;

  ElementView<NodeId, DAGView<DS>> AppendNode() const;
  ElementView<EdgeId, DAGView<DS>> AppendEdge() const;

  ElementView<NodeId, DAGView<DS>> AddNode(NodeId id);
  ElementView<EdgeId, DAGView<DS>> AddEdge(EdgeId id, NodeId parent, NodeId child,
                                           CladeIdx clade);  // TODO
  ElementView<EdgeId, DAGView<DS>> AppendEdge(NodeId parent, NodeId child,
                                              CladeIdx clade) const;  // TODO

  size_t GetNodesCount() const;
  size_t GetEdgesCount() const;

  auto GetNodes() const;
  auto GetEdges() const;

  void InitializeNodes(size_t size) const;

  template <typename F>
  auto& GetFeatureStorage() const;
  template <typename F>
  auto& GetFeatureStorage(NodeId id) const;
  template <typename F>
  auto& GetFeatureStorage(EdgeId id) const;
  template <typename Id, typename F>
  auto& GetFeatureExtraStorage() const;

 private:
  DS& dag_storage_;
};