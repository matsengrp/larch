#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename Storage>
class DAGView : public Storage::template DAGViewBaseType<Storage> {
 public:
  constexpr static const bool is_mutable = not std::is_const_v<Storage>;
  using Immutable = DAGView<const Storage>;
  using StorageType = Storage;
  using Node = NodeView<DAGView<Storage>>;
  using Edge = EdgeView<DAGView<Storage>>;
  using NodesContainerFeatures = typename Storage::NodesContainerType::FeaturesType;
  using EdgesContainerFeatures = typename Storage::EdgesContainerType::FeaturesType;

  explicit DAGView(Storage& storage);
  operator Immutable() const;
  template <typename Feature>
  auto& Get();
  template <typename Feature>
  void Set(Feature&& feature);
  Node AddNode(NodeId id);
  Node AppendNode();
  Edge AddEdge(EdgeId id, NodeId parent, NodeId child, CladeIdx clade);
  Edge AppendEdge(NodeId parent, NodeId child, CladeIdx clade);
  void InitializeNodes(size_t nodes_count);
  void BuildConnections();
  void BuildConnectionsRaw();
  auto GetNodes();
  auto GetEdges();
  auto Get(NodeId id);
  auto Get(EdgeId id);
  size_t GetNodesCount();
  size_t GetEdgesCount();
  bool IsTree();
  bool HaveRoot();
  auto GetRoot();
  auto GetLeafs();

 private:
  DAG_FEATURE_FRIENDS;
  DAG_VIEW_FRIENDS;
  auto& GetStorage() const;
  Storage& storage_;
};