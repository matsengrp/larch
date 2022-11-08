#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename Storage, typename... Features>
class DAGView : public std::conditional_t<
                    std::is_const_v<Storage>,
                    FeatureReader<Features, DAGView<Storage, Features...>>,
                    FeatureWriter<Features, DAGView<Storage, Features...>>>... {
 public:
  constexpr static const bool is_mutable = not std::is_const_v<Storage>;
  using Immutable = DAGView<const Storage, Features...>;
  using Node = typename Storage::NodesContainerType::StorageType::template ViewType<
      DAGView<Storage, Features...>>;
  using Edge = typename Storage::EdgesContainerType::StorageType::template ViewType<
      DAGView<Storage, Features...>>;

  explicit DAGView(Storage& storage);
  operator Immutable() const;
  auto AddNode(NodeId id) -> Node;
  auto AppendNode() -> Node;
  auto AddEdge(EdgeId id, NodeId parent, NodeId child, CladeIdx clade) -> Edge;
  auto AppendEdge(NodeId parent, NodeId child, CladeIdx clade) -> Edge;
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