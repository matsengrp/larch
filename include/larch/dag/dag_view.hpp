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
  auto AddNode(NodeId id) const -> Node;
  auto AppendNode() const -> Node;
  auto AddEdge(EdgeId id, NodeId parent, NodeId child, CladeIdx clade) const -> Edge;
  auto AppendEdge(NodeId parent, NodeId child, CladeIdx clade) const -> Edge;
  void InitializeNodes(size_t nodes_count) const;
  void BuildConnections() const;
  void BuildConnectionsRaw() const;
  auto GetNodes() const;
  auto GetEdges() const;
  auto Get(NodeId id) const;
  auto Get(EdgeId id) const;
  size_t GetNodesCount() const;
  size_t GetEdgesCount() const;
  bool IsTree() const;
  bool HaveRoot() const;
  auto GetRoot() const;
  auto GetLeafs() const;

 private:
  DAG_FEATURE_FRIENDS;
  DAG_VIEW_FRIENDS;
  auto& GetStorage() const;
  Storage& storage_;
};