#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename DAG>
Fragment<DAG>::Fragment(DAG dag, std::vector<NodeId>&& nodes,
                        std::vector<EdgeId>&& edges, NodeId root_node_id)
    : dag_{dag},
      nodes_{std::forward<std::vector<NodeId>>(nodes)},
      edges_{std::forward<std::vector<EdgeId>>(edges)},
      root_node_id_{root_node_id} {}

template <typename DAG>
void Fragment<DAG>::AssertUA() const {
  dag_.AssertUA();
}

template <typename DAG>
size_t Fragment<DAG>::GetNodesCount() const {
  return nodes_.size();
}

template <typename DAG>
size_t Fragment<DAG>::GetEdgesCount() const {
  return edges_.size();
}

template <typename DAG>
auto Fragment<DAG>::Get(NodeId id) const {
  return dag_.Get(id);
}

template <typename DAG>
auto Fragment<DAG>::Get(EdgeId id) const {
  return dag_.Get(id);
}

template <typename DAG>
auto Fragment<DAG>::GetNodes() const {
  return ranges::views::all(nodes_) | Transform::ToNodes(dag_);
}

template <typename DAG>
auto Fragment<DAG>::GetEdges() const {
  return ranges::views::all(edges_) | Transform::ToEdges(dag_);
}

template <typename DAG>
auto Fragment<DAG>::GetRoot() const {
  return dag_.Get(root_node_id_);
}

template <typename DAG>
bool Fragment<DAG>::IsTree() const {
  return GetNodesCount() == GetEdgesCount() + 1;
}
