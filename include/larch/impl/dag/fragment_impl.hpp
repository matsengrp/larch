#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename DAG>
FragmentView<DAG>::FragmentView(FragmentStorage<DAG>& storage) : storage_{storage} {}

template <typename DAG>
void FragmentView<DAG>::AssertUA() const {
  storage_.dag_.AssertUA();
}

template <typename DAG>
size_t FragmentView<DAG>::GetNodesCount() const {
  return storage_.nodes_.size();
}

template <typename DAG>
size_t FragmentView<DAG>::GetEdgesCount() const {
  return storage_.edges_.size();
}

template <typename DAG>
auto FragmentView<DAG>::Get(NodeId id) const {
  return storage_.dag_.Get(id);
}

template <typename DAG>
auto FragmentView<DAG>::Get(EdgeId id) const {
  return storage_.dag_.Get(id);
}

template <typename DAG>
auto FragmentView<DAG>::GetNodes() const {
  return ranges::views::all(storage_.nodes_) | Transform::ToNodes(storage_.dag_);
}

template <typename DAG>
auto FragmentView<DAG>::GetEdges() const {
  return ranges::views::all(storage_.edges_) | Transform::ToEdges(storage_.dag_);
}

template <typename DAG>
auto FragmentView<DAG>::GetRoot() const {
  return storage_.dag_.GetRoot();
}
