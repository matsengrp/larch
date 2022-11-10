#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Storage, typename... Features>
Storage& DefaultNodesContainer<Storage, Features...>::AddNode(NodeId id) {
  return GetOrInsert(nodes_, id);
}

template <typename Storage, typename... Features>
void DefaultNodesContainer<Storage, Features...>::InitializeNodes(size_t nodes_count) {
  nodes_.resize(nodes_count);
}

template <typename Storage, typename... Features>
auto DefaultNodesContainer<Storage, Features...>::View() {
  return nodes_ | ranges::views::all;
}

template <typename Storage, typename... Features>
auto DefaultNodesContainer<Storage, Features...>::View() const {
  return nodes_ | ranges::views::all;
}

template <typename Storage, typename... Features>
Storage& DefaultNodesContainer<Storage, Features...>::NodeAt(NodeId id) {
  return nodes_.at(id.value);
}

template <typename Storage, typename... Features>
const Storage& DefaultNodesContainer<Storage, Features...>::NodeAt(NodeId id) const {
  return nodes_.at(id.value);
}

template <typename Storage, typename... Features>
size_t DefaultNodesContainer<Storage, Features...>::Count() const {
  return nodes_.size();
}