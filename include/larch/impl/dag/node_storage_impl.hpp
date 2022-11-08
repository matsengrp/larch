#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename... Features>
void DefaultNodeStorage<Features...>::ClearConnections() {
  parents_.clear();
  clades_.clear();
}

template <typename... Features>
void DefaultNodeStorage<Features...>::AddEdge(CladeIdx clade, EdgeId id,
                                              bool this_node_is_parent) {
  if (this_node_is_parent) {
    GetOrInsert(clades_, clade).push_back(id);
  } else {
    parents_.push_back(id);
  }
}