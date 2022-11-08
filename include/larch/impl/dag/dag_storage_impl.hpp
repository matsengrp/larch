#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename NodesContainer, typename EdgesContainer, typename... Features>
auto DefaultDAGStorage<NodesContainer, EdgesContainer, Features...>::View() {
  return DAGView<DefaultDAGStorage<NodesContainer, EdgesContainer, Features...>,
                 Features...>{*this};
}