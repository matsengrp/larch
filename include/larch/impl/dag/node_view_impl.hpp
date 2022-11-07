#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename DAGType, typename... Features>
NodeView<DAGType, Features...>::NodeView(DAGType dag, NodeId id) : dag_{dag}, id_{id} {}

template <typename DAGType, typename... Features>
NodeView<DAGType, Features...>::operator NodeId() {
  return id_;
}

template <typename DAGType, typename... Features>
auto& NodeView<DAGType, Features...>::GetDAG() {
  return dag_;
}

template <typename DAGType, typename... Features>
NodeId NodeView<DAGType, Features...>::GetId() {
  return id_;
}

template <typename DAGType, typename... Features>
auto NodeView<DAGType, Features...>::GetParents() {
  return GetStorage().parents_ | Transform::ToEdges(dag_);
}

template <typename DAGType, typename... Features>
auto NodeView<DAGType, Features...>::GetClades() {
  return GetStorage().clades_ |
         ranges::views::transform([*this](const std::vector<EdgeId>& clade) {
           return clade | Transform::ToEdges(dag_);
         });
}

template <typename DAGType, typename... Features>
auto NodeView<DAGType, Features...>::GetClade(CladeIdx clade) {
  return GetStorage().clades_.at(clade.value) | Transform::ToEdges(dag_);
}

template <typename DAGType, typename... Features>
size_t NodeView<DAGType, Features...>::GetCladesCount() {
  return GetStorage().clades_.size();
}

template <typename DAGType, typename... Features>
auto NodeView<DAGType, Features...>::GetChildren() {
  return GetClades() | ranges::views::join;
}

template <typename DAGType, typename... Features>
auto NodeView<DAGType, Features...>::GetSingleParent() {
  Assert(GetStorage().parents_.size() == 1);
  return *GetParents().begin();
}

template <typename DAGType, typename... Features>
auto NodeView<DAGType, Features...>::GetFirstChild() {
  Assert(not IsLeaf());
  return *GetChildren().begin();
}

template <typename DAGType, typename... Features>
auto NodeView<DAGType, Features...>::GetFirstClade() {
  return GetClade({0});
}

template <typename DAGType, typename... Features>
bool NodeView<DAGType, Features...>::IsRoot() {
  return GetStorage().parents_.empty();
}

template <typename DAGType, typename... Features>
bool NodeView<DAGType, Features...>::IsLeaf() {
  return ranges::all_of(GetStorage().clades_,
                        [](const auto& clade) { return clade.empty(); });
}

template <typename DAGType, typename... Features>
auto& NodeView<DAGType, Features...>::GetStorage() const {
  return dag_.storage_.nodes_.nodes_.at(id_.value);
}
