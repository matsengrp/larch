#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename DAGType, typename... Features>
NodeView<DAGType, Features...>::NodeView(const DAGType& dag, NodeId id)
    : dag_{dag}, id_{id} {}

template <typename DAGType, typename... Features>
NodeView<DAGType, Features...>::operator NodeId() const {
  return id_;
}

template <typename DAGType, typename... Features>
auto& NodeView<DAGType, Features...>::GetDAG() const {
  return dag_;
}

template <typename DAGType, typename... Features>
NodeId NodeView<DAGType, Features...>::GetId() const {
  return id_;
}

template <typename DAGType, typename... Features>
auto NodeView<DAGType, Features...>::GetParents() const {
  return GetStorage().parents_ | Transform::ToEdges(dag_);
}

template <typename DAGType, typename... Features>
auto NodeView<DAGType, Features...>::GetClades() const {
  return GetStorage().clades_ |
         ranges::views::transform([*this](const std::vector<EdgeId>& clade) {
           return clade | Transform::ToEdges(dag_);
         });
}

template <typename DAGType, typename... Features>
auto NodeView<DAGType, Features...>::GetClade(CladeIdx clade) const {
  return GetStorage().clades_.at(clade.value) | Transform::ToEdges(dag_);
}

template <typename DAGType, typename... Features>
size_t NodeView<DAGType, Features...>::GetCladesCount() const {
  return GetStorage().clades_.size();
}

template <typename DAGType, typename... Features>
auto NodeView<DAGType, Features...>::GetChildren() const {
  return GetClades() | ranges::views::join;
}

// template <typename DAGType, typename... Features>
// auto NodeView<DAGType, Features...>::GetSingleParent() const;

template <typename DAGType, typename... Features>
auto NodeView<DAGType, Features...>::GetFirstChild() const {
  Assert(not IsLeaf());
  return *GetChildren().begin();
}

template <typename DAGType, typename... Features>
auto NodeView<DAGType, Features...>::GetFirstClade() const {
  return GetClade({0});
}

template <typename DAGType, typename... Features>
bool NodeView<DAGType, Features...>::IsRoot() const {
  return GetStorage().parents_.empty();
}

template <typename DAGType, typename... Features>
bool NodeView<DAGType, Features...>::IsLeaf() const {
  for (auto& clade : GetStorage().clades_) {
    if (not clade.empty()) {
      return false;
    }
  }
  return true;
}

template <typename DAGType, typename... Features>
auto& NodeView<DAGType, Features...>::GetStorage() const {
  return dag_.storage_.nodes_.nodes_.at(id_.value);
}