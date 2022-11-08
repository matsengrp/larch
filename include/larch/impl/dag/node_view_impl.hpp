#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename DAG>
NodeView<DAG>::NodeView(DAG dag, NodeId id) : dag_{dag}, id_{id} {}

template <typename DAG>
NodeView<DAG>::operator NodeId() {
  return id_;
}

template <typename DAG>
auto& NodeView<DAG>::GetDAG() {
  return dag_;
}

template <typename DAG>
NodeId NodeView<DAG>::GetId() {
  return id_;
}

template <typename DAG>
auto NodeView<DAG>::GetParents() {
  return GetStorage().parents_ | Transform::ToEdges(dag_);
}

template <typename DAG>
auto NodeView<DAG>::GetClades() {
  return GetStorage().clades_ |
         ranges::views::transform([*this](const std::vector<EdgeId>& clade) {
           return clade | Transform::ToEdges(dag_);
         });
}

template <typename DAG>
auto NodeView<DAG>::GetClade(CladeIdx clade) {
  return GetStorage().clades_.at(clade.value) | Transform::ToEdges(dag_);
}

template <typename DAG>
size_t NodeView<DAG>::GetCladesCount() {
  return GetStorage().clades_.size();
}

template <typename DAG>
auto NodeView<DAG>::GetChildren() {
  return GetClades() | ranges::views::join;
}

template <typename DAG>
auto NodeView<DAG>::GetSingleParent() {
  Assert(GetStorage().parents_.size() == 1);
  return *GetParents().begin();
}

template <typename DAG>
auto NodeView<DAG>::GetFirstChild() {
  Assert(not IsLeaf());
  return *GetChildren().begin();
}

template <typename DAG>
auto NodeView<DAG>::GetFirstClade() {
  return GetClade({0});
}

template <typename DAG>
bool NodeView<DAG>::IsRoot() {
  return GetStorage().parents_.empty();
}

template <typename DAG>
bool NodeView<DAG>::IsLeaf() {
  return ranges::all_of(GetStorage().clades_,
                        [](const auto& clade) { return clade.empty(); });
}

template <typename DAG>
auto& NodeView<DAG>::GetStorage() const {
  return dag_.storage_.nodes_.nodes_.at(id_.value);
}

template <typename DAG>
template <typename Feature>
auto& NodeView<DAG>::GetFeatureStorage() const {
  if constexpr (tuple_contians_v<typename DAG::NodesContainerFeatures, Feature>) {
    return std::get<Feature>(dag_.storage_.nodes_.features_).at(id_.value);
  } else {
    return std::get<Feature>(dag_.storage_.nodes_.nodes_.at(id_.value).features_);
  }
}