#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Id, typename DAGViewType>
template <typename Feature>
inline constexpr bool ElementView<Id, DAGViewType>::contains_feature =
    DAGViewType::template contains_element_feature<Id, Feature>;

template <typename Id, typename DAGViewType>
ElementView<Id, DAGViewType>::ElementView(DAGViewType dag_view, Id id)
    : dag_view_{dag_view}, id_{id} {
  Assert(id.value not_eq NoId);
}

template <typename Id, typename DAGViewType>
ElementView<Id, DAGViewType>::operator Id() const {
  return GetId();
}

template <typename Id, typename DAGViewType>
auto ElementView<Id, DAGViewType>::Const() const {
  return dag_view_.Const().Get(id_);
}

template <typename Id, typename DAGViewType>
DAGViewType ElementView<Id, DAGViewType>::GetDAG() const {
  return dag_view_;
}

template <typename Id, typename DAGViewType>
Id ElementView<Id, DAGViewType>::GetId() const {
  return id_;
}

template <typename Id, typename DAGViewType>
template <typename Feature>
auto& ElementView<Id, DAGViewType>::GetFeatureStorage() const {
  return dag_view_.template GetFeatureStorage<Feature>(id_);
}

template <typename Id, typename DAGViewType>
template <typename Feature>
auto& ElementView<Id, DAGViewType>::GetFeatureExtraStorage() const {
  return dag_view_.template GetFeatureExtraStorage<Id, Feature>();
}
