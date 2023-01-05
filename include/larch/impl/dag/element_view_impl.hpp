#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Id, typename DV>
ElementView<Id, DV>::ElementView(DV dag_view, Id id) : dag_view_{dag_view}, id_{id} {
  Assert(id.value not_eq NoId);
}

template <typename Id, typename DV>
ElementView<Id, DV>::operator Id() const {
  return GetId();
}

template <typename Id, typename DV>
DV ElementView<Id, DV>::GetDAG() const {
  return dag_view_;
}

template <typename Id, typename DV>
Id ElementView<Id, DV>::GetId() const {
  return id_;
}

template <typename Id, typename DV>
template <typename F>
auto& ElementView<Id, DV>::GetFeatureStorage() const {
  return dag_view_.template GetFeatureStorage<F>(id_);
}

template <typename Id, typename DV>
template <typename F>
auto& ElementView<Id, DV>::GetFeatureExtraStorage() const {
  return dag_view_.template GetFeatureExtraStorage<Id, F>();
}