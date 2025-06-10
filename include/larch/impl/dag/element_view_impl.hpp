#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <Component C, typename DAGViewType>
template <typename Feature>
inline constexpr bool ElementView<C, DAGViewType>::contains_feature =
    DAGViewType::template contains_element_feature<C, Feature>;

template <Component C, typename DAGViewType>
ElementView<C, DAGViewType>::ElementView(DAGViewType dag_view, Id<C> id)
    : dag_view_{dag_view}, id_{id} {
  Assert(id.value not_eq NoId);
  // Assert(id.value < dag_view.template GetElementsCount<C>());
}

template <Component C, typename DAGViewType>
ElementView<C, DAGViewType>::operator Id<C>() const {
  // LARCH_DEBUG_USE;
  return GetId();
}

template <Component C, typename DAGViewType>
auto ElementView<C, DAGViewType>::Const() const {
  // LARCH_DEBUG_USE;
  return dag_view_.Const().Get(id_);
}

template <Component C, typename DAGViewType>
DAGViewType ElementView<C, DAGViewType>::GetDAG() const {
  // LARCH_DEBUG_USE;
  return dag_view_;
}

template <Component C, typename DAGViewType>
Id<C> ElementView<C, DAGViewType>::GetId() const {
  // LARCH_DEBUG_USE;
  return id_;
}

template <Component C, typename DAGViewType>
template <typename Feature>
auto ElementView<C, DAGViewType>::GetFeatureStorage() const {
  // LARCH_DEBUG_USE;
  return dag_view_.template GetFeatureStorage<Feature>(id_);
}

template <Component C, typename DAGViewType>
template <typename Feature>
auto ElementView<C, DAGViewType>::GetFeatureExtraStorage() const {
  // LARCH_DEBUG_USE;
  return dag_view_.template GetFeatureExtraStorage<C, Feature>();
}
