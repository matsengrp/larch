#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
template <Component C, typename Feature>
inline constexpr bool ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2,
                                       ViewBase>::contains_element_feature = [] {
  // NOLINTBEGIN
  if constexpr (TargetView::StorageType::template contains_element_feature<C,
                                                                           Feature>) {
    return true;
  } else {
    if constexpr (C == Component::Node) {
      return OnNodes::template contains_element_feature<Feature>;
    } else {
      return OnEdges::template contains_element_feature<Feature>;
    }
  }
  // NOLINTEND
}();

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
template <template <typename, typename> typename Base>
DAGView<typename ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::Self,
        Base>
ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::View() {
  return DAGView<Self, Base>{static_cast<Self&>(*this)};
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
template <template <typename, typename> typename Base>
DAGView<const typename ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2,
                                        ViewBase>::Self,
        Base>
ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::View() const {
  return DAGView<const Self, Base>{static_cast<const Self&>(*this)};
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
NodeId ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::AppendNode() {
  additional_node_features_storage_.push_back({});
  return GetTarget().AppendNode().GetId();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
EdgeId ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::AppendEdge() {
  additional_edge_features_storage_.push_back({});
  return GetTarget().AppendEdge().GetId();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::AddNode(
    NodeId id) {
  std::ignore = GetOrInsert(additional_node_features_storage_, id);
  GetTarget().AddNode(id);
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::AddEdge(
    EdgeId id) {
  std::ignore = GetOrInsert(additional_edge_features_storage_, id);
  GetTarget().AddEdge(id);
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
template <typename VT>
size_t ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::GetNodesCount()
    const {
  return GetTarget().GetStorage().template GetNodesCount<VT>();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
template <typename VT>
size_t ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::GetEdgesCount()
    const {
  return GetTarget().GetStorage().template GetEdgesCount<VT>();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
template <typename VT>
auto ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::GetNodes() const {
  return GetTarget().GetStorage().template GetNodes<VT>();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
template <typename VT>
auto ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::GetEdges() const {
  return GetTarget().GetStorage().template GetEdges<VT>();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::InitializeNodes(
    size_t size) {
  GetTarget().InitializeNodes(size);
  additional_node_features_storage_.resize(size);
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::InitializeEdges(
    size_t size) {
  GetTarget().InitializeEdges(size);
  additional_edge_features_storage_.resize(size);
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::ClearNodes() {
  GetTarget().ClearNodes();
  additional_node_features_storage_.clear();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::ClearEdges() {
  GetTarget().GetStorage().ClearEdges();
  additional_edge_features_storage_.clear();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
template <typename Feature>
auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2,
                       ViewBase>::GetFeatureStorage() {
  if constexpr (tuple_contains_v<decltype(additional_dag_features_storage_), Feature>) {
    return std::get<Feature>(additional_dag_features_storage_);
  } else {
    return GetTarget().template GetFeatureStorage<Feature>();
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
template <typename Feature>
const auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2,
                             ViewBase>::GetFeatureStorage() const {
  if constexpr (tuple_contains_v<decltype(additional_dag_features_storage_), Feature>) {
    return std::get<Feature>(additional_dag_features_storage_);
  } else {
    return GetTarget().template GetFeatureStorage<Feature>();
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
template <typename F>
auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2,
                       ViewBase>::GetFeatureStorage(NodeId id) {
  if constexpr (tuple_contains_v<std::remove_reference_t<
                                     decltype(additional_node_features_storage_.at(0))>,
                                 F>) {
    return std::get<F>(additional_node_features_storage_.at(id.value));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
template <typename F>
const auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2,
                             ViewBase>::GetFeatureStorage(NodeId id) const {
  if constexpr (tuple_contains_v<
                    typename decltype(additional_node_features_storage_)::value_type,
                    F>) {
    return std::get<F>(additional_node_features_storage_.at(id.value));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
template <typename F>
auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2,
                       ViewBase>::GetFeatureStorage(EdgeId id) {
  if constexpr (tuple_contains_v<
                    typename decltype(additional_edge_features_storage_)::value_type,
                    F>) {
    return std::get<F>(additional_edge_features_storage_.at(id.value));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
template <typename F>
const auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2,
                             ViewBase>::GetFeatureStorage(EdgeId id) const {
  if constexpr (tuple_contains_v<
                    typename decltype(additional_edge_features_storage_)::value_type,
                    F>) {
    return std::get<F>(additional_edge_features_storage_.at(id.value));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
template <Component C, typename F>
auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2,
                       ViewBase>::GetFeatureExtraStorage() {
  if constexpr (std::remove_reference_t<Target>::template contains_element_feature<C,
                                                                                   F>) {
    return GetTarget().template GetFeatureExtraStorage<C, F>();
  } else {
    if constexpr (C == Component::Node) {
      return std::get<ExtraFeatureStorage<F>>(additional_node_extra_features_storage_);
    } else {
      return std::get<ExtraFeatureStorage<F>>(additional_edge_extra_features_storage_);
    }
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
template <Component C, typename F>
const auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2,
                             ViewBase>::GetFeatureExtraStorage() const {
  if constexpr (Target::template contains_element_feature<C, F>) {
    return GetTarget().template GetFeatureExtraStorage<C, F>();
  } else {
    if constexpr (C == Component::Node) {
      return std::get<ExtraFeatureStorage<F>>(additional_node_extra_features_storage_);
    } else {
      return std::get<ExtraFeatureStorage<F>>(additional_edge_extra_features_storage_);
    }
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::ExtendDAGStorage(
    Target&& target)
    : target_{std::forward<Target>(target)} {
  additional_node_features_storage_.resize(GetTarget().GetNodesCount());
  additional_edge_features_storage_.resize(GetTarget().GetEdgesCount());
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
auto ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::GetTarget() {
  return ViewOf(target_);
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase>
auto ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase>::GetTarget()
    const {
  return ViewOf(target_);
}
