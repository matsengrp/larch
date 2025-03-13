#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
template <Component C, typename Feature>
inline constexpr bool ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                                       Cont>::contains_element_feature = [] {
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
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
template <template <typename, typename> typename Base>
DAGView<typename ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                                  Cont>::Self,
        Base>
ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase, Cont>::View() {
  return DAGView<Self, Base>{static_cast<Self&>(*this)};
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
template <template <typename, typename> typename Base>
DAGView<const typename ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                                        Cont>::Self,
        Base>
ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase, Cont>::View() const {
  return DAGView<const Self, Base>{static_cast<const Self&>(*this)};
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
NodeId
ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase, Cont>::AppendNode() {
  if constexpr (Cont == IdContinuity::Dense) {
    additional_node_features_storage_.push_back({});
  } else {
    std::ignore = additional_node_features_storage_[GetNextAvailableId<Component::Node,
                                                                       TargetView>()];
  }
  return GetTarget().AppendNode().GetId();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
EdgeId
ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase, Cont>::AppendEdge() {
  if constexpr (Cont == IdContinuity::Dense) {
    additional_edge_features_storage_.push_back({});
  } else {
    std::ignore = additional_edge_features_storage_[GetNextAvailableId<Component::Edge,
                                                                       TargetView>()];
  }
  return GetTarget().AppendEdge().GetId();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase, Cont>::AddNode(
    NodeId id) {
  std::ignore = additional_node_features_storage_[id];
  GetTarget().AddNode(id);
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase, Cont>::AddEdge(
    EdgeId id) {
  std::ignore = additional_edge_features_storage_[id];
  GetTarget().AddEdge(id);
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
template <typename VT>
size_t ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                        Cont>::GetNodesCount() const {
  return GetTarget().GetStorage().template GetNodesCount<VT>();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
template <typename VT>
size_t ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                        Cont>::GetEdgesCount() const {
  return GetTarget().GetStorage().template GetEdgesCount<VT>();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
template <typename VT>
auto ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase, Cont>::GetNodes()
    const {
  return GetTarget().GetStorage().template GetNodes<VT>();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
template <typename VT>
auto ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase, Cont>::GetEdges()
    const {
  return GetTarget().GetStorage().template GetEdges<VT>();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                      Cont>::InitializeNodes(size_t size) {
  GetTarget().InitializeNodes(size);
  additional_node_features_storage_.resize(size);
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                      Cont>::InitializeEdges(size_t size) {
  GetTarget().InitializeEdges(size);
  additional_edge_features_storage_.resize(size);
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                      Cont>::ClearNodes() {
  GetTarget().ClearNodes();
  additional_node_features_storage_.clear();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                      Cont>::ClearEdges() {
  GetTarget().GetStorage().ClearEdges();
  additional_edge_features_storage_.clear();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
template <typename Feature>
auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                       Cont>::GetFeatureStorage() {
  if constexpr (tuple_contains_v<decltype(additional_dag_features_storage_), Feature,
                                 FeatureEquivalent>) {
    return tuple_get<Feature, FeatureEquivalent>(additional_dag_features_storage_);
  } else {
    return GetTarget().template GetFeatureStorage<Feature>();
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
template <typename Feature>
const auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                             Cont>::GetFeatureStorage() const {
  if constexpr (tuple_contains_v<decltype(additional_dag_features_storage_), Feature,
                                 FeatureEquivalent>) {
    return tuple_get<Feature, FeatureEquivalent>(additional_dag_features_storage_);
  } else {
    return GetTarget().template GetFeatureStorage<Feature>();
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
template <typename F>
auto&& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                        Cont>::GetFeatureStorage(NodeId id) {
  if constexpr (tuple_contains_v<
                    std::remove_reference_t<
                        decltype(additional_node_features_storage_.at(NodeId{0}))>,
                    F, FeatureEquivalent>) {
    return tuple_get<F, FeatureEquivalent>(additional_node_features_storage_.at(id));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
template <typename F>
const auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                             Cont>::GetFeatureStorage(NodeId id) const {
  if constexpr (tuple_contains_v<typename decltype(additional_node_features_storage_)::
                                     value_type::second_type,
                                 F, FeatureEquivalent>) {
    return tuple_get<F, FeatureEquivalent>(additional_node_features_storage_.at(id));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
template <typename F>
auto&& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                        Cont>::GetFeatureStorage(EdgeId id) {
  if constexpr (tuple_contains_v<typename decltype(additional_edge_features_storage_)::
                                     value_type::second_type,
                                 F, FeatureEquivalent>) {
    return tuple_get<F, FeatureEquivalent>(additional_edge_features_storage_.at(id));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
template <typename F>
const auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                             Cont>::GetFeatureStorage(EdgeId id) const {
  if constexpr (tuple_contains_v<typename decltype(additional_edge_features_storage_)::
                                     value_type::second_type,
                                 F, FeatureEquivalent>) {
    return tuple_get<F, FeatureEquivalent>(additional_edge_features_storage_.at(id));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
template <Component C, typename F>
auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                       Cont>::GetFeatureExtraStorage() {
  if constexpr (std::remove_reference_t<Target>::template contains_element_feature<C,
                                                                                   F>) {
    return GetTarget().template GetFeatureExtraStorage<C, F>();
  } else {
    if constexpr (C == Component::Node) {
      return tuple_get<ExtraFeatureStorage<F>, FeatureEquivalent>(
          additional_node_extra_features_storage_);
    } else {
      return tuple_get<ExtraFeatureStorage<F>, FeatureEquivalent>(
          additional_edge_extra_features_storage_);
    }
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
template <Component C, typename F>
const auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                             Cont>::GetFeatureExtraStorage() const {
  if constexpr (Target::template contains_element_feature<C, F>) {
    return GetTarget().template GetFeatureExtraStorage<C, F>();
  } else {
    if constexpr (C == Component::Node) {
      return tuple_get<ExtraFeatureStorage<F>, FeatureEquivalent>(
          additional_node_extra_features_storage_);
    } else {
      return tuple_get<ExtraFeatureStorage<F>, FeatureEquivalent>(
          additional_edge_extra_features_storage_);
    }
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase, Cont>::ExtendDAGStorage(
    Target&& target)
    : target_{std::forward<Target>(target)} {
  additional_node_features_storage_.resize(GetTarget().GetNodesCount());
  additional_edge_features_storage_.resize(GetTarget().GetEdgesCount());
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
auto ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase,
                      Cont>::GetTarget() {
  return ViewOf(target_);
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2, template <typename, typename> typename ViewBase,
          IdContinuity Cont>
auto ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2, ViewBase, Cont>::GetTarget()
    const {
  return ViewOf(target_);
}
