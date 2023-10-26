#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
template <Component C, typename Feature>
inline constexpr bool ExtendDAGStorage<ShortName, Target, Arg0, Arg1,
                                       Arg2>::contains_element_feature = [] {
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
          typename Arg2>
DAGView<typename ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::Self>
ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::View() {
  return DAGView<Self>{static_cast<Self&>(*this)};
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
DAGView<const typename ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::Self>
ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::View() const {
  return DAGView<const Self>{static_cast<const Self&>(*this)};
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
NodeId ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::AppendNode() {
  additional_node_features_storage_.push_back({});
  return GetTarget().AppendNode().GetId();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
EdgeId ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::AppendEdge() {
  additional_edge_features_storage_.push_back({});
  return GetTarget().AppendEdge().GetId();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::AddNode(NodeId id) {
  std::ignore = GetOrInsert(additional_node_features_storage_, id);
  GetTarget().AddNode(id);
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::AddEdge(EdgeId id) {
  std::ignore = GetOrInsert(additional_edge_features_storage_, id);
  GetTarget().AddEdge(id);
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
size_t ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetNodesCount() const {
  return GetTarget().GetNodesCount();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
size_t ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetEdgesCount() const {
  return GetTarget().GetEdgesCount();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
auto ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetNodes() const {
  return ranges::views::indices(GetNodesCount()) |
         ranges::views::transform([this](size_t i) {
           return ElementView{this->View(), NodeId{i}};
         });
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
auto ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetEdges() const {
  return ranges::views::indices(GetEdgesCount()) |
         ranges::views::transform([this](size_t i) {
           return ElementView{this->View(), EdgeId{i}};
         });
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::InitializeNodes(
    size_t size) {
  GetTarget().InitializeNodes(size);
  additional_node_features_storage_.resize(size);
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::InitializeEdges(
    size_t size) {
  GetTarget().InitializeEdges(size);
  additional_edge_features_storage_.resize(size);
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::ClearNodes() {
  GetTarget().ClearNodes();
  additional_node_features_storage_.clear();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
void ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::ClearEdges() {
  GetTarget().GetStorage().ClearEdges();
  additional_edge_features_storage_.clear();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
template <typename F>
auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetFeatureStorage() {
  return additional_dag_features_storage_.template GetFeatureStorage<F>();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
template <typename F>
const auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetFeatureStorage()
    const {
  return additional_dag_features_storage_.template GetFeatureStorage<F>();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
template <typename F>
auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetFeatureStorage(
    NodeId id) {
  if constexpr (tuple_contains_v<std::remove_reference_t<
                                     decltype(additional_node_features_storage_.at(0))>,
                                 F>) {
    return std::get<F>(additional_node_features_storage_.at(id.value));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
template <typename F>
const auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetFeatureStorage(
    NodeId id) const {
  if constexpr (tuple_contains_v<
                    typename decltype(additional_node_features_storage_)::value_type,
                    F>) {
    return std::get<F>(additional_node_features_storage_.at(id.value));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
template <typename F>
auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetFeatureStorage(
    EdgeId id) {
  if constexpr (tuple_contains_v<
                    typename decltype(additional_edge_features_storage_)::value_type,
                    F>) {
    return std::get<F>(additional_edge_features_storage_.at(id.value));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
template <typename F>
const auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetFeatureStorage(
    EdgeId id) const {
  if constexpr (tuple_contains_v<
                    typename decltype(additional_edge_features_storage_)::value_type,
                    F>) {
    return std::get<F>(additional_edge_features_storage_.at(id.value));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
template <Component C, typename F>
auto& ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetFeatureExtraStorage() {
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
          typename Arg2>
template <Component C, typename F>
const auto&
ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetFeatureExtraStorage() const {
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
          typename Arg2>
template <Component C>
auto ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetContainer() ->
    typename TargetView::StorageType::template Container<C>& {
  return GetTarget().GetStorage().template GetContainer<C>();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
template <Component C>
auto ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetContainer() const
    -> const typename TargetView::StorageType::template Container<C>& {
  return GetTarget().GetStorage().template GetContainer<C>();
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::ExtendDAGStorage(Target&& target)
    : target_{std::forward<Target>(target)},
      additional_dag_features_storage_{ViewOf(target_)} {
  additional_node_features_storage_.resize(GetTarget().GetNodesCount());
  additional_edge_features_storage_.resize(GetTarget().GetEdgesCount());
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
auto ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetTarget() {
  return ViewOf(target_);
}

template <typename ShortName, typename Target, typename Arg0, typename Arg1,
          typename Arg2>
auto ExtendDAGStorage<ShortName, Target, Arg0, Arg1, Arg2>::GetTarget() const {
  return ViewOf(target_);
}
