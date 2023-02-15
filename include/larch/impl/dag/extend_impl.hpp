#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
template <typename Id, typename Feature>
inline constexpr bool ExtendDAGStorage<Target, Arg0, Arg1,
                                       Arg2>::contains_element_feature = [] {
  if constexpr (TargetView::StorageType::template contains_element_feature<Id,
                                                                           Feature>) {
    return true;
  } else {
    if constexpr (std::is_same_v<Id, NodeId>) {
      return OnNodes::template contains_element_feature<Feature>;
    } else {
      return OnEdges::template contains_element_feature<Feature>;
    }
  }
}();

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::ExtendDAGStorage() = default;

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::ExtendDAGStorage(Target&& target)
    : target_{std::forward<Target>(target)} {
  additional_node_features_storage_.resize(GetTarget().GetNodesCount());
  additional_edge_features_storage_.resize(GetTarget().GetEdgesCount());
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
auto ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::View() {
  return DAGView<ExtendDAGStorage<Target, Arg0, Arg1, Arg2>>{*this};
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
auto ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::View() const {
  return DAGView<const ExtendDAGStorage<Target, Arg0, Arg1, Arg2>>{*this};
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
NodeId ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::AppendNode() {
  additional_node_features_storage_.push_back({});
  return GetTarget().AppendNode().GetId();
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
EdgeId ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::AppendEdge() {
  additional_edge_features_storage_.push_back({});
  return GetTarget().AppendEdge().GetId();
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
void ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::AddNode(NodeId id) {
  std::ignore = GetOrInsert(additional_node_features_storage_, id);
  GetTarget().AddNode(id);
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
void ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::AddEdge(EdgeId id) {
  std::ignore = GetOrInsert(additional_edge_features_storage_, id);
  GetTarget().AddEdge(id);
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
size_t ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::GetNodesCount() const {
  return GetTarget().GetNodesCount();
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
size_t ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::GetEdgesCount() const {
  return GetTarget().GetEdgesCount();
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
auto ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::GetNodes() const {
  return ranges::views::indices(GetNodesCount()) |
         ranges::views::transform([this](size_t i) {
           return ElementView{this->View(), NodeId{i}};
         });
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
auto ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::GetEdges() const {
  return ranges::views::indices(GetEdgesCount()) |
         ranges::views::transform([this](size_t i) {
           return ElementView{this->View(), EdgeId{i}};
         });
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
void ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::InitializeNodes(size_t size) {
  GetTarget().InitializeNodes(size);
  additional_node_features_storage_.resize(size);
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
template <typename F>
auto& ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::GetFeatureStorage() {
  if constexpr (tuple_contains_v<decltype(additional_dag_features_storage_), F>) {
    return std::get<F>(additional_dag_features_storage_);
  } else {
    return GetTarget().template GetFeatureStorage<F>();
  }
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
template <typename F>
const auto& ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::GetFeatureStorage() const {
  if constexpr (tuple_contains_v<decltype(additional_dag_features_storage_), F>) {
    return std::get<F>(additional_dag_features_storage_);
  } else {
    return GetTarget().template GetFeatureStorage<F>();
  }
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
template <typename F>
auto& ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::GetFeatureStorage(NodeId id) {
  if constexpr (tuple_contains_v<
                    std::decay_t<decltype(additional_node_features_storage_.at(0))>,
                    F>) {
    return std::get<F>(additional_node_features_storage_.at(id.value));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
template <typename F>
const auto& ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::GetFeatureStorage(
    NodeId id) const {
  if constexpr (tuple_contains_v<
                    std::decay_t<decltype(additional_node_features_storage_.at(0))>,
                    F>) {
    return std::get<F>(additional_node_features_storage_.at(id.value));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
template <typename F>
auto& ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::GetFeatureStorage(EdgeId id) {
  if constexpr (tuple_contains_v<
                    std::decay_t<decltype(additional_edge_features_storage_.at(0))>,
                    F>) {
    return std::get<F>(additional_edge_features_storage_.at(id.value));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
template <typename F>
const auto& ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::GetFeatureStorage(
    EdgeId id) const {
  if constexpr (tuple_contains_v<
                    std::decay_t<decltype(additional_edge_features_storage_.at(0))>,
                    F>) {
    return std::get<F>(additional_edge_features_storage_.at(id.value));
  } else {
    return GetTarget().template GetFeatureStorage<F>(id);
  }
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
template <typename Id, typename F>
auto& ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::GetFeatureExtraStorage() {
  if constexpr (std::decay_t<Target>::template contains_element_feature<Id, F>) {
    return GetTarget().template GetFeatureExtraStorage<Id, F>();
  } else {
    if constexpr (std::is_same_v<Id, NodeId>) {
      return std::get<ExtraFeatureStorage<F>>(additional_node_extra_features_storage_);
    } else {
      return std::get<ExtraFeatureStorage<F>>(additional_edge_extra_features_storage_);
    }
  }
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
template <typename Id, typename F>
const auto& ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::GetFeatureExtraStorage() const {
  if constexpr (Target::template contains_element_feature<Id, F>) {
    return GetTarget().template GetFeatureExtraStorage<Id, F>();
  } else {
    if constexpr (std::is_same_v<Id, NodeId>) {
      return std::get<ExtraFeatureStorage<F>>(additional_node_extra_features_storage_);
    } else {
      return std::get<ExtraFeatureStorage<F>>(additional_edge_extra_features_storage_);
    }
  }
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
auto ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::GetTarget() {
  return ViewOf(target_);
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
auto ExtendDAGStorage<Target, Arg0, Arg1, Arg2>::GetTarget() const {
  return ViewOf(target_);
}

template <typename Target, typename Arg0, typename Arg1, typename Arg2>
auto ExtendStorage(Target&& target, Arg0, Arg1, Arg2) {
  return ExtendDAGStorage<Target, Arg0, Arg1, Arg2>{std::forward<Target>(target)};
}