#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::ExtendDAGStorage(DV dv, Arg0, Arg1, Arg2)
    : target_dag_view_{dv} {
  additional_node_features_storage_.resize(target_dag_view_.NodesCount());
  additional_edge_features_storage_.resize(target_dag_view_.EdgesCount());
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
auto ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::View() {
  return DAGView<ExtendDAGStorage<DV, Arg0, Arg1, Arg2>>{*this};
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
NodeId ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::AppendNode() {
  additional_node_features_storage_.push_back({});
  return target_dag_view_.AppendNode().GetId();
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
EdgeId ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::AppendEdge() {
  additional_edge_features_storage_.push_back({});
  return target_dag_view_.AppendEdge().GetId();
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
void ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::AddNode(NodeId id) {
  GetOrInsert(additional_node_features_storage_, id);
  target_dag_view_.AddNode(id);
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
void ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::AddEdge(EdgeId id) {
  GetOrInsert(additional_edge_features_storage_, id);
  target_dag_view_.AddEdge(id);
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
auto ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::GetNodes() const {
  return target_dag_view_.GetNodes() |
         ranges::views::transform([*this, idx = size_t{}](auto&) mutable {
           return ElementView<NodeId, DV>{*this, {idx++}};
         });
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
auto ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::GetEdges() const {
  return target_dag_view_.GetEdges() |
         ranges::views::transform([*this, idx = size_t{}](auto&) mutable {
           return ElementView<EdgeId, DV>{*this, {idx++}};
         });
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
void ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::InitializeNodes(size_t size) const {
  target_dag_view_.InitializeNodes(size);
  additional_node_features_storage_.resize(size);
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
template <typename F>
auto& ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::GetFeatureStorage() {
  if constexpr (tuple_contains_v<decltype(additional_dag_features_storage_), F>) {
    return std::get<F>(additional_dag_features_storage_);
  } else {
    return target_dag_view_.template GetFeatureStorage<F>();
  }
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
template <typename F>
const auto& ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::GetFeatureStorage() const {
  if constexpr (tuple_contains_v<decltype(additional_dag_features_storage_), F>) {
    return std::get<F>(additional_dag_features_storage_);
  } else {
    return target_dag_view_.template GetFeatureStorage<F>();
  }
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
template <typename F>
auto& ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::GetFeatureStorage(NodeId id) {
  if constexpr (tuple_contains_v<
                    std::decay_t<decltype(additional_node_features_storage_.at(0))>,
                    F>) {
    return std::get<F>(additional_node_features_storage_.at(id.value));
  } else {
    return target_dag_view_.template GetFeatureStorage<F>(id);
  }
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
template <typename F>
const auto& ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::GetFeatureStorage(NodeId id) const {
  if constexpr (tuple_contains_v<
                    std::decay_t<decltype(additional_node_features_storage_.at(0))>,
                    F>) {
    return std::get<F>(additional_node_features_storage_.at(id.value));
  } else {
    return target_dag_view_.template GetFeatureStorage<F>(id);
  }
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
template <typename F>
auto& ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::GetFeatureStorage(EdgeId id) {
  if constexpr (tuple_contains_v<
                    std::decay_t<decltype(additional_edge_features_storage_.at(0))>,
                    F>) {
    return std::get<F>(additional_edge_features_storage_.at(id.value));
  } else {
    return target_dag_view_.template GetFeatureStorage<F>(id);
  }
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
template <typename F>
const auto& ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::GetFeatureStorage(EdgeId id) const {
  if constexpr (tuple_contains_v<
                    std::decay_t<decltype(additional_edge_features_storage_.at(0))>,
                    F>) {
    return std::get<F>(additional_edge_features_storage_.at(id.value));
  } else {
    return target_dag_view_.template GetFeatureStorage<F>(id);
  }
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
template <typename Id, typename F>
auto& ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::GetFeatureExtraStorage() {
  if constexpr (DV::template contains_element_feature<Id, F>) {
    return target_dag_view_.template GetFeatureExtraStorage<Id, F>();
  } else {
    if constexpr (std::is_same_v<Id, NodeId>) {
      return std::get<ExtraFeatureStorage<F>>(additional_node_extra_features_storage_);
    } else {
      return std::get<ExtraFeatureStorage<F>>(additional_edge_extra_features_storage_);
    }
  }
}

template <typename DV, typename Arg0, typename Arg1, typename Arg2>
template <typename Id, typename F>
const auto& ExtendDAGStorage<DV, Arg0, Arg1, Arg2>::GetFeatureExtraStorage() const {
  if constexpr (DV::template contains_element_feature<Id, F>) {
    return target_dag_view_.template GetFeatureExtraStorage<Id, F>();
  } else {
    if constexpr (std::is_same_v<Id, NodeId>) {
      return std::get<ExtraFeatureStorage<F>>(additional_node_extra_features_storage_);
    } else {
      return std::get<ExtraFeatureStorage<F>>(additional_edge_extra_features_storage_);
    }
  }
}