#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename CRTP, typename Feature, typename Tag>
auto& GetFeatureStorage(const FeatureMutableView<Feature, CRTP, Tag>* feature) {
  return static_cast<const CRTP&>(*feature).template GetFeatureStorage<Tag>();
}

template <typename CRTP, typename Feature, typename Tag>
const auto& GetFeatureStorage(const FeatureConstView<Feature, CRTP, Tag>* feature) {
  return static_cast<const CRTP&>(*feature).Const().template GetFeatureStorage<Tag>();
}

std::ostream& operator<<(std::ostream& os, const NodeId node_id) {
  os << "NodeId::" << node_id.value;
  return os;
}
bool operator==(NodeId lhs, NodeId rhs) { return lhs.value == rhs.value; }
bool operator!=(NodeId lhs, NodeId rhs) { return lhs.value != rhs.value; }
bool operator<(NodeId lhs, NodeId rhs) { return lhs.value < rhs.value; }

size_t std::hash<NodeId>::operator()(NodeId id) const noexcept { return id.value; }

std::ostream& operator<<(std::ostream& os, const EdgeId edge_id) {
  os << "EdgeId::" << edge_id.value;
  return os;
}
bool operator==(EdgeId lhs, EdgeId rhs) { return lhs.value == rhs.value; }
bool operator!=(EdgeId lhs, EdgeId rhs) { return lhs.value != rhs.value; }
bool operator<(EdgeId lhs, EdgeId rhs) { return lhs.value < rhs.value; }
size_t std::hash<EdgeId>::operator()(EdgeId id) const noexcept { return id.value; }

std::ostream& operator<<(std::ostream& os, const CladeIdx clade_id) {
  os << "CladeId::" << clade_id.value;
  return os;
}
bool operator==(CladeIdx lhs, CladeIdx rhs) { return lhs.value == rhs.value; }
bool operator!=(CladeIdx lhs, CladeIdx rhs) { return lhs.value != rhs.value; }
bool operator<(CladeIdx lhs, CladeIdx rhs) { return lhs.value < rhs.value; }

namespace Transform {

auto GetParent() {
  return ranges::views::transform([](auto&& i) { return i.GetParent(); });
}
auto GetChild() {
  return ranges::views::transform([](auto&& i) { return i.GetChild(); });
}
auto GetId() {
  return ranges::views::transform([](auto&& i) { return i.GetId(); });
}
template <typename DAG>
auto ToNodes(DAG dag) {
  return ranges::views::transform([dag](auto&& i) {
    return typename DAG::NodeView{dag, i};
  });
}
template <typename DAG>
auto ToEdges(DAG dag) {
  return ranges::views::transform([dag](auto&& i) {
    return typename DAG::EdgeView{dag, i};
  });
}
auto ToConst() {
  return ranges::views::transform([](auto&& i) { return i.Const(); });
}

auto ToView() {
  return ranges::views::transform([](auto&& i) { return i.View(); });
}

}  // namespace Transform

template <template <typename...> typename Template, size_t I, typename... Ts>
static constexpr auto select_argument() {
  if constexpr (I < std::tuple_size_v<std::tuple<Ts...>>) {
    if constexpr (is_specialization_v<std::tuple_element_t<I, std::tuple<Ts...>>,
                                      Template>) {
      return std::get<I>(std::tuple<Ts...>{});
    } else {
      return select_argument<Template, I + 1, Ts...>();
    }
  } else {
    return Template<>{};
  }
}

template <typename T>
auto ViewOf(T&& dag) {
  if constexpr (std::decay_t<T>::role == Role::View) {
    return dag;
  } else {
    return dag.View();
  }
}
