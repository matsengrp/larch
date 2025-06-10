#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Target, Component C>
FragmentElementsContainer<Target, C>::FragmentElementsContainer(
    Target target, std::vector<Id<C>>&& ids)
    : target_{target}, ids_{std::move(ids)} {
  std::vector<Id<C>> unique{ids_};
  unique |= ranges::actions::sort(std::less<Id<C>>{}) |
            ranges::actions::unique(std::equal_to<Id<C>>{});
  Assert(unique.size() == ids_.size());
  if constexpr (C == Component::Node) {
    for (auto i : ids_) {
      fragment_element_features_.insert({i, DAGNeighbors{}});
    }
  }
}

template <typename Target, Component C>
template <typename VT>
size_t FragmentElementsContainer<Target, C>::GetCount() const {
  return ids_.size();
}

template <typename Target, Component C>
template <typename Feature, typename E>
auto FragmentElementsContainer<Target, C>::GetFeatureStorage(Id<C> id, E) {
  Assert(ranges::contains(ids_, id));
  if constexpr (FeatureEquivalent<Feature, Neighbors>::value) {
    return std::ref(fragment_element_features_.at(id));
  } else {
    return target_.GetStorage().template GetFeatureStorage<Feature>(id);
  }
}

template <typename Target, Component C>
template <typename Feature, typename E>
auto FragmentElementsContainer<Target, C>::GetFeatureStorage(Id<C> id, E) const {
  Assert(ranges::contains(ids_, id));
  if constexpr (FeatureEquivalent<Feature, Neighbors>::value) {
    return std::cref(fragment_element_features_.at(id));
  } else {
    return target_.Const().GetStorage().template GetFeatureStorage<Feature>(id);
  }
}

template <typename Target, Component C>
template <typename Feature>
auto FragmentElementsContainer<Target, C>::GetFeatureExtraStorage() {
  return target_.GetStorage().template GetFeatureExtraStorage<C, Feature>();
}

template <typename Target, Component C>
template <typename Feature>
auto FragmentElementsContainer<Target, C>::GetFeatureExtraStorage() const {
  return target_.Const().GetStorage().template GetFeatureExtraStorage<C, Feature>();
}
