#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Target, Component C>
FragmentElementsContainer<Target, C>::FragmentElementsContainer(
    Target target, std::vector<Id<C>>&& ids, NodeId root_node_id)
    : target_{target}, ids_{ids}, root_node_id_{root_node_id} {}

template <typename Target, Component C>
size_t FragmentElementsContainer<Target, C>::GetCount() const {
  return ids_.size();
}

template <typename Target, Component C>
template <typename Feature>
auto& FragmentElementsContainer<Target, C>::GetFeatureStorage(Id<C> id) {
  return target_.GetStorage().template GetFeatureStorage<Feature>(id);
}

template <typename Target, Component C>
template <typename Feature>
const auto& FragmentElementsContainer<Target, C>::GetFeatureStorage(Id<C> id) const {
  return target_.Const().GetStorage().template GetFeatureStorage<Feature>(id);
}

template <typename Target, Component C>
template <typename Feature>
auto& FragmentElementsContainer<Target, C>::GetFeatureExtraStorage() {
  return target_.GetStorage().template GetFeatureExtraStorage<C, Feature>();
}

template <typename Target, Component C>
template <typename Feature>
const auto& FragmentElementsContainer<Target, C>::GetFeatureExtraStorage() const {
  return target_.Const().GetStorage().template GetFeatureExtraStorage<C, Feature>();
}
