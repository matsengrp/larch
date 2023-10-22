#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename Target, Component C>
struct FragmentElementsContainer {
  static_assert(Target::role == Role::View);
  static_assert(Target::component == Component::DAG);

  using TargetStorage = typename Target::StorageType;
  using TargetContainer = typename TargetStorage::template Container<C>;
  using FeatureTypes = typename TargetStorage::FeatureTypes;
  using AllFeatureTypes =
      std::conditional_t<C == Component::Node, typename TargetStorage::AllNodeFeatures,
                         typename TargetStorage::AllEdgeFeatures>;

  template <typename Feature>
  static const bool contains_element_feature =
      TargetStorage::template contains_element_feature<C, Feature>;

  template <typename CRTP>
  struct ConstElementViewBase : TargetStorage::ConstElementViewBase<C, CRTP> {};
  template <typename CRTP>
  struct MutableElementViewBase : TargetStorage::MutableElementViewBase<C, CRTP> {
    using TargetStorage::template MutableElementViewBase<C, CRTP>::operator=;
  };

  template <typename CRTP>
  struct ExtraConstElementViewBase : TargetContainer::ExtraConstElementViewBase<CRTP> {
  };
  template <typename CRTP>
  struct ExtraMutableElementViewBase
      : TargetContainer::ExtraMutableElementViewBase<CRTP> {};

  FragmentElementsContainer(Target target, std::vector<Id<C>>&& ids,
                            NodeId root_node_id);
  MOVE_ONLY(FragmentElementsContainer);

  size_t GetCount() const;

  // Id<C> Append();
  // void Add(Id<C> id);
  // void Initialize(size_t size);

  template <typename Feature>
  auto& GetFeatureStorage(Id<C> id);
  template <typename Feature>
  const auto& GetFeatureStorage(Id<C> id) const;

  template <typename Feature>
  auto& GetFeatureExtraStorage();
  template <typename Feature>
  const auto& GetFeatureExtraStorage() const;

  auto All() const { return ids_ | ranges::views::all; }

 private:
  Target target_;
  const std::vector<Id<C>> ids_;
  const NodeId root_node_id_;
};

template <typename Target>
struct FragmentExtraStorage {
  constexpr static const Component component = Component::DAG;
  constexpr static const Role role = Role::Storage;

  static_assert(Target::role == Role::View);
  static_assert(Target::component == Component::DAG);

  explicit FragmentExtraStorage(Target target) : target_{target} {}
  MOVE_ONLY(FragmentExtraStorage);

  using FeatureTypes = typename Target::StorageType::ExtraStorageType::FeatureTypes;

  template <template <typename, typename> typename T, typename CRTP>
  using Base = typename Target::StorageType::ExtraStorageType::template Base<T, CRTP>;

  template <typename Feature>
  auto& GetFeatureStorage() {
    return target_.GetStorage().template GetFeatureStorage<Feature>();
  }

  template <typename Feature>
  const auto& GetFeatureStorage() const {
    return target_.GetStorage().template GetFeatureStorage<Feature>();
  }

  template <typename Storage>
  auto& GetTargetStorage(Storage&) const {
    return target_.GetStorage();
  }

  template <typename Storage>
  auto& GetTargetStorage(Storage&) {
    return target_.GetStorage();
  }

 private:
  Target target_;
};

template <typename Target>
using FragmentStorageFor =
    DAGStorage<FragmentElementsContainer<Target, Component::Node>,
               FragmentElementsContainer<Target, Component::Edge>,
               FragmentExtraStorage<Target>>;

template <typename Target>
FragmentStorageFor<Target> FragmentStorage(Target target, std::vector<NodeId>&& nodes,
                                           std::vector<EdgeId>&& edges,
                                           NodeId root_node_id) {
  return {FragmentElementsContainer<Target, Component::Node>{target, std::move(nodes),
                                                             root_node_id},
          FragmentElementsContainer<Target, Component::Edge>{target, std::move(edges),
                                                             root_node_id},
          FragmentExtraStorage<Target>{target}};
}
