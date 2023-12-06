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

  static constexpr IdContinuity id_continuity = IdContinuity::Sparse;

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

  FragmentElementsContainer(Target target, std::vector<Id<C>>&& ids);
  MOVE_ONLY(FragmentElementsContainer);

  size_t GetCount() const;

  Id<C> GetNextAvailableId() const { return target_.template GetNextAvailableId<C>(); }

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
  std::conditional_t<
      C == Component::Node,
      IdContainer<NodeId, Neighbors, IdContinuity::Sparse, Ordering::Ordered>,
      std::tuple<>>
      fragment_element_features_;
};

template <typename Target>
struct FragmentExtraStorage {
  constexpr static const Component component = Component::DAG;
  constexpr static const Role role = Role::Storage;

  static_assert(Target::role == Role::View);
  static_assert(Target::component == Component::DAG);

  FragmentExtraStorage(Target target, NodeId root_node_id)
      : target_{target}, connections_{root_node_id} {}
  MOVE_ONLY(FragmentExtraStorage);

  using FeatureTypes = typename Target::StorageType::ExtraStorageType::FeatureTypes;

  template <template <typename, typename> typename T, typename CRTP>
  using Base = typename Target::StorageType::ExtraStorageType::template Base<T, CRTP>;

  template <typename Feature>
  auto& GetFeatureStorage() {
    if constexpr (std::is_same_v<Feature, Connections>) {
      return connections_;
    } else {
      return target_.GetStorage().template GetFeatureStorage<Feature>();
    }
  }

  template <typename Feature>
  const auto& GetFeatureStorage() const {
    if constexpr (std::is_same_v<Feature, Connections>) {
      return connections_;
    } else {
      return target_.GetStorage().template GetFeatureStorage<Feature>();
    }
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
  Connections connections_;
};

template <typename Target>
struct FragmentStorage;

template <typename Target>
struct LongNameOf<FragmentStorage<Target>> {
  using type = DAGStorage<
      FragmentStorage<Target>, FragmentElementsContainer<Target, Component::Node>,
      FragmentElementsContainer<Target, Component::Edge>, FragmentExtraStorage<Target>>;
};

template <typename Target>
struct FragmentStorage : LongNameOf<FragmentStorage<Target>>::type {
  MOVE_ONLY(FragmentStorage);

  using LongNameType = typename LongNameOf<FragmentStorage>::type;
  using LongNameType::LongNameType;

  static FragmentStorage FromView(const Target& target, std::vector<NodeId>&& nodes,
                                  std::vector<EdgeId>&& edges, NodeId root_node_id) {
    static_assert(Target::role == Role::View);

    FragmentStorage result{
        FragmentElementsContainer<Target, Component::Node>{target, std::move(nodes)},
        FragmentElementsContainer<Target, Component::Edge>{target, std::move(edges)},
        FragmentExtraStorage<Target>{target, root_node_id}};

    auto view = result.View();
    view.BuildConnections();
    Assert(view.GetRoot().GetId() == root_node_id);
    return result;
  }
};

template <typename Target, typename = std::enable_if_t<Target::role == Role::View>>
FragmentStorage<Target> AddFragmentStorage(const Target& target,
                                           std::vector<NodeId>&& nodes,
                                           std::vector<EdgeId>&& edges,
                                           NodeId root_node_id) {
  return FragmentStorage<Target>::FromView(target, std::move(nodes), std::move(edges),
                                           root_node_id);
}
