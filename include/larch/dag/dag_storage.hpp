#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename... Features>
struct ExtraStorage {
  MOVE_ONLY(ExtraStorage);
  ExtraStorage() = default;

  constexpr static const Component component = Component::DAG;
  constexpr static const Role role = Role::Storage;

  using FeatureTypes = std::tuple<Features...>;

  template <template <typename, typename> typename T, typename CRTP>
  struct Base : T<Features, CRTP>... {};

  template <typename Feature>
  auto GetFeatureStorage() {
    return std::ref(tuple_get<Feature, FeatureEquivalent>(features_storage_));
  }

  template <typename Feature>
  auto GetFeatureStorage() const {
    return std::cref(tuple_get<Feature, FeatureEquivalent>(features_storage_));
  }

  template <typename Storage>
  auto& GetTargetStorage(Storage& storage) const {
    return storage;
  }

  template <typename Storage>
  auto& GetTargetStorage(Storage& storage) {
    return storage;
  }

 private:
  FeatureTypes features_storage_;
};

/**
 * Owning storage for an entire DAG. All views (DAG, node or edge) are
 * internally pointing to an instance of DAGStorage or classes providing
 * the same interface, like ExtendDAGStorage.
 */
template <typename ShortName, typename NodesContainerT, typename EdgesContainerT,
          typename ExtraStorageT, template <typename, typename> typename ViewBase>
struct DAGStorage {
  static_assert(ExtraStorageT::role == Role::Storage);
  static_assert(ExtraStorageT::component == Component::DAG);

  constexpr static const Component component = Component::DAG;
  constexpr static const Role role = Role::Storage;

  using Self =
      std::conditional_t<std::is_same_v<ShortName, void>, DAGStorage, ShortName>;

  using ViewType = DAGView<Self, ViewBase>;
  using ConstViewType = DAGView<const Self, ViewBase>;

  using ExtraStorageType = ExtraStorageT;
  using FeatureTypes = typename ExtraStorageT::FeatureTypes;
  template <Component C>
  using Container =
      std::conditional_t<C == Component::Node, NodesContainerT, EdgesContainerT>;
  using AllNodeFeatures = typename NodesContainerT::AllFeatureTypes;
  using AllEdgeFeatures = typename EdgesContainerT::AllFeatureTypes;

  template <Component C, typename CRTP>
  using ConstElementViewBase =
      std::conditional_t<C == Component::Node,
                         typename NodesContainerT::template ConstElementViewBase<CRTP>,
                         typename EdgesContainerT::template ConstElementViewBase<CRTP>>;
  template <Component C, typename CRTP>
  using MutableElementViewBase = std::conditional_t<
      C == Component::Node,
      typename NodesContainerT::template MutableElementViewBase<CRTP>,
      typename EdgesContainerT::template MutableElementViewBase<CRTP>>;

  template <Component C, typename Feature>
  static const bool contains_element_feature;

  template <typename CRTP>
  struct ConstDAGViewBase : ExtraStorageT::template Base<FeatureConstView, CRTP>,
                            ExtraStorageT::template Base<ExtraFeatureConstView, CRTP>,
                            NodesContainerT::template ExtraConstElementViewBase<CRTP>,
                            EdgesContainerT::template ExtraConstElementViewBase<CRTP> {
  };

  template <typename CRTP>
  struct MutableDAGViewBase
      : ExtraStorageT::template Base<FeatureMutableView, CRTP>,
        ExtraStorageT::template Base<ExtraFeatureMutableView, CRTP>,
        NodesContainerT::template ExtraMutableElementViewBase<CRTP>,
        EdgesContainerT::template ExtraMutableElementViewBase<CRTP> {
    template <typename F>
    constexpr auto& AsFeature() const noexcept {
      return static_cast<const ExtraFeatureMutableView<F, CRTP>&>(*this);
    }
  };

  DAGStorage() = default;
  DAGStorage(NodesContainerT&& nodes_container, EdgesContainerT&& edges_container,
             ExtraStorageT&& features_storage);
  MOVE_ONLY(DAGStorage);

  template <template <typename, typename> typename Base = ViewBase>
  DAGView<Self, Base> View();
  template <template <typename, typename> typename Base = ViewBase>
  DAGView<const Self, Base> View() const;

  NodeId AppendNode();
  EdgeId AppendEdge();

  void AddNode(NodeId id);
  void AddEdge(EdgeId id);

  template <typename VT>
  auto GetNodes() const {
    return nodes_container_.template All<VT>();
  }
  template <typename VT>
  auto GetEdges() const {
    return edges_container_.template All<VT>();
  }

  template <typename VT>
  size_t GetNodesCount() const;
  template <typename VT>
  size_t GetEdgesCount() const;

  template <Component C, typename VT>
  Id<C> GetNextAvailableId() const {
    if constexpr (C == Component::Node) {
      return nodes_container_.template GetNextAvailableId<VT>();
    } else {
      return edges_container_.template GetNextAvailableId<VT>();
    }
  }

  template <typename VT>
  bool ContainsId(NodeId id) const {
    return nodes_container_.template ContainsId<VT>(id);
  }

  template <typename VT>
  bool ContainsId(EdgeId id) const {
    return edges_container_.template ContainsId<VT>(id);
  }

  void InitializeNodes(size_t size);
  void InitializeEdges(size_t size);

  void ClearNodes();
  void ClearEdges();

  template <typename Feature>
  auto GetFeatureStorage(NodeId id);

  template <typename Feature>
  auto GetFeatureStorage(NodeId id) const;

  template <typename Feature>
  auto GetFeatureStorage(EdgeId id);

  template <typename Feature>
  auto GetFeatureStorage(EdgeId id) const;

  template <Component C, typename Feature>
  auto GetFeatureExtraStorage();

  template <Component C, typename Feature>
  auto GetFeatureExtraStorage() const;

  template <typename Feature>
  auto GetFeatureStorage();

  template <typename Feature>
  auto GetFeatureStorage() const;

  auto& GetTargetStorage() { return features_storage_.GetTargetStorage(*this); }
  auto& GetTargetStorage() const { return features_storage_.GetTargetStorage(*this); }

  NodesContainerT& GetNodesContainer() { return nodes_container_; }
  EdgesContainerT& GetEdgesContainer() { return edges_container_; }

 private:
  NodesContainerT nodes_container_;
  EdgesContainerT edges_container_;
  ExtraStorageT features_storage_;
};
