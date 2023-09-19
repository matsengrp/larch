#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * Owning storage for an entire DAG. All views (DAG, node or edge) are
 * internally pointing to an instance of DAGStorage or classes providing
 * the same interface, like ExtendDAGStorage.
 */
template <typename NodesContainerT, typename EdgesContainerT, typename... Features>
struct DAGStorage {
 public:
  constexpr static const Component component = Component::DAG;
  constexpr static const Role role = Role::Storage;

  using FeatureTypes = std::tuple<Features...>;
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
  struct ConstDAGViewBase : FeatureConstView<Features, CRTP>...,
                            ExtraFeatureConstView<Features, CRTP>...,
                            NodesContainerT::template ExtraConstElementViewBase<CRTP>,
                            EdgesContainerT::template ExtraConstElementViewBase<CRTP> {
  };
  template <typename CRTP>
  struct MutableDAGViewBase
      : ConstDAGViewBase<CRTP>,
        FeatureMutableView<Features, CRTP>...,
        ExtraFeatureMutableView<Features, CRTP>...,
        NodesContainerT::template ExtraMutableElementViewBase<CRTP>,
        EdgesContainerT::template ExtraMutableElementViewBase<CRTP> {
    template <typename F>
    constexpr auto& AsFeature() const noexcept {
      return static_cast<const ExtraFeatureMutableView<F, CRTP>&>(*this);
    }
  };

  DAGStorage() = default;
  MOVE_ONLY(DAGStorage);

  auto View();
  auto View() const;

  NodeId AppendNode();
  EdgeId AppendEdge();

  void AddNode(NodeId id);
  void AddEdge(EdgeId id);

  size_t GetNodesCount() const;
  size_t GetEdgesCount() const;

  void InitializeNodes(size_t size);
  void InitializeEdges(size_t size);

  void ClearNodes();
  void ClearEdges();

  template <typename Feature>
  auto& GetFeatureStorage(NodeId id);

  template <typename Feature>
  const auto& GetFeatureStorage(NodeId id) const;

  template <typename Feature>
  auto& GetFeatureStorage(EdgeId id);

  template <typename Feature>
  const auto& GetFeatureStorage(EdgeId id) const;

  template <Component C, typename Feature>
  auto& GetFeatureExtraStorage();

  template <Component C, typename Feature>
  const auto& GetFeatureExtraStorage() const;

  template <typename Feature>
  auto& GetFeatureStorage();

  template <typename Feature>
  const auto& GetFeatureStorage() const;

 private:
  NodesContainerT nodes_container_;
  EdgesContainerT edges_container_;
  std::tuple<Features...> features_storage_;
};
