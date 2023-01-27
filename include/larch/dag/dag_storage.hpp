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
  using FeatureTypes = std::tuple<Features...>;

  template <typename Id, typename CRTP>
  using ConstElementViewBase =
      std::conditional_t<std::is_same_v<Id, NodeId>,
                         typename NodesContainerT::template ConstElementViewBase<CRTP>,
                         typename EdgesContainerT::template ConstElementViewBase<CRTP>>;
  template <typename Id, typename CRTP>
  using MutableElementViewBase = std::conditional_t<
      std::is_same_v<Id, NodeId>,
      typename NodesContainerT::template MutableElementViewBase<CRTP>,
      typename EdgesContainerT::template MutableElementViewBase<CRTP>>;

  template <typename Id, typename Feature>
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
        EdgesContainerT::template ExtraMutableElementViewBase<CRTP> {};

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

  template <typename Feature>
  auto& GetFeatureStorage(NodeId id);

  template <typename Feature>
  const auto& GetFeatureStorage(NodeId id) const;

  template <typename Feature>
  auto& GetFeatureStorage(EdgeId id);

  template <typename Feature>
  const auto& GetFeatureStorage(EdgeId id) const;

  template <typename Id, typename Feature>
  auto& GetFeatureExtraStorage();

  template <typename Id, typename Feature>
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