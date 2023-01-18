#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * Owning storage for an entire DAG. All views (DAG, node or edge) are
 * internally pointing to an instance of DAGStorage or classes providing
 * the same interface, like ExtendDAGStorage.
 */
template <typename NC, typename EC, typename... Fs>
struct DAGStorage {
 public:
  using FeatureTypes = std::tuple<Fs...>;

  template <typename Id, typename CRTP>
  using ConstElementViewBase =
      std::conditional_t<std::is_same_v<Id, NodeId>,
                         typename NC::template ConstElementViewBase<CRTP>,
                         typename EC::template ConstElementViewBase<CRTP>>;
  template <typename Id, typename CRTP>
  using MutableElementViewBase =
      std::conditional_t<std::is_same_v<Id, NodeId>,
                         typename NC::template MutableElementViewBase<CRTP>,
                         typename EC::template MutableElementViewBase<CRTP>>;

  template <typename Id, typename Feature>
  static const bool contains_element_feature;

  template <typename CRTP>
  struct ConstDAGViewBase : FeatureConstView<Fs, CRTP>...,
                            ExtraFeatureConstView<Fs, CRTP>...,
                            NC::template ExtraConstElementViewBase<CRTP>,
                            EC::template ExtraConstElementViewBase<CRTP> {};
  template <typename CRTP>
  struct MutableDAGViewBase : ConstDAGViewBase<CRTP>,
                              FeatureMutableView<Fs, CRTP>...,
                              ExtraFeatureMutableView<Fs, CRTP>...,
                              NC::template ExtraMutableElementViewBase<CRTP>,
                              EC::template ExtraMutableElementViewBase<CRTP> {};

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

  template <typename F>
  auto& GetFeatureStorage(NodeId id);

  template <typename F>
  const auto& GetFeatureStorage(NodeId id) const;

  template <typename F>
  auto& GetFeatureStorage(EdgeId id);

  template <typename F>
  const auto& GetFeatureStorage(EdgeId id) const;

  template <typename Id, typename F>
  auto& GetFeatureExtraStorage();

  template <typename Id, typename F>
  const auto& GetFeatureExtraStorage() const;

  template <typename F>
  auto& GetFeatureStorage();

  template <typename F>
  const auto& GetFeatureStorage() const;

 private:
  NC nodes_container_;
  EC edges_container_;
  std::tuple<Fs...> features_storage_;
};