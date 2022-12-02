#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

namespace Extend {

template <typename... Fs>
struct Nodes {
  using Storage = std::vector<std::tuple<Fs...>>;
  using ExtraStorage = std::tuple<ExtraFeatureStorage<Fs>...>;
  template <typename CRTP>
  struct ConstView : FeatureConstView<Fs, CRTP>... {};
  template <typename CRTP>
  struct MutableView : ConstView<CRTP>, FeatureMutableView<Fs, CRTP>... {
    using FeatureMutableView<Fs, CRTP>::operator=...;
  };
};

template <typename... Fs>
struct Edges {
  using Storage = std::vector<std::tuple<Fs...>>;
  using ExtraStorage = std::tuple<ExtraFeatureStorage<Fs>...>;
  template <typename CRTP>
  struct ConstView : FeatureConstView<Fs, CRTP>... {};
  template <typename CRTP>
  struct MutableView : ConstView<CRTP>, FeatureMutableView<Fs, CRTP>... {
    using FeatureMutableView<Fs, CRTP>::operator=...;
  };
};

template <typename... Fs>
struct DAG {
  using Storage = std::tuple<Fs...>;
  template <typename CRTP>
  struct ConstView : FeatureConstView<Fs, CRTP>... {};
  template <typename CRTP>
  struct MutableView : ConstView<CRTP>, FeatureMutableView<Fs, CRTP>... {};
};

template <typename...>
struct Empty {
  using Storage = std::tuple<>;
  template <typename, typename>
  struct ConstView {};
  template <typename, typename>
  struct MutableView {};
};
}  // namespace Extend

template <typename DV, typename Arg0 = Extend::Empty<>, typename Arg1 = Extend::Empty<>,
          typename Arg2 = Extend::Empty<>>
struct ExtendDAGStorage {
 public:
  using OnNodes = select_argument_t<Extend::Nodes, Arg0, Arg1, Arg2>;
  using OnEdges = select_argument_t<Extend::Edges, Arg0, Arg1, Arg2>;
  using OnDAG = select_argument_t<Extend::DAG, Arg0, Arg1, Arg2>;

  template <typename Id, typename CRTP>
  struct ConstElementViewBase
      : DV::StorageType::template ConstElementViewBase<Id, CRTP>,
        OnNodes::template ConstView<CRTP> {};
  template <typename Id, typename CRTP>
  struct MutableElementViewBase;

  template <typename CRTP>
  struct MutableElementViewBase<NodeId, CRTP>
      : DV::StorageType::template MutableElementViewBase<NodeId, CRTP>,
        OnNodes::template MutableView<CRTP> {
    using OnNodes::template MutableView<CRTP>::operator=;
  };

  template <typename CRTP>
  struct MutableElementViewBase<EdgeId, CRTP>
      : DV::StorageType::template MutableElementViewBase<EdgeId, CRTP>,
        OnEdges::template MutableView<CRTP> {
    using OnEdges::template MutableView<CRTP>::operator=;
  };

  template <typename CRTP>
  struct ConstDAGViewBase : DV::StorageType::template ConstDAGViewBase<CRTP>,
                            OnDAG::template ConstView<CRTP> {};
  template <typename CRTP>
  struct MutableDAGViewBase : DV::StorageType::template MutableDAGViewBase<CRTP>,
                              OnDAG::template MutableView<CRTP> {};

  MOVE_ONLY(ExtendDAGStorage);

  explicit ExtendDAGStorage(DV dv, Arg0 = Extend::Empty<>{}, Arg1 = Extend::Empty<>{},
                            Arg2 = Extend::Empty<>{});

  auto View();

  NodeId AppendNode();

  EdgeId AppendEdge();

  void AddNode(NodeId id);

  void AddEdge(EdgeId id);

  auto GetNodes() const;

  auto GetEdges() const;

  void InitializeNodes(size_t size) const;

  template <typename F>
  auto& GetFeatureStorage();

  template <typename F>
  const auto& GetFeatureStorage() const;

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

 private:
  std::decay_t<DV> target_dag_view_;
  typename OnNodes::Storage additional_node_features_storage_;
  typename OnEdges::Storage additional_edge_features_storage_;
  typename OnDAG::Storage additional_dag_features_storage_;
  typename OnNodes::ExtraStorage additional_node_extra_features_storage_;
  typename OnEdges::ExtraStorage additional_edge_extra_features_storage_;
};