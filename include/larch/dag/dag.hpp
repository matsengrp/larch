#pragma once

#include <vector>
#include <tuple>
#include <type_traits>

#include <tbb/concurrent_unordered_set.h>

#include "larch/common.hpp"

template <typename T>
using ConcurrentUnorderedSet =
    tbb::concurrent_unordered_set<T, std::hash<T>, std::equal_to<T>>;

//////////////////////////////////////////////////////////////////////////////////////

template <typename Feature, typename CRTP, typename Tag = Feature>
struct FeatureConstView;
template <typename Feature, typename CRTP, typename Tag = Feature>
struct FeatureMutableView;

template <typename Feature>
struct ExtraFeatureStorage {};

template <typename CRTP, typename Feature, typename Tag>
auto& GetFeatureStorage(const FeatureMutableView<Feature, CRTP, Tag>* feature) {
  return static_cast<const CRTP&>(*feature).template GetFeatureStorage<Tag>();
}

template <typename CRTP, typename Feature, typename Tag>
const auto& GetFeatureStorage(const FeatureConstView<Feature, CRTP, Tag>* feature) {
  return static_cast<const CRTP&>(*feature).template GetFeatureStorage<Tag>();
}

//////////////////////////////////////////////////////////////////////////////////////

template <typename... Fs>
struct ElementStorage {
  std::tuple<Fs...> features_storage_;

  ElementStorage() = default;
  MOVE_ONLY(ElementStorage);

  template <typename Feature>
  constexpr static const bool contains_element_feature =
      tuple_contains_v<std::tuple<Fs...>, Feature>;

  using ExtraStorage = std::tuple<ExtraFeatureStorage<Fs>...>;

  template <typename CRTP>
  struct ConstElementViewBase : FeatureConstView<Fs, CRTP>... {};
  template <typename CRTP>
  struct MutableElementViewBase : ConstElementViewBase<CRTP>,
                                  FeatureMutableView<Fs, CRTP>... {
    using FeatureMutableView<Fs, CRTP>::operator=...;
  };

  template <typename F>
  auto& GetFeatureStorage() {
    return std::get<F>(features_storage_);
  }

  template <typename F>
  const auto& GetFeatureStorage() const {
    return std::get<F>(features_storage_);
  }
};

//////////////////////////////////////////////////////////////////////////////////////

template <typename Id, typename ES, typename... Fs>
struct ElementsContainer {
  std::vector<ES> elements_storage_;
  std::vector<std::tuple<Fs...>> features_storage_;
  std::tuple<ExtraFeatureStorage<Fs>...> extra_features_storage_;
  typename ES::ExtraStorage elements_extra_features_storage_;

  ElementsContainer() = default;
  MOVE_ONLY(ElementsContainer);

  template <typename CRTP>
  struct ConstElementViewBase : ES::template ConstElementViewBase<CRTP>,
                                FeatureConstView<Fs, CRTP>... {};
  template <typename CRTP>
  struct MutableElementViewBase : ES::template MutableElementViewBase<CRTP>,
                                  FeatureMutableView<Fs, CRTP>... {
    using ES::template MutableElementViewBase<CRTP>::operator=;
    using FeatureMutableView<Fs, CRTP>::operator=...;
  };

  template <typename Feature>
  constexpr static const bool contains_element_feature =
      tuple_contains_v<std::tuple<Fs...>, Feature> or
      ES::template contains_element_feature<Feature>;

  template <typename F>
  auto& GetFeatureStorage(Id id) {
    if constexpr (tuple_contains_v<std::tuple<Fs...>, F>) {
      return std::get<F>(features_storage_.at(id.value));
    } else {
      return elements_storage_.at(id.value).template GetFeatureStorage<F>();
    }
  }

  template <typename F>
  const auto& GetFeatureStorage(Id id) const {
    if constexpr (tuple_contains_v<std::tuple<Fs...>, F>) {
      return std::get<F>(features_storage_.at(id.value));
    } else {
      return elements_storage_.at(id.value).template GetFeatureStorage<F>();
    }
  }

  template <typename F>
  auto& GetFeatureExtraStorage() {
    if constexpr (tuple_contains_v<decltype(extra_features_storage_), F>) {
      return std::get<ExtraFeatureStorage<F>>(extra_features_storage_);
    } else {
      return std::get<ExtraFeatureStorage<F>>(elements_extra_features_storage_);
    }
  }

  template <typename F>
  const auto& GetFeatureExtraStorage() const {
    if constexpr (tuple_contains_v<decltype(extra_features_storage_), F>) {
      return std::get<ExtraFeatureStorage<F>>(extra_features_storage_);
    } else {
      return std::get<ExtraFeatureStorage<F>>(elements_extra_features_storage_);
    }
  }

  size_t GetCount() const { return elements_storage_.size(); }

  Id Append() {
    Id result{GetCount()};
    elements_storage_.push_back({});
    features_storage_.push_back({});
    return result;
  }

  void Add(Id id) {
    std::ignore = GetOrInsert(elements_storage_, id);
    std::ignore = GetOrInsert(features_storage_, id);
  }

  void Initialize(size_t size) {
    elements_storage_.resize(size);
    features_storage_.resize(size);
  }
};

//////////////////////////////////////////////////////////////////////////////////////

struct NodeId {
  size_t value = NoId;
};

inline bool operator==(NodeId lhs, NodeId rhs) { return lhs == rhs; }
inline bool operator<(NodeId lhs, NodeId rhs) { return lhs < rhs; }

template <>
struct std::hash<NodeId> {
  inline size_t operator()(NodeId id) const noexcept { return id.value; }
};

struct EdgeId {
  size_t value = NoId;
};

inline bool operator==(EdgeId lhs, EdgeId rhs) { return lhs == rhs; }
inline bool operator<(EdgeId lhs, EdgeId rhs) { return lhs < rhs; }

struct CladeIdx {
  size_t value = NoId;
};

inline bool operator==(CladeIdx lhs, CladeIdx rhs) { return lhs == rhs; }
inline bool operator<(CladeIdx lhs, CladeIdx rhs) { return lhs < rhs; }

template <typename DS>
struct DAGView;

//////////////////////////////////////////////////////////////////////////////////////

template <typename NC, typename EC, typename... Fs>
struct DAGStorage {
  NC nodes_container_;
  EC edges_container_;
  std::tuple<Fs...> features_storage_;

  DAGStorage() = default;
  MOVE_ONLY(DAGStorage);

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

  template <typename CRTP>
  struct ConstDAGViewBase : FeatureConstView<Fs, CRTP>... {};
  template <typename CRTP>
  struct MutableDAGViewBase : ConstDAGViewBase<CRTP>,
                              FeatureMutableView<Fs, CRTP>... {};

  auto View() { return DAGView<DAGStorage<NC, EC, Fs...>>{*this}; }
  auto View() const { return DAGView<DAGStorage<NC, EC, Fs...>>{*this}; }

  template <typename F>
  auto& GetFeatureStorage(NodeId id) {
    return nodes_container_.template GetFeatureStorage<F>(id);
  }

  template <typename F>
  const auto& GetFeatureStorage(NodeId id) const {
    return nodes_container_.template GetFeatureStorage<F>(id);
  }

  template <typename F>
  auto& GetFeatureStorage(EdgeId id) {
    return edges_container_.template GetFeatureStorage<F>(id);
  }

  template <typename F>
  const auto& GetFeatureStorage(EdgeId id) const {
    return edges_container_.template GetFeatureStorage<F>(id);
  }

  template <typename Id, typename Feature>
  constexpr static const bool contains_element_feature = [] {
    if constexpr (std::is_same_v<Id, NodeId>) {
      return NC::template contains_element_feature<Feature>;
    } else {
      return EC::template contains_element_feature<Feature>;
    }
  }();

  template <typename Id, typename F>
  auto& GetFeatureExtraStorage() {
    if constexpr (std::is_same_v<Id, NodeId>) {
      return nodes_container_.template GetFeatureExtraStorage<F>();
    } else {
      return edges_container_.template GetFeatureExtraStorage<F>();
    }
  }

  template <typename Id, typename F>
  const auto& GetFeatureExtraStorage() const {
    if constexpr (std::is_same_v<Id, NodeId>) {
      return nodes_container_.template GetFeatureExtraStorage<F>();
    } else {
      return edges_container_.template GetFeatureExtraStorage<F>();
    }
  }

  template <typename F>
  auto& GetFeatureStorage() {
    return std::get<F>(features_storage_);
  }

  template <typename F>
  const auto& GetFeatureStorage() const {
    return std::get<F>(features_storage_);
  }

  NodeId AppendNode() { return nodes_container_.Append(); }

  EdgeId AppendEdge() { return edges_container_.Append(); }

  void AddNode(NodeId id) { nodes_container_.Add(id); }

  void AddEdge(EdgeId id) { edges_container_.Add(id); }

  size_t GetNodesCount() const { return nodes_container_.GetCount(); }

  size_t GetEdgesCount() const { return edges_container_.GetCount(); }

  const auto& GetNodes() const { return nodes_container_.elements_storage_; }  // TODO
  const auto& GetEdges() const { return edges_container_.elements_storage_; }  // TODO

  void InitializeNodes(size_t size) { nodes_container_.Initialize(size); }
};

//////////////////////////////////////////////////////////////////////////////////////

template <typename Id, typename DV>
struct ElementView;

template <typename Id, typename DV>
struct element_view_base {
  using type = std::conditional_t<
      DV::is_mutable,
      typename DV::template MutableElementViewBase<Id, ElementView<Id, DV>>,
      typename DV::template ConstElementViewBase<Id, ElementView<Id, DV>>>;
};
template <typename Id, typename DV>
using element_view_base_t = typename element_view_base<Id, DV>::type;

template <typename Id, typename DV>
struct ElementView : element_view_base_t<Id, std::decay_t<DV>> {
  std::decay_t<DV> dag_view_;
  Id id_;

  using element_view_base_t<Id, DV>::operator=;

  ElementView(DV dag_view, Id id) : dag_view_{dag_view}, id_{id} {}

  operator Id() const { return GetId(); }

  DV GetDAG() const { return dag_view_; }
  Id GetId() const { return id_; }

  template <typename F>
  auto& GetFeatureStorage() const {
    return dag_view_.template GetFeatureStorage<F>(id_);
  }

  template <typename F>
  auto& GetFeatureExtraStorage() const {
    return dag_view_.template GetFeatureExtraStorage<Id, F>();
  }
};

//////////////////////////////////////////////////////////////////////////////////////

template <typename DS>
struct DAGView
    : std::conditional_t<std::is_const_v<DS>,
                         typename DS::template ConstDAGViewBase<DAGView<DS>>,
                         typename DS::template MutableDAGViewBase<DAGView<DS>>> {
  DS& dag_storage_;

  template <typename Id, typename Feature>
  constexpr static const bool contains_element_feature =
      DS::template contains_element_feature<Id, Feature>;

  using NodeView = ElementView<NodeId, DAGView<DS>>;
  using EdgeView = ElementView<EdgeId, DAGView<DS>>;
  using StorageType = DS;

  explicit DAGView(DS& dag_storage) : dag_storage_{dag_storage} {}

  operator DAGView<const DS>() const { return DAGView<const DS>{dag_storage_}; }

  constexpr static const bool is_mutable = not std::is_const_v<DS>;
  template <typename Id, typename CRTP>
  using ConstElementViewBase = typename DS::template ConstElementViewBase<Id, CRTP>;
  template <typename Id, typename CRTP>
  using MutableElementViewBase = typename DS::template MutableElementViewBase<Id, CRTP>;

  ElementView<NodeId, DAGView<DS>> Get(NodeId id) const { return {*this, id}; }
  ElementView<EdgeId, DAGView<DS>> Get(EdgeId id) const { return {*this, id}; }

  ElementView<NodeId, DAGView<DS>> AppendNode() const {
    NodeId result = dag_storage_.AppendNode();
    return {*this, result};
  }

  ElementView<EdgeId, DAGView<DS>> AppendEdge() const {
    EdgeId result = dag_storage_.AppendEdge();
    return {*this, result};
  }

  ElementView<NodeId, DAGView<DS>> AddNode(NodeId id) {
    dag_storage_.AddNode(id);
    return {*this, id};
  }

  ElementView<EdgeId, DAGView<DS>> AddEdge(EdgeId id, NodeId parent, NodeId child,
                                           CladeIdx clade) {  // TODO
    dag_storage_.AddEdge(id);
    auto result = Get(id);
    result.Set(parent, child, clade);
    return result;
  }

  ElementView<EdgeId, DAGView<DS>> AppendEdge(NodeId parent, NodeId child,
                                              CladeIdx clade) const {  // TODO
    auto result = AppendEdge();
    result.Set(parent, child, clade);
    return result;
  }

  template <typename F>
  auto& GetFeatureStorage() const {
    return dag_storage_.template GetFeatureStorage<F>();
  }
  template <typename F>
  auto& GetFeatureStorage(NodeId id) const {
    return dag_storage_.template GetFeatureStorage<F>(id);
  }
  template <typename F>
  auto& GetFeatureStorage(EdgeId id) const {
    return dag_storage_.template GetFeatureStorage<F>(id);
  }
  template <typename Id, typename F>
  auto& GetFeatureExtraStorage() const {
    return dag_storage_.template GetFeatureExtraStorage<Id, F>();
  }

  size_t GetNodesCount() const { return dag_storage_.GetNodesCount(); }

  size_t GetEdgesCount() const { return dag_storage_.GetEdgesCount(); }

  auto GetNodes() const {
    return dag_storage_.GetNodes() |
           ranges::views::transform([*this, idx = size_t{}](auto&) mutable {
             return ElementView<NodeId, DAGView<DS>>{*this, {idx++}};
           });
  }

  auto GetEdges() const {
    return dag_storage_.GetEdges() |
           ranges::views::transform([*this, idx = size_t{}](auto&) mutable {
             return ElementView<EdgeId, DAGView<DS>>{*this, {idx++}};
           });
  }

  void InitializeNodes(size_t size) const { dag_storage_.InitializeNodes(size); }
};

//////////////////////////////////////////////////////////////////////////////////////

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
  template <template <typename...> typename Template, size_t I, typename... Ts>
  static constexpr auto SelectArgument() {
    if constexpr (I < std::tuple_size_v<std::tuple<Ts...>>) {
      if constexpr (is_specialization_v<std::tuple_element_t<I, std::tuple<Ts...>>,
                                        Template>) {
        return std::get<I>(std::tuple<Ts...>{});
      } else {
        return SelectArgument<Template, I + 1, Ts...>();
      }
    } else {
      return Template<>{};
    }
  }

  using OnNodes = decltype(SelectArgument<Extend::Nodes, 0, Arg0, Arg1, Arg2>());
  using OnEdges = decltype(SelectArgument<Extend::Edges, 0, Arg0, Arg1, Arg2>());
  using OnDAG = decltype(SelectArgument<Extend::DAG, 0, Arg0, Arg1, Arg2>());

  DV target_dag_view_;
  typename OnNodes::Storage additional_node_features_storage_;
  typename OnEdges::Storage additional_edge_features_storage_;
  typename OnDAG::Storage additional_dag_features_storage_;
  typename OnNodes::ExtraStorage additional_node_extra_features_storage_;
  typename OnEdges::ExtraStorage additional_edge_extra_features_storage_;

  ExtendDAGStorage() = default;
  MOVE_ONLY(ExtendDAGStorage);

  explicit ExtendDAGStorage(DV dv, Arg0 = Extend::Empty<>{}, Arg1 = Extend::Empty<>{},
                            Arg2 = Extend::Empty<>{})
      : target_dag_view_{dv} {
    additional_node_features_storage_.resize(target_dag_view_.NodesCount());
    additional_edge_features_storage_.resize(target_dag_view_.EdgesCount());
  }

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

  auto View() { return DAGView<ExtendDAGStorage<DV, Arg0, Arg1, Arg2>>{*this}; }

  template <typename F>
  auto& GetFeatureStorage() {
    if constexpr (tuple_contains_v<decltype(additional_dag_features_storage_), F>) {
      return std::get<F>(additional_dag_features_storage_);
    } else {
      return target_dag_view_.template GetFeatureStorage<F>();
    }
  }

  template <typename F>
  const auto& GetFeatureStorage() const {
    if constexpr (tuple_contains_v<decltype(additional_dag_features_storage_), F>) {
      return std::get<F>(additional_dag_features_storage_);
    } else {
      return target_dag_view_.template GetFeatureStorage<F>();
    }
  }

  template <typename F>
  auto& GetFeatureStorage(NodeId id) {
    if constexpr (tuple_contains_v<
                      std::decay_t<decltype(additional_node_features_storage_.at(0))>,
                      F>) {
      return std::get<F>(additional_node_features_storage_.at(id.value));
    } else {
      return target_dag_view_.template GetFeatureStorage<F>(id);
    }
  }

  template <typename F>
  const auto& GetFeatureStorage(NodeId id) const {
    if constexpr (tuple_contains_v<
                      std::decay_t<decltype(additional_node_features_storage_.at(0))>,
                      F>) {
      return std::get<F>(additional_node_features_storage_.at(id.value));
    } else {
      return target_dag_view_.template GetFeatureStorage<F>(id);
    }
  }

  template <typename F>
  auto& GetFeatureStorage(EdgeId id) {
    if constexpr (tuple_contains_v<
                      std::decay_t<decltype(additional_edge_features_storage_.at(0))>,
                      F>) {
      return std::get<F>(additional_edge_features_storage_.at(id.value));
    } else {
      return target_dag_view_.template GetFeatureStorage<F>(id);
    }
  }

  template <typename F>
  const auto& GetFeatureStorage(EdgeId id) const {
    if constexpr (tuple_contains_v<
                      std::decay_t<decltype(additional_edge_features_storage_.at(0))>,
                      F>) {
      return std::get<F>(additional_edge_features_storage_.at(id.value));
    } else {
      return target_dag_view_.template GetFeatureStorage<F>(id);
    }
  }

  template <typename Id, typename F>
  auto& GetFeatureExtraStorage() {
    if constexpr (DV::template contains_element_feature<Id, F>) {
      return target_dag_view_.template GetFeatureExtraStorage<Id, F>();
    } else {
      if constexpr (std::is_same_v<Id, NodeId>) {
        return std::get<ExtraFeatureStorage<F>>(
            additional_node_extra_features_storage_);
      } else {
        return std::get<ExtraFeatureStorage<F>>(
            additional_edge_extra_features_storage_);
      }
    }
  }

  template <typename Id, typename F>
  const auto& GetFeatureExtraStorage() const {
    if constexpr (DV::template contains_element_feature<Id, F>) {
      return target_dag_view_.template GetFeatureExtraStorage<Id, F>();
    } else {
      if constexpr (std::is_same_v<Id, NodeId>) {
        return std::get<ExtraFeatureStorage<F>>(
            additional_node_extra_features_storage_);
      } else {
        return std::get<ExtraFeatureStorage<F>>(
            additional_edge_extra_features_storage_);
      }
    }
  }

  NodeId AppendNode() {
    additional_node_features_storage_.push_back({});
    return target_dag_view_.AppendNode().GetId();
  }

  EdgeId AppendEdge() {
    additional_edge_features_storage_.push_back({});
    return target_dag_view_.AppendEdge().GetId();
  }

  void AddNode(NodeId id) {
    GetOrInsert(additional_node_features_storage_, id);
    target_dag_view_.AddNode(id);
  }

  void AddEdge(EdgeId id) {
    GetOrInsert(additional_edge_features_storage_, id);
    target_dag_view_.AddEdge(id);
  }

  auto GetNodes() const {
    return target_dag_view_.GetNodes() |
           ranges::views::transform([*this, idx = size_t{}](auto&) mutable {
             return ElementView<NodeId, DV>{*this, {idx++}};
           });
  }

  auto GetEdges() const {
    return target_dag_view_.GetEdges() |
           ranges::views::transform([*this, idx = size_t{}](auto&) mutable {
             return ElementView<EdgeId, DV>{*this, {idx++}};
           });
  }

  void InitializeNodes(size_t size) const {
    target_dag_view_.InitializeNodes(size);
    additional_node_features_storage_.resize(size);
  }
};

//////////////////////////////////////////////////////////////////////////////////////

template <typename Feature>
struct Deduplicate {
  const Feature* feature_ = nullptr;
};

template <typename Feature>
struct ExtraFeatureStorage<Deduplicate<Feature>> {
  ConcurrentUnorderedSet<Feature> deduplicated_;
};

template <typename CRTP, typename Feature>
auto& GetFeatureStorage(
    const FeatureConstView<Feature, CRTP, Deduplicate<Feature>>* feature) {
  const Feature* result = static_cast<CRTP&>(*feature)
                              .template GetFeatureStorage<Deduplicate<Feature>>()
                              .feature_;
  Assert(result);
  return *result;
}

template <typename Feature, typename CRTP>
struct FeatureConstView<Deduplicate<Feature>, CRTP>
    : FeatureConstView<Feature, CRTP, Deduplicate<Feature>> {};

template <typename Feature, typename CRTP>
struct FeatureMutableView<Deduplicate<Feature>, CRTP> {
  auto& operator=(Feature&& feature) {
    auto& deduplicated = static_cast<CRTP&>(*this)
                             .template GetFeatureExtraStorage<Deduplicate<Feature>>()
                             .deduplicated_;
    const Feature* result =
        std::addressof(*deduplicated.insert(std::forward<Feature>(feature)).first);
    static_cast<CRTP&>(*this)
        .template GetFeatureStorage<Deduplicate<Feature>>()
        .feature_ = result;
    return *this;
  }
};

//////////////////////////////////////////////////////////////////////////////////////

namespace Transform {

inline auto GetParent() {
  return ranges::views::transform([](auto&& i) { return i.GetParent(); });
}
inline auto GetChild() {
  return ranges::views::transform([](auto&& i) { return i.GetChild(); });
}
inline auto GetId() {
  return ranges::views::transform([](auto&& i) { return i.GetId(); });
}
template <typename DAG>
inline auto ToNodes(DAG dag) {
  return ranges::views::transform([dag](auto&& i) {
    return typename DAG::NodeView{dag, i};
  });
}
template <typename DAG>
inline auto ToEdges(DAG dag) {
  return ranges::views::transform([dag](auto&& i) {
    return typename DAG::EdgeView{dag, i};
  });
}

}  // namespace Transform

//////////////////////////////////////////////////////////////////////////////////////

struct Neighbors {
  std::vector<EdgeId> parents_;
  std::vector<std::vector<EdgeId>> clades_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<Neighbors, CRTP, Tag> {
  auto GetParents() const {
    auto dag = static_cast<const CRTP&>(*this).GetDAG();
    return GetFeatureStorage(this).parents_ | Transform::ToEdges(dag);
  }
  auto GetClades() const {
    auto dag = static_cast<const CRTP&>(*this).GetDAG();
    return GetFeatureStorage(this).clades_ |
           ranges::views::transform([*this, dag](const std::vector<EdgeId>& clade) {
             return clade | Transform::ToEdges(dag);
           });
  }
  auto GetClade(CladeIdx clade) const {
    auto dag = static_cast<const CRTP&>(*this).GetDAG();
    return GetFeatureStorage(this).clades_.at(clade.value) | Transform::ToEdges(dag);
  }
  size_t GetParentsCount() const { return GetFeatureStorage(this).parents_.size(); }
  size_t GetCladesCount() const { return GetFeatureStorage(this).clades_.size(); }
  auto GetChildren() const { return GetClades() | ranges::views::join; }
  auto GetSingleParent() const {
    Assert(GetParentsCount() == 1);
    return *GetParents().begin();
  }
  auto GetFirstChild() const {
    Assert(not IsLeaf());
    return *GetChildren().begin();
  }
  auto GetFirstClade() const { return GetClade({0}); }
  bool IsRoot() const { return GetFeatureStorage(this).parents_.empty(); }
  bool IsLeaf() const {
    return ranges::all_of(GetFeatureStorage(this).clades_,
                          [](const auto& clade) { return clade.empty(); });
  }
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<Neighbors, CRTP, Tag> {
  void ClearConnections() {
    GetFeatureStorage(this).parents_.clear();
    GetFeatureStorage(this).clades_.clear();
  }
  void AddEdge(CladeIdx clade, EdgeId id, bool this_node_is_parent) {
    if (this_node_is_parent) {
      GetOrInsert(GetFeatureStorage(this).clades_, clade).push_back(id);
    } else {
      GetFeatureStorage(this).parents_.push_back(id);
    }
  }
};

//////////////////////////////////////////////////////////////////////////////////////

struct Endpoints {
  NodeId parent_;
  NodeId child_;
  CladeIdx clade_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<Endpoints, CRTP, Tag> {
  auto GetParent() const {
    auto dag = static_cast<const CRTP&>(*this).GetDAG();
    return typename decltype(dag)::NodeView{dag, GetParentId()};
  }
  auto GetChild() const {
    auto dag = static_cast<const CRTP&>(*this).GetDAG();
    return typename decltype(dag)::NodeView{dag, GetChildId()};
  }
  CladeIdx GetClade() const { return GetFeatureStorage(this).clade_; }
  NodeId GetParentId() const { return GetFeatureStorage(this).parent_; }
  NodeId GetChildId() const { return GetFeatureStorage(this).child_; }
  std::pair<NodeId, NodeId> GetNodeIds() const { return {GetParentId(), GetChildId()}; }
  bool IsRoot() const;
  bool IsLeaf() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<Endpoints, CRTP, Tag> {
  void Set(NodeId parent, NodeId child, CladeIdx clade) {
    GetFeatureStorage(this).parent_ = parent;
    GetFeatureStorage(this).child_ = child;
    GetFeatureStorage(this).clade_ = clade;
  }
};

//////////////////////////////////////////////////////////////////////////////////////

struct Connections {
  NodeId root_ = {NoId};
  std::vector<NodeId> leafs_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<Connections, CRTP, Tag> {
  bool IsTree() const {
    auto& dag = static_cast<const CRTP&>(*this);
    return dag.GetNodesCount() == dag.GetEdgesCount() + 1;
  }
  bool HaveRoot() const { return GetFeatureStorage(this).root_.value != NoId; }
  auto GetRoot() const {
    auto& dag = static_cast<const CRTP&>(*this);
    using Node = typename CRTP::NodeView;
    return Node{dag, GetFeatureStorage(this).root_};
  }
  auto GetLeafs() const {
    auto& dag = static_cast<const CRTP&>(*this);
    return GetFeatureStorage(this).leafs_ |
           ranges::views::transform([*this, dag, idx = size_t{}](auto&) mutable {
             using Node = typename CRTP::NodeView;
             return Node{dag, {idx++}};
           });
  }
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<Connections, CRTP, Tag> {
  void BuildConnections() const {
    auto& storage = GetFeatureStorage(this);
    auto& dag = static_cast<const CRTP&>(*this);
    storage.root_ = {NoId};
    storage.leafs_ = {};
    BuildConnectionsRaw();
    for (auto node : dag.GetNodes()) {
      for (auto clade : node.GetClades()) {
        Assert(not clade.empty() && "Empty clade");
      }
      if (node.IsRoot()) {
        Assert(storage.root_.value == NoId && "Duplicate root");
        storage.root_ = node;
      }
      if (node.IsLeaf()) {
        storage.leafs_.push_back(node);
      }
    }
  }
  void BuildConnectionsRaw() const {
    auto& dag = static_cast<const CRTP&>(*this);
    for (auto node : dag.GetNodes()) {
      node.ClearConnections();
    }
    EdgeId edge_id = {0};
    for (auto edge : dag.GetEdges()) {
      Assert(edge.GetParentId().value != NoId && "Edge has no parent");
      Assert(edge.GetChildId().value != NoId && "Edge has no child");
      Assert(edge.GetClade().value != NoId && "Edge has no clade index");
      edge.GetParent().AddEdge(edge.GetClade(), edge, true);
      edge.GetChild().AddEdge(edge.GetClade(), edge, false);
      ++edge_id.value;
    }
  }
};

//////////////////////////////////////////////////////////////////////////////////////