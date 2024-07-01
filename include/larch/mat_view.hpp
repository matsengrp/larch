#pragma once

#include "larch/mat_conversion.hpp"

template <typename DAGStorageType, typename DAGViewType>
struct CondensedViewBase : DefaultViewBase<DAGStorageType, DAGViewType> {
  static constexpr inline bool is_condensed = true;
};

template <typename DAGStorageType, typename DAGViewType>
struct UncondensedViewBase : DefaultViewBase<DAGStorageType, DAGViewType> {
  static constexpr inline bool is_condensed = false;
};

template <typename T, typename = void>
struct CheckIsCondensed : std::false_type {};

template <typename T>
struct CheckIsCondensed<T, std::void_t<decltype(T::BaseType::is_condensed)>>
    : std::bool_constant<T::BaseType::is_condensed> {};

struct MATNodeStorage {
  MOVE_ONLY(MATNodeStorage);
  MATNodeStorage() = default;

  inline MATNodeStorage Copy() const { return {}; }
};

template <>
struct ExtraFeatureStorage<MATNodeStorage> {
  ExtraFeatureStorage() = default;
  MOVE_ONLY(ExtraFeatureStorage);
  MAT::Tree* mat_tree_ = nullptr;
  NodeId ua_node_id_;
  std::map<NodeId, std::vector<std::string>> condensed_nodes_;
  std::map<std::string, MAT::Node*> reversed_condensed_nodes_;
  size_t condensed_nodes_count_ = 0;
};

template <typename CRTP>
struct ExtraFeatureConstView<MATNodeStorage, CRTP> {
  const MAT::Tree& GetMAT() const {
    auto& dag = static_cast<const CRTP&>(*this);
    auto* mat = dag.template GetFeatureExtraStorage<Component::Node, MATNodeStorage>()
                    .mat_tree_;
    Assert(mat != nullptr);
    return *mat;
  }

  constexpr static bool IsCondensed() { return CheckIsCondensed<CRTP>::value; }

  auto GetUncondensed() const {
    static_assert(IsCondensed());
    auto& dag = static_cast<const CRTP&>(*this);
    return dag.GetStorage().template View<UncondensedViewBase>();
  }
};

struct MATEdgeStorage;

template <typename CRTP>
struct ExtraFeatureMutableView<MATNodeStorage, CRTP> {
  void SetMAT(MAT::Tree* mat) const {
    Assert(mat != nullptr);
    auto& dag = static_cast<const CRTP&>(*this);
    auto& node_storage =
        dag.template GetFeatureExtraStorage<Component::Node, MATNodeStorage>();
    auto& edge_storage =
        dag.template GetFeatureExtraStorage<Component::Edge, MATEdgeStorage>();

    node_storage.mat_tree_ = mat;
    edge_storage.mat_tree_ = mat;
    size_t ua_node_id = 0;
    Assert(mat->get_node(ua_node_id) == nullptr);
    node_storage.ua_node_id_ = NodeId{ua_node_id};
    edge_storage.ua_node_id_ = NodeId{ua_node_id};
    auto& cn = node_storage.condensed_nodes_;
    for (auto& [i, j] : mat->condensed_nodes) {
      auto& nodes = cn[NodeId{i}];
      for (auto& k : j) {
        nodes.push_back(k);
        node_storage.reversed_condensed_nodes_.insert_or_assign(k, mat->get_node(i));
      }
      node_storage.condensed_nodes_count_ += nodes.size();
    }
    edge_storage.condensed_nodes_ = node_storage.condensed_nodes_;
    edge_storage.reversed_condensed_nodes_ = node_storage.reversed_condensed_nodes_;
    edge_storage.condensed_nodes_count_ = node_storage.condensed_nodes_count_;
  }
};

template <typename T>
struct MATValidator {
  MATValidator(const T& storage) : storage_{storage} {}

  auto CladesRange() const { return storage_.GetClades(); }

  auto ParentsRange() const { return storage_.GetParents(); }

 private:
  const T& storage_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<MATNodeStorage, CRTP, Tag> {
  auto GetParents() const {
    return ranges::views::iota(size_t{0}, GetParentsCount()) |
           ranges::views::transform([this](size_t) { return GetSingleParent(); });
  }
  auto GetClades() const {
    return ranges::views::iota(size_t{0}, GetCladesCount()) |
           ranges::views::transform([this](size_t i) { return GetClade(CladeIdx{i}); });
  }

  auto GetClade(CladeIdx clade) const {
    auto [dag_node, mat, mat_node, is_ua] = access();
    size_t node_id =
        is_ua ? mat.root->node_id : mat_node->children.at(clade.value)->node_id;
    return ranges::views::iota(node_id, node_id + 1) |
           Transform::ToId<Component::Edge>() | Transform::ToEdges(dag_node.GetDAG());
  }

  size_t GetParentsCount() const {
    auto [dag_node, mat, mat_node, is_ua] = access();
    if (is_ua) {
      return 0;
    } else {
      return 1;
    }
  }

  size_t GetCladesCount() const {
    auto [dag_node, mat, mat_node, is_ua] = access();
    if (is_ua) {
      return 1;
    }
    if constexpr (CheckIsCondensed<decltype(dag_node.GetDAG())>::value) {
    }
    return mat_node->children.size();
  }

  auto GetChildren() const {
    auto [dag_node, mat, mat_node, is_ua] = access();
    auto dag = dag_node.GetDAG();

    using Edge = typename decltype(dag)::EdgeView;

    auto child_edge = [dag, dag_node, is_ua](const MAT::Node* i) {
      if (is_ua) {
        return Edge{dag, EdgeId{dag_node.GetFirstChild().GetId().value}};
      } else {
        Assert(i != nullptr);
        return Edge{dag, EdgeId{i->node_id}};
      }
    };
    if (is_ua) {
      return empty_node | ranges::views::transform(child_edge);
    } else {
      return mat_node->children | ranges::views::transform(child_edge);
    }
  }

  auto GetSingleParent() const {
    auto [dag_node, mat, mat_node, is_ua] = access();
    Assert(not is_ua);
    EdgeId parent{mat_node->parent == nullptr ? mat.root->node_id : mat_node->node_id};
    auto dag = dag_node.GetDAG();
    return typename decltype(dag)::EdgeView{dag, parent};
  }

  auto GetFirstParent() const {
    Assert(GetParentsCount() == 1);
    return GetSingleParent();
  }

  auto GetFirstChild() const {
    Assert(not IsLeaf());
    return (*GetFirstClade().begin()).GetChild();
  }

  auto GetFirstClade() const {
    Assert(not IsLeaf());
    return GetClade({0});
  }

  bool IsUA() const {
    auto [dag_node, mat, mat_node, is_ua] = access();
    return is_ua;
  }

  bool IsTreeRoot() const;

  bool IsLeaf() const { return GetCladesCount() == 0; }

  auto GetLeafsBelow() const;

  void Validate(bool recursive = false, bool allow_dag = false) const {
    auto node = static_cast<const CRTP&>(*this).Const();
    ValidateImpl(node, MATValidator{*this}, recursive, allow_dag);
  }

  auto GetParentNodes() const;
  auto GetChildNodes() const;

  bool ContainsParent(NodeId node) const {
    auto [dag_node, mat, mat_node, is_ua] = access();
    if (is_ua) {
      return false;
    }
    Assert(mat_node != nullptr);
    if (mat_node->parent == nullptr) {
      return node == GetUA();
    }
    return mat_node->parent->node_id == node.value;
  }

  bool ContainsChild(NodeId node) const {
    auto [dag_node, mat, mat_node, is_ua] = access();
    if (is_ua) {
      return node.value == mat.root->node_id;
    }
    for (auto* i : mat_node->children) {
      if (i->node_id == node.value) {
        return true;
      }
    }
    return false;
  }

  std::string ParentsToString() const;
  std::string ChildrenToString() const;

 private:
  static inline std::vector<MAT::Node*> empty_node{nullptr};

  NodeId GetUA() const {
    auto dag_node = static_cast<const CRTP&>(*this);
    auto dag = dag_node.GetDAG();
    return dag.template GetFeatureExtraStorage<Component::Node, MATNodeStorage>()
        .ua_node_id_;
  }

  auto access() const {
    auto dag_node = static_cast<const CRTP&>(*this);
    NodeId id = dag_node.GetId();
    auto dag = dag_node.GetDAG();
    auto& mat = dag.GetMAT();
    auto* mat_node = mat.get_node(id.value);
    bool is_ua = dag.template GetFeatureExtraStorage<Component::Node, MATNodeStorage>()
                     .ua_node_id_ == id;
    if (not is_ua) {
      Assert(mat_node != nullptr);
    }
    return std::make_tuple(dag_node, std::ref(mat), mat_node, is_ua);
  }
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<MATNodeStorage, CRTP, Tag> {};

struct MATEdgeStorage {
  MOVE_ONLY(MATEdgeStorage);
  MATEdgeStorage() = default;

  inline MATEdgeStorage Copy() const { return {}; }

  mutable EdgeMutations mutations_;
};

template <>
struct ExtraFeatureStorage<MATEdgeStorage> {
  ExtraFeatureStorage() = default;
  MOVE_ONLY(ExtraFeatureStorage);
  MAT::Tree* mat_tree_ = nullptr;
  NodeId ua_node_id_;
  std::map<NodeId, std::vector<std::string>> condensed_nodes_;
  std::map<std::string, MAT::Node*> reversed_condensed_nodes_;
  size_t condensed_nodes_count_ = 0;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<MATEdgeStorage, CRTP, Tag> {
  auto GetParent() const {
    auto [dag_edge, mat, mat_node, is_ua] = access();
    if (is_ua) {
      return dag_edge.GetDAG().GetRoot();
    }
    auto& storage = dag_edge.template GetFeatureExtraStorage<MATEdgeStorage>();
    if (mat_node == nullptr) {
      // SEGFAULT ERROR HERE: there is a segfault with dag_edge.GetChild() call
      auto cn_id_iter = storage.condensed_nodes_.find(dag_edge.GetChild().GetId());
      Assert(cn_id_iter != storage.condensed_nodes_.end());
      auto condensed_mat_node_str = storage.condensed_nodes_.at(dag_edge.GetChild().GetId());
      auto condensed_mat_node = storage.reversed_condensed_nodes_.at(condensed_mat_node_str.front());
      Assert(condensed_mat_node->parent != nullptr);
      return dag_edge.GetDAG().Get(NodeId{condensed_mat_node->parent->node_id});
      // TODO return parent
    }
    Assert(mat_node->parent != nullptr);
    return dag_edge.GetDAG().Get(NodeId{mat_node->parent->node_id});
  }

  auto GetChild() const {
    auto [dag_edge, mat, mat_node, is_ua] = access();
    if (is_ua) {
      return dag_edge.GetDAG().Get(NodeId{mat.root->node_id});
    }
    return dag_edge.GetDAG().Get(NodeId{mat_node->node_id});
  }

  CladeIdx GetClade() const {
    auto [dag_edge, mat, mat_node, is_ua] = access();
    CladeIdx result{0};
    if (is_ua) {
      return result;
    }
    Assert(mat_node->parent != nullptr);
    for (auto* i : mat_node->parent->children) {
      if (i == mat_node) {
        return result;
      }
      ++result.value;
    }
    Fail("Clade not found");
  }

  auto GetParentId() const {
    auto [dag_edge, mat, mat_node, is_ua] = access();
    if (is_ua) {
      return dag_edge.GetDAG().GetRoot().GetId();
    }
    Assert(mat_node->parent != nullptr);
    return NodeId{mat_node->parent->node_id};
  }

  auto GetChildId() const {
    auto [dag_edge, mat, mat_node, is_ua] = access();
    return NodeId{mat_node->node_id};
  }

  std::pair<NodeId, NodeId> GetNodeIds() const {
    auto [dag_edge, mat, mat_node, is_ua] = access();
    if (is_ua) {
      return {dag_edge.GetDAG().GetRoot(), {mat.root->node_id}};
    }
    Assert(mat_node->parent != nullptr);
    return {{mat_node->parent->node_id}, {mat_node->node_id}};
  }

  bool IsUA() const {
    auto [dag_edge, mat, mat_node, is_ua] = access();
    return is_ua;
  }

  bool IsTreeRoot() const { return GetParent().IsTreeRoot(); }

  bool IsLeaf() const { return GetChild().IsLeaf(); }

  const EdgeMutations& GetEdgeMutations() const {
    auto [dag_edge, mat, mat_node, is_ua] = access();
    auto& storage = dag_edge.template GetFeatureStorage<MATEdgeStorage>();
    if (storage.mutations_.empty()) {
      storage.mutations_ = EdgeMutations{
          mat_node->mutations |
          ranges::views::transform(
              [](const MAT::Mutation& mut)
                  -> std::pair<MutationPosition, std::pair<char, char>> {
                static const std::array<char, 4> decode = {'A', 'C', 'G', 'T'};
                return {{static_cast<size_t>(mut.get_position())},
                        {decode.at(one_hot_to_two_bit(mut.get_par_one_hot())),
                         decode.at(one_hot_to_two_bit(mut.get_mut_one_hot()))}};
              })};
    }
    return storage.mutations_;
  }

 private:
  auto access() const {
    auto dag_edge = static_cast<const CRTP&>(*this);
    EdgeId id = dag_edge.GetId();
    auto dag = dag_edge.GetDAG();
    auto& mat = dag.GetMAT();
    Assert(id.value < mat.get_size_upper());
    MAT::Node* mat_node = mat.get_node(id.value);
    bool is_ua = (mat.root == mat_node);
    return std::make_tuple(dag_edge, std::ref(mat), mat_node, is_ua);
  }
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<MATEdgeStorage, CRTP, Tag> {};

class MATStorageImpl;

template <Component C, bool Condensed>
struct MATElementsContainerBase {
  using ElementStorageT =
      std::conditional_t<C == Component::Node,
                         std::conditional_t<Condensed, MATNodeStorage, MATNodeStorage>,
                         std::conditional_t<Condensed, MATEdgeStorage, MATEdgeStorage>>;
  using FeatureTypes = std::tuple<ElementStorageT>;
  using AllFeatureTypes = FeatureTypes;

  static constexpr IdContinuity id_continuity = IdContinuity::Dense;

  template <typename Feature>
  static const bool contains_element_feature =
      tuple_contains_v<AllFeatureTypes, Feature>;

  template <typename CRTP>
  struct ConstElementViewBase : FeatureConstView<ElementStorageT, CRTP> {};
  template <typename CRTP>
  struct MutableElementViewBase : FeatureMutableView<ElementStorageT, CRTP> {
    using FeatureMutableView<ElementStorageT, CRTP>::operator=;
  };

  template <typename CRTP>
  struct ExtraConstElementViewBase : ExtraFeatureConstView<ElementStorageT, CRTP> {};
  template <typename CRTP>
  struct ExtraMutableElementViewBase : ExtraFeatureMutableView<ElementStorageT, CRTP> {
  };

  MATElementsContainerBase(MATStorageImpl& impl) : impl_{impl} {}

  size_t GetCount() const {
    size_t count = GetMAT().get_size_upper();
    if constexpr (C == Component::Edge) {
      Assert(count > 0);
      count -= 1;
    }
    if constexpr (not Condensed) {
      count += extra_storage_.condensed_nodes_count_;
    }
    return count;
  }

  Id<C> GetNextAvailableId() const { return {GetCount()}; }

  template <typename Feature>
  const auto& GetFeatureStorage(Id<C> id) const {
    if constexpr (C == Component::Node) {
      // TODO
    } else {
      if (features_storage_.empty()) {
        features_storage_.resize(GetMAT().get_size_upper());
      }
      return std::get<Feature>(features_storage_.at(id));
    }
  }

  template <typename Feature>
  auto& GetFeatureStorage(Id<C> id) {
    if constexpr (C == Component::Node) {
      // TODO
    } else {
      if (features_storage_.empty()) {
        features_storage_.resize(GetMAT().get_size_upper());
      }
      return std::get<Feature>(features_storage_.at(id));
    }
  }

  template <typename Feature>
  auto& GetFeatureExtraStorage() {
    static_assert(std::is_same_v<Feature, ElementStorageT>);
    return extra_storage_;
  }

  template <typename Feature>
  const auto& GetFeatureExtraStorage() const {
    static_assert(std::is_same_v<Feature, ElementStorageT>);
    return extra_storage_;
  }

  auto All() const {
    if constexpr (C == Component::Node) {
      return ranges::views::iota(size_t{0}, GetCount()) |
             ranges::views::filter([this](size_t i) {
               return GetMAT().get_node(i) != nullptr or
                      i == extra_storage_.ua_node_id_.value;
             }) |
             ranges::views::transform([](size_t i) -> NodeId { return {i}; });
    } else {
      return ranges::views::iota(size_t{0}, GetCount()) |
             ranges::views::filter([this](size_t i) {
               return GetMAT().get_node(i) != nullptr or
                      i != extra_storage_.ua_node_id_.value;
             }) |
             ranges::views::transform([](size_t i) -> EdgeId { return {i}; });
    }
  }

 private:
  MAT::Tree& GetMAT() {
    Assert(extra_storage_.mat_tree_ != nullptr);
    return *extra_storage_.mat_tree_;
  }

  const MAT::Tree& GetMAT() const {
    Assert(extra_storage_.mat_tree_ != nullptr);
    return *extra_storage_.mat_tree_;
  }

  MATStorageImpl& impl_;
  mutable IdContainer<Id<C>, AllFeatureTypes, id_continuity> features_storage_ = {};
  ExtraFeatureStorage<ElementStorageT> extra_storage_ = {};
};

class UncondensedNodesContainer
    : public MATElementsContainerBase<Component::Node, false> {
 public:
  UncondensedNodesContainer(MATStorageImpl& impl)
      : MATElementsContainerBase<Component::Node, false>{impl} {}
};

class UncondensedEdgesContainer
    : public MATElementsContainerBase<Component::Edge, false> {
 public:
  UncondensedEdgesContainer(MATStorageImpl& impl)
      : MATElementsContainerBase<Component::Edge, false>{impl} {}
};

class CondensedNodesContainer : public MATElementsContainerBase<Component::Node, true> {
 public:
  CondensedNodesContainer(MATStorageImpl& impl)
      : MATElementsContainerBase<Component::Node, true>{impl} {}
};

class CondensedEdgesContainer : public MATElementsContainerBase<Component::Edge, true> {
 public:
  CondensedEdgesContainer(MATStorageImpl& impl)
      : MATElementsContainerBase<Component::Edge, true>{impl} {}
};
;

class MATStorageImpl {
 private:
  DAGStorage<void, UncondensedNodesContainer, UncondensedEdgesContainer,
             ExtraStorage<Connections>, UncondensedViewBase>
      uncondensed_;
  DAGStorage<void, CondensedNodesContainer, CondensedEdgesContainer,
             ExtraStorage<Connections>, CondensedViewBase>
      condensed_;

 public:
  using UncondensedT = decltype(uncondensed_);
  using CondensedT = decltype(condensed_);

  MATStorageImpl()
      : uncondensed_{{*this}, {*this}, {}}, condensed_{{*this}, {*this}, {}} {}

  UncondensedT&& GetUncondensed() { return std::move(uncondensed_); }
  CondensedT&& GetCondensed() { return std::move(condensed_); }
};

using MATViewStorage = MATStorageImpl::CondensedT;
