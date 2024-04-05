#pragma once

#include "larch/mat_conversion.hpp"

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
};

struct MATEdgeStorage;

template <typename CRTP>
struct ExtraFeatureMutableView<MATNodeStorage, CRTP> {
  void SetMAT(MAT::Tree* mat) const {
    Assert(mat != nullptr);
    auto& dag = static_cast<const CRTP&>(*this);
    dag.template GetFeatureExtraStorage<Component::Node, MATNodeStorage>().mat_tree_ =
        mat;
    dag.template GetFeatureExtraStorage<Component::Edge, MATEdgeStorage>().mat_tree_ =
        mat;
    size_t ua_node_id = 0;  // mat->get_size_upper();
    Assert(mat->get_node(ua_node_id) == nullptr);
    dag.template GetFeatureExtraStorage<Component::Node, MATNodeStorage>().ua_node_id_ =
        NodeId{ua_node_id};
    dag.template GetFeatureExtraStorage<Component::Edge, MATEdgeStorage>().ua_node_id_ =
        NodeId{ua_node_id};
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
    if (mat_node->parent == nullptr) {
      return false;
    }
    return mat_node->parent->node_id == node.value;
  }

  bool ContainsChild(NodeId node) const {
    auto [dag_node, mat, mat_node, is_ua] = access();
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

struct MATNodesContainer {
  using FeatureTypes = std::tuple<MATNodeStorage>;
  using AllFeatureTypes = FeatureTypes;

  static constexpr IdContinuity id_continuity = IdContinuity::Dense;

  template <typename Feature>
  static const bool contains_element_feature =
      tuple_contains_v<AllFeatureTypes, Feature>;

  template <typename CRTP>
  struct ConstElementViewBase : FeatureConstView<MATNodeStorage, CRTP> {};
  template <typename CRTP>
  struct MutableElementViewBase : FeatureMutableView<MATNodeStorage, CRTP> {
    using FeatureMutableView<MATNodeStorage, CRTP>::operator=;
  };

  template <typename CRTP>
  struct ExtraConstElementViewBase : ExtraFeatureConstView<MATNodeStorage, CRTP> {};
  template <typename CRTP>
  struct ExtraMutableElementViewBase : ExtraFeatureMutableView<MATNodeStorage, CRTP> {};

  MATNodesContainer() = default;
  MOVE_ONLY(MATNodesContainer);

  size_t GetCount() const { return GetMAT().get_size_upper(); }

  NodeId GetNextAvailableId() const { return {GetCount()}; }

  template <typename Feature>
  auto& GetFeatureExtraStorage() {
    static_assert(std::is_same_v<Feature, MATNodeStorage>);
    return extra_node_storage_;
  }

  template <typename Feature>
  const auto& GetFeatureExtraStorage() const {
    static_assert(std::is_same_v<Feature, MATNodeStorage>);
    return extra_node_storage_;
  }

  auto All() const {
    return ranges::views::iota(size_t{0}, GetCount()) |
           ranges::views::filter([this](size_t i) {
             return GetMAT().get_node(i) != nullptr or
                    i == extra_node_storage_.ua_node_id_.value;
           }) |
           ranges::views::transform([](size_t i) -> NodeId { return {i}; });
  }

 private:
  MAT::Tree& GetMAT() {
    Assert(extra_node_storage_.mat_tree_ != nullptr);
    Assert(extra_node_storage_.ua_node_id_.value != NoId);
    return *extra_node_storage_.mat_tree_;
  }

  const MAT::Tree& GetMAT() const {
    Assert(extra_node_storage_.mat_tree_ != nullptr);
    Assert(extra_node_storage_.ua_node_id_.value != NoId);
    return *extra_node_storage_.mat_tree_;
  }

  ExtraFeatureStorage<MATNodeStorage> extra_node_storage_;
};

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
};

template <typename CRTP, typename Tag>
struct FeatureConstView<MATEdgeStorage, CRTP, Tag> {
  auto GetParent() const {
    auto [dag_edge, mat, mat_node, is_ua] = access();
    if (is_ua) {
      return dag_edge.GetDAG().GetRoot();
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
    bool is_ua = mat.root == mat_node;
    return std::make_tuple(dag_edge, std::ref(mat), mat_node, is_ua);
  }
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<MATEdgeStorage, CRTP, Tag> {};

struct MATEdgesContainer {
  using FeatureTypes = std::tuple<MATEdgeStorage>;
  using AllFeatureTypes = FeatureTypes;

  static constexpr IdContinuity id_continuity = IdContinuity::Dense;

  template <typename Feature>
  static const bool contains_element_feature =
      tuple_contains_v<AllFeatureTypes, Feature>;

  template <typename CRTP>
  struct ConstElementViewBase : FeatureConstView<MATEdgeStorage, CRTP> {};
  template <typename CRTP>
  struct MutableElementViewBase : FeatureMutableView<MATEdgeStorage, CRTP> {
    using FeatureMutableView<MATEdgeStorage, CRTP>::operator=;
  };

  template <typename CRTP>
  struct ExtraConstElementViewBase : ExtraFeatureConstView<MATEdgeStorage, CRTP> {};
  template <typename CRTP>
  struct ExtraMutableElementViewBase : ExtraFeatureMutableView<MATEdgeStorage, CRTP> {};

  MATEdgesContainer() = default;
  MOVE_ONLY(MATEdgesContainer);

  size_t GetCount() const {
    Assert(GetMAT().get_size_upper() > 0);
    return GetMAT().get_size_upper() - 1;
  }

  EdgeId GetNextAvailableId() const { return {GetCount()}; }

  template <typename Feature>
  const auto& GetFeatureStorage(EdgeId id) const {
    if (features_storage_.empty()) {
      features_storage_.resize(GetMAT().get_size_upper());
    }
    return std::get<Feature>(features_storage_.at(id));
  }

  template <typename Feature>
  auto& GetFeatureStorage(EdgeId id) {
    if (features_storage_.empty()) {
      features_storage_.resize(GetMAT().get_size_upper());
    }
    return std::get<Feature>(features_storage_.at(id));
  }

  template <typename Feature>
  auto& GetFeatureExtraStorage() {
    static_assert(std::is_same_v<Feature, MATEdgeStorage>);
    return extra_edge_storage_;
  }

  template <typename Feature>
  const auto& GetFeatureExtraStorage() const {
    static_assert(std::is_same_v<Feature, MATEdgeStorage>);
    return extra_edge_storage_;
  }

  auto All() const {
    return ranges::views::iota(size_t{0}, GetCount()) |
           ranges::views::filter([this](size_t i) {
             return i != extra_edge_storage_.ua_node_id_.value and
                    GetMAT().get_node(i) != nullptr;
           }) |
           ranges::views::transform([](size_t i) -> EdgeId { return {i}; });
  }

 private:
  MAT::Tree& GetMAT() {
    Assert(extra_edge_storage_.mat_tree_ != nullptr);
    return *extra_edge_storage_.mat_tree_;
  }

  const MAT::Tree& GetMAT() const {
    Assert(extra_edge_storage_.mat_tree_ != nullptr);
    return *extra_edge_storage_.mat_tree_;
  }

  mutable IdContainer<EdgeId, AllFeatureTypes, id_continuity> features_storage_;
  ExtraFeatureStorage<MATEdgeStorage> extra_edge_storage_;
};

using MATViewStorage =
    DAGStorage<void, MATNodesContainer, MATEdgesContainer, ExtraStorage<Connections>>;
