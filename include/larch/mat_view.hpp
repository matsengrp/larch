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
    auto& dag = static_cast<const CRTP&>(*this);
    dag.template GetFeatureExtraStorage<Component::Node, MATNodeStorage>().mat_tree_ =
        mat;
    dag.template GetFeatureExtraStorage<Component::Edge, MATEdgeStorage>().mat_tree_ =
        mat;
  }
};

template <typename CRTP, typename Tag>
struct FeatureConstView<MATNodeStorage, CRTP, Tag> {
  auto GetParents() const {
    auto [dag_node, mat, mat_node] = access();
    return ranges::views::iota(size_t{0}, GetParentsCount()) |
           ranges::views::transform([this](size_t i) { return GetSingleParent(); });
  }
  auto GetClades() const {
    auto [dag_node, mat, mat_node] = access();
    return ranges::views::iota(size_t{0}, mat_node->children.size()) |
           ranges::views::transform([this](size_t i) { return GetClade(CladeIdx{i}); });
  }

  auto GetClade(CladeIdx clade) const {
    auto [dag_node, mat, mat_node] = access();
    size_t node_id = mat_node->children.at(clade.value)->node_id;
    return ranges::views::iota(node_id, node_id + 1) |
           Transform::ToId<Component::Edge>() | Transform::ToEdges(dag_node.GetDAG());
  }

  size_t GetParentsCount() const {
    auto [dag_node, mat, mat_node] = access();
    if (mat_node->parent == nullptr) {
      return 0;
    } else {
      return 1;
    }
  }

  size_t GetCladesCount() const {
    auto [dag_node, mat, mat_node] = access();
    return mat_node->children.size();
  }

  auto GetChildren() const {
    auto [dag_node, mat, mat_node] = access();
    auto dag = dag_node.GetDAG();
    return mat_node->children | ranges::views::transform([dag](MAT::Node* i) {
             return typename decltype(dag)::EdgeView{dag, EdgeId{i->node_id}};
           });
  }

  auto GetSingleParent() const {
    auto [dag_node, mat, mat_node] = access();
    Assert(mat_node->parent != nullptr);
    return Transform::ToEdges(dag_node.GetDAG())(EdgeId{mat_node->parent->node_id});
  }

  auto GetFirstParent() const {
    Assert(GetParentsCount() == 1);
    return GetSingleParent();
  }

  auto GetFirstChild() const;
  auto GetFirstClade() const;

  bool IsUA() const { return GetParentsCount() == 0; }

  bool IsTreeRoot() const;

  bool IsLeaf() const { return GetCladesCount() == 0; }

  auto GetLeafsBelow() const;
  void Validate(bool recursive = false, bool allow_dag = false) const;
  auto GetParentNodes() const;
  auto GetChildNodes() const;
  bool ContainsParent(NodeId node) const;
  bool ContainsChild(NodeId node) const;
  std::string ParentsToString() const;
  std::string ChildrenToString() const;

 private:
  auto access() const {
    auto dag_node = static_cast<const CRTP&>(*this);
    NodeId id = dag_node.GetId();
    auto dag = dag_node.GetDAG();
    auto& mat = dag.GetMAT();
    auto* mat_node = mat.get_node(id.value);
    Assert(mat_node != nullptr);
    return std::make_tuple(dag_node, std::ref(mat), mat_node);
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
           ranges::views::filter(
               [this](size_t i) { return GetMAT().get_node(i) != nullptr; }) |
           ranges::views::transform([](size_t i) -> NodeId { return {i}; });
  }

 private:
  MAT::Tree& GetMAT() {
    Assert(extra_node_storage_.mat_tree_ != nullptr);
    return *extra_node_storage_.mat_tree_;
  }

  const MAT::Tree& GetMAT() const {
    Assert(extra_node_storage_.mat_tree_ != nullptr);
    return *extra_node_storage_.mat_tree_;
  }

  ExtraFeatureStorage<MATNodeStorage> extra_node_storage_;
};

struct MATEdgeStorage {
  MOVE_ONLY(MATEdgeStorage);
  MATEdgeStorage() = default;

  inline MATEdgeStorage Copy() const { return {}; }
};

template <>
struct ExtraFeatureStorage<MATEdgeStorage> {
  ExtraFeatureStorage() = default;
  MOVE_ONLY(ExtraFeatureStorage);
  MAT::Tree* mat_tree_ = nullptr;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<MATEdgeStorage, CRTP, Tag> {
  auto GetParent() const {
    auto [dag_edge, mat, mat_node] = access();
    Assert(mat_node->parent != nullptr);
    return dag_edge.GetDAG().Get(NodeId{mat_node->parent->node_id});
  }

  auto GetChild() const {
    auto [dag_edge, mat, mat_node] = access();
    return dag_edge.GetDAG().Get(NodeId{mat_node->node_id});
  }

  CladeIdx GetClade() const {
    auto [dag_edge, mat, mat_node] = access();
    Assert(mat_node->parent != nullptr);
    CladeIdx result{};
    for (auto* i : mat_node->parent->children) {
      if (i == mat_node) {
        return result;
      }
      ++result.value;
    }
    Fail("Clade not found");
  }

  auto GetParentId() const {
    auto [dag_edge, mat, mat_node] = access();
    Assert(mat_node->parent != nullptr);
    return NodeId{mat_node->parent->node_id};
  }

  auto GetChildId() const {
    auto [dag_edge, mat, mat_node] = access();
    return NodeId{mat_node->node_id};
  }

  std::pair<NodeId, NodeId> GetNodeIds() const {
    auto [dag_edge, mat, mat_node] = access();
    Assert(mat_node->parent != nullptr);
    return {{mat_node->parent->node_id}, {mat_node->node_id}};
  }

  bool IsUA() const { return GetParent().IsUA(); }

  bool IsTreeRoot() const { return GetParent().IsTreeRoot(); }

  bool IsLeaf() const { return GetChild().IsLeaf(); }

  EdgeMutations GetEdgeMutations() const {
    // TODO implment
    return EdgeMutations{};
  }

 private:
  auto access() const {
    auto dag_edge = static_cast<const CRTP&>(*this);
    EdgeId id = dag_edge.GetId();
    auto dag = dag_edge.GetDAG();
    auto& mat = dag.GetMAT();
    auto* mat_node = mat.get_node(id.value);
    Assert(mat_node != nullptr);
    return std::make_tuple(dag_edge, std::ref(mat), mat_node);
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
    auto nodes_count = GetMAT().get_size_upper();
    if (nodes_count == 0) {
      return 0;
    }
    return nodes_count - 1;
  }

  EdgeId GetNextAvailableId() const { return {GetCount()}; }

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
           ranges::views::filter(
               [this](size_t i) { return GetMAT().get_node(i) != nullptr; }) |
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

  ExtraFeatureStorage<MATEdgeStorage> extra_edge_storage_;
};
