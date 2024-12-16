#pragma once

#include "larch/mat_conversion.hpp"

#define MV_UA_NODE_ID 100000

template <typename DAGStorageType, typename DAGViewType>
struct CondensedViewBase : DefaultViewBase<DAGStorageType, DAGViewType> {
  static constexpr inline bool is_condensed = true;
};

template <typename DAGStorageType, typename DAGViewType>
struct UncondensedViewBase : DefaultViewBase<DAGStorageType, DAGViewType> {
  static constexpr inline bool is_condensed = false;
};

namespace {

template <typename T, typename = void>
struct CheckIsCondensedHelper : std::false_type {};

template <typename Storage, template <typename, typename> typename Base>
struct CheckIsCondensedHelper<
    DAGView<Storage, Base>,
    std::void_t<decltype(DAGView<Storage, Base>::BaseType::is_condensed)>>
    : std::true_type {};

template <typename T, typename = void>
struct CheckIsCondensedExtendedHelper : std::false_type {};

template <typename Storage, template <typename, typename> typename Base>
struct CheckIsCondensedExtendedHelper<
    DAGView<Storage, Base>,
    std::void_t<decltype(DAGView<Storage, Base>::StorageType::TargetView::BaseType::
                             is_condensed)>> : std::true_type {};

}  // namespace

template <typename T, typename = void>
struct CheckIsCondensed : std::false_type {};  // TODO USE_MAT_VIEW: incomplete

template <typename T>
struct CheckIsCondensed<T, std::enable_if_t<CheckIsCondensedHelper<T>::value>>
    : std::bool_constant<T::BaseType::is_condensed> {};

template <typename T>
struct CheckIsCondensed<T, std::enable_if_t<CheckIsCondensedExtendedHelper<T>::value and
                                            not CheckIsCondensedHelper<T>::value>>
    : std::bool_constant<T::StorageType::TargetView::BaseType::is_condensed> {};

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
  std::map<NodeId, std::string> node_id_to_sampleid_map_;
  std::map<std::string, size_t> sampleid_to_mat_node_id_map_;
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

  auto GetNodeFromMAT(MAT::Node* mat_node) const {
    Assert(mat_node != nullptr);
    auto& dag = static_cast<const CRTP&>(*this);
    return typename CRTP::NodeView{dag, NodeId{mat_node->node_id}};
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
    size_t ua_node_id = MV_UA_NODE_ID;
    size_t num_nodes = mat->get_size_upper();
    Assert(mat->get_node(ua_node_id) == nullptr);
    node_storage.ua_node_id_ = NodeId{ua_node_id};
    edge_storage.ua_node_id_ = NodeId{ua_node_id};
    auto& cn = node_storage.condensed_nodes_;

    size_t last_uncondensed_node_id = 1;
    for (size_t ict = 1; ict < num_nodes + 1; ict++) {
      if (mat->get_node(ict) == nullptr) {
        last_uncondensed_node_id = ict;
        break;
      }
    }
    for (auto& [i, j] : mat->condensed_nodes) {
      auto& nodes = cn[NodeId{i}];
      auto uncondensed_dag_nodes_id = NodeId{i};
      for (auto& k : j) {
        nodes.push_back(k);
        node_storage.reversed_condensed_nodes_.insert_or_assign(k, mat->get_node(i));
        if (k != j.front()) {
          uncondensed_dag_nodes_id = NodeId{last_uncondensed_node_id};
          last_uncondensed_node_id++;
          for (size_t ict = last_uncondensed_node_id; ict < num_nodes + 1; ict++) {
            if (mat->get_node(ict) == nullptr) {
              break;
            } else {
              last_uncondensed_node_id++;
            }
          }
        }
        node_storage.node_id_to_sampleid_map_.insert_or_assign(uncondensed_dag_nodes_id,
                                                               k);
        node_storage.sampleid_to_mat_node_id_map_.insert_or_assign(
            k, uncondensed_dag_nodes_id.value);
      }
      node_storage.condensed_nodes_count_ += nodes.size();
    }
    edge_storage.condensed_nodes_ = node_storage.condensed_nodes_;
    edge_storage.reversed_condensed_nodes_ = node_storage.reversed_condensed_nodes_;
    edge_storage.node_id_to_sampleid_map_ = node_storage.node_id_to_sampleid_map_;
    edge_storage.sampleid_to_mat_node_id_map_ =
        node_storage.sampleid_to_mat_node_id_map_;
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
    auto [dag_node, mat, mat_node, is_ua] = access();
    return GetChildren() | Transform::GetId() |
           ranges::views::transform([dag_node](EdgeId child) {
             return ranges::views::iota(child.value, child.value + 1) |
                    Transform::ToId<Component::Edge>() |
                    Transform::ToEdges(dag_node.GetDAG());
           });
  }

  auto GetClade(CladeIdx clade) const {
    auto [dag_node, mat, mat_node, is_ua] = access();
    size_t node_id = MV_UA_NODE_ID;  // mat.root->node_id;
    if constexpr (CheckIsCondensed<decltype(dag_node.GetDAG())>::value) {
      if (not is_ua) {
        Assert(mat_node != nullptr);
        node_id = mat_node->children.at(clade.value)->node_id;
      }
    } else {
      if (not is_ua) {
        size_t ctr = 0;
        for (auto&& i : GetChildren() | Transform::GetChild()) {
          if (ctr == clade.value) {
            node_id = i.GetId().value;
            break;
          }
          ++ctr;
        }
        Assert(node_id != mat.root->node_id);
      }
    }
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
      Assert(mat_node != nullptr);
      return mat_node->children.size();
    }
    // if current node is condensed, then it's a leaf node;
    if (mat_node == nullptr) {
      return 0;
    }
    // TODO: add case for uncondensed nodes
    //  something like this, only better:
    auto& storage = dag_node.template GetFeatureExtraStorage<MATNodeStorage>();
    size_t ct = 0;
    for (auto* c : mat_node->children) {
      auto c_node_id = c->node_id;
      auto cn_id_iter = storage.condensed_nodes_.find(NodeId{c_node_id});
      if (cn_id_iter != storage.condensed_nodes_.end()) {
        ct += storage.condensed_nodes_.at(NodeId{c_node_id}).size();
      } else {
        ++ct;
      }
    }
    return ct;
  }

  auto GetChildren() const {
    auto [dag_node, mat, mat_node, is_ua] = access();
    auto dag = dag_node.GetDAG();

    using Node = typename decltype(dag)::NodeView;
    using Edge = typename decltype(dag)::EdgeView;

    struct Result : ranges::view_facade<Result> {
      Result() : done(true) {}
      Result(std::tuple<Node, std::reference_wrapper<const MAT::Tree>, MAT::Node*, bool>
                 access,
             size_t count)
          : access_{access}, clades_count{count} {
        if (mat_node() == nullptr or mat_node()->children.empty()) {
          if (is_ua()) {
            c_node_id = dag_node()
                            .template GetFeatureExtraStorage<MATNodeStorage>()
                            .mat_tree_->root->node_id;
          } else {
            done = true;
          }
        } else {
          child_iter = mat_node()->children.begin();
          c_node_id = (*child_iter)->node_id;
          if (not condensed()) {
            // check if this child node needs uncondensing
            cn_id_iter = storage().condensed_nodes_.find(NodeId{c_node_id});
            if (cn_id_iter != storage().condensed_nodes_.end()) {
              cn_str = cn_id_iter->second.begin();
              c_node_id = storage().sampleid_to_mat_node_id_map_.at(*cn_str);
            }
          }
        }
      }

      using Storage = decltype(std::declval<Node>()
                                   .template GetFeatureExtraStorage<MATNodeStorage>());

      bool empty() const { return clades_count == 0; }
      size_t size() const { return clades_count; }

     private:
      friend ranges::range_access;

      Edge read() const {
        if (is_ua()) {
          return Edge{dag_node().GetDAG(), EdgeId{MV_UA_NODE_ID}};
        } else {
          return Edge{dag_node().GetDAG(), EdgeId{c_node_id}};
        }
      }

      bool equal(ranges::default_sentinel_t) const { return is_done(); }

      void next() {
        if (is_done()) {
          return;
        }
        if (is_ua()) {
          done = true;
          return;
        }
        if (condensed()) {
          if (++child_iter != mat_node()->children.end()) {
            c_node_id = (*child_iter)->node_id;
          } else {
            done = true;
          }
        } else {
          if (not advance_condensed()) {
            if (not advance_child()) {
              done = true;
            } else {
              c_node_id = (*child_iter)->node_id;
              cn_id_iter = storage().condensed_nodes_.find(NodeId{c_node_id});
              if (cn_id_iter != storage().condensed_nodes_.end()) {
                cn_str = cn_id_iter->second.begin();
                c_node_id = storage().sampleid_to_mat_node_id_map_.at(*cn_str);
              }
            }
          }
          if (not done) {
            if (cn_id_iter == storage().condensed_nodes_.end() or
                cn_str == cn_id_iter->second.end()) {
              c_node_id = (*child_iter)->node_id;
            } else {
              c_node_id = storage().sampleid_to_mat_node_id_map_.at(*cn_str);
            }
          }
        }
      }

      bool advance_condensed() {
        if (cn_id_iter == storage().condensed_nodes_.end()) {
          return false;
        }
        if (cn_str == cn_id_iter->second.end()) {
          return false;
        }
        return ++cn_str != cn_id_iter->second.end();
      }

      bool advance_child() {
        if (child_iter == mat_node()->children.end()) {
          return false;
        }
        return ++child_iter != mat_node()->children.end();
      }

      static constexpr bool condensed() {
        return CheckIsCondensed<decltype(std::get<0>(access_.value()).GetDAG())>::value;
      }

      Node dag_node() const { return std::get<0>(access_.value()); }

      MAT::Node* mat_node() const { return std::get<2>(access_.value()); }

      Storage& storage() const {
        return dag_node().template GetFeatureExtraStorage<MATNodeStorage>();
      }

      bool is_ua() const { return std::get<3>(access_.value()); }

      bool is_done() const { return done or not access_.has_value(); }

      std::optional<
          std::tuple<Node, std::reference_wrapper<const MAT::Tree>, MAT::Node*, bool>>
          access_;
      decltype(std::declval<MAT::Node>().children.begin()) child_iter = {};
      decltype(std::declval<Storage>().condensed_nodes_.begin()) cn_id_iter = {};
      decltype(cn_id_iter->second.begin()) cn_str = {};
      size_t c_node_id = NoId;
      bool done = false;
      size_t clades_count = 0;
    };

    return Result{access(), GetCladesCount()};
  }

  auto GetSingleParent() const {
    auto [dag_node, mat, mat_node, is_ua] = access();
    Assert(not is_ua);
    auto dag = dag_node.GetDAG();
    if constexpr (not CheckIsCondensed<decltype(dag_node.GetDAG())>::value) {
      if (mat_node == nullptr) {
        auto& storage = dag_node.template GetFeatureExtraStorage<MATNodeStorage>();
        auto cn_id_str = storage.node_id_to_sampleid_map_.find(dag_node.GetId());
        Assert(cn_id_str != storage.node_id_to_sampleid_map_.end());
        auto condensed_mat_node_str =
            storage.node_id_to_sampleid_map_.at(dag_node.GetId());
        auto condensed_mat_node_iter =
            storage.reversed_condensed_nodes_.find(condensed_mat_node_str);
        Assert(condensed_mat_node_iter != storage.reversed_condensed_nodes_.end());
        auto& condensed_mat_node =
            storage.reversed_condensed_nodes_.at(condensed_mat_node_str);
        Assert(condensed_mat_node->parent != nullptr);
        EdgeId parent{condensed_mat_node->parent == nullptr ? mat.root->node_id
                                                            : dag_node.GetId().value};
        return typename decltype(dag)::EdgeView{dag, parent};
      }
    }
    Assert(mat_node != nullptr);
    EdgeId parent{mat_node->parent == nullptr ? mat.root->node_id : mat_node->node_id};
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

  // bool IsTreeRoot() const;

  bool IsTreeRoot() const {
    auto [dag_node, mat, mat_node, is_ua] = access();
    if (mat_node != nullptr) {
      return mat_node->node_id == mat.root->node_id;
    }
    return false;
  }

  MAT::Node* GetMATNode() const {
    auto [dag_node, mat, mat_node, is_ua] = access();
    return mat_node;
  }

  bool HaveMATNode() const { return GetMATNode() != nullptr; }

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
    if (mat_node == nullptr) {
      auto& storage = dag_node.template GetFeatureExtraStorage<MATNodeStorage>();
      auto cn_id_str = storage.node_id_to_sampleid_map_.find(dag_node.GetId());
      Assert(cn_id_str != storage.node_id_to_sampleid_map_.end());
      auto condensed_mat_node_str =
          storage.node_id_to_sampleid_map_.at(dag_node.GetId());
      auto condensed_mat_node_iter =
          storage.reversed_condensed_nodes_.find(condensed_mat_node_str);
      Assert(condensed_mat_node_iter != storage.reversed_condensed_nodes_.end());
      auto& condensed_mat_node =
          storage.reversed_condensed_nodes_.at(condensed_mat_node_str);
      Assert(condensed_mat_node->parent != nullptr);
      return condensed_mat_node->parent->node_id == node.value;
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
    if (mat_node == nullptr) {
      return false;
    }
    if constexpr (CheckIsCondensed<decltype(dag_node.GetDAG())>::value) {
      for (auto* i : mat_node->children) {
        if (i->node_id == node.value) {
          return true;
        }
      }
      return false;
    }
    size_t node_value = node.value;
    // if the NodeId in question is one that has been condensed, we need to match the
    // value of its reference condensed node value to the MAT's condensed counterpart
    auto& storage = dag_node.template GetFeatureExtraStorage<MATNodeStorage>();
    auto cn_id_str = storage.node_id_to_sampleid_map_.find(node);
    if (cn_id_str != storage.node_id_to_sampleid_map_.end()) {
      auto condensed_mat_node_str = storage.node_id_to_sampleid_map_.at(node);
      auto condensed_mat_node =
          storage.reversed_condensed_nodes_.at(condensed_mat_node_str);
      node_value = condensed_mat_node->node_id;
    }
    for (auto* i : mat_node->children) {
      if (i->node_id == node_value) {
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
    auto& mat = *dag.template GetFeatureExtraStorage<Component::Node, MATNodeStorage>()
                     .mat_tree_;
    auto* mat_node = mat.get_node(id.value);
    bool is_ua = dag.template GetFeatureExtraStorage<Component::Node, MATNodeStorage>()
                     .ua_node_id_ == id;
    if ((not is_ua) and (mat_node == nullptr)) {
      auto& storage = dag_node.template GetFeatureExtraStorage<MATNodeStorage>();
      if (storage.node_id_to_sampleid_map_.find(id) ==
          storage.node_id_to_sampleid_map_.end()) {
        is_ua = true;
      }
    }
    return std::make_tuple(dag_node, std::ref(mat), mat_node, is_ua);
  }
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<MATNodeStorage, CRTP, Tag> {
  void ClearConnections() const {
    // TODO USE_MAT_VIEW
  }
  void SetSingleParent(EdgeId /*parent*/) const {
    // TODO USE_MAT_VIEW
  }
  void AddEdge(CladeIdx /*clade*/, EdgeId /*id*/, bool /*this_node_is_parent*/) const {
    // TODO USE_MAT_VIEW
  }
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
  std::map<NodeId, std::vector<std::string>> condensed_nodes_;
  std::map<std::string, MAT::Node*> reversed_condensed_nodes_;
  std::map<NodeId, std::string> node_id_to_sampleid_map_;
  std::map<std::string, size_t> sampleid_to_mat_node_id_map_;
  size_t condensed_nodes_count_ = 0;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<MATEdgeStorage, CRTP, Tag> {
  auto GetParent() const {
    auto [dag_edge, mat, mat_node, is_ua] = access();
    if (is_ua) {
      return dag_edge.GetDAG().GetRoot();
    }
    // if it's an edge above a condensed node
    auto& storage = dag_edge.template GetFeatureExtraStorage<MATEdgeStorage>();
    if (mat_node == nullptr) {
      [[maybe_unused]] auto cn_id_str =
          storage.node_id_to_sampleid_map_.find(dag_edge.GetChildId());
      Assert(cn_id_str != storage.node_id_to_sampleid_map_.end());
      [[maybe_unused]] auto condensed_mat_node_str =
          storage.node_id_to_sampleid_map_.at(dag_edge.GetChildId());
      [[maybe_unused]] auto condensed_mat_node_iter =
          storage.reversed_condensed_nodes_.find(condensed_mat_node_str);
      Assert(condensed_mat_node_iter != storage.reversed_condensed_nodes_.end());
      [[maybe_unused]] auto& condensed_mat_node =
          storage.reversed_condensed_nodes_.at(condensed_mat_node_str);
      Assert(condensed_mat_node->parent != nullptr);
      return dag_edge.GetDAG().Get(NodeId{condensed_mat_node->parent->node_id});
    }
    Assert(mat_node != nullptr);
    Assert(mat_node->parent != nullptr);
    return dag_edge.GetDAG().Get(NodeId{mat_node->parent->node_id});
  }

  auto GetChild() const {
    auto [dag_edge, mat, mat_node, is_ua] = access();
    if (is_ua) {
      return dag_edge.GetDAG().Get(NodeId{mat.root->node_id});
    }
    if (mat_node != nullptr) {
      return dag_edge.GetDAG().Get(NodeId{mat_node->node_id});
    }
    return dag_edge.GetDAG().Get(dag_edge.GetChildId());
  }

  CladeIdx GetClade() const {
    auto [dag_edge, mat, mat_node, is_ua] = access();
    CladeIdx result{0};
    if (is_ua) {
      return result;
    }
    if constexpr (CheckIsCondensed<decltype(dag_edge.GetDAG())>::value) {
      Assert(mat_node != nullptr);
      Assert(mat_node->parent != nullptr);
      for (auto* i : mat_node->parent->children) {
        if (i == mat_node) {
          return result;
        }
        ++result.value;
      }
    } else {
      // if it's an edge above a condensed node
      auto& storage = dag_edge.template GetFeatureExtraStorage<MATEdgeStorage>();
      if (mat_node == nullptr) {
        auto cn_id_str = storage.node_id_to_sampleid_map_.find(dag_edge.GetChildId());
        Assert(cn_id_str != storage.node_id_to_sampleid_map_.end());
        auto condensed_mat_node_str =
            storage.node_id_to_sampleid_map_.at(dag_edge.GetChildId());
        auto condensed_mat_node_iter =
            storage.reversed_condensed_nodes_.find(condensed_mat_node_str);
        Assert(condensed_mat_node_iter != storage.reversed_condensed_nodes_.end());
        mat_node = storage.reversed_condensed_nodes_.at(condensed_mat_node_str);
      }
      Assert(mat_node != nullptr);
      Assert(mat_node->parent != nullptr);

      for (auto* c : mat_node->parent->children) {
        auto c_node_id = c->node_id;
        if (c_node_id == dag_edge.GetId().value) {
          return result;
        }
        auto cn_id_iter = storage.condensed_nodes_.find(NodeId{c_node_id});
        if (cn_id_iter != storage.condensed_nodes_.end()) {
          for (auto cn_str : storage.condensed_nodes_.at(NodeId{c_node_id})) {
            auto current_id = storage.sampleid_to_mat_node_id_map_.at(cn_str);
            if (current_id == dag_edge.GetChildId().value) {
              return result;
            }
            ++result.value;
          }
        } else {
          ++result.value;
        }
      }
    }
    Fail("Clade not found");
  }

  auto GetParentId() const {
    auto [dag_edge, mat, mat_node, is_ua] = access();
    if (is_ua) {
      return dag_edge.GetDAG().GetRoot().GetId();
    }
    // if it's an edge above a condensed node
    if (mat_node == nullptr) {
      [[maybe_unused]] auto& storage =
          dag_edge.template GetFeatureExtraStorage<MATEdgeStorage>();
      [[maybe_unused]] auto cn_id_str =
          storage.node_id_to_sampleid_map_.find(dag_edge.GetChildId());
      Assert(cn_id_str != storage.node_id_to_sampleid_map_.end());
      [[maybe_unused]] auto condensed_mat_node_str =
          storage.node_id_to_sampleid_map_.at(dag_edge.GetChildId());
      [[maybe_unused]] auto condensed_mat_node_iter =
          storage.reversed_condensed_nodes_.find(condensed_mat_node_str);
      Assert(condensed_mat_node_iter != storage.reversed_condensed_nodes_.end());
      [[maybe_unused]] auto& condensed_mat_node =
          storage.reversed_condensed_nodes_.at(condensed_mat_node_str);
      Assert(condensed_mat_node->parent != nullptr);
      return NodeId{condensed_mat_node->parent->node_id};
    }
    Assert(mat_node->parent != nullptr);
    return NodeId{mat_node->parent->node_id};
  }

  auto GetChildId() const {
    auto [dag_edge, mat, mat_node, is_ua] = access();
    return NodeId{dag_edge.GetId().value};
  }

  std::pair<NodeId, NodeId> GetNodeIds() const {
    auto [dag_edge, mat, mat_node, is_ua] = access();
    if (is_ua) {
      return {dag_edge.GetDAG().GetRoot(), {mat.root->node_id}};
    }
    // if it's an edge above a condensed node
    if (mat_node == nullptr) {
      return {GetParentId(), GetChildId()};
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
    auto& id_storage = dag_edge.template GetFeatureExtraStorage<MATEdgeStorage>();

    if (mat_node == nullptr) {
      [[maybe_unused]] auto cn_id_str =
          id_storage.node_id_to_sampleid_map_.find(dag_edge.GetChildId());
      Assert(cn_id_str != id_storage.node_id_to_sampleid_map_.end());
      [[maybe_unused]] auto condensed_mat_node_str =
          id_storage.node_id_to_sampleid_map_.at(dag_edge.GetChildId());
      [[maybe_unused]] auto condensed_mat_node_iter =
          id_storage.reversed_condensed_nodes_.find(condensed_mat_node_str);
      Assert(condensed_mat_node_iter != id_storage.reversed_condensed_nodes_.end());
      mat_node = id_storage.reversed_condensed_nodes_.at(condensed_mat_node_str);
    }
    Assert(mat_node != nullptr);
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
    auto& mat = *dag.template GetFeatureExtraStorage<Component::Edge, MATEdgeStorage>()
                     .mat_tree_;
    Assert(id.value == MV_UA_NODE_ID or id.value < mat.get_size_upper());
    MAT::Node* mat_node = mat.get_node(id.value);

    bool is_ua =
        mat_node != nullptr ? ((mat.root->node_id == mat_node->node_id)) : false;
    if ((not is_ua) and (mat_node == nullptr)) {
      auto& storage = dag_edge.template GetFeatureExtraStorage<MATEdgeStorage>();
      // if (storage.condensed_nodes_.find(NodeId{id.value}) ==
      //     storage.condensed_nodes_.end()) {
      if (storage.node_id_to_sampleid_map_.find(NodeId{id.value}) ==
          storage.node_id_to_sampleid_map_.end()) {
        is_ua = true;
      }
    }

    return std::make_tuple(dag_edge, std::ref(mat), mat_node, is_ua);
  }
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<MATEdgeStorage, CRTP, Tag> {
  void Set(NodeId /*parent*/, NodeId /*child*/, CladeIdx /*clade*/) const {
    // TODO USE_MAT_VIEW
  }
};

template <Component C, bool>
struct MATElementsContainerBase {
  using ElementStorageT =
      std::conditional_t<C == Component::Node, MATNodeStorage, MATEdgeStorage>;
  using FeatureTypes =
      std::conditional_t<C == Component::Node, std::tuple<ElementStorageT, Neighbors>,
                         std::tuple<ElementStorageT, Endpoints>>;
  using AllFeatureTypes = FeatureTypes;

  static constexpr IdContinuity id_continuity = IdContinuity::Sparse;

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

  MATElementsContainerBase() = default;

  template <typename VT>
  size_t GetCount() const {
    size_t count = GetMAT().get_node_idx() + 1;

    if constexpr (C == Component::Edge) {
      Assert(count > 0);
      count -= 1;
    }
    if constexpr (not CheckIsCondensed<VT>::value) {
      count -= extra_storage_.condensed_nodes_.size();
      count += extra_storage_.condensed_nodes_count_;
    }
    return count;
  }

  template <typename VT>
  Id<C> GetNextAvailableId() const {
    return {GetCount<VT>()};
  }

  template <typename VT>
  bool ContainsId(Id<C> id) const {
    return id.value < GetCount<VT>();
  }

  template <typename Feature>
  const auto& GetFeatureStorage(Id<C> id) const {
    if constexpr (C == Component::Node) {
      if constexpr (std::is_same_v<Feature, Neighbors>) {
        return Neighbors{};  // TODO USE_MAT_VIEW
      }
      // TODO
    } else {
      if (features_storage_.empty()) {
        features_storage_.resize(GetMAT().get_size_upper());
      }
      return std::get<Feature>(features_storage_.at(id));
    }
  }

  template <typename Feature>
  auto&& GetFeatureStorage(Id<C> id) {
    if constexpr (C == Component::Node) {
      if constexpr (std::is_same_v<Feature, Neighbors>) {
        static Neighbors neighbors;
        return neighbors;  // TODO USE_MAT_VIEW
      }
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

  template <typename VT>
  auto All() const {
    size_t iota_max = GetCount<VT>() + 1;
    // getcount returns the number of nodes/edges in the tree.
    // if view is condensed, we need to check all possible node ids,
    // which means checking the max number of nodes, condensed or not.

    static constexpr const bool is_condensed = CheckIsCondensed<VT>::value;
    if constexpr (is_condensed) {
      iota_max -= extra_storage_.condensed_nodes_.size();
      iota_max += extra_storage_.condensed_nodes_count_;
    }
    if constexpr (C == Component::Node) {
      return ranges::views::iota(size_t{0}, iota_max) |
             ranges::views::transform([](size_t i) -> size_t {
               if (i == 0) {
                 return MV_UA_NODE_ID;
               } else {
                 return i - 1;
               }
             }) |
             ranges::views::filter([this](size_t i) {
               if constexpr (is_condensed) {
                 return GetMAT().get_node(i) != nullptr or
                        i == extra_storage_.ua_node_id_.value;
               }
               return GetMAT().get_node(i) != nullptr or
                      i == extra_storage_.ua_node_id_.value or
                      (extra_storage_.node_id_to_sampleid_map_.find(NodeId{i}) !=
                       extra_storage_.node_id_to_sampleid_map_.end());
             }) |
             ranges::views::transform([](size_t i) -> NodeId { return {i}; });
    } else {
      return ranges::views::iota(size_t{0}, iota_max + 1) |
             ranges::views::filter([this](size_t i) {
               if constexpr (is_condensed) {
                 return GetMAT().get_node(i) != nullptr and
                        i != extra_storage_.ua_node_id_.value;
               }
               return (GetMAT().get_node(i) != nullptr and
                       i != extra_storage_.ua_node_id_.value) or
                      (extra_storage_.node_id_to_sampleid_map_.find(NodeId{i}) !=
                       extra_storage_.node_id_to_sampleid_map_.end());
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

  mutable IdContainer<Id<C>, AllFeatureTypes, id_continuity> features_storage_ = {};
  ExtraFeatureStorage<ElementStorageT> extra_storage_ = {};
};

class UncondensedNodesContainer
    : public MATElementsContainerBase<Component::Node, false> {
 public:
  UncondensedNodesContainer() : MATElementsContainerBase<Component::Node, false>{} {}
};

class UncondensedEdgesContainer
    : public MATElementsContainerBase<Component::Edge, false> {
 public:
  UncondensedEdgesContainer() : MATElementsContainerBase<Component::Edge, false>{} {}
};

class CondensedNodesContainer : public MATElementsContainerBase<Component::Node, true> {
 public:
  CondensedNodesContainer() : MATElementsContainerBase<Component::Node, true>{} {}
};

class CondensedEdgesContainer : public MATElementsContainerBase<Component::Edge, true> {
 public:
  CondensedEdgesContainer() : MATElementsContainerBase<Component::Edge, true>{} {}
};

using UncondensedMATViewStorage =
    DAGStorage<void, UncondensedNodesContainer, UncondensedEdgesContainer,
               ExtraStorage<Connections>, UncondensedViewBase>;

using CondensedMATViewStorage =
    DAGStorage<void, CondensedNodesContainer, CondensedEdgesContainer,
               ExtraStorage<Connections>, CondensedViewBase>;

using CondensedMADAGStorage =
    ExtendStorageType<void, CondensedMATViewStorage, Extend::Nodes<CompactGenome>,
                      Extend::DAG<ReferenceSequence>, Extend::Empty<>,
                      CondensedViewBase>;

using UncondensedMADAGStorage =
    ExtendStorageType<void, UncondensedMATViewStorage,
                      Extend::Nodes<CompactGenome, Deduplicate<SampleId>>,
                      Extend::DAG<ReferenceSequence>, Extend::Empty<>,
                      UncondensedViewBase>;

using CondensedMergeDAGStorage = ExtendStorageType<
    void, CondensedMATViewStorage, Extend::Nodes<Deduplicate<CompactGenome>>,
    Extend::DAG<ReferenceSequence>, Extend::Empty<>, CondensedViewBase>;

using UncondensedMergeDAGStorage =
    ExtendStorageType<void, UncondensedMATViewStorage,
                      Extend::Nodes<Deduplicate<CompactGenome>, Deduplicate<SampleId>>,
                      Extend::DAG<ReferenceSequence>, Extend::Empty<>,
                      UncondensedViewBase>;
