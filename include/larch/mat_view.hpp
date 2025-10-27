#pragma once

#include "larch/mat_conversion.hpp"
#include "larch/debug.hpp"

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

/**
 * @brief Empty storage tag for MAT nodes in the view system.
 *
 * MATNodeStorage is an empty marker type used to tag DAG nodes as part of a MAT view.
 * It doesn't store any data itself since all node information is accessed directly
 * from the underlying MAT structure. The view system uses this tag to identify and
 * provide MAT-specific operations on DAG nodes.
 */
struct MATNodeStorage {
  MOVE_ONLY(MATNodeStorage);
  MATNodeStorage() = default;

  template <typename CRTP>
  inline MATNodeStorage Copy(const CRTP*) const {
    return {};
  }
};

/**
 * @brief Empty storage tag for MAT edges in the view system.
 *
 * MATEdgeStorage is an empty marker type used to tag DAG edges as part of a MAT view.
 * It doesn't store any data itself since all edge information is accessed directly
 * from the underlying MAT structure. The view system uses this tag to identify and
 * provide MAT-specific operations on DAG edges.
 */
struct MATEdgeStorage {
  MOVE_ONLY(MATEdgeStorage);
  MATEdgeStorage() = default;

  template <typename CRTP>
  inline MATEdgeStorage Copy(const CRTP*) const {
    return {};
  }
};

/**
 * @brief Range adaptor for iterating over children of a MAT node in a DAG view.
 *
 * MATChildrenRange provides a lazy range interface for accessing the children of a node
 * in a MAT view. It handles both condensed and uncondensed nodes, automatically
 * expanding condensed nodes (which represent multiple samples) into their individual
 * children when iterating in an uncondensed view. This class implements the
 * ranges::view_facade interface for seamless integration with range-based algorithms.
 *
 * @tparam DAG The DAG view type being iterated over
 */
template <typename DAG>
struct MATChildrenRange : ranges::view_facade<MATChildrenRange<DAG>> {
  using Node = typename DAG::NodeView;
  using Edge = typename DAG::EdgeView;

  MATChildrenRange() : done(true) {}
  MATChildrenRange(
      std::tuple<Node, std::reference_wrapper<const MAT::Tree>, MAT::Node*, bool>
          access,
      size_t count)
      : access_{access}, clades_count{count} {
    if (mat_node() == nullptr or mat_node()->children.empty()) {
      if (is_ua()) {
        c_node_id = dag_node()
                        .template GetFeatureExtraStorage<MATNodeStorage>()
                        .get()
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
                               .template GetFeatureExtraStorage<MATNodeStorage>()
                               .get());

  bool empty() const { return clades_count == 0; }

  size_t size() const { return clades_count; }

 private:
  friend ranges::range_access;

  EdgeId read() const { return EdgeId{c_node_id}; }

  bool equal(ranges::default_sentinel_t) const { return is_done(); }

  bool equal(const MATChildrenRange& other) const {
    if (is_done()) {
      return other.is_done();
    }
    if (other.is_done()) {
      return is_done();
    }
    return c_node_id == other.c_node_id;
  }

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
    return dag_node().template GetFeatureExtraStorage<MATNodeStorage>().get();
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

/**
 * @brief Provides neighbor relationship navigation for nodes in a MAT view.
 *
 * MATNeighbors implements the Neighbors interface for MAT views, providing access to
 * parent and child relationships of nodes. It handles the translation between MAT node
 * structures and DAG neighbor concepts, including support for both condensed and
 * uncondensed views. The class provides methods to query parents, clades (groups of
 * children), and leaf nodes below a given node in the tree hierarchy.
 */
struct MATNeighbors : Neighbors {
  template <typename CRTP>
  inline DAGNeighbors Copy(const CRTP* crtp) const {
    DAGNeighbors result;
    result.GetParentsMutable(crtp) = ranges::to_vector(GetParents(crtp));
    for (auto i : GetClades(crtp)) {
      result.GetCladesMutable(crtp).push_back(std::move(ranges::to_vector(i)));
    }
    return result;
  }

  template <typename CRTP>
  auto GetParents(const CRTP* crtp) const {
    static_assert(CRTP::role == Role::View);
    const size_t count = GetParentsCount(crtp);
    const EdgeId result = count == 0 ? EdgeId{NoId} : GetSingleParent(crtp);
    return ranges::views::iota(size_t{0}, count) |
           ranges::views::transform([result](size_t) { return result; });
  }

  template <typename CRTP>
  auto GetClades(const CRTP* crtp) const {
    return GetChildren(crtp) | ranges::views::transform([](EdgeId child) {
             return ranges::views::iota(child.value, child.value + 1) |
                    Transform::ToId<Component::Edge>();
           });
  }

  template <typename CRTP>
  auto GetLeafsBelow(const CRTP*) const {
    Fail("Leafs below not available on MATView");
    static const std::array<NodeId, 1> empty = {{{NoId}}};
    return empty | ranges::views::take(0);
  }

  template <typename CRTP>
  auto& GetParentsMutable(const CRTP*) {
    Fail("Can't modify MATNeighbors");
    return *Unreachable<std::vector<EdgeId>>();
  }

  template <typename CRTP>
  auto& GetCladesMutable(const CRTP*) {
    Fail("Can't modify MATNeighbors");
    return *Unreachable<std::vector<std::vector<EdgeId>>>();
  }

  template <typename CRTP>
  auto& GetLeafsBelowMutable(const CRTP*) {
    Fail("Can't modify MATNeighbors");
    return *Unreachable<std::vector<std::vector<NodeId>>>();
  }

 private:
  template <typename CRTP>
  EdgeId GetSingleParent(const CRTP* crtp) const {
    // DebugUse(crtp);
    auto [dag_node, mat, mat_node, is_ua] = access_const(crtp);
    Assert(not is_ua);
    if constexpr (not CheckIsCondensed<decltype(dag_node.GetDAG())>::value) {
      if (mat_node == nullptr) {
        auto& storage =
            dag_node.template GetFeatureExtraStorage<MATNodeStorage>().get();
        [[maybe_unused]] auto cn_id_str =
            storage.node_id_to_sampleid_map_.find(dag_node.GetId());
        Assert(cn_id_str != storage.node_id_to_sampleid_map_.end());
        auto condensed_mat_node_str =
            storage.node_id_to_sampleid_map_.at(dag_node.GetId());
        [[maybe_unused]] auto condensed_mat_node_iter =
            storage.reversed_condensed_nodes_.find(condensed_mat_node_str);
        Assert(condensed_mat_node_iter != storage.reversed_condensed_nodes_.end());
        auto& condensed_mat_node =
            storage.reversed_condensed_nodes_.at(condensed_mat_node_str);
        Assert(condensed_mat_node->parent != nullptr);
        EdgeId parent{condensed_mat_node->parent == nullptr ? mat.root->node_id
                                                            : dag_node.GetId().value};
        return parent;
      }
    }
    Assert(mat_node != nullptr);
    EdgeId parent{mat_node->parent == nullptr ? mat.root->node_id : mat_node->node_id};
    return parent;
  }

  template <typename CRTP>
  size_t GetParentsCount(const CRTP* crtp) const {
    auto [dag_node, mat, mat_node, is_ua] = access_const(crtp);
    if (is_ua) {
      return 0;
    } else {
      return 1;
    }
  }

  template <typename CRTP>
  auto GetChildren(const CRTP* crtp) const {
    auto [dag_node, mat, mat_node, is_ua] = access_const(crtp);
    auto dag = dag_node.GetDAG();
    return MATChildrenRange<decltype(dag)>{access_const(crtp), GetCladesCount(crtp)} |
           ranges::views::common;
  }

  template <typename CRTP>
  size_t GetCladesCount(const CRTP* crtp) const {
    auto [dag_node, mat, mat_node, is_ua] = access_const(crtp);
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
    auto& storage = dag_node.template GetFeatureExtraStorage<MATNodeStorage>().get();
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

  template <typename CRTP>
  auto access_const(const CRTP* crtp) const {
    return static_cast<const FeatureConstView<MATNodeStorage, CRTP, MATNodeStorage>&>(
               *crtp)
        .access();
  }

  template <typename CRTP>
  auto access_mutable(const CRTP* crtp) {
    return static_cast<const FeatureMutableView<MATNodeStorage, CRTP, MATNodeStorage>&>(
               *crtp)
        .access();
  }
};

template <>
struct OverlayFeatureType<MATNeighbors> {
  using store_type = DAGNeighbors;
  using const_view_type = std::variant<std::reference_wrapper<const DAGNeighbors>,
                                       std::reference_wrapper<const MATNeighbors>>;
  using mutable_view_type = std::variant<std::reference_wrapper<DAGNeighbors>,
                                         std::reference_wrapper<MATNeighbors>>;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<MATNeighbors, CRTP, Tag>
    : FeatureConstView<Neighbors, CRTP, Tag> {
  bool IsCondensedInMAT() { return CheckIsCondensed<CRTP>::value; }
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<MATNeighbors, CRTP, Tag>
    : FeatureMutableView<Neighbors, CRTP, Tag> {};

/**
 * @brief Provides endpoint information for edges in a MAT view.
 *
 * MATEndpoints implements the Endpoints interface for MAT views, providing access to
 * the parent node, child node, and clade index of edges. It handles the mapping between
 * MAT edge representations and DAG edge endpoints, including special handling for
 * condensed nodes that represent multiple samples. The class correctly identifies edge
 * endpoints even when nodes have been uncondensed from the original MAT structure.
 */
struct MATEndpoints : Endpoints {
  template <typename CRTP>
  inline DAGEndpoints Copy(const CRTP* crtp) const {
    DAGEndpoints result;
    result.SetParent(crtp, GetParent(crtp));
    result.SetChild(crtp, GetChild(crtp));
    result.SetClade(crtp, GetClade(crtp));
    return result;
  }

  template <typename CRTP>
  NodeId GetParent(const CRTP* crtp) const {
    auto [dag_edge, mat, mat_node, is_ua] = access_const(crtp);
    if (is_ua) {
      return dag_edge.GetDAG().GetRoot().GetId();
    }
    // if it's an edge above a condensed node
    if (mat_node == nullptr) {
      [[maybe_unused]] auto& storage =
          dag_edge.template GetFeatureExtraStorage<MATEdgeStorage>().get();
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

  template <typename CRTP>
  NodeId GetChild(const CRTP* crtp) const {
    auto [dag_edge, mat, mat_node, is_ua] = access_const(crtp);
    return NodeId{dag_edge.GetId().value};
  }

  template <typename CRTP>
  CladeIdx GetClade(const CRTP* crtp) const {
    auto [dag_edge, mat, mat_node, is_ua] = access_const(crtp);
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
      auto& storage = dag_edge.template GetFeatureExtraStorage<MATEdgeStorage>().get();
      if (mat_node == nullptr) {
#ifdef KEEP_ASSERTS
        auto cn_id_str = storage.node_id_to_sampleid_map_.find(dag_edge.GetChildId());
        Assert(cn_id_str != storage.node_id_to_sampleid_map_.end());
#endif
        auto condensed_mat_node_str =
            storage.node_id_to_sampleid_map_.at(dag_edge.GetChildId());
#ifdef KEEP_ASSERTS
        auto condensed_mat_node_iter =
            storage.reversed_condensed_nodes_.find(condensed_mat_node_str);
        Assert(condensed_mat_node_iter != storage.reversed_condensed_nodes_.end());
#endif
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

  template <typename CRTP>
  void SetParent(const CRTP*, NodeId) {
    Fail("Can't modify MATEndpoints");
  }
  template <typename CRTP>
  void SetChild(const CRTP*, NodeId) {
    Fail("Can't modify MATEndpoints");
  }
  template <typename CRTP>
  void SetClade(const CRTP*, CladeIdx) {
    Fail("Can't modify MATEndpoints");
  }

 private:
  template <typename CRTP>
  auto access_const(const CRTP* crtp) const {
    return static_cast<const FeatureConstView<MATEdgeStorage, CRTP, MATEdgeStorage>&>(
               *crtp)
        .access();
  }

  template <typename CRTP>
  auto access_mutable(const CRTP* crtp) {
    return static_cast<const FeatureMutableView<MATEdgeStorage, CRTP, MATEdgeStorage>&>(
               *crtp)
        .access();
  }
};

template <>
struct OverlayFeatureType<MATEndpoints> {
  using store_type = DAGEndpoints;
  using const_view_type = std::variant<std::reference_wrapper<const DAGEndpoints>,
                                       std::reference_wrapper<const MATEndpoints>>;
  using mutable_view_type = std::variant<std::reference_wrapper<DAGEndpoints>,
                                         std::reference_wrapper<MATEndpoints>>;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<MATEndpoints, CRTP, Tag>
    : FeatureConstView<Endpoints, CRTP, Tag> {};

template <typename CRTP, typename Tag>
struct FeatureMutableView<MATEndpoints, CRTP, Tag>
    : FeatureMutableView<Endpoints, CRTP, Tag> {};

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
  size_t max_id_ = 0;
};

template <typename CRTP>
struct ExtraFeatureConstView<MATNodeStorage, CRTP> {
  const MAT::Tree& GetMAT() const {
    auto& dag = static_cast<const CRTP&>(*this);
    auto* mat = dag.template GetFeatureExtraStorage<Component::Node, MATNodeStorage>()
                    .get()
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

  auto GetUncondensedNodeFromMAT(MATNodePtr node) const {
    auto& dag = static_cast<const CRTP&>(*this);
    auto dag_node = GetNodeFromMAT(node);
    auto& storage =
        dag.template GetFeatureExtraStorage<Component::Node, MATNodeStorage>().get();

    std::vector<NodeId> to_ret;
    auto it = storage.condensed_nodes_.find(dag_node.GetId());
    if (it != storage.condensed_nodes_.end()) {
      to_ret = ranges::to_vector(
          it->second | ranges::views::transform([&](const std::string& x) -> NodeId {
            return NodeId{storage.sampleid_to_mat_node_id_map_.at(x)};
          }));
    } else {
      to_ret.push_back(dag_node.GetId());
    }
    return to_ret;
  }
};

template <typename CRTP>
struct ExtraFeatureMutableView<MATNodeStorage, CRTP> {
  void SetMAT(MAT::Tree* mat) const {
    Assert(mat != nullptr);
    auto& dag = static_cast<const CRTP&>(*this);
    auto& node_storage =
        dag.template GetFeatureExtraStorage<Component::Node, MATNodeStorage>().get();
    auto& edge_storage =
        dag.template GetFeatureExtraStorage<Component::Edge, MATEdgeStorage>().get();

    node_storage.mat_tree_ = mat;
    edge_storage.mat_tree_ = mat;
    size_t num_nodes = mat->get_size_upper();
    auto& cn = node_storage.condensed_nodes_;

    size_t last_uncondensed_node_id = 1;
    for (size_t ict = 1; ict < num_nodes + 1; ict++) {
      if (mat->get_node(ict) != nullptr) {
        last_uncondensed_node_id = ict + 1;
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

    size_t max_id = 0;
    for (auto* i : mat->depth_first_expansion()) {
      if (i->node_id > max_id) {
        max_id = i->node_id;
      }
    }
    if (max_id < last_uncondensed_node_id) {
      max_id = last_uncondensed_node_id;
    }
    node_storage.max_id_ = max_id;
    edge_storage.max_id_ = max_id;

    size_t ua_node_id = 0;
    mat->fix_node_idx();
    if (mat->get_node(0) != nullptr) {
      size_t count = mat->get_node_idx() + 1;
      count -= node_storage.condensed_nodes_.size();
      count += node_storage.condensed_nodes_count_;
      ua_node_id = count;
    }
    if (max_id > ua_node_id) {
      ua_node_id = max_id + 1;
    }

    Assert(mat->get_node(ua_node_id) == nullptr);

    node_storage.ua_node_id_ = NodeId{ua_node_id};
    edge_storage.ua_node_id_ = NodeId{ua_node_id};
  }
};

template <typename CRTP, typename Tag>
struct FeatureConstView<MATNodeStorage, CRTP, Tag> {
  MAT::Node* GetMATNode() const {
    auto [dag_node, mat, mat_node, is_ua] = access();
    return mat_node;
  }

  bool HaveMATNode() const { return GetMATNode() != nullptr; }

 private:
  friend struct MATNeighbors;

  static inline std::vector<MAT::Node*> empty_node{nullptr};

  NodeId GetUA() const {
    auto dag_node = static_cast<const CRTP&>(*this);
    auto dag = dag_node.GetDAG();
    return dag.template GetFeatureExtraStorage<Component::Node, MATNodeStorage>()
        .get()
        .ua_node_id_;
  }

  auto access() const {
    auto dag_node = static_cast<const CRTP&>(*this);
    NodeId id = dag_node.GetId();
    auto dag = dag_node.GetDAG();
    auto& mat = *dag.template GetFeatureExtraStorage<Component::Node, MATNodeStorage>()
                     .get()
                     .mat_tree_;
    auto* mat_node = mat.get_node(id.value);
    bool is_ua = dag.template GetFeatureExtraStorage<Component::Node, MATNodeStorage>()
                     .get()
                     .ua_node_id_ == id;
    return std::make_tuple(dag_node, std::ref(mat), mat_node, is_ua);
  }

 private:
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<MATNodeStorage, CRTP, Tag> {};

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
  size_t max_id_ = 0;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<MATEdgeStorage, CRTP, Tag> {
 private:
  friend struct MATEndpoints;

  auto access() const {
    auto& dag_edge = static_cast<const CRTP&>(*this);
    EdgeId id = dag_edge.GetId();
    auto dag = dag_edge.GetDAG();
    auto& mat = *dag.template GetFeatureExtraStorage<Component::Edge, MATEdgeStorage>()
                     .get()
                     .mat_tree_;
    // Assert((id.value == dag_edge.template GetFeatureExtraStorage<MATEdgeStorage>()
    //                         .get()
    //                         .ua_node_id_.value) or
    //        (id.value < mat.get_size_upper()));
    MAT::Node* mat_node = mat.get_node(id.value);

    bool is_ua =
        mat_node != nullptr ? ((mat.root->node_id == mat_node->node_id)) : false;
    if ((not is_ua) and (mat_node == nullptr)) {
      auto& storage = dag_edge.template GetFeatureExtraStorage<MATEdgeStorage>().get();
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
struct FeatureMutableView<MATEdgeStorage, CRTP, Tag> {};

/**
 * @brief Base container for managing MAT elements (nodes or edges) in a DAG view.
 *
 * MATElementsContainerBase provides the core infrastructure for storing and accessing
 * MAT nodes or edges within a DAG view system. It manages the mapping between MAT
 * element IDs and their corresponding DAG representations, handling both condensed and
 * uncondensed views. The container supports efficient lookup, iteration, and feature
 * storage for MAT elements, with special handling for condensed nodes that represent
 * multiple samples and the universal ancestor (UA) node.
 *
 * @tparam C The component type (Node or Edge)
 */
template <Component C, bool>
struct MATElementsContainerBase {
  using ElementStorageT =
      std::conditional_t<C == Component::Node, MATNodeStorage, MATEdgeStorage>;
  using FeatureTypes =
      std::conditional_t<C == Component::Node,
                         std::tuple<ElementStorageT, MATNeighbors>,
                         std::tuple<ElementStorageT, MATEndpoints, EdgeMutations>>;
  using AllFeatureTypes = FeatureTypes;

  static constexpr IdContinuity id_continuity = IdContinuity::Sparse;

  template <typename Feature>
  static const bool contains_element_feature =
      tuple_contains_v<AllFeatureTypes, Feature, FeatureEquivalent>;

  template <typename CRTP>
  struct ConstElementViewBase : ToFeatureConstBase<CRTP, AllFeatureTypes> {};
  template <typename CRTP>
  struct MutableElementViewBase : ToFeatureMutableBase<CRTP, AllFeatureTypes> {
    using FeatureMutableView<ElementStorageT, CRTP>::operator=;
  };

  template <typename CRTP>
  struct ExtraConstElementViewBase : ExtraFeatureConstView<ElementStorageT, CRTP> {};
  template <typename CRTP>
  struct ExtraMutableElementViewBase : ExtraFeatureMutableView<ElementStorageT, CRTP> {
  };

  MOVE_ONLY_DEF_CTOR(MATElementsContainerBase);

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
    std::size_t result = GetCount<VT>();
    while (ContainsId<VT>({result})) {  // TODO linear search
      ++result;
    }
    return {result};
  }

  template <typename VT>
  bool ContainsId(Id<C> id) const {
    if (id.value == extra_storage_.ua_node_id_.value) {
      if constexpr (C == Component::Node) {
        return true;
      } else {
        return false;
      }
    }
    if constexpr (not CheckIsCondensed<VT>::value) {
      if (GetMAT().get_node(id.value) != nullptr) {
        return true;
      }
      if (extra_storage_.node_id_to_sampleid_map_.find(NodeId{id.value}) !=
          extra_storage_.node_id_to_sampleid_map_.end()) {
        return true;
      }
    }
    return GetMAT().get_node(id.value) != nullptr;
  }

  template <typename Feature, typename E>
  auto GetFeatureStorage(Id<C> id, E elem) const {
    auto& feat_storage = [this](Id<C> i) -> auto& {
      return features_storage_[i.value];
    }(id);
    if constexpr (std::is_same_v<Feature, EdgeMutations>) {
      EdgeMutations& storage = std::get<EdgeMutations>(feat_storage);
      if (storage.empty()) {
        MAT::Node* mat_node = elem.GetChild().GetMATNode();
        // check if the edge and its child node are uncondensed from a condensed node
        // in the MAT, since in this case we need to extract the condensed node/edge's
        // mutations.
        if (mat_node == nullptr) {
          if (extra_storage_.node_id_to_sampleid_map_.find(NodeId{id.value}) !=
              extra_storage_.node_id_to_sampleid_map_.end()) {
            // find the ID of the condensed node in the MAT
            auto condensed_mat_node_str =
                extra_storage_.node_id_to_sampleid_map_.at(NodeId{id.value});
            // set mat_node to point to the condensed node
            mat_node =
                extra_storage_.reversed_condensed_nodes_.at(condensed_mat_node_str);
          }
        }
        storage = EdgeMutations{
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
      return std::cref(storage);
    } else {
      return std::cref(tuple_get<Feature, FeatureEquivalent>(feat_storage));
    }
  }

  template <typename Feature, typename E>
  auto GetFeatureStorage(Id<C> id, E elem) {
    if constexpr (std::is_same_v<Feature, EdgeMutations>) {
      static_cast<const MATElementsContainerBase*>(this)
          ->template GetFeatureStorage<Feature>(id, elem);
    }
    // TODO static_assert(not std::is_same_v<Feature, EdgeMutations>);
    return std::ref(tuple_get<Feature, FeatureEquivalent>(features_storage_[id.value]));
  }

  template <typename Feature>
  auto GetFeatureExtraStorage() {
    static_assert(std::is_same_v<Feature, ElementStorageT>);
    return std::ref(extra_storage_);
  }

  template <typename Feature>
  auto GetFeatureExtraStorage() const {
    static_assert(std::is_same_v<Feature, ElementStorageT>);
    return std::cref(extra_storage_);
  }

  template <typename VT>
  auto All() const {
    size_t iota_max = extra_storage_.max_id_ + 2;
    // getcount returns the number of nodes/edges in the tree.
    // if view is condensed, we need to check all possible node ids,
    // which means checking the max number of nodes, condensed or not.

    static constexpr const bool is_condensed = CheckIsCondensed<VT>::value;
    if constexpr (is_condensed) {
      iota_max -= extra_storage_.condensed_nodes_.size();
      iota_max += extra_storage_.condensed_nodes_count_;
    }
    if constexpr (C == Component::Node) {
      return ranges::views::iota(size_t{0}, iota_max + 2) |
             ranges::views::transform([iota_max, this](size_t i) {
               if (i > iota_max + 1) {
                 if (extra_storage_.ua_node_id_.value > iota_max + 1) {
                   return extra_storage_.ua_node_id_.value;
                 }
                 return iota_max + 2;
               }
               return i;
             }) |
             ranges::views::filter([this](size_t i) {
               if constexpr (is_condensed) {
                 return GetMAT().get_node(i) != nullptr or
                        i == extra_storage_.ua_node_id_.value;
               }
               bool result = GetMAT().get_node(i) != nullptr or
                             i == extra_storage_.ua_node_id_.value or
                             (extra_storage_.node_id_to_sampleid_map_.find(NodeId{i}) !=
                              extra_storage_.node_id_to_sampleid_map_.end());
               return result;
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
  template <typename, typename>
  friend struct ExtraFeatureMutableView;
  MAT::Tree& GetMAT() {
    Assert(extra_storage_.mat_tree_ != nullptr);
    return *extra_storage_.mat_tree_;
  }

  const MAT::Tree& GetMAT() const {
    Assert(extra_storage_.mat_tree_ != nullptr);
    return *extra_storage_.mat_tree_;
  }

  mutable ConcurrentSparseIdMap<AllFeatureTypes> features_storage_;
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
