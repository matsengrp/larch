#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

#include "larch/debug.hpp"

/**
 * Basic per-node feature.
 */
struct Neighbors {
  MOVE_ONLY(Neighbors);
  Neighbors() = default;
};

template <typename CRTP>
struct FeatureConstView<Neighbors, CRTP, Neighbors> {
  static_assert(sizeof(CRTP) == NoId, "Neighbors is an abstract feature");
};

template <typename CRTP>
struct FeatureMutableView<Neighbors, CRTP, Neighbors> {
  static_assert(sizeof(CRTP) == NoId, "Neighbors is an abstract feature");
};

struct DAGNeighbors : Neighbors {
  template <typename CRTP>
  inline DAGNeighbors Copy(const CRTP*) const {
    DAGNeighbors result;
    result.parents_ = parents_;
    result.clades_ = clades_;
    result.leafs_below_ = leafs_below_;
    return result;
  }

  template <typename CRTP>
  auto GetParents(const CRTP*) const {
    return parents_ | ranges::view::all;
  }
  template <typename CRTP>
  auto GetClades(const CRTP*) const {
    return clades_ |
           ranges::views::transform([](auto& i) { return i | ranges::view::all; });
  }
  template <typename CRTP>
  auto GetLeafsBelow(const CRTP*) const {
    return leafs_below_ |
           ranges::views::transform([](auto& i) { return i | ranges::view::all; });
  }

  template <typename CRTP>
  auto& GetParentsMutable(const CRTP*) {
    return parents_;
  }
  template <typename CRTP>
  auto& GetCladesMutable(const CRTP*) {
    return clades_;
  }
  template <typename CRTP>
  auto& GetLeafsBelowMutable(const CRTP*) {
    return leafs_below_;
  }

 private:
  std::vector<EdgeId> parents_;
  std::vector<std::vector<EdgeId>> clades_;
  std::vector<std::vector<NodeId>> leafs_below_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<Neighbors, CRTP, Tag> {
  static_assert(std::is_base_of_v<Neighbors, std::decay_t<Tag>>);
  auto GetParents() const;
  auto GetClades() const;
  auto GetClade(CladeIdx clade) const;
  size_t GetParentsCount() const;
  size_t GetCladesCount() const;
  auto GetChildren() const;
  auto GetSingleParent() const;
  auto GetFirstParent() const;
  auto GetFirstChild() const;
  auto GetFirstClade() const;
  bool IsUA() const;
  bool IsTreeRoot() const;
  bool IsLeaf() const;
  auto GetLeafsBelow() const;
  void Validate(bool recursive = false, bool allow_dag = false) const;
  auto GetParentNodes() const;
  auto GetChildNodes() const;
  bool ContainsParent(NodeId node) const;
  bool ContainsChild(NodeId node) const;
  std::string ParentsToString() const;
  std::string ChildrenToString() const;

 private:
  auto GetStorageParents() const {
    auto storage = GetFeatureStorage(this);
    auto* self = static_cast<const CRTP*>(this);
    if constexpr (is_variant_v<decltype(storage)>) {
      using Var =
          variant_of_views<decltype(std::get<0>(storage).get().GetParents(self)),
                           decltype(std::get<1>(storage).get().GetParents(self))>;
      return std::visit([self](auto& x) { return Var{x.get().GetParents(self)}; },
                        storage);
    } else {
      return storage.get().GetParents(self);
    }
  }

  auto GetStorageClades() const {
    auto storage = GetFeatureStorage(this);
    auto* self = static_cast<const CRTP*>(this);
    if constexpr (is_variant_v<decltype(storage)>) {
      using Var =
          variant_of_views<decltype(std::get<0>(storage).get().GetClades(self)),
                           decltype(std::get<1>(storage).get().GetClades(self))>;
      return std::visit([self](auto& x) { return Var{x.get().GetClades(self)}; },
                        storage);
    } else {
      return storage.get().GetClades(self);
    }
  }

  auto GetStorageLeafsBelow() const {
    auto storage = GetFeatureStorage(this);
    auto* self = static_cast<const CRTP*>(this);
    if constexpr (is_variant_v<decltype(storage)>) {
      using Var =
          variant_of_views<decltype(std::get<0>(storage).get().GetLeafsBelow(self)),
                           decltype(std::get<1>(storage).get().GetLeafsBelow(self))>;
      return std::visit([self](auto& x) { return Var{x.get().GetLeafsBelow(self)}; },
                        storage);
    } else {
      return storage.get().GetLeafsBelow(self);
    }
  }
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<Neighbors, CRTP, Tag> {
  static_assert(std::is_base_of_v<Neighbors, std::decay_t<Tag>>);
  void ClearConnections() const;
  void AddEdge(CladeIdx clade, EdgeId id, bool this_node_is_parent) const;
  void RemoveParent(EdgeId edge) const;
  void ChangeParent(EdgeId from, EdgeId to) const;
  void SetSingleParent(EdgeId parent) const;
  void RemoveChild(CladeIdx clade, EdgeId child) const;
  void ChangeChild(CladeIdx clade, EdgeId from, EdgeId to) const;
  void CalculateLeafsBelow() const;

 private:
  auto& GetStorage() const {
    auto storage = GetFeatureStorage(this);
    if constexpr (is_variant_v<decltype(storage)>) {
      if (not std::holds_alternative<std::reference_wrapper<DAGNeighbors>>(storage)) {
        Fail("Only DAGNeighbors can be modified");
      }
      return std::get<std::reference_wrapper<DAGNeighbors>>(storage).get();
    } else {
      if constexpr (not std::is_same_v<decltype(storage),
                                       std::reference_wrapper<DAGNeighbors>>) {
        Fail("Only DAGNeighbors can be modified");
      } else {
        return storage.get();
      }
    }
  }
};

template <typename CRTP, typename Tag>
struct FeatureConstView<DAGNeighbors, CRTP, Tag>
    : FeatureConstView<Neighbors, CRTP, Tag> {};

template <typename CRTP, typename Tag>
struct FeatureMutableView<DAGNeighbors, CRTP, Tag>
    : FeatureMutableView<Neighbors, CRTP, Tag> {};
