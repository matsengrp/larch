#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * Basic per-edge feature.
 */
struct Endpoints {
  MOVE_ONLY(Endpoints);
  Endpoints() = default;
};

template <typename CRTP>
struct FeatureConstView<Endpoints, CRTP, Endpoints> {
  static_assert(false, "Endpoints is an abstract feature");
};

template <typename CRTP>
struct FeatureMutableView<Endpoints, CRTP, Endpoints> {
  static_assert(false, "Endpoints is an abstract feature");
};

struct DAGEndpoints : Endpoints {
  template <typename CRTP>
  inline DAGEndpoints Copy(const CRTP*) const {
    DAGEndpoints result;
    result.parent_ = parent_;
    result.child_ = child_;
    result.clade_ = clade_;
    return result;
  }

  template <typename CRTP>
  NodeId GetParent(const CRTP*) const {
    return parent_;
  }
  template <typename CRTP>
  NodeId GetChild(const CRTP*) const {
    return child_;
  }
  template <typename CRTP>
  CladeIdx GetClade(const CRTP*) const {
    return clade_;
  }

  template <typename CRTP>
  void SetParent(const CRTP*, NodeId parent) {
    parent_ = parent;
  }
  template <typename CRTP>
  void SetChild(const CRTP*, NodeId child) {
    child_ = child;
  }
  template <typename CRTP>
  void SetClade(const CRTP*, CladeIdx clade) {
    clade_ = clade;
  }

 private:
  NodeId parent_;
  NodeId child_;
  CladeIdx clade_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<Endpoints, CRTP, Tag> {
  static_assert(std::is_base_of_v<Endpoints, std::decay_t<Tag>>);
  auto GetParent() const;
  auto GetChild() const;
  CladeIdx GetClade() const;
  NodeId GetParentId() const;
  NodeId GetChildId() const;
  std::pair<NodeId, NodeId> GetNodeIds() const;
  bool IsUA() const;
  bool IsTreeRoot() const;
  bool IsLeaf() const;

 private:
  NodeId GetStorageParent() const {
    auto storage = GetFeatureStorage(this);
    auto* self = static_cast<const CRTP*>(this);
    if constexpr (is_variant_v<decltype(storage)>) {
      return std::visit([self](auto& x) { return x.get().GetParent(self); }, storage);
    } else {
      return storage.get().GetParent(self);
    }
  }

  NodeId GetStorageChild() const {
    auto storage = GetFeatureStorage(this);
    auto* self = static_cast<const CRTP*>(this);
    if constexpr (is_variant_v<decltype(storage)>) {
      return std::visit([self](auto& x) { return x.get().GetChild(self); }, storage);
    } else {
      return storage.get().GetChild(self);
    }
  }

  CladeIdx GetStorageClade() const {
    auto storage = GetFeatureStorage(this);
    auto* self = static_cast<const CRTP*>(this);
    if constexpr (is_variant_v<decltype(storage)>) {
      return std::visit([self](auto& x) { return x.get().GetClade(self); }, storage);
    } else {
      return storage.get().GetClade(self);
    }
  }
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<Endpoints, CRTP, Tag> {
  static_assert(std::is_base_of_v<Endpoints, std::decay_t<Tag>>);
  void Set(NodeId parent, NodeId child, CladeIdx clade) const;
  void SetParent(NodeId parent) const;
  void SetChild(NodeId child) const;
  void SetClade(CladeIdx clade) const;

 private:
  auto& GetStorage() const {
    auto storage = GetFeatureStorage(this);
    if constexpr (is_variant_v<decltype(storage)>) {
      if (not std::holds_alternative<std::reference_wrapper<DAGEndpoints>>(storage)) {
        Fail("Only DAGEndpoints can be modified");
      }
      return std::get<std::reference_wrapper<DAGEndpoints>>(storage).get();
    } else {
      if constexpr (not std::is_same_v<decltype(storage),
                                       std::reference_wrapper<DAGEndpoints>>) {
        Fail("Only DAGEndpoints can be modified");
      } else {
        return GetFeatureStorage(this).get();
      }
    }
  }
};

template <typename CRTP, typename Tag>
struct FeatureConstView<DAGEndpoints, CRTP, Tag>
    : FeatureConstView<Endpoints, CRTP, Tag> {};

template <typename CRTP, typename Tag>
struct FeatureMutableView<DAGEndpoints, CRTP, Tag>
    : FeatureMutableView<Endpoints, CRTP, Tag> {};
