#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename CRTP, typename Tag>
auto FeatureConstView<Neighbors, CRTP, Tag>::GetParents() const {
  auto dag = static_cast<const CRTP&>(*this).GetDAG();
  return GetFeatureStorage(this).parents_ | Transform::ToEdges(dag);
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Neighbors, CRTP, Tag>::GetClades() const {
  auto dag = static_cast<const CRTP&>(*this).GetDAG();
  return GetFeatureStorage(this).clades_ |
         ranges::views::transform([dag](const std::vector<EdgeId>& clade) {
           return clade | Transform::ToEdges(dag);
         });
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Neighbors, CRTP, Tag>::GetClade(CladeIdx clade) const {
  auto dag = static_cast<const CRTP&>(*this).GetDAG();
  return GetFeatureStorage(this).clades_.at(clade.value) | Transform::ToEdges(dag);
}

template <typename CRTP, typename Tag>
size_t FeatureConstView<Neighbors, CRTP, Tag>::GetParentsCount() const {
  return GetFeatureStorage(this).parents_.size();
}

template <typename CRTP, typename Tag>
size_t FeatureConstView<Neighbors, CRTP, Tag>::GetCladesCount() const {
  return GetFeatureStorage(this).clades_.size();
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Neighbors, CRTP, Tag>::GetChildren() const {
  return GetClades() | ranges::views::join;
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Neighbors, CRTP, Tag>::GetSingleParent() const {
  Assert(GetParentsCount() == 1);
  return *GetParents().begin();
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Neighbors, CRTP, Tag>::GetFirstParent() const {
  Assert(not IsRoot());
  return *GetParents().begin();
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Neighbors, CRTP, Tag>::GetFirstChild() const {
  Assert(not IsLeaf());
  return *GetChildren().begin();
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Neighbors, CRTP, Tag>::GetFirstClade() const {
  return GetClade({0});
}

template <typename CRTP, typename Tag>
bool FeatureConstView<Neighbors, CRTP, Tag>::IsRoot() const {
  return GetFeatureStorage(this).parents_.empty();
}

template <typename CRTP, typename Tag>
bool FeatureConstView<Neighbors, CRTP, Tag>::IsLeaf() const {
  return ranges::all_of(GetFeatureStorage(this).clades_,
                        [](const auto& clade) { return clade.empty(); });
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Neighbors, CRTP, Tag>::ClearConnections() const {
  GetFeatureStorage(this).parents_.clear();
  GetFeatureStorage(this).clades_.clear();
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Neighbors, CRTP, Tag>::AddEdge(CladeIdx clade, EdgeId id,
                                                       bool this_node_is_parent) const {
  if (this_node_is_parent) {
    GetOrInsert(GetFeatureStorage(this).clades_, clade).push_back(id);
  } else {
    GetFeatureStorage(this).parents_.push_back(id);
  }
}