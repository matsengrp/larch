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
auto FeatureConstView<Neighbors, CRTP, Tag>::GetLeafsBelow() const {
  auto dag = static_cast<const CRTP&>(*this).GetDAG();
  return GetFeatureStorage(this).leafs_below_ |
         ranges::views::transform([dag](const std::vector<NodeId>& i) {
           return i | Transform::ToNodes(dag);
         });
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

template <typename CRTP, typename Tag>
void FeatureMutableView<Neighbors, CRTP, Tag>::RemoveParent(EdgeId edge) const {
  auto& parents = GetFeatureStorage(this).parents_;
  auto it = ranges::find(parents, edge);
  Assert(it != parents.end());
  parents.erase(it);
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Neighbors, CRTP, Tag>::RemoveChild(CladeIdx clade,
                                                           EdgeId child) const {
  auto node = static_cast<const CRTP&>(*this);
  auto& clades = GetFeatureStorage(this).clades_;
  auto& children = clades.at(clade.value);
  auto it = ranges::find(children, child);
  Assert(it != children.end());
  children.erase(it);
  if (children.empty()) {
    clades.erase(clades.begin() + clade.value);
    for (size_t i = clade.value; i < clades.size(); ++i) {
      for (EdgeId edge : clades.at(i)) {
        node.GetDAG().Get(edge).SetClade({i});
      }
    }
  }
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Neighbors, CRTP, Tag>::CalculateLeafsBelow() const {
  auto& self = static_cast<const CRTP&>(*this);
  std::vector<std::vector<NodeId>> result;
  result.reserve(self.GetCladesCount());
  for (auto clade : self.GetClades()) {
    std::vector<NodeId> clade_leafs;
    clade_leafs.reserve(clade.size());
    for (auto child : clade | Transform::GetChild()) {
      if (child.IsLeaf()) {
        clade_leafs.push_back(child);
      } else {
        child.CalculateLeafsBelow();
        for (auto i : child.GetLeafsBelow() | ranges::views::join) {
          clade_leafs.push_back(i);
        }
      }
    }
    clade_leafs |= ranges::actions::sort(
        [](NodeId lhs, NodeId rhs) { return lhs.value < rhs.value; });
    result.push_back(std::move(clade_leafs));
  }
  GetFeatureStorage(this).leafs_below_ = std::move(result);
}
