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
  auto& children = GetFeatureStorage(this).clades_.at(clade.value);
  auto it = ranges::find(children, child);
  Assert(it != children.end());
  children.erase(it);
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Neighbors, CRTP, Tag>::GetParentNodes() const {
  return GetParents() | Transform::GetParent();
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Neighbors, CRTP, Tag>::GetChildNodes() const {
  return GetChildren() | Transform::GetChild();
}

template <typename CRTP, typename Tag>
bool FeatureConstView<Neighbors, CRTP, Tag>::ContainsParent(NodeId node) const {
  return ranges::any_of(GetParents(), [node](auto parent_edge) {
    return parent_edge.GetParent().GetId() == node;
  });
}

template <typename CRTP, typename Tag>
bool FeatureConstView<Neighbors, CRTP, Tag>::ContainsChild(NodeId node) const {
  return ranges::any_of(GetChildren(), [node](auto child_edge) {
    return child_edge.GetChild().GetId() == node;
  });
}

template <typename CRTP, typename Tag>
std::string FeatureConstView<Neighbors, CRTP, Tag>::ParentsToString() const {
  std::stringstream os;
  os << "[";
  for (auto parent : GetParentNodes()) {
    os << parent.GetId().value << ", ";
  }
  os << "]";
  return os.str();
}

template <typename CRTP, typename Tag>
std::string FeatureConstView<Neighbors, CRTP, Tag>::ChildrenToString() const {
  std::stringstream os;
  os << "[";
  for (auto child : GetChildNodes()) {
    os << child.GetId().value << ", ";
  }
  os << "]";
  return os.str();
}
