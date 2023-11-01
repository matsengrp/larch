#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

#include <set>
#include "larch/contiguous_set.hpp"

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
  Assert(not IsUA());
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
bool FeatureConstView<Neighbors, CRTP, Tag>::IsUA() const {
  return GetFeatureStorage(this).parents_.empty();
}

template <typename CRTP, typename Tag>
bool FeatureConstView<Neighbors, CRTP, Tag>::IsTreeRoot() const {
  if (GetParentsCount() != 1) {
    return false;
  }
  return GetSingleParent().IsUA();
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

namespace {

template <typename Edge>
void CheckCycle(Edge edge, std::vector<std::pair<bool, bool>>& visited_finished) {
  auto& [visited, finished] = visited_finished.at(edge.GetId().value);
  if (finished) {
    return;
  }
  if (visited) {
    Assert(false && "Cycle detected");
  }
  visited = true;
  for (auto child : edge.GetChild().GetChildren()) {
    CheckCycle(child, visited_finished);
  }
  finished = true;
}

}  // namespace

struct SampleId;

template <typename CRTP, typename Tag>
void FeatureConstView<Neighbors, CRTP, Tag>::Validate(bool recursive,
                                                      bool allow_dag) const {
  auto node = static_cast<const CRTP&>(*this).Const();
  auto dag = node.GetDAG();
  auto& storage = GetFeatureStorage(this);
  if (node.IsUA()) {
    Assert(dag.HaveUA());
    Assert(node.GetId() == dag.GetRoot());
  } else {
    size_t children_count = 0;
    for ([[maybe_unused]] auto child : node.GetChildren()) {
      ++children_count;
    }
    Assert(children_count != 1);
  }
  std::set<std::string> sample_ids;
  if (node.IsLeaf()) {
    Assert(storage.clades_.empty());
    if constexpr (std::remove_reference_t<decltype(node)>::template contains_feature<
                      SampleId> or
                  std::remove_reference_t<decltype(node)>::template contains_feature<
                      Deduplicate<SampleId>>) {
      Assert(node.HaveSampleId());
      Assert(sample_ids.insert(node.GetSampleId().value()).second);
    }
  }
  for (auto& i : storage.clades_) {
    Assert(not i.empty());
  }
  if (not allow_dag) {
    if (storage.parents_.size() > 1) {
      throw std::runtime_error{std::string{"Mulptiple parents at node "} +
                               std::to_string(node.GetId().value)};
    }
  }
  if (not storage.parents_.empty()) {
    auto edge = dag.Get(storage.parents_.at(0));
    if (edge.GetChild().GetId() != node.GetId()) {
      throw std::runtime_error{std::string{"Mismatch parent at node "} +
                               std::to_string(node.GetId().value) + " : " +
                               std::to_string(edge.GetChildId().value) +
                               ", should be " + std::to_string(node.GetId().value)};
    }
  }
  for (CladeIdx i{0}; i.value < storage.clades_.size(); ++i.value) {
    auto& clade = storage.clades_.at(i.value);
    if (not allow_dag) {
      if (clade.size() != 1) {
        std::string children;
        for (auto j : clade) {
          children += std::to_string(dag.Get(j).GetChild().GetId().value) + ", ";
        }
        throw std::runtime_error{std::string{"Mulptiple children at node "} +
                                 std::to_string(node.GetId().value) + ", clade " +
                                 std::to_string(i.value) + " : " + children};
      }
    }
    auto edge = dag.Get(clade.at(0));
    if (edge.GetParent().GetId() != node.GetId()) {
      throw std::runtime_error{std::string{"Mismatch child edge id at node "} +
                               std::to_string(node.GetId().value) + ", clade " +
                               std::to_string(i.value) + " : " +
                               std::to_string(edge.GetParent().GetId().value) +
                               ", should be " + std::to_string(node.GetId().value)};
    }
    if (edge.GetClade() != i) {
      throw std::runtime_error{std::string{"Mismatch child edge clade at node "} +
                               std::to_string(node.GetId().value) + ", clade " +
                               std::to_string(i.value)};
    }
    if (recursive) {
      edge.GetChild().Validate(true, allow_dag);
    }
  }

  if (recursive and node.IsUA()) {
    std::vector<std::pair<bool, bool>> visited_finished;
    visited_finished.resize(dag.GetEdgesCount());
    for (auto i : node.GetChildren()) {
      CheckCycle(i, visited_finished);
    }

    for (auto i : dag.GetEdges()) {
      Assert(i.GetChild().ContainsParent(i.GetParent()));
      Assert(i.GetParent().ContainsChild(i.GetChild()));
    }

    ContiguousSet<NodeId> leafs1;
    ContiguousSet<NodeId> leafs2;
    ranges::actions::insert(leafs1, dag.GetLeafs());
    ranges::actions::insert(leafs2, dag.GetNodes() | ranges::view::filter([](auto i) {
                                      return i.IsLeaf();
                                    }));
    Assert(ranges::equal(leafs1, leafs2));
  }
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
void FeatureMutableView<Neighbors, CRTP, Tag>::ChangeParent(EdgeId from,
                                                            EdgeId to) const {
  auto& parents = GetFeatureStorage(this).parents_;
  auto it = ranges::find(parents, from);
  Assert(it != parents.end());
  *it = to;
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Neighbors, CRTP, Tag>::SetSingleParent(EdgeId parent) const {
  auto& parents = GetFeatureStorage(this).parents_;
  parents.clear();
  parents.push_back(parent);
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
    clades.erase(clades.begin() + static_cast<ssize_t>(clade.value));
    for (size_t i = clade.value; i < clades.size(); ++i) {
      for (EdgeId edge : clades.at(i)) {
        node.GetDAG().Get(edge).SetClade({i});
      }
    }
  }
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Neighbors, CRTP, Tag>::ChangeChild(CladeIdx clade, EdgeId from,
                                                           EdgeId to) const {
  auto& clades = GetFeatureStorage(this).clades_;
  auto& children = clades.at(clade.value);
  auto it = ranges::find(children, from);
  Assert(it != children.end());
  *it = to;
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
