#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

#include <iostream>
#include <set>
#include "larch/contiguous_set.hpp"
#include "larch/id_container.hpp"

template <typename CRTP, typename Tag>
auto FeatureConstView<Neighbors, CRTP, Tag>::GetParents() const {
  auto dag = static_cast<const CRTP&>(*this).GetDAG();
  return GetStorageParents() | Transform::ToEdges(dag);
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Neighbors, CRTP, Tag>::GetClades() const {
  auto dag = static_cast<const CRTP&>(*this).GetDAG();
  return GetStorageClades() | ranges::views::transform([dag](auto&& clade) {
           return clade | Transform::ToEdges(dag);
         });
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Neighbors, CRTP, Tag>::GetClade(CladeIdx clade) const {
  auto dag = static_cast<const CRTP&>(*this).GetDAG();
  return GetStorageClades().at(clade.value) | Transform::ToEdges(dag);
}

template <typename CRTP, typename Tag>
size_t FeatureConstView<Neighbors, CRTP, Tag>::GetParentsCount() const {
  return static_cast<size_t>(GetStorageParents().size());
}

template <typename CRTP, typename Tag>
size_t FeatureConstView<Neighbors, CRTP, Tag>::GetCladesCount() const {
  return static_cast<size_t>(GetStorageClades().size());
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
  return GetStorageParents().empty();
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
  Assert(ranges::all_of(GetStorageClades(),
                        [](auto&& clade) { return not clade.empty(); }));
  return GetStorageClades().empty();
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Neighbors, CRTP, Tag>::GetLeafsBelow() const {
  auto dag = static_cast<const CRTP&>(*this).GetDAG();
  return GetStorageLeafsBelow() | ranges::views::transform([dag](auto&& i) {
           return i | Transform::ToNodes(dag);
         });
}

namespace {

template <typename DAG>
void MADAGToDOTCycle(DAG dag, std::ostream& out, const std::deque<EdgeId>& path1,
                     const std::deque<EdgeId>& path2) {
  out << "digraph G {\n";
  out << "  forcelabels=true\n";
  out << "  nodesep=1.0\n";
  out << "  ranksep=2.0\n";
  out << "  ratio=1.0\n";
  out << "  node [color=azure4,fontcolor=black,penwidth=4]\n";
  out << "  edge [color=azure3,fontcolor=black,penwidth=4]\n";
  for (auto edge : dag.Const().GetEdges()) {
    const bool in_path1 = std::find(path1.begin(), path1.end(), edge) != path1.end();
    const bool in_path2 = std::find(path2.begin(), path2.end(), edge) != path2.end();
    out << "  \"" << CompactGenomeToString(edge.GetParent()) << "\" -> \""
        << CompactGenomeToString(edge.GetChild()) << "\"";
    out << "[ headlabel=\"";
    out << EdgeMutationsToString(edge);
    out << "\"";
    if (in_path1) {
      out << "color=red,penwidth=40,arrowsize=2";
    }
    if (in_path2) {
      out << "color=blue,penwidth=40,arrowsize=2";
    }
    out << "]\n";
  }
  out << "}\n";
}

template <typename Edge, typename VisitedType>
void CheckCycle(Edge edge, VisitedType& visited_finished, std::deque<EdgeId>& path) {
  auto& [visited, finished, first_path] = visited_finished[edge.GetId()];
  if (finished) {
    return;
  }
  path.push_back(edge);
  if (visited) {
    MADAGToDOTCycle(edge.GetDAG(), std::cout, path, first_path);
    Assert(false && "Cycle detected");
  }
  first_path = path;
  visited = true;
  for (auto child : edge.GetChild().GetChildren()) {
    CheckCycle(child, visited_finished, path);
  }
  path.pop_back();
  finished = true;
}

}  // namespace

struct SampleId;

template <typename CRTP, typename Tag>
void FeatureConstView<Neighbors, CRTP, Tag>::Validate(
    [[maybe_unused]] bool recursive, [[maybe_unused]] bool allow_dag) const {
#ifndef NDEBUG
  auto node = static_cast<const CRTP&>(*this).Const();
  auto dag = node.GetDAG();
  if (node.IsUA()) {
    // Assert(dag.HaveUA());
    Assert(node.GetId() == dag.GetRoot());
  } else if (allow_dag) {
    size_t children_count = 0;
    for ([[maybe_unused]] auto child : node.GetChildren()) {
      ++children_count;
    }
    Assert(children_count != 1);
  }
  ContiguousSet<std::string> sample_ids;
  if (node.IsLeaf()) {
    Assert(GetStorageClades().empty());
    using NodeT = std::remove_reference_t<decltype(node)>;
    if constexpr (NodeT::template contains_feature<SampleId> or
                  NodeT::template contains_feature<Deduplicate<SampleId>>) {
      if constexpr (not is_specialization_v<std::remove_const_t<std::remove_reference_t<
                                                decltype(dag.GetStorage())>>,
                                            FragmentStorage>) {
        Assert(node.HaveSampleId());
        Assert(sample_ids.insert(std::string{node.GetSampleId().value()}).second);
      }
    }
  }
  for (auto&& i : GetStorageClades()) {
    Assert(not i.empty());
  }
  if (not allow_dag) {
    if (GetStorageParents().size() > 1) {
      throw std::runtime_error{std::string{"Multiple parents at node "} +
                               std::to_string(node.GetId().value)};
    }
  }
  if (not GetStorageParents().empty()) {
    auto edge = dag.Get(*GetStorageParents().begin());
    if (edge.GetChild().GetId() != node.GetId()) {
      throw std::runtime_error{std::string{"Mismatch parent at node "} +
                               std::to_string(node.GetId().value) + " : " +
                               std::to_string(edge.GetChildId().value) +
                               ", should be " + std::to_string(node.GetId().value)};
    }
  } else {
    Assert(node.IsUA());
  }
  CladeIdx cidx{0};
  for (auto&& clade : GetStorageClades()) {
    CladeIdx i{cidx.value++};
    if (not allow_dag) {
      if (clade.size() != 1) {
        std::string children;
        for (auto j : clade) {
          children += std::to_string(dag.Get(j).GetChild().GetId().value) + ", ";
        }
        throw std::runtime_error{std::string{"Multiple children at node "} +
                                 std::to_string(node.GetId().value) + ", clade " +
                                 std::to_string(i.value) + " : " + children};
      }
    }
    auto edge = dag.Get(*clade.begin());
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
    if (not allow_dag) {
      Assert(dag.IsTree());
    }
    size_t node_count = 0;
    for ([[maybe_unused]] auto i : dag.GetNodes()) {
      ++node_count;
    }

    if (node_count != dag.GetNodesCount()) {
      std::cout << "node_count: " << node_count
                << "  dag.GetNodesCount(): " << dag.GetNodesCount() << "\n";
    }
    Assert(node_count == dag.GetNodesCount());
    size_t edge_count = 0;
    for ([[maybe_unused]] auto i : dag.GetEdges()) {
      ++edge_count;
    }
    if (edge_count != dag.GetEdgesCount()) {
      std::cout << "edge_count: " << edge_count
                << "  dag.GetEdgesCount(): " << dag.GetEdgesCount() << "\n"
                << std::flush;
    }
    Assert(edge_count == dag.GetEdgesCount());
    if (not allow_dag) {
      constexpr IdContinuity id_continuity =
          decltype(dag)::template id_continuity<Component::Node>;
      IdContainer<EdgeId, std::tuple<bool, bool, std::deque<EdgeId>>, id_continuity,
                  Ordering::Ordered>
          visited_finished;
      visited_finished.reserve(dag.GetEdgesCount());
      std::deque<EdgeId> path;
      for (auto i : node.GetChildren()) {
        CheckCycle(i, visited_finished, path);
      }
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
#endif
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Neighbors, CRTP, Tag>::ClearConnections() const {
  GetStorage().GetParentsMutable(static_cast<const CRTP*>(this)).clear();
  GetStorage().GetCladesMutable(static_cast<const CRTP*>(this)).clear();
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Neighbors, CRTP, Tag>::AddEdge(CladeIdx clade, EdgeId id,
                                                       bool this_node_is_parent) const {
  if (this_node_is_parent) {
    GetOrInsert(GetStorage().GetCladesMutable(static_cast<const CRTP*>(this)), clade)
        .push_back(id);
  } else {
    GetStorage().GetParentsMutable(static_cast<const CRTP*>(this)).push_back(id);
  }
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Neighbors, CRTP, Tag>::RemoveParent(EdgeId edge) const {
  auto& parents = GetStorage().GetParentsMutable(static_cast<const CRTP*>(this));
  auto it = ranges::find(parents, edge);
  Assert(it != parents.end());
  parents.erase(it);
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Neighbors, CRTP, Tag>::ChangeParent(EdgeId from,
                                                            EdgeId to) const {
  auto& parents = GetStorage().GetParentsMutable(static_cast<const CRTP*>(this));
  auto it = ranges::find(parents, from);
  Assert(it != parents.end());
  *it = to;
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Neighbors, CRTP, Tag>::SetSingleParent(EdgeId parent) const {
  auto& parents = GetStorage().GetParentsMutable(static_cast<const CRTP*>(this));
  parents.clear();
  parents.push_back(parent);
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Neighbors, CRTP, Tag>::RemoveChild(CladeIdx clade,
                                                           EdgeId child) const {
  auto node = static_cast<const CRTP&>(*this);
  auto& clades = GetStorage().GetCladesMutable(static_cast<const CRTP*>(this));
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
  auto& clades = GetStorage().GetCladesMutable(static_cast<const CRTP*>(this));
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
  GetStorage().GetLeafsBelowMutable(static_cast<const CRTP*>(this)) = std::move(result);
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
