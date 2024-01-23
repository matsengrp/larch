#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

#include <atomic>
#include <iostream>

#include "larch/parallel/reduction.hpp"

template <typename CRTP, typename Tag>
bool FeatureConstView<Connections, CRTP, Tag>::IsTree() const {
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.GetNodesCount() == dag.GetEdgesCount() + 1;
}
template <typename CRTP, typename Tag>
bool FeatureConstView<Connections, CRTP, Tag>::HaveRoot() const {
  return GetFeatureStorage(this).root_.value != NoId;
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Connections, CRTP, Tag>::GetRoot() const {
  auto& dag = static_cast<const CRTP&>(*this);
  using Node = typename CRTP::NodeView;
  return Node{dag, GetFeatureStorage(this).root_};
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Connections, CRTP, Tag>::GetLeafs() const {
  auto& dag = static_cast<const CRTP&>(*this);
  return GetFeatureStorage(this).leafs_ | Transform::ToNodes(dag);
}

template <typename CRTP, typename Tag>
auto FeatureConstView<Connections, CRTP, Tag>::GetLeafsCount() const {
  return GetFeatureStorage(this).leafs_.size();
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Connections, CRTP, Tag>::BuildConnections() const {
  BuildConnectionsRaw();
  BuildRootAndLeafs();
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Connections, CRTP, Tag>::BuildConnectionsRaw() const {
  auto& dag = static_cast<const CRTP&>(*this);
  ParallelForEach(dag.GetNodes(), [](auto node) { node.ClearConnections(); });
  for (auto edge : dag.GetEdges()) {
    Assert(edge.GetParentId().value != NoId && "Edge has no parent");
    Assert(edge.GetChildId().value != NoId && "Edge has no child");
    Assert(edge.GetClade().value != NoId && "Edge has no clade index");
    Assert(edge.GetParentId() != edge.GetChildId() && "Edge is looped");
    edge.GetParent().AddEdge(edge.GetClade(), edge, true);
    edge.GetChild().AddEdge(edge.GetClade(), edge, false);
  }
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Connections, CRTP, Tag>::BuildRootAndLeafs() const {
  auto& storage = GetFeatureStorage(this);
  auto& dag = static_cast<const CRTP&>(*this);
  storage.root_ = {NoId};
  storage.leafs_ = {};
  std::atomic<size_t> root_id{NoId};
  Reduction<std::vector<NodeId>> leafs{32};
  ParallelForEach(dag.GetNodes(), [&](auto node) {
    for (auto clade : node.GetClades()) {
      Assert(not clade.empty() && "Empty clade");
    }
    if (node.IsUA()) {
      const size_t previous = root_id.exchange(node.GetId().value);
      if (previous != NoId) {
        std::cout << "Duplicate root: " << previous << " and " << node.GetId().value
                  << "\n";
      }
      Assert(previous == NoId);
    }
    if (node.IsLeaf()) {
      leafs.AddElement([](std::vector<NodeId>& ls, NodeId id) { ls.push_back(id); },
                       node);
    }
  });
  storage.root_.value = root_id.load();
  leafs.GatherAndClear(
      [&leafs](auto buckets, auto& stor) {
        for (auto&& bucket : buckets) {
          stor.leafs_.reserve(leafs.size_approx());
          stor.leafs_.insert(stor.leafs_.end(), bucket.begin(), bucket.end());
        }
      },
      storage);
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Connections, CRTP, Tag>::AddLeaf(NodeId id) const {
  auto& storage = GetFeatureStorage(this);
  // TODO make leafs ContiguousSet ?
  if (std::find(storage.leafs_.begin(), storage.leafs_.end(), id) ==
      storage.leafs_.end()) {
    storage.leafs_.push_back(id);
  }
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Connections, CRTP, Tag>::ClearConnections() const {
  auto& dag = static_cast<const CRTP&>(*this);
  for (auto node : dag.GetNodes()) {
    node.ClearConnections();
  }
  dag.GetStorage().ClearEdges();
}

template <typename CRTP, typename Tag>
std::map<std::set<NodeId>, std::set<NodeId>>
FeatureMutableView<Connections, CRTP, Tag>::BuildCladeUnionMap() const {
  auto& dag = static_cast<const CRTP&>(*this);
  std::map<std::set<NodeId>, std::set<NodeId>> clade_union_map;
  dag.GetRoot().CalculateLeafsBelow();
  for (auto node : dag.GetNodes()) {
    std::set<NodeId> full_leafset;
    if (node.GetLeafsBelow().size() > 0) {
      for (const auto clade_leafset : node.GetLeafsBelow()) {
        full_leafset.insert(clade_leafset.begin(), clade_leafset.end());
      }
    } else {
      full_leafset.insert(node.GetId());
    }
    if (clade_union_map.find(full_leafset) == clade_union_map.end()) {
      // TODO Copy()
      clade_union_map[full_leafset] = std::set<NodeId>();
    }
    clade_union_map[full_leafset].insert(node.GetId());
  }
  return clade_union_map;
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Connections, CRTP, Tag>::MakeComplete() const {
  auto& dag = static_cast<const CRTP&>(*this);
  auto clade_union_map = dag.BuildCladeUnionMap();
  size_t taxon_count = dag.GetLeafsCount();
  dag.ClearConnections();
  // Connect rootsplit nodes.
  bool rootsplits_found = false;
  for (auto& [clade_union, node_ids] : clade_union_map) {
    if (clade_union.size() == taxon_count) {
      for (auto node_id : node_ids) {
        dag.AppendEdge(dag.GetRoot().GetId(), node_id, {0});
      }
      rootsplits_found = true;
      break;
    }
  }
  Assert(rootsplits_found);
  // Connect all other nodes.
  for (auto parent_node : dag.GetNodes()) {
    auto leaf_sets = parent_node.GetLeafsBelow();
    for (size_t clade_idx = 0; clade_idx < leaf_sets.size(); clade_idx++) {
      auto leaf_clade = leaf_sets[clade_idx];
      std::set<NodeId> clade_set(leaf_clade.begin(), leaf_clade.end());
      auto possible_children = clade_union_map.find(clade_set);
      if (possible_children != clade_union_map.end()) {
        for (auto child_node_id : possible_children->second) {
          if (parent_node.GetId() == child_node_id) continue;
          dag.AppendEdge(parent_node.GetId(), child_node_id, {clade_idx});
        }
      }
    }
  }
}
