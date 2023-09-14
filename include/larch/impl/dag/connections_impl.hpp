#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

#include <atomic>
#include <tbb/parallel_for_each.h>
#include <tbb/concurrent_vector.h>
#include <iostream>

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
void FeatureMutableView<Connections, CRTP, Tag>::BuildConnections() const {
  auto& storage = GetFeatureStorage(this);
  auto& dag = static_cast<const CRTP&>(*this);
  storage.root_ = {NoId};
  storage.leafs_ = {};
  BuildConnectionsRaw();
  std::atomic<size_t> root_id{NoId};
  tbb::concurrent_vector<NodeId> leafs;
  tbb::parallel_for_each(dag.GetNodes(), [&](auto node) {
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
      leafs.push_back(node);
    }
  });
  storage.root_.value = root_id.load();
  storage.leafs_.insert(storage.leafs_.end(), leafs.begin(), leafs.end());
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Connections, CRTP, Tag>::BuildConnectionsRaw() const {
  auto& dag = static_cast<const CRTP&>(*this);
  tbb::parallel_for_each(dag.GetNodes(), [](auto node) { node.ClearConnections(); });
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
void FeatureMutableView<Connections, CRTP, Tag>::AddLeaf(NodeId id) const {
  auto& storage = GetFeatureStorage(this);
  storage.leafs_.push_back(id);
}
