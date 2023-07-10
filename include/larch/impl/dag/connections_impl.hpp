#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

#include "larch/parallel/reduction.hpp"
#include "larch/parallel/for_loop.hpp"

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
  return GetFeatureStorage(this).leafs_ |
         ranges::views::transform([dag, idx = size_t{}](auto&) mutable {
           using Node = typename CRTP::NodeView;
           return Node{dag, {idx++}};
         });
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Connections, CRTP, Tag>::BuildConnections() const {
  auto& storage = GetFeatureStorage(this);
  auto& dag = static_cast<const CRTP&>(*this);
  storage.root_ = {NoId};
  storage.leafs_ = {};
  BuildConnectionsRaw();
  std::atomic<size_t> root_id{NoId};
  Reduction<NodeId> leafs{DefaultScheduler().WorkersCount()};
  seq_for_each(dag.Const().GetNodesCount(), [&](size_t i, size_t worker) {
    auto node = dag.Get(NodeId{i});
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
      leafs.Emplace(worker, node);
    }
  });
  storage.root_.value = root_id.load();
  Assert(storage.root_.value != NoId);
  leafs.Consume([&](auto all_leafs) {
    storage.leafs_.insert(storage.leafs_.end(), all_leafs.begin(), all_leafs.end());
  });
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Connections, CRTP, Tag>::BuildConnectionsRaw() const {
  auto& dag = static_cast<const CRTP&>(*this);
  seq_for_each(dag.GetNodesCount(),
               [&](size_t i, size_t) { dag.Get(NodeId{i}).ClearConnections(); });
  for (auto edge : dag.GetEdges()) {
    Assert(edge.GetParentId().value != NoId && "Edge has no parent");
    Assert(edge.GetChildId().value != NoId && "Edge has no child");
    Assert(edge.GetClade().value != NoId && "Edge has no clade index");
    Assert(edge.GetParentId() != edge.GetChildId() && "Edge is looped");
    edge.GetParent().AddEdge(edge.GetClade(), edge, true);
    edge.GetChild().AddEdge(edge.GetClade(), edge, false);
  };
}
