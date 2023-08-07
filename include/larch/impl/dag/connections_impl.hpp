#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

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
  return GetFeatureStorage(this).leafs_ |
         ranges::views::transform([*this, dag, idx = size_t{}](auto&) mutable {
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
  for (auto node : dag.GetNodes()) {
    for (auto clade : node.GetClades()) {
      Assert(not clade.empty() && "Empty clade");
    }
    if (node.IsUA()) {
      if (storage.root_.value != NoId) {
        std::cout << "Duplicate root: " << storage.root_.value << " and "
                  << node.GetId().value << "\n";
      }
      Assert(storage.root_.value == NoId && "Duplicate root");
      storage.root_ = node;
    }
    if (node.IsLeaf()) {
      storage.leafs_.push_back(node);
    }
  }
}

template <typename CRTP, typename Tag>
void FeatureMutableView<Connections, CRTP, Tag>::BuildConnectionsRaw() const {
  auto& dag = static_cast<const CRTP&>(*this);
  for (auto node : dag.GetNodes()) {
    node.ClearConnections();
  }
  for (auto edge : dag.GetEdges()) {
    Assert(edge.GetParentId().value != NoId && "Edge has no parent");
    Assert(edge.GetChildId().value != NoId && "Edge has no child");
    Assert(edge.GetClade().value != NoId && "Edge has no clade index");
    Assert(edge.GetParentId() != edge.GetChildId() && "Edge is looped");
    edge.GetParent().AddEdge(edge.GetClade(), edge, true);
    edge.GetChild().AddEdge(edge.GetClade(), edge, false);
  }
}
