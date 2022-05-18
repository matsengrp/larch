#include <unordered_map>
#include <numeric>

auto HistoryDAG::GetNodes() const {
  return nodes_ |
         ranges::views::transform([this, idx = size_t{}](const NodeStorage&) mutable {
           return Node{*this, {idx++}};
         });
}

auto HistoryDAG::GetNodes() {
  return nodes_ |
         ranges::views::transform([this, idx = size_t{}](NodeStorage&) mutable {
           return MutableNode{*this, {idx++}};
         });
}

auto HistoryDAG::GetEdges() const {
  return edges_ | ranges::views::transform(
                      [this, idx = size_t{}](const EdgeStorage<Weight>&) mutable {
                        return Edge{*this, {idx++}};
                      });
}

auto HistoryDAG::GetEdges() {
  return edges_ |
         ranges::views::transform([this, idx = size_t{}](EdgeStorage<Weight>&) mutable {
           return MutableEdge{*this, {idx++}};
         });
}

auto HistoryDAG::GetLeafs() const {
  return leafs_ | ranges::views::transform([this](NodeId node_id) {
           return Node{*this, node_id};
         });
}

auto HistoryDAG::GetLeafs() {
  return leafs_ | ranges::views::transform([this](NodeId node_id) {
           return MutableNode{*this, node_id};
         });
}

auto HistoryDAG::TraversePreOrder() const {
  return ranges::subrange(PreOrderIterator{GetRoot()}, PreOrderIterator<Node>{});
}

auto HistoryDAG::TraversePreOrder() {
  return ranges::subrange(PreOrderIterator{GetRoot()}, PreOrderIterator<MutableNode>{});
}

auto HistoryDAG::TraversePostOrder() const {
  return ranges::subrange(PostOrderIterator{GetRoot()}, PostOrderIterator<Node>{});
}

auto HistoryDAG::TraversePostOrder() {
  return ranges::subrange(PostOrderIterator{GetRoot()},
                          PostOrderIterator<MutableNode>{});
}
