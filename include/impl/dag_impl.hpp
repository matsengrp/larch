// The functions in this file are documented where declared in `include/dag.hpp`
#include <unordered_map>
#include <numeric>

#include <range/v3/view/subrange.hpp>

auto DAG::GetNodes() const {
  return nodes_ |
         ranges::views::transform([this, idx = size_t{}](const NodeStorage&) mutable {
           return Node{*this, {idx++}};
         });
}

auto DAG::GetNodes() {
  return nodes_ |
         ranges::views::transform([this, idx = size_t{}](NodeStorage&) mutable {
           return MutableNode{*this, {idx++}};
         });
}

auto DAG::GetEdges() const {
  return edges_ |
         ranges::views::transform([this, idx = size_t{}](const EdgeStorage&) mutable {
           return Edge{*this, {idx++}};
         });
}

auto DAG::GetEdges() {
  return edges_ |
         ranges::views::transform([this, idx = size_t{}](EdgeStorage&) mutable {
           return MutableEdge{*this, {idx++}};
         });
}

auto DAG::GetLeafs() const { return leafs_ | Transform::ToNodes(*this); }

auto DAG::GetLeafs() { return leafs_ | Transform::ToNodes(*this); }

#if 0
auto DAG::TraversePreOrder(NodeId below_node) const {
  return ranges::subrange(PreOrderIterator{Get(below_node)}, PreOrderIterator<Node>{});
}

auto DAG::TraversePreOrder(NodeId below_node) {
  return ranges::subrange(PreOrderIterator{Get(below_node)},
                          PreOrderIterator<MutableNode>{});
}

auto DAG::TraversePostOrder(NodeId below_node) const {
  return ranges::subrange(PostOrderIterator{Get(below_node)},
                          PostOrderIterator<Node>{});
}

auto DAG::TraversePostOrder(NodeId below_node) {
  return ranges::subrange(PostOrderIterator{Get(below_node)},
                          PostOrderIterator<MutableNode>{});
}

auto DAG::TraversePreOrder() const { return TraversePreOrder(GetRoot()); }

auto DAG::TraversePreOrder() { return TraversePreOrder(GetRoot()); }

auto DAG::TraversePostOrder() const { return TraversePostOrder(GetRoot()); }

auto DAG::TraversePostOrder() { return TraversePostOrder(GetRoot()); }
#endif
