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
  return edges_ | ranges::views::transform(
                      [this, idx = size_t{}](const EdgeStorage<Weight>&) mutable {
                        return Edge{*this, {idx++}};
                      });
}

auto DAG::GetEdges() {
  return edges_ |
         ranges::views::transform([this, idx = size_t{}](EdgeStorage<Weight>&) mutable {
           return MutableEdge{*this, {idx++}};
         });
}

auto DAG::GetLeafs() const { return leafs_ | Transform::ToNodes(*this); }

auto DAG::GetLeafs() { return leafs_ | Transform::ToNodes(*this); }

auto DAG::TraversePreOrder() const {
  return ranges::subrange(PreOrderIterator{GetRoot()}, PreOrderIterator<Node>{});
}

auto DAG::TraversePreOrder() {
  return ranges::subrange(PreOrderIterator{GetRoot()}, PreOrderIterator<MutableNode>{});
}

auto DAG::TraversePostOrder() const {
  return ranges::subrange(PostOrderIterator{GetRoot()}, PostOrderIterator<Node>{});
}

auto DAG::TraversePostOrder() {
  return ranges::subrange(PostOrderIterator{GetRoot()},
                          PostOrderIterator<MutableNode>{});
}
