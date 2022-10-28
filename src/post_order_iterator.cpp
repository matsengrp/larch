#include "post_order_iterator.hpp"

#include "dag.hpp"

#if 0
template <typename NodeType>
PostOrderIterator<NodeType>::PostOrderIterator(NodeType node) {
  stack_.push(*node.GetChildren().begin());
  PushToNextLeaf();
}

template <typename NodeType>
typename PostOrderIterator<NodeType>::value_type
PostOrderIterator<NodeType>::operator*() {
  Assert(not stack_.empty());
  EdgeType edge = stack_.top();
  NodeType node = visit_root_ ? edge.GetParent() : edge.GetChild();
  if (not visited_nodes_.insert(node).second) {
    this->operator++();
    return this->operator*();
  }
  return {edge.GetDAG(), node, edge};
}

template <typename NodeType>
PostOrderIterator<NodeType>& PostOrderIterator<NodeType>::operator++() {
  Assert(not end_sentinel_);
  Assert(not stack_.empty());

  if (visit_root_) {
    end_sentinel_ = true;
    return *this;
  }

  EdgeType top = stack_.top();
  auto sibling = top.FindNextSibling();

  if (top.IsRoot() && not sibling) {
    visit_root_ = true;
    return *this;
  }

  stack_.pop();

  if (sibling) {
    stack_.push(*sibling);
    PushToNextLeaf();
  }
  return *this;
}

template <typename NodeType>
PostOrderIterator<NodeType> PostOrderIterator<NodeType>::operator++(int) {
  PostOrderIterator result = *this;
  this->operator++();
  return result;
}

template <typename NodeType>
bool PostOrderIterator<NodeType>::operator==(const PostOrderIterator& other) const {
  if (other.stack_.empty() && end_sentinel_) {
    return true;
  }
  return stack_ == other.stack_;
}

template <typename NodeType>
bool PostOrderIterator<NodeType>::operator!=(const PostOrderIterator& other) const {
  return not(*this == other);
}

template <typename NodeType>
void PostOrderIterator<NodeType>::PushToNextLeaf() {
  while (not stack_.top().IsLeaf()) {
    stack_.push(*stack_.top().GetChild().GetChildren().begin());
  }
}

template class PostOrderIterator<Node>;
template class PostOrderIterator<MutableNode>;
#endif