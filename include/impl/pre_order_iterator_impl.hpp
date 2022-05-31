template <typename NodeType>
PreOrderIterator<NodeType>::PreOrderIterator(NodeType node) {
  stack_.push(*node.GetChildren().begin());
}

template <typename NodeType>
typename PreOrderIterator<NodeType>::value_type PreOrderIterator<NodeType>::operator*()
    const {
  Assert(not stack_.empty());
  EdgeType top = stack_.top();
  return {top.GetDAG(), root_visited_ ? top.GetChild() : top.GetParent(), top};
}

template <typename NodeType>
PreOrderIterator<NodeType>& PreOrderIterator<NodeType>::operator++() {
  Assert(not stack_.empty());
  auto top = stack_.top();
  if (not root_visited_) {
    root_visited_ = true;
    return *this;
  }
  auto top_first_child = GetFirstChild(top);
  if (top_first_child.has_value()) {
    stack_.push(*top_first_child);
  } else {
    stack_.pop();
    auto next = top.FindNextSibling();
    if (next.has_value()) {
      stack_.push(*next);
    } else {
      std::optional<EdgeType> next;
      auto NextFound = [&]() {
        auto top_next = stack_.top().FindNextSibling();
        if (top_next.has_value()) {
          next.emplace(*top_next);
          return true;
        }
        return false;
      };
      while (not stack_.empty() and not NextFound()) stack_.pop();
      if (not stack_.empty()) {
        stack_.pop();
        stack_.push(*next);
      }
    }
  }
  return *this;
}

template <typename NodeType>
PreOrderIterator<NodeType> PreOrderIterator<NodeType>::operator++(int) {
  PreOrderIterator result = *this;
  this->operator++();
  return result;
}

template <typename NodeType>
bool PreOrderIterator<NodeType>::operator==(const PreOrderIterator& other) const {
  return stack_ == other.stack_;
}

template <typename NodeType>
bool PreOrderIterator<NodeType>::operator!=(const PreOrderIterator& other) const {
  return not(*this == other);
}

template <typename NodeType>
std::optional<typename PreOrderIterator<NodeType>::EdgeType>
PreOrderIterator<NodeType>::GetFirstChild(EdgeType edge) {
  if (edge.IsLeaf()) return std::nullopt;
  return *edge.GetChild().GetChildren().begin();
}
