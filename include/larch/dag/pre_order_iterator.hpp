#pragma once

#include <utility>
#include <cstddef>
#include <iterator>
#include <stack>
#include <set>

#include "larch/dag/traverse_value.hpp"

/*
 * PreOrderIterator provides a pre-order node iterator, given a DAG node.
 */
#if 0
template <typename NodeType>
class PreOrderIterator {
 public:
  using EdgeType = decltype(*std::declval<NodeType>().GetChildren().begin());

  using iterator_category = std::forward_iterator_tag;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using value_type =
      TraverseValue<std::conditional_t<NodeType::is_mutable, DAG&, const DAG&>>;
  using pointer = value_type*;
  using reference = value_type&;
  using const_pointer = const pointer;
  using const_reference = const value_type&;

  explicit PreOrderIterator(NodeType node);
  PreOrderIterator() = default;
  value_type operator*();
  PreOrderIterator& operator++();
  PreOrderIterator operator++(int);
  bool operator==(const PreOrderIterator& other) const;
  bool operator!=(const PreOrderIterator& other) const;

 private:
  static std::optional<EdgeType> GetFirstChild(EdgeType edge);

  std::stack<EdgeType> stack_;
  std::set<NodeId> visited_nodes_;
  bool root_visited_ = false;
};
#endif