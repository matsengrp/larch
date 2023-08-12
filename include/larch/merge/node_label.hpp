#pragma once

#include <cstddef>

class CompactGenome;
class LeafSet;

/*
 * NodeLabel stores the data which formally defines a node: a CompactGenome,
 * and a set of child clades, stored in a LeafSet. This allows appropriate
 * equality testing between node objects during MAD operations, like merging.
 */
class NodeLabel {
 public:
  inline NodeLabel();
  inline NodeLabel(const CompactGenome* cg, const LeafSet* ls);

  inline const CompactGenome* GetCompactGenome() const;
  inline const LeafSet* GetLeafSet() const;

  inline void SetCompactGenome(const CompactGenome* cg);
  inline void SetLeafSet(const LeafSet* ls);

  inline bool operator==(const NodeLabel& rhs) const noexcept;

  [[nodiscard]] inline size_t Hash() const noexcept;

  inline bool Empty() const {
    if (compact_genome_ == nullptr or leaf_set_ == nullptr) {
      return true;
    }
    if (compact_genome_ == CompactGenome::Empty() or leaf_set_ == LeafSet::Empty()) {
      return true;
    }
    if (compact_genome_->empty() and leaf_set_->empty()) {
      // return true;
    }
    return false;
  }

  inline std::string ToString() const {
    std::string result = "[";
    result += compact_genome_->ToString();
    result += "(";
    result += leaf_set_->ToString();
    result += ")";
    return result;
  }

 private:
  const CompactGenome* compact_genome_;
  const LeafSet* leaf_set_;
};

template <>
struct std::hash<NodeLabel> {
  inline std::size_t operator()(const NodeLabel& nl) const noexcept;
};

template <>
struct std::equal_to<NodeLabel> {
  inline bool operator()(const NodeLabel& lhs, const NodeLabel& rhs) const noexcept;
};

#include "larch/impl/merge/node_label_impl.hpp"
