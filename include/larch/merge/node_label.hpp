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
  inline NodeLabel(const CompactGenome* cg, const LeafSet* ls, const SampleId* id);

  inline const CompactGenome* GetCompactGenome() const;
  inline const LeafSet* GetLeafSet() const;
  inline const SampleId* GetSampleId() const;

  inline void SetCompactGenome(const CompactGenome* cg);
  inline void SetLeafSet(const LeafSet* ls);
  inline void SetSampleId(const SampleId* id);

  inline bool operator==(const NodeLabel& rhs) const noexcept;

  [[nodiscard]] inline size_t Hash() const noexcept;

  inline bool Empty() const {
    Assert(compact_genome_);
    Assert(leaf_set_);
    Assert(sample_id_);
    if (compact_genome_ == CompactGenome::Empty() and sample_id_ == SampleId::Empty()) {
      return true;
    }
    if (leaf_set_ == LeafSet::Empty()) {
      return true;
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
  const SampleId* sample_id_;
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
