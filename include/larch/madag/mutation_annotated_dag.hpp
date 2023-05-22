/**
 * MADAG extends DAG, containing additional node data (CompactGenomes) and edge
 * data (EdgeMutations) in vectors ordered by node/edge IDs.
 *
 * When edge_mutations_ is populated, MADAG implements the Mutation Annotated
 * DAG object, in which nucleotide sequences on nodes are determined by the
 * nucleotide sequence on the UA node (reference_sequence_) and a sequence of
 * mutations associated to any path of edges from the UA node to the node of
 * interest.
 *
 * In order to merge mutation annotated DAGs, nodes must be compared using
 * their sequence and child clade sets. Node sequences are represented by
 * CompactGenome objects, which may be computed from edge mutations using
 * GetCompactGenomes().
 *
 */

#pragma once

#include <string_view>
#include <vector>
#include <unordered_map>
#include <utility>
#include <optional>

#include "larch/dag/dag.hpp"
#include "larch/madag/mutation_base.hpp"
#include "larch/madag/compact_genome.hpp"
#include "larch/madag/edge_mutations.hpp"

struct ReferenceSequence {
  std::string reference_sequence_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<ReferenceSequence, CRTP, Tag> {
  const std::string &GetReferenceSequence() const;
  void AssertUA() const;
  bool HaveUA() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<ReferenceSequence, CRTP, Tag> {
  void SetReferenceSequence(std::string_view reference_sequence) const;
  void SetLeafCompactGenomesFromSequenceMap(
      const std::unordered_map<NodeId, std::string> &leaf_sequence_map) const;
  void AddUA(const EdgeMutations &mutations_at_root) const;
  void RecomputeCompactGenomes(bool recompute_leaves = false) const;
  void RecomputeEdgeMutations() const;
};

struct SampleId {
  std::optional<std::string> sample_id_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<SampleId, CRTP, Tag> {
  const std::optional<std::string> &GetSampleId() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<SampleId, CRTP, Tag> {
  void SetSampleId(const std::optional<std::string> &sample_id) const;
};

#include "larch/impl/madag/mutation_annotated_dag_impl.hpp"

using MADAGStorage =
    ExtendDAGStorage<DefaultDAGStorage, Extend::Nodes<CompactGenome, SampleId>,
                     Extend::Edges<EdgeMutations>, Extend::DAG<ReferenceSequence>>;

using MADAG = DAGView<const MADAGStorage>;
using MutableMADAG = DAGView<MADAGStorage>;
