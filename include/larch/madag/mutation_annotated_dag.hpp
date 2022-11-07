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
#include "larch/madag/compact_genome.hpp"
#include "larch/madag/edge_mutations.hpp"

namespace Mutation_Annotated_Tree {
class Tree;
class Node;
}  // namespace Mutation_Annotated_Tree

class ReferenceSequence {
 public:
 private:
  DAG_FEATURE_FRIENDS;
  std::string reference_sequence_;
};

template <typename View>
class FeatureReader<ReferenceSequence, View> {
 public:
  const std::string& GetReferenceSequence();
  void AssertUA();
  bool HaveUA();
};

template <typename View>
class FeatureWriter<ReferenceSequence, View>
    : public FeatureReader<ReferenceSequence, View> {
 public:
  void SetReferenceSequence(std::string_view reference_sequence);
  void AddUA(const EdgeMutations& mutations_at_root);
  void RecomputeCompactGenomes();
  void RecomputeEdgeMutations();
};

class SampleId {
 public:
 private:
  DAG_FEATURE_FRIENDS;
  std::optional<std::string> sample_id_;
};

template <typename View>
class FeatureReader<SampleId, View> {
 public:
  const std::optional<std::string>& GetSampleId();
};

template <typename View>
class FeatureWriter<SampleId, View> : public FeatureReader<SampleId, View> {
 public:
  void SetSampleId(const std::optional<std::string>& sample_id);
};

#include "larch/impl/madag/mutation_annotated_dag_impl.hpp"

using MADAGStorage = DefaultDAGStorage<
    DefaultNodesContainer<DefaultNodeStorage<CompactGenome, SampleId>>,
    DefaultEdgesContainer<DefaultEdgeStorage<EdgeMutations>>, ReferenceSequence>;

using MADAG = DAGView<const MADAGStorage, ReferenceSequence>;
using MutableMADAG = DAGView<MADAGStorage, ReferenceSequence>;
