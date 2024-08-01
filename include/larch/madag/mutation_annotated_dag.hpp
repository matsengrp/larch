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
#include <utility>
#include <optional>

#include "larch/dag/dag.hpp"
#include "larch/madag/mutation_base.hpp"
#include "larch/madag/compact_genome.hpp"
#include "larch/madag/edge_mutations.hpp"
#include "larch/madag/sample_id.hpp"

struct ReferenceSequence {
  MOVE_ONLY_DEF_CTOR(ReferenceSequence);
  std::string reference_sequence_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<ReferenceSequence, CRTP, Tag> {
  const std::string& GetReferenceSequence() const;
  void AssertUA() const;
  bool HaveUA() const;
};

using NodeSeqMap =
    IdContainer<NodeId, std::string, IdContinuity::Sparse, Ordering::Unordered>;
using NodeMutMap = IdContainer<NodeId, ContiguousMap<MutationPosition, MutationBase>,
                               IdContinuity::Sparse, Ordering::Unordered>;

template <typename CRTP, typename Tag>
struct FeatureMutableView<ReferenceSequence, CRTP, Tag> {
  void SetReferenceSequence(std::string_view reference_sequence) const;
  void SetCompactGenomesFromNodeSequenceMap(const NodeSeqMap& sequence_map) const;
  void SetCompactGenomesFromNodeMutationMap(NodeMutMap&& node_mutation_map) const;
  void UpdateCompactGenomesFromNodeMutationMap(NodeMutMap&& node_mutation_map) const;
  void AddUA(const EdgeMutations& mutations_at_root) const;
  template <IdContinuity Cont = IdContinuity::Dense>
  void RecomputeCompactGenomes(bool recompute_leaves = true) const;
  void SampleIdsFromCG(bool coerce = false) const;
  void RecomputeEdgeMutations() const;
};

#include "larch/impl/madag/mutation_annotated_dag_impl.hpp"

template <typename Target = DefaultDAGStorage,
          template <typename, typename> typename ViewBase = DefaultViewBase>
struct MADAGStorage;

template <typename Target, template <typename, typename> typename ViewBase>
struct LongNameOf<MADAGStorage<Target, ViewBase>> {
  using type =
      ExtendDAGStorage<MADAGStorage<Target>, Target,
                       Extend::Nodes<CompactGenome, Deduplicate<SampleId>>,
                       Extend::Edges<EdgeMutations>, Extend::DAG<ReferenceSequence>,
                       ViewBase, IdContinuity::Sparse>; // TODO change to Dense
};

template <typename Target, template <typename, typename> typename ViewBase>
struct MADAGStorage : LongNameOf<MADAGStorage<Target, ViewBase>>::type {
  SHORT_NAME(MADAGStorage);
};

template <typename Target, typename = std::enable_if_t<Target::role == Role::Storage>>
MADAGStorage<Target> AddMADAG(Target&& dag) {
  return MADAGStorage<Target>::Consume(std::move(dag));
}

template <typename Target, typename = std::enable_if_t<Target::role == Role::View>>
MADAGStorage<Target> AddMADAG(const Target& dag) {
  return MADAGStorage<Target>::FromView(dag);
}

using MADAG = typename MADAGStorage<DefaultDAGStorage>::ConstViewType;
using MutableMADAG = typename MADAGStorage<DefaultDAGStorage>::ViewType;
