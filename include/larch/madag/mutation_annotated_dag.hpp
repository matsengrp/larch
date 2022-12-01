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

struct ReferenceSequence {
  std::string reference_sequence_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<ReferenceSequence, CRTP, Tag> {
  const std::string& GetReferenceSequence() const {
    return GetFeatureStorage(this).reference_sequence_;
  }

  void AssertUA() const {
    auto& dag = static_cast<const CRTP&>(*this);
    Assert(dag.HaveRoot());
    auto ua = dag.GetRoot();
    Assert(ua.GetCladesCount() == 1);
  }

  bool HaveUA() const {
    auto& dag = static_cast<const CRTP&>(*this);
    Assert(dag.HaveRoot());
    auto ua = dag.GetRoot();
    return ua.GetCladesCount() == 1;
  }
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<ReferenceSequence, CRTP, Tag> {
  void SetReferenceSequence(std::string_view reference_sequence) const {
    GetFeatureStorage(this).reference_sequence_ = reference_sequence;
  }

  void AddUA(const EdgeMutations& mutations_at_root) const {
    auto dag = static_cast<const CRTP&>(*this);
    using Node = typename decltype(dag)::NodeView;
    using Edge = typename decltype(dag)::EdgeView;
    Assert(not dag.HaveUA());
    Node root = dag.GetRoot();
    Node ua_node = dag.AppendNode();
    Edge ua_edge = dag.AppendEdge(ua_node, root, {0});
    ua_edge.SetEdgeMutations(mutations_at_root.Copy());
    dag.BuildConnections();
    dag.AssertUA();
  }

  void RecomputeCompactGenomes() const {
    auto dag = static_cast<const CRTP&>(*this);
    using Node = typename decltype(dag)::NodeView;
    using Edge = typename decltype(dag)::EdgeView;

    std::vector<CompactGenome> new_cgs;
    new_cgs.resize(dag.GetNodesCount());
    auto ComputeCG = [&new_cgs, dag](auto& self, Node for_node) {
      CompactGenome& compact_genome = new_cgs.at(for_node.GetId().value);
      if (for_node.IsRoot()) {
        compact_genome = {};
        return;
      }
      if (not compact_genome.empty()) {
        return;
      }
      Edge edge = *(for_node.GetParents().begin());
      self(self, edge.GetParent());
      const EdgeMutations& mutations = edge.GetEdgeMutations();
      const CompactGenome& parent = new_cgs.at(edge.GetParentId().value);
      compact_genome.AddParentEdge(mutations, parent, dag.GetReferenceSequence());
    };
    for (Node node : dag.GetNodes()) {
      ComputeCG(ComputeCG, node);
    }
    for (Node node : dag.GetNodes()) {
      node.SetCompactGenome(std::move(new_cgs.at(node.GetId().value)));
    }
    std::unordered_map<CompactGenome, NodeId> leaf_cgs;
    for (Node node : dag.GetNodes()) {
      if (node.IsLeaf()) {
        bool success =
            leaf_cgs.emplace(node.GetCompactGenome().Copy(), node.GetId()).second;
        if (not success) {
          /*
          std::cout << "Error in ComputeCompactGenomes: had a non-unique leaf node at "
                    << node.GetId().value << " also seen at "
                    << leaf_cgs.at(node.GetCompactGenome().Copy()).value
                    << "\nCompact Genome is\n"
                    << node.GetCompactGenome().ToString() << "\n"
                    << std::flush;
          */
          Fail("Error in ComputeCompactGenomes: had a non-unique leaf node");
        }
      }
    }
  }

  void RecomputeEdgeMutations() const {
    auto dag = static_cast<const CRTP&>(*this);
    using Edge = typename decltype(dag)::EdgeView;

    for (Edge edge : dag.GetEdges()) {
      edge.SetEdgeMutations(CompactGenome::ToEdgeMutations(
          dag.GetReferenceSequence(), edge.GetParent().GetCompactGenome(),
          edge.GetChild().GetCompactGenome()));
    }
  }
};

struct SampleId {
  std::optional<std::string> sample_id_;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<SampleId, CRTP, Tag> {
  const std::optional<std::string>& GetSampleId() const {
    return GetFeatureStorage(this).sample_id_;
  }
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<SampleId, CRTP, Tag> {
  void SetSampleId(const std::optional<std::string>& sample_id) const {
    GetFeatureStorage(this).sample_id_ = sample_id;
  }
};

#include "larch/impl/madag/mutation_annotated_dag_impl.hpp"

using MADAGStorage = DAGStorage<
    ElementsContainer<NodeId, ElementStorage<Neighbors, CompactGenome, SampleId>>,
    ElementsContainer<EdgeId, ElementStorage<Endpoints, EdgeMutations>>, Connections,
    ReferenceSequence>;

using MADAG = DAGView<const MADAGStorage>;
using MutableMADAG = DAGView<MADAGStorage>;
