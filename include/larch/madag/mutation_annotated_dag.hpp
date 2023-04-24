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

// Strongly typed class to represent a base
struct MutationBase {
  using BitArray = std::array<bool, 2>;

  MutationBase(const BitArray m_value) { value = m_value; };
  MutationBase(const char m_char_in) {
    for (const auto &[m_base, m_char] : mut_to_char_map) {
      if (m_char_in == m_char) {
        value = m_base.value;
        return;
      }
    }
    // ASSERT(false, "ERROR: Invalid char given for MutationBase constructor.");
  }

  // NOTE: Due to mapping, we can simply invert the bool vector to get a base's
  // complement.
  MutationBase GetComplement() const {
    MutationBase m_out{{!value[0], !value[1]}};
    return m_out;
  }

  char ToChar() const { return mut_to_char_map.find(value)->second; }

  friend std::ostream &operator<<(std::ostream &os, const MutationBase m_in) {
    os << m_in.ToChar();
    return os;
  }

  static std::string ToString(std::vector<MutationBase> m_in) {
    std::string str_out = "";
    for (size_t i = 0; i < m_in.size(); i++) {
      str_out[i] += m_in[i].ToChar();
    }
    return str_out;
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const std::vector<MutationBase> m_in) {
    os << ToString(m_in);
    return os;
  }

  friend bool operator==(const MutationBase &lhs, const MutationBase &rhs) {
    return lhs.value == rhs.value;
  }
  friend bool operator<(const MutationBase &lhs, const MutationBase &rhs) {
    return lhs.value < rhs.value;
  }
  friend bool operator==(const MutationBase &lhs, const BitArray &rhs) {
    return lhs.value == rhs;
  }
  friend bool operator==(const MutationBase &lhs, const char &rhs) {
    return lhs.ToChar() == rhs;
  }

  // Two-bit array represents the four possible DNA bases.
  BitArray value;

  struct DNA {
    static const MutationBase A, C, G, T;
  };
  static const std::array<MutationBase, 4> bases;
  static const std::map<MutationBase, char> mut_to_char_map;
};
const MutationBase MutationBase::DNA::A = MutationBase({0, 0});
const MutationBase MutationBase::DNA::C = MutationBase({0, 1});
const MutationBase MutationBase::DNA::G = MutationBase({1, 0});
const MutationBase MutationBase::DNA::T = MutationBase({1, 1});
const std::array<MutationBase, 4> MutationBase::bases = {
    MutationBase::DNA::A, MutationBase::DNA::C, MutationBase::DNA::G,
    MutationBase::DNA::T};
const std::map<MutationBase, char> MutationBase::mut_to_char_map = {
    {MutationBase::DNA::A, 'A'},
    {MutationBase::DNA::C, 'C'},
    {MutationBase::DNA::G, 'G'},
    {MutationBase::DNA::T, 'T'}};

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
  void AddUA(const EdgeMutations &mutations_at_root) const;
  void RecomputeCompactGenomes() const;
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
