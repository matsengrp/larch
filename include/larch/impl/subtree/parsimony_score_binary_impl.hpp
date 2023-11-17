#include <algorithm>

template <typename DAG>
ParsimonyScore_::Weight ParsimonyScore_::ComputeLeaf(DAG, NodeId) const {
  return 0;
}

template <typename DAG>
ParsimonyScore_::Weight ParsimonyScore_::ComputeEdge(DAG dag, EdgeId edge_id) const {
  return dag.Get(edge_id).GetEdgeMutations().size();
}

bool ParsimonyScore_::Compare(ParsimonyScore_::Weight lhs,
                              ParsimonyScore_::Weight rhs) const {
  return lhs < rhs;
}

bool ParsimonyScore_::CompareEqual(ParsimonyScore_::Weight lhs,
                                   ParsimonyScore_::Weight rhs) const {
  return lhs == rhs;
}

inline bool ParsimonyScore_::IsIdentity(ParsimonyScore_::Weight weight) const {
  return weight == Identity;
}

ParsimonyScore_::Weight ParsimonyScore_::Combine(ParsimonyScore_::Weight lhs,
                                                 ParsimonyScore_::Weight rhs) const {
  return lhs + rhs;
}

bool MaxParsimonyScore_::Compare(MaxParsimonyScore_::Weight lhs,
                                 MaxParsimonyScore_::Weight rhs) const {
  return lhs > rhs;
}
