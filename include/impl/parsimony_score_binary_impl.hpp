#include <algorithm>

ParsimonyScore_::Weight ParsimonyScore_::ComputeLeaf(const MADAG& dag, NodeId node_id) {
  return 0;
}

ParsimonyScore_::Weight ParsimonyScore_::ComputeEdge(const MADAG& dag, EdgeId edge_id) {
  return dag.GetEdgeMutations(edge_id).size();
}

bool ParsimonyScore_::Compare(ParsimonyScore_::Weight lhs,
                              ParsimonyScore_::Weight rhs) {
  return lhs < rhs;
}

bool ParsimonyScore_::CompareEqual(ParsimonyScore_::Weight lhs,
                                   ParsimonyScore_::Weight rhs) {
  return lhs == rhs;
}

inline bool ParsimonyScore_::IsIdentity(ParsimonyScore_::Weight weight) {
  return weight == Identity;
}

ParsimonyScore_::Weight ParsimonyScore_::Combine(ParsimonyScore_::Weight lhs,
                                                 ParsimonyScore_::Weight rhs) {
  return lhs + rhs;
}
