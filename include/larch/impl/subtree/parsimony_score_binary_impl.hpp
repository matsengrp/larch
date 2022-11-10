#include <algorithm>

template <typename Node>
ParsimonyScore_::Weight ParsimonyScore_::ComputeLeaf(Node) {
  return 0;
}

template <typename Edge>
ParsimonyScore_::Weight ParsimonyScore_::ComputeEdge(Edge edge) {
  return edge.GetEdgeMutations().size();
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
