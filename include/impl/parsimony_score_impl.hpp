#include <algorithm>

ParsimonyScore::Weight ParsimonyScore::ComputeLeaf(const MADAG& dag, NodeId node_id) {
  return 0;
}

ParsimonyScore::Weight ParsimonyScore::ComputeEdge(const MADAG& dag, EdgeId edge_id) {
  return dag.GetEdgeMutations(edge_id).size();
}

bool ParsimonyScore::Compare(Weight lhs, Weight rhs) {
  return lhs < rhs;
}

ParsimonyScore::Weight ParsimonyScore::Combine(Weight lhs, Weight rhs) {
  return lhs + rhs;
}