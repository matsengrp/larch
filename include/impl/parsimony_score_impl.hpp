#include <algorithm>

ParsimonyScore::Weight ParsimonyScore::ComputeLeaf(Node node) { return 0; }

ParsimonyScore::Weight ParsimonyScore::ComputeEdge(Edge edge) {
  return dag_.GetEdgeMutations(edge).size();
}

bool ParsimonyScore::Compare(Weight lhs, Weight rhs) { return lhs < rhs; }

ParsimonyScore::Weight ParsimonyScore::Combine(Weight lhs, Weight rhs) {
  return lhs + rhs;
}

void ParsimonyScore::MinWeightEdge(Edge edge) {}

void ParsimonyScore::VisitNode(Node node, Weight&& weight) {}