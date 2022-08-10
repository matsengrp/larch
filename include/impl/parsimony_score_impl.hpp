#include <algorithm>

template <typename MinWeightEdgeCallback, typename VisitNodeCallback>
typename ParsimonyScore<MinWeightEdgeCallback, VisitNodeCallback>::Weight
ParsimonyScore<MinWeightEdgeCallback, VisitNodeCallback>::ComputeLeaf(Node node) {
  return 0;
}

template <typename MinWeightEdgeCallback, typename VisitNodeCallback>
typename ParsimonyScore<MinWeightEdgeCallback, VisitNodeCallback>::Weight
ParsimonyScore<MinWeightEdgeCallback, VisitNodeCallback>::ComputeEdge(Edge edge) {
  return dag_.GetEdgeMutations(edge).size();
}

template <typename MinWeightEdgeCallback, typename VisitNodeCallback>
bool ParsimonyScore<MinWeightEdgeCallback, VisitNodeCallback>::Compare(Weight lhs,
                                                                       Weight rhs) {
  return lhs < rhs;
}

template <typename MinWeightEdgeCallback, typename VisitNodeCallback>
typename ParsimonyScore<MinWeightEdgeCallback, VisitNodeCallback>::Weight
ParsimonyScore<MinWeightEdgeCallback, VisitNodeCallback>::Combine(Weight lhs,
                                                                  Weight rhs) {
  return lhs + rhs;
}

template <typename MinWeightEdgeCallback, typename VisitNodeCallback>
void ParsimonyScore<MinWeightEdgeCallback, VisitNodeCallback>::MinWeightEdge(
    Edge edge) {
  min_weight_edge_callback_(edge);
}

template <typename MinWeightEdgeCallback, typename VisitNodeCallback>
void ParsimonyScore<MinWeightEdgeCallback, VisitNodeCallback>::VisitNode(
    Node node, Weight&& weight) {
  visit_node_callback_(node, std::forward<Weight>(weight));
}