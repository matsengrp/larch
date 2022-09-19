#include <algorithm>

template <typename WeightOps>
WeightAccumulator<WeightOps>::Weight WeightAccumulator<WeightOps>::ComputeLeaf(const MADAG& dag, NodeId node_id) {
  return WeightCounter<WeightOps>({weight_ops.ComputeLeaf(dag, node_id)},
          std::forward<WeightOps>(weight_ops));
}

template <typename WeightOps>
WeightAccumulator<WeightOps>::Weight WeightAccumulator<WeightOps>::ComputeEdge(const MADAG& dag, EdgeId edge_id) {
  return WeightCounter<WeightOps>({weight_ops.ComputeEdge(dag, edge_id)},
          std::forward<WeightOps>(weight_ops));
}

template <typename WeightOps>
std::pair<WeightAccumulator<WeightOps>::Weight, std::vector<size_t>> WithinCladeAccumOptimum(std::vector<WeightAccumulator<WeightOps>::Weight> inweights) {
    std::vector<size_t> optimal_indices;
    std::iota(optimal_indices.begin(), optimal_indices.end(), 0);
    return {std::accumulate(inweights.begin(), inweights.end(), WeightAccumulator<WeightOps>(std::forward<WeightOps>(weight_ops_)), WeightCounter<WeightOps>::operator+()), optimal_indices};
}

WeightAccumulator<WeightOps>::Weight BetweenClades(std::vector<WeightAccumulator<WeightOps>::Weight> inweights) {
    return std::accumulate(inweights.begin(), inweights.end(), WeightAccumulator<WeightOps>(std::forward<WeightOps>(weight_ops_)), WeightCounter<WeightOps>::operator*());
}

WeightAccumulator<WeightOps>::Weight AboveNode(WeightAccumulator<WeightOps>::Weight edgeweight, WeightAccumulator<WeightOps>::Weight childnodeweight) {
    // because edgeweight should have come from ComputeEdge:
    assert(edgeweight.size() == 1);
    auto edgepair = edgeweight.begin();
    std::map<typename WeightOps::Weight, size_t> result;
    for (auto const& childitem : childnodeweight) {
        result[weight_ops_.AboveNode(edgepair->first, childitem.first)] += childitem.second;
    }
    return WeightCounter<WeightOps>(result, std::forward<WeightOps>(weight_ops_));
}



template <typename WeightOps>
WeightCounter<WeightOps>::WeightCounter(const WeightOps&& weight_ops) : weight_ops_{weight_ops} {}

template <typename WeightOps>
WeightCounter<WeightOps>::WeightCounter(std::vector<typename WeightOps::Weight> inweights, const WeightOps&& weight_ops) : weight_ops_{weight_ops} {
    for (auto weight : inweights) {
        weights_[weight] ++;
    }
}

template <typename WeightOps>
WeightCounter<WeightOps>::WeightCounter(std::map<typename WeightOps::Weight, size_t> inweights, const WeightOps&& weight_ops) : weights_{inweights}, weight_ops_{weight_ops} {
    for (auto weight : inweights) {
        weights_[weight] ++;
    }
}

template <typename WeightOps>
std::map<typename WeightOps::Weight, size_t> WeightCounter<WeightOps>::GetWeights() {
    return weights_;
}

template<typename WeightOps>
WeightCounter<WeightOps> WeightCounter<WeightOps>::operator+(
        WeightCounter<WeightOps> lhs,
        WeightCounter<WeightOps> rhs) {
    std::map<typename WeightOps::Weight, size_t> result = lhs.GetWeights();
    for (auto const& map_pair : rhs.GetWeights()) {
        result[map_pair.first] += map_pair.second;
    }
    return WeightCounter<WeightOps>(result, std::forward<WeightOps>(weight_ops_));
}
        
template<typename WeightOps>
WeightCounter<WeightOps> WeightCounter<WeightOps>::operator*(
        WeightCounter<WeightOps> lhs,
        WeightCounter<WeightOps> rhs) {
    std::map<typename WeightOps::Weight, size_t> result;
    for (auto const& lpair : lhs.GetWeights()) {
        for (auto const& rpair : rhs.GetWeights()) {
            auto value = weight_ops_.BetweenClades({lpair.first, rpair.first});
            result[value] += lpair.second * rpair.second;
        }
    }
    return WeightCounter<WeightOps>(result, std::forward<WeightOps>(weight_ops_));
}
