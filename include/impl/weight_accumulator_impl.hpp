#include <algorithm>

template <typename WeightOps>
typename WeightAccumulator<WeightOps>::Weight WeightAccumulator<WeightOps>::ComputeLeaf(const MADAG& dag, NodeId node_id) {
  return WeightCounter<WeightOps>({weight_ops_.ComputeLeaf(dag, node_id)},
          std::forward<WeightOps>(weight_ops_));
}

template <typename WeightOps>
typename WeightAccumulator<WeightOps>::Weight WeightAccumulator<WeightOps>::ComputeEdge(const MADAG& dag, EdgeId edge_id) {
  return WeightCounter<WeightOps>({weight_ops_.ComputeEdge(dag, edge_id)},
          std::forward<WeightOps>(weight_ops_));
}

template <typename WeightOps>
std::pair<typename WeightAccumulator<WeightOps>::Weight, std::vector<size_t>> WeightAccumulator<WeightOps>::WithinCladeAccumOptimum(std::vector<typename WeightAccumulator<WeightOps>::Weight> inweights) {
    std::vector<size_t> optimal_indices;
    std::iota(optimal_indices.begin(), optimal_indices.end(), 0);
    return {std::accumulate(inweights.begin(), inweights.end(), WeightAccumulator<WeightOps>(std::forward<WeightOps>(weight_ops_)), WeightCounter<WeightOps>::operator+()), optimal_indices};
}

template <typename WeightOps>
typename WeightAccumulator<WeightOps>::Weight WeightAccumulator<WeightOps>::BetweenClades(std::vector<typename WeightAccumulator<WeightOps>::Weight> inweights) {
    return std::accumulate(inweights.begin(), inweights.end(), WeightAccumulator<WeightOps>(std::forward<WeightOps>(weight_ops_)), WeightCounter<WeightOps>::operator*());
}

template <typename WeightOps>
typename WeightAccumulator<WeightOps>::Weight WeightAccumulator<WeightOps>::AboveNode(typename WeightAccumulator<WeightOps>::Weight edgeweight, typename WeightAccumulator<WeightOps>::Weight childnodeweight) {
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
WeightCounter<WeightOps>::WeightCounter(WeightOps&& weight_ops) : weight_ops_{weight_ops} {}

template <typename WeightOps>
WeightCounter<WeightOps>::WeightCounter(WeightCounter<WeightOps>& incounter) : WeightCounter(incounter.GetWeights(), incounter.GetWeightOps()) {}

template <typename WeightOps>
WeightCounter<WeightOps>::WeightCounter(std::vector<typename WeightOps::Weight> inweights, WeightOps&& weight_ops) : weight_ops_{weight_ops} {
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

template <typename WeightOps>
WeightOps&& WeightCounter<WeightOps>::GetWeightOps() {
    return weight_ops_;
}

template<typename WeightOps>
WeightCounter<WeightOps> WeightCounter<WeightOps>::operator+(WeightCounter<WeightOps> rhs) {
    std::map<typename WeightOps::Weight, size_t> result = weights_;
    for (auto const& map_pair : rhs.GetWeights()) {
        result[map_pair.first] += map_pair.second;
    }
    return WeightCounter<WeightOps>(result, std::forward<WeightOps>(weight_ops_));
}
        
template<typename WeightOps>
WeightCounter<WeightOps> WeightCounter<WeightOps>::operator*(WeightCounter<WeightOps> rhs) {
    std::map<typename WeightOps::Weight, size_t> result;
    for (auto const& lpair : weights_) {
        for (auto const& rpair : rhs.GetWeights()) {
            auto value = weight_ops_.BetweenClades({lpair.first, rpair.first});
            result[value] += lpair.second * rpair.second;
        }
    }
    return WeightCounter<WeightOps>(result, std::forward<WeightOps>(weight_ops_));
}

template<typename WeightOps>
WeightCounter<WeightOps> WeightCounter<WeightOps>::operator=(WeightCounter<WeightOps> rhs) { return WeightCounter(&rhs); }

template<typename WeightOps>
WeightCounter<WeightOps>& WeightCounter<WeightOps>::operator=(WeightCounter<WeightOps>& rhs) { return WeightCounter(rhs); }

template<typename WeightOps>
WeightCounter<WeightOps>&& WeightCounter<WeightOps>::operator=(WeightCounter<WeightOps>&& rhs) { return rhs; }
