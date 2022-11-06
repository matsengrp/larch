#include <algorithm>

template <typename WeightOps>
WeightAccumulator<WeightOps>::WeightAccumulator(const WeightOps& ops)
    : weight_ops_{ops} {}

template <typename WeightOps>
typename WeightAccumulator<WeightOps>::Weight WeightAccumulator<WeightOps>::ComputeLeaf(
    MADAG dag, NodeId node_id) {
  return WeightCounter<WeightOps>({weight_ops_.ComputeLeaf(dag, node_id)}, weight_ops_);
}

template <typename WeightOps>
typename WeightAccumulator<WeightOps>::Weight WeightAccumulator<WeightOps>::ComputeEdge(
    MADAG dag, EdgeId edge_id) {
  return WeightCounter<WeightOps>({weight_ops_.ComputeEdge(dag, edge_id)}, weight_ops_);
}

template <typename WeightOps>
std::pair<typename WeightAccumulator<WeightOps>::Weight, std::vector<size_t>>
WeightAccumulator<WeightOps>::WithinCladeAccumOptimum(
    const std::vector<typename WeightAccumulator<WeightOps>::Weight>& inweights) {
  std::vector<size_t> optimal_indices;
  std::iota(optimal_indices.begin(), optimal_indices.end(), 0);
  return {std::accumulate(inweights.begin(), inweights.end(), Weight{weight_ops_},
                          [](auto& lhs, auto& rhs) { return lhs + rhs; }),
          optimal_indices};
}

template <typename WeightOps>
typename WeightAccumulator<WeightOps>::Weight
WeightAccumulator<WeightOps>::BetweenClades(
    const std::vector<typename WeightAccumulator<WeightOps>::Weight>& inweights) const {
  if (inweights.size() == 1) {
    return inweights[0];
  }
  return std::accumulate(std::next(inweights.begin()), inweights.end(),
                         *inweights.begin(),
                         [](auto& lhs, auto& rhs) { return lhs * rhs; });
}

template <typename WeightOps>
typename WeightAccumulator<WeightOps>::Weight WeightAccumulator<WeightOps>::AboveNode(
    typename WeightAccumulator<WeightOps>::Weight edgeweight,
    typename WeightAccumulator<WeightOps>::Weight childnodeweight) {
  // because edgeweight should have come from ComputeEdge:
  Assert(edgeweight.GetWeights().size() == 1);
  auto edgepair = edgeweight.GetWeights().begin();
  std::map<typename WeightOps::Weight, Count> result;
  for (auto const& childitem : childnodeweight.GetWeights()) {
    result[weight_ops_.AboveNode(edgepair->first, childitem.first)] += childitem.second;
  }
  return WeightCounter<WeightOps>(std::move(result),
                                  std::forward<WeightOps>(weight_ops_));
}

template <typename WeightOps>
WeightCounter<WeightOps>::WeightCounter(WeightOps&& weight_ops)
    : weight_ops_{weight_ops} {}

template <typename WeightOps>
WeightCounter<WeightOps>::WeightCounter(
    const std::vector<typename WeightOps::Weight>& inweights,
    const WeightOps& weight_ops)
    : weight_ops_{weight_ops} {
  for (const auto& weight : inweights) {
    weights_[weight]++;
  }
}

template <typename WeightOps>
WeightCounter<WeightOps>::WeightCounter(
    std::map<typename WeightOps::Weight, Count>&& inweights,
    const WeightOps& weight_ops)
    : weights_{inweights}, weight_ops_{weight_ops} {}

template <typename WeightOps>
const std::map<typename WeightOps::Weight, Count>&
WeightCounter<WeightOps>::GetWeights() const {
  return weights_;
}

template <typename WeightOps>
const WeightOps& WeightCounter<WeightOps>::GetWeightOps() const {
  return weight_ops_;
}

template <typename WeightOps>
WeightCounter<WeightOps> WeightCounter<WeightOps>::operator+(
    const WeightCounter<WeightOps>& rhs) const {
  std::map<typename WeightOps::Weight, Count> result = weights_;
  for (auto const& map_pair : rhs.GetWeights()) {
    result[map_pair.first] += map_pair.second;
  }
  return WeightCounter<WeightOps>(std::move(result), weight_ops_);
}

template <typename WeightOps>
WeightCounter<WeightOps> WeightCounter<WeightOps>::operator*(
    const WeightCounter<WeightOps>& rhs) const {
  std::map<typename WeightOps::Weight, Count> result;
  for (auto const& lpair : weights_) {
    for (auto const& rpair : rhs.GetWeights()) {
      auto value = weight_ops_.BetweenClades({lpair.first, rpair.first});
      result[value] += lpair.second * rpair.second;
    }
  }
  return WeightCounter<WeightOps>(std::move(result), weight_ops_);
}

template <typename WeightOps>
WeightCounter<WeightOps>& WeightCounter<WeightOps>::operator=(
    const WeightCounter<WeightOps>& rhs) {
  weights_ = rhs.weights_;
  weight_ops_ = rhs.weight_ops_;
  return *this;
}

template <typename WeightOps>
WeightCounter<WeightOps>& WeightCounter<WeightOps>::operator=(
    WeightCounter<WeightOps>&& rhs) noexcept {
  weights_ = rhs.weights_;
  weight_ops_ = rhs.weight_ops_;
  return *this;
}

template <typename WeightOps>
bool WeightCounter<WeightOps>::operator==(const WeightCounter<WeightOps>& rhs) {
  return weights_ == rhs.GetWeights();
}

template <typename WeightOps>
bool WeightCounter<WeightOps>::operator!=(const WeightCounter<WeightOps>& rhs) {
  return weights_ != rhs.GetWeights();
}

template <typename WeightOps>
std::ostream& operator<<(std::ostream& os,
                         const WeightCounter<WeightOps>& weight_counter) {
  auto weights = weight_counter.GetWeights();
  if (weights.empty()) {
    os << "{}";
    return os;
  }
  auto addpair = [&os](auto const& mappair) {
    os << mappair.first << ": " << mappair.second;
  };
  os << "{";
  for (auto i = weights.begin(); i != --weights.end(); i++) {
    addpair(*i);
    os << ", ";
  }
  addpair(*(--weights.end()));
  os << "}";
  return os;
}
