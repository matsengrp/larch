// Implementation of Sankoff algorithm
// Included from larch/subtree/sankoff.hpp

#include <algorithm>
#include <stdexcept>

inline size_t BaseToIndex(char base) {
  switch (base) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
    default:
      throw std::invalid_argument(std::string("Invalid base: ") + base);
  }
}

inline char IndexToBase(size_t index) {
  static const char bases[] = {'A', 'C', 'G', 'T'};
  if (index >= 4) {
    throw std::out_of_range("Base index must be 0-3");
  }
  return bases[index];
}

inline std::vector<size_t> GetCompatibleIndices(MutationBase base) {
  std::vector<size_t> result;
  char c = base.ToChar();
  if (c == 'N') {
    // Ambiguous: all bases are compatible
    result = {0, 1, 2, 3};
  } else {
    result.push_back(BaseToIndex(c));
  }
  return result;
}

// ============================================================================
// SankoffScorer implementation
// ============================================================================

template <typename DAG>
SankoffScorer<DAG>::SankoffScorer(DAG dag) : dag_{dag} {
  CollectVariableSites();
  InitializeCostMatrices();
}

template <typename DAG>
SankoffScorer<DAG>::SankoffScorer(
    DAG dag, const std::unordered_map<size_t, SiteCostMatrix>& cost_matrices)
    : dag_{dag}, custom_cost_matrices_{cost_matrices} {
  CollectVariableSites();
  InitializeCostMatrices();
}

template <typename DAG>
SankoffScorer<DAG>::SankoffScorer(DAG dag, const SiteCostMatrix& uniform_matrix)
    : dag_{dag}, uniform_cost_matrix_{uniform_matrix} {
  CollectVariableSites();
  InitializeCostMatrices();
}

template <typename DAG>
SiteCostMatrix SankoffScorer<DAG>::MakeUniformCostMatrix() {
  SiteCostMatrix matrix;
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      matrix[i][j] = (i == j) ? 0.0 : 1.0;
    }
  }
  return matrix;
}

template <typename DAG>
void SankoffScorer<DAG>::CollectVariableSites() {
  // Collect all positions that have mutations on any edge
  std::set<size_t> variable_positions;

  for (auto edge : dag_.GetEdges()) {
    for (const auto& [pos, mut] : edge.GetEdgeMutations()) {
      variable_positions.insert(pos.value);
    }
  }

  // Build variable_sites vector and site_to_index map
  dp_table_.variable_sites.clear();
  dp_table_.site_to_index.clear();

  size_t idx = 0;
  for (size_t pos : variable_positions) {
    dp_table_.variable_sites.push_back(MutationPosition{pos});
    dp_table_.site_to_index[pos] = idx++;
  }

  // Allocate DP tables
  size_t num_nodes = dag_.GetNodesCount();
  size_t num_edges = dag_.GetEdgesCount();
  size_t num_sites = dp_table_.variable_sites.size();

  dp_table_.dp_costs.resize(num_nodes);
  for (auto& node_costs : dp_table_.dp_costs) {
    node_costs.resize(num_sites);
    for (auto& site_costs : node_costs) {
      site_costs.fill(kSankoffInfinity);
    }
  }

  dp_table_.traceback.resize(num_edges);
  for (auto& edge_traceback : dp_table_.traceback) {
    edge_traceback.resize(num_sites);
    for (auto& site_traceback : edge_traceback) {
      site_traceback.fill(0);
    }
  }

  // Allocate ancestral bases storage
  ancestral_bases_.resize(num_nodes);
  for (auto& node_bases : ancestral_bases_) {
    node_bases.resize(num_sites, 0);
  }
}

template <typename DAG>
void SankoffScorer<DAG>::InitializeCostMatrices() {
  size_t num_sites = dp_table_.variable_sites.size();
  dp_table_.cost_matrices.resize(num_sites);

  SiteCostMatrix default_matrix = MakeUniformCostMatrix();

  for (size_t site_idx = 0; site_idx < num_sites; ++site_idx) {
    size_t pos = dp_table_.variable_sites[site_idx].value;

    if (uniform_cost_matrix_.has_value()) {
      // Use the single uniform matrix for all sites
      dp_table_.cost_matrices[site_idx] = *uniform_cost_matrix_;
    } else if (custom_cost_matrices_.count(pos) > 0) {
      // Use custom matrix for this position
      dp_table_.cost_matrices[site_idx] = custom_cost_matrices_.at(pos);
    } else {
      // Use default uniform matrix
      dp_table_.cost_matrices[site_idx] = default_matrix;
    }
  }
}

template <typename DAG>
MutationBase SankoffScorer<DAG>::GetLeafBase(typename DAG::NodeView leaf,
                                             MutationPosition pos) const {
  const auto& compact_genome = leaf.GetCompactGenome();
  auto maybe_base = compact_genome[pos];
  if (maybe_base.has_value()) {
    return *maybe_base;
  }
  // Not in compact genome means it matches the reference
  const std::string& ref_seq = dag_.GetReferenceSequence();
  // Position is 1-based
  return MutationBase{ref_seq.at(pos.value - 1)};
}

template <typename DAG>
void SankoffScorer<DAG>::ComputeNodeCosts(typename DAG::NodeView node) {
  size_t node_idx = node.GetId().value;
  size_t num_sites = dp_table_.variable_sites.size();

  if (node.IsLeaf()) {
    // Leaf node: cost 0 for observed base, infinity for others
    for (size_t site_idx = 0; site_idx < num_sites; ++site_idx) {
      MutationPosition pos = dp_table_.variable_sites[site_idx];
      MutationBase observed = GetLeafBase(node, pos);
      std::vector<size_t> compatible = GetCompatibleIndices(observed);

      // Initialize all to infinity
      dp_table_.dp_costs[node_idx][site_idx].fill(kSankoffInfinity);

      // Set compatible bases to 0
      for (size_t base_idx : compatible) {
        dp_table_.dp_costs[node_idx][site_idx][base_idx] = 0.0;
      }
    }
    return;
  }

  // Internal node: recursively compute children first, then aggregate
  // Initialize costs to 0 (will accumulate from children)
  for (size_t site_idx = 0; site_idx < num_sites; ++site_idx) {
    dp_table_.dp_costs[node_idx][site_idx].fill(0.0);
  }

  // Process each child edge
  for (auto clade : node.GetClades()) {
    for (EdgeId edge_id : clade) {
      auto edge = dag_.Get(edge_id);
      auto child = edge.GetChild();
      size_t edge_idx = edge_id.value;

      // Recursively compute child costs first
      ComputeNodeCosts(child);

      size_t child_idx = child.GetId().value;

      // For each site, compute contribution from this child
      for (size_t site_idx = 0; site_idx < num_sites; ++site_idx) {
        const SiteCostMatrix& cost_matrix = dp_table_.cost_matrices[site_idx];
        const SiteCosts& child_costs = dp_table_.dp_costs[child_idx][site_idx];

        // For each possible parent base, find optimal child base
        for (size_t parent_base = 0; parent_base < 4; ++parent_base) {
          double best_cost = kSankoffInfinity;
          uint8_t best_child_base = 0;

          for (size_t child_base = 0; child_base < 4; ++child_base) {
            double cost =
                cost_matrix[parent_base][child_base] + child_costs[child_base];
            if (cost < best_cost) {
              best_cost = cost;
              best_child_base = static_cast<uint8_t>(child_base);
            }
          }

          // Store traceback
          dp_table_.traceback[edge_idx][site_idx][parent_base] = best_child_base;

          // Add to node cost
          dp_table_.dp_costs[node_idx][site_idx][parent_base] += best_cost;
        }
      }
    }
  }
}

template <typename DAG>
double SankoffScorer<DAG>::ComputeScoreBelow(typename DAG::NodeView node) {
  // Run bottom-up DP
  ComputeNodeCosts(node);

  // Sum minimum costs across all sites at root
  size_t node_idx = node.GetId().value;
  size_t num_sites = dp_table_.variable_sites.size();

  total_score_ = 0.0;
  for (size_t site_idx = 0; site_idx < num_sites; ++site_idx) {
    const SiteCosts& costs = dp_table_.dp_costs[node_idx][site_idx];
    double min_cost = *std::min_element(costs.begin(), costs.end());
    if (min_cost < kSankoffInfinity) {
      total_score_ += min_cost;
    }
  }

  score_computed_ = true;
  return total_score_;
}

template <typename DAG>
void SankoffScorer<DAG>::TracebackFromNode(
    typename DAG::NodeView node, const std::vector<uint8_t>& assigned_bases) {
  size_t node_idx = node.GetId().value;
  size_t num_sites = dp_table_.variable_sites.size();

  // Store assigned bases for this node
  ancestral_bases_[node_idx] = assigned_bases;

  if (node.IsLeaf()) {
    return;
  }

  // Propagate to children using traceback
  for (auto clade : node.GetClades()) {
    for (EdgeId edge_id : clade) {
      auto edge = dag_.Get(edge_id);
      auto child = edge.GetChild();
      size_t edge_idx = edge_id.value;

      std::vector<uint8_t> child_bases(num_sites);
      for (size_t site_idx = 0; site_idx < num_sites; ++site_idx) {
        uint8_t parent_base = assigned_bases[site_idx];
        child_bases[site_idx] = dp_table_.traceback[edge_idx][site_idx][parent_base];
      }

      TracebackFromNode(child, child_bases);
    }
  }
}

template <typename DAG>
void SankoffScorer<DAG>::ReconstructAncestralSequences(typename DAG::NodeView root) {
  if (!score_computed_) {
    throw std::runtime_error(
        "Must call ComputeScoreBelow() before ReconstructAncestralSequences()");
  }

  size_t root_idx = root.GetId().value;
  size_t num_sites = dp_table_.variable_sites.size();

  // Choose optimal base at root for each site
  std::vector<uint8_t> root_bases(num_sites);
  for (size_t site_idx = 0; site_idx < num_sites; ++site_idx) {
    const SiteCosts& costs = dp_table_.dp_costs[root_idx][site_idx];
    size_t best_base = 0;
    double best_cost = costs[0];
    for (size_t base = 1; base < 4; ++base) {
      if (costs[base] < best_cost) {
        best_cost = costs[base];
        best_base = base;
      }
    }
    root_bases[site_idx] = static_cast<uint8_t>(best_base);
  }

  // Traceback from root
  TracebackFromNode(root, root_bases);

  ancestors_reconstructed_ = true;
}

template <typename DAG>
char SankoffScorer<DAG>::GetReconstructedBase(NodeId node_id,
                                              size_t site_index) const {
  if (!ancestors_reconstructed_) {
    throw std::runtime_error(
        "Must call ReconstructAncestralSequences() before GetReconstructedBase()");
  }
  return IndexToBase(ancestral_bases_[node_id.value][site_index]);
}

template <typename DAG>
CompactGenome SankoffScorer<DAG>::GetReconstructedGenome(
    typename DAG::NodeView node) const {
  if (!ancestors_reconstructed_) {
    throw std::runtime_error(
        "Must call ReconstructAncestralSequences() before GetReconstructedGenome()");
  }

  const std::string& ref_seq = dag_.GetReferenceSequence();
  size_t node_idx = node.GetId().value;
  size_t num_sites = dp_table_.variable_sites.size();

  ContiguousMap<MutationPosition, MutationBase> mutations;

  for (size_t site_idx = 0; site_idx < num_sites; ++site_idx) {
    MutationPosition pos = dp_table_.variable_sites[site_idx];
    char reconstructed = IndexToBase(ancestral_bases_[node_idx][site_idx]);
    char reference = ref_seq.at(pos.value - 1);  // 1-based position

    if (reconstructed != reference) {
      mutations[pos] = MutationBase{reconstructed};
    }
  }

  return CompactGenome{std::move(mutations)};
}
