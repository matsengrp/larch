// MatOptimizeScoringBackend implementation

template <typename DAG>
template <typename DAGView>
bool MatOptimizeScoringBackend<DAG>::Initialize(
    const DAGView& dag, const Profitable_Moves& move,
    const std::vector<Node_With_Major_Allele_Set_Change>& changes) {
  // Convert MAT::Node* to NodeId
  src_ = dag.GetNodeFromMAT(move.src);
  dst_ = dag.GetNodeFromMAT(move.dst);
  lca_ = dag.GetNodeFromMAT(move.LCA);
  move_score_change_ = move.score_change;

  // Build the changed Fitch set map, keyed by NodeId
  constexpr int end_sentinel = 2147483647;
  for (const auto& node_with_allele_set_change : changes) {
    if (not node_with_allele_set_change.node->is_leaf()) {
      Assert(node_with_allele_set_change.node != nullptr);

      ContiguousMap<MutationPosition, Mutation_Count_Change> node_map;
      for (const auto& mutation_count_change :
           node_with_allele_set_change.major_allele_set_change) {
        if (mutation_count_change.get_position() >= end_sentinel) {
          continue;
        }
        MutationPosition pos = {
            static_cast<size_t>(mutation_count_change.get_position())};
        node_map.insert({pos, mutation_count_change});
      }

      // Store by NodeId
      NodeId node_id = dag.GetNodeFromMAT(node_with_allele_set_change.node);
      changed_fitch_set_map_.insert({node_id, node_map.Copy()});

      // Also store by MATNodePtr for legacy access
      mat_keyed_fitch_set_map_.insert(
          {node_with_allele_set_change.node, std::move(node_map)});
    }
  }

  return true;
}

template <typename DAG>
ContiguousSet<MutationPosition>
MatOptimizeScoringBackend<DAG>::GetSitesWithScoringChanges(NodeId node) const {
  ContiguousSet<MutationPosition> result;
  auto it = changed_fitch_set_map_.find(node);
  if (it != changed_fitch_set_map_.end()) {
    for (auto& map_pair : it->second) {
      result.insert(map_pair.first);
    }
  }
  return result;
}

template <typename DAG>
bool MatOptimizeScoringBackend<DAG>::HasScoringChanges(NodeId node) const {
  return changed_fitch_set_map_.find(node) != changed_fitch_set_map_.end();
}

template <typename DAG>
template <typename DAGView>
std::pair<MAT::Mutations_Collection,
          std::optional<ContiguousMap<MutationPosition, Mutation_Count_Change>>>
MatOptimizeScoringBackend<DAG>::GetFitchSetParts(const DAGView& dag, NodeId node,
                                                 bool is_leaf, bool is_move_target,
                                                 bool is_move_new) const {
  auto mat_node = dag.Get(node).GetMATNode();
  if (is_leaf or is_move_target) {
    // if it's a leaf node, then the fitch sets don't change.
    // if it's the target node, then fitch sets can't have changed,
    // but the fitch sets recorded in tree_'s changed fitch set map
    // relative to this node are meant for this node's new parent!
    return {mat_node->mutations, std::nullopt};
  } else if (is_move_new) {
    // Then fitch set changes are relative to the target node
    auto result = mat_keyed_fitch_set_map_.find(dag.Get(dst_).GetMATNode());
    return {dag.Get(dst_).GetMATNode()->mutations,
            result == mat_keyed_fitch_set_map_.end()
                ? std::nullopt
                : std::make_optional(result->second.Copy())};
  } else {
    auto result = mat_keyed_fitch_set_map_.find(mat_node);
    return {mat_node->mutations, result == mat_keyed_fitch_set_map_.end()
                                     ? std::nullopt
                                     : std::make_optional(result->second.Copy())};
  }
}

template <typename DAG>
template <typename DAGView>
FitchSet MatOptimizeScoringBackend<DAG>::GetFitchSetAtSite(
    const DAGView& dag, NodeId node, MutationPosition site, bool is_leaf,
    bool is_move_target, bool is_move_new) const {
  Assert(site.value <= dag.GetReferenceSequence().size());
  auto [old_fitch_sets, changes] =
      GetFitchSetParts(dag, node, is_leaf, is_move_target, is_move_new);

  auto build_new_fitch_set = [](nuc_one_hot old_major_alleles,
                                nuc_one_hot old_boundary_alleles, nuc_one_hot rem_set,
                                nuc_one_hot add_set) -> FitchSet {
    Assert(old_major_alleles <= 16);
    Assert(old_boundary_alleles <= 16);
    Assert(rem_set <= 16);
    Assert(add_set <= 16);
    if ((old_major_alleles & add_set) != 0) {
      auto result = old_major_alleles & add_set;
      Assert(0 < result);
      Assert(result <= 16);
      return result;
    } else {
      auto result = (old_major_alleles & ~rem_set) | (old_boundary_alleles & add_set);
      Assert(result <= 16);
      if (result > 0) {
        return result;
      } else {
        result = old_major_alleles | (old_boundary_alleles & ~rem_set);
        return result;
      }
    }
  };

  if (old_fitch_sets.find(static_cast<int>(site.value)) ==
      old_fitch_sets.mutations.end()) {
    // if no fitch set is recorded on the corresponding MAT node, we can use
    // a singleton set containing the base at this site in the parent compact genome
    // In the special case that the node is the new move, its Fitch sets (and Fitch set
    // changes) should be calculated relative to the old parent of GetMoveTarget, and so
    // we will retrieve the base of that parent for the Fitch set in this case.

    nuc_one_hot boundary_allele =
        0;  // from the 3 conditions for empty fitch set (see Cheng's comment)

    // Get the parent base - need to handle is_move_new specially
    nuc_one_hot old_parent_base;
    if (is_move_new) {
      old_parent_base =
          base_to_singleton(dag.Get(dst_)
                                .GetOld()
                                .GetSingleParent()
                                .GetParent()
                                .GetCompactGenome()
                                .GetBase(site, dag.GetReferenceSequence()));
    } else {
      old_parent_base =
          base_to_singleton(dag.Get(node)
                                .GetOld()
                                .GetSingleParent()
                                .GetParent()
                                .GetCompactGenome()
                                .GetBase(site, dag.GetReferenceSequence()));
    }

    if (changes.has_value() and changes.value().Contains(site)) {
      return build_new_fitch_set(old_parent_base, boundary_allele,
                                 changes.value().at(site).get_decremented(),
                                 changes.value().at(site).get_incremented());
    } else {
      return FitchSet(old_parent_base);
    }
  } else if (changes.has_value() and changes.value().Contains(site)) {
    return FitchSet(build_new_fitch_set(
        old_fitch_sets.find(static_cast<int>(site.value))->get_all_major_allele(),
        old_fitch_sets.find(static_cast<int>(site.value))->get_boundary1_one_hot(),
        changes.value().at(site).get_decremented(),
        changes.value().at(site).get_incremented()));
  } else {
    return FitchSet(
        old_fitch_sets.find(static_cast<int>(site.value))->get_all_major_allele());
  }
}

template <typename DAG>
MutationBase MatOptimizeScoringBackend<DAG>::SelectBase(const FitchSet& fitch_set,
                                                        MutationBase old_base,
                                                        MutationBase parent_base) {
  // If parent base is in the Fitch set, prefer it
  if (fitch_set.find(parent_base.ToChar())) {
    return parent_base;
  }
  // Otherwise, if old base is in the Fitch set, keep it
  if (fitch_set.find(old_base.ToChar())) {
    return old_base;
  }
  // Otherwise, pick the first available base from the Fitch set
  return MutationBase{fitch_set.at(0)};
}

// MLScoringBackend implementation

#ifdef USE_NETAM

template <typename DAG>
template <typename DAGView>
std::string MLScoringBackend<DAG>::ExpandSequence(const DAGView& dag, NodeId node) {
  const auto& ref_seq = dag.GetReferenceSequence();
  std::string result{ref_seq.begin(), ref_seq.end()};

  const auto& compact_genome = dag.Get(node).GetCompactGenome();
  for (const auto& [pos, base] : compact_genome) {
    // MutationPosition is 1-indexed
    result[pos.value - 1] = base.ToChar();
  }

  return result;
}

template <typename DAG>
template <typename DAGView>
double MLScoringBackend<DAG>::ComputeEdgeLogLikelihood(const DAGView& dag,
                                                       NodeId parent,
                                                       NodeId child) const {
  if (model_ == nullptr) {
    return 0.0;
  }

  torch::NoGradGuard no_grad;

  // Expand sequences from compact genomes
  std::string parent_seq = ExpandSequence(dag, parent);
  std::string child_seq = ExpandSequence(dag, child);

  // Encode parent sequence for model input
  // encode_sequence returns (encoded_1d, wt_mod_2d) without batch dimension
  auto [encoded_1d, wt_mod_2d] = model_->encoder().encode_sequence(parent_seq);

  // Add batch dimension for model input
  auto encoded = encoded_1d.unsqueeze(0);  // [1, site_count]
  auto wt_mod = wt_mod_2d.unsqueeze(0);    // [1, site_count, 4]

  // Create mask (all true for valid positions)
  auto mask = torch::ones({1, encoded.size(1)}, torch::kBool);

  // Run model forward pass
  auto [rates, csp_logits] = (*model_)->forward(encoded, mask, wt_mod);

  // Apply softmax to get CSP probabilities
  auto csp = torch::softmax(csp_logits, -1);

  // Encode sequences as base indices
  auto parent_indices = netam::kmer_sequence_encoder::encode_bases(parent_seq);
  auto child_indices = netam::kmer_sequence_encoder::encode_bases(child_seq);

  // Compute log-likelihood
  auto log_likelihood =
      netam::poisson_context_log_likelihood(rates, csp, parent_indices, child_indices);

  return log_likelihood.item<double>();
}

#endif  // USE_NETAM

template <typename DAG>
template <typename DAGView>
bool MLScoringBackend<DAG>::Initialize([[maybe_unused]] const DAGView& dag, NodeId src,
                                       NodeId dst, NodeId lca) {
  src_ = src;
  dst_ = dst;
  lca_ = lca;
  ml_score_change_ = 0.0;

#ifdef USE_NETAM
  if (model_ == nullptr) {
    return true;
  }

  // Compute total log-likelihood for all edges in the affected subtree
  // Starting from LCA, traverse down and sum edge log-likelihoods
  //
  // For the initial implementation, we compute the log-likelihood of the
  // new tree configuration. A proper score "change" would require comparing
  // old vs new, but that requires access to the old tree state.

  double total_ll = 0.0;

  // Use a worklist to traverse from LCA down
  std::vector<NodeId> worklist;
  worklist.push_back(lca);

  while (!worklist.empty()) {
    NodeId current = worklist.back();
    worklist.pop_back();

    auto node = dag.Get(current);
    if (node.IsLeaf()) {
      continue;
    }

    // Process all children of this node
    for (auto clade : node.GetClades()) {
      for (auto child_edge : clade) {
        NodeId child_id = child_edge.GetChild().GetId();

        // Skip UA node
        if (child_edge.GetChild().IsUA()) {
          continue;
        }

        // Compute log-likelihood for this edge
        double edge_ll = ComputeEdgeLogLikelihood(dag, current, child_id);
        total_ll += edge_ll;

        // Add child to worklist for further traversal
        worklist.push_back(child_id);
      }
    }
  }

  // Store negated log-likelihood (lower is better, consistent with parsimony)
  ml_score_change_ = -total_ll;

#endif  // USE_NETAM

  return true;
}

template <typename DAG>
template <typename DAGView>
typename MLScoringBackend<DAG>::NucleotideSet
MLScoringBackend<DAG>::GetNucleotideSetAtSite(const DAGView& dag, NodeId node,
                                              MutationPosition site, bool, bool,
                                              bool) const {
  // For ML backend, just return the base from the compact genome
  auto base =
      dag.Get(node).GetCompactGenome().GetBase(site, dag.GetReferenceSequence());
  return NucleotideSet{base.ToChar()};
}

// ParsimonyOnlyScoringBackend implementation

template <typename DAG>
template <typename TreeView>
NodeId ParsimonyOnlyScoringBackend<DAG>::FindTreeRoot(const TreeView& tree) const {
  auto ua = tree.GetRoot();
  for (auto clade : ua.GetClades()) {
    for (auto edge : clade) {
      return edge.GetChild().GetId();
    }
  }
  Fail("No tree root found");
}

template <typename DAG>
template <typename TreeView>
void ParsimonyOnlyScoringBackend<DAG>::FitchVisit(
    const TreeView& tree, NodeId node_id,
    std::unordered_map<size_t, std::vector<uint8_t>>& fitch_sets, int& score) const {
  auto node = tree.Get(node_id);
  auto& node_sets = fitch_sets[node_id.value];
  node_sets.resize(variable_sites_.size());

  if (node.IsLeaf()) {
    for (size_t i = 0; i < variable_sites_.size(); i++) {
      auto base =
          node.GetCompactGenome().GetBase(variable_sites_[i], tree.GetReferenceSequence());
      node_sets[i] = static_cast<uint8_t>(base_to_singleton(base));
    }
  } else {
    // Visit children first (postorder)
    std::vector<NodeId> children;
    for (auto clade : node.GetClades()) {
      for (auto edge : clade) {
        auto child_id = edge.GetChild().GetId();
        FitchVisit(tree, child_id, fitch_sets, score);
        children.push_back(child_id);
      }
    }

    // Compute Fitch sets from children
    for (size_t i = 0; i < variable_sites_.size(); i++) {
      uint8_t intersection = fitch_sets[children[0].value][i];
      uint8_t union_set = fitch_sets[children[0].value][i];
      for (size_t c = 1; c < children.size(); c++) {
        uint8_t child_set = fitch_sets[children[c].value][i];
        intersection &= child_set;
        union_set |= child_set;
      }
      if (intersection != 0) {
        node_sets[i] = intersection;
      } else {
        node_sets[i] = union_set;
        score++;
      }
    }
  }
}

template <typename DAG>
template <typename DAGView>
bool ParsimonyOnlyScoringBackend<DAG>::Initialize(const DAGView& dag, NodeId src,
                                                   NodeId dst, NodeId lca) {
  src_ = src;
  dst_ = dst;
  lca_ = lca;

  // Collect variable sites from the original tree's edge mutations
  auto original = dag.GetOriginal();
  std::set<size_t> var_sites_set;
  for (auto edge : original.GetEdges()) {
    for (auto& [pos, mut] : edge.GetEdgeMutations()) {
      var_sites_set.insert(pos.value);
    }
  }
  variable_sites_.clear();
  site_to_index_.clear();
  size_t idx = 0;
  for (auto site : var_sites_set) {
    variable_sites_.push_back({site});
    site_to_index_[site] = idx++;
  }

  if (variable_sites_.empty()) {
    score_change_ = 0;
    return true;
  }

  // Run Fitch on original tree (pre-move)
  int old_score = 0;
  std::unordered_map<size_t, std::vector<uint8_t>> old_fitch_sets;
  NodeId old_root = FindTreeRoot(original);
  FitchVisit(original, old_root, old_fitch_sets, old_score);

  // Run Fitch on overlay tree (post-move)
  // Use .Const() to read through the const overlay path, which falls through
  // to the target for non-overlaid features (avoids "Can't modify non-overlaid
  // node" error on mutable overlays).
  int new_score = 0;
  new_fitch_sets_.clear();
  auto const_dag = dag.Const();
  NodeId new_root = FindTreeRoot(const_dag);
  FitchVisit(const_dag, new_root, new_fitch_sets_, new_score);

  // Score change = new parsimony - old parsimony (negative = improvement)
  score_change_ = new_score - old_score;

  // Find changed sites per node
  changed_sites_map_ = {};
  for (auto& [node_val, new_sets] : new_fitch_sets_) {
    NodeId node_id{node_val};
    auto old_it = old_fitch_sets.find(node_val);
    if (old_it == old_fitch_sets.end()) {
      // New node (appended by SPR): all sites are "changed"
      ContiguousSet<MutationPosition> sites;
      for (auto& site : variable_sites_) {
        sites.insert(site);
      }
      changed_sites_map_.insert({node_id, std::move(sites)});
    } else {
      ContiguousSet<MutationPosition> sites;
      for (size_t i = 0; i < variable_sites_.size(); i++) {
        if (new_sets[i] != old_it->second[i]) {
          sites.insert(variable_sites_[i]);
        }
      }
      if (not sites.empty()) {
        changed_sites_map_.insert({node_id, std::move(sites)});
      }
    }
  }

  return true;
}

template <typename DAG>
ContiguousSet<MutationPosition>
ParsimonyOnlyScoringBackend<DAG>::GetSitesWithScoringChanges(NodeId node) const {
  auto it = changed_sites_map_.find(node);
  if (it != changed_sites_map_.end()) {
    return it->second.Copy();
  }
  return {};
}

template <typename DAG>
bool ParsimonyOnlyScoringBackend<DAG>::HasScoringChanges(NodeId node) const {
  auto it = changed_sites_map_.find(node);
  return it != changed_sites_map_.end();
}

template <typename DAG>
template <typename DAGView>
std::pair<MAT::Mutations_Collection,
          std::optional<ContiguousMap<MutationPosition, Mutation_Count_Change>>>
ParsimonyOnlyScoringBackend<DAG>::GetFitchSetParts(const DAGView&, NodeId node, bool,
                                                    bool, bool) const {
  // Return empty Mutations_Collection and a changes map with keys = changed positions
  ContiguousMap<MutationPosition, Mutation_Count_Change> changes;
  auto it = changed_sites_map_.find(node);
  if (it != changed_sites_map_.end()) {
    for (auto site : it->second) {
      changes.insert({site, Mutation_Count_Change{}});
    }
  }
  return {MAT::Mutations_Collection{},
          changes.empty() ? std::nullopt : std::make_optional(std::move(changes))};
}

template <typename DAG>
template <typename DAGView>
FitchSet ParsimonyOnlyScoringBackend<DAG>::GetFitchSetAtSite(
    const DAGView& dag, NodeId node, MutationPosition site, bool is_leaf, bool,
    bool) const {
  if (is_leaf) {
    // For leaves, return the observed base as a singleton
    auto base =
        dag.GetOriginal().Get(node).GetCompactGenome().GetBase(
            site, dag.GetReferenceSequence());
    return FitchSet{static_cast<int>(static_cast<uint8_t>(base_to_singleton(base)))};
  }

  // Look up in the pre-computed Fitch sets
  auto node_it = new_fitch_sets_.find(node.value);
  if (node_it != new_fitch_sets_.end()) {
    auto site_it = site_to_index_.find(site.value);
    if (site_it != site_to_index_.end() && site_it->second < node_it->second.size()) {
      return FitchSet{static_cast<int>(node_it->second[site_it->second])};
    }
  }

  // Site not in variable sites: return the reference base
  const auto& ref = dag.GetReferenceSequence();
  Assert(site.value > 0 and site.value <= ref.size());
  return FitchSet{static_cast<int>(static_cast<uint8_t>(
      base_to_singleton(MutationBase{ref[site.value - 1]})))};
}

template <typename DAG>
MutationBase ParsimonyOnlyScoringBackend<DAG>::SelectBase(const FitchSet& fitch_set,
                                                           MutationBase old_base,
                                                           MutationBase parent_base) {
  // Same logic as MatOptimizeScoringBackend
  if (fitch_set.find(parent_base.ToChar())) {
    return parent_base;
  }
  if (fitch_set.find(old_base.ToChar())) {
    return old_base;
  }
  return MutationBase{fitch_set.at(0)};
}

template <typename DAG>
const ContiguousMap<MATNodePtr, ContiguousMap<MutationPosition, Mutation_Count_Change>>&
ParsimonyOnlyScoringBackend<DAG>::GetChangedFitchSetMap() const {
  static ContiguousMap<MATNodePtr, ContiguousMap<MutationPosition, Mutation_Count_Change>>
      empty;
  return empty;
}
