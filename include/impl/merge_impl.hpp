const CompactGenome* CompactGenome::Empty() {
  static CompactGenome empty = {};
  return &empty;
}

CompactGenome::CompactGenome(const Mutations& mutations, const CompactGenome& parent,
                             std::string_view reference_sequence)
    : mutations_{[&] {
        std::vector<std::pair<MutationPosition, char>> result{parent.mutations_.begin(),
                                                              parent.mutations_.end()};
        for (auto [pos, base] : mutations) {
          const bool is_valid = base != reference_sequence.at(pos.value - 1);
          auto it =
              std::lower_bound(result.begin(), result.end(), pos,
                               [](std::pair<MutationPosition, char> lhs,
                                  MutationPosition rhs) { return lhs.first < rhs; });
          if (it != result.end() and it->first == pos) {
            if (is_valid) {
              it->second = base;
            } else {
              result.erase(it);
            }
          } else {
            if (is_valid) {
              result.insert(it, {pos, base});
            }
          }
        }
        return result;
      }()},
      hash_{ComputeHash(mutations_)} {}

CompactGenome::CompactGenome(std::vector<std::pair<MutationPosition, char>>&& mutations)
    : mutations_{mutations}, hash_{ComputeHash(mutations_)} {}

bool CompactGenome::operator==(const CompactGenome& rhs) const noexcept {
  if (hash_ != rhs.hash_) return false;
  return mutations_ == rhs.mutations_;
}

size_t CompactGenome::Hash() const noexcept { return hash_; }

std::optional<char> CompactGenome::operator[](MutationPosition pos) const {
  auto it = std::lower_bound(mutations_.begin(), mutations_.end(), pos,
                             [](std::pair<MutationPosition, char> lhs,
                                MutationPosition rhs) { return lhs.first < rhs; });
  if (it != mutations_.end() and it->first == pos) {
    return it->second;
  } else {
    return std::nullopt;
  }
}

auto CompactGenome::begin() const { return mutations_.begin(); }

auto CompactGenome::end() const { return mutations_.end(); }

bool CompactGenome::empty() const { return mutations_.empty(); }

Mutations CompactGenome::ToEdgeMutations(std::string_view reference_sequence,
                                         const CompactGenome& parent,
                                         const CompactGenome& child) {
  Mutations result;
  for (auto [pos, child_base] : child) {
    char parent_base = reference_sequence.at(pos.value - 1);
    auto opt_parent_base = parent[pos];
    if (opt_parent_base.has_value()) {
      parent_base = opt_parent_base.value();
    }
    if (parent_base != child_base) {
      result[pos] = child_base;
    }
  }

  for (auto [pos, parent_base] : parent) {
    char child_base = reference_sequence.at(pos.value - 1);
    auto opt_child_base = child[pos];
    if (opt_child_base.has_value()) {
      child_base = opt_child_base.value();
    }
    if (child_base != parent_base) {
      result[pos] = child_base;
    }
  }
  return result;
}

size_t CompactGenome::ComputeHash(
    const std::vector<std::pair<MutationPosition, char>>& mutations) {
  size_t result = 0;
  for (auto [pos, base] : mutations) {
    result = HashCombine(result, pos.value);
    result = HashCombine(result, base);
  }
  return result;
}

const LeafSet* LeafSet::Empty() {
  static LeafSet empty = {};
  return &empty;
}

LeafSet::LeafSet(Node node, const std::vector<NodeLabel>& labels,
                 std::vector<LeafSet>& computed_leafsets)
    : clades_{[&] {
        std::vector<std::vector<const CompactGenome*>> clades;
        clades.reserve(node.GetClades().size());
        for (auto clade : node.GetClades()) {
          std::vector<const CompactGenome*> clade_leafs;
          clade_leafs.reserve(clade.size());
          for (Node child : clade | Transform::GetChild()) {
            const LeafSet& child_leaf_set = computed_leafsets.at(child.GetId().value);
            if (child.IsLeaf()) {
              clade_leafs.push_back(labels.at(child.GetId().value).compact_genome);
            } else {
              for (auto& child_leafs :
                   computed_leafsets.at(child.GetId().value).clades_) {
                clade_leafs.insert(clade_leafs.end(), child_leafs.begin(),
                                   child_leafs.end());
              }
            }
          }
          clade_leafs |= ranges::actions::sort | ranges::actions::unique;
          clades.emplace_back(std::move(clade_leafs));
        }
        clades |= ranges::actions::sort;
        return clades;
      }()},
      hash_{ComputeHash(clades_)} {}

LeafSet::LeafSet(std::vector<std::vector<const CompactGenome*>>&& clades)
    : clades_{clades}, hash_{ComputeHash(clades_)} {}

bool LeafSet::operator==(const LeafSet& rhs) const noexcept {
  return clades_ == rhs.clades_;
}

size_t LeafSet::Hash() const noexcept { return hash_; }

auto LeafSet::begin() const { return clades_.begin(); }

auto LeafSet::end() const { return clades_.end(); }

bool LeafSet::empty() const { return clades_.empty(); }

size_t LeafSet::ComputeHash(
    const std::vector<std::vector<const CompactGenome*>>& clades) {
  size_t hash = 0;
  for (auto& clade : clades) {
    for (auto leaf : clade) {
      hash = HashCombine(hash, leaf->Hash());
    }
  }
  return hash;
}

NodeLabel::NodeLabel()
    : compact_genome{CompactGenome::Empty()}, leaf_set{LeafSet::Empty()} {}

NodeLabel::NodeLabel(const CompactGenome* cg, const LeafSet* ls)
    : compact_genome{cg}, leaf_set{ls} {}

bool NodeLabel::operator==(const NodeLabel& rhs) const noexcept {
  return compact_genome == rhs.compact_genome && leaf_set == rhs.leaf_set;
}

size_t NodeLabel::Hash() const noexcept {
  return HashCombine(reinterpret_cast<std::uintptr_t>(compact_genome),
                     reinterpret_cast<std::uintptr_t>(leaf_set));
}

bool EdgeLabel::operator==(const EdgeLabel& rhs) const noexcept {
  return parent_compact_genome == rhs.parent_compact_genome &&
         parent_leaf_set == rhs.parent_leaf_set &&
         child_compact_genome == rhs.child_compact_genome &&
         child_leaf_set == rhs.child_leaf_set;
}

size_t EdgeLabel::Hash() const noexcept {
  size_t hash = reinterpret_cast<std::uintptr_t>(parent_compact_genome);
  hash = HashCombine(hash, reinterpret_cast<std::uintptr_t>(parent_leaf_set));
  hash = HashCombine(hash, reinterpret_cast<std::uintptr_t>(child_compact_genome));
  hash = HashCombine(hash, reinterpret_cast<std::uintptr_t>(child_leaf_set));
  return hash;
}

Merge::Merge(std::string_view reference_sequence)
    : reference_sequence_{reference_sequence} {}

void Merge::AddTrees(const std::vector<std::reference_wrapper<const DAG>>& trees,
                     const std::vector<std::vector<Mutations>>& mutations,
                     bool show_progress) {
  std::vector<size_t> tree_idxs;
  tree_idxs.resize(trees.size());
  std::iota(tree_idxs.begin(), tree_idxs.end(), trees_.size());

  trees_.insert(trees_.end(), trees.begin(), trees.end());
  mutations_.insert(mutations_.end(), mutations.begin(), mutations.end());
  tree_labels_.resize(trees_.size());

  ComputeCompactGenomes(tree_idxs, show_progress);
  ComputeLeafSets(tree_idxs, show_progress);
  MergeTrees(tree_idxs);
  Assert(result_nodes_.size() == result_dag_.GetNodes().size());
  Assert(result_edges_.size() == result_dag_.GetEdges().size());
  result_dag_.BuildConnections();
}

void Merge::AddDAGs(const std::vector<std::reference_wrapper<const DAG>>& dags,
                    std::vector<std::vector<CompactGenome>>&& compact_genomes,
                    bool show_progress) {
  std::vector<size_t> tree_idxs;
  tree_idxs.resize(dags.size());
  std::iota(tree_idxs.begin(), tree_idxs.end(), trees_.size());

  trees_.insert(trees_.end(), dags.begin(), dags.end());
  tree_labels_.resize(trees_.size());

  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const DAG& tree = trees_.at(tree_idx).get();
    std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    labels.resize(tree.GetNodes().size());
    for (size_t node_idx = 0; node_idx < tree.GetNodes().size(); ++node_idx) {
      auto cg_iter = all_compact_genomes_.insert(
          std::move(compact_genomes.at(tree_idx).at(node_idx)));
      labels.at(node_idx).compact_genome = std::addressof(*cg_iter.first);
    }
  });

  ComputeLeafSets(tree_idxs, show_progress);
  MergeTrees(tree_idxs);
  Assert(result_nodes_.size() == result_dag_.GetNodes().size());
  Assert(result_edges_.size() == result_dag_.GetEdges().size());
  result_dag_.BuildConnections();
}

DAG& Merge::GetResult() { return result_dag_; }

const DAG& Merge::GetResult() const { return result_dag_; }

const std::unordered_map<NodeLabel, NodeId>& Merge::GetResultNodes() const {
  return result_nodes_;
}

std::vector<Mutations> Merge::ComputeResultEdgeMutations() const {
  std::vector<Mutations> result;
  result.resize(result_dag_.GetEdges().size());
  for (auto& [label, edge_id] : result_edges_) {
    Assert(label.parent_compact_genome);
    Assert(label.child_compact_genome);
    const CompactGenome& parent = *label.parent_compact_genome;
    const CompactGenome& child = *label.child_compact_genome;
    Mutations& muts = result.at(edge_id.value);
    muts = CompactGenome::ToEdgeMutations(reference_sequence_, parent, child);
  }
  return result;
}

void Merge::ComputeCompactGenomes(const std::vector<size_t>& tree_idxs,
                                  bool show_progress) {
  if (show_progress) {
    std::cout << "Computing compact genomes " << std::flush;
  }
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const DAG& tree = trees_.at(tree_idx).get();
    const std::vector<Mutations>& edge_mutations = mutations_.at(tree_idx);
    std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    labels.resize(tree.GetNodes().size());
    std::vector<CompactGenome> computed_cgs =
        ComputeCompactGenomesDAG(tree, edge_mutations, reference_sequence_);

    // TODO move to a test case
    for (Edge edge : tree.GetEdges()) {
      const CompactGenome& parent = computed_cgs.at(edge.GetParentId().value);
      const CompactGenome& child = computed_cgs.at(edge.GetChildId().value);
      Mutations mutations =
          CompactGenome::ToEdgeMutations(reference_sequence_, parent, child);
      Assert(mutations == edge_mutations.at(edge.GetId().value));
    }

    for (size_t node_idx = 0; node_idx < tree.GetNodes().size(); ++node_idx) {
      auto cg_iter = all_compact_genomes_.insert(std::move(computed_cgs.at(node_idx)));
      labels.at(node_idx).compact_genome = std::addressof(*cg_iter.first);
    }
    if (show_progress) {
      std::cout << "." << std::flush;
    }
  });
  if (show_progress) {
    std::cout << " done.\n";
  }
}

void Merge::ComputeLeafSets(const std::vector<size_t>& tree_idxs, bool show_progress) {
  if (show_progress) {
    std::cout << "Computing leaf sets " << std::flush;
  }
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const DAG& tree = trees_.at(tree_idx).get();
    std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    std::vector<LeafSet> computed_ls = ComputeLeafSetsDAG(tree, labels);
    for (size_t node_idx = 0; node_idx < tree.GetNodes().size(); ++node_idx) {
      auto ls_iter = all_leaf_sets_.insert(std::move(computed_ls.at(node_idx)));
      labels.at(node_idx).leaf_set = std::addressof(*ls_iter.first);
    }
    if (show_progress) {
      std::cout << "." << std::flush;
    }
  });
  if (show_progress) {
    std::cout << " done.\n";
  }
}

void Merge::MergeTrees(const std::vector<size_t>& tree_idxs) {
  NodeId node_id{result_dag_.GetNodes().size()};
  std::mutex mtx;
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    for (auto label : labels) {
      std::unique_lock<std::mutex> lock{mtx};
      if (result_nodes_.try_emplace(label, node_id).second) {
        ++node_id.value;
      }
    }
  });
  tbb::concurrent_vector<std::pair<EdgeLabel, EdgeId>> added_edges;
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const DAG& tree = trees_.at(tree_idx).get();
    const std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    EdgeId edge_id{result_dag_.GetEdges().size()};
    for (Edge edge : tree.GetEdges()) {
      auto& parent_label = labels.at(edge.GetParentId().value);
      auto& child_label = labels.at(edge.GetChildId().value);
      auto ins =
          result_edges_.insert({{parent_label.compact_genome, parent_label.leaf_set,
                                 child_label.compact_genome, child_label.leaf_set},
                                {}});
      if (ins.second) {
        ins.first->second = edge_id;
        edge_id.value++;
        added_edges.push_back(*ins.first);
      }
    }
  });
  result_dag_.InitializeNodes(result_nodes_.size());
  EdgeId edge_id{result_dag_.GetEdges().size()};
  for (auto& [edge, id] : added_edges) {
    auto parent =
        result_nodes_.find(NodeLabel{edge.parent_compact_genome, edge.parent_leaf_set});
    auto child =
        result_nodes_.find(NodeLabel{edge.child_compact_genome, edge.child_leaf_set});
    Assert(parent != result_nodes_.end());
    Assert(child != result_nodes_.end());
    Assert(parent->second.value < result_dag_.GetNodes().size());
    Assert(child->second.value < result_dag_.GetNodes().size());
    id = edge_id;
    result_dag_.AddEdge(edge_id, parent->second, child->second, {0});
    edge_id.value++;
  }
}

std::vector<CompactGenome> Merge::ComputeCompactGenomes(
    const DAG& tree, const std::vector<Mutations>& edge_mutations,
    std::string_view reference_sequence) {
  std::vector<CompactGenome> result;
  result.resize(tree.GetNodes().size());
  for (auto [node, edge] : tree.TraversePreOrder()) {
    const Mutations& mutations = edge_mutations.at(edge.GetId().value);
    const CompactGenome& parent = result.at(edge.GetParentId().value);
    CompactGenome& compact_genome = result.at(node.GetId().value);
    compact_genome = CompactGenome{mutations, parent, reference_sequence};
  }
  return result;
}

std::vector<CompactGenome> Merge::ComputeCompactGenomesDAG(
    const DAG& tree, const std::vector<Mutations>& edge_mutations,
    std::string_view reference_sequence) {
  std::vector<CompactGenome> result;
  result.resize(tree.GetNodes().size());
  auto ComputeCG = [&](auto& self, Node node) {
    if (node.IsRoot()) {
      return;
    }
    CompactGenome& compact_genome = result.at(node.GetId().value);
    if (not compact_genome.empty()) {
      return;
    }
    Edge edge = *node.GetParents().begin();
    self(self, edge.GetParent());
    const Mutations& mutations = edge_mutations.at(edge.GetId().value);
    const CompactGenome& parent = result.at(edge.GetParentId().value);
    compact_genome = CompactGenome{mutations, parent, reference_sequence};
  };
  for (Node node : tree.GetNodes()) {
    ComputeCG(ComputeCG, node);
  }
  return result;
}

std::vector<LeafSet> Merge::ComputeLeafSets(const DAG& tree,
                                            const std::vector<NodeLabel>& labels) {
  std::vector<LeafSet> result;
  result.resize(tree.GetNodes().size());
  for (Node node : tree.TraversePostOrder()) {
    result.at(node.GetId().value) = LeafSet{node, labels, result};
  }
  return result;
}

std::vector<LeafSet> Merge::ComputeLeafSetsDAG(const DAG& tree,
                                               const std::vector<NodeLabel>& labels) {
  std::vector<LeafSet> result;
  result.resize(tree.GetNodes().size());
  auto ComputeLS = [&](auto& self, Node node) {
    LeafSet& leaf_set = result.at(node.GetId().value);
    if (not leaf_set.empty()) {
      return;
    }
    for (Node child : node.GetChildren() | Transform::GetChild()) {
      self(self, child);
    }
    result.at(node.GetId().value) = LeafSet{node, labels, result};
  };
  for (Node node : tree.GetNodes()) {
    ComputeLS(ComputeLS, node);
  }
  return result;
}
