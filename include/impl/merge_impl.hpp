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
          for (Node child : clade | ranges::views::transform(Transform::GetChild)) {
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
        clades |= ranges::actions::sort | ranges::actions::unique;
        return clades;
      }()},
      hash_{ComputeHash(clades_)} {}

LeafSet::LeafSet(std::vector<std::vector<const CompactGenome*>>&& clades)
    : clades_{clades}, hash_{ComputeHash(clades_)} {}

bool LeafSet::operator==(const LeafSet& rhs) const noexcept {
  if (hash_ != rhs.hash_) return false;
  return clades_ == rhs.clades_;
}

size_t LeafSet::Hash() const noexcept { return hash_; }

size_t LeafSet::ComputeHash(
    const std::vector<std::vector<const CompactGenome*>>& clades) {
  size_t hash = 0;
  for (auto& clade : clades) {
    for (auto& leaf : clade) {
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
  size_t hash = 0;
  hash = HashCombine(hash, reinterpret_cast<std::uintptr_t>(parent_compact_genome));
  hash = HashCombine(hash, reinterpret_cast<std::uintptr_t>(parent_leaf_set));
  hash = HashCombine(hash, reinterpret_cast<std::uintptr_t>(child_compact_genome));
  hash = HashCombine(hash, reinterpret_cast<std::uintptr_t>(child_leaf_set));
  return hash;
}

Merge::Merge(std::string_view reference_sequence, const std::vector<HistoryDAG>& trees,
             const std::vector<std::vector<Mutations>>& mutations, bool show_progress)
    : reference_sequence_{reference_sequence},
      trees_{trees},
      mutations_{mutations},
      show_progress_{show_progress} {
  tree_labels_.resize(trees_.size());
}

void Merge::Run() {
  ComputeCompactGenomes();
  ComputeLeafSets();
  MergeTrees();
  result_.BuildConnections();
}

HistoryDAG& Merge::GetResult() { return result_; }

const HistoryDAG& Merge::GetResult() const { return result_; }

const std::vector<std::vector<NodeLabel>>& Merge::GetTreeLabels() const {
  return tree_labels_;
}

const std::unordered_map<NodeLabel, NodeId>& Merge::GetResultNodes() const {
  return result_nodes_;
}

const ConcurrentUnorderedSet<EdgeLabel>& Merge::GetResultEdges() const {
  return result_edges_;
}

void Merge::ComputeCompactGenomes() {
  std::vector<size_t> tree_idxs;
  tree_idxs.resize(trees_.size());
  std::iota(tree_idxs.begin(), tree_idxs.end(), 0);

  if (show_progress_) {
    std::cout << "Computing compact genomes " << std::flush;
  }
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const HistoryDAG& tree = trees_.at(tree_idx);
    const std::vector<Mutations>& edge_mutations = mutations_.at(tree_idx);
    std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    labels.resize(tree.GetNodes().size());
    std::vector<CompactGenome> computed_cgs =
        ComputeCompactGenomes(tree, edge_mutations, reference_sequence_);
    for (size_t node_idx = 0; node_idx < tree.GetNodes().size(); ++node_idx) {
      auto cg_iter = all_compact_genomes_.insert(std::move(computed_cgs.at(node_idx)));
      labels.at(node_idx).compact_genome = std::addressof(*cg_iter.first);
    }
    if (show_progress_) {
      std::cout << "." << std::flush;
    }
  });
  if (show_progress_) {
    std::cout << " done.\n";
  }
}

void Merge::ComputeLeafSets() {
  std::vector<size_t> tree_idxs;
  tree_idxs.resize(trees_.size());
  std::iota(tree_idxs.begin(), tree_idxs.end(), 0);

  if (show_progress_) {
    std::cout << "Computing leaf sets " << std::flush;
  }
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const HistoryDAG& tree = trees_.at(tree_idx);
    std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    std::vector<LeafSet> computed_ls = ComputeLeafSets(tree, labels);
    for (size_t node_idx = 0; node_idx < tree.GetNodes().size(); ++node_idx) {
      auto ls_iter = all_leaf_sets_.insert(std::move(computed_ls.at(node_idx)));
      labels.at(node_idx).leaf_set = std::addressof(*ls_iter.first);
    }
    if (show_progress_) {
      std::cout << "." << std::flush;
    }
  });
  if (show_progress_) {
    std::cout << " done.\n";
  }
}

void Merge::MergeTrees() {
  std::vector<size_t> tree_idxs;
  tree_idxs.resize(trees_.size());
  std::iota(tree_idxs.begin(), tree_idxs.end(), 0);

  size_t node_id = 0;
  std::mutex mtx;
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    for (auto label : labels) {
      std::unique_lock<std::mutex> lock{mtx};
      auto i = result_nodes_.find(label);
      if (i == result_nodes_.end()) {
        result_nodes_.insert({label, {node_id++}});
      }
    }
  });
  tbb::parallel_for_each(tree_idxs.begin(), tree_idxs.end(), [&](size_t tree_idx) {
    const HistoryDAG& tree = trees_.at(tree_idx);
    const std::vector<NodeLabel>& labels = tree_labels_.at(tree_idx);
    for (Edge edge : tree.GetEdges()) {
      auto& parent_label = labels.at(edge.GetParent().GetId().value);
      auto& child_label = labels.at(edge.GetChild().GetId().value);
      result_edges_.insert({parent_label.compact_genome, parent_label.leaf_set,
                            child_label.compact_genome, child_label.leaf_set});
    }
  });
  result_.InitializeComponents(result_nodes_.size(), result_edges_.size());
  size_t edge_id = 0;
  for (auto& edge : result_edges_) {
    auto parent =
        result_nodes_.find(NodeLabel{edge.parent_compact_genome, edge.parent_leaf_set});
    auto child =
        result_nodes_.find(NodeLabel{edge.child_compact_genome, edge.child_leaf_set});
    Assert(parent != result_nodes_.end());
    Assert(child != result_nodes_.end());
    Assert(parent->second.value != NoId);
    Assert(child->second.value != NoId);
    result_.AddEdge({edge_id++}, parent->second, child->second, {0});
  }
}

std::vector<CompactGenome> Merge::ComputeCompactGenomes(
    const HistoryDAG& tree, const std::vector<Mutations>& edge_mutations,
    std::string_view reference_sequence) {
  std::vector<CompactGenome> result;
  result.resize(tree.GetNodes().size());
  for (auto iter : tree.TraversePreOrder()) {
    const Mutations& mutations = edge_mutations.at(iter.GetEdge().GetId().value);
    const CompactGenome& parent = result.at(iter.GetEdge().GetParent().GetId().value);
    CompactGenome& compact_genome = result.at(iter.GetNode().GetId().value);
    compact_genome = CompactGenome{mutations, parent, reference_sequence};
  }
  return result;
}

std::vector<LeafSet> Merge::ComputeLeafSets(const HistoryDAG& tree,
                                            const std::vector<NodeLabel>& labels) {
  std::vector<LeafSet> result;
  result.resize(tree.GetNodes().size());
  for (Node node : tree.TraversePostOrder()) {
    result.at(node.GetId().value) = LeafSet{node, labels, result};
  }
  return result;
}
