const LeafSet* LeafSet::Empty() {
  static const LeafSet empty = {};
  return &empty;
}

template <typename Node>
LeafSet::LeafSet(Node node, const std::vector<NodeLabel>& labels,
                 std::vector<LeafSet>& computed_leafsets)
    : clades_{[&] {
        std::vector<std::vector<UniqueData>> clades;
        clades.reserve(node.GetCladesCount());
        for (auto clade : node.GetClades()) {
          std::vector<UniqueData> clade_leafs;
          clade_leafs.reserve(clade.size());
          for (Node child : clade | Transform::GetChild()) {
            if (child.IsLeaf()) {
              Assert(child.Const().HaveSampleId());
              UniqueData id = labels.at(child.GetId().value).GetSampleId();
              clade_leafs.push_back(id);
            } else {
              for (auto& child_leafs :
                   computed_leafsets.at(child.GetId().value).clades_) {
                clade_leafs.insert(clade_leafs.end(), child_leafs.begin(),
                                   child_leafs.end());
              }
            }
          }
          clade_leafs |= ranges::actions::sort | ranges::actions::unique;
          Assert(not clade_leafs.empty());
          clades.emplace_back(std::move(clade_leafs));
        }
        clades |= ranges::actions::sort;
        return clades;
      }()},
      hash_{ComputeHash(clades_)} {}

template <typename Node>
LeafSet::LeafSet(Node node, const ContiguousMap<NodeId, NodeLabel>& labels,
                 ContiguousMap<NodeId, LeafSet>& computed_leafsets)
    : clades_{[&] {
        std::vector<std::vector<UniqueData>> clades;
        clades.reserve(node.GetCladesCount());
        for (auto clade : node.GetClades()) {
          std::vector<UniqueData> clade_leafs;
          clade_leafs.reserve(clade.size());
          for (Node child : clade | Transform::GetChild()) {
            if (child.IsLeaf()) {
              Assert(child.Const().HaveSampleId());
              UniqueData id = labels.at(child.GetId()).GetSampleId();
              clade_leafs.push_back(id);
            } else {
              for (auto& child_leafs : computed_leafsets.at(child.GetId()).clades_) {
                clade_leafs.insert(clade_leafs.end(), child_leafs.begin(),
                                   child_leafs.end());
              }
            }
          }
          clade_leafs |= ranges::actions::sort | ranges::actions::unique;
          Assert(not clade_leafs.empty());
          clades.emplace_back(std::move(clade_leafs));
        }
        clades |= ranges::actions::sort;
        return clades;
      }()},
      hash_{ComputeHash(clades_)} {}

LeafSet::LeafSet(std::vector<std::vector<UniqueData>>&& clades)
    : clades_{std::forward<std::vector<std::vector<UniqueData>>>(clades)},
      hash_{ComputeHash(clades_)} {}

bool LeafSet::operator==(const LeafSet& rhs) const noexcept {
  return clades_ == rhs.clades_;
}

size_t LeafSet::Hash() const noexcept { return hash_; }

auto LeafSet::begin() const -> decltype(clades_.begin()) { return clades_.begin(); }

auto LeafSet::end() const -> decltype(clades_.end()) { return clades_.end(); }

bool LeafSet::empty() const { return clades_.empty(); }

size_t LeafSet::size() const { return clades_.size(); }

std::vector<LeafSet::UniqueData> LeafSet::ToParentClade() const {
  std::vector<UniqueData> result = ranges::to_vector(clades_ | ranges::views::join);
  result |= ranges::actions::sort | ranges::actions::unique;
  return result;
}

size_t LeafSet::ParentCladeSize() const {
  size_t result = 0;
  for (const auto& clade : clades_) {
    result += clade.size();
  }
  return result;
}

const std::vector<std::vector<LeafSet::UniqueData>>& LeafSet::GetClades() const {
  return clades_;
}

std::string LeafSet::ToString() const {
  std::string result = "{";
  for (const auto& clade : GetClades()) {
    for (const auto* cg : clade) {
      Assert(cg != nullptr);
      result += cg->ToString();
      result += ", ";
    }
    result += "; ";
  }
  return result + "}";
}

template <typename DAGType>
ContiguousMap<NodeId, LeafSet> LeafSet::ComputeLeafSets(
    DAGType dag, const ContiguousMap<NodeId, NodeLabel>& labels) {
  ContiguousMap<NodeId, LeafSet> result;
  auto ComputeLS = [&](auto& self, auto for_node) {
    const LeafSet& leaf_set = result[for_node.GetId()];
    if (not leaf_set.empty()) {
      return;
    }
    for (auto child : for_node.GetChildren() | Transform::GetChild()) {
      self(self, child);
    }
    result.insert({for_node, LeafSet{for_node, labels, result}});
  };
  for (auto node : dag.GetNodes()) {
    ComputeLS(ComputeLS, node);
  }
  return result;
}

size_t LeafSet::ComputeHash(
    const std::vector<std::vector<UniqueData>>& clades) noexcept {
  size_t hash = 0;
  for (const auto& clade : clades) {
    for (const auto* leaf : clade) {
      hash = HashCombine(hash, leaf->Hash());
    }
  }
  return hash;
}

std::size_t std::hash<LeafSet>::operator()(const LeafSet& ls) const noexcept {
  return ls.Hash();
}

bool std::equal_to<LeafSet>::operator()(const LeafSet& lhs,
                                        const LeafSet& rhs) const noexcept {
  return lhs == rhs;
}
