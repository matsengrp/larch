const CompactGenome* CompactGenome::Empty() {
  static CompactGenome empty = {};
  return &empty;
}

CompactGenome::CompactGenome(const Mutations& mutations, const CompactGenome& parent,
                             std::string_view reference_sequence)
    : mutations_{[&] {
        std::map<MutationPosition, char> result{parent.mutations_};
        for (auto [pos, base] : mutations) {
          if (base != reference_sequence.at(pos.value - 1)) {
            result[pos] = base;
          } else {
            result.erase(pos);
          }
        }
        return result;
      }()},
      hash_{[this] {
        size_t result = 0;
        for (auto [pos, base] : mutations_) {
          result = HashCombine(result, pos.value);
          result = HashCombine(result, base);
        }
        return result;
      }()} {}

bool CompactGenome::operator==(const CompactGenome& rhs) const {
  if (hash_ != rhs.hash_) return false;
  return mutations_ == rhs.mutations_;
}

bool CompactGenome::operator<(const CompactGenome& rhs) const {
  return mutations_ < rhs.mutations_;
}

size_t CompactGenome::Hash() const { return hash_; }

const LeafSet* LeafSet::Empty() {
  static LeafSet empty = {};
  return &empty;
}

LeafSet::LeafSet(Node node, const std::vector<NodeLabel>& labels)
    : clades_{[&] {
        std::set<std::set<const CompactGenome*>> clades;
        for (auto clade : node.GetClades()) {
          std::set<const CompactGenome*> clade_leafs;
          for (Node child : clade | ranges::views::transform(Transform::GetChild)) {
            const LeafSet& child_leaf_set = *labels.at(child.GetId().value).leaf_set;
            if (child.IsLeaf()) {
              clade_leafs.insert(labels.at(child.GetId().value).compact_genome);
            } else {
              for (auto& child_leafs :
                   (*labels.at(child.GetId().value).leaf_set).clades_) {
                clade_leafs.insert(child_leafs.begin(), child_leafs.end());
              }
            }
          }
          clades.emplace(std::move(clade_leafs));
        }
        return clades;
      }()},
      hash_{[this] {
        size_t hash = 0;
        return hash;
      }()} {}

bool LeafSet::operator==(const LeafSet& rhs) const {
  if (hash_ != rhs.hash_) return false;
  return clades_ == rhs.clades_;
}

bool LeafSet::operator<(const LeafSet& rhs) const { return clades_ < rhs.clades_; }

size_t LeafSet::Hash() const { return hash_; }

NodeLabel::NodeLabel()
    : compact_genome{CompactGenome::Empty()}, leaf_set{LeafSet::Empty()} {}

bool NodeLabel::operator==(const NodeLabel& rhs) const {
  return compact_genome == rhs.compact_genome && leaf_set == rhs.leaf_set;
}

size_t NodeLabel::Hash() const {
  return HashCombine(reinterpret_cast<std::uintptr_t>(compact_genome),
                     reinterpret_cast<std::uintptr_t>(leaf_set));
}

bool EdgeLabel::operator==(const EdgeLabel& rhs) const {
  return parent_compact_genome == rhs.parent_compact_genome &&
         parent_leaf_set == rhs.parent_leaf_set &&
         child_compact_genome == rhs.child_compact_genome &&
         child_leaf_set == rhs.child_leaf_set;
}

size_t EdgeLabel::Hash() const {
  size_t hash = 0;
  hash = HashCombine(hash, reinterpret_cast<std::uintptr_t>(parent_compact_genome));
  hash = HashCombine(hash, reinterpret_cast<std::uintptr_t>(parent_leaf_set));
  hash = HashCombine(hash, reinterpret_cast<std::uintptr_t>(child_compact_genome));
  hash = HashCombine(hash, reinterpret_cast<std::uintptr_t>(child_leaf_set));
  return hash;
}
