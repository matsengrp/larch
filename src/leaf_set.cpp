#include "leaf_set.hpp"

#include <range/v3/action/sort.hpp>
#include <range/v3/action/unique.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/range/conversion.hpp>

#include "dag.hpp"
#include "compact_genome.hpp"
#include "node_label.hpp"

const LeafSet* LeafSet::Empty() {
  static LeafSet empty = {};
  return &empty;
}

LeafSet::LeafSet(Node node, const std::vector<NodeLabel>& labels,
                 std::vector<LeafSet>& computed_leafsets)
    : clades_{[&] {
        std::vector<std::vector<const CompactGenome*>> clades;
        clades.reserve(node.GetCladesCount());
        for (auto clade : node.GetClades()) {
          std::vector<const CompactGenome*> clade_leafs;
          clade_leafs.reserve(clade.size());
          for (Node child : clade | Transform::GetChild()) {
            const LeafSet& child_leaf_set = computed_leafsets.at(child.GetId().value);
            if (child.IsLeaf()) {
              clade_leafs.push_back(labels.at(child.GetId().value).GetCompactGenome());
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

auto LeafSet::begin() const -> decltype(clades_.begin()) { return clades_.begin(); }

auto LeafSet::end() const -> decltype(clades_.end()) { return clades_.end(); }

bool LeafSet::empty() const { return clades_.empty(); }

size_t LeafSet::size() const { return clades_.size(); }

std::vector<const CompactGenome*> LeafSet::ToParentClade() const {
  std::vector<const CompactGenome*> result =
      ranges::to_vector(clades_ | ranges::views::join);
  result |= ranges::actions::sort | ranges::actions::unique;
  return result;
}

const std::vector<std::vector<const CompactGenome*>>& LeafSet::GetClades() const {
  return clades_;
}

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
