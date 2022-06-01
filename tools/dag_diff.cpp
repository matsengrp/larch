#include <cstdlib>
#include <iostream>
#include <set>
#include <unordered_set>

#include <range/v3/view/enumerate.hpp>

#include "arguments.hpp"
#include "merge.hpp"
#include "dag_loader.hpp"

[[noreturn]] static void Usage() {
  std::cout << "Usage:\n";
  std::cout << "dag_diff -p,--proto file -j,--json file"
               "[-o,--output filename]\n";
  std::cout << "  -p,--proto     Input protobuf DAG filename\n";
  std::cout << "  -j,--json      Input json DAG filename\n";

  std::exit(EXIT_SUCCESS);
}

[[noreturn]] static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
}

class DeepNodeLabel {
 public:
  explicit DeepNodeLabel(const NodeLabel& label)
      : compact_genome{label.compact_genome->Copy()}, leaf_set{[&] {
          std::set<std::set<CompactGenome>> result;
          for (auto& clade : *label.leaf_set) {
            std::set<CompactGenome> leafs;
            for (auto* leaf : clade) {
              leafs.insert(leaf->Copy());
            }
            result.insert(std::move(leafs));
          }
          return result;
        }()} {}

  CompactGenome compact_genome;
  std::set<std::set<CompactGenome>> leaf_set;
};

struct DeepNodeLabelHash {
  size_t operator()(const DeepNodeLabel& label) const {
    size_t hash = label.compact_genome.Hash();
    for (auto& clade : label.leaf_set) {
      for (auto& leaf : clade) {
        hash = HashCombine(hash, leaf.Hash());
      }
    }
    return hash;
  }
};

struct DeepNodeLabelEq {
  size_t operator()(const DeepNodeLabel& lhs, const DeepNodeLabel& rhs) const {
    if (not(lhs.compact_genome == rhs.compact_genome)) {
      return false;
    }
    if (not(lhs.leaf_set == rhs.leaf_set)) {
      return false;
    }
    return true;
  }
};

static void Print(const CompactGenome& cg) {
  for (auto [pos, base] : cg) {
    std::cout << pos.value << base << " ";
  }
}

static void Print(const DeepNodeLabel& label) {
  std::cout << "{";
  Print(label.compact_genome);
  std::cout << "} [";
  for (auto i : label.leaf_set | ranges::view::enumerate) {
    for (auto j : i.second | ranges::view::enumerate) {
      Print(j.second);
      if (j.first + 1 < i.second.size()) {
        std::cout << ", ";
      }
    }
    if (i.first + 1 < label.leaf_set.size()) {
      std::cout << " | ";
    }
  }
  std::cout << "]";
}

static int TakeDiff(std::string_view proto_filename, std::string_view json_filename) {
  std::string lhs_reference_sequence;
  std::vector<std::vector<Mutations>> lhs_mutations;
  std::vector<DAG> lhs_trees;
  lhs_mutations.push_back({});
  lhs_trees.push_back(LoadDAGFromProtobufGZ(proto_filename, lhs_reference_sequence,
                                            lhs_mutations.at(0)));
  Merge lhs_merge{lhs_reference_sequence};
  std::vector<std::reference_wrapper<const DAG>> lhs_tree_refs{lhs_trees.begin(),
                                                               lhs_trees.end()};
  lhs_merge.AddTrees(lhs_tree_refs, lhs_mutations, true);

  std::string rhs_reference_sequence;
  std::vector<std::vector<CompactGenome>> rhs_compact_genomes;
  std::vector<DAG> rhs_trees;
  rhs_trees.push_back(LoadDAGFromJsonGZ(json_filename, rhs_reference_sequence));
  rhs_compact_genomes.push_back(LoadCompactGenomesJsonGZ(json_filename));
  Merge rhs_merge{rhs_reference_sequence};
  std::vector<std::reference_wrapper<const DAG>> rhs_tree_refs{rhs_trees.begin(),
                                                               rhs_trees.end()};
  rhs_merge.AddDAGs(rhs_tree_refs, std::move(rhs_compact_genomes), true);

  // std::string rhs_reference_sequence;
  // std::vector<std::vector<Mutations>> rhs_mutations;
  // std::vector<DAG> rhs_trees;
  // rhs_mutations.push_back({});
  // rhs_trees.push_back(LoadDAGFromProtobufGZ(json_filename, rhs_reference_sequence,
  //                                           rhs_mutations.at(0)));
  // Merge rhs_merge{rhs_reference_sequence};
  // std::vector<std::reference_wrapper<const DAG>> rhs_tree_refs{rhs_trees.begin(),
  //                                                              rhs_trees.end()};
  // rhs_merge.AddTrees(rhs_tree_refs, rhs_mutations, true);

  size_t not_found_in_lhs = 0, not_found_in_rhs = 0;

  std::unordered_set<DeepNodeLabel, DeepNodeLabelHash, DeepNodeLabelEq> lhs_nodes;
  std::unordered_set<DeepNodeLabel, DeepNodeLabelHash, DeepNodeLabelEq> rhs_nodes;

  for (auto lhs : lhs_merge.GetResultNodes()) {
    lhs_nodes.insert(DeepNodeLabel{lhs.first});
  }

  for (auto rhs : rhs_merge.GetResultNodes()) {
    rhs_nodes.insert(DeepNodeLabel{rhs.first});
  }

  for (auto& lhs : lhs_nodes) {
    if (rhs_nodes.find(lhs) == rhs_nodes.end()) {
      std::cout << "Not found in rhs: ";
      Print(lhs);
      std::cout << "\n";
      ++not_found_in_rhs;
    }
  }

  for (auto& rhs : rhs_nodes) {
    if (lhs_nodes.find(rhs) == lhs_nodes.end()) {
      // std::cout << "Not found in lhs: ";
      // Print(rhs.first);
      // std::cout << "\n";
      ++not_found_in_lhs;
    }
  }

  std::cout << "Not found in rhs: " << not_found_in_rhs << "\n";
  std::cout << "Not found in lhs: " << not_found_in_lhs << "\n";

  return EXIT_SUCCESS;
}

int main(int argc, char** argv) {
  Arguments args = GetArguments(argc, argv);

  std::string_view proto_filename;
  std::string_view json_filename;

  for (auto [name, params] : args) {
    if (name == "-h" or name == "--help") {
      Usage();
    } else if (name == "-p" or name == "--proto") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      proto_filename = *params.begin();
    } else if (name == "-j" or name == "--json") {
      if (params.empty()) {
        std::cerr << "Filename not specified.\n";
        Fail();
      }
      json_filename = *params.begin();
    } else {
      std::cerr << "Unknown argument.\n";
      Fail();
    }
  }

  if (proto_filename.empty() or json_filename.empty()) {
    std::cerr << "Specify two input file names.\n";
    Fail();
  }

  return TakeDiff(proto_filename, json_filename);
}
