#pragma once

#include <vector>
#include <cstdint>
#include <cstddef>
#include <string>
#include <unordered_set>
#include <chrono>
#include <mutex>

#include <tbb/concurrent_unordered_map.h>

#include "larch/common.hpp"

static uint8_t one_hot_to_two_bit(uint8_t arg) {
  return static_cast<uint8_t>(31 - __builtin_clz(static_cast<unsigned int>(arg)));
}

class nuc_one_hot {
 public:
  nuc_one_hot() : nuc{0xff} {};
  nuc_one_hot(uint8_t n) : nuc{n} {}
  operator uint8_t() const { return nuc; }

 private:
  uint8_t nuc;
};

namespace Mutation_Annotated_Tree {
class Mutation {
 public:
  Mutation(const std::string& chromosome, int pos, nuc_one_hot mut, nuc_one_hot par,
           nuc_one_hot tie, nuc_one_hot ref = 0)
      : position{pos},
        par_mut_nuc{static_cast<uint8_t>((par << 4) | (mut))},
        boundary1_all_major_allele{tie} {
    auto ins_result = chromosome_map.emplace(chromosome, chromosome_map.size());
    if (ins_result.second) {
      std::lock_guard<std::mutex> lk(ref_lock);
      chromosomes.push_back(chromosome);
    }
    chrom_idx = ins_result.first->second;
    if (ref) {
      std::lock_guard<std::mutex> lk(ref_lock);
      refs.resize(std::max(static_cast<int>(refs.size()), position + 1), 0);
      refs[position] = ref;
    }
  }
  inline int get_position() const { return position; }
  inline nuc_one_hot get_par_one_hot() const { return par_mut_nuc >> 4; }
  inline nuc_one_hot get_mut_one_hot() const { return par_mut_nuc & 0xf; }
  nuc_one_hot get_all_major_allele() const { return boundary1_all_major_allele & 0xf; }
  nuc_one_hot get_boundary1_one_hot() const { return boundary1_all_major_allele >> 4; }

  static inline std::vector<nuc_one_hot> refs{};
  static inline std::mutex ref_lock{};
  static inline tbb::concurrent_unordered_map<std::string, uint8_t> chromosome_map{};
  static inline std::vector<std::string> chromosomes{};

 private:
  int position;
  uint8_t chrom_idx;
  uint8_t par_mut_nuc;
  uint8_t boundary1_all_major_allele;
};

class Mutations_Collection {
 public:
  std::vector<Mutation> mutations;
  typedef std::vector<Mutation>::iterator iterator;
  typedef std::vector<Mutation>::const_iterator const_iterator;
  iterator begin() { return mutations.begin(); }
  iterator end() { return mutations.end(); }
  const_iterator begin() const { return mutations.begin(); }
  const_iterator end() const { return mutations.end(); }
  void reserve(size_t n) { mutations.reserve(n); }
  void push_back(const Mutation& m) {
    if (m.get_position() >= static_cast<int>((Mutation::refs.size() + 1)) &&
        m.get_position() != INT_MAX) {
      Assert(false and "strange size");
    }
    if (!mutations.empty()) {
      if (m.get_position() <= mutations.back().get_position()) {
        Assert(false and "Adding out of order");
      }
    }
    mutations.push_back(m);
  }
  iterator find_next(int pos) {
    auto iter = mutations.begin();
    for (; iter < mutations.end(); iter++) {
      if (iter->get_position() >= pos) {
        break;
      }
    }
    assert(iter == mutations.begin() || (iter - 1)->get_position() < pos);
    assert(iter == mutations.end() || (iter)->get_position() >= pos);
    return iter;
  }
  iterator find(int position) {
    auto iter = find_next(position);
    if (iter != mutations.end() && iter->get_position() > position) {
      return mutations.end();
    }
    return iter;
  }
};

class Node {
 public:
  Node(size_t id) : node_id{id}, parent{nullptr} {};
  bool is_leaf() const { return children.empty(); }
  size_t node_id;
  Node* parent;
  size_t level;
  size_t dfs_index;
  size_t dfs_end_index;
  Mutations_Collection mutations;
  std::vector<Node*> children;
  std::vector<std::string> clade_annotations;
};

static size_t level_helper(const Node* node) {
  size_t level = 0;
  for (auto child : node->children) {
    level = std::max(level, level_helper(child));
  }
  return level + 1;
}

static void depth_first_expansion_helper(
    Mutation_Annotated_Tree::Node* node,
    std::vector<Mutation_Annotated_Tree::Node*>& vec, size_t& index, size_t level) {
  vec.push_back(node);
  node->level = level;
  node->dfs_index = index;
  index++;
  for (auto c : node->children) {
    depth_first_expansion_helper(c, vec, index, level + 1);
  }
  node->dfs_end_index = index - 1;
}

class Tree {
 public:
  Node* root;
  std::vector<Node*> all_nodes;
  void register_node_serial(Node* node) {
    all_nodes.resize(std::max(all_nodes.size(), node->node_id + 1), nullptr);
    all_nodes[node->node_id] = node;
  }
  size_t get_num_annotations() const {
    size_t ret = 0;
    if (root != NULL) {
      ret = root->clade_annotations.size();
    }
    return ret;
  }
  size_t get_max_level() const {
    size_t max_level = level_helper(root);
    return max_level;
  }
  std::vector<Node*> depth_first_expansion(Node* node = nullptr) const {
    std::vector<Node*> traversal;
    if (node == NULL) {
      node = root;
    }
    size_t index = 0;
    if (node == NULL) {
      return traversal;
    }
    depth_first_expansion_helper(node, traversal, index, 0);
    return traversal;
  }
};

inline void save_mutation_annotated_tree(Mutation_Annotated_Tree::Tree tree,
                                         std::string filename) {
  // TODO
  std::ignore = tree;
  std::ignore = filename;
}

}  // namespace Mutation_Annotated_Tree
namespace MAT = Mutation_Annotated_Tree;

struct Profitable_Moves {
  int score_change;
  MAT::Node* src;
  MAT::Node* dst;
  MAT::Node* LCA;
};

class Mutation_Count_Change {
 public:
  Mutation_Count_Change()
      : position{INT_MAX}, decremented_allele{0}, incremented_allele{0} {}
  int get_position() const { return position; }
  nuc_one_hot get_decremented() const { return decremented_allele; }
  nuc_one_hot get_incremented() const { return incremented_allele; }

 private:
  int position;
  nuc_one_hot decremented_allele;
  nuc_one_hot incremented_allele;
};

typedef std::vector<Mutation_Count_Change> Mutation_Count_Change_Collection;

struct Node_With_Major_Allele_Set_Change {
  MAT::Node* node;
  Mutation_Count_Change_Collection major_allele_set_change;
};

struct Move_Found_Callback {
  virtual bool operator()(
      Profitable_Moves& move, int best_score_change,
      [[maybe_unused]] std::vector<Node_With_Major_Allele_Set_Change>&
          node_with_major_allele_set_change) {
    return move.score_change <= best_score_change;
  }
};

struct Mutation_Pos_Only_Comparator {
  bool operator()(const Mutation_Annotated_Tree::Mutation& first,
                  const Mutation_Annotated_Tree::Mutation& second) const {
    return (first.get_position() == second.get_position());
  }
};
struct Mutation_Pos_Only_Hash {
  size_t operator()(const Mutation_Annotated_Tree::Mutation& in) const {
    return in.get_position();
  }
};

typedef std::unordered_set<Mutation_Annotated_Tree::Mutation, Mutation_Pos_Only_Hash,
                           Mutation_Pos_Only_Comparator>
    Mutation_Set;
typedef tbb::concurrent_unordered_map<size_t, Mutation_Set> Original_State_t;

inline void check_samples(const Mutation_Annotated_Tree::Node* root,
                          Original_State_t& samples, const MAT::Tree* tree,
                          bool ignore_missed_samples = false) {
  // TODO
  std::ignore = root;
  std::ignore = samples;
  std::ignore = tree;
  std::ignore = ignore_missed_samples;
}

inline void reassign_states(MAT::Tree& t, Original_State_t& origin_states) {
  // TODO
  std::ignore = t;
  std::ignore = origin_states;
}

inline size_t optimize_inner_loop(
    std::vector<MAT::Node*>& nodes_to_search, MAT::Tree& t, int radius,
    Move_Found_Callback& callback, bool allow_drift = false, bool search_all_dir = true,
    int minutes_between_save = 0, bool no_write_intermediate = true,
    std::chrono::steady_clock::time_point search_end_time =
        std::chrono::steady_clock::time_point::max(),
    std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now(),
    bool log_moves = false, int iteration = 1, std::string intermediate_template = "",
    std::string intermediate_pb_base_name = "", std::string intermediate_nwk_out = "") {
  // TODO
  std::ignore = nodes_to_search;
  std::ignore = t;
  std::ignore = radius;
  std::ignore = callback;
  std::ignore = allow_drift;
  std::ignore = search_all_dir;
  std::ignore = minutes_between_save;
  std::ignore = no_write_intermediate;
  std::ignore = search_end_time;
  std::ignore = start_time;
  std::ignore = log_moves;
  std::ignore = iteration;
  std::ignore = intermediate_template;
  std::ignore = intermediate_pb_base_name;
  std::ignore = intermediate_nwk_out;
  return 0;
}
