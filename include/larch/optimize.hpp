#pragma once

#ifdef USE_USHER
#error Dont include this header when optimizing with Usher
#else

#include <vector>
#include <cstdint>
#include <cstddef>
#include <string>
#include <unordered_set>
#include <chrono>
#include <mutex>
#include <random>
#include <set>
#include <queue>

#include "larch/parallel/growable_hash_map.hpp"

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

// Convert nuc_id back to IUPAC base
inline char get_nuc(int8_t nuc_id) {
  char ret = 'N';
  // assert ((nuc_id >= 1) && (nuc_id <= 15));
  switch (nuc_id) {
    case 1:
      ret = 'A';
      break;
    case 2:
      ret = 'C';
      break;
    case 3:
      ret = 'M';
      break;
    case 4:
      ret = 'G';
      break;
    case 5:
      ret = 'R';
      break;
    case 6:
      ret = 'S';
      break;
    case 7:
      ret = 'V';
      break;
    case 8:
      ret = 'T';
      break;
    case 9:
      ret = 'W';
      break;
    case 10:
      ret = 'Y';
      break;
    case 11:
      ret = 'H';
      break;
    case 12:
      ret = 'K';
      break;
    case 13:
      ret = 'D';
      break;
    case 14:
      ret = 'B';
      break;
    default:
      ret = 'N';
      break;
  }
  return ret;
}

namespace Mutation_Annotated_Tree {
class Mutation {
 public:
  Mutation(const std::string& chromosome, int pos, nuc_one_hot mut, nuc_one_hot par,
           nuc_one_hot tie, nuc_one_hot ref = 0)
      : position{pos},
        par_mut_nuc{static_cast<uint8_t>((par << 4) | (mut))},
        boundary1_all_major_allele{tie} {
    auto ins_result = chromosome_map.insert({chromosome, chromosome_map.size()});
    if (ins_result.second) {
      std::lock_guard<std::mutex> lk(ref_lock);
      chromosomes.push_back(chromosome);
    }
    chrom_idx = ins_result.first;
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
  static inline GrowableHashMap<std::string, uint8_t> chromosome_map{32};
  static inline std::vector<std::string> chromosomes{};

  inline bool is_masked() const { return (position < 0); }
  inline std::string get_string() const {
    if (is_masked()) {
      return "MASKED";
    } else {
      return get_nuc(get_par_one_hot()) + std::to_string(position) +
             get_nuc(get_mut_one_hot());
    }
  }

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
  void clear() { mutations.clear(); }
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
  size_t bfs_index;  // index in bfs
  void delete_this() {
    for (Node* n : children) {
      n->delete_this();
    }
    delete this;
  }
  int branch_length = {};
  bool have_masked = {};
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
  // TODO ConcurrentUnorderedMap
  typedef std::unordered_map<size_t, std::vector<std::string> > condensed_node_t;
  Node* root = nullptr;
  std::vector<Node*> all_nodes;
  size_t node_idx = 0;
  size_t num_nodes = 0;
  condensed_node_t condensed_nodes;
  std::unordered_map<size_t, std::string> node_names;
  std::unordered_map<std::string, size_t> node_name_to_idx_map;
  void register_node_serial(Node* node) {
    all_nodes.resize(std::max(all_nodes.size(), node->node_id + 1), nullptr);
    all_nodes[node->node_id] = node;
  }
  void register_node_serial(Node* node, std::string& name) {
    register_node_serial(node);
    node_names.emplace(node->node_id, name);
    node_name_to_idx_map.emplace(name, node->node_id);
  }
  size_t node_name_to_node_idx(const std::string& in) const {
    auto iter = node_name_to_idx_map.find(in);
    if (iter == node_name_to_idx_map.end()) {
      return -1;
    }
    return iter->second;
  }
  size_t get_num_annotations() const {
    size_t ret = 0;
    if (root != NULL) {
      ret = root->clade_annotations.size();
    }
    return ret;
  }
  size_t get_max_level() const { return level_helper(root); }
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
  void delete_nodes() {
    if (!root) {
      return;
    }
    root->delete_this();
    root = nullptr;
  }
  void fix_node_idx() {
    auto dfs = depth_first_expansion();
    for (auto node : dfs) {
      node_idx = std::max(node_idx, node->node_id);
      ++num_nodes;
    }
    node_idx++;
  }

  size_t get_node_idx() const { return num_nodes; }
  size_t get_size_upper() const { return all_nodes.size(); }

  Node* get_node(const std::string& identifier) const {
    auto iter = node_name_to_idx_map.find(identifier);
    if (iter != node_name_to_idx_map.end()) {
      return all_nodes[iter->second];
    }
    return NULL;
  }
  Node* get_node(size_t idx) const {
    if (idx >= all_nodes.size()) {
      return nullptr;
    }

    return all_nodes[idx];
  }
  std::string get_node_name(size_t idx) const {
    auto node_name_iter = node_names.find(idx);
    if (node_name_iter == node_names.end()) {
      return "";
    }
    return node_name_iter->second;
  }
  static void get_leaves_helper(const Node* root, std::vector<Node*>& out) {
    for (auto child : root->children) {
      if (child->is_leaf()) {
        out.push_back(child);
      } else {
        get_leaves_helper(child, out);
      }
    }
  }
  std::vector<Node*> get_leaves(const Node* r = nullptr) const {
    std::vector<Node*> out;
    if (r == nullptr) {
      r = this->root;
    }
    get_leaves_helper(r, out);
    return out;
  }

  // TODO implement
  std::string get_newick_string(bool print_internal, bool print_branch_len,
                                bool retain_original_branch_len,
                                bool uncondense_leaves) const {
    std::ignore = print_internal;
    std::ignore = print_branch_len;
    std::ignore = retain_original_branch_len;
    std::ignore = uncondense_leaves;
    return "TODO";
  }

  std::vector<Mutation_Annotated_Tree::Node*> breadth_first_expansion(
      std::string nid = "") {
    std::vector<Node*> traversal;
    Node* node;
    if (nid == "") {
      if (root == NULL) {
        return traversal;
      }
      node = root;
    } else {
      node = get_node(nid);
    }
    size_t idx = 0;
    std::queue<Node*> remaining_nodes;
    remaining_nodes.push(node);
    while (remaining_nodes.size() > 0) {
      Node* curr_node = remaining_nodes.front();
      curr_node->bfs_index = idx++;
      traversal.push_back(curr_node);
      remaining_nodes.pop();
      for (auto c : curr_node->children) {
        remaining_nodes.push(c);
      }
    }

    return traversal;
  }
  void condense_leaves(std::vector<std::string> = std::vector<std::string>{}) {}
  void uncondense_leaves() {}
  size_t root_ident = {};
  size_t max_level = {};
  size_t curr_internal_node = {};
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
  static Move_Found_Callback& default_instance() {
      static Move_Found_Callback instance;
      return instance;
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
// TODO ConcurrentUnorderedMap
typedef std::unordered_map<size_t, Mutation_Set> Original_State_t;

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
    bool allow_drift = false, bool search_all_dir = true,
    int minutes_between_save = 0, bool no_write_intermediate = true,
    std::chrono::steady_clock::time_point search_end_time =
        std::chrono::steady_clock::time_point::max(),
    std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now(),
    bool log_moves = false, int iteration = 1, std::string intermediate_template = "",
    std::string intermediate_pb_base_name = "", std::string intermediate_nwk_out = "",
    // std::optional<uint32_t> user_seed = std::nullopt,
    Move_Found_Callback& callback = Move_Found_Callback::default_instance()) {
  Assert(not nodes_to_search.empty());
  std::random_device random_device;
  auto rand = random_device();//user_seed.value_or(random_device());
  std::mt19937 gen(rand);
  std::uniform_int_distribution<size_t> dist(0, nodes_to_search.size() - 1);
  std::mutex random_mtx;
  auto random_node = [&] {
    std::unique_lock lock{random_mtx};
    return dist(gen);
  };

  std::vector<size_t> idxs;
  idxs.resize(nodes_to_search.size());
  std::iota(idxs.begin(), idxs.end(), 0);

  ParallelForEach(idxs, [&](size_t) {
    retry:
      MAT::Node* src = nodes_to_search.at(random_node());
      if (src->parent == nullptr) {
        goto retry;
      }
      MAT::Node* dst = [&] {
        while (true) {
        outer:
          MAT::Node* result = nodes_to_search.at(random_node());
          if (result != src and result->parent != nullptr and result != src->parent) {
            MAT::Node* parent = result->parent;
            while (parent != nullptr) {
              if (parent == src) {
                goto outer;
              }
              parent = parent->parent;
            }
            return result;
          }
        }
      }();
      std::set<MAT::Node*> parents;
      MAT::Node* lca = nullptr;
      MAT::Node* src_parent = src;
      MAT::Node* dst_parent = dst;
      while (lca == nullptr and (src_parent != nullptr or dst_parent != nullptr)) {
        if (src_parent != nullptr and not parents.insert(src_parent).second) {
          lca = src_parent;
          break;
        } else {
          if (src_parent != nullptr) {
            src_parent = src_parent->parent;
          }
        }
        if (dst_parent != nullptr and not parents.insert(dst_parent).second) {
          lca = dst_parent;
          break;
        } else {
          if (dst_parent != nullptr) {
            dst_parent = dst_parent->parent;
          }
        }
      }
      Assert(lca != nullptr);
      Assert(src != dst);
      if (lca->parent == nullptr) {
        goto retry;
      }
      Profitable_Moves move{-1, src, dst, lca};
      int best_score_change = 0;
      std::vector<Node_With_Major_Allele_Set_Change> node_with_major_allele_set_change;

      callback(move, best_score_change, node_with_major_allele_set_change);
  });

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

#endif
