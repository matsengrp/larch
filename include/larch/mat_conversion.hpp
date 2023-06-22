#pragma once

#if 0
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wdeprecated-copy-with-user-provided-copy"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wconversion"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#include "src/matOptimize/check_samples.hpp"
#include "src/matOptimize/tree_rearrangement_internal.hpp"
#include "src/matOptimize/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators.hpp"
#pragma GCC diagnostic pop
#endif

#include "larch/optimize/optimize.hpp"

#include "larch/madag/mutation_annotated_dag.hpp"

// namespace Mutation_Annotated_Tree {
// class Tree;
// class Node;
// }  // namespace Mutation_Annotated_Tree
// namespace MAT = Mutation_Annotated_Tree;

using MATNodePtr = MAT::Node*;

struct MATConversion {
  MATNodePtr mat_node_ptr_ = nullptr;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<MATConversion, CRTP, Tag> {
  bool HaveMATNode() const;
  MATNodePtr GetMATNode() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<MATConversion, CRTP, Tag> {
  MATNodePtr GetMutableMATNode() const;

 private:
  template <typename, typename>
  friend struct ExtraFeatureMutableView;
  void SetMATNode(MATNodePtr ptr) const;
};

template <>
struct ExtraFeatureStorage<MATConversion> {
  ExtraFeatureStorage() = default;
  MOVE_ONLY(ExtraFeatureStorage);
  MAT::Tree* mat_tree_ = nullptr;
  ContiguousMap<MATNodePtr, NodeId> reverse_map_;
};

template <typename CRTP>
struct ExtraFeatureConstView<MATConversion, CRTP> {
  const MAT::Tree& GetMAT() const;
  auto GetNodeFromMAT(MATNodePtr ptr) const;
};

template <typename CRTP>
struct ExtraFeatureMutableView<MATConversion, CRTP> {
  MAT::Tree& GetMutableMAT() const;
  auto GetMutableNodeFromMAT(MATNodePtr ptr) const;
  void BuildMAT(MAT::Tree& tree) const;
  void BuildFromMAT(MAT::Tree& mat, std::string_view reference_sequence) const;

 private:
  template <typename Node>
  static void BuildHelper(Node dag_node, MATNodePtr mat_par_node, MAT::Tree& new_tree);
  template <typename Node, typename DAG>
  static void BuildHelper(MATNodePtr par_node, Node node, DAG dag);
};

template <typename DAG>
auto AddMATConversion(DAG&& dag) {
  return ExtendStorage(std::forward<DAG>(dag), Extend::Nodes<MATConversion>{});
}

#include "larch/impl/mat_conversion_impl.hpp"
