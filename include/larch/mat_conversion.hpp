#pragma once

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

#include "larch/madag/mutation_annotated_dag.hpp"

namespace Mutation_Annotated_Tree {
class Tree;
class Node;
}  // namespace Mutation_Annotated_Tree
namespace MAT = Mutation_Annotated_Tree;

struct MATConversion {
  size_t mat_node_id_ = NoId;
};

template <typename CRTP, typename Tag>
struct FeatureConstView<MATConversion, CRTP, Tag> {
  bool HaveMATNode() const;
  const MAT::Node& GetMATNode() const;
  size_t GetMATNodeId() const;
};

template <typename CRTP, typename Tag>
struct FeatureMutableView<MATConversion, CRTP, Tag> {
  MAT::Node& GetMutableMATNode() const;

 private:
  template <typename, typename>
  friend struct ExtraFeatureMutableView;
  void SetMATNodeId(size_t id) const;
};

template <>
struct ExtraFeatureStorage<MATConversion> {
  ExtraFeatureStorage() = default;
  MOVE_ONLY(ExtraFeatureStorage);
  MAT::Tree mat_tree_;
};

template <typename CRTP>
struct ExtraFeatureConstView<MATConversion, CRTP> {
  const MAT::Tree& GetMAT() const;
};

template <typename CRTP>
struct ExtraFeatureMutableView<MATConversion, CRTP> {
  MAT::Tree& GetMutableMAT() const;
  MAT::Tree& BuildMAT() const;
  void BuildFromMAT(MAT::Tree&& mat, std::string_view reference_sequence) const;

 private:
  template <typename Node>
  static void BuildHelper(Node dag_node, MAT::Node* mat_par_node, MAT::Tree& new_tree);
  template <typename Node, typename DAG>
  static void BuildHelper(MAT::Node* par_node, Node node, DAG dag);
};

template <typename DAG>
auto AddMATConversion(DAG&& dag) {
  return ExtendStorage(std::forward<DAG>(dag), Extend::Nodes<MATConversion>{});
}

#include "larch/impl/mat_conversion_impl.hpp"