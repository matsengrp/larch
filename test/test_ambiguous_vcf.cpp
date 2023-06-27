
#include "test_common.hpp"

#include "larch/dag_loader.hpp"
#include "larch/spr/spr_view.hpp"

#include "larch/usher_glue.hpp"

[[maybe_unused]] static auto LoadMATFromProtobuf(std::string_view path,
                                                 std::string_view ref_seq) {
  auto dag_storage = LoadTreeFromProtobuf(path, ref_seq);
  auto dag = dag_storage.View();
  Assert(dag.IsTree());
  auto mat_conv = AddMATConversion(dag);
  mat_conv.View().BuildMAT(mat);
  return mat;
}

std::ostream &operator<<(std::ostream &os, const MAT::Mutation mut) {
  os << mut.get_string();
  return os;
}

std::ostream &operator<<(std::ostream &os, const MAT::Tree &mat) {
  os << "newick: " << mat.get_newick_string(true, false, false, true) << std::endl;

  std::vector<MAT::Node *> nodes;
  nodes.push_back(mat.root);
  while (!nodes.empty()) {
    auto node = nodes.back();
    nodes.pop_back();
    for (auto child : node->children) {
      nodes.push_back(child);
    }
    os << "node: " << node->node_id << " " << mat.get_node_name(node->node_id)
       << " -> children [ ";
    for (auto child : node->children) {
      os << child->node_id << " ";
    }
    os << "]" << std::endl;
    os << "muts: [ ";
    for (const auto mut : node->mutations) {
      os << mut.get_string() << " ";
    }
    os << "]" << std::endl;
  }
  return os;
}

std::ostream &operator<<(std::ostream &os, const Original_State_t &og_state) {}

[[maybe_unused]] static auto test_ambiguous_vcf() {
  std::string mat_pb_filename = "data/test_ambiguous_vcf/unambiguous_mat.pb";
  std::string vcf_filename = "data/test_ambiguous_vcf/ambiguous_sample_dag.vcf";
  MAT::Tree mat_from_pb;

  // - construct a DAG from protobuf (tree shaped, since it's from a MAT protobuf)
  auto dag_storage = LoadTreeFromProtobuf(mat_pb_filename, "GAA");
  auto dag = dag_storage.View();

  std::string file_out = "_ignore/unambiguous_mat.dot";
  std::ofstream out;
  out.open(file_out);
  MADAGToDOT(dag, out);
  out.close();

  // - create a MAT conversion object for the DAG
  Assert(dag.IsTree());
  auto mat_conv = AddMATConversion(dag);
  mat_conv.View().BuildMAT(mat_from_pb);
  std::cout << "MAT_FROM_PB: " << mat_from_pb << std::endl;

  // - call VCF_input() on the MAT
  std::cout << "VCF [before]" << std::endl;
  VCF_input(vcf_filename.c_str(), mat_from_pb);
  std::cout << "MAT_FROM_VCF: " << mat_from_pb << std::endl;
  std::cout << "VCF [after]" << std::endl;

  // - create an Original_State_t object to hold the ambiguous leaf data
  Original_State_t og_state;

  // - call check_samples() on the MAT and pass in the Original_State_t object to
  // populate it.
  check_samples(mat_from_pb.root, og_state, &mat_from_pb);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[]() { test_ambiguous_vcf(); }, "Loading VCFs with Ambiguities Test"});
