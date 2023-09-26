#include "test_common.hpp"
#include "sample_dag.hpp"
#include "larch/dag_loader.hpp"

using Storage = MergeDAGStorage;

template <typename DAG>
[[maybe_unused]] bool check_leaf_sample_ids(DAG dag, MAT::Tree& tree) {
  std::set<std::string> dag_leaf_sample_ids;
  std::set<std::string> mat_leaf_node_names;

  for (auto dag_leaf: dag.GetLeafs()) {
    Assert(dag_leaf.HaveSampleId());
    dag_leaf_sample_ids.insert(dag_leaf.GetSampleId().value_or("NA"));
  }
  for (auto mat_leaf: tree.get_leaves()) {
    mat_leaf_node_names.insert(tree.get_node_name(mat_leaf->node_id).c_str());
  }
  return dag_leaf_sample_ids == mat_leaf_node_names;
}

[[maybe_unused]] void test_sample_id_conversion(std::string_view input_dag_path, std::string_view refseq_path, std::string vcf_path) {
  MADAGStorage dag_storage =
      refseq_path.empty()
          ? LoadDAGFromProtobuf(input_dag_path)
          : LoadTreeFromProtobuf(input_dag_path, LoadReferenceSequence(refseq_path));

  dag_storage.View().RecomputeCompactGenomes(true);
  LoadVCFData(dag_storage, vcf_path);
  auto dag = dag_storage.View();
  // if the DAG is from a DAG protobuf file, then it needs to be equipped with SampleIds
  if (vcf_path.empty()) {
    dag.SampleIdsFromCG();
  }

  Assert(dag.IsTree());
  auto mat_conv = AddMATConversion(dag);
  MAT::Tree mat;
  mat_conv.View().BuildMAT(mat);

  // check BuildMAT
  check_MAT_MADAG_Eq(mat, dag);
  Assert(check_leaf_sample_ids(dag, mat));

  // check BuildFromMAT
  auto dag_from_mat = AddMATConversion(Storage{{}});
  dag_from_mat.View().BuildFromMAT(mat, dag.GetReferenceSequence());
  check_MAT_MADAG_Eq(mat, dag_from_mat.View());
  Assert(check_leaf_sample_ids(dag_from_mat.View(), mat));
}


[[maybe_unused]] static const auto test_added1 =
    add_test({[]() { test_sample_id_conversion("data/startmat/startmat_no_ancestral.pb.gz",
                                               "data/startmat/refseq.txt.gz",
                                               ""); },
                     "Check SampleIds on startmat"});

[[maybe_unused]] static const auto test_added2 =
    add_test({[]() { test_sample_id_conversion("data/test_ambiguous_vcf/amb_mat.pb",
                                               "data/test_ambiguous_vcf/sample_reference_sequence.fasta",
                                               "data/test_ambiguous_vcf/amb.vcf"); },
                     "Check SampleIds on MakeSampleDAG with ambiguous VCF input"});
