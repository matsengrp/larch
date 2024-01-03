#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/mat_view.hpp"

using Storage = MergeDAGStorage<>;

using MVStorage =
    DAGStorage<void, MATNodesContainer, MATEdgesContainer, ExtraStorage<Connections>>;

void test_mat_view(std::string_view input_dag_path, std::string_view refseq_path,
                   std::string vcf_path) {
  MADAGStorage dag_storage =
      refseq_path.empty()
          ? LoadDAGFromProtobuf(input_dag_path)
          : LoadTreeFromProtobuf(input_dag_path, LoadReferenceSequence(refseq_path));

  dag_storage.View().RecomputeCompactGenomes(true);
  dag_storage.View().GetRoot().Validate(true);
  LoadVCFData(dag_storage, vcf_path);
  auto dag = dag_storage.View();
  // if the DAG is from a DAG protobuf file, then it needs to be equipped with SampleIds
  if (vcf_path.empty()) {
    dag.SampleIdsFromCG();
  }

  TestAssert(dag.IsTree());
  auto mat_conv = AddMATConversion(dag);
  MAT::Tree mat;
  mat_conv.View().BuildMAT(mat);
  mat_conv.View().GetRoot().Validate(true);

  // check BuildMAT
  check_MAT_MADAG_Eq(mat, dag);

  MVStorage mv_storage;
  auto mv = mv_storage.View();
  mv.SetMAT(std::addressof(mat));
  mv.BuildRootAndLeafs();

  // check BuildFromMAT
  auto dag_from_mat = AddMATConversion(Storage::EmptyDefault());
  dag_from_mat.View().BuildFromMAT(mat, dag.GetReferenceSequence());
  dag_from_mat.View().GetRoot().Validate(true);
  check_MAT_MADAG_Eq(mat, dag_from_mat.View());
}

[[maybe_unused]] static const auto test_added1 =
    add_test({[]() {
                test_mat_view("data/startmat/startmat_no_ancestral.pb.gz",
                              "data/startmat/refseq.txt.gz", "");
              },
              "MATView"});
