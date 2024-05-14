#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/mat_view.hpp"
#include "sample_dag.hpp"

using Storage = ExtendStorageType<void, MATViewStorage, Extend::Nodes<CompactGenome>,
                                  Extend::DAG<ReferenceSequence>, Extend::Empty<>,
                                  CondensedViewBase>;

template <typename DAGView>
void test_mat_view_impl(DAGView dag) {
  // Create MAT Conversion
  TestAssert(dag.IsTree());
  dag.GetRoot().Validate(true);
  auto mat_conv = AddMATConversion(dag);
  MAT::Tree mat;
  mat_conv.View().BuildMAT(mat);
  mat_conv.View().GetRoot().Validate(true);

  // check BuildMAT
  check_MAT_MADAG_Eq(mat, dag);

  std::cout << "\nDAG view\n";
  MADAGToDOT(dag, std::cout);

  // Create MAT View
  MATStorageImpl matview_storage_impl;
  auto matview_storage = matview_storage_impl.GetCondensed();
  matview_storage.View().SetMAT(std::addressof(mat));
  auto storage = Storage::Consume(std::move(matview_storage));
  auto mv = storage.View();
  static_assert(mv.IsCondensed());
  mv.SetReferenceSequence(dag.GetReferenceSequence());
  mv.BuildRootAndLeafs();
  mv.RecomputeCompactGenomes();
  mv.GetRoot().Validate(true, false);

  auto umv = mv.GetUncondensed();
  static_assert(not umv.IsCondensed());

  // check BuildFromMAT
  auto dag_from_mat = AddMATConversion(MergeDAGStorage<>::EmptyDefault());
  dag_from_mat.View().BuildFromMAT(mat, dag.GetReferenceSequence());
  dag_from_mat.View().GetRoot().Validate(true);
  check_MAT_MADAG_Eq(mat, dag_from_mat.View());

  auto merge_mv =
      ExtendDAGStorage<void, decltype(mv), Extend::Nodes<SampleId>, Extend::Empty<>,
                       Extend::Empty<>, DefaultViewBase>::FromView(mv);
  merge_mv.View().SampleIdsFromCG();
  Merge merge(mv.GetReferenceSequence());

  std::cout << "\n\nMAT view\n";
  MADAGToDOT(merge_mv.View(), std::cout);

  merge.AddDAG(merge_mv.View());
  merge.GetResult().GetRoot().Validate(true, true);

  std::cout << "\n\nMerge result\n";
  MADAGToDOT(merge.GetResult(), std::cout);
}

void test_mat_view(std::string_view input_dag_path, std::string_view refseq_path,
                   std::string vcf_path) {
  MADAGStorage dag_storage =
      refseq_path.empty()
          ? LoadDAGFromProtobuf(input_dag_path)
          : LoadTreeFromProtobuf(input_dag_path, LoadReferenceSequence(refseq_path));

  auto dag = dag_storage.View();
  dag.RecomputeCompactGenomes(true);
  LoadVCFData(dag_storage, vcf_path);
  // if the DAG is from a DAG protobuf file, then it needs to be equipped with SampleIds
  if (vcf_path.empty()) {
    dag.SampleIdsFromCG();
  }
  test_mat_view_impl(dag);
}

void test_sample_dag() {
  MADAGStorage dag_storage = make_sample_dag();
  test_mat_view_impl(dag_storage.View());
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[]() { test_sample_dag(); }, "MATView: sample dag"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[]() {
                test_mat_view("data/startmat/startmat_no_ancestral.pb.gz",
                              "data/startmat/refseq.txt.gz", "");
              },
              "MATView: startmat"});
