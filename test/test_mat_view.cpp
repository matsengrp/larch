#if USE_MAT_VIEW
#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/mat_view.hpp"
#include "sample_dag.hpp"

using Storage = CondensedMADAGStorage;

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

  // condense MAT prior to building a MAT View
  Original_State_t origin_states;
  check_samples(mat.root, origin_states, &mat);
  reassign_states(mat, origin_states);

  // condense MAT leaves
  std::vector<std::string> condense_arg{};
  mat.condense_leaves(condense_arg);
  mat.fix_node_idx();

  // std::cout << "\nDAG view\n";
  // MADAGToDOT(dag, std::cout);

  // Create MAT View
  CondensedMATViewStorage matview_storage;
  matview_storage.View().SetMAT(std::addressof(mat));
  auto storage = Storage::Consume(std::move(matview_storage));
  auto mv = storage.View();
  static_assert(mv.IsCondensed());
  mv.SetReferenceSequence(dag.GetReferenceSequence());
  mv.BuildRootAndLeafs();
  mv.RecomputeCompactGenomes<IdContinuity::Sparse>();
  mv.GetRoot().Validate(true, false);

  // std::cout << "\nCondensed MAT view\n";
  // MADAGToDOT(mv, std::cout);

  auto umv = mv.GetUncondensed();
  static_assert(not umv.IsCondensed());

  // std::cout << "\nUnndensed MAT view\n";
  // MADAGToDOT(umv, std::cout);

  // check BuildFromMAT
  auto dag_from_mat = AddMATConversion(MergeDAGStorage<>::EmptyDefault());
  dag_from_mat.View().BuildFromMAT(mat, dag.GetReferenceSequence());
  dag_from_mat.View().GetRoot().Validate(true);
  check_MAT_MADAG_Eq(mat, dag_from_mat.View());

  auto merge_mv = ExtendDAGStorage<void, decltype(mv), Extend::Nodes<SampleId>,
                                   Extend::Empty<>, Extend::Empty<>, DefaultViewBase,
                                   IdContinuity::Sparse>::FromView(mv);
  merge_mv.View().RecomputeCompactGenomes<IdContinuity::Sparse>(true);
  merge_mv.View().SampleIdsFromCG();

  Merge merge(mv.GetReferenceSequence());

  // std::cout << "\n\nMAT view\n";
  // MADAGToDOT(merge_mv.View(), std::cout);

  merge.AddDAG(merge_mv.View());
  merge.GetResult().GetRoot().Validate(true, true);

  // std::cout << "\n\nMerge result\n";
  // MADAGToDOT(merge.GetResult(), std::cout);
}

void test_condensed_mat_view() {
  // create a sample dag with ambiguous leaves such that some of the
  // leaves are condensed when condense_leaves() is called.
  auto dag_storage = make_condensing_ambiguous_sample_dag();
  auto dag = dag_storage.View();
  dag.SampleIdsFromCG();

  // create a MAT that is copied from the dag
  auto mat_conv = AddMATConversion(dag);
  MAT::Tree mat;
  mat_conv.View().BuildMAT(mat);
  mat_conv.View().GetRoot().Validate(true);
  check_MAT_MADAG_Eq(mat, dag);

  // std::cout << "\n\nDAG view\n";
  // MADAGToDOT(dag_storage.View(), std::cout);

  // optimal labelings for MAT (this disambiguates the leaf nodes so they can be
  // condensed) NOTE: ordinarily, reassign_states optimizes internal node labels and
  // collapses any empty edges, which means that the topology and internal node labels
  // can also change, but this sample dag is handpicked so that only the leaf cgs will
  // change with reassign_states. This means that we can therefore compare the
  // uncondensed to the condensed MATs
  Original_State_t origin_states;
  check_samples(mat.root, origin_states, &mat);
  reassign_states(mat, origin_states);

  // condense MAT leaves
  std::vector<std::string> condense_arg{};
  mat.condense_leaves(condense_arg);
  mat.fix_node_idx();

  // create a MATView from the condensed MAT
  CondensedMATViewStorage matview_storage;
  matview_storage.View().SetMAT(std::addressof(mat));
  auto storage = Storage::Consume(std::move(matview_storage));
  auto mv = storage.View();
  static_assert(mv.IsCondensed());
  mv.SetReferenceSequence(dag.GetReferenceSequence());
  mv.BuildRootAndLeafs();

  mv.GetRoot().Validate(true, false);
  mv.RecomputeCompactGenomes<IdContinuity::Sparse>(true);

  auto umv = mv.GetUncondensed();
  static_assert(not umv.IsCondensed());
  auto merge_umv = ExtendDAGStorage<void, decltype(umv), Extend::Nodes<SampleId>,
                                    Extend::Empty<>, Extend::Empty<>, DefaultViewBase,
                                    IdContinuity::Sparse>::FromView(umv);
  umv.RecomputeCompactGenomes<IdContinuity::Sparse>(true);

  merge_umv.View().SampleIdsFromCG();
  // std::cout << "\n\nUncondensed MAT view\n";
  // MADAGToDOT(merge_umv.View(), std::cout);

  // FIXME umv.GetRoot().Validate(true, false);

  // TODO: check that dag and mat_view(condensed=False) are equal
  auto merge_mv = ExtendDAGStorage<void, decltype(mv), Extend::Nodes<SampleId>,
                                   Extend::Empty<>, Extend::Empty<>, DefaultViewBase,
                                   IdContinuity::Sparse>::FromView(mv);
  merge_mv.View().SampleIdsFromCG();
  // std::cout << "\n\nCondensed MAT view\n";
  // MADAGToDOT(merge_mv.View(), std::cout);

  // TODO: check that mat and mat_view(condensed=True) are equal
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

void test_sample_dag1() {
  MADAGStorage dag_storage = make_sample_dag();
  test_mat_view_impl(dag_storage.View());
}

void test_sample_dag2() {
  MADAGStorage dag_storage = make_sample_dag_non_extremal_ua_node_id();
  test_mat_view_impl(dag_storage.View());
}

void test_sample_dag3() {
  MADAGStorage dag_storage = make_sample_dag_maximal_node_id();
  test_mat_view_impl(dag_storage.View());
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[]() { test_sample_dag1(); }, "MATView: sample dag"});

[[maybe_unused]] static const auto test_added1 = add_test(
    {[]() { test_sample_dag2(); }, "MATView: sample dag non-extremal UA NodeId"});

[[maybe_unused]] static const auto test_added2 =
    add_test({[]() { test_sample_dag3(); }, "MATView: sample dag maximal UA NodeId"});

[[maybe_unused]] static const auto test_added4 = add_test(
    {[]() {
       test_mat_view("data/seedtree/seedtree.pb.gz", "data/seedtree/refseq.txt.gz", "");
     },
     "MATView: seedtree"});

[[maybe_unused]] static const auto test_added5 =
    add_test({[]() { test_condensed_mat_view(); }, "MATView: condensing"});
#endif
