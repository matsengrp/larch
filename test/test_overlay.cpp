#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#include "larch/mat_view.hpp"
#include "sample_dag.hpp"

using Storage = CondensedMADAGStorage;

static void test_overlay_dag(std::string_view input_dag_path,
                         std::string_view refseq_path) {
  std::string reference_sequence = LoadReferenceSequence(refseq_path);
  MADAGStorage input_dag_storage =
      LoadTreeFromProtobuf(input_dag_path, reference_sequence);
  auto input_dag = input_dag_storage.View();
  input_dag.RecomputeCompactGenomes(true);
  auto overlay_dag_storage = AddOverlay<void>(input_dag);
  auto overlay_dag = overlay_dag_storage.View();

  [[maybe_unused]] auto input_node = input_dag.Get(NodeId{10});
  auto overlay_node = overlay_dag.Get(NodeId{10});

  TestAssert(not input_node.GetCompactGenome().empty());
  TestAssert(not overlay_node.GetCompactGenome().empty());
  TestAssert(input_node.GetCompactGenome() == overlay_node.GetCompactGenome());

  overlay_node.SetOverlay<CompactGenome>();

  TestAssert(overlay_node.IsOverlaid<CompactGenome>());
  TestAssert(not overlay_node.GetCompactGenome().empty());
  TestAssert(input_node.GetCompactGenome() == overlay_node.GetCompactGenome());

  ContiguousMap<MutationPosition, MutationBase> new_cg;
  new_cg.insert({{5}, {'A'}});
  new_cg.insert({{10}, {'G'}});
  new_cg.insert({{15}, {'T'}});
  new_cg.insert({{20}, {'C'}});
  overlay_node = CompactGenome{std::move(new_cg)};

  TestAssert(not overlay_node.GetCompactGenome().empty());
  TestAssert(input_node.GetCompactGenome() != overlay_node.GetCompactGenome());
}

static void test_overlay_mat_view() {
  // create a MAT from which we can build a MATView
  auto dag_storage = make_sample_dag();
  auto dag = dag_storage.View();
  dag.SampleIdsFromCG();
  auto mat_conv = AddMATConversion(dag);
  MAT::Tree mat;
  mat_conv.View().BuildMAT(mat);
  mat_conv.View().GetRoot().Validate(true);
  std::vector<std::string> condense_arg{};
  mat.condense_leaves(condense_arg);
  mat.fix_node_idx();

  // Create MAT View
  CondensedMATViewStorage matview_storage;
  matview_storage.View().SetMAT(std::addressof(mat));
  auto storage = Storage::Consume(std::move(matview_storage));
  auto mv = storage.View();
  mv.SetReferenceSequence(dag.GetReferenceSequence());
  mv.BuildRootAndLeafs();
  mv.RecomputeCompactGenomes<IdContinuity::Sparse>();
  mv.GetRoot().Validate(true, false);

  auto overlay_dag_storage = AddOverlay<void>(mv);
  auto overlay_dag = overlay_dag_storage.View();

  // test overlay compact genome
  [[maybe_unused]] auto input_node = mv.Get(NodeId{2});
  auto overlay_node = overlay_dag.Get(NodeId{2});

  TestAssert(not input_node.GetCompactGenome().empty());
  TestAssert(not overlay_node.GetCompactGenome().empty());
  TestAssert(input_node.GetCompactGenome() == overlay_node.GetCompactGenome());

  overlay_node.SetOverlay<CompactGenome>();

  TestAssert(overlay_node.IsOverlaid<CompactGenome>());
  TestAssert(not overlay_node.GetCompactGenome().empty());
  TestAssert(input_node.GetCompactGenome() == overlay_node.GetCompactGenome());

  ContiguousMap<MutationPosition, MutationBase> new_cg;
  new_cg.insert({{1}, {'C'}});
  new_cg.insert({{2}, {'C'}});
  overlay_node = CompactGenome{std::move(new_cg)};

  TestAssert(not overlay_node.GetCompactGenome().empty());
  TestAssert(input_node.GetCompactGenome() != overlay_node.GetCompactGenome());

  // test overlay connectivity (edge Endpoints, node Neighbors)
  auto overlay_edge = overlay_dag.Get(EdgeId{3});
  auto overlay_child = overlay_edge.GetChild();
  auto overlay_new_parent = overlay_dag.Get(NodeId{10});

  overlay_edge.SetOverlay<Endpoints>();
  overlay_new_parent.SetOverlay<Neighbors>();

  TestAssert(overlay_edge.IsOverlaid<Endpoints>());
  TestAssert(overlay_new_parent.IsOverlaid<Neighbors>());

  overlay_edge.Set(overlay_new_parent, overlay_child, {2});
  overlay_new_parent.AddEdge({2}, overlay_edge, true);
  TestAssert(overlay_edge.GetParent() == overlay_new_parent);
  TestAssert(overlay_child.GetSingleParent().GetParent() == overlay_new_parent);
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] {
                test_overlay_dag("data/20D_from_fasta/1final-tree-1.nh1.pb.gz",
                             "data/20D_from_fasta/refseq.txt.gz");
              },
              "Overlay: DAG"});

[[maybe_unused]] static const auto test_added1 =
    add_test({[] {
                test_overlay_mat_view();
              },
              "Overlay: MATView"});
