#include "test_common.hpp"
#include "larch/dag_loader.hpp"
#if USE_MAT_VIEW
#include "larch/mat_view.hpp"
#endif
#include "sample_dag.hpp"

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

  auto overlay_connections_node = overlay_dag.Get(NodeId{7});
  auto overlay_connections_edge = overlay_connections_node.GetSingleParent();
  auto old_parent_node = overlay_connections_node.GetSingleParent().GetParent().GetId();
  auto new_parent_node = overlay_dag.Get(NodeId{10});
  auto new_parent_num_clades = new_parent_node.GetCladesCount();
  new_parent_node.SetOverlay<Neighbors>();
  overlay_connections_node.SetOverlay<Neighbors>();
  overlay_connections_edge.SetOverlay<Endpoints>();
  overlay_connections_edge.Set(new_parent_node, overlay_connections_node, {new_parent_num_clades});
  new_parent_node.AddEdge({new_parent_num_clades}, overlay_connections_edge, true);
  overlay_connections_node.SetSingleParent(overlay_connections_edge);

  TestAssert(overlay_connections_node.GetDAG().GetOriginal().Get(overlay_connections_node).GetSingleParent().GetParent().GetId() == old_parent_node)
}

#if USE_MAT_VIEW
using Storage = CondensedMADAGStorage;

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

  // test accessing children after SetOverlay
  auto node_7 = overlay_dag.Get(NodeId{7});
  TestAssert(node_7.GetCladesCount() == 2);
  node_7.SetOverlay<Neighbors>();
  TestAssert(node_7.GetCladesCount() == 2);
  TestAssert(node_7.ContainsChild(overlay_dag.Get(NodeId{1})));
  TestAssert(node_7.ContainsChild(overlay_dag.Get(NodeId{2})));

  // test overlay connectivity (edge Endpoints, node Neighbors)
  auto overlay_edge = overlay_dag.Get(EdgeId{3});
  auto overlay_child = overlay_edge.GetChild();
  auto overlay_new_parent = overlay_dag.Get(NodeId{10});

  auto overlay_old_parent = overlay_child.GetSingleParent().GetParent();
  auto overlay_old_parent_edge = overlay_old_parent.GetSingleParent();
  auto overlay_old_parent_edge_parent = overlay_old_parent_edge.GetParent();
  auto overlay_old_parent_clade = overlay_old_parent_edge.GetClade().value;
  auto overlay_child_sib_1 = overlay_dag.Get(NodeId{4});
  auto overlay_child_sib_2 = overlay_dag.Get(NodeId{7});
  auto overlay_child_sib_edge_1 = overlay_dag.Get(EdgeId{4});
  auto overlay_child_sib_edge_2 = overlay_dag.Get(EdgeId{7});

  overlay_edge.SetOverlay<Endpoints>();
  overlay_new_parent.SetOverlay<Neighbors>();
  overlay_child.SetOverlay<Neighbors>();

  overlay_old_parent.SetOverlay<Neighbors>();
  overlay_child_sib_1.SetOverlay<Neighbors>();
  // overlay_child_sib_2.SetOverlay<Neighbors>();
  overlay_old_parent_edge.SetOverlay<Endpoints>();
  overlay_child_sib_edge_1.SetOverlay<Endpoints>();
  overlay_child_sib_edge_2.SetOverlay<Endpoints>();

  TestAssert(overlay_edge.IsOverlaid<Endpoints>());
  TestAssert(overlay_new_parent.IsOverlaid<Neighbors>());

  overlay_child.ClearConnections();
  overlay_edge.Set(overlay_new_parent, overlay_child, {2});
  overlay_new_parent.AddEdge({2}, overlay_edge, true);
  overlay_child.SetSingleParent(overlay_edge);

  TestAssert(overlay_edge.GetParent() == overlay_new_parent);
  TestAssert(overlay_child.GetSingleParent().GetParent() == overlay_new_parent);

  TestAssert(overlay_edge.GetParent().GetParentsCount() == 1);
  TestAssert(overlay_edge.GetChild().GetParentsCount() == 1);
  TestAssert(overlay_edge.GetChild().template IsOverlaid<Neighbors>());
  TestAssert(overlay_edge.GetChild().ContainsParent(overlay_edge.GetParent()));
  TestAssert(overlay_edge.GetParent().ContainsChild(overlay_edge.GetChild()));

  // make the overlay_dag a valid tree-shaped DAG by removing the old parent's stored
  // connection to overlay_child
  overlay_child_sib_edge_1.Set(overlay_old_parent, overlay_child_sib_1, {0});
  overlay_child_sib_edge_2.Set(overlay_old_parent, overlay_child_sib_2, {1});
  overlay_old_parent.ClearConnections();
  overlay_old_parent_edge.Set(overlay_old_parent_edge_parent, overlay_old_parent,
                              {overlay_old_parent_clade});
  overlay_old_parent.SetSingleParent(overlay_old_parent_edge);
  overlay_old_parent.AddEdge({0}, overlay_child_sib_edge_1, true);
  overlay_old_parent.AddEdge({1}, overlay_child_sib_edge_2, true);
  overlay_child_sib_1.SetSingleParent(overlay_child_sib_edge_1);
  overlay_child_sib_2.SetSingleParent(overlay_child_sib_edge_2);
  // make sure the current overlay is a valid one
  overlay_dag.GetRoot().Validate(true, false);

  // test appending node/edge to overlaid dag.
  auto new_node = overlay_dag.AppendNode();
  auto new_edge = overlay_dag.AppendEdge();
  auto child_node_1 = overlay_dag.Get(NodeId({8}));
  auto child_node_2 = overlay_dag.Get(NodeId({9}));
  auto child_node_3 = overlay_dag.Get(NodeId({3}));
  auto parent_node = overlay_dag.Get(NodeId({10}));
  auto child_edge_1 = child_node_1.GetSingleParent();
  auto child_edge_2 = child_node_2.GetSingleParent();
  auto child_edge_3 = child_node_3.GetSingleParent();
  auto parent_edge = parent_node.GetSingleParent();

  if (not child_node_1.IsOverlaid<Neighbors>()) {
    child_node_1.SetOverlay<Neighbors>();
  }
  if (not child_node_2.IsOverlaid<Neighbors>()) {
    child_node_2.SetOverlay<Neighbors>();
  }
  if (not child_node_3.IsOverlaid<Neighbors>()) {
    child_node_3.SetOverlay<Neighbors>();
  }
  if (not parent_node.IsOverlaid<Neighbors>()) {
    parent_node.SetOverlay<Neighbors>();
  }
  if (not child_edge_1.IsOverlaid<Endpoints>()) {
    child_edge_1.SetOverlay<Endpoints>();
  }
  if (not child_edge_2.IsOverlaid<Endpoints>()) {
    child_edge_2.SetOverlay<Endpoints>();
  }
  if (not child_edge_3.IsOverlaid<Endpoints>()) {
    child_edge_3.SetOverlay<Endpoints>();
  }
  if (not parent_edge.IsOverlaid<Endpoints>()) {
    parent_edge.SetOverlay<Endpoints>();
  }

  // OPERATIONS:
  // - create a cherry of child_nodes 1 and 2 with new_node as the parent of this cherry
  // - and set new_node and child_node_3 as children of parent_node.
  // step 1: clear parent node
  parent_node.ClearConnections();
  // step 2: connect parent_node to its parent
  parent_node.SetSingleParent(parent_edge);
  // step 3: add child_node_3 as child of parent_node
  child_edge_3.Set(parent_node, child_node_3, {0});
  parent_node.AddEdge({0}, child_edge_3, true);
  child_node_3.SetSingleParent(child_edge_3);
  // step 4: add new_node as child of parent_node
  new_edge.Set(parent_node, new_node, {1});
  parent_node.AddEdge({1}, new_edge, true);
  new_node.SetSingleParent(new_edge);
  // step 5: add child_node_1 as child of new_node
  child_edge_1.Set(new_node, child_node_1, {0});
  new_node.AddEdge({0}, child_edge_1, true);
  child_node_1.SetSingleParent(child_edge_1);
  // step 6: add child_node_2 as child of new_node
  child_edge_2.Set(new_node, child_node_2, {1});
  new_node.AddEdge({1}, child_edge_2, true);
  child_node_2.SetSingleParent(child_edge_2);

  TestAssert(parent_node.GetParentsCount() == 1);
  TestAssert(new_node.GetParentsCount() == 1);
  TestAssert(child_node_1.GetParentsCount() == 1);
  TestAssert(child_node_2.GetParentsCount() == 1);
  TestAssert(child_node_3.GetParentsCount() == 1);
  TestAssert(not new_node.IsUA());
  TestAssert(not new_node.IsLeaf());
  TestAssert(new_node.GetCladesCount() == 2);
  TestAssert(new_node.ContainsChild(child_node_1));
  TestAssert(new_node.ContainsChild(child_node_2));
  TestAssert(child_node_1.ContainsParent(new_node));
  TestAssert(child_node_2.ContainsParent(new_node));
  TestAssert(new_edge.GetParent() == parent_node);
  TestAssert(new_edge.GetChild() == new_node);

  TestAssert(child_node_1.GetDAG().GetOriginal().Get(child_node_1).GetSingleParent().GetParent().GetId() == parent_node.GetId())

  // make sure the current overlay is a valid one
  overlay_dag.GetRoot().Validate(true, false);
}
#endif

[[maybe_unused]] static const auto test_added0 =
    add_test({[] {
                test_overlay_dag("data/20D_from_fasta/1final-tree-1.nh1.pb.gz",
                                 "data/20D_from_fasta/refseq.txt.gz");
              },
              "Overlay: DAG"});
#if USE_MAT_VIEW
[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_overlay_mat_view(); }, "Overlay: MATView"});
#endif
