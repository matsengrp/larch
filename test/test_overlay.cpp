#include "test_common.hpp"
#include "larch/dag_loader.hpp"

static void test_overlay(std::string_view input_dag_path,
                         std::string_view refseq_path) {
  std::string reference_sequence = LoadReferenceSequence(refseq_path);
  MADAGStorage input_dag_storage =
      LoadTreeFromProtobuf(input_dag_path, reference_sequence);
  auto input_dag = input_dag_storage.View();
  input_dag.RecomputeCompactGenomes();
  auto overlay_dag_storage = AddOverlay(input_dag);
  auto overlay_dag = overlay_dag_storage.View();

  auto input_node = input_dag.Get(NodeId{10});
  auto overlay_node = overlay_dag.Get(NodeId{10});

  Assert(not input_node.GetCompactGenome().empty());
  Assert(not overlay_node.GetCompactGenome().empty());
  Assert(input_node.GetCompactGenome() == overlay_node.GetCompactGenome());

  overlay_node.SetOverlay<CompactGenome>();

  Assert(overlay_node.IsOverlaid<CompactGenome>());
  Assert(not overlay_node.GetCompactGenome().empty());
  Assert(input_node.GetCompactGenome() == overlay_node.GetCompactGenome());

  ContiguousMap<MutationPosition, char> new_cg;
  new_cg.insert({{5}, 'A'});
  new_cg.insert({{10}, 'G'});
  new_cg.insert({{15}, 'T'});
  new_cg.insert({{20}, 'C'});
  overlay_node = CompactGenome{std::move(new_cg)};

  Assert(not overlay_node.GetCompactGenome().empty());
  Assert(input_node.GetCompactGenome() != overlay_node.GetCompactGenome());
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] {
                test_overlay("data/20D_from_fasta/1final-tree-1.nh1.pb.gz",
                             "data/20D_from_fasta/refseq.txt.gz");
              },
              "Overlay"});
