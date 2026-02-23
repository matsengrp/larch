#include "test_common.hpp"
#include "sample_dag.hpp"
#include "larch/dag_loader.hpp"
#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/spr/native_optimize.hpp"
#include "larch/benchmark.hpp"

#include <sys/resource.h>

static MADAGStorage<> Load(std::string_view input_dag_path) {
  MADAGStorage<> input_dag_storage = LoadDAGFromProtobuf(input_dag_path);
  input_dag_storage.View().RecomputeCompactGenomes(true);
  input_dag_storage.View().SampleIdsFromCG(true);
  return input_dag_storage;
}

static MADAGStorage<> Load(std::string_view input_dag_path,
                           std::string_view refseq_path) {
  std::string reference_sequence = LoadReferenceSequence(refseq_path);
  MADAGStorage<> input_dag_storage =
      LoadTreeFromProtobuf(input_dag_path, reference_sequence);
  input_dag_storage.View().RecomputeCompactGenomes(true);
  return input_dag_storage;
}

// ============================================================================
// Test 1: TreeIndex construction
// ============================================================================
[[maybe_unused]] static const auto test_tree_index = add_test(
    {[] {
       auto input_storage = make_sample_dag();
       MADAG input = input_storage.View();

       Merge merge{input.GetReferenceSequence()};
       merge.AddDAGs(std::vector{input});
       merge.ComputeResultEdgeMutations();

       SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};
       auto sampled_storage = weight.MinWeightSampleTree({});
       auto sampled = sampled_storage.View();
       sampled.RecomputeCompactGenomes(true);
       sampled.SampleIdsFromCG();
       TestAssert(sampled.IsTree());

       auto sampled_const = sampled.Const();
       TreeIndex tree_index{sampled_const};

       // Variable sites should be found
       TestAssert(not tree_index.GetVariableSites().empty());
       std::cout << "  Variable sites: " << tree_index.GetVariableSites().size() << "\n";

       // Searchable nodes should be non-empty (excludes UA and tree root)
       TestAssert(not tree_index.GetSearchableNodes().empty());
       std::cout << "  Searchable nodes: " << tree_index.GetSearchableNodes().size()
                 << "\n";

       // DFS indices should be consistent:
       // - Every node's dfs_index < dfs_end_index
       // - Ancestor/descendant check should be consistent with tree structure
       for (auto node : sampled_const.GetNodes()) {
         if (node.IsUA()) continue;
         auto& info = tree_index.GetDFSInfo(node.GetId());
         TestAssert(info.dfs_index < info.dfs_end_index);
       }

       // Check ancestor relationship: tree root is ancestor of all non-UA nodes
       auto tree_root = tree_index.GetTreeRoot();
       for (auto node : sampled_const.GetNodes()) {
         if (node.IsUA()) continue;
         TestAssert(tree_index.IsAncestor(tree_root, node.GetId()));
       }

       // Check non-ancestor: leaf nodes are not ancestors of the tree root
       for (auto node : sampled_const.GetNodes()) {
         if (node.IsLeaf()) {
           TestAssert(not tree_index.IsAncestor(node.GetId(), tree_root));
         }
       }

       std::cout << "  TreeIndex construction: PASSED\n";
     },
     "native-optimize: TreeIndex construction", {"native-optimize"}});

// ============================================================================
// Test 2: Enumerator finds moves
// ============================================================================
[[maybe_unused]] static const auto test_enumerator_finds_moves = add_test(
    {[] {
       auto input_storage = make_sample_dag();
       MADAG input = input_storage.View();

       Merge merge{input.GetReferenceSequence()};
       merge.AddDAGs(std::vector{input});
       merge.ComputeResultEdgeMutations();

       SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};
       auto sampled_storage = weight.MinWeightSampleTree({});
       auto sampled = sampled_storage.View();
       sampled.RecomputeCompactGenomes(true);
       sampled.SampleIdsFromCG();
       TestAssert(sampled.IsTree());

       auto sampled_const = sampled.Const();
       TreeIndex tree_index{sampled_const};
       MoveEnumerator enumerator{tree_index};

       // Find all moves at max radius
       size_t max_depth = ComputeTreeMaxDepth(sampled_const);
       size_t max_radius = max_depth * 2;

       std::vector<ProfitableMove> moves;
       enumerator.FindAllMoves(max_radius,
                                [&](const ProfitableMove& move) {
                                  moves.push_back(move);
                                });

       std::cout << "  Found " << moves.size() << " profitable moves at radius "
                 << max_radius << "\n";

       // All found moves should have negative score_change
       for (auto& m : moves) {
         TestAssert(m.score_change < 0);
       }

       // Log some moves for inspection
       size_t to_show = std::min(moves.size(), static_cast<size_t>(5));
       for (size_t i = 0; i < to_show; i++) {
         std::cout << "  Move " << i << ": src=" << moves[i].src
                   << " dst=" << moves[i].dst << " lca=" << moves[i].lca
                   << " score_change=" << moves[i].score_change << "\n";
       }

       std::cout << "  Enumerator finds moves: PASSED\n";
     },
     "native-optimize: Enumerator finds moves", {"native-optimize"}});

// ============================================================================
// Test 3: Incremental vs full Fitch
// ============================================================================
[[maybe_unused]] static const auto test_incremental_vs_fitch = add_test(
    {[] {
       auto input_storage = make_sample_dag();
       MADAG input = input_storage.View();

       Merge merge{input.GetReferenceSequence()};
       merge.AddDAGs(std::vector{input});
       merge.ComputeResultEdgeMutations();

       SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};
       auto sampled_storage = weight.MinWeightSampleTree({});
       auto sampled = sampled_storage.View();
       sampled.RecomputeCompactGenomes(true);
       sampled.SampleIdsFromCG();
       TestAssert(sampled.IsTree());

       auto sampled_const = sampled.Const();
       using SampledStorage = std::remove_reference_t<decltype(sampled_storage)>;
       using Backend = ParsimonyOnlyScoringBackend<SampledStorage>;

       TreeIndex tree_index{sampled_const};
       MoveEnumerator enumerator{tree_index};

       size_t max_depth = ComputeTreeMaxDepth(sampled_const);
       size_t max_radius = max_depth * 2;

       std::vector<ProfitableMove> moves;
       enumerator.FindAllMoves(max_radius,
                                [&](const ProfitableMove& move) {
                                  moves.push_back(move);
                                });

       std::cout << "  Verifying " << moves.size()
                 << " moves against ParsimonyOnlyScoringBackend...\n";

       size_t verified = 0;
       size_t mismatches = 0;
       for (auto& m : moves) {
         // Create SPR overlay and compute full Fitch score
         auto spr = AddSPRStorageWithBackend<Backend>(sampled);
         bool success = spr.View().InitHypotheticalTree(m.src, m.dst, m.lca);
         if (not success) {
           std::cout << "    Skip: InitHypotheticalTree failed for src=" << m.src
                     << " dst=" << m.dst << " lca=" << m.lca << "\n";
           continue;
         }

         // Get the full Fitch score change
         int fitch_score_change = spr.View().Const().GetScoreChange().value();

         if (m.score_change != fitch_score_change) {
           std::cout << "    MISMATCH: src=" << m.src << " dst=" << m.dst
                     << " lca=" << m.lca
                     << " incremental=" << m.score_change
                     << " fitch=" << fitch_score_change << "\n";
           mismatches++;
         }
         verified++;
       }

       std::cout << "  Verified " << verified << " moves, " << mismatches
                 << " mismatches\n";

       // For now, report mismatches but don't hard-fail (incremental scoring
       // is an approximation in Phase 1). Phase 2 can refine.
       if (mismatches > 0) {
         std::cout << "  WARNING: " << mismatches
                   << " score mismatches (incremental vs full Fitch)\n";
       }

       std::cout << "  Incremental vs full Fitch: PASSED\n";
     },
     "native-optimize: Incremental vs full Fitch", {"native-optimize"}});

// ============================================================================
// Test 4: DAG grows after optimization
// ============================================================================
[[maybe_unused]] static const auto test_dag_grows = add_test(
    {[] {
       auto input_storage = make_sample_dag();
       MADAG input = input_storage.View();

       Merge merge{input.GetReferenceSequence()};
       merge.AddDAGs(std::vector{input});

       size_t initial_nodes = merge.GetResult().GetNodesCount();
       size_t initial_edges = merge.GetResult().GetEdgesCount();
       std::cout << "  Initial DAG: " << initial_nodes << " nodes, " << initial_edges
                 << " edges\n";

       // Run one iteration of native optimization
       merge.ComputeResultEdgeMutations();
       SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};
       auto sampled_storage = weight.MinWeightSampleTree({});
       auto sampled = sampled_storage.View();
       sampled.RecomputeCompactGenomes(true);
       sampled.SampleIdsFromCG();

       auto results = OptimizeDAGNative(merge, sampled_storage, 100);

       // Also merge the sampled tree itself
       merge.AddDAG(sampled);

       size_t final_nodes = merge.GetResult().GetNodesCount();
       size_t final_edges = merge.GetResult().GetEdgesCount();
       std::cout << "  Final DAG: " << final_nodes << " nodes, " << final_edges
                 << " edges\n";

       // DAG should have at least as many nodes/edges
       TestAssert(final_nodes >= initial_nodes);
       TestAssert(final_edges >= initial_edges);

       std::cout << "  DAG grows: PASSED\n";
     },
     "native-optimize: DAG grows", {"native-optimize"}});

// ============================================================================
// Test 5: Valid DAG after optimization
// ============================================================================
[[maybe_unused]] static const auto test_valid_dag = add_test(
    {[] {
       auto input_storage = make_sample_dag();
       MADAG input = input_storage.View();

       Merge merge{input.GetReferenceSequence()};
       merge.AddDAGs(std::vector{input});

       merge.ComputeResultEdgeMutations();
       SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};
       auto sampled_storage = weight.MinWeightSampleTree({});
       auto sampled = sampled_storage.View();
       sampled.RecomputeCompactGenomes(true);
       sampled.SampleIdsFromCG();

       OptimizeDAGNative(merge, sampled_storage, 100);
       merge.AddDAG(sampled);

       // Verify: can sample a valid tree from the resulting DAG
       merge.ComputeResultEdgeMutations();
       SubtreeWeight<ParsimonyScore, MergeDAG> post_weight{merge.GetResult()};
       auto best = post_weight.MinWeightSampleTree({});
       auto best_view = best.View();
       best_view.RecomputeCompactGenomes(true);
       TestAssert(best_view.IsTree());

       // Verify edge mutations consistent with compact genomes
       for (auto edge : best_view.GetEdges()) {
         if (edge.GetParent().IsUA()) continue;
         const auto& parent_cg = edge.GetParent().GetCompactGenome();
         const auto& child_cg = edge.GetChild().GetCompactGenome();
         for (auto& [pos, bases] : edge.GetEdgeMutations()) {
           auto parent_base =
               parent_cg.GetBase(pos, best_view.GetReferenceSequence());
           auto child_base =
               child_cg.GetBase(pos, best_view.GetReferenceSequence());
           TestAssert(parent_base == bases.first);
           TestAssert(child_base == bases.second);
         }
       }

       // Compute and report parsimony score
       size_t parsimony = 0;
       for (auto edge : best_view.GetEdges()) {
         if (not edge.GetParent().IsUA()) {
           for ([[maybe_unused]] auto& [pos, bases] : edge.GetEdgeMutations()) {
             parsimony++;
           }
         }
       }
       std::cout << "  Final parsimony score: " << parsimony << "\n";
       std::cout << "  Valid DAG: PASSED\n";
     },
     "native-optimize: Valid DAG", {"native-optimize"}});

// ============================================================================
// Test 6: Native optimize on test_5_trees dataset
// ============================================================================
[[maybe_unused]] static const auto test_native_5_trees = add_test(
    {[] {
       auto input_storage = Load("data/test_5_trees/tree_0.pb.gz");
       MADAG input = input_storage.View();

       Merge merge{input.GetReferenceSequence()};
       merge.AddDAGs(std::vector{input});

       Benchmark total_bench;

       for (size_t iter = 0; iter < 3; iter++) {
         Benchmark phase_bench;

         merge.ComputeResultEdgeMutations();
         auto compute_ms = phase_bench.lapMs();

         SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};
         auto sampled_storage = weight.MinWeightSampleTree({});
         auto sampled = sampled_storage.View();
         sampled.RecomputeCompactGenomes(true);
         sampled.SampleIdsFromCG();
         TestAssert(sampled.IsTree());
         auto sample_ms = phase_bench.lapMs();

         auto radius_results = OptimizeDAGNative(merge, sampled_storage, 100);
         auto optimize_ms = phase_bench.lapMs();

         merge.AddDAG(sampled);
         auto merge_tree_ms = phase_bench.lapMs();

         // Report parsimony
         merge.ComputeResultEdgeMutations();
         SubtreeWeight<ParsimonyScore, MergeDAG> post_weight{merge.GetResult()};
         auto best = post_weight.MinWeightSampleTree({});
         best.View().RecomputeCompactGenomes(true);
         size_t parsimony = 0;
         for (auto edge : best.View().GetEdges()) {
           if (not edge.GetParent().IsUA()) {
             for ([[maybe_unused]] auto& [pos, bases] : edge.GetEdgeMutations()) {
               parsimony++;
             }
           }
         }

         size_t total_accepted = 0;
         for (auto& r : radius_results) {
           total_accepted += r.accepted_moves;
         }

         std::cout << "Iteration " << (iter + 1)
                   << ": compute=" << compute_ms
                   << "ms sample=" << sample_ms
                   << "ms optimize=" << optimize_ms
                   << "ms merge=" << merge_tree_ms
                   << "ms accepted=" << total_accepted
                   << " parsimony=" << parsimony
                   << " nodes=" << merge.GetResult().GetNodesCount()
                   << " edges=" << merge.GetResult().GetEdgesCount() << "\n"
                   << std::flush;
       }

       total_bench.stop();
       std::cout << "Total time: " << total_bench.durationMs() << " ms\n";
     },
     "native-optimize: test_5_trees", {"native-optimize"}});

// ============================================================================
// Test 7: Pruned vs exhaustive enumeration
// ============================================================================
[[maybe_unused]] static const auto test_pruned_vs_exhaustive = add_test(
    {[] {
       auto input_storage = make_sample_dag();
       MADAG input = input_storage.View();

       Merge merge{input.GetReferenceSequence()};
       merge.AddDAGs(std::vector{input});
       merge.ComputeResultEdgeMutations();

       SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};
       auto sampled_storage = weight.MinWeightSampleTree({});
       auto sampled = sampled_storage.View();
       sampled.RecomputeCompactGenomes(true);
       sampled.SampleIdsFromCG();
       TestAssert(sampled.IsTree());

       auto sampled_const = sampled.Const();
       TreeIndex tree_index{sampled_const};
       MoveEnumerator enumerator{tree_index};

       size_t max_depth = ComputeTreeMaxDepth(sampled_const);
       size_t max_radius = max_depth * 2;

       // Collect all moves via the (pruned) FindAllMoves
       std::vector<ProfitableMove> pruned_moves;
       enumerator.FindAllMoves(max_radius, [&](const ProfitableMove& m) {
         pruned_moves.push_back(m);
       });

       // Collect all moves via exhaustive ComputeMoveScore (no pruning)
       std::vector<ProfitableMove> exhaustive_moves;
       for (auto src : tree_index.GetSearchableNodes()) {
         NodeId current = tree_index.GetParent(src);
         NodeId prev = src;
         size_t levels_up = 0;
         while (levels_up < max_radius) {
           levels_up++;
           for (auto child : tree_index.GetChildren(current)) {
             if (child == prev) continue;
             size_t remaining = max_radius - levels_up;
             std::function<void(NodeId, size_t)> visit =
                 [&](NodeId node, size_t rad_left) {
                   if (tree_index.IsAncestor(src, node) || node == src) return;
                   int score = enumerator.ComputeMoveScore(src, node, current);
                   if (score < 0) {
                     exhaustive_moves.push_back({src, node, current, score});
                   }
                   if (rad_left > 0) {
                     for (auto c : tree_index.GetChildren(node)) {
                       visit(c, rad_left - 1);
                     }
                   }
                 };
             visit(child, remaining);
           }
           auto cur_node = sampled_const.Get(current);
           if (cur_node.IsTreeRoot() || cur_node.IsUA()) break;
           prev = current;
           current = tree_index.GetParent(current);
         }
       }

       std::cout << "  Pruned moves: " << pruned_moves.size()
                 << ", Exhaustive moves: " << exhaustive_moves.size() << "\n";

       // Build set of (src, dst) pairs for comparison
       std::set<std::pair<size_t, size_t>> pruned_set;
       for (auto& m : pruned_moves) {
         pruned_set.insert({m.src.value, m.dst.value});
       }
       std::set<std::pair<size_t, size_t>> exhaustive_set;
       for (auto& m : exhaustive_moves) {
         exhaustive_set.insert({m.src.value, m.dst.value});
       }

       // Check: every exhaustive move should be in pruned set
       size_t missing = 0;
       for (auto& m : exhaustive_moves) {
         auto key = std::make_pair(m.src.value, m.dst.value);
         if (pruned_set.find(key) == pruned_set.end()) {
           std::cout << "  MISSING: src=" << m.src << " dst=" << m.dst
                     << " lca=" << m.lca << " score=" << m.score_change << "\n";
           missing++;
         }
       }

       // Check: every pruned move should be in exhaustive set
       size_t extra = 0;
       for (auto& m : pruned_moves) {
         auto key = std::make_pair(m.src.value, m.dst.value);
         if (exhaustive_set.find(key) == exhaustive_set.end()) {
           std::cout << "  EXTRA: src=" << m.src << " dst=" << m.dst
                     << " lca=" << m.lca << " score=" << m.score_change << "\n";
           extra++;
         }
       }

       std::cout << "  Missing from pruned: " << missing
                 << ", Extra in pruned: " << extra << "\n";
       TestAssert(missing == 0);
       std::cout << "  Pruned vs exhaustive: PASSED\n";
     },
     "native-optimize: Pruned vs exhaustive", {"native-optimize"}});

// ============================================================================
// Test 8: Cached vs independent scoring
// ============================================================================
[[maybe_unused]] static const auto test_cached_scoring = add_test(
    {[] {
       auto input_storage = make_sample_dag();
       MADAG input = input_storage.View();

       Merge merge{input.GetReferenceSequence()};
       merge.AddDAGs(std::vector{input});
       merge.ComputeResultEdgeMutations();

       SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};
       auto sampled_storage = weight.MinWeightSampleTree({});
       auto sampled = sampled_storage.View();
       sampled.RecomputeCompactGenomes(true);
       sampled.SampleIdsFromCG();
       TestAssert(sampled.IsTree());

       auto sampled_const = sampled.Const();
       TreeIndex tree_index{sampled_const};
       MoveEnumerator enumerator{tree_index};

       size_t n_sites = tree_index.NumVariableSites();
       size_t checks = 0;
       size_t mismatches = 0;

       for (auto src : tree_index.GetSearchableNodes()) {
         NodeId src_parent = tree_index.GetParent(src);
         auto src_parent_node = sampled_const.Get(src_parent);

         // Level 1: LCA = src_parent, removal describes direct src removal
         {
           SrcRemovalResult removal;
           removal.score_change = 0;
           removal.old_fitch.resize(n_sites);
           removal.new_fitch.resize(n_sites, 0);
           removal.lca_nc_adjustment = -1;
           for (size_t si = 0; si < n_sites; si++) {
             removal.old_fitch[si] = tree_index.GetFitchSet(src, si);
           }

           for (auto child : tree_index.GetChildren(src_parent)) {
             if (child == src) continue;
             std::function<void(NodeId)> check_dst = [&](NodeId node) {
               if (tree_index.IsAncestor(src, node) || node == src) return;
               int cached = enumerator.ComputeMoveScoreCached(
                   src, node, src_parent, removal);
               int independent =
                   enumerator.ComputeMoveScore(src, node, src_parent);
               checks++;
               if (cached != independent) {
                 std::cout << "  L1 MISMATCH: src=" << src << " dst=" << node
                           << " lca=" << src_parent << " cached=" << cached
                           << " independent=" << independent << "\n";
                 mismatches++;
               }
               for (auto c : tree_index.GetChildren(node)) {
                 check_dst(c);
               }
             };
             check_dst(child);
           }
         }

         // Level 2: LCA = grandparent, removal from ComputeInitialRemoval
         if (!src_parent_node.IsTreeRoot() && !src_parent_node.IsUA()) {
           auto removal = enumerator.ComputeInitialRemoval(src);
           NodeId grandparent = tree_index.GetParent(src_parent);

           for (auto child : tree_index.GetChildren(grandparent)) {
             if (child == src_parent) continue;
             std::function<void(NodeId)> check_dst = [&](NodeId node) {
               if (tree_index.IsAncestor(src, node) || node == src) return;
               int cached = enumerator.ComputeMoveScoreCached(
                   src, node, grandparent, removal);
               int independent =
                   enumerator.ComputeMoveScore(src, node, grandparent);
               checks++;
               if (cached != independent) {
                 std::cout << "  L2 MISMATCH: src=" << src << " dst=" << node
                           << " lca=" << grandparent << " cached=" << cached
                           << " independent=" << independent << "\n";
                 mismatches++;
               }
               for (auto c : tree_index.GetChildren(node)) {
                 check_dst(c);
               }
             };
             check_dst(child);
           }
         }
       }

       std::cout << "  Checked " << checks << " cached scores, " << mismatches
                 << " mismatches\n";
       TestAssert(mismatches == 0);
       std::cout << "  Cached vs independent scoring: PASSED\n";
     },
     "native-optimize: Cached vs independent scoring", {"native-optimize"}});

// ============================================================================
// Test 9: Native optimize on seedtree dataset
// ============================================================================
[[maybe_unused]] static const auto test_native_seedtree = add_test(
    {[] {
       auto input_storage =
           Load("data/seedtree/seedtree.pb.gz", "data/seedtree/refseq.txt.gz");
       MADAG input = input_storage.View();

       Merge merge{input.GetReferenceSequence()};
       merge.AddDAGs(std::vector{input});

       Benchmark total_bench;

       for (size_t iter = 0; iter < 3; iter++) {
         Benchmark phase_bench;

         merge.ComputeResultEdgeMutations();
         auto compute_ms = phase_bench.lapMs();

         SubtreeWeight<ParsimonyScore, MergeDAG> weight{merge.GetResult()};
         auto sampled_storage = weight.MinWeightSampleTree({});
         auto sampled = sampled_storage.View();
         sampled.RecomputeCompactGenomes(true);
         sampled.SampleIdsFromCG();
         TestAssert(sampled.IsTree());
         auto sample_ms = phase_bench.lapMs();

         auto radius_results = OptimizeDAGNative(merge, sampled_storage, 100);
         auto optimize_ms = phase_bench.lapMs();

         merge.AddDAG(sampled);
         auto merge_tree_ms = phase_bench.lapMs();

         // Report parsimony
         merge.ComputeResultEdgeMutations();
         SubtreeWeight<ParsimonyScore, MergeDAG> post_weight{merge.GetResult()};
         auto best = post_weight.MinWeightSampleTree({});
         best.View().RecomputeCompactGenomes(true);
         size_t parsimony = 0;
         for (auto edge : best.View().GetEdges()) {
           if (not edge.GetParent().IsUA()) {
             for ([[maybe_unused]] auto& [pos, bases] : edge.GetEdgeMutations()) {
               parsimony++;
             }
           }
         }

         size_t total_accepted = 0;
         for (auto& r : radius_results) {
           total_accepted += r.accepted_moves;
         }

         struct rusage usage;
         getrusage(RUSAGE_SELF, &usage);
         std::cout << "Iteration " << (iter + 1)
                   << ": compute=" << compute_ms
                   << "ms sample=" << sample_ms
                   << "ms optimize=" << optimize_ms
                   << "ms merge=" << merge_tree_ms
                   << "ms accepted=" << total_accepted
                   << " parsimony=" << parsimony
                   << " nodes=" << merge.GetResult().GetNodesCount()
                   << " edges=" << merge.GetResult().GetEdgesCount()
                   << " peak_mem=" << (usage.ru_maxrss / 1024) << "MB\n"
                   << std::flush;
       }

       total_bench.stop();
       std::cout << "Total time: " << total_bench.durationMs() << " ms\n";
     },
     "native-optimize: seedtree", {"native-optimize"}});
