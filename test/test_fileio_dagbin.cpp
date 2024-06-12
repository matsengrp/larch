#include "test_common.hpp"
#include "test_common_dag.hpp"

#include "larch/benchmark.hpp"
#include "larch/dag_loader.hpp"

[[maybe_unused]] bool compare_files(const std::string& file1,
                                    const std::string& file2) {
  std::ifstream stream1(file1, std::ios::binary);
  std::ifstream stream2(file2, std::ios::binary);

  if (!stream1 || !stream2) {
    return false;
  }

  char byte1, byte2;

  while (stream1.get(byte1) && stream2.get(byte2)) {
    if (byte1 != byte2) {
      return false;
    }
  }

  if (!stream1.eof() || !stream2.eof()) {
    return false;
  }

  return true;
}

[[maybe_unused]] static void test_dagbin_sample_dag() {
  std::string protobuf_path = test_output_folder + "/sample_dag.pb";
  std::string dagbin_path = test_output_folder + "/sample_dag.dagbin";

  auto sample_dag = make_sample_dag();
  dag_info(sample_dag);

  // Compare read/write dag protobuf vs dagbin
  std::cout << "store dags..." << std::endl;
  StoreDAGToProtobuf(sample_dag.View(), protobuf_path);
  StoreDAGToDagbin(sample_dag.View(), dagbin_path);
  std::cout << "load dags..." << std::endl;
  auto sample_dag_from_protobuf = LoadDAGFromProtobuf(protobuf_path);
  auto sample_dag_from_dagbin = LoadDAGFromDagbin(dagbin_path);
  std::cout << "compare dags..." << std::endl;
  TestAssert(compare_treedags(sample_dag_from_dagbin.View(),
                              sample_dag_from_protobuf.View()) &&
             "Loading Sample DAG via Protobuf and Dagbin do not have the same result.");

  // Compare dag to sample dag.
  sample_dag_from_protobuf.View().RecomputeCompactGenomes();
  TestAssert(compare_treedags(sample_dag_from_protobuf.View(), sample_dag.View()) &&
             "Loaded Protobuf DAG does not match Sample DAG.");
  sample_dag_from_dagbin.View().RecomputeCompactGenomes();
  TestAssert(compare_treedags(sample_dag_from_dagbin.View(), sample_dag.View()) &&
             "Loaded Dagbin DAG does not match Sample DAG.");

  // Append test.
  auto big_sample_dag = make_big_sample_dag_topology(true);
  StoreDAGToDagbin(big_sample_dag.View(), dagbin_path);
  auto big_sample_dag_incomplete = LoadDAGFromDagbin(dagbin_path);
  big_sample_dag_incomplete.View().RecomputeCompactGenomes();
  TestAssert(
      compare_treedags(big_sample_dag_incomplete.View(), big_sample_dag.View()) &&
      "Loaded Dagbin DAG does not match Incomplete Big Sample DAG.");

  auto edge_id = big_sample_dag.View().GetEdgesCount();
  big_sample_dag.View().AddEdge({edge_id++}, {15}, {13}, {1});
  big_sample_dag.View().AddEdge({edge_id++}, {14}, {12}, {1});
  big_sample_dag.View().AddEdge({edge_id++}, {11}, {9}, {1});
  big_sample_dag.View().AddEdge({edge_id++}, {10}, {8}, {1});

  StoreDAGToDagbin(big_sample_dag.View(), dagbin_path, true);
  auto big_sample_dag_complete = LoadDAGFromDagbin(dagbin_path);
  big_sample_dag_complete.View().RecomputeCompactGenomes();
  TestAssert(compare_treedags(big_sample_dag_complete.View(), big_sample_dag.View()) &&
             "Using append, loaded Dagbin DAG does not match Complete Big Sample DAG.");
}

[[maybe_unused]] static void test_dagbin_via_larchusher(std::string_view input_dag_path,
                                                        int iter, bool use_seed,
                                                        bool save_both) {
  std::string output_dag_path_no_ext = test_output_folder + "/temp";
  std::string output_dag_path_protobuf = test_output_folder + "/temp.pb";
  std::string output_dag_path_dagbin = test_output_folder + "/temp.dagbin";
  std::string inter_dag_path_protobuf = output_dag_path_protobuf + ".intermediate";
  std::string inter_dag_path_dagbin = output_dag_path_dagbin + ".intermediate";
  std::string other_options = "--thread 1 ";
  if (use_seed) {
    other_options += "--seed 42 ";
  }

  if (save_both) {
    other_options += " --output-format debug-all";
    auto [command, result] = run_larch_usher(input_dag_path, output_dag_path_no_ext,
                                             std::nullopt, iter, other_options);
    TestAssert((result == 0) && "larch-usher debug-all run failed.");
  } else {
    auto [command1, result1] = run_larch_usher(input_dag_path, output_dag_path_protobuf,
                                               std::nullopt, iter, other_options);
    TestAssert((result1 == 0) && "larch-usher protobuf run failed.");
    auto [command2, result2] = run_larch_usher(input_dag_path, output_dag_path_dagbin,
                                               std::nullopt, iter, other_options);
    TestAssert((result2 == 0) && "larch-usher dagbin run failed.");
  }

  // Check files formats differ.
  TestAssert(!compare_files(output_dag_path_protobuf, output_dag_path_dagbin) &&
             "final dag files both incorrectly written in the same format.");
  TestAssert(!compare_files(inter_dag_path_protobuf, inter_dag_path_dagbin) &&
             "intermediate dag files both incorrectly written in the same format.");

  // Check final DAGs
  auto dag_protobuf = LoadDAGFromProtobuf(output_dag_path_protobuf);
  auto dag_dagbin = LoadDAGFromDagbin(output_dag_path_dagbin);
  std::cout << "\nFINAL DAGS: " << std::endl;
  dag_info(dag_protobuf);
  dag_info(dag_dagbin);

  // Check intermediate DAGs
  auto inter_dag_protobuf = LoadDAGFromProtobuf(inter_dag_path_protobuf);
  auto inter_dag_dagbin = LoadDAGFromDagbin(inter_dag_path_dagbin);
  std::cout << "INTER DAGS: " << std::endl;
  dag_info(inter_dag_protobuf);
  dag_info(inter_dag_dagbin);

  TestAssert(compare_treedags(dag_protobuf.View(), dag_dagbin.View()) &&
             "larch-usher final DAG: Protobuf and Dagbin do not have same result.");
  TestAssert(
      compare_treedags(inter_dag_protobuf.View(), inter_dag_dagbin.View()) &&
      "larch-usher intermediate DAG: Protobuf and Dagbin do not have same result.");
}

[[maybe_unused]] static void test_dagbin(
    std::string_view input_dag_path, FileFormat input_format,
    std::optional<std::string> refseq_path = std::nullopt) {
  std::string output_dag_path_protobuf = test_output_folder + "/temp.pb";
  std::string output_dag_path_dagbin = test_output_folder + "/temp.dagbin";

  auto dag = LoadDAG(input_dag_path, input_format, refseq_path);

  Benchmark timer;
  timer.start();
  StoreDAG(dag.View(), output_dag_path_protobuf, FileFormat::ProtobufDAG);
  std::cout << "DAG stored as protobuf in: " << timer.lapFormatUs() << std::endl;
  timer.start();
  StoreDAG(dag.View(), output_dag_path_dagbin, FileFormat::Dagbin);
  std::cout << "DAG stored as dagbin in: " << timer.lapFormatUs() << std::endl;

  timer.start();
  auto dag_protobuf = LoadDAG(output_dag_path_protobuf, FileFormat::ProtobufDAG);
  std::cout << "DAG loaded from protobuf in: " << timer.lapFormatUs() << std::endl;
  timer.start();
  auto dag_dagbin = LoadDAG(output_dag_path_dagbin, FileFormat::Dagbin);
  std::cout << "DAG loaded from dagbin in: " << timer.lapFormatUs() << std::endl;

  std::cout << "\nDAG_INFO: " << std::endl;
  dag_info(dag);
  dag_info(dag_protobuf);
  dag_info(dag_dagbin);

  TestAssert(compare_treedags(dag.View(), dag_protobuf.View()) &&
             "Test dagbin: Original DAG and Protobuf do not have same result.");
  TestAssert(compare_treedags(dag.View(), dag_dagbin.View()) &&
             "Test dagbin: Original DAG and Dagbin do not have same result.");
}

[[maybe_unused]] static const auto test_added0 =
    add_test({test_dagbin_sample_dag, "Load dagbin: Sample DAG"});

const std::string input_dag_path = "data/test_5_trees/full_dag.pb.gz";

[[maybe_unused]] static const auto test_added1 =
    add_test({[] { test_dagbin(input_dag_path, FileFormat::ProtobufDAG); },
              "Load dagbin: test_5_trees"});
[[maybe_unused]] static const auto test_added2 =
    add_test({[] { test_dagbin_via_larchusher(input_dag_path, 3, true, true); },
              "Load dagbin: test_5_trees, 3 iters, simultaneous write"});
[[maybe_unused]] static const auto test_added3 =
    add_test({[] { test_dagbin_via_larchusher(input_dag_path, 3, true, false); },
              "Load dagbin: test_5_trees, 3 iters, seeded"});

const std::string big_input_dag_path = "data/big_test/big_test.pb.gz";

[[maybe_unused]] static const auto test_added4 =
    add_test({[] { test_dagbin(big_input_dag_path, FileFormat::ProtobufDAG); },
              "Load dagbin: big_test",
              {"slow"}});
