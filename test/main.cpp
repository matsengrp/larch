#include <iostream>
#include <set>
#include <regex>
#include <chrono>

#ifdef USE_USHER
#include <mpi.h>
#endif

#include "test_common.hpp"

static std::vector<Test>& get_all_tests() {
  static std::vector<Test> all_tests{};
  return all_tests;
}

bool add_test(const Test& test) noexcept {
  get_all_tests().push_back(test);
  return true;
}

int main(int argc, char* argv[]) {
#ifdef USE_USHER
  int ignored{};
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &ignored);
#endif
  bool no_catch = false;
  bool opt_list_names = false;
  bool opt_test_range = false;

  std::pair<int, int> range;
  std::regex regex{".*"};
  std::set<std::string> include_tags;
  std::set<std::string> exclude_tags;
  for (int i = 1; i < argc; ++i) {
    if (std::string("nocatch") == argv[i]) {
      no_catch = true;
    } else if (std::string("--list") == argv[i]) {
      opt_list_names = true;
    } else if (std::string("--range") == argv[i]) {
      opt_test_range = true;
      range.first = atoi(argv[++i]);
      range.second = atoi(argv[++i]);
    } else if (std::string("+tag") == argv[i]) {
      include_tags.insert(argv[++i]);
    } else if (std::string("-tag") == argv[i]) {
      exclude_tags.insert(argv[++i]);
    } else {
      regex = argv[i];
    }
  }

  std::vector<std::pair<int, Test>> tests;
  {
    std::cout << "LIST TESTS:" << std::endl;
    int test_id = 0;
    for (auto& test : get_all_tests()) {
      test_id++;
      bool included = false;
      bool excluded = false;
      for (auto& tag : test.tags) {
        if (include_tags.find(tag) != include_tags.end()) {
          included = true;
        }
        if (exclude_tags.find(tag) != exclude_tags.end()) {
          excluded = true;
        }
      }
      if (excluded or (not include_tags.empty() and not included)) {
        continue;
      }
      std::smatch match;
      if (std::regex_match(test.name, match, regex)) {
        if ((!opt_test_range) ||
            ((test_id >= range.first) && (test_id <= range.second))) {
          tests.push_back({test_id, test});
          std::cout << "  [" << test_id << "] " << test.name;
          if (not test.tags.empty()) {
            std::cout << " | Tags: ";
            for (size_t i = 0; i < test.tags.size(); ++i) {
              std::cout << test.tags.at(i);
              if (i + 1 < test.tags.size()) {
                std::cout << ", ";
              }
            }
          }
          std::cout << std::endl;
        }
      }
    }
    std::cout << std::endl;
  }

  if (opt_list_names) {
    return EXIT_SUCCESS;
  }

  {
    std::vector<bool> test_results;
    std::vector<double> test_runtimes;
    std::vector<std::pair<int, Test>> failed;

    size_t ran = 0;
    const auto num_tests = tests.size();
    std::cout << "RUNNING " << num_tests << " TEST(S) ..." << std::endl;
    for (auto& [test_id, test] : tests) {
      ++ran;
      auto start_time = std::chrono::steady_clock::now();

      std::string run_number_str =
          "  (" + std::to_string(ran) + " of " + std::to_string(num_tests) + ")";
      std::string test_id_str = "[" + std::to_string(test_id) + "]";
      std::string test_name_str = "'" + test.name + "'";
      std::string test_header_str = run_number_str + " " + test_id_str;

      std::cout << test_header_str << " TEST RUN: " << test_name_str << " Begins... "
                << std::endl
                << std::flush;

      if (no_catch) {
        test.entry();
        std::cout << test_header_str << " TEST RESULT: " << test_name_str << " Passed."
                  << std::endl;
        test_results.push_back(true);
      } else {
        try {
          test.entry();
          std::cout << test_header_str << " TEST RESULT: " << test_name_str << " Passed."
                    << std::endl;
          test_results.push_back(true);
        } catch (const std::exception& e) {
          failed.push_back({test_id, test});
          std::cerr << test_header_str << " TEST RESULT: " << test_name_str
                    << " failed with '" << e.what() << "'" << std::endl;
          test_results.push_back(false);
        }
      }

      auto diff_time = std::chrono::steady_clock::now() - start_time;
      auto run_time = double(std::chrono::duration <double, std::milli>(diff_time).count()) / 1000.0;
      test_runtimes.push_back(run_time);
    }
    std::cout << std::endl;

    std::cout << "SUMMARY:" << std::endl;
    for (size_t i = 0; i < tests.size(); i++) {
      auto [test_id, test] = tests[i];
      std::string run_number_str =
          "(" + std::to_string(ran) + " of " + std::to_string(num_tests) + ")";
      std::string test_id_str = "[" + std::to_string(test_id) + "]";
      std::string test_name_str = "'" + test.name + "'";
      std::string test_header_str = run_number_str + " " + test_id_str;
      std::string test_result_str = (test_results[i] ? "PASS" : "FAIL");
      std::string test_runtime_str = std::to_string(test_runtimes[i]) + "s";

      std::cout << "  " << test_id_str << " " << test_result_str
                << " | " << test_runtime_str << " | " << test_name_str << std::endl;
    }
    std::cout << std::endl;

    if (not failed.empty()) {
      std::cerr << "TESTS FAILED: " << failed.size() << "/" << num_tests << std::endl;
      for (auto& [test_id, test] : failed) {
        std::string test_id_str = "[" + std::to_string(test_id) + "]";
        std::cerr << "  " << test_id_str << " " << test.name << "\n";
      }
      print_peak_mem();
      return EXIT_FAILURE;
    } else {
      std::cout << "ALL TESTS PASSED." << std::endl;
      print_peak_mem();
      return EXIT_SUCCESS;
    }
  }
}
