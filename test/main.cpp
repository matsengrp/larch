#include <iostream>
#include <set>
#include <regex>

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
  std::pair<size_t, size_t> range;
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
      range.first = static_cast<size_t>(atoi(argv[++i]));
      range.second = static_cast<size_t>(atoi(argv[++i]));
    } else if (std::string("+tag") == argv[i]) {
      include_tags.insert(argv[++i]);
    } else if (std::string("-tag") == argv[i]) {
      exclude_tags.insert(argv[++i]);
    } else {
      regex = argv[i];
    }
  }

  std::cout << "LIST TESTS:" << std::endl;
  std::vector<Test> tests;
  size_t test_counter = 1;
  for (auto& test : get_all_tests()) {
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
          ((test_counter >= range.first) && (test_counter <= range.second))) {
        tests.push_back(test);
        std::cout << "  [" << test_counter << "] " << test.name;
        if (not test.tags.empty()) {
          std::cout << " Tags: ";
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
    test_counter++;
  }

  if (opt_list_names) {
    return EXIT_SUCCESS;
  }

  std::vector<Test> failed;
  size_t ran = 0;
  const auto num_tests = tests.size();
  std::cout << "RUNNING " << num_tests << " TEST(S) ..." << std::endl;
  for (auto& test : tests) {
    ++ran;
    std::string run_number =
        "  (" + std::to_string(ran) + "/" + std::to_string(num_tests) + ")";

    std::cout << run_number << " TEST RUN: Begin '" << test.name << "' ... "
              << std::endl
              << std::flush;

    if (no_catch) {
      test.entry();
      std::cout << run_number << " TEST RESULT: Passed." << std::endl;
    } else {
      try {
        test.entry();
        std::cout << run_number << " TEST RESULT: Passed." << std::endl;
      } catch (const std::exception& e) {
        failed.push_back(test);
        std::cerr << run_number << " TEST RESULT: Test '" << test.name
                  << "' failed with '" << e.what() << "'" << std::endl;
      }
    }
  }
  if (not failed.empty()) {
    std::cerr << "TESTS FAILED: " << failed.size() << "/" << num_tests << std::endl;
    for (auto& test : failed) {
      std::cerr << "  " << test.name << "\n";
    }
    print_peak_mem();
    return EXIT_FAILURE;
  }

  std::cout << "ALL TESTS PASSED." << std::endl;
  print_peak_mem();
  return EXIT_SUCCESS;
}
