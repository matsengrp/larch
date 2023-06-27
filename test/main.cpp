#include <iostream>
#include <vector>
#include <regex>

#include <mpi.h>

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
  int ignored{};
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &ignored);
  bool no_catch = false;
  bool opt_list_names = false;
  bool opt_test_range = false;
  std::pair<size_t, size_t> range;
  std::regex regex{".*"};
  for (int i = 1; i < argc; ++i) {
    if (std::string("nocatch") == argv[i]) {
      no_catch = true;
    } else if (std::string("--list") == argv[i]) {
      opt_list_names = true;
    } else if (std::string("--range") == argv[i]) {
      opt_test_range = true;
      range.first = atoi(argv[++i]);
      range.second = atoi(argv[++i]);
    } else {
      regex = argv[i];
    }
  }

  std::cout << "Run tests:" << std::endl;
  std::vector<Test> tests;
  size_t test_counter = 1;
  for (auto& test : get_all_tests()) {
    std::smatch match;
    if (std::regex_match(test.name, match, regex)) {
      if ((!opt_test_range) ||
          ((test_counter >= range.first) && (test_counter <= range.second))) {
        tests.push_back(test);
        std::cout << "[" << test_counter << "] " << test.name << std::endl;
      }
    }
    test_counter++;
  }

  if (opt_list_names) {
    return EXIT_SUCCESS;
  }

  size_t failed = 0, ran = 0;
  const auto num_tests = tests.size();
  std::cout << "Running " << num_tests << " tests" << std::endl;
  for (auto& test : tests) {
    ++ran;
    std::cout << "Running test: " << test.name << " (" << ran << "/" << num_tests
              << ") ..." << std::endl
              << std::flush;

    if (no_catch) {
      test.entry();
      std::cout << " passed." << std::endl;

    } else {
      try {
        test.entry();
        std::cout << " passed." << std::endl;
      } catch (const std::exception& e) {
        ++failed;
        std::cerr << std::endl
                  << "Test '" << test.name << "' failed with '" << e.what() << "'"
                  << std::endl;
      }
    }
  }
  if (failed > 0) {
    std::cerr << "Failed tests: " << failed << "/" << num_tests << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "All tests passed" << std::endl;
  return EXIT_SUCCESS;
}
