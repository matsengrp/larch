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
  std::regex regex{".*"};
  for (int i = 1; i < argc; ++i) {
    if (std::string("nocatch") == argv[i]) {
      no_catch = true;
    } else {
      regex = argv[i];
    }
  }

  std::vector<Test> tests;
  for (auto& test : get_all_tests()) {
    std::smatch match;
    if (std::regex_match(test.name, match, regex)) {
      tests.push_back(test);
    }
  }

  size_t failed = 0, ran = 0;
  const auto num_tests = tests.size();
  std::cout << "Running " << num_tests << " tests" << std::endl;
  for (auto& test : tests) {
    ++ran;
    std::cout << "Running test: " << test.name << " (" << ran << "/" << num_tests
              << ") ..." << std::flush;

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
