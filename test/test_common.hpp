#pragma once

#include <functional>
#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>

#include "larch/dag_loader.hpp"
#include "larch/spr/spr_view.hpp"

struct Test {
  Test(std::function<void()> e, std::string n, std::vector<std::string> ts)
      : entry{e}, name{n}, tags{ts} {}
  Test(std::function<void()> e, std::string n) : entry{e}, name{n} {}
  std::function<void()> entry;
  std::string name;
  std::vector<std::string> tags;
};

bool add_test(const Test& test) noexcept;

#define TestAssert(x)                                                       \
  {                                                                         \
    if (not(x)) {                                                           \
      throw std::runtime_error("TestAssert failed: \"" #x "\" in " __FILE__ \
                               ":" TOSTRING(__LINE__));                     \
    }                                                                       \
  }

inline void print_peak_mem() {
  std::ifstream str{"/proc/self/status"};
  std::string line;
  while (std::getline(str, line)) {
    if (line.find("VmPeak:") == 0) {
      std::cout << line << "\n";
      break;
    }
  }
}

inline const std::string test_output_folder = "data/_ignore/";

inline std::pair<std::string, int> run_larch_usher(
    std::string_view input_dag_path, std::string_view output_dag_path,
    std::optional<std::string_view> refseq_path = std::nullopt,
    std::optional<int> iter = std::nullopt,
    std::optional<std::string_view> other_options = std::nullopt,
    bool do_print_stdout = true, bool do_print_stderr = false,
    bool do_print_summary = true) {
  std::string log_path = test_output_folder + "/optimization_log";

  std::stringstream ss;
  // ss << "/usr/bin/time ";
  ss << "./larch-usher ";
  ss << "-i " << input_dag_path << " ";
  ss << "-o " << output_dag_path << " ";
  ss << "-l " << log_path << " ";
  if (refseq_path.has_value()) {
    ss << "-r " << refseq_path.value() << " ";
  }
  if (iter.has_value()) {
    ss << "-c " << iter.value() << " ";
  }
  if (other_options.has_value()) {
    ss << other_options.value() << " ";
  }
  if (!do_print_stdout) {
    ss << "> /dev/null ";
    if (!do_print_stderr) {
      ss << "2>&1 ";
    }
  }

  std::string command = ss.str();
  if (do_print_summary) {
    std::cout << ">COMMAND_EXECUTE: \"" << command << "\"" << std::endl;
  }
  auto result = std::system(command.c_str());
  if (do_print_summary) {
    std::cout << ">COMMAND_RESULT: " << result << " "
              << (result ? "FAILURE" : "SUCCESS") << std::endl;
  }
  return std::make_pair(command, result);
}
