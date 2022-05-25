#include <cstdlib>
#include <iostream>

#include "arguments.hpp"
#include "merge.hpp"
#include "history_dag_loader.hpp"

int main(int argc, char** argv) {
  Arguments args = GetArguments(argc, argv);

  std::string_view lhs_filename;
  std::string_view rhs_filename;

  for (auto [name, params] : args) {
    std::cout << "Name: " << name << "\n";
    for (auto par : params) {
        std::cout << "  Param: " << par << "\n";
    }
  }

  return EXIT_SUCCESS;
}
