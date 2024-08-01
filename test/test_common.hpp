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

#define TestThrow(x)                                                       \
  {                                                                        \
    bool caught = false;                                                   \
    try {                                                                  \
      (x);                                                                 \
    } catch (...) {                                                        \
      caught = true;                                                       \
    }                                                                      \
    if (not caught) {                                                      \
      throw std::runtime_error("TestThrow failed: \"" #x "\" in " __FILE__ \
                               ":" TOSTRING(__LINE__));                    \
    }                                                                      \
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
