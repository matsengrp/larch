#pragma once

#include <string_view>
#include <initializer_list>
#include <iostream>
#include <iomanip>

#include "larch/common.hpp"

inline auto GetArguments(int argc, char** argv) {
  return ranges::views::counted(argv, argc) | ranges::views::drop(1) |
         ranges::views::transform([](auto i) { return std::string_view{i}; }) |
         ranges::views::chunk_by([](auto lhs, auto rhs) {
           return lhs.find('-') == 0 and rhs.find('-') != 0;
         }) |
         ranges::views::transform([](auto i) {
           return std::make_pair(*i.begin(), i | ranges::views::drop(1));
         });
}

using Arguments = decltype(GetArguments(0, nullptr));

inline auto FormatUsage(
    const std::string& program_desc, const std::vector<std::string>& usage_examples,
    const std::vector<std::pair<std::string, std::string>>& flag_desc_pairs) {
  std::stringstream os;

  if (!program_desc.empty()) {
    os << program_desc << "\n\n";
  }

  os << "Usage:\n";
  for (const auto& usage_example : usage_examples) {
    os << "  " << usage_example << "\n";
  }
  os << "\n";

  size_t max_flag_width = 0;
  for (const auto& [flag, desc] : flag_desc_pairs) {
    std::ignore = desc;
    max_flag_width = std::max(max_flag_width, flag.size());
  }

  os << "Options:\n";
  for (const auto& [flag, desc] : flag_desc_pairs) {
    const auto lines = SplitString(desc, '\n');
    for (size_t i = 0; i < lines.size(); i++) {
      os << "  " << std::left << std::setw(static_cast<int>(max_flag_width) + 2)
         << ((i == 0) ? flag : "");
      os << ((i == 0) ? "" : "  ") << lines[i] << "\n";
    }
  }

  return os.str();
}
