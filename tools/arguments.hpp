#pragma once

#include <string_view>
#include <initializer_list>

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
