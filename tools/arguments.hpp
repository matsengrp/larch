#pragma once

#include <string_view>
#include <initializer_list>

#include "common.hpp"

inline auto GetArguments(int argc, char** argv) {
  return ranges::views::counted(argv, argc) | ranges::views::drop(1) |
         ranges::views::transform([](auto i) { return std::string_view{i}; }) |
         ranges::views::group_by([](auto _, auto i) {
           std::ignore = _;
           return i.find('-') != 0;
         }) |
         ranges::views::transform([](auto i) {
           return std::make_pair(*i.begin(), i | ranges::views::drop(1));
         });
}

using Arguments = decltype(GetArguments(0, nullptr));
