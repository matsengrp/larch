#pragma once

#include <string_view>
#include <initializer_list>

#include <range/v3/view/counted.hpp>
#include <range/v3/view/group_by.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/drop.hpp>

inline auto GetArguments(int argc, char** argv) {
  return ranges::views::counted(argv, argc) | ranges::views::drop(1) |
         ranges::views::transform([](auto i) { return std::string_view{i}; }) |
         ranges::views::group_by([](auto _, auto i) { return i.find('-') != 0; }) |
         ranges::views::transform([](auto i) {
           return std::make_pair(*i.begin(), i | ranges::views::drop(1));
         });
}

using Arguments = decltype(GetArguments(0, nullptr));
