#pragma once

#include <string_view>
#include <initializer_list>

#include <range/v3/view/counted.hpp>
#include <range/v3/view/group_by.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/drop.hpp>

inline auto GetArguments(int argc, char** argv) {
  return ranges::view::counted(argv, argc) | ranges::view::drop(1) |
         ranges::view::transform([](auto i) { return std::string_view{i}; }) |
         ranges::view::group_by([](auto _, auto i) { return i.find('-') != 0; }) |
         ranges::view::transform([](auto i) {
           return std::make_pair(*i.begin(), i | ranges::view::drop(1));
         });
}

using Arguments = decltype(GetArguments(0, nullptr));
