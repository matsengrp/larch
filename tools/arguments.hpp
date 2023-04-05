#pragma once

#include <string_view>
#include <initializer_list>

#include "larch/common.hpp"

inline auto GetArguments(int argc, char** argv) {
  return ranges::views::counted(argv, argc) | ranges::views::drop(1) |
         ranges::views::transform([](auto i) { return std::string_view{i}; }) |
         ranges::views::group_by([](auto _, auto i) {
           // should return true if this is not a named argument (e.g. a value)
           // only arguments of length 1 are allowed to have a single -,
           // otherwise, they're interpreted as non-arguments (e.g. -02 is
           // a number, -2 would be interpreted as an argument)
           std::ignore = _;
           if (i.length() < 2){
             return true;
           } else if (i.length() == 2) {
             return i.find('-') != 0;
           } else {
             return i[0] != '-' or i[1] != '-';
           }
         }) |
         ranges::views::transform([](auto i) {
           return std::make_pair(*i.begin(), i | ranges::views::drop(1));
         });
}

using Arguments = decltype(GetArguments(0, nullptr));
