#pragma once

#include <string_view>
#include <initializer_list>
#include <iostream>
#include <iomanip>

#include "version.hpp"
#include "larch/common.hpp"

[[noreturn]] inline static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
}

inline static void Version(std::string program_name = "larch-usher") {
  std::cout << program_name << "\n";
  std::cout << "Build version: " << VERSION_NUMBER << "\n";
  std::cout << "Build date: " << GIT_COMMIT_DATE << "\n";
  std::cout << "Build commit: " << GIT_COMMIT_HASH << "\n";
  std::cout << "Build type: " << BUILD_TYPE << "\n";

  exit(EXIT_SUCCESS);
}

//////////////////////////////////////////////////////////////////////////////////////

// Commandline Parsers

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

// Helper to detect vector types
template <typename T>
struct is_vector : std::false_type {};
template <typename T, typename Alloc>
struct is_vector<std::vector<T, Alloc>> : std::true_type {};

template <typename T>
T parse(const std::string& str) {
  std::istringstream iss(str);
  T value;
  iss >> value;
  return value;
}

// If number of params can vary but must be non-zero, use req_num_params = -1.
// If number of params can vary including zero, use req_num_params = -2.
template <bool do_parse = true, typename Params, typename ParamType>
inline void ParseOption(std::string_view name, const Params& params_in,
                        ParamType& params_out, const int req_num_params) {
  auto num_params = std::distance(params_in.begin(), params_in.end());
  if (req_num_params >= 0) {
    if (num_params != req_num_params) {
      std::cerr << "ERROR: Incorrect number of params for `" << name
                << "` option (expected " << req_num_params << ", got " << num_params
                << ").\n";
      std::exit(EXIT_FAILURE);
    }
  } else if (req_num_params == -1) {
    if (num_params == 0) {
      std::cerr << "ERROR: Incorrect number of params for `" << name
                << "` options. Must be non-zero.\n";
      std::exit(EXIT_FAILURE);
    }
  }

  if constexpr (do_parse) {
    // flag types -- flip the initial value
    if constexpr (std::is_same<bool, ParamType>::value and req_num_params == 0) {
      params_out = !params_out;
      return;
    }
    // strings, ints, and real numbers
    else if constexpr (std::is_same<std::string, ParamType>::value or
                       std::is_integral<ParamType>::value or
                       std::is_floating_point<ParamType>::value) {
      params_out = parse<ParamType>(std::string{*params_in.begin()});
      return;
    }
    // vectors
    else if constexpr (is_vector<ParamType>::value) {
      using VectorType = typename ParamType::value_type;
      static_assert(std::is_same<std::string, VectorType>::value or
                        std::is_integral<VectorType>::value or
                        std::is_floating_point<VectorType>::value,
                    "Vector type not supported");
      params_out.clear();
      for (const auto& param : params_in) {
        params_out.push_back(parse<VectorType>(std::string{param}));
      }
    } else {
      static_assert(!std::is_same<ParamType, ParamType>::value,
                    "Param type is not supported.");
    }
  }
}
