#pragma once

#include <string_view>
#include <initializer_list>
#include <iostream>
#include <iomanip>

#include "larch/common.hpp"

[[noreturn]] inline static void Fail() {
  std::cerr << "Run with -h or --help to see usage.\n";

  std::exit(EXIT_FAILURE);
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
                        ParamType& params_out, const int req_num_params = 1) {
  auto num_params = std::distance(params_in.begin(), params_in.end());
  if (req_num_params >= 0) {
    if (num_params != req_num_params) {
      std::cerr << "Incorrect number of params for `" << name << "` option (expected "
                << req_num_params << ", got " << num_params << ").\n";
      std::exit(EXIT_FAILURE);
    }
  } else if (req_num_params == -1) {
    if (num_params == 0) {
      std::cerr << "Incorrect number of params for `" << name
                << "` options.  Must be nonzero.\n";
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
      for (const auto& param : params_in) {
        params_out = parse<ParamType>(std::string{param});
        return;
      }
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

//////////////////////////////////////////////////////////////////////////////////////

// DAG Info (Parsimony scores, RF Distance scores, etc.)

// enum class DAGScoreType {
//   MinParsimony,
//   MaxParsimony,
//   MinSumRFDistance,
//   MaxSumRFDistance
// };

// template <typename DAG, typename WeightOps>
// auto DAG_GetAllScores(DAG dag, WeightOps weight_ops = {}) {
//   SumRFDistance sum_rf_dist_weight_ops{merge, merge};
//   SubtreeWeight<WeightAccumulator<SumRFDistance>, DAG> sum_rf_dist_counter{
//       merge.GetResult()};
//   auto sum_rf_dist_counts = sum_rf_dist_counter.ComputeWeightBelow(
//       merge.GetResult().GetRoot(), WeightAccumulator{sum_rf_dist_weight_ops});
//   using Weights = decltype(sum_rf_dist_counts.GetWeights());
//   Weights all_sum_rf_dist_data{sum_rf_dist_counts.GetWeights()};

//   auto shift_sum = sum_rf_dist_weight_ops.GetOps().GetShiftSum();
//   for (auto& score_count : all_sum_rf_dist_data) {
//     score_count.first += shift_sum;
//   }

//   auto min_sum_rf_dist_data =
//       *std::min_element(sum_rf_dist_counts.GetWeights().begin(),
//                         sum_rf_dist_counts.GetWeights().end(), scorecount_compare);
//   auto max_sum_rf_dist_data =
//       *std::max_element(sum_rf_dist_counts.GetWeights().begin(),
//                         sum_rf_dist_counts.GetWeights().end(), scorecount_compare);

//   return {all_sum_rf_dist_data, min_sum_rf_dist_data, max_sum_rf_dist_data};
// }

// template <typename DAG>
// auto DAG_GetScore(DAG dag, WeightOps weight_ops = {}) {
//   SubtreeWeight<SumRFDistance, DAG> min_sum_rf_dist{dag.GetResult()};
//   SumRFDistance min_rf_weight_ops{dag, dag};
//   SubtreeWeight<MaxSumRFDistance, DAG> max_sum_rf_dist{dag.GetResult()};
//   MaxSumRFDistance max_rf_weight_ops{merge, merge};
//   auto shiftsum = min_rf_weight_ops.GetOps().GetShiftSum();

//   auto min_rf_distance =
//   min_sum_rf_dist.ComputeWeightBelow(merge.GetResult().GetRoot(),
//                                                             min_rf_weight_ops) +
//                          shiftsum;
//   auto min_rf_count =
//       min_sum_rf_dist.MinWeightCount(merge.GetResult().GetRoot(),
//       min_rf_weight_ops);
//   auto max_rf_distance =
//   max_sum_rf_dist.ComputeWeightBelow(merge.GetResult().GetRoot(),
//                                                             max_rf_weight_ops) +
//                          shiftsum;
//   auto max_rf_count =
//       max_sum_rf_dist.MinWeightCount(merge.GetResult().GetRoot(),
//       max_rf_weight_ops);

//   return {min_rf_distance, min_rf_count, max_rf_distance, max_rf_count};
// }
