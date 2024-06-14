#include "larch/subtree/subtree_weight.hpp"
#include "larch/subtree/parsimony_score.hpp"
#include "larch/subtree/weight_accumulator.hpp"

#include <iostream>
#include <string_view>
#include <string>

#include "test_common.hpp"

#include "larch/dag_loader.hpp"

using Weight = typename WeightAccumulator<ParsimonyScore>::Weight;
using Counter = WeightCounter<ParsimonyScore>;

Counter make_counter(const std::vector<ParsimonyScore::Weight> inweights) {
  return Counter(inweights, {});
}

// template <typename T>
// std::string to_string(T& inobject) {
//   std::ostringstream result;
//   result << inobject;
//   return result.str();
// }

// template std::string to_string(Count&);
// template std::string to_string(size_t&);

static void test_counter_add(Counter lhs, Counter rhs, Counter expected_counter) {
  Counter result = lhs + rhs;

  std::cout << "lhs = " << lhs << "\n";
  std::cout << "rhs = " << rhs << "\n";
  std::cout << "(lhs + rhs) = " << result << "\n";

  TestAssert(result == expected_counter);
}

static void test_counter_multiply(Counter lhs, Counter rhs, Counter expected_counter) {
  Counter result = lhs * rhs;

  std::cout << "lhs = " << lhs << "\n";
  std::cout << "rhs = " << rhs << "\n";
  std::cout << "(lhs * rhs) = " << result << "\n";

  TestAssert(result == expected_counter);
}

[[maybe_unused]] static const auto test_added0 = add_test(
    {[] {
       test_counter_multiply(make_counter({}), make_counter({2, 2, 3, 3, 3, 4}),
                             make_counter({}));
     },
     "Counter: multiply1"});

[[maybe_unused]] static const auto test_added1 = add_test(
    {[] {
       test_counter_multiply(make_counter({0}), make_counter({2, 2, 3, 3, 3, 4}),
                             make_counter({2, 2, 3, 3, 3, 4}));
     },
     "Counter: multiply2"});

[[maybe_unused]] static const auto test_added2 = add_test(
    {[] {
       test_counter_multiply(
           make_counter({2, 2, 3}), make_counter({2, 2, 3, 3, 3, 4}),
           make_counter({4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7}));
     },
     "Counter: multiply3"});

[[maybe_unused]] static const auto test_added3 = add_test(
    {[] {
       test_counter_multiply(make_counter({0}), make_counter({2}), make_counter({2}));
     },
     "Counter: multiply4"});

[[maybe_unused]] static const auto test_added4 = add_test(
    {[] { test_counter_add(make_counter({}), make_counter({2}), make_counter({2})); },
     "Counter: add1"});

[[maybe_unused]] static const auto test_added5 = add_test(
    {[] {
       test_counter_add(make_counter({0, 1, 2, 2, 3}), make_counter({2, 2, 2, 3, 4}),
                        make_counter({0, 1, 2, 2, 3, 2, 2, 2, 3, 4}));
     },
     "Counter: add2"});
