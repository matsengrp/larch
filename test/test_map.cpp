#include "test_common.hpp"

#include "larch/contiguous_map.hpp"
#include <map>

template <typename K, typename V>
bool map_eq(const std::pair<K, V>& lhs, const std::pair<K, V>& rhs) {
  return lhs == rhs;
}

template <typename K, typename V>
struct TestMap {
  template <typename Fn, typename CheckFn>
  decltype(auto) test(Fn fn, CheckFn check_fn) {
    finally check_eq{[this] { TestAssert(ranges::equal(map1, map2, map_eq<K, V>)); }};
    return check_fn(fn(map1), fn(map2));
  }

  ContiguousMap<K, V> map1;
  std::map<K, V> map2;
};

static void test_map() {
  auto eq = [](auto&& lhs, auto&& rhs) {
    TestAssert(lhs == rhs);
    return lhs;
  };
  auto inspair_eq = [](auto&& lhs, auto&& rhs) {
    TestAssert(lhs.second == rhs.second);
    TestAssert(lhs.first->first == rhs.first->first);
    TestAssert(lhs.first->second == rhs.first->second);
    return lhs;
  };
  auto find_eq = [](auto&& lhs, auto&& rhs) {
    if (lhs.first == lhs.second.end()) {
      TestAssert(rhs.first == rhs.second.end());
    } else {
      TestAssert(lhs.first->first == rhs.first->first);
      TestAssert(lhs.first->second == rhs.first->second);
    }
    return lhs.first;
  };

  TestMap<int, int> m;

  m.test([](auto&& x) { return x.size(); }, eq);
  m.test([](auto&& x) { return x.empty(); }, eq);

  m.test([](auto&& x) { return x.insert({5, 10}); }, inspair_eq);
  m.test([](auto&& x) { return x.insert({5, 10}); }, inspair_eq);
  m.test([](auto&& x) { return x.insert({6, 12}); }, inspair_eq);
  m.test([](auto&& x) { return x.insert({6, 12}); }, inspair_eq);
  m.test([](auto&& x) { return x.insert({5, 10}); }, inspair_eq);
  m.test([](auto&& x) { return x.insert({5, 15}); }, inspair_eq);
  m.test([](auto&& x) { return x.insert({5, 10}); }, inspair_eq);

  m.test([](auto&& x) { return std::make_pair(x.find(5), std::ref(x)); }, find_eq);
  m.test([](auto&& x) { return std::make_pair(x.find(6), std::ref(x)); }, find_eq);
  m.test([](auto&& x) { return std::make_pair(x.find(7), std::ref(x)); }, find_eq);

  m.test(
      [](auto&& x) {
        x[7] = 14;
        return true;
      },
      eq);
  m.test(
      [](auto&& x) {
        x[6] = 18;
        return true;
      },
      eq);

  m.test([](auto&& x) { return x.insert_or_assign(5, 15); }, inspair_eq);
  m.test([](auto&& x) { return x.insert_or_assign(8, 16); }, inspair_eq);
  m.test([](auto&& x) { return x.insert_or_assign(5, 20); }, inspair_eq);
  m.test([](auto&& x) { return x.insert_or_assign(10, 20); }, inspair_eq);
}

[[maybe_unused]] static const auto test_added = add_test({[] { test_map(); }, "Map"});
