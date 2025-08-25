#pragma once

#include <cstddef>
#include <limits>
#include <vector>
#include <set>
#include <unordered_set>
#include <tuple>
#include <typeinfo>
#include <random>
#include <optional>
#include <variant>
#include <mutex>
#include <shared_mutex>
#include <execution>
#include <thread>
#include <atomic>
#include <iostream>

//////////////////////////////////////////////////////////////////////////////////////

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wstack-usage="
#include <range/v3/all.hpp>
#pragma GCC diagnostic pop

#include "larch/debug.hpp"

static constexpr const size_t NoId = std::numeric_limits<size_t>::max();

template <typename T>
struct remove_cvref {
  using type = std::remove_cv_t<std::remove_reference_t<T>>;
};

template <typename T>
using remove_cvref_t = typename remove_cvref<T>::type;

template <typename T>
struct type_identity {
  using type = T;
};

template <typename Fn>
struct finally {
  finally(Fn&& fn) : fn_{std::forward<Fn>(fn)} {}
  ~finally() { std::invoke(fn_); }

 private:
  Fn fn_;
};

///////////////////////////////////////////////////////

template <typename T>
struct is_variant : std::false_type {};

template <typename... Args>
struct is_variant<std::variant<Args...>> : std::true_type {};

template <typename T>
static constexpr const bool is_variant_v = is_variant<T>::value;

///////////////////////////////////////////////////////

template <typename T>
using value_t_helper = ranges::iter_value_t<ranges::iterator_t<T>>;

template <typename... Views>
struct variant_of_views_helper {
  static_assert(sizeof...(Views) > 0);
  static_assert((ranges::view_<std::remove_reference_t<Views>> and ...));
  static_assert((ranges::common_range<std::remove_reference_t<Views>> and ...));

  using value_t = value_t_helper<
      std::tuple_element_t<0, std::tuple<std::remove_reference_t<Views>...>>>;

  constexpr static bool has_nested = ranges::range<value_t>;

  static_assert(has_nested or
                (std::is_convertible_v<value_t_helper<std::remove_reference_t<Views>>,
                                       value_t> and
                 ...));

  using views_variant_t = std::variant<remove_cvref_t<Views>...>;

  using iter_variant_t = std::variant<ranges::iterator_t<const Views&>...>;

  static constexpr auto index_seq = std::index_sequence_for<Views...>{};

  template <size_t... I>
  static bool iters_equal(const iter_variant_t& lhs, const iter_variant_t& rhs,
                          std::index_sequence<I...>) {
    if (lhs.index() != rhs.index()) {
      return false;
    }
    return ((I == lhs.index() and std::get<I>(lhs) == std::get<I>(rhs)) or ...);
  }
};

template <typename...>
class variant_of_views;

template <typename... Views>
struct variant_of_views_nested_helper {
  static_assert(variant_of_views_helper<Views...>::has_nested);
  using variant_of_views_t =
      variant_of_views<value_t_helper<std::remove_reference_t<Views>>...>;
};

template <typename... Views>
class variant_of_views_iterator {
  using helper = variant_of_views_helper<Views...>;

 public:
  using iterator_category = std::input_iterator_tag;
  using value_type = typename helper::value_t;
  using difference_type = std::ptrdiff_t;

  variant_of_views_iterator() = default;
  variant_of_views_iterator(const variant_of_views_iterator&) = default;
  variant_of_views_iterator& operator=(const variant_of_views_iterator&) = default;

  variant_of_views_iterator(const typename helper::iter_variant_t& iter)
      : iter_{iter} {}

  bool operator==(const variant_of_views_iterator& other) const {
    return helper::iters_equal(iter_, other.iter_, helper::index_seq);
  }

  bool operator!=(const variant_of_views_iterator& other) const {
    return not(*this == other);
  }

  decltype(auto) operator*() const {
    if constexpr (helper::has_nested) {
      return std::visit(
          [](auto& x) {
            using nested =
                typename variant_of_views_nested_helper<Views...>::variant_of_views_t;
            return nested{*x};
          },
          iter_);
    } else {
      return std::visit([](auto& x) { return *x; }, iter_);
    }
  }

  variant_of_views_iterator& operator++() {
    std::visit([](auto& x) { ++x; }, iter_);
    return *this;
  }

  void operator++(int) { ++*this; }

 private:
  template <
      typename Sentinel,
      typename = std::enable_if_t<
          (ranges::sentinel_for<Sentinel, ranges::iterator_t<remove_cvref_t<Views>>> or
           ...)>>
  friend bool operator==(const variant_of_views_iterator<Views...>& iter,
                         Sentinel&& sent) {
    return iter.iter_ == sent;
  }

  template <
      typename Sentinel,
      typename = std::enable_if_t<
          (ranges::sentinel_for<Sentinel, ranges::iterator_t<remove_cvref_t<Views>>> or
           ...)>>
  friend bool operator==(Sentinel&& sent,
                         const variant_of_views_iterator<Views...>& iter) {
    return iter.iter_ == sent;
  }

  typename helper::iter_variant_t iter_;
};

template <typename... Views>
class variant_of_views {
  using helper = variant_of_views_helper<Views...>;
  using iterator = variant_of_views_iterator<Views...>;

 public:
  variant_of_views() = default;
  variant_of_views(const variant_of_views<Views...>&) = default;
  variant_of_views(variant_of_views<Views...>&&) = default;
  variant_of_views(variant_of_views<Views...>&) = default;
  variant_of_views& operator=(const variant_of_views<Views...>&) = default;
  variant_of_views& operator=(variant_of_views<Views...>&&) = default;

  template <typename T>
  variant_of_views(T&& view) : view_{std::forward<T>(view)} {}

  iterator begin() const {
    return std::visit([](auto& x) { return iterator{x.begin()}; }, view_);
  }

  iterator end() const {
    return std::visit([](auto& x) { return iterator{x.end()}; }, view_);
  }

  size_t size() const {
    return std::visit([](auto& x) { return static_cast<size_t>(x.size()); }, view_);
  }

  bool empty() const {
    return std::visit([](auto& x) { return x.empty(); }, view_);
  }

 private:
  typename helper::views_variant_t view_;
};

template <typename... Views>
constexpr bool ranges::enable_view<variant_of_views<Views...>> = true;

///////////////////////////////////////////////////////

enum class IdContinuity { Dense, Sparse };

enum class Ordering { Ordered, Unordered };

struct NodeId;
struct EdgeId;
struct CladeIdx;

template <typename T, typename Id>
[[nodiscard]] static T& GetOrInsert(std::vector<T>& data, Id id) {
  size_t idx{};
  if constexpr (std::is_same_v<Id, size_t>) {
    idx = id;
  } else {
    idx = id.value;
  }

  if (idx >= data.size()) {
    data.resize(idx + 1);
  }
  return data[idx];
}

template <typename T>
void MoveElements(std::vector<T>&& source, std::vector<T>& destination) {
  destination.resize(0);
  destination.insert(destination.end(), std::make_move_iterator(source.begin()),
                     std::make_move_iterator(source.end()));
  source.resize(0);
}

namespace Transform {
template <typename T>
inline auto To() {
  return ranges::views::transform([&](auto&& i) { return T{i}; });
}
}  // namespace Transform

inline constexpr const auto HashCombine = [](size_t lhs, size_t rhs) noexcept {
  constexpr const size_t GoldenRatioFractional = 0x9e3779b97f4a7c15;
  constexpr const size_t LeftShift64 = 12;
  constexpr const size_t RightShift64 = 4;
  lhs ^= rhs + GoldenRatioFractional + (lhs << LeftShift64) + (lhs >> RightShift64);
  return lhs;
};

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifdef KEEP_ASSERTS
#ifdef USE_CPPTRACE
#define Assert(x)                                                       \
  {                                                                     \
    if (not(x)) {                                                       \
      std::cerr << "Assert failed: \"" #x "\" in " __FILE__             \
                   ":" TOSTRING(__LINE__) "\n";                         \
      DebugItem::print_current_trace();                                 \
      throw std::runtime_error("Assert failed: \"" #x "\" in " __FILE__ \
                               ":" TOSTRING(__LINE__));                 \
    }                                                                   \
  }
#else
#define Assert(x)                                                       \
  {                                                                     \
    if (not(x)) {                                                       \
      std::cerr << "Assert failed: \"" #x "\" in " __FILE__             \
                   ":" TOSTRING(__LINE__) "\n";                         \
      throw std::runtime_error("Assert failed: \"" #x "\" in " __FILE__ \
                               ":" TOSTRING(__LINE__));                 \
    }                                                                   \
  }
#endif
#else
#define Assert(x) \
  {}
#endif

[[noreturn]] inline void Fail(const char* msg) {
#ifdef KEEP_ASSERTS
  std::cerr << msg << "\n";
#endif
#ifdef USE_CPPTRACE
  DebugItem::print_current_trace();
#endif
  throw std::runtime_error(msg);
}

template <typename ReturnType = void>
[[noreturn]] ReturnType* Unreachable() {
  Fail("Unreachable");
}

#define MOVE_ONLY(x)                    \
  x(x&&) noexcept = default;            \
  x(const x&) = delete;                 \
  x& operator=(x&&) noexcept = default; \
  x& operator=(const x&) = delete

#define MOVE_ONLY_VIRT_DTOR(x)          \
  x(x&&) noexcept = default;            \
  x(const x&) = delete;                 \
  x& operator=(x&&) noexcept = default; \
  x& operator=(const x&) = delete;      \
  virtual ~x() = default

#define NO_COPY(x)      \
  x(const x&) = delete; \
  x& operator=(const x&) = delete

#define MOVE_ONLY_DEF_CTOR(x) \
  MOVE_ONLY(x);               \
  x() = default

//////////////////////////////////////////////////////////////////////////////////////

template <typename, template <typename...> typename>
struct is_specialization : std::false_type {};

template <template <typename...> typename Template, typename... Args>
struct is_specialization<Template<Args...>, Template> : std::true_type {};

template <typename L, template <typename...> typename R>
constexpr bool is_specialization_v = is_specialization<L, R>::value;

//////////////////////////////////////////////////////////////////////////////////////

template <typename, typename, template <typename, typename> typename>
struct tuple_contains {};

template <typename... Types, typename Type, template <typename, typename> typename Comp>
struct tuple_contains<std::tuple<Types...>, Type, Comp>
    : std::bool_constant<(Comp<Types, Type>::value || ...)> {};
template <typename... Types, typename Type, template <typename, typename> typename Comp>
struct tuple_contains<const std::tuple<Types...>, Type, Comp>
    : std::bool_constant<(Comp<Types, Type>::value || ...)> {};

template <typename Tuple, typename Type, template <typename, typename> typename Comp>
inline constexpr bool tuple_contains_v = tuple_contains<Tuple, Type, Comp>::value;

//////////////////////////////////////////////////////////////////////////////////////

template <typename T, template <typename, typename> typename Comp, typename Tuple,
          size_t I = 0>
constexpr auto& tuple_get(Tuple& x) {
  if constexpr (I < std::tuple_size_v<Tuple>) {
    using E = std::remove_const_t<std::tuple_element_t<I, Tuple>>;
    if constexpr (Comp<T, E>::value) {
      return std::get<I>(x);
    } else {
      if constexpr (I + 1 < std::tuple_size_v<Tuple>) {
        return tuple_get<T, Comp, Tuple, I + 1>(x);
      } else {
        static_assert(sizeof(T) == NoId, "Not found");
      }
    }
  } else {
    static_assert(sizeof(T) == NoId, "Not found");
  }
}

//////////////////////////////////////////////////////////////////////////////////////

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wc++20-extensions"
template <typename T>
inline const auto tuple_to_string_impl =
    []<std::size_t... I>(std::index_sequence<I...>) {
      std::string result = "std::tuple<";
      result += (... + (std::string{typeid(std::tuple_element_t<I, T>).name()} +
                        std::string{", "}));
      result.pop_back();
      result.pop_back();
      result += ">";
      return result;
    };

template <typename T>
inline std::string tuple_to_string() {
  return tuple_to_string_impl<T>(std::make_index_sequence<std::tuple_size_v<T>>());
}
#pragma GCC diagnostic pop

//////////////////////////////////////////////////////////////////////////////////////

// Generic containers to stream.

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& vector) {
  os << "[ ";
  for (const auto& element : vector) {
    os << element << ", ";
  }
  os << "]";
  return os;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const std::set<T>& set) {
  os << "{ ";
  for (const auto& element : set) {
    os << element << ", ";
  }
  os << "}";
  return os;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const std::unordered_set<T>& set) {
  os << "{ ";
  for (const auto& element : set) {
    os << element << ", ";
  }
  os << "}";
  return os;
}

template <typename K, typename V>
inline std::ostream& operator<<(std::ostream& os, const std::map<K, V>& map) {
  os << "{ ";
  for (const auto& [key, value] : map) {
    os << key << ": " << value << ", ";
  }
  os << "}";
  return os;
}

template <typename K, typename V>
inline std::ostream& operator<<(std::ostream& os, const std::unordered_map<K, V>& map) {
  os << "{ ";
  for (const auto& [key, value] : map) {
    os << key << ": " << value << ", ";
  }
  os << "}";
  return os;
}

template <typename T1, typename T2>
inline std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2>& tup) {
  os << "[ " << tup.first << ", " << tup.second << " ], ";
  return os;
}

template <typename T1, typename T2, typename T3>
inline std::ostream& operator<<(std::ostream& os, const std::tuple<T1, T2, T3>& tup) {
  os << "[ " << std::get<0>(tup) << ", " << std::get<1>(tup) << ", " << std::get<2>(tup)
     << " ], ";
  return os;
}

//////////////////////////////////////////////////////////////////////////////////////

// Common utilities

class RandomNumberGenerator {
 public:
  RandomNumberGenerator(std::optional<uint32_t> user_seed) {
    std::random_device rd;
    seed_ = user_seed.value_or(rd());
    generator_ = std::mt19937(seed_);
  }

  template <typename T>
  T Generate() {
    if constexpr (std::is_integral<T>::value) {
      std::uniform_int_distribution<T> dist;
      return dist(generator_);
    } else if constexpr (std::is_floating_point<T>::value) {
      std::uniform_real_distribution<T> dist(0.0, 1.0);
      return dist(generator_);
    } else {
      static_assert(!sizeof(T),
                    "Only integral and floating point types are supported.");
    }
  }

  uint32_t GenerateSeed() { return Generate<uint32_t>(); }

  std::mt19937 generator_;
  uint32_t seed_;
};

inline std::vector<std::string> SplitString(const std::string& str, char delim) {
  std::vector<std::string> tokens;
  size_t start = 0;
  size_t end = str.find(delim);

  while (end != std::string::npos) {
    tokens.push_back(str.substr(start, end - start));
    start = end + 1;
    end = str.find(delim, start);
  }
  tokens.push_back(str.substr(start));

  return tokens;
}

inline int ParseNumber(std::string_view str) {
  int result{};
  std::istringstream stream{std::string{str}};
  stream >> result;
  if (stream.fail()) {
    throw std::runtime_error("Invalid number");
  }
  return result;
}

//////////////////////////////////////////////////////////////////////////////////////
