#pragma once

#ifdef USE_CPPTRACE

#include <cpptrace/cpptrace.hpp>
#include <cpptrace/formatting.hpp>

enum class DebugState { Undef, Constructed, MovedFrom, Destructed };

enum class DebugConstructType { Undef, Default, Move, Copy };
enum class DebugAssignType { Undef, Move, Copy };
enum class DebugMoveFromType { Undef, Construct, Assign };

struct Debug;

struct DebugItem {
  DebugState state = DebugState::Undef;
  DebugConstructType construct_type = DebugConstructType::Undef;
  cpptrace::stacktrace construct_trace = {};
  DebugAssignType last_assign_to_type = DebugAssignType::Undef;
  cpptrace::stacktrace last_assign_to_trace = {};
  DebugMoveFromType move_from_type = DebugMoveFromType::Undef;
  cpptrace::stacktrace move_from_trace = {};
  cpptrace::stacktrace destruct_trace = {};

  void constructor() {
    switch (state) {
      case DebugState::Undef:
      case DebugState::Destructed:
        *this = {};
        construct_type = DebugConstructType::Default;
        construct_trace = generate_trace();
        state = DebugState::Constructed;
        break;
      default:
        fail("Double construct");
    }
  }

  void destructor() {
    switch (state) {
      case DebugState::Undef:
        fail("Early destruct");
        break;
      case DebugState::Destructed:
        fail("Double destruct");
        break;
      default:
        destruct_trace = generate_trace();
        state = DebugState::Destructed;
        break;
    }
  }

  void move_constructor_from() {
    switch (state) {
      case DebugState::Undef:
        fail("Early move");
        break;
      case DebugState::MovedFrom:
        fail("Double move");
        break;
      case DebugState::Destructed:
        fail("Move from destructed");
        break;
      default:
        move_from_type = DebugMoveFromType::Construct;
        move_from_trace = generate_trace();
        state = DebugState::MovedFrom;
        break;
    }
  }

  void move_constructor_to() {
    switch (state) {
      case DebugState::Undef:
      case DebugState::Destructed:
        *this = {};
        construct_type = DebugConstructType::Move;
        construct_trace = generate_trace();
        state = DebugState::Constructed;
        break;
      default:
        fail("Double construct");
    }
  }

  void move_assignment_from() {
    switch (state) {
      case DebugState::Undef:
        fail("Early move");
        break;
      case DebugState::MovedFrom:
        fail("Double move");
        break;
      case DebugState::Destructed:
        fail("Move from destructed");
        break;
      default:
        move_from_type = DebugMoveFromType::Assign;
        move_from_trace = generate_trace();
        state = DebugState::MovedFrom;
        break;
    }
  }

  void move_assignment_to() {
    switch (state) {
      case DebugState::Undef:
        fail("Early assign");
        break;
      default:
        *this = {};
        last_assign_to_type = DebugAssignType::Move;
        last_assign_to_trace = generate_trace();
        state = DebugState::Constructed;
        break;
    }
  }

  void copy_constructor_from() {
    switch (state) {
      case DebugState::Undef:
        fail("Early copy");
        break;
      case DebugState::MovedFrom:
        fail("Moved from");
        break;
      case DebugState::Destructed:
        fail("Copy from destructed");
        break;
      default:
        break;
    }
  }

  void copy_constructor_to() {
    switch (state) {
      case DebugState::Undef:
      case DebugState::Destructed:
        *this = {};
        construct_type = DebugConstructType::Copy;
        construct_trace = generate_trace();
        state = DebugState::Constructed;
        break;
      default:
        fail("Double construct");
    }
  }

  void copy_assignment_from() {
    switch (state) {
      case DebugState::Undef:
        fail("Early copy");
        break;
      case DebugState::MovedFrom:
        fail("Moved from");
        break;
      case DebugState::Destructed:
        fail("Copy from destructed");
        break;
      default:
        break;
    }
  }

  void copy_assignment_to() {
    switch (state) {
      case DebugState::Undef:
        fail("Early assign");
        break;
      default:
        *this = {};
        last_assign_to_type = DebugAssignType::Copy;
        last_assign_to_trace = generate_trace();
        state = DebugState::Constructed;
        break;
    }
  }

  void use() const {
    if (state != DebugState::Constructed) {
      fail("Using not constructed");
    }
  }

  static void print_current_trace() {
    cpptrace::formatter formatter = cpptrace::get_default_formatter();
    formatter.paths(cpptrace::formatter::path_mode::basename);
    print(generate_trace(), formatter);
  }

 private:
  [[noreturn]] void fail(const char* msg) const {
    cpptrace::formatter formatter = cpptrace::get_default_formatter();
    formatter.paths(cpptrace::formatter::path_mode::basename);
    std::cout << "\n\nLarchDebug: " << msg << "\n";
    print(generate_trace(), formatter);
    std::cout << "\nConstructed at trace:\n";
    print(construct_trace, formatter);
    std::cout << "\nMoved from trace:\n";
    print(move_from_trace, formatter);
    std::cout << "\nAssign to trace:\n";
    print(last_assign_to_trace, formatter);
    std::cout << "\nDestructed trace:\n";
    print(destruct_trace, formatter);
    throw std::runtime_error(msg);
  }

  static void print(const cpptrace::stacktrace& trace,
                    const cpptrace::formatter& formatter) {
    size_t frame_no = 0;
    for (auto& i : trace.frames) {
      if ((i.filename.find("include/larch") != std::string::npos or
           i.filename.find("test/test_") != std::string::npos) and
          i.symbol.find("Debug::") != 0 and i.symbol.find("DebugItem::") != 0) {
        std::cout << "#" << frame_no << "  ";
        formatter.print(std::cout, i);
        std::cout << "\n" << cpptrace::get_snippet(i.filename, i.line.value(), 5, true);
        std::cout << "\n\n";
      }
      ++frame_no;
    }
  }

  static cpptrace::stacktrace generate_trace() {
    return cpptrace::generate_trace(0, 200);
  }
};

struct DebugBucket {
  mutable std::recursive_mutex mutex;
  std::unordered_map<const Debug*, DebugItem> items;
};

static inline DebugBucket& GetDebugBucket(const Debug* ptr) {
  static std::array<DebugBucket, 32> buckets_{};
  return buckets_.at(std::hash<const Debug*>()(ptr) % buckets_.size());
}

struct Debug {
  Debug() {
    DebugBucket& bucket = GetDebugBucket(this);
    std::unique_lock lock{bucket.mutex};
    bucket.items[this].constructor();
  }

  ~Debug() {
    DebugBucket& bucket = GetDebugBucket(this);
    std::unique_lock lock{bucket.mutex};
    bucket.items[this].destructor();
  }

  Debug(Debug&& other) {
    DebugBucket& bucket = GetDebugBucket(this);
    DebugBucket& bucket_other = GetDebugBucket(&other);
    std::scoped_lock lock{bucket.mutex, bucket_other.mutex};
    bucket_other.items[&other].move_constructor_from();
    bucket.items[this].move_constructor_to();
  }

  Debug& operator=(Debug&& other) {
    DebugBucket& bucket = GetDebugBucket(this);
    DebugBucket& bucket_other = GetDebugBucket(&other);
    std::scoped_lock lock{bucket.mutex, bucket_other.mutex};
    bucket_other.items[&other].move_assignment_from();
    bucket.items[this].move_assignment_to();
    return *this;
  }

  Debug(const Debug& other) {
    DebugBucket& bucket = GetDebugBucket(this);
    DebugBucket& bucket_other = GetDebugBucket(&other);
    std::scoped_lock lock{bucket.mutex, bucket_other.mutex};
    bucket_other.items[&other].copy_constructor_from();
    bucket.items[this].copy_constructor_to();
  }

  Debug& operator=(const Debug& other) {
    DebugBucket& bucket = GetDebugBucket(this);
    DebugBucket& bucket_other = GetDebugBucket(&other);
    std::scoped_lock lock{bucket.mutex, bucket_other.mutex};
    bucket_other.items[&other].copy_assignment_from();
    bucket.items[this].copy_assignment_to();
    return *this;
  }

  void use() const {
    const DebugBucket& bucket = GetDebugBucket(this);
    std::unique_lock lock{bucket.mutex};
    auto it = bucket.items.find(this);
    if (it == bucket.items.end()) {
      throw std::runtime_error("Item not found");
    }
    it->second.use();
  }

 private:
};

template <typename T>
void DebugUse(const T* x) {
  x->debug_.use();
}

#define LARCH_DEBUG_THIS                 \
  template <typename DEBARGT_>           \
  friend void DebugUse(const DEBARGT_*); \
  Debug debug_
#define LARCH_DEBUG_USE debug_.use()

#else

#define LARCH_DEBUG_THIS static_assert(true)
#define LARCH_DEBUG_USE
#define DebugUse(x)

#endif
