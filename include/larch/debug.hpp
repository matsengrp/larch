#pragma once

#include <cpptrace/cpptrace.hpp>
#include <cpptrace/formatting.hpp>

#include "larch/parallel/parallel_common.hpp"

enum class DebugState { Undef, Constructed, MovedFrom, Destructed };

enum class DebugConstructType { Undef, Default, Move, Copy };
enum class DebugAssignType { Undef, Move, Copy };
enum class DebugMoveFromType { Undef, Construct, Assign };

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

 private:
  [[noreturn]] void fail(const char* msg) const {
    cpptrace::formatter formatter = cpptrace::get_default_formatter();
    formatter.paths(cpptrace::formatter::path_mode::basename);
    formatter.snippets(true);
    formatter.snippet_context(5);
    std::cout << "\n\nLarchDebug: " << msg << "\n";
    print(generate_trace(), formatter);
    std::cout << "\nConstructed at trace:\n";
    print(construct_trace, formatter);
    std::cout << "\nMoved from trace:\n";
    print(move_from_trace, formatter);
    std::cout << "\nDestructed trace:\n";
    print(destruct_trace, formatter);
    Fail(msg);
  }

  void print(const cpptrace::stacktrace& trace,
             const cpptrace::formatter& formatter) const {
    for (auto& i : trace.frames) {
      if (i.filename.find("include/larch") != std::string::npos) {
        formatter.print(std::cout, i);
        std::cout << "\n\n";
      }
    }
  }

  cpptrace::stacktrace generate_trace() const {
    // return {};
    return cpptrace::generate_trace();
  }
};

struct DebugBucket {
  mutable std::mutex mutex;
  std::unordered_map<const void*, DebugItem> items;
};

static inline DebugBucket& GetDebugBucket(const void* ptr) {
  static std::array<DebugBucket, 32> buckets_{};
  return buckets_.at(std::hash<const void*>()(ptr) % buckets_.size());
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
    std::unique_lock lock{bucket.mutex};
    if (&bucket == &bucket_other) {
      bucket.items[this].move_constructor_to();
      bucket_other.items[&other].move_constructor_from();
    } else {
      std::unique_lock lock_other{bucket_other.mutex};
      bucket.items[this].move_constructor_to();
      bucket_other.items[&other].move_constructor_from();
    }
  }

  Debug& operator=(Debug&& other) {
    DebugBucket& bucket = GetDebugBucket(this);
    DebugBucket& bucket_other = GetDebugBucket(&other);
    std::unique_lock lock{bucket.mutex};
    if (&bucket == &bucket_other) {
      bucket.items[this].move_assignment_to();
      bucket_other.items[&other].move_assignment_from();
    } else {
      std::unique_lock lock_other{bucket_other.mutex};
      bucket.items[this].move_assignment_to();
      bucket_other.items[&other].move_assignment_from();
    }
    return *this;
  }

  Debug(const Debug& other) {
    DebugBucket& bucket = GetDebugBucket(this);
    DebugBucket& bucket_other = GetDebugBucket(&other);
    std::unique_lock lock{bucket.mutex};
    if (&bucket == &bucket_other) {
      bucket.items[this].copy_constructor_to();
      bucket_other.items[&other].copy_constructor_from();
    } else {
      std::unique_lock lock_other{bucket_other.mutex};
      bucket.items[this].copy_constructor_to();
      bucket_other.items[&other].copy_constructor_from();
    }
  }

  Debug& operator=(const Debug& other) {
    DebugBucket& bucket = GetDebugBucket(this);
    DebugBucket& bucket_other = GetDebugBucket(&other);
    std::unique_lock lock{bucket.mutex};
    if (&bucket == &bucket_other) {
      bucket.items[this].copy_assignment_to();
      bucket_other.items[&other].copy_assignment_from();
    } else {
      std::unique_lock lock_other{bucket_other.mutex};
      bucket.items[this].copy_assignment_to();
      bucket_other.items[&other].copy_assignment_from();
    }
    return *this;
  }

  void use() const {
    const DebugBucket& bucket = GetDebugBucket(this);
    std::unique_lock lock{bucket.mutex};
    auto it = bucket.items.find(this);
    if (it == bucket.items.end()) {
      Fail("Item not found");
    }
    it->second.use();
  }

 private:
};
