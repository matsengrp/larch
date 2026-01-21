#pragma once

#include <compare>

/**
 * @brief Abstract score type for SPR move evaluation.
 *
 * This class wraps a score value that can represent different underlying
 * scoring models:
 * - For matOptimize/Fitch: parsimony score (1 point per mutation)
 * - For ML/Sankoff: log-likelihood based score (currently stub returning 0)
 *
 * The class provides comparison and arithmetic operations to allow unified
 * treatment of scores regardless of the underlying scoring model.
 */
class Score {
 public:
  Score() : value_{0} {}
  explicit Score(int v) : value_{v} {}

  bool operator<(const Score& o) const { return value_ < o.value_; }
  bool operator>(const Score& o) const { return value_ > o.value_; }
  bool operator<=(const Score& o) const { return value_ <= o.value_; }
  bool operator>=(const Score& o) const { return value_ >= o.value_; }
  bool operator==(const Score& o) const { return value_ == o.value_; }
  bool operator!=(const Score& o) const { return value_ != o.value_; }

  Score operator+(const Score& o) const { return Score{value_ + o.value_}; }
  Score operator-(const Score& o) const { return Score{value_ - o.value_}; }
  Score& operator+=(const Score& o) {
    value_ += o.value_;
    return *this;
  }
  Score& operator-=(const Score& o) {
    value_ -= o.value_;
    return *this;
  }

  Score operator-() const { return Score{-value_}; }

  int value() const { return value_; }

 private:
  int value_;
};
