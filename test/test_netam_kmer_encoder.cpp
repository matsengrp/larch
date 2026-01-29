// Tests for kmer_sequence_encoder
#include "test_common.hpp"

#ifdef USE_NETAM
#include <netam/kmer_sequence_encoder.hpp>

#include <torch/torch.h>

using namespace netam;

namespace {

// Helper to create a minimal YAML config for testing
YAML::Node make_config(std::size_t kmer_length, std::size_t site_count) {
  YAML::Node yaml;
  yaml["kmer_length"] = kmer_length;
  yaml["site_count"] = site_count;
  return yaml;
}

// ============================================================================
// encode_bases tests (static method)
// ============================================================================

void test_encode_bases_standard() {
  auto result = kmer_sequence_encoder::encode_bases("ACGT");
  TestAssert(result.size(0) == 4);
  TestAssert(result[0].item<int64_t>() == 0);  // A
  TestAssert(result[1].item<int64_t>() == 1);  // C
  TestAssert(result[2].item<int64_t>() == 2);  // G
  TestAssert(result[3].item<int64_t>() == 3);  // T
}

void test_encode_bases_lowercase() {
  auto result = kmer_sequence_encoder::encode_bases("acgt");
  TestAssert(result.size(0) == 4);
  TestAssert(result[0].item<int64_t>() == 0);  // a -> A
  TestAssert(result[1].item<int64_t>() == 1);  // c -> C
  TestAssert(result[2].item<int64_t>() == 2);  // g -> G
  TestAssert(result[3].item<int64_t>() == 3);  // t -> T
}

void test_encode_bases_unknown() {
  auto result = kmer_sequence_encoder::encode_bases("N");
  TestAssert(result.size(0) == 1);
  TestAssert(result[0].item<int64_t>() == 4);  // N -> 4 (unknown)
}

void test_encode_bases_mixed_with_unknown() {
  auto result = kmer_sequence_encoder::encode_bases("ANCT");
  TestAssert(result.size(0) == 4);
  TestAssert(result[0].item<int64_t>() == 0);  // A
  TestAssert(result[1].item<int64_t>() == 4);  // N -> unknown
  TestAssert(result[2].item<int64_t>() == 1);  // C
  TestAssert(result[3].item<int64_t>() == 3);  // T
}

void test_encode_bases_empty() {
  auto result = kmer_sequence_encoder::encode_bases("");
  TestAssert(result.size(0) == 0);
}

// ============================================================================
// Constructor and property tests
// ============================================================================

void test_constructor_kmer_length_3() {
  auto yaml = make_config(3, 100);
  kmer_sequence_encoder encoder{yaml};

  TestAssert(encoder.kmer_length() == 3);
  TestAssert(encoder.site_count() == 100);
  // 4^3 + 1 (for N placeholder) = 65
  TestAssert(encoder.kmer_count() == 65);
}

void test_constructor_kmer_length_5() {
  auto yaml = make_config(5, 200);
  kmer_sequence_encoder encoder{yaml};

  TestAssert(encoder.kmer_length() == 5);
  TestAssert(encoder.site_count() == 200);
  // 4^5 + 1 = 1025
  TestAssert(encoder.kmer_count() == 1025);
}

// ============================================================================
// encode_sequence tests
// ============================================================================

void test_encode_sequence_output_shape() {
  auto yaml = make_config(3, 10);
  kmer_sequence_encoder encoder{yaml};

  auto [encoded, wt_modifier] = encoder.encode_sequence("ACGTACGT");

  // encoded should have site_count elements
  TestAssert(encoded.size(0) == 10);
  // wt_modifier should be [site_count, 4]
  TestAssert(wt_modifier.size(0) == 10);
  TestAssert(wt_modifier.size(1) == 4);
}

void test_encode_sequence_short_sequence() {
  // Sequence shorter than site_count
  auto yaml = make_config(3, 100);
  kmer_sequence_encoder encoder{yaml};

  auto [encoded, wt_modifier] = encoder.encode_sequence("ACG");

  TestAssert(encoded.size(0) == 100);
  TestAssert(wt_modifier.size(0) == 100);
}

void test_encode_sequence_lowercase() {
  auto yaml = make_config(3, 10);
  kmer_sequence_encoder encoder{yaml};

  auto [encoded_lower, wt_lower] = encoder.encode_sequence("acgtacgt");
  auto [encoded_upper, wt_upper] = encoder.encode_sequence("ACGTACGT");

  // Should produce identical results
  TestAssert(torch::equal(encoded_lower, encoded_upper));
  TestAssert(torch::equal(wt_lower, wt_upper));
}

void test_encode_sequence_wt_modifier_values() {
  auto yaml = make_config(3, 5);
  kmer_sequence_encoder encoder{yaml};

  auto [encoded, wt_modifier] = encoder.encode_sequence("ACG");

  // Position 0 has 'A' (index 0), so wt_modifier[0][0] should be -BIG
  TestAssert(wt_modifier[0][0].item<float>() < -1e8f);
  TestAssert(wt_modifier[0][1].item<float>() == 0.0f);
  TestAssert(wt_modifier[0][2].item<float>() == 0.0f);
  TestAssert(wt_modifier[0][3].item<float>() == 0.0f);

  // Position 1 has 'C' (index 1)
  TestAssert(wt_modifier[1][0].item<float>() == 0.0f);
  TestAssert(wt_modifier[1][1].item<float>() < -1e8f);
  TestAssert(wt_modifier[1][2].item<float>() == 0.0f);
  TestAssert(wt_modifier[1][3].item<float>() == 0.0f);

  // Position 2 has 'G' (index 2)
  TestAssert(wt_modifier[2][0].item<float>() == 0.0f);
  TestAssert(wt_modifier[2][1].item<float>() == 0.0f);
  TestAssert(wt_modifier[2][2].item<float>() < -1e8f);
  TestAssert(wt_modifier[2][3].item<float>() == 0.0f);

  // Positions beyond sequence length should be all zeros
  TestAssert(wt_modifier[3][0].item<float>() == 0.0f);
  TestAssert(wt_modifier[3][1].item<float>() == 0.0f);
  TestAssert(wt_modifier[3][2].item<float>() == 0.0f);
  TestAssert(wt_modifier[3][3].item<float>() == 0.0f);
}

void test_encode_sequence_with_n_in_middle() {
  auto yaml = make_config(3, 10);
  kmer_sequence_encoder encoder{yaml};

  auto [encoded, wt_modifier] = encoder.encode_sequence("ACNGT");

  // wt_modifier at position 2 (N) should be all zeros (no base to mask)
  TestAssert(wt_modifier[2][0].item<float>() == 0.0f);
  TestAssert(wt_modifier[2][1].item<float>() == 0.0f);
  TestAssert(wt_modifier[2][2].item<float>() == 0.0f);
  TestAssert(wt_modifier[2][3].item<float>() == 0.0f);
}

void test_encode_sequence_empty() {
  auto yaml = make_config(3, 10);
  kmer_sequence_encoder encoder{yaml};

  auto [encoded, wt_modifier] = encoder.encode_sequence("");

  // Should still produce tensors of correct shape
  TestAssert(encoded.size(0) == 10);
  TestAssert(wt_modifier.size(0) == 10);
  TestAssert(wt_modifier.size(1) == 4);

  // All wt_modifier should be zeros (no bases to mask)
  TestAssert(torch::all(wt_modifier == torch::tensor(0)).item<bool>());
}

void test_encode_sequence_single_base() {
  auto yaml = make_config(3, 5);
  kmer_sequence_encoder encoder{yaml};

  auto [encoded, wt_modifier] = encoder.encode_sequence("T");

  TestAssert(encoded.size(0) == 5);
  // Position 0 has 'T' (index 3)
  TestAssert(wt_modifier[0][3].item<float>() < -1e8f);
}

void test_encode_sequence_all_same_base() {
  auto yaml = make_config(3, 5);
  kmer_sequence_encoder encoder{yaml};

  auto [encoded, wt_modifier] = encoder.encode_sequence("AAAAA");

  // All valid k-mers should be "AAA" -> same index
  // First position is NNA (due to padding), last is AAN
  // But middle positions should all encode AAA
  auto encoded_accessor = encoded.accessor<int32_t, 1>();

  // Positions 1, 2, 3 should have the same k-mer index (AAA)
  TestAssert(encoded_accessor[1] == encoded_accessor[2]);
  TestAssert(encoded_accessor[2] == encoded_accessor[3]);
}

void test_encode_sequence_kmer_indices_valid() {
  auto yaml = make_config(3, 10);
  kmer_sequence_encoder encoder{yaml};

  auto [encoded, wt_modifier] = encoder.encode_sequence("ACGTACGTAC");

  auto encoded_accessor = encoded.accessor<int32_t, 1>();

  // All indices should be within valid range [0, kmer_count)
  for (int64_t i = 0; i < encoded.size(0); ++i) {
    TestAssert(encoded_accessor[i] >= 0);
    TestAssert(static_cast<std::size_t>(encoded_accessor[i]) <
               encoder.kmer_count());
  }
}

void test_encode_sequence_padding_creates_zero_index() {
  // K-mers containing N should map to index 0
  auto yaml = make_config(3, 5);
  kmer_sequence_encoder encoder{yaml};

  auto [encoded, wt_modifier] = encoder.encode_sequence("ACG");

  auto encoded_accessor = encoded.accessor<int32_t, 1>();

  // Position 0: k-mer is "NAC" (contains N from padding) -> should be 0
  TestAssert(encoded_accessor[0] == 0);

  // Position 2: k-mer is "CGN" (contains N from padding) -> should be 0
  TestAssert(encoded_accessor[2] == 0);
}

void test_encode_sequence_long_sequence_truncated() {
  auto yaml = make_config(3, 5);
  kmer_sequence_encoder encoder{yaml};

  // Sequence longer than site_count
  auto [encoded, wt_modifier] = encoder.encode_sequence("ACGTACGTACGT");

  // Output should be limited to site_count
  TestAssert(encoded.size(0) == 5);
  TestAssert(wt_modifier.size(0) == 5);
}

}  // namespace

[[maybe_unused]] static bool reg_test_kmer_encoder =
    add_test({test_encode_bases_standard,
              "Netam KmerEncoder: encode_bases standard", {"netam"}}) &&
    add_test({test_encode_bases_lowercase,
              "Netam KmerEncoder: encode_bases lowercase", {"netam"}}) &&
    add_test({test_encode_bases_unknown,
              "Netam KmerEncoder: encode_bases unknown", {"netam"}}) &&
    add_test({test_encode_bases_mixed_with_unknown,
              "Netam KmerEncoder: encode_bases mixed with unknown",
              {"netam"}}) &&
    add_test({test_encode_bases_empty, "Netam KmerEncoder: encode_bases empty",
              {"netam"}}) &&
    add_test({test_constructor_kmer_length_3,
              "Netam KmerEncoder: constructor kmer_length 3", {"netam"}}) &&
    add_test({test_constructor_kmer_length_5,
              "Netam KmerEncoder: constructor kmer_length 5", {"netam"}}) &&
    add_test({test_encode_sequence_output_shape,
              "Netam KmerEncoder: encode_sequence output shape", {"netam"}}) &&
    add_test({test_encode_sequence_short_sequence,
              "Netam KmerEncoder: encode_sequence short sequence",
              {"netam"}}) &&
    add_test({test_encode_sequence_lowercase,
              "Netam KmerEncoder: encode_sequence lowercase", {"netam"}}) &&
    add_test({test_encode_sequence_wt_modifier_values,
              "Netam KmerEncoder: encode_sequence wt_modifier values",
              {"netam"}}) &&
    add_test({test_encode_sequence_with_n_in_middle,
              "Netam KmerEncoder: encode_sequence with N in middle",
              {"netam"}}) &&
    add_test({test_encode_sequence_empty,
              "Netam KmerEncoder: encode_sequence empty", {"netam"}}) &&
    add_test({test_encode_sequence_single_base,
              "Netam KmerEncoder: encode_sequence single base", {"netam"}}) &&
    add_test({test_encode_sequence_all_same_base,
              "Netam KmerEncoder: encode_sequence all same base", {"netam"}}) &&
    add_test({test_encode_sequence_kmer_indices_valid,
              "Netam KmerEncoder: encode_sequence kmer indices valid",
              {"netam"}}) &&
    add_test({test_encode_sequence_padding_creates_zero_index,
              "Netam KmerEncoder: encode_sequence padding creates zero index",
              {"netam"}}) &&
    add_test({test_encode_sequence_long_sequence_truncated,
              "Netam KmerEncoder: encode_sequence long sequence truncated",
              {"netam"}});

#endif  // USE_NETAM
