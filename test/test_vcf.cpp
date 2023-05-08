#include "test_common.hpp"

#include "larch/vcf/vcf_file.hpp"

static void test_vcf(std::string_view path) {
  VcfFile file{path};
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_vcf("data/new_samples.vcf"); },
              "VCF"});
