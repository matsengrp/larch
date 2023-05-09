#include "test_common.hpp"

#include "larch/vcf/vcf_file.hpp"

static void test_vcf(std::string_view path) {
  VcfFile file{path};
  while (auto record = file.read()) {
    std::cout << " --- Record ---\n";
    std::cout << "Pos: " << record->GetPos() << "\n";
    std::cout << "RefLength: " << record->GetRefLength() << "\n";
    std::cout << "ID: " << record->GetId() << "\n";
    std::cout << "Ref: " << record->GetRef() << "\n";
    std::cout << "Qual: " << record->GetQual() << "\n";
  }
}

[[maybe_unused]] static const auto test_added0 =
    add_test({[] { test_vcf("data/new_samples.vcf"); }, "VCF"});
