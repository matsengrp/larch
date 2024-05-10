#include "test_common.hpp"
#include "test_common_dag.hpp"

#include "larch/vcf/vcf_original_state.hpp"

static void test_vcf() {
  VcfFile file{"data/MakeSampleDAG_vcf.vcf"};
  auto orig_state = MakeOriginalState(file);
  auto dag_storage = AddMATConversion(MakeSampleDAG());
  MAT::Tree tree;
  dag_storage.View().BuildMAT(tree);
  reassign_states(tree, orig_state);

  // for (auto i : file.GetSamples()) {
  //   std::cout << "Sample: " << i << "\n";
  // }
  // while (auto record = file.read()) {
  //   std::cout << " --- Record ---\n";
  //   std::cout << "Pos: " << record->GetPos() << "\n";
  //   std::cout << "ChromId: " << record->GetChromId() << "\n";
  //   std::cout << "Chrom: " << record->GetChrom(file.GetHeader()) << "\n";
  //   std::cout << "RefLength: " << record->GetRefLength() << "\n";
  //   std::cout << "ID: " << record->GetId() << "\n";
  //   std::cout << "Ref: " << record->GetRef() << "\n";
  //   std::cout << "Alt: " << record->GetAlt() << "\n";
  //   std::cout << "Qual: " << record->GetQual() << "\n";
  //   for (auto i : record->GetAlleles()) {
  //     std::cout << "Allele: " << i << "\n";
  //   }
  // }
}

[[maybe_unused]] static const auto test_added0 = add_test({[] { test_vcf(); }, "VCF"});
