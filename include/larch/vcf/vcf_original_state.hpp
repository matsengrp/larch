#pragma once

#include "larch/vcf/vcf_file.hpp"
#include "larch/usher_glue.hpp"

Original_State_t MakeOriginalState(VcfFile& vcf) {
  Original_State_t result;
  while (auto record = vcf.read()) {
    auto ref_nuc = static_cast<uint8_t>(MAT::get_nuc_id(record->GetRef().at(0)));
    MAT::Mutation mut_template{std::string{record->GetChrom(vcf.GetHeader())},
                               static_cast<int>(record->GetPos()),
                               0,
                               ref_nuc,
                               0,
                               ref_nuc};
    std::vector<nuc_one_hot> allele_translated;
    for (char i : record->GetAlt()) {
      if (i == ',') {
        continue;
      }
      allele_translated.push_back(static_cast<uint8_t>(MAT::get_nuc_id(i)));
    }
    std::vector<MAT::Mutation> mutations;
    for (auto i : record->GetAlleles()) {
      MAT::Mutation this_mut(mut_template);
      this_mut.set_mut_one_hot(allele_translated[allele_idx - 1]);
      this_mut.set_descendant_mut(allele_translated[allele_idx - 1]);
      this_mut.set_auxillary(allele_translated[allele_idx - 1], 0);
      mutations.push_back(this_mut);
    }
  }
  return result;
}
