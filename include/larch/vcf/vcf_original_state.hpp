#pragma once

#include "larch/vcf/vcf_file.hpp"
#include "larch/usher_glue.hpp"

Original_State_t MakeOriginalState(VcfFile& vcf) {
  Original_State_t result;
  size_t index = 0;
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
    Mutation_Set mutations;
    for (nuc_one_hot i :
         record->GetAlleles() |
             ranges::views::transform([](auto sample) -> nuc_one_hot {
               return static_cast<uint8_t>(MAT::get_nuc_id(sample.at(0)));
             })) {
      MAT::Mutation this_mut(mut_template);
      this_mut.set_mut_one_hot(i);
      this_mut.set_descendant_mut(i);
      this_mut.set_auxillary(i, 0);
      mutations.insert(std::move(this_mut));
    }
    result.insert({index++, std::move(mutations)});
  }
  return result;
}
