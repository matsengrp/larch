#pragma once

#include <string_view>

#include "larch/common.hpp"
#include "larch/vcf/include_hts.hpp"

class VcfFile {
public:
  explicit VcfFile(std::string_view path) {
    hts_file_ = hts_open(std::string{path}.c_str(), "r");
    Assert(hts_file_ != nullptr);
    rec_ = bcf_init1();
    Assert(rec_ != nullptr);
    hdr_ = bcf_hdr_read(hts_file_);
    Assert(hdr_ != nullptr);
  }
private:
  htsFile* hts_file_ = nullptr;
  bcf1_t* rec_ = nullptr;
  bcf_hdr_t* hdr_ = nullptr;
};
