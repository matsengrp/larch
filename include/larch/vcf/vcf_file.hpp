#pragma once

#include <string_view>
#include <memory>
#include <type_traits>

#include "larch/common.hpp"
#include "larch/vcf/include_hts.hpp"

template <typename T, auto* FnDestroy>
struct CObject {
  struct Deleter {
    void operator()(T* ptr) const noexcept { FnDestroy(ptr); };
  };
  using ptr = std::unique_ptr<T, Deleter>;
  static ptr null() { return {nullptr, {}}; }
  template <auto* FnCreate, typename... Args>
  static ptr make(Args&&... args) {
    T* result = FnCreate(std::forward<Args>(args)...);
    Assert(result != nullptr);
    return {result, {}};
  }
};

namespace libhts {
using file = CObject<::htsFile, ::hts_close>;
using header = CObject<::bcf_hdr_t, ::bcf_hdr_destroy>;
using record = CObject<::bcf1_t, ::bcf_destroy>;

file::ptr make_file(std::string_view path, std::string_view mode) {
  return file::make<::hts_open>(std::string{path}.c_str(), std::string{mode}.c_str());
}

header::ptr make_header(const file::ptr& file) {
  Assert(file);
  return header::make<::bcf_hdr_read>(file.get());
}

record::ptr make_record() { return record::make<::bcf_init>(); }

};  // namespace libhts

class VcfFile {
 public:
  explicit VcfFile(std::string_view path) {
    file_ = libhts::make_file(path, "r");
    header_ = libhts::make_header(file_);
    record_ = libhts::make_record();
  }

 private:
  libhts::file::ptr file_ = libhts::file::null();
  libhts::header::ptr header_ = libhts::header::null();
  libhts::record::ptr record_ = libhts::record::null();
};
