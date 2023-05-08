#pragma once

#include <string_view>
#include <memory>
#include <type_traits>

#include "larch/common.hpp"
#include "larch/vcf/include_hts.hpp"

template <typename FnCreate, typename FnDestroy, typename... Args>
auto MakeCObject(FnCreate&& create, FnDestroy&& destroy, Args&&... args) {
  using type = std::invoke_result_t<FnCreate, Args...>;
  return std::unique_ptr<std::remove_pointer_t<type>, FnDestroy>(
      std::invoke(std::forward<FnCreate>(create), std::forward<Args>(args)...),
      std::forward<FnDestroy>(destroy));
}

namespace libhts {
using file = decltype(MakeCObject(hts_open, hts_close, "", ""));
using record = decltype(MakeCObject(bcf_init, bcf_destroy));
using header = decltype(
    MakeCObject(bcf_hdr_read, bcf_hdr_destroy, static_cast<htsFile*>(nullptr)));
};  // namespace libhts

class VcfFile {
 public:
  explicit VcfFile(std::string_view path)
      : hts_file_{MakeCObject(hts_open, hts_close, std::string{path}.c_str(), "r")},
        rec_{MakeCObject(bcf_init, bcf_destroy)},
        hdr_{MakeCObject(bcf_hdr_read, bcf_hdr_destroy, hts_file_.get())} {
    Assert(hts_file_ != nullptr);
    Assert(rec_ != nullptr);
    Assert(hdr_ != nullptr);
  }

 private:
  libhts::file hts_file_;
  libhts::record rec_;
  libhts::header hdr_;
};
