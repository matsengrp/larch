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

class VcfRecord {
 public:
  VcfRecord() : record_{libhts::make_record()} {}

  hts_pos_t GetPos() const { return record_->pos; }

  int32_t GetChromId() const { return record_->rid; }

  std::string_view GetChrom(const libhts::header::ptr& header) const {
    return ::bcf_hdr_id2name(header.get(), record_->rid);
  }

  hts_pos_t GetRefLength() const { return record_->rlen; }

  std::string_view GetId() const { return record_->d.id; }

  std::string_view GetRef() const { return record_->d.als; }

  std::string_view GetAlt() const { return record_->d.als + GetRef().size() + 1; }

  float GetQual() const { return record_->qual; }

  auto GetAlleles() {
    size_t count = static_cast<size_t>(record_->n_allele);
    return ranges::subrange(record_->d.allele, record_->d.allele + count) |
           ranges::views::transform([](char* i) { return std::string_view(i); });
  }

 private:
  friend class VcfFile;
  libhts::record::ptr record_;
};

class VcfFile {
 public:
  explicit VcfFile(std::string_view path) {
    file_ = libhts::make_file(path, "r");
    header_ = libhts::make_header(file_);
  }

  const libhts::header::ptr& GetHeader() { return header_; }

  auto GetSamples() {
    auto* hdr = header_.get();
    size_t count = static_cast<size_t>(bcf_hdr_nsamples(hdr));
    return ranges::subrange(hdr->samples, hdr->samples + count) |
           ranges::views::transform([](char* i) { return std::string_view(i); });
  }

  std::optional<VcfRecord> read() {
    VcfRecord result{};
    int error = ::bcf_read(file_.get(), header_.get(), result.record_.get());
    if (error < -1) {
      std::cerr << "bcf_read failed with " << error << "\n";
      Fail("bcf_read");
    } else if (error == 0) {
      Assert(::bcf_unpack(result.record_.get(), BCF_UN_ALL) == 0);
      if (result.GetPos() == 0) {
        return std::nullopt;
      }
      return result;
    } else {
      return std::nullopt;
    }
  }

 private:
  libhts::file::ptr file_ = libhts::file::null();
  libhts::header::ptr header_ = libhts::header::null();
};
