#include <netam/pcp_dataframe.hpp>
#include <netam/common.hpp>

#include <zlib.h>

#include <array>
#include <memory>

namespace netam {

pcp_dataframe::pcp_dataframe(const std::filesystem::path& csv_gz_path)
    : path_{csv_gz_path} {}

generator<std::string> pcp_dataframe::lines(const std::filesystem::path& path) {
  auto close = [](::gzFile x) static {
    if (x != nullptr) {
      ::gzclose(x);
    }
  };
  std::unique_ptr<std::remove_pointer_t<gzFile>, decltype(close)> file{
      ::gzopen(path.c_str(), "rb"), close};

  if (file == nullptr) {
    fail("Can't open gzip file");
  }

  std::string buffer;
  std::array<char, 4096> chunk;

  for (int bytes = ::gzread(file.get(), chunk.data(), chunk.size()); bytes > 0;
       bytes = ::gzread(file.get(), chunk.data(), chunk.size())) {
    buffer.append(chunk.data(), unsigned_cast(bytes));

    std::size_t pos;
    while ((pos = buffer.find('\n')) != std::string::npos) {
      co_yield buffer.substr(0, pos);
      buffer.erase(0, pos + 1);
    }
  }

  if (not buffer.empty()) {
    co_yield buffer;
  }
}
}  // namespace netam
