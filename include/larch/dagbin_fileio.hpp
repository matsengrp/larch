#pragma once

#include <iostream>
#include <fstream>
#include <vector>

class DagbinFileIO {
 public:
  enum class SectionId : char {
    Header = 'H',
    RefSeq = 'R',
    Nodes = 'N',
    Edges = 'E',
    CompactGenomes = 'C',
    Test = 'X'
  };

  friend std::ostream &operator<<(std::ostream &os, SectionId section_id) {
    std::cout << "SectionId::" << static_cast<char>(section_id);
    return os;
  }

  // Magic number is `dagbin`.
  inline static const std::vector<unsigned char> MAGIC_NUMBER = {0x44, 0x41, 0x47,
                                                                 0x42, 0x49, 0x4E};

  struct Header {
    size_t node_count = 0;
    size_t edge_count = 0;
    size_t leaf_count = 0;
  };

  static const size_t batch_size = 250;

  /* Uses magic number to check if file is in dagbin format */
  inline static bool IsFileDagbinFormat(std::string_view path);

  inline static MADAGStorage<> ReadDAG(std::string_view path);

  template <typename DAG>
  inline static void WriteDAG(DAG dag, std::string_view path);

  template <typename DAG>
  inline static void AppendDAG(DAG dag, std::string_view path);

 private:
  template <typename T, typename iostream>
  inline static T ReadData(iostream &infile);

  template <typename T, typename iostream>
  inline static void WriteData(iostream &outfile, const T data);

  template <typename iostream>
  inline static std::string ReadString(iostream &infile);

  template <typename iostream>
  inline static void WriteString(iostream &outfile, const std::string &str);

  template <typename iostream>
  inline static std::vector<std::streampos> ReadLinkedList(iostream &infile);

  template <typename iostream>
  inline static std::vector<std::pair<std::streampos, SectionId>> ReadLabeledLinkedList(
      iostream &infile);

  template <typename iostream>
  inline static void WriteLinkedList(iostream &outfile,
                                     const std::vector<std::streampos> &offsets);

  template <typename iostream>
  inline static bool CheckMagicNumber(iostream &infile, bool do_assert = true);

  template <typename iostream>
  inline static void WriteMagicNumber(iostream &outfile);

  template <typename iostream, typename DAG>
  inline static Header ReadHeader(iostream &infile, DAG dag);

  template <typename iostream, typename DAG>
  inline static void WriteHeader(iostream &outfile, const DAG dag);

  template <typename iostream, typename DAG>
  inline static void ReadReferenceSequence(iostream &infile, DAG dag);

  template <typename iostream, typename DAG>
  inline static void WriteReferenceSequence(iostream &outfile, const DAG dag);

  template <typename iostream, typename DAG>
  inline static void ReadNodes(iostream &infile, DAG dag);

  template <typename iostream, typename DAG>
  inline static void WriteNodes(iostream &outfile, const DAG dag,
                                std::optional<size_t> min_id_opt = std::nullopt,
                                std::optional<size_t> max_id_opt = std::nullopt);

  template <typename iostream, typename DAG>
  inline static void ReadEdges(iostream &infile, DAG dag);

  template <typename iostream, typename DAG>
  inline static void WriteEdges(iostream &outfile, const DAG dag,
                                std::optional<size_t> min_id_opt = std::nullopt,
                                std::optional<size_t> max_id_opt = std::nullopt);
};

#pragma GCC push_options
#pragma GCC optimize ("O0")
#pragma GCC visibility push(default)

#include "larch/impl/dagbin_fileio_impl.hpp"

#pragma GCC visibility pop
#pragma GCC pop_options