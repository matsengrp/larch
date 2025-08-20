#pragma once

#include <iostream>
#include <fstream>
#include <vector>

/**
 * @brief Handles binary serialization and deserialization of phylogenetic DAGs.
 * 
 * DagbinFileIO provides a custom binary format for efficiently storing and loading
 * phylogenetic DAG structures to/from disk. The format uses a section-based layout
 * with a magic number header for file validation. Sections include header metadata,
 * reference sequence, nodes, edges, and compact genomes. The format supports both
 * complete DAG writing and incremental appending operations.
 * 
 * The binary format is structured as:
 * - Magic number (6 bytes: "DAGBIN")
 * - Linked list of section offsets
 * - Individual sections (Header, RefSeq, Nodes, Edges, CompactGenomes)
 */
class DagbinFileIO {
 public:
  /**
   * @brief Identifiers for different sections in the dagbin file format.
   */
  enum class SectionId : char {
    Header = 'H',         ///< Header section containing metadata
    RefSeq = 'R',         ///< Reference sequence section
    Nodes = 'N',          ///< Nodes section
    Edges = 'E',          ///< Edges section  
    CompactGenomes = 'C', ///< Compact genomes section
    Test = 'X'            ///< Test section for debugging
  };

  friend std::ostream &operator<<(std::ostream &os, SectionId section_id) {
    std::cout << "SectionId::" << static_cast<char>(section_id);
    return os;
  }

  /// Magic number bytes spelling "DAGBIN" for file format identification
  inline static const std::vector<unsigned char> MAGIC_NUMBER = {0x44, 0x41, 0x47,
                                                                 0x42, 0x49, 0x4E};

  /**
   * @brief Metadata header for dagbin files.
   * 
   * Contains counts of various DAG components for validation and pre-allocation
   * during deserialization.
   */
  struct Header {
    size_t node_count = 0;  ///< Total number of nodes in the DAG
    size_t edge_count = 0;  ///< Total number of edges in the DAG
    size_t leaf_count = 0;  ///< Number of leaf nodes
  };

  /// Batch size for processing nodes/edges in chunks
  static const size_t batch_size = 250;

  /**
   * @brief Checks if a file is in dagbin format by verifying the magic number.
   * @param path Path to the file to check
   * @return true if file has valid dagbin magic number, false otherwise
   */
  inline static bool IsFileDagbinFormat(std::string_view path);

  /**
   * @brief Reads a complete DAG from a dagbin file.
   * @param path Path to the dagbin file
   * @return Deserialized DAG storage object
   */
  inline static MADAGStorage<> ReadDAG(std::string_view path);

  /**
   * @brief Writes a complete DAG to a new dagbin file.
   * @tparam DAG Type of the DAG to serialize
   * @param dag The DAG to write
   * @param path Output file path
   */
  template <typename DAG>
  inline static void WriteDAG(DAG dag, std::string_view path);

  /**
   * @brief Appends DAG data to an existing dagbin file.
   * @tparam DAG Type of the DAG to append
   * @param dag The DAG data to append
   * @param path Path to the existing dagbin file
   */
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

#include "larch/impl/dagbin_fileio_impl.hpp"
