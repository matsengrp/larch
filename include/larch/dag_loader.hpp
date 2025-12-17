/**
 * Functions for loading and writing mutation annotated dags.
 */
#pragma once

#include <string_view>
#include <unordered_map>

#include "larch/merge/merge.hpp"
#include "larch/mat_conversion.hpp"

enum class FileFormat {
  Infer,
  Dagbin,
  Protobuf,
  ProtobufDAG,
  ProtobufTree,
  JsonDAG,
  DebugAll
};

const std::vector<std::pair<std::string, FileFormat>> file_extension_names = {
    {"dagbin", FileFormat::Dagbin},        {"protobuf", FileFormat::Protobuf},
    {"pb", FileFormat::Protobuf},          {"dag-pb", FileFormat::ProtobufDAG},
    {"pb_dag", FileFormat::ProtobufDAG},   {"tree-pb", FileFormat::ProtobufTree},
    {"pb_tree", FileFormat::ProtobufTree}, {"json", FileFormat::JsonDAG},
    {"dag-json", FileFormat::JsonDAG},     {"json_dag", FileFormat::JsonDAG},
    {"debug-all", FileFormat::DebugAll}};

inline FileFormat InferFileFormat(std::string_view path);

[[nodiscard]] MADAGStorage<> LoadDAG(
    std::string_view input_dag_path, FileFormat file_format = FileFormat::Infer,
    std::optional<std::string> refseq_path = std::nullopt);

[[nodiscard]] MADAGStorage<> LoadDAGFromProtobuf(std::string_view path);

[[nodiscard]] MADAGStorage<> LoadDAGFromDagbin(std::string_view path);

[[nodiscard]] MADAGStorage<> LoadTreeFromProtobuf(std::string_view path,
                                                  std::string_view reference_sequence);

/*
[[nodiscard]] MADAGStorage LoadTreeFromProtobuf(std::string_view path,
                                                std::string_view reference_sequence,
                                                std::string_view vcf_path);
*/

[[nodiscard]] MADAGStorage<> LoadDAGFromJson(std::string_view path);

[[nodiscard]] MADAGStorage<> LoadTreeFromFastaNewick(std::string_view fasta_path,
                                                     std::string_view newick_path,
                                                     std::string_view reference_path);

[[nodiscard]] std::string LoadReferenceSequence(std::string_view path);

[[nodiscard]] std::unordered_map<std::string, std::string> LoadFasta(
    std::string_view path);

using CompactGenomeData = ContiguousMap<MutationPosition, MutationBase>;

std::unordered_map<std::string, CompactGenomeData> LoadCompactGenomeDataFromVCF(
    const std::string& path, const std::string& ref_seq);

void MADAGApplyCompactGenomeData(
    MADAGStorage<>& dag_storage,
    const std::unordered_map<std::string, CompactGenomeData>& mut_map,
    bool silence_warnings = true);

void LoadVCFData(MADAGStorage<>& dag_storage, std::string& vcf_path,
                 bool silence_warnings = true);

template <typename DAG>
void StoreDAG(DAG dag, std::string_view output_dag_path,
              FileFormat file_format = FileFormat::Infer, bool append_changes = false);

template <typename DAG>
void StoreDAGToProtobuf(DAG dag, std::string_view path);

template <typename DAG>
void StoreDAGToDagbin(DAG dag, std::string_view path, bool append_changes = false);

template <typename DAG>
void StoreTreeToProtobuf(DAG dag, std::string_view path);

template <typename DAG, typename iostream>
void MADAGToDOT(DAG dag, iostream& out);

template <typename DAG, typename iostream>
void FragmentToDOT(DAG dag, const std::vector<EdgeId>& edges, iostream& out);

template <typename iostream>
void NewickToDOT(const std::string& newick, iostream& out);

std::string ToEdgeMutationsString(const MAT::Node* node);

template <typename iostream>
void MATToDOT(const MAT::Tree& mat, iostream& out);

#include "larch/impl/dag_loader_impl.hpp"
