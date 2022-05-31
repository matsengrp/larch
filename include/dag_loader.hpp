#pragma once

#include <string_view>
#include <unordered_map>

#include "merge.hpp"

[[nodiscard]] DAG LoadDAGFromProtobufGZ(std::string_view path, std::string& ref_seq,
                                        std::vector<Mutations>& mutations);

[[nodiscard]] DAG LoadTreeFromProtobufGZ(std::string_view path,
                                         std::vector<Mutations>& mutations);

[[nodiscard]] std::string LoadRefseqFromJsonGZ(std::string_view path);

[[nodiscard]] DAG LoadDAGFromJsonGZ(std::string_view path, std::string& refseq);

[[nodiscard]] std::vector<CompactGenome> LoadCompactGenomesJsonGZ(
    std::string_view path);

void StoreDAGToProtobuf(const DAG& dag, std::string_view reference_sequence,
                        const std::vector<Mutations>& edge_parent_mutations,
                        std::string_view path);
