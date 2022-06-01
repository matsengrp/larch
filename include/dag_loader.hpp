#pragma once

#include <string_view>
#include <unordered_map>

#include "merge.hpp"

[[nodiscard]] DAG LoadDAGFromProtobuf(std::string_view path,
                                      std::string& reference_sequence,
                                      std::vector<Mutations>& mutations);

[[nodiscard]] DAG LoadTreeFromProtobuf(std::string_view path,
                                       std::vector<Mutations>& mutations);

[[nodiscard]] std::string LoadRefseqFromJson(std::string_view path);

[[nodiscard]] DAG LoadDAGFromJson(std::string_view path,
                                  std::string& reference_sequence);

[[nodiscard]] std::vector<CompactGenome> LoadCompactGenomesJson(std::string_view path);

void StoreDAGToProtobuf(const DAG& dag, std::string_view reference_sequence,
                        const std::vector<Mutations>& edge_parent_mutations,
                        std::string_view path);
