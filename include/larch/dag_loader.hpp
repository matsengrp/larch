/**
 * Functions for loading and writing mutation annotated dags.
 */
#pragma once

#include <string_view>
#include <unordered_map>

#include "larch/merge/merge.hpp"

[[nodiscard]] MADAGStorage LoadDAGFromProtobuf(std::string_view path);

[[nodiscard]] MADAGStorage LoadTreeFromProtobuf(std::string_view path,
                                                std::string_view reference_sequence);

[[nodiscard]] MADAGStorage LoadDAGFromJson(std::string_view path);

[[nodiscard]] std::string LoadReferenceSequence(std::string_view path);

template <typename DAG>
void StoreDAGToProtobuf(DAG dag, std::string_view path);

void StoreTreeToProtobuf(MADAG dag, std::string_view path);

void MADAGToDOT(MADAG dag, std::ostream& out);

#include "larch/impl/dag_loader_impl.hpp"