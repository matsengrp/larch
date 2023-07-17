/**
 * Functions for loading and writing mutation annotated dags.
 */
#pragma once

#include <string_view>
#include <unordered_map>

#include "larch/merge/merge.hpp"
#include "larch/mat_conversion.hpp"

[[nodiscard]] MADAGStorage LoadDAGFromProtobuf(std::string_view path);

[[nodiscard]] MADAGStorage LoadTreeFromProtobuf(std::string_view path,
                                                std::string_view reference_sequence);

/*
[[nodiscard]] MADAGStorage LoadTreeFromProtobuf(std::string_view path,
                                                std::string_view reference_sequence,
                                                std::string_view vcf_path);
*/

[[nodiscard]] MADAGStorage LoadDAGFromJson(std::string_view path);

[[nodiscard]] std::string LoadReferenceSequence(std::string_view path);

template <typename DAG>
void StoreDAGToProtobuf(DAG dag, std::string_view path);

template <typename DAG>
void StoreTreeToProtobuf(DAG dag, std::string_view path);

template <typename DAG>
void MADAGToDOT(DAG dag, std::ostream& out);

template <typename DAG>
void FragmentToDOT(DAG dag, const std::vector<EdgeId>& edges, std::ostream& out);

std::string ToEdgeMutationsString(const MAT::Node* node);

void MATToDOT(const MAT::Tree& mat, std::ostream& out);

#include "larch/impl/dag_loader_impl.hpp"
