#pragma once

#include "larch/mat_conversion.hpp"

template <typename DAG>
void check_edge_mutations(DAG dag);

template <typename DAG, typename RadiusCallback, typename ReassignCallback>
auto optimize_dag_direct(DAG dag, Move_Found_Callback& callback,
                         RadiusCallback&& radius_callback,
                         ReassignCallback&& reassign_callback,
                         std::optional<uint32_t> user_seed = std::nullopt);

#include "larch/impl/usher_glue_impl.hpp"
#include "larch/impl/produce_mat_impl.hpp"
