#pragma once

#include "larch/mat_conversion.hpp"

inline void fill_static_reference_sequence(std::string_view dag_ref);

template <typename DAG>
void check_edge_mutations(DAG dag);

template <typename DAG, typename RadiusCallback>
auto optimize_dag_direct(DAG dag, Move_Found_Callback& callback,
                                 RadiusCallback&& radius_callback);

#include "larch/impl/usher_glue_impl.hpp"
#include "larch/impl/produce_mat_impl.hpp"