
add_library(usher-sampled-lib
src/usher-sampled/fix_par_mut.cpp
src/usher-sampled/sampler.cpp
src/usher-sampled/static_tree_mapper/build_idx.cpp
src/usher-sampled/static_tree_mapper/search.cpp
src/usher-sampled/Min_back_FS.cpp
src/usher-sampled/import_vcf.cpp
src/usher-sampled/wait_debug.cpp
src/usher-sampled/place_sample_follower.cpp
src/usher-sampled/place_sample.cpp
src/usher-sampled/main_mapper.cpp
src/usher-sampled/testers.cpp
src/usher-sampled/utils.cpp
src/usher-sampled/place_sample_shared.cpp
src/usher-sampled/multiple_placement.cpp
src/usher-sampled/usher_only_utils.cpp
src/usher-sampled/static_tree_mapper/build_idx.cpp
src/usher-sampled/static_tree_mapper/search.cpp
src/matOptimize/main_helper.cpp
src/matOptimize/detailed_mutations_store.cpp
src/matOptimize/detailed_mutations_load.cpp
src/matOptimize/mutation_annotated_tree.cpp
src/matOptimize/mutation_annotated_tree_node.cpp
src/matOptimize/mutation_annotated_tree_load_store.cpp
src/matOptimize/mutation_annotated_tree_nuc_util.cpp
src/matOptimize/Mutation_Collection.cpp
src/matOptimize/check_samples.cpp
src/matOptimize/optimize_tree.cpp
src/matOptimize/condense.cpp
src/matOptimize/Fitch_Sankoff.cpp
src/matOptimize/main_load_tree.cpp
src/matOptimize/priority_conflict_resolver.cpp
src/matOptimize/apply_move/move_node_no_reassign.cpp
src/matOptimize/apply_move/forward_pass.cpp
src/matOptimize/apply_move/one_level_fitch_sankoff.cpp
src/matOptimize/apply_move/one_level_fitch_sankoff_binary.cpp
src/matOptimize/apply_move/backward_pass.cpp
src/matOptimize/apply_move/apply_move.cpp
src/matOptimize/apply_move/debug.cpp
src/matOptimize/Profitable_Moves_Enumerators/process_terminal_nodes.cpp
src/matOptimize/Profitable_Moves_Enumerators/Profitable_Moves_Enumerators-debug.cpp
src/matOptimize/Profitable_Moves_Enumerators/possible_change.cpp
src/matOptimize/Profitable_Moves_Enumerators/process_LCA.cpp
src/matOptimize/Profitable_Moves_Enumerators/Profitable_Moves_Enumerator.cpp
src/matOptimize/Profitable_Moves_Enumerators/process_inidividual_mutation.cpp
src/matOptimize/Profitable_Moves_Enumerators/downward_integrated.cpp
src/matOptimize/Profitable_Moves_Enumerators/process_intermediate_node.cpp
src/matOptimize/Profitable_Moves_Enumerators/check_move_profitable_dst_branch.cpp
src/matOptimize/Profitable_Moves_Enumerators/upward_integrated.cpp
src/matOptimize/Profitable_Moves_Enumerators/range_tree.cpp
src/matOptimize/Profitable_Moves_Enumerators/check_move_profitable_LCA.cpp
src/matOptimize/Profitable_Moves_Enumerators/check_move_profitable.cpp
${PROTO_SRCS}
${PROTO_HDRS}
${DETAILED_MUTATIONS_PROTO_SRCS}
${DETAILED_MUTATIONS_PROTO_HDRS}
)

add_dependencies(usher-sampled-lib usher-sampled)
target_include_directories(usher-sampled-lib PUBLIC "${PROJECT_BINARY_DIR}")
target_compile_options(usher-sampled-lib PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
target_link_libraries(usher-sampled-lib PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB  ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${ISAL_LIB})

add_library(matOptimize-lib
src/matOptimize/matOptimize-main.cpp
src/matOptimize/mutation_annotated_tree.cpp
src/matOptimize/mutation_annotated_tree_node.cpp
src/matOptimize/mutation_annotated_tree_load_store.cpp
src/matOptimize/detailed_mutations_store.cpp
src/matOptimize/detailed_mutations_load.cpp
src/matOptimize/mutation_annotated_tree_nuc_util.cpp
src/matOptimize/optimize_tree.cpp
src/matOptimize/import_vcf_fast.cpp
src/matOptimize/condense.cpp
src/matOptimize/VCF_load_tree.cpp
src/matOptimize/main_load_tree.cpp
src/matOptimize/main_helper.cpp
src/matOptimize/Mutation_Collection.cpp
src/matOptimize/Fitch_Sankoff.cpp
src/matOptimize/check_samples.cpp
src/matOptimize/priority_conflict_resolver.cpp
src/matOptimize/transpose_vcf/transposed_vcf_patch.cpp
${patch_tree}
${New_Profitable_Moves_Enumerators}
${PROTO_SRCS}
${PROTO_HDRS}
${DETAILED_MUTATIONS_PROTO_SRCS}
${DETAILED_MUTATIONS_PROTO_HDRS}
)

add_dependencies(matOptimize-lib matOptimize)
target_include_directories(matOptimize-lib PUBLIC "${PROJECT_BINARY_DIR}")
target_compile_options(matOptimize-lib PRIVATE -DTBB_SUPPRESS_DEPRECATED_MESSAGES)
target_link_libraries(matOptimize-lib PRIVATE stdc++  ${Boost_LIBRARIES} ${TBB_IMPORTED_TARGETS} ${Protobuf_LIBRARIES} ZLIB::ZLIB  ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} ${ISAL_LIB})
