# set(XT_INCLUDES "/opt/homebrew/Cellar/xtensor/0.25.0/include")
# add_executable (test test.cpp)
# target_include_directories(test PRIVATE ${XT_INCLUDES} ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
# target_link_libraries(test ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

add_executable (estimate_error estimate_error.cpp)
target_include_directories(estimate_error PRIVATE ${XT_INCLUDES} ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(estimate_error ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

# add_executable (test_weight_propagation test_weight_propagation.cpp)
# target_include_directories(test_weight_propagation PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
# target_link_libraries(test_weight_propagation ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

add_executable (test_refactor_2 test_refactor_2.cpp)
target_include_directories(test_refactor_2 PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_refactor_2 ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

add_executable (test_ge test_ge.cpp)
target_include_directories(test_ge PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_ge ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

# add_executable (test_encoder test_encoder.cpp)
# target_include_directories(test_encoder PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
# target_link_libraries(test_encoder ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

# add_executable (test_decomposer test_decomposer.cpp)
# target_include_directories(test_decomposer PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
# target_link_libraries(test_decomposer ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

# add_executable (test_interleaver test_interleaver.cpp)
# target_include_directories(test_interleaver PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
# target_link_libraries(test_interleaver ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

# add_executable (test_error_collector test_error_collector.cpp)
# target_include_directories(test_error_collector PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
# target_link_libraries(test_error_collector ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

add_executable (test_refactor test_refactor.cpp)
target_include_directories(test_refactor PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_refactor ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

add_executable (test_reconstructor test_reconstructor.cpp)
target_include_directories(test_reconstructor PRIVATE ${EVA_INCLUDES} ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_reconstructor ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})
