find_package(MPI REQUIRED)

add_executable(refactorGE refactorGE.cpp)
target_include_directories(refactorGE PRIVATE ${MPI_INCLUDE_PATH} ${MGARDx_INCLUDES})
target_link_libraries(refactorGE ${PROJECT_NAME} ${ZSTD_LIB} MPI::MPI_CXX adios2::adios2)

add_executable(refactorGE_SZ3 refactorGE_SZ3.cpp)
target_include_directories(refactorGE_SZ3 PRIVATE ${MPI_INCLUDE_PATH} ${MGARDx_INCLUDES})
target_link_libraries(refactorGE_SZ3 ${PROJECT_NAME} ${ZSTD_LIB} MPI::MPI_CXX adios2::adios2)

add_executable(refactorGE_SZ3_delta refactorGE_SZ3_delta.cpp)
target_include_directories(refactorGE_SZ3_delta PRIVATE ${MPI_INCLUDE_PATH} ${MGARDx_INCLUDES})
target_link_libraries(refactorGE_SZ3_delta ${PROJECT_NAME} ${ZSTD_LIB} MPI::MPI_CXX adios2::adios2)

add_executable(testGE_VTOT testGE_VTOT.cpp)
target_include_directories(testGE_VTOT PRIVATE ${MPI_INCLUDE_PATH} ${MGARDx_INCLUDES})
target_link_libraries(testGE_VTOT ${PROJECT_NAME} ${ZSTD_LIB} MPI::MPI_CXX adios2::adios2)

add_executable(testGE_VTOT_SZ3 testGE_VTOT_SZ3.cpp)
target_include_directories(testGE_VTOT_SZ3 PRIVATE ${MPI_INCLUDE_PATH} ${MGARDx_INCLUDES})
target_link_libraries(testGE_VTOT_SZ3 ${PROJECT_NAME} ${ZSTD_LIB} MPI::MPI_CXX adios2::adios2)

add_executable(testGE_VTOT_SZ3_delta testGE_VTOT_SZ3_delta.cpp)
target_include_directories(testGE_VTOT_SZ3_delta PRIVATE ${MPI_INCLUDE_PATH} ${MGARDx_INCLUDES})
target_link_libraries(testGE_VTOT_SZ3_delta ${PROJECT_NAME} ${ZSTD_LIB} MPI::MPI_CXX adios2::adios2)