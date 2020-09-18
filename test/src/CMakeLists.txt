
add_executable(compress_dictionary compress_dictionary.cpp)
target_link_libraries(compress_dictionary common sdsl malloc_count)

FetchContent_GetProperties(r-index)
FetchContent_GetProperties(shaped_slp)

set(FOLCA_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/folca)
set(SUX_SOURCE_DIR ${shaped_slp_SOURCE_DIR}/external/sux/sux)


add_executable(ms matching_statistics.cpp)
target_link_libraries(ms common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(ms PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(ms PUBLIC "-std=c++17")

add_executable(shapedslp_test shapedslp_test.cpp)
target_link_libraries(shapedslp_test common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(shapedslp_test PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(shapedslp_test PUBLIC "-std=c++17")