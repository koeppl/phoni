
add_executable(compress_dictionary compress_dictionary.cpp)
target_link_libraries(compress_dictionary common sdsl malloc_count)

add_executable(rlbwt rlbwt.cpp)
target_link_libraries(rlbwt common sdsl malloc_count)

FetchContent_GetProperties(r-index)
FetchContent_GetProperties(shaped_slp)
FetchContent_GetProperties(ssw)

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

add_executable(rlems rle_matching_statistics.cpp)
target_link_libraries(rlems common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(rlems PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(rlems PUBLIC "-std=c++17")


add_executable(phoni phoni.cpp)
target_link_libraries(phoni common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(phoni PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(phoni PUBLIC "-std=c++17")

add_executable(phoni_compatibility phoni_compatibility_driver.cpp)
target_link_libraries(phoni_compatibility common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(phoni_compatibility PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(phoni_compatibility PUBLIC "-std=c++17")

add_executable(build_phoni build_phoni.cpp)
target_link_libraries(build_phoni common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(build_phoni PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(build_phoni PUBLIC "-std=c++17")


add_executable(bwt2rlbwt bwt2rlbwt.cpp)
target_compile_options(bwt2rlbwt PUBLIC "-std=c++17")

add_executable(patternstats patternstats.cpp)
target_compile_options(patternstats PUBLIC "-std=c++17")

add_executable(rlebwt_ms_build rlebwt_ms_build.cpp)
target_link_libraries(rlebwt_ms_build common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(rlebwt_ms_build PUBLIC    "../../include/ms" 
                                        "../../include/common" 
                                        "${shaped_slp_SOURCE_DIR}" 
                                        "${FOLCA_SOURCE_DIR}" 
                                        "${SUX_SOURCE_DIR}/function" 
                                        "${SUX_SOURCE_DIR}/support"
                                        )
target_compile_options(rlebwt_ms_build PUBLIC "-std=c++17")




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

add_executable(align align.cpp)
target_link_libraries(align common sdsl divsufsort divsufsort64 malloc_count ri ssw)
target_include_directories(align PUBLIC    "../../include/ms" 
                                            "../../include/common" 
                                            "${shaped_slp_SOURCE_DIR}" 
                                            "${FOLCA_SOURCE_DIR}" 
                                            "${SUX_SOURCE_DIR}/function" 
                                            "${SUX_SOURCE_DIR}/support"
                                            "${ssw_SOURCE_DIR}/src"
                                        )
target_compile_options(align PUBLIC "-std=c++17")
