set(TWOB_SOURCES energy2b.cpp 
                 x2b-v9x.cpp 
                 poly-2b-v6x-nogrd.cpp 
                 poly-2b-v6x.cpp 
# BEGIN SECTION CMAKELISTS
                 poly_2b_A1B2_A1B2_v1x.cpp 
                 poly_2b_A1B2_A1B2_v1.cpp
                 x2b_A1B2_A1B2_v1x.cpp 
                 x2b_A1B2Z2_C1D2_v1x.cpp 
                 poly_2b_A1B2Z2_C1D2_v1x.cpp 
                 poly_2b_A1B2Z2_C1D2_v1.cpp
# END SECTION CMAKELISTS
)

#SET_SOURCE_FILES_PROPERTIES(poly-2b-v6x.cpp PROPERTIES COMPILE_FLAGS -O0)

add_library(2b OBJECT ${TWOB_SOURCES})
target_include_directories(2b PRIVATE ${CMAKE_SOURCE_DIR}) 
