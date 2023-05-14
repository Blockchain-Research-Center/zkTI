include_directories(.)

aux_source_directory(. DIR_SRCS)

add_executable(
        zkTI ${DIR_SRCS}
)

target_link_libraries(
        zkTI snark
)
target_include_directories( 
        zkTI PUBLIC include ${DEPENDS_DIR}/libsnark ${DEPENDS_DIR}/libsnark/depends/libfqfft
)

# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -march=native -mavx2 -faligned-new -lstdc++ -DNDEBUG")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -ggdb -O0 -march=native -mavx2 -faligned-new -lstdc++")