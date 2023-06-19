include_directories(.)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mavx2 -faligned-new -lstdc++ -DNDEBUG")

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

message(STATUS "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
message(STATUS "OPT_FLAGS = ${OPT_FLAGS}")
