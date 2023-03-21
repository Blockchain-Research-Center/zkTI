include_directories(.)

add_executable(
        zkTI
                main.cpp   
)
target_link_libraries(
        zkTI
                snark
)
target_include_directories( 
        zkTI
                PUBLIC
                include
                ${DEPENDS_DIR}/libsnark
                ${DEPENDS_DIR}/libsnark/depends/libfqfft
)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -march=native -mavx2 -faligned-new -lstdc++ -DNDEBUG")
