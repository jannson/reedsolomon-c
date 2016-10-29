cmake_minimum_required (VERSION 2.6)
project (examples)
#set(CMAKE_VERBOSE_MAKEFILE ON)
#add_definitions( -O0 -g -rdynamic -fno-inline -Wstrict-prototypes -std=gnu99 -Wall)
add_definitions( -Os -Wstrict-prototypes -Wstrict-prototypes -std=gnu99 -Wall)

set(CMAKE_SYSTEM_NAME Linux)

SET(CMAKE_BUILD_TYPE debug)
#SET(CMAKE_EXE_LINKER_FLAGS "-static")

include_directories(${PROJECT_SOURCE_DIR}/src)

add_executable(test_rs ../test_rs.c)
target_link_libraries(test_rs)

#add_executable(fec_test fec_test.c ../rs.c)
#target_link_libraries(fec_test)

add_executable(simple-encoder simple-encoder.c ../rs.c)
target_link_libraries(simple-encoder)

