# A minimal consumer of the installed obdeect stream package. Built by
# python/tests/test_cmake_package.py with -Dobdeect_stream_DIR=<directory of
# obdeect_stream-config.cmake>.
cmake_minimum_required(VERSION 3.20)
project(obdeect_stream_consumer LANGUAGES C CXX)

find_package(obdeect_stream CONFIG REQUIRED)

add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE obdeect::stream)

add_executable(consumer_c main_c.c)
target_link_libraries(consumer_c PRIVATE obdeect::stream_c)
