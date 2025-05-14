cmake_minimum_required(VERSION 3.10)
project(SimpleImageProcessor)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(SimpleImage SHARED Image.cpp)

target_include_directories(SimpleImage PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(demo demo.cpp)
target_link_libraries(demo SimpleImage)
