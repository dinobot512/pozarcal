cmake_minimum_required(VERSION 3.16)
project(fire_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(fire_sim
    main.cpp
)

target_include_directories(fire_sim PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(fire_sim PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -O2>
    $<$<CXX_COMPILER_ID:MSVC>:/W4 /O2>
)
