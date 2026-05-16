cmake_minimum_required(VERSION 3.20)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address -fno-omit-frame-pointer -Wframe-larger-than=4096 -Wall -Wextra -pedantic -Werror")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -pedantic -Werror")
set(CMAKE_CXX_FLAGS_RELEASE "-g -Wall -Wextra -Werror -O2")

project(liberqs)

add_library(erqs
  src/quantumstate.cpp
  src/operations.cpp
  src/optimization.cpp
)

target_include_directories(erqs PUBLIC include)
