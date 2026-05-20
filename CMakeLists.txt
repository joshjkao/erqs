cmake_minimum_required(VERSION 4.0)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(erqs LANGUAGES CXX)

add_compile_options(
  -Wall
  -Wextra
  -Wpedantic
  -Wshadow
  -Wnon-virtual-dtor
  -Wold-style-cast
  -Wunused
  -Woverloaded-virtual
  -Wconversion
  -Wsign-conversion
  # -Wnull-dereference
  -Wformat=2
  $<$<NOT:$<CONFIG:Debug>>:-Werror>
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address)
endif()

add_subdirectory(liberqs)

find_package(NLopt)

# executable for benchmarking stuff
add_executable(benchmark benchmark.cpp)
target_link_libraries(benchmark PRIVATE erqs NLopt::nlopt)

# executable for testing stuff
add_executable(main main.cpp)
target_link_libraries(main PRIVATE erqs NLopt::nlopt)

include(FetchContent)
FetchContent_Declare(
  cxxopts
  GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
  GIT_TAG v3.3.1
)
FetchContent_MakeAvailable(cxxopts)

# executable for testing stuff
add_executable(randomtest_inner randomtest_inner.cpp)
target_link_libraries(randomtest_inner PRIVATE erqs NLopt::nlopt cxxopts::cxxopts)

# executable for testing stuff
add_executable(randomtest_simplify randomtest_simplify.cpp)
target_link_libraries(randomtest_simplify PRIVATE erqs NLopt::nlopt cxxopts::cxxopts)

# executable for testing stuff
add_executable(randomtest_expectation randomtest_expectation.cpp)
target_link_libraries(randomtest_expectation PRIVATE erqs NLopt::nlopt cxxopts::cxxopts)

# executable for optimizing
# add_executable(opt opt.cpp)
# target_link_libraries(opt PRIVATE erqs NLopt::nlopt)

# tests
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.17.0
)
FetchContent_MakeAvailable(googletest)
add_executable(tests tests.cpp)
target_link_libraries(tests PRIVATE GTest::gtest_main erqs NLopt::nlopt)
include(GoogleTest)
gtest_discover_tests(tests)
