cmake_minimum_required(VERSION 3.22)

project(vMF)

add_library(vMF INTERFACE)
target_include_directories(vMF INTERFACE include)
target_compile_features(vMF INTERFACE cxx_std_17)

add_executable(
        test-vMF
        test-vMF.cpp
        )
target_link_libraries(test-vMF vMF)

target_include_directories(
        test-vMF PUBLIC
        /home/tin/git-repos/cpp-soraxas_toolbox/include/
)