cmake_minimum_required(VERSION 3.16)

project(ACSP VERSION 0.2.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#if (CMAKE_COMPILER_IS_GNUCXX)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")
#elseif (MSVC)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  /std:c++17")
#else ()
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++17 -O3")
#endif (CMAKE_COMPILER_IS_GNUCXX)


set(PROJECT_TOP_DIR ${PROJECT_SOURCE_DIR})
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/bin)


#add_library(ACSP INTERFACE)
#target_include_directories(ACSP PUBLIC src)

option(ENABLE_UNIT_TESTS "Enable unit tests" ON)
option(ENABLE_DEVELOPMENT "Enable development" ON)

if(${ENABLE_UNIT_TESTS})
    add_definitions(-DENABLE_UNIT_TESTS)
    add_subdirectory(test)
endif()

if(${ENABLE_DEVELOPMENT})
    add_executable(develop devlop.cpp)
    target_include_directories(develop PUBLIC src)
endif()


