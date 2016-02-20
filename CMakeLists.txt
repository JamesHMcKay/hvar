cmake_minimum_required(VERSION 2.8)
project(main)
message(status "${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(src include)
set(HEADER_FILES include/boost_offset.hpp include/get_data.hpp include/hvar.hpp include/Figures.hpp)
set(FIGURES src/Figure_1_a.cpp src/Figure_1_b.cpp src/Figure_2.cpp src/Figure_3.cpp src/Figure_4_6.cpp src/Figure_7.cpp src/Figure_8.cpp src/Figure_9.cpp src/Figure_13.cpp src/Figure_14.cpp src/Figure_15.cpp src/Figure_23.cpp)



set(Boost_INCLUDE_DIR /usr/local/Cellar/boost/1.59.0/include)
set(Boost_LIBRARY_DIR /usr/local/Cellar/boost/1.59.0/lib)
find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})



add_executable(main src/main.cpp src/boost_offset.cpp src/hvar.cpp src/get_data.cpp src/scanner.cpp ${HEADER_FILES})

add_executable(Figures src/Figures.cpp src/boost_offset.cpp src/hvar.cpp src/get_data.cpp src/scanner.cpp ${FIGURES} ${HEADER_FILES})

#add_executable(gamma_test src/gamma_test.cpp ${HEADER_FILES})

target_link_libraries( main ${Boost_LIBRARIES} )


if (CMAKE_COMPILER_IS_GNUCXX)
    message(status "Setting GCC flags")

    # Remove the -g flag when not debugging
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -Wall -Wextra -W -pedantic -std=c++0x -O3")
else()
    message(status "Setting MSVC flags")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHc-")
endif()
install(TARGETS main DESTINATION ${CMAKE_BINARY_DIR}/lib)