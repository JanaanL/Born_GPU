project(QTCUBE)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)




  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3")
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
        message(SEND_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

enable_language(CXX)
#find_library(LIBTBB tbb REQUIRED)

FIND_PACKAGE(genericIO NO_MODULE REQUIRED)
INCLUDE_DIRECTORIES(${genericIO_DIRNAME_include})




install(PROGRAMS ${PROJECT_BINARY_DIR}/Sep_cube

 DESTINATION bin
)

add_subdirectory(base)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/base)

ADD_EXECUTABLE(RTM3D main.cc)
TARGET_LINK_LIBRARIES(RTM3D base genericCpp)
#add_subdirectory(gui)




#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/base)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/gui)
#link_directories(${CURRENT_BINARY_DIR}/base ${CURRENT_BINARY_DIR}/gui)



#ADD_EXECUTABLE(qt_cube main.cpp ${QtApp_RCC_SRCS} )
#set_property(TARGET qt_cube PROPERTY LINKER_LANGUAGE CXX)
#TARGET_LINK_LIBRARIES(qt_cube gui)
#TARGET_LINK_LIBRARIES(qt_cube base)
#TARGET_LINK_LIBRARIES(qt_cube genericCpp ${GENERIC_EXPORT_LIBS} ${LIBTBB})

#TARGET_LINK_LIBRARIES(qt_cube ${LIBTBB})

#install(TARGETS qt_cube

# DESTINATION bin
#)
