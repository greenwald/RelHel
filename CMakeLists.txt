cmake_minimum_required(VERSION 2.8.12)
project(RelHel)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmakeModules)

# Default C++ flags
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if (COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
	message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

CHECK_CXX_COMPILER_FLAG("-fdiagnostics-color=always" COMPILER_SUPPORTS_DIAGNOSTIC_COLORS)
if (COMPILER_SUPPORTS_DIAGNOSTIC_COLORS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wno-long-long -Wno-missing-braces -Werror=overloaded-virtual -pthread" )

# Default build type is 'Release' (cmake default one is '')
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
		FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# build with cmake -DCMAKE_BUILD_TYPE=Debug
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g2" )

# build with cmake -DCMAKE_BUILD_TYPE=Release
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -D ELPP_DISABLE_DEBUG_LOGS" )

# Add subdirectories
add_subdirectory(src)
add_subdirectory(programs)

