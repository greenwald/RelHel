set(INCLUDE_DIRECTORIES
	${RelHel_SOURCE_DIR}/include)
include_directories(${INCLUDE_DIRECTORIES})

set(RelHel_SOURCES
	ClebschGordanBox.cc
	MathUtils.cc
	TFhh.cc
	TFracNum.cc
	TJSS.cc
	TJwfTensor.cc
	TLSAmpl.cc
	TLSContrib.cc
	TLSNonRel.cc
	TSpinWaveFunction.cc
	generatePrimeNumberCache.cc
	primeNumbers.cc
	)

add_library(RelHel SHARED ${RelHel_SOURCES})

# install destinations can be passed via the command line:
# cmake -DLIBRARY_OUTPUT_DIRECTORY:PATH=<lib_path> <path_to_CMakeLists.tex>
# otherwise, default LD_LIBRARY_PATH
if(NOT DEFINED LIBRARY_OUTPUT_DIRECTORY)
	set(LIBRARY_OUTPUT_DIRECTORY ${RelHel_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE})
endif()

if(NOT DEFINED INCLUDE_OUTPUT_DIRECTORY)
	set(INCLUDE_OUTPUT_DIRECTORY ${RelHel_SOURCE_DIR}/include/RelHel)
endif()

install(TARGETS RelHel LIBRARY DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})

# Matches all the headers in ${RelHelDID}/include and its subdirs
file(GLOB_RECURSE
	 INSTALL_INCLUDES ${RelHel_SOURCE_DIR}/include/*.h)

# message(STATUS "${INSTALL_INCLUDES}")

# install(FILE       ${INSTALL_INCLUDES}
#	    DESTINATION ${INCLUDE_OUTPUT_DIRECTORY}
#	)
