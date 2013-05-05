# ---------------------------- Checking the OS ------------------------------
IF("${CMAKE_SYSTEM}" MATCHES "Linux")
	MESSAGE(STATUS "Compiling on Linux")
	SET(LINUX 1)
	EXEC_PROGRAM("uname" ARGS -m OUTPUT_VARIABLE myArch)
	IF("${myArch}" MATCHES "x86_64")
		MESSAGE(STATUS "  64-bit architecture detected")
		SET(LINUX64 1)
	ENDIF("${myArch}" MATCHES "x86_64")
ELSEIF(APPLE)
	SET(LINUX 0)
	MESSAGE(STATUS "Compiling on MacOSX")
ELSEIF(WIN32)
	SET(LINUX 0)
	IF(CYGWIN)
		MESSAGE(STATUS "Compiling on CygWin")
		SET(CYGWIN 1)
	ELSE(CYGWIN)
		MESSAGE(FATAL_ERROR "Compiling on Windows is not possible at the moment")
	ENDIF(CYGWIN)
ENDIF("${CMAKE_SYSTEM}" MATCHES "Linux")

include(ExternalProject)

# Set variables for external libraries
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
SET(PROJECT_EXTERNAL_DIR ${PROJECT_SOURCE_DIR}/external)
INCLUDE_DIRECTORIES(${PROJECT_EXTERNAL_DIR})


# Library: FLIRT 
set(FLIRT_INSTALL_DIR ${PROJECT_EXTERNAL_DIR}/flirtlib)
ExternalProject_Add(flirt
	PREFIX ${FLIRT_INSTALL_DIR}
	SVN_REPOSITORY https://svn.openslam.org/data/svn/flirtlib/trunk
	SVN_TRUST_CERT 1
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${FLIRT_INSTALL_DIR} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DBUILD_DOC=OFF -DBUILD_GUI=OFF
	INSTALL_DIR ${FLIRT_INSTALL_DIR}
)
INCLUDE_DIRECTORIES(${FLIRT_INSTALL_DIR}/include/flirtlib)
LINK_DIRECTORIES(${FLIRT_INSTALL_DIR}/lib/flirtlib)
install(DIRECTORY ${FLIRT_INSTALL_DIR}/share/flirtlib/data DESTINATION share/${PROJECT_NAME} PATTERN ".svn" EXCLUDE)

# Library: Boost
set(Boost_USE_MULTITHREADED ON)
FIND_PACKAGE(Boost 1.36.0 COMPONENTS program_options serialization)
#FIND_PACKAGE(Boost COMPONENTS program_options math graph)
IF(Boost_FOUND)
	MESSAGE(STATUS "Boost found")
	IF($ENV{VERBOSE})
 		MESSAGE(STATUS "   Boost_INCLUDE_DIR = ${Boost_INCLUDE_DIR}")
 		MESSAGE(STATUS "   Boost_LINK_DIRECTORIES = ${Boost_LINK_DIRECTORIES}")
		MESSAGE(STATUS "   Boost Version = m:${Boost_MINOR_VERSION} M:${Boost_MAJOR_VERSION} v:${Boost_VERSION}")
	ENDIF($ENV{VERBOSE})
	INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
	LINK_DIRECTORIES(${Boost_LINK_DIRECTORIES})
ELSE(Boost_FOUND)
	MESSAGE(FATAL_ERROR "	Boost not found, it is REQUIRED to build the FLIRT.\n Boost can be found on http://www.boost.org")
ENDIF(Boost_FOUND)

# Library: Eigen3
INCLUDE(./build_tools/FindEigen3.cmake REQUIRED)
IF(EIGEN3_FOUND)
        MESSAGE(STATUS "Eigen3 found")
        IF($ENV{VERBOSE})
                MESSAGE(STATUS "   EIGEN3_INCLUDE_DIR = ${EIGEN3_INCLUDE_DIR}")
        ENDIF($ENV{VERBOSE})
        INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
ELSE(EIGEN3_FOUND)
        MESSAGE(FATAL_ERROR "   Eigen3 not found, but it is REQUIRED to build FLIRT.")
ENDIF(EIGEN3_FOUND)

# Optimization
SET(LOCAL_OPTIMIZATION_FLAGS "-O3 -march=nocona -msse4.1 -funroll-loops -ftree-vectorize -DTBB -mfpmath=both")
# SET(LOCAL_OPTIMIZATION_FLAGS "-O3 -march=native -msse4.2 -msse4.1 -funroll-loops -ftree-vectorize -DTBB -mfpmath=both")
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${LOCAL_OPTIMIZATION_FLAGS}")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=gnu++0x ${LOCAL_OPTIMIZATION_FLAGS}")

#Debug (it can clash with optimization)
SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g3 -ggdb3 -fno-inline -Wall")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=gnu++0x -O0 -g3 -ggdb3 -fno-inline -Wall")

#Profiling
#ADD_DEFINITIONS(-pg)
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")


# RPATH Handling during building and installation phases
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(CMAKE_INSTALL_RPATH "${FLIRT_INSTALL_DIR}/lib/flirtlib:$ORIGIN/../lib/${CMAKE_PROJECT_NAME}")
# 
# # use, i.e. don't skip the full RPATH for the build tree
# SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# # when building, don't use the install RPATH already
# # (but later on when installing)
# SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
# # set(INSTALL_RPATH_USE_LINK_PATH true)
# # Set the install RPATH relative to the binaries

# # add the automatically determined parts of the RPATH
# # which point to directories outside the build tree to the install RPATH
# SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# # the RPATH to be used when installing, but only if it's not a system directory
# # LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_PROJECT_NAME}" isSystemDir)
# # IF("${isSystemDir}" STREQUAL "-1")
#  #  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_PROJECT_NAME}")
# # ENDIF("${isSystemDir}" STREQUAL "-1")


add_subdirectory(gflip)
ADD_SUBDIRECTORY(vocabulary)
ADD_SUBDIRECTORY(applications)

message(STATUS "INSTALL_RPATH: ${INSTALL_RPATH}")
message(STATUS "CMAKE_INSTALL_RPATH: ${CMAKE_INSTALL_RPATH}")


