#----------------------------------------------------------------------------
# Setup the project
# cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()


project(FieldCorrection C CXX)


add_definitions("-std=c++11")
#----------------------------------------------------------------------------
# Find Geant4 package, activating all available Vis drivers by default
# You can set WITH_GEANT4_VIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
# option(WITH_GEANT4_VIS "Build example with Geant4 Vis drivers" ON)
# if(WITH_GEANT4_VIS)
#   find_package(Geant4 REQUIRED gdml vis_all)
# else()
#   find_package(Geant4 REQUIRED gdml)
# endif()

find_package(CGAL QUIET)
if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
else()
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
# include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
#list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
#list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS)
add_definitions(${ROOT_CXX_FLAGS})
include(${ROOT_USE_FILE})
# find_package(Boost REQUIRED)

#----------------------------------------------------------------------------
# Find package for multi-threading
find_package (Threads)

#----------------------------------------------------------------------------
# Set CRY Loacation
#
# set(CRY_HOME_DIR /opt/cry_v1.7)

#----------------------------------------------------------------------------
# Find CRY libraries
#
# find_library(CRY_LIBRARIES CRY ${CRY_HOME_DIR}/lib REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
#include_directories(${ROOT_INCLUDE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})
# include_directories(${Boost_INCLUDE_DIR})
# file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cpp)
# file(GLOB headers ${PROJECT_SOURCE_DIR}/*.hpp)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(FieldCorrection FieldCorrection.cpp Laser.cpp Laser.hpp LaserTrack.cpp LaserTrack.hpp ThreeVector.hpp Matrix3x3.cpp Matrix3x3.hpp TPCVolumeHandler.cpp TPCVolumeHandler.hpp Interpolation3D.cpp Interpolation3D.hpp)
target_link_libraries(FieldCorrection ${ROOT_LIBRARIES} ${CGAL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build G01. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
#### FOLLOWING LINES COMMENTED OUT BY MRM ####
#
#set(FieldCorrection_Files Field.root)
#
#foreach(_File ${FieldCorrection_Files})
#  configure_file(
#    ${PROJECT_SOURCE_DIR}/${_File}
#    ${PROJECT_BINARY_DIR}/${_File}
#    COPYONLY
#    )
#endforeach()

#----------------------------------------------------------------------------
# Copy GDMLSchema folder

# file(COPY ${PROJECT_SOURCE_DIR}/GDMLSchema DESTINATION ${PROJECT_BINARY_DIR})

# set(GDML)
# make_directory(${PROJECT_BINARY_DIR}/GDMLSchema)

# file(GLOB GDMLSchema RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/*)
# foreach(templateFile ${GDMLSchema})
#     set(GDMLPath ${PROJECT_SOURCE_DIR}/${templateFile})
#     if(NOT IS_DIRECTORY ${GDMLPath})
# 	message(STATUS "Configuring file ${templateFile}")
# 	configure_file(${GDMLPath}
# 		       ${PROJECT_BINARY_DIR}/${templateFile}
#                        @ONLY)
#     endif(NOT IS_DIRECTORY ${GDMLPath})
# endforeach(templateFile)

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
# add_custom_target(gamma-uboone DEPENDS load_gdml)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS FieldCorrection DESTINATION bin)

