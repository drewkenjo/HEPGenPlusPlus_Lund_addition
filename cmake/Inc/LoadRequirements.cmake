message(STATUS "Looking for dependencies..")

#
# Find ROOT
#
find_package (ROOT REQUIRED)
include(${ROOT_USE_FILE})
include_directories (${ROOT_INCLUDE_DIRS})

if(NOT ROOT_mathmore_FOUND)
     message(FATAL_ERROR "ROOT component mathmore not found")
endif()

#
# Find GSL
#
#find_package(GSL REQUIRED)
#include_directories (${GSL_INCLUDE_DIRS})

#
# Find GEANT
#
find_package(Geant4 REQUIRED ui_all vis_all)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

#
# Curses
#
find_package(Curses)
