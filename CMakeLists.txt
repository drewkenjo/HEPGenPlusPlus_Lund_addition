project (hepgen_cpp)
cmake_minimum_required(VERSION 3.0)


# Add module to cmake
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Inc/  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)

message("-------------------------------------")
message("--- HEPGen++ Package installation ---")
message("-------------------------------------")

if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    add_compile_options(-Wall -Wextra)
else()
    add_compile_options(-w)
endif()

#
# Requirements
#

include(AutoDetectCompiler)
include(LoadRequirements)


#### Check input variables
if (${HEPGEN_STANDALONE}   MATCHES "YES")
message("*** WANT STANDALONE COMPILATION ***")
set (USE_STANDALONE 1)
set (HEPGEN_DISABLE_ROOT "YES")
endif (${HEPGEN_STANDALONE}   MATCHES "YES")

if (${HEPGEN_DISABLE_ROOT} MATCHES "YES")
message("*** WANT DISABLED ROOT")
set (ROOT_FOUND 0)
else (${HEPGEN_DISABLE_ROOT} MATCHES "YES")
find_package (ROOT)
set (USE_ROOT 1)
endif (${HEPGEN_DISABLE_ROOT} MATCHES "YES")

if (${HEPGEN_ENABLE_PYTHON} MATCHES "YES")
message("*** WANT PYTHON BINDINGS")
set (PYTHON_WANTED 1)
else (${HEPGEN_ENABLE_PYTHON} MATCHES "YES")
set (PYTHON_WANTED 0)
endif (${HEPGEN_ENABLE_PYTHON} MATCHES "YES")




if (${HEPGEN_EXPERIMENTAL_BUILD} MATCHES "YES")
message("*** WANT EXPERIMENTAL BUILD")
set (USE_EXPERIMENTAL 1)
endif (${HEPGEN_EXPERIMENTAL_BUILD} MATCHES "YES")

if (${HEPGEN_EVDISGL} MATCHES "YES")
message("*** WANT OPENGL EVENT DISPLAY")
set (USE_GL 1)
endif (${HEPGEN_EVDISGL} MATCHES "YES")

#### Set config.h file input

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/libhepgen/config.h )


#### Set include directories to add
include_directories(hepgen_c)
include_directories(libhepgen)


#### Only can build mosse code or GK Pi0 if root is enabled
include_directories(libVGG)
add_subdirectory(libVGG)

include_directories(libGKPi0)
add_subdirectory(libGKPi0)

#### add opengl subdirectories
if (${USE_GL})
  add_subdirectory(evDisGL)
endif(${USE_GL})


#### add standard subdirectories
add_subdirectory(libhepgen)
add_subdirectory(hepgen_c)
add_subdirectory(tools)

### add python bindings if wanted
if (${PYTHON_WANTED})
  add_subdirectory(share/swig)
endif (${PYTHON_WANTED})


#### search datacards
FILE(GLOB datacards resources/datacards/*.data)


#### install misc stuff
install (FILES "share/beamfile.dat" DESTINATION bin)
install (DIRECTORY share DESTINATION .)
install (DIRECTORY "share/steering/" DESTINATION share/hepgen/share/steering)
install (CODE "EXECUTE_PROCESS (COMMAND echo -- Untar pi0 tables into ${CMAKE_INSTALL_PREFIX}/share/pi0_cache/)")
install (CODE "EXECUTE_PROCESS (COMMAND tar -xzf output.tar.gz WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/share/pi0_cache/ )")
