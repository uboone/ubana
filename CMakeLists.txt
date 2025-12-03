# ======================================================================
#  ubana main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/lardataobj/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules REQUIRED)
project(ubana VERSION 10.04.07.19 LANGUAGES CXX)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  ALLOW_DEPRECATIONS
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-expansion-to-defined 
  #-Wno-error=unused-variable 
)


cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

cet_cmake_module_directories(Modules BINARY)

find_package( ubreco REQUIRED EXPORT )
find_package( ubcv REQUIRED EXPORT )
find_package( larana REQUIRED EXPORT )
find_package( larpandora REQUIRED EXPORT )
find_package( swtrigger REQUIRED EXPORT )
find_package(Eigen3 REQUIRED)

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(ubana)

# Cmake modules.
add_subdirectory(Modules)

# Unit tests.
add_subdirectory(test)

# ups - table and config files

# packaging utility
cet_cmake_config()
