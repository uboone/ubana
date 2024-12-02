#[================================================================[.rst:

Name:    Findswtrigger.cmake

Purpose: find_package module for ups product swtrigger.

Created: 31-Aug-2023  H. Greenlee

------------------------------------------------------------------

The swtrigger ups product defines the following environment variables,
which are used by this module.

SWTRIGGER_INCDIR - Include path
SWTRIGGER_LIBDIR - Library path

This module creates the following targets, corresponding to the libraries
in the library directory.

swtrigger::SWTriggerBase - libSWTriggerBase.so
swtrigger::FEMBeamTrigger - libFEMBeamTrigger.so


#]================================================================]

# Don't do anything of this package has already been found.

if(NOT swtrigger_FOUND)

  # First hunt for the swtrigger include directory.

  message("Finding package swtrigger")
  find_file(_swtrigger_h NAMES SWTriggerBase HINTS ENV SWTRIGGER_INCDIR NO_CACHE)
  if(_swtrigger_h)
    get_filename_component(_swtrigger_include_dir ${_swtrigger_h} DIRECTORY)
    message("Found swtrigger include directory ${_swtrigger_include_dir}")
    set(swtrigger_FOUND TRUE)
  else()
    message("Could not find swtrigger include directory")
  endif()

  # Next hunt for the swtrigger libraries.

  if(swtrigger_FOUND)

    if(NOT TARGET swtrigger::SWTriggerBase)

      # Hunt for this library.

      find_library(_swtrigger_lib_path LIBRARY NAMES SWTriggerBase HINTS ENV SWTRIGGER_LIBDIR REQUIRED NO_CACHE)
      message("Found swtrigger library ${_swtrigger_lib_path}")

      # Make target.

      message("Making target swtrigger::SWTriggerBase")
      add_library(swtrigger::SWTriggerBase SHARED IMPORTED)
      set_target_properties(swtrigger::SWTriggerBase PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${_swtrigger_include_dir}"
        IMPORTED_LOCATION "${_swtrigger_lib_path}"
      )
      unset(_swtrigger_lib_path)
    endif()

    if(NOT TARGET swtrigger::FEMBeamTrigger)

      # Hunt for this library.

      find_library(_swtrigger_lib_path LIBRARY NAMES FEMBeamTrigger HINTS ENV SWTRIGGER_LIBDIR REQUIRED NO_CACHE)
      message("Found swtrigger library ${_swtrigger_lib_path}")

      # Make target.

      message("Making target swtrigger::FEMBeamTrigger")
      add_library(swtrigger::FEMBeamTrigger SHARED IMPORTED)
      set_target_properties(swtrigger::FEMBeamTrigger PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${_swtrigger_include_dir}"
        IMPORTED_LOCATION "${_swtrigger_lib_path}"
        INTERFACE_LINK_LIBRARIES "swtrigger::SWTriggerBase"
      )
      unset(_swtrigger_lib_path)
    endif()
  endif()
endif()
