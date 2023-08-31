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

swtrigger::Swtrigger_Base - libSwtrigger_Base.so
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

    # Internal transitive dependencies.

    set(_swtrigger_tdep_FEMBeamTrigger "SWTriggerBase")

    # Loop over libraries.

    foreach(_swtrigger_lib_name IN ITEMS SWTriggerBase FEMBeamTrigger )
      if(NOT TARGET swtrigger::${_swtrigger_lib_name})

        # Hunt for this library.

        find_library(_swtrigger_lib_path LIBRARY NAMES ${_swtrigger_lib_name} HINTS ENV SWTRIGGER_LIBDIR REQUIRED NO_CACHE)
        message("Found swtrigger library ${_swtrigger_lib_path}")

        # Make target.

        message("Making target swtrigger::${_swtrigger_lib_name}")
        add_library(swtrigger::${_swtrigger_lib_name} SHARED IMPORTED)

        # Calculate internal transitive dependencies.

        set(_swtrigger_tdep)
        if(_swtrigger_tdep_${_swtrigger_lib_name})
          set(_swtrigger_tdep "swtrigger::${_swtrigger_tdep_${_swtrigger_lib_name}}")
        endif()

        set_target_properties(swtrigger::${_swtrigger_lib_name} PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${_swtrigger_include_dir}"
          IMPORTED_LOCATION "${_swtrigger_lib_path}"
          INTERFACE_LINK_LIBRARIES "${_swtrigger_tdep}"
        )
      endif()

      # End of loop over libraries.

      unset(_swtrigger_lib_path)
    endforeach()
  endif()
endif()
