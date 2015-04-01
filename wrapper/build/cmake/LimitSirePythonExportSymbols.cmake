#
# The aim of this file is to make the changes necessary so that 
# *only* the init_PyModule function is exported from each python
# library. This will allow "strip" to remove most of the 
# unnecessary symbols from the library, thereby reducings
# its size.
#

macro ( EXPORT_THIS_SYMBOL_ONLY  _symbol _mangled_symbol )

  if ( CMAKE_COMPILER_IS_GNUCC )
    # use gcc/ld
    #set( CMAKE_SHARED_LINKER_FLAGS 
    #     "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-exported_symbol,${_mangled_symbol}" )

  elseif ( ${CMAKE_CXX_COMPILER} MATCHES "xlC" )
    # we need to create an 'exported_symbols' file that just contains
    # this symbol
    set( _filename "${CMAKE_BINARY_DIR}/${_symbol}_export_file" )
    file( WRITE "${_filename}" "${_symbol}\n" )
    set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -bE:${_filename}" )

  elseif ( ${CMAKE_CXX_COMPILER} MATCHES "icpc" )

  endif()

endmacro()
