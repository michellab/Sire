########################################
#
# CMake file for library:SireMaths
#
########################################

if (SIRE_DISABLE_FORTRAN)
  add_definitions( "-DSIRE_DISABLE_FORTRAN" )
  message( STATUS "Disabling SireMaths BLAS/LAPACK/LINPACK interfaces.")
  message( STATUS "Have to do this as Fortran code has been disabled.")

else()
  # We first create a library (SireBLASPACK) that provides an interface
  # to BLAS, LAPACK and LINPACK. This library either links
  # to those libraries, or it provides the necessary functions
  # from those libraries if they are not available.

  set( SIRE_BLAS_FUNCTIONS
       DGEMM
       DGEMV
     )

  set( SIRE_LAPACK_FUNCTIONS
       DSYEV
     )

  set( SIRE_LINPACK_FUNCTIONS
       DGECO
       DGEDI
     )

  if ( SIRE_USE_BLAS )
    list(APPEND SIREBLASPACK_LIBRARIES ${BLAS_LIBRARIES})
    set( SIREBLASPACK_HAS_LIBRARIES TRUE )
  else()
    set( SIRE_BLAS_SOURCES
         ThirdParty/BLAS/idamax.f
         ThirdParty/BLAS/dasum.f
         ThirdParty/BLAS/daxpy.f
         ThirdParty/BLAS/dcopy.f
         ThirdParty/BLAS/ddot.f
         ThirdParty/BLAS/dgemm.f
         ThirdParty/BLAS/dgemv.f
         ThirdParty/BLAS/dger.f
         ThirdParty/BLAS/dnrm2.f
         ThirdParty/BLAS/dscal.f
         ThirdParty/BLAS/dswap.f
         ThirdParty/BLAS/dsymv.f
         ThirdParty/BLAS/dsyr2.f
         ThirdParty/BLAS/dsyr2k.f
         ThirdParty/BLAS/dtrmm.f
         ThirdParty/BLAS/dtrmv.f
         ThirdParty/BLAS/lsame.f
         ThirdParty/BLAS/xerbla.f
       )
  endif()

  if ( SIRE_USE_LINPACK )
    list(APPEND SIREBLASPACK_LIBRARIES ${LINPACK_LIBRARIES})
    set( SIREBLASPACK_HAS_LIBRARIES TRUE )
  else()
    set( SIRE_LINPACK_SOURCES
         ThirdParty/LINPACK/dgeco.f
         ThirdParty/LINPACK/dgedi.f
         ThirdParty/LINPACK/dgefa.f
       )
  endif()

  if ( SIRE_USE_LAPACK )
    list(APPEND SIREBLASPACK_LIBRARIES ${LAPACK_LIBRARIES})
    set( SIREBLASPACK_HAS_LIBRARIES TRUE )
  else()
    set( SIRE_LAPACK_SOURCES
         ThirdParty/LAPACK/disnan.f
         ThirdParty/LAPACK/dlae2.f
         ThirdParty/LAPACK/dlaev2.f
         ThirdParty/LAPACK/dlaisnan.f
         ThirdParty/LAPACK/dlamch.f
         ThirdParty/LAPACK/dlanst.f
         ThirdParty/LAPACK/dlansy.f
         ThirdParty/LAPACK/dlapy2.f
         ThirdParty/LAPACK/dlascl.f
         ThirdParty/LAPACK/dlarf.f
         ThirdParty/LAPACK/dlarfb.f
         ThirdParty/LAPACK/dlartg.f
         ThirdParty/LAPACK/dlarfg.f
         ThirdParty/LAPACK/dlarft.f
         ThirdParty/LAPACK/dlaset.f
         ThirdParty/LAPACK/dlasr.f
         ThirdParty/LAPACK/dlasrt.f
         ThirdParty/LAPACK/dlassq.f
         ThirdParty/LAPACK/dlatrd.f
         ThirdParty/LAPACK/dorg2l.f
         ThirdParty/LAPACK/dorg2r.f
         ThirdParty/LAPACK/dorgql.f
         ThirdParty/LAPACK/dorgqr.f
         ThirdParty/LAPACK/dorgtr.f
         ThirdParty/LAPACK/dsyev.f
         ThirdParty/LAPACK/dsytd2.f
         ThirdParty/LAPACK/dsytrd.f
         ThirdParty/LAPACK/dsteqr.f
         ThirdParty/LAPACK/dsterf.f
         ThirdParty/LAPACK/ilaenv.f
         ThirdParty/LAPACK/iladlc.f
         ThirdParty/LAPACK/iladlr.f
         ThirdParty/LAPACK/ieeeck.f
         ThirdParty/LAPACK/iparmq.f
       )
  endif()

  if (NOT CMAKE_GENERATOR STREQUAL "Xcode" )
    # We need to sort out the naming convention for the
    # fortran blas, lapack and linpack libraries
    if ( USING_CMAKE_2_6 )
      message( STATUS "Detecting C++/Fortran interface using cmake-2.6 code..." )
      include( FortranCInterface )
      create_fortran_c_interface( Sire SIRE_BLAS_FUNCTIONS sire_blas_f.h )
      create_fortran_c_interface( Sire SIRE_LAPACK_FUNCTIONS sire_lapack_f.h )
      create_fortran_c_interface( Sire SIRE_LINPACK_FUNCTIONS sire_linpack_f.h )
    else()
      message( STATUS "Detecting C++/Fortran interface using cmake-2.8 code..." )
      include( FortranCInterface )

      FortranCInterface_VERIFY( CXX QUIET )

      if (FortranCInterface_VERIFIED_CXX AND FortranCInterface_GLOBAL_FOUND)
        FortranCInterface_HEADER( sire_blas_f.h 
                                  MACRO_NAMESPACE "Sire"
                                  SYMBOL_NAMESPACE "Sire"
                                  SYMBOLS ${SIRE_BLAS_FUNCTIONS} )

        FortranCInterface_HEADER( sire_lapack_f.h
                                  MACRO_NAMESPACE "Sire"
                                  SYMBOL_NAMESPACE "Sire"
                                  SYMBOLS ${SIRE_LAPACK_FUNCTIONS} )

        FortranCInterface_HEADER( sire_linpack_f.h
                                  MACRO_NAMESPACE "Sire"
                                  SYMBOL_NAMESPACE "Sire"
                                  SYMBOLS ${SIRE_LINPACK_FUNCTIONS} )

        include_directories( BEFORE ${CMAKE_BINARY_DIR}/src/libs/SireMaths )
      else()
        message( STATUS "Fortran/C++ compilers are incompatible. Disabling Fortran code!" )
        set( SIRE_DISABLE_FORTRAN ON )
        add_definitions( "-DSIRE_DISABLE_FORTRAN" )        
      endif()
    endif()
  endif()
endif()

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

if ( NOT SIRE_DISABLE_FORTRAN )
  # Define the sources in SireBLASPACK
  set( SIREBLASPACK_SOURCES
   
       ${SIRE_BLAS_SOURCES}
       ${SIRE_LAPACK_SOURCES}
       ${SIRE_LINPACK_SOURCES}

       sire_blaspack.f
     )
endif()

# Define the headers in SireMaths
set ( SIREMATHS_HEADERS
      accumulator.h
      align.h
      axisset.h
      boys.h
      complex.h
      constants.h
      distvector.h
      errors.h
      freeenergyaverage.h
      gamma.h
      histogram.h
      line.h
      linearap.h
      maths.h
      matrix.h
      multidouble.h
      multifixed.h
      multifloat.h
      multiint.h
      multiuint.h
      nmatrix.h
      nvector.h
      n4matrix.h
      plane.h
      quaternion.h
      rangenerator.h
      rational.h
      rotate.h
      sphere.h
      torsion.h
      triangle.h
      trigmatrix.h
      vector.h
      vector3d.hpp
      vectorproperty.h

      sire_blas.h
      sire_lapack.h
      sire_linpack.h

      third_party/eig3/eig3.cpp
    )

# Define the sources in SireMaths
set ( SIREMATHS_SOURCES

      register_siremaths.cpp

      accumulator.cpp
      align.cpp
      axisset.cpp
      boys.cpp
      complex.cpp
      distvector.cpp
      errors.cpp
      freeenergyaverage.cpp
      gamma.cpp
      histogram.cpp
      line.cpp
      linearap.cpp
      matrix.cpp
      multidouble.cpp
      multifixed.cpp
      multifloat.cpp
      multiint.cpp
      multiuint.cpp
      nmatrix.cpp
      nvector.cpp
      n4matrix.cpp
      plane.cpp
      quaternion.cpp
      rangenerator.cpp
      rational.cpp
      sphere.cpp
      torsion.cpp
      triangle.cpp
      trigmatrix.cpp
      vector.cpp
      vectorproperty.cpp

      sire_blas.cpp
      sire_lapack.cpp
      sire_linpack.cpp

      ${SIREMATHS_HEADERS}
    )

if (NOT SIRE_DISABLE_FORTRAN)
  add_library (SireBLASPACK ${SIREBLASPACK_SOURCES})

  set_target_properties (SireBLASPACK PROPERTIES
                         VERSION ${SIRE_VERSION}
                         SOVERSION ${S_VERSION_MAJOR}
                        )

  if ( SIREBLASPACK_HAS_LIBRARIES )
    target_link_libraries (SireBLASPACK
                           ${SIREBLASPACK_LIBRARIES}
                          )
  endif()
endif()

add_library (SireMaths ${SIREMATHS_SOURCES})

set_target_properties (SireMaths PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

# Link to other Sire libraries, and GSL
target_link_libraries (SireMaths
                       SireID
                       SireBase
                       SireUnits
                       SireStream
                       ${SIREMATHS_EXTRA_LIBRARIES}
                      )

if (NOT SIRE_DISABLE_FORTRAN)
  target_link_libraries (SireMaths SireBLASPACK )
endif()

save_sire_variable( "SIREMATHS_EXTRA_LIBRARIES" "${SIREMATHS_EXTRA_LIBRARIES}" )

# installation
install( TARGETS SireMaths EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

if (NOT SIRE_DISABLE_FORTRAN)
  install( TARGETS SireBLASPACK EXPORT SireLibraries
           RUNTIME DESTINATION ${SIRE_BIN}
           LIBRARY DESTINATION ${SIRE_LIBS}
           ARCHIVE DESTINATION ${SIRE_ARCHIVES}
         )
endif()

install( FILES ${SIREMATHS_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireMaths )
