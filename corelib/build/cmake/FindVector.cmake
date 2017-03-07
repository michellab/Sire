# Find and test vector support (i.e. SSE, AVX, AVX512)
# used within Sire
#
#  SIRE_VECTOR_FLAGS    - flags passed to the compiler to support vectorisation,
#                         and to inform the application of what is available

function( GET_SIRE_VECTOR_FLAGS SSE2_FLAG SSE4_FLAG AVX_FLAG AVX512F_FLAG )
  # look for the header files needed to support different vector schemes
  check_include_files( emmintrin.h HAVE_EMMINTRIN_H ) # SSE2
  check_include_files( smmintrin.h HAVE_SMMINTRIN_H ) # SSE4
  check_include_files( immintrin.h HAVE_IMMINTRIN_H ) # AVX
  check_include_files( zmmintrin.h HAVE_ZMMINTRIN_H ) # AVX512

  # check that the compiler has the required flags
  CHECK_CXX_COMPILER_FLAG( "${SSE2_FLAG}" HAVE_SSE2_FLAG )
  CHECK_CXX_COMPILER_FLAG( "${SSE4_FLAG}" HAVE_SSE4_FLAG )
  CHECK_CXX_COMPILER_FLAG( "${AVX_FLAG}" HAVE_AVX_FLAG )
  CHECK_CXX_COMPILER_FLAG( "${AVX512F_FLAG}" HAVE_AVX512F_FLAG )

  if ( HAVE_EMMINTRIN_H AND HAVE_SSE2_FLAG )
    set( HAVE_SSE2_SYSTEM 1 )
    message( STATUS "We can compile and run a program with SSE2 support" )
  else()
    set( HAVE_SSE2_SYSTEM 0 )
    message( STATUS "Not possible to compile and run a program with SSE2 support" )
  endif()

  if ( HAVE_SMMINTRIN_H AND HAVE_SSE4_FLAG )
    set( HAVE_SSE4_SYSTEM 1 )
    message( STATUS "We can compile and run a program with SSE4 support" )
  else()
    set( HAVE_SSE4_SYSTEM 0 )
    message( STATUS "Not possible to compile and run a program with SSE4 support" )
  endif()

  # Need this to check that C source files can run
  INCLUDE(CheckCSourceRuns)

  # check if we have a computer that can compile and run AVX
  set( HAVE_AVX_SYSTEM 0 )

  if ( HAVE_IMMINTRIN_H AND HAVE_AVX_FLAG )
    set(CMAKE_OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_FLAGS "${AVX_FLAG}")
    CHECK_C_SOURCE_RUNS("
       #include <immintrin.h>
       int main()
       {
          __m256i a = _mm256_setzero_si256();
          return 0;
       }"
       CAN_RUN_AVX_PROGRAM)
    set(CMAKE_REQUIRED_FLAGS ${CMAKE_OLD_REQUIRED_FLAGS})

    if (CAN_RUN_AVX_PROGRAM)
      set( HAVE_AVX_SYSTEM 1 )
    endif()
  endif()

  if (HAVE_AVX_SYSTEM)
    message( STATUS "We can compile and run a program with AVX support" )
  else()
    message( STATUS "Not possible to compile and run a program with AVX support" )
  endif()

  # check if we have a computer that can compile and run AVX512F
  set( HAVE_AVX512F_SYSTEM 0 )

  if ( HAVE_ZMMINTRIN_H AND HAVE_AVX512F_FLAG )
    set(CMAKE_OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_FLAGS "${AVX512F_FLAG}")
    CHECK_C_SOURCE_RUNS("
       #include <zmmintrin.h>
       int main()
       {
          __m512d a = _mm512_setzero_pd();
          return 0;
       }"
       CAN_RUN_AVX512F_PROGRAM)
    set(CMAKE_REQUIRED_FLAGS ${CMAKE_OLD_REQUIRED_FLAGS})

    if (CAN_RUN_AVX512F_PROGRAM)
      set( HAVE_AVX512F_SYSTEM 1 )
    endif()
  endif()  

  if (HAVE_AVX512F_SYSTEM)
    message( STATUS "We can compile and run a program with AVX512-F support" )
  else()
    message( STATUS "Not possible to compile and run a program with AVX512-F support" )
  endif()

  message( STATUS "CAN WE READ GLOBAL VARIABLES? ${SIRE_DISABLE_AVX}" )

endfunction( GET_SIRE_VECTOR_FLAGS )
