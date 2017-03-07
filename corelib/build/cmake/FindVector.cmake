# Find and test vector support (i.e. SSE, AVX, AVX512)
# used within Sire
#
#  SIRE_VECTOR_FLAGS    - flags passed to the compiler to support vectorisation,
#                         and to inform the application of what is available
# Call this function using
#
# GET_SIRE_VECTOR_FLAGS( omp-simd-flag sse2-flag avx-flag avx512f-flag )
#
# This will test what is supported and will return the value as SIRE_VECTOR_FLAGS
#
# This will read SIRE_DISABLE_SSE, SIRE_DISABLE_AVX and SIRE_DISABLE_AVX512F
# to decide whether or not these levels of vectorisation should be disabled

function( GET_SIRE_VECTOR_FLAGS OMP_SIMD_FLAG SSE2_FLAG AVX_FLAG AVX512F_FLAG )
  #message( STATUS "Checking flags ${OMP_SIMD_FLAG} | ${SSE2_FLAG} | ${AVX_FLAG} | ${AVX512F_FLAG}" )

  if (NOT SIRE_VECTORISE)
    message( STATUS "Disabling all vectorisation!" )
    set( SIRE_VECTOR_FLAGS "" )
    return()
  endif()

  set( ENABLE_AVX512F 1 )
  set( ENABLE_AVX 1 )
  set( ENABLE_SSE2 1 )

  if ( SIRE_DISABLE_SSE )
    set( ENABLE_SSE2 0 )
  endif()

  if ( SIRE_DISABLE_AVX )
    set( ENABLE_AVX 0 )
  endif()

  if ( SIRE_DISABLE_AVX512F )
    set( ENABLE_AVX512F 0 )
  endif()

  # check that the compiler has the required flags
  CHECK_CXX_COMPILER_FLAG( "${OMP_SIMD_FLAG}" HAVE_OMP_SIMD_FLAG )
  CHECK_CXX_COMPILER_FLAG( "${SSE2_FLAG}" HAVE_SSE2_FLAG )
  CHECK_CXX_COMPILER_FLAG( "${SSE4_FLAG}" HAVE_SSE4_FLAG )
  CHECK_CXX_COMPILER_FLAG( "${AVX_FLAG}" HAVE_AVX_FLAG )
  CHECK_CXX_COMPILER_FLAG( "${AVX512F_FLAG}" HAVE_AVX512F_FLAG )

  # Need this to check that C source files can run
  INCLUDE(CheckCSourceRuns)

  if ( HAVE_OMP_SIMD_FLAG )
    set( HAVE_OMP_SIMD_SYSTEM 1 )
  else()
    set( HAVE_OMP_SIMD_SYSTEM 0 )
  endif()

  # check if we have a computer that can compile and run SSE2
  set( HAVE_SSE2_SYSTEM 0 )

  if ( HAVE_SSE2_FLAG )
    set(CMAKE_OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_FLAGS "${SSE2_FLAG}")

    CHECK_C_SOURCE_RUNS("
       #include <emmintrin.h>
       int main()
       {
          __m128d a = _mm_set_pd1(1.0);
          return 0;
       }"
       CAN_RUN_SSE2_PROGRAM)
    set(CMAKE_REQUIRED_FLAGS ${CMAKE_OLD_REQUIRED_FLAGS})

    if (CAN_RUN_SSE2_PROGRAM)
      set( HAVE_SSE2_SYSTEM 1 )
    endif()
  endif()

  # check if we have a computer that can compile and run AVX
  set( HAVE_AVX_SYSTEM 0 )

  if ( HAVE_AVX_FLAG )
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

  # check if we have a computer that can compile and run AVX512F
  set( HAVE_AVX512F_SYSTEM 0 )

  if ( HAVE_AVX512F_FLAG )
    set(CMAKE_OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_FLAGS "${AVX512F_FLAG}")
    CHECK_C_SOURCE_RUNS("
       #include <immintrin.h>
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

  if (HAVE_OMP_SIMD_SYSTEM)
    message( STATUS "We can compile and run a program using pragma omp simd" )
  else()
    message( STATUS "Not possible to compile using pragma omp simd" )
  endif()

  if (HAVE_SSE2_SYSTEM)
    message( STATUS "We can compile and run a program with SSE2 support" )
  else()
    message( STATUS "Not possible to compile and run a program with SSE2 support" )
  endif()

  if (HAVE_AVX_SYSTEM)
    message( STATUS "We can compile and run a program with AVX support" )
  else()
    message( STATUS "Not possible to compile and run a program with AVX support" )
  endif()

  if (HAVE_AVX512F_SYSTEM)
    message( STATUS "We can compile and run a program with AVX512-F support" )
  else()
    message( STATUS "Not possible to compile and run a program with AVX512-F support" )
  endif()

  if (HAVE_AVX512F_SYSTEM AND ENABLE_AVX512F)
    set( VECTOR_FLAGS "${AVX512F_FLAG} -DSIRE_USE_AVX -DSIRE_USE_AVX2 -DSIRE_USE_AVX512F" )
    message( STATUS "Compiling using AVX-512F. Note that the resulting binary will only work on really new Intel processors" )
  elseif ( HAVE_AVX_SYSTEM AND ENABLE_AVX )
    set( VECTOR_FLAGS "${AVX_FLAG} -DSIRE_USE_AVX" )
  elseif ( HAVE_SSE2_SYSTEM AND ENABLE_SSE2 )
    set( VECTOR_FLAGS "${SSE2_FLAG} -DSIRE_USE_SSE" )
  endif()

  if (HAVE_OMP_SIMD_SYSTEM)
    set( VECTOR_FLAGS "${OMP_SIMD_FLAG} ${VECTOR_FLAGS}" )
  endif()

  message( STATUS "Using vectorisation flags ${VECTOR_FLAGS}" )

  # return this value as SIRE_VECTOR_FLAGS to the parent caller
  set( SIRE_VECTOR_FLAGS "${VECTOR_FLAGS}" PARENT_SCOPE )

endfunction( GET_SIRE_VECTOR_FLAGS )
