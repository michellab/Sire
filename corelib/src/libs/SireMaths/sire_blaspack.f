
c
c   Stub that provides information about the BLAS, LAPACK and LINPACK
c   libraries - are they linked here, or are machine optimised
c   versions being used?
c

       integer function sire_using_internal_blas()
       implicit none
       
       sire_using_internal_blas = 1
       return
       end
       
       integer function sire_using_internal_lapack()
       implicit none
       
       sire_using_internal_lapack = 1
       return
       end
       
       integer function sire_using_internal_linpack()
       implicit none
       
       sire_using_internal_linpack = 1
       return
       end
       