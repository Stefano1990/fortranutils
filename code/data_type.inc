C
C RJ: - Modified data type constants. Necessary for the
C       values of integer and real types to be different,
c       else the code will not work on platforms where
C       the sizes of integer and real types are the same
C       (like the Cray T3E). Also taking the opportunity
C       to define similar constants for the benefit of
C       MPI, so as to make it easier to port the code 
C       (for instance, REAL_TYPE_MPI must be defined
C       MPI_REAL on the T3E and MPI_DOUBLE_PRECISION on
C       SGIs)
C
C       In case it should be necessary to extract the 
C       size of the data types from the values of these
C       constants (and this is not done at the moment),
C       I have defined INTEGER_TYPE equal to its size in
C       bytes and REAL_TYPE to its size + 64. Similar
C       schemes could be adopted to squeeze more types
C       into one integer.
C
      INTEGER INTEGER_TYPE
      INTEGER REAL_TYPE
      INTEGER INTEGER_TYPE_MPI
      INTEGER REAL_TYPE_MPI
      INTEGER REAL_TYPE_MPI_WRITE
      INTEGER REAL_TYPE_MPI_READ

#if defined SGI || defined DEC
      DATA INTEGER_TYPE     / 4 /
      DATA REAL_TYPE        / 72 /
      DATA INTEGER_TYPE_MPI / MPI_INTEGER /
      DATA REAL_TYPE_MPI    / MPI_DOUBLE_PRECISION /
#elif defined T3E
      DATA INTEGER_TYPE     / 8 /
      DATA REAL_TYPE        / 72 /
      DATA INTEGER_TYPE_MPI / MPI_INTEGER /
      DATA REAL_TYPE_MPI    / MPI_REAL /
#endif
#ifdef WRITE_DOUBLE
      DATA REAL_TYPE_MPI_WRITE / MPI_REAL8 /
#else
      DATA REAL_TYPE_MPI_WRITE / MPI_REAL4 /
#endif
#ifdef READ_DOUBLE
      DATA REAL_TYPE_MPI_READ  / MPI_REAL8 /
#else
      DATA REAL_TYPE_MPI_READ  / MPI_REAL4 /
#endif
