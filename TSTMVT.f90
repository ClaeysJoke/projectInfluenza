PROGRAM TSTMVT
  !
  !     Simple Test program for MVT Methods
  !
  USE PRECISION_MODEL
  USE MVSTAT
  IMPLICIT NONE
  INTEGER                                   :: I, IVLS, INF, NU
  INTEGER,                        PARAMETER :: IP = 1, MX = 100000
  REAL(KIND=STND),                PARAMETER :: ABSEPS = 1E-5_STND
  !
  INTEGER,                        PARAMETER :: M = 5, N = 4
  REAL(KIND=STND), DIMENSION(N,N)         :: COV
  REAL(KIND=STND), DIMENSION(M,N)         :: CNS
  REAL(KIND=STND), DIMENSION(M)            :: LW
  REAL(KIND=STND), DIMENSION(M)            :: UP
  REAL(KIND=STND), DIMENSION(M), PARAMETER :: DL = 0
  INTEGER,         DIMENSION(M), PARAMETER :: FN = 2
  REAL(KIND=STND) :: VAL, ERR
  COV = 0; CNS = 0; LW = 0; UP = 1; NU = 8;
  DO I = 1, N
     COV(I,I) = 1; CNS(I,I) = 1
  END DO
  CNS(M,:) = 1
  CALL MVPRNT( N, COV, NU, M, LW, CNS, UP, FN, DL,      &
       MX, ABSEPS, "      4-d Simplex", IP )
  COV = 1E0_STND/N;
  DO I = 1, N; COV(I,I) = 1; END DO
     CALL MVPRNT( N, COV, NU, M, LW, CNS, UP, FN, DL,      &
          MX, ABSEPS, "      4-d Simplex", IP )
     UP(M) = -1    ! Infeasible Constraint
     CALL MVPRNT( N, COV, NU, M, LW, CNS, UP, FN, DL,      &
          MX, ABSEPS, "   Bad 4-d Simplex", IP )
   END PROGRAM TSTMVT
!
SUBROUTINE MVPRNT( N, COVRNC, NU, M, LOWER, CONSTR, UPPER, INFIN, DELTA,  &
                MAXPTS, ABSEPS, LABEL, IP )
USE PRECISION_MODEL
USE MVSTAT
IMPLICIT NONE
INTEGER,                         INTENT(IN) :: N, NU, M, MAXPTS, IP
REAL(KIND=STND), DIMENSION(M),   INTENT(IN) :: LOWER, UPPER, DELTA
REAL(KIND=STND), DIMENSION(N,N), INTENT(IN) :: COVRNC
REAL(KIND=STND), DIMENSION(M,N), INTENT(IN) :: CONSTR
INTEGER,         DIMENSION(M),   INTENT(IN) :: INFIN
CHARACTER(LEN=*),                INTENT(IN) :: LABEL
REAL(KIND=STND),                 INTENT(IN) :: ABSEPS
!
INTEGER                                     :: I, INF, IVLS
REAL(KIND=STND)                             :: LWR, UPR, ALPHA, TALPHA
REAL(KIND=STND)                             :: ERROR, VALUE
!
PRINT "(/10X,A)", LABEL
PRINT "(""           Number of Dimensions is "",I3)", N
PRINT "(""          Number of Constraints is "",I3)", M
PRINT "(""      Number of Degrees of Freedom is "",I3)", NU
PRINT "(""     Maximum # of Function Evaluations is "",I7)", MAXPTS
IF ( IP .GT. 0 ) THEN
   PRINT "(""    Delta  Lower  Upper     Constraints "")"
   DO I = 1, M
      IF ( INFIN(I) < 0 ) THEN
         PRINT "(I2, "" -00     00   "", 7F7.3)",  CONSTR(I,1:N)
      ELSE IF ( INFIN(I) == 0 ) THEN
         PRINT "(I2, "" -00   "", 10F7.3)", I, UPPER(I), CONSTR(I,1:N)
      ELSE IF ( INFIN(I) == 1 ) THEN
         PRINT "(I2,F7.3,""  00   "",9F7.3)", I, LOWER(I), CONSTR(I,1:N)
      ELSE
         PRINT "(I2, 11F7.3)",                                       &
           I, DELTA(I), LOWER(I), UPPER(I), CONSTR(I,1:N)
      ENDIF
   END DO
   PRINT "(""     Lower Left Covariance Matrix "")"
   DO I = 1, N
      PRINT "(I2, 11F7.3)", I, COVRNC(I,1:I)
   END DO
END IF
CALL MVDIST( N, COVRNC, NU, M, LOWER, CONSTR, UPPER, INFIN, DELTA,      &
             MAXPTS, ABSEPS, 0E0_STND, ERROR, VALUE, IVLS, INF )
PRINT "(5X, ""Value(Error): "",F9.7,""("",F9.7,"")"")", VALUE, ERROR
PRINT "(5X,""Evaluations(Inform): "",I8,""("",I2,"")"")", IVLS, INF
END SUBROUTINE MVPRNT
