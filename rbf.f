!                    *******************
                     SUBROUTINE RBF_2D
!                    *******************
     &(STAVAL,POINTS,POINVAL,NPOIN,NSTA)
!
!***********************************************************************
! TELEMAC3D   V8P1
!***********************************************************************
!
!brief    USES RADIAL BASIS FUNCTION INTERPOLATION TO COMPUTE A WIND FIELD
!+               THAT VARIES IN TIME AND SPACE
!
!history M. GAMBARINI (POLIMI)
!+        22/12/2020
!+        V8P1
!+
!+
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_SPECIAL
!
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     VARIABLES THAT ARE USED BY THE SUBROUTINE TO PASS DATA IN AND OUT
      INTEGER,INTENT(IN) :: NPOIN, NSTA
      DOUBLE PRECISION, INTENT(IN) :: POINTS(NPOIN,2)
      DOUBLE PRECISION, INTENT(IN) :: STAVAL(NSTA,3)
      DOUBLE PRECISION, INTENT(INOUT) :: POINVAL(NPOIN)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     VARIABLES INTERNAL TO THE SUBROUTINE
      INTEGER :: I, J, L
      DOUBLE PRECISION :: A(NSTA, NSTA)
      DOUBLE PRECISION :: X(NSTA)
      DOUBLE PRECISION :: B(NSTA)
      DOUBLE PRECISION :: DISTSQUARED
      DOUBLE PRECISION :: K

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: COL
      DOUBLE PRECISION :: PIVOTRATIO, TEMPSUM, TEMP

      INTEGER :: MAXPIVOTROW
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
  ! EXPONENTIAL COEFFICIENT (1/TYPICAL DISTANCE BETWEEN STATIONS^2)
      K = 2.E-5

      DO I=1,NSTA
        DO J=1,NSTA
          DISTSQUARED = (STAVAL(I,1)-STAVAL(J,1))**2 +
     &                  (STAVAL(I,2)-STAVAL(J,2))**2
          A(I,J) = EXP(-K*DISTSQUARED)
        ENDDO
        B(I) = STAVAL(I,3)
      ENDDO

      PRINT *, A

      ! SOLVE THE LINEAR SYSTEM HERE
      ! GAUSS ELIMINATION
      DO I=1,NSTA-1
        !FIND MAXIMUM PIVOT
        ALLOCATE(COL(NSTA-I+1))
        DO J=I,NSTA
          COL(J-I+1) = ABS(A(J,I))
        ENDDO
        MAXPIVOTROW = I -1 + MAXLOC(COL, 1)
        DEALLOCATE(COL)
        !ROW SUBSTITUTION
        DO J=I,NSTA
          TEMP = A(I,J)
          A(I,J) = A(MAXPIVOTROW,J)
          A(MAXPIVOTROW, J) = TEMP
        ENDDO
        TEMP = B(I)
        B(I) = B(MAXPIVOTROW)
        B(MAXPIVOTROW) = TEMP

        DO J=I+1,NSTA
          IF (A(J,I).NE.0.) THEN
            PIVOTRATIO = A(J,I)/A(I,I)
            DO L=I,NSTA
              A(J,L) = A(J,L) - PIVOTRATIO*A(I,L)
            ENDDO
            B(J) = B(J) - PIVOTRATIO*B(I)
          ENDIF
        ENDDO

      ENDDO


      ! BACK SUBSTITUTION
      X(NSTA) = B(NSTA)/A(NSTA,NSTA)

      DO I=NSTA-1, 1, -1
        TEMPSUM = 0
        DO J = I+1,NSTA
          TEMPSUM = TEMPSUM + A(I,J)*X(J)
        ENDDO
        X(I) = (B(I) - TEMPSUM)/A(I,I)
      ENDDO

    !  PRINT*, A

      DO I = 1,NPOIN
        POINVAL(I) = 0
        DO J = 1,NSTA
          DISTSQUARED = (POINTS(I,1)-STAVAL(J,1))**2 +
     &                  (POINTS(I,2)-STAVAL(J,2))**2
          POINVAL(I) = POINVAL(I) + X(J) * EXP(-K*DISTSQUARED)
        ENDDO
      ENDDO

!
!-----------------------------------------------------------------------
!
      RETURN
      END
