!                    *******************
                     SUBROUTINE IDWM_T2D
!                    *******************
     &(STAVAL,POINTS,POINVAL,NPOIN,NSTA)
!
!***********************************************************************
! TELEMAC2D   V7P1                                   29/07/2016
!***********************************************************************
!
!brief    USES INVERSE DISTANCE WEIGHTING METHOD TO COMPUTE A WIND FIELD
!+               THAT VARIES IN TIME AND SPACE
!
!history  P. PRODANOVIC (RIGGS ENGINEERING LTD)
!+        23/04/2014
!+        V6P3
!+   Initial version.
!
!history  P. PRODANOVIC (RIGGS ENGINEERING LTD)
!+        29/05/2016
!+        V7P0
!+   Improvement to manage divisions by 0 when station and mesh nodes
!+   are on the same point. Also, now allows proper interpolation for
!+   directions between 359 and 1 degrees azimuth.
!
!history M. GAMBARINI (POLIMI)
!+        3/12/2020
!+        V8P1
!+   Corrected cases in which station and mesh nodes are the same point
!+   or have the same ordinate.
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
      DOUBLE PRECISION :: DISTS(NSTA), WEIGHTS(NSTA)
      INTEGER :: I, J
      DOUBLE PRECISION :: DISTTHRES
      DOUBLE PRECISION :: TOTWEIGHTS
      INTEGER :: INVDISTEXP
      INTEGER :: MININDEX
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     BELOW THIS DISTANCE POINTS ARE CONSIDERED COINCIDENT
      DISTTHRES = 1.E-2
      INVDISTEXP = 2

      DO I = 1, NPOIN
        DO J = 1, NSTA
          DISTS(J) = SQRT( (POINTS(I,1) - STAVAL(J,1))**2 +
     &      (POINTS(I,2) - STAVAL(J,2))**2 )
        ENDDO
        MININDEX = MINLOC(DISTS, 1)
        IF (MINVAL(DISTS).LT.DISTTHRES) THEN
          WRITE(LU, *) "FOUND POINT CLOSE TO STATION"
          POINVAL(I) = STAVAL(MININDEX, 3)
        ELSE
          TOTWEIGHTS = 0.
!     COMPUTE WEIGHTS AND THEIR TOTAL FOR NORMALIZATION
          DO J = 1, NSTA
            WEIGHTS(J) = 1/DISTS(J)**INVDISTEXP
            TOTWEIGHTS = TOTWEIGHTS + WEIGHTS(J)
          ENDDO
!     NORMALIZE WEIGHTS AND COMBINE THE VALUES
          POINVAL(I) = 0.
          DO J = 1, NSTA
            POINVAL(I) = POINVAL(I) +
     &        WEIGHTS(J)/TOTWEIGHTS*STAVAL(J,3)
          ENDDO

        ENDIF

      ENDDO
!
!-----------------------------------------------------------------------
!
      RETURN
      END
