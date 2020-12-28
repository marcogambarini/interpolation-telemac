!                    ****************
                     SUBROUTINE METEO
!                    ****************
!
     &(PATMOS,WINDX,WINDY,FUAIR,FVAIR,AT,LT,NPOIN,VENT,ATMOS,
     & ATMFILEA,ATMFILEB,FILES,LISTIN,
     & PATMOS_VALUE,AWATER_QUALITY,PLUIE,AOPTWIND,AWIND_SPD)
!
!***********************************************************************
! TELEMAC2D   V7P2
!***********************************************************************
!
!brief    COMPUTES ATMOSPHERIC PRESSURE AND WIND VELOCITY FIELDS
!+               (IN GENERAL FROM INPUT DATA FILES).
!
!warning  CAN BE ADAPTED BY USER
!
!history  J-M HERVOUET (LNHE)
!+        02/01/2004
!+        V5P4
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!history  J-M HERVOUET (EDF R&D, LNHE)
!+        30/01/2013
!+        V6P3
!+   Now 2 options with an example for reading a file. Extra arguments.
!
!history  C.-T. PHAM (LNHE)
!+        09/07/2014
!+        V7P0
!+   Reading a file of meteo data for exchange with atmosphere
!+   Only the wind is used here
!
!history R.ATA (LNHE)
!+        09/11/2014
!+        V7P0
!+  introducion of water quality option + pluie is introduced as
!+   an optional parameter + remove of my_option which is replaced
!+   by a new keyword + value of patmos managed also with a new keyword
!
!history  J-M HERVOUET (EDF R&D, LNHE)
!+        07/01/2015
!+        V7P0
!+  Adding optional arguments to remove USE DECLARATIONS_TELEMAC2D.
!
!history R.ATA (LNHE)
!+        16/11/2015
!+        V7P0
!+  Adding USE WAQTEL...
!
!history A. LEROY (LNHE)
!+        25/11/2015
!+        V7P1
!+  INTERPMETEO now writes directly in variables of WAQTEL which
!+  can be used by the other modules. This makes it possible to
!+  remove subsequent calls to INTERPMETEO in TELEMAC3D
!
!history  P. PRODANOVIC (RIGGS ENGINEERING LTD)
!+        15/06/2016
!+        V7P0
!+   Converts the wind data to cartesian form, then interpolates. This
!+   eliminates errors when interpolating direction between 359 and 1
!+   degrees azimuth.
!
!history J.-M. HERVOUET (RETIRED)
!+        01/07/2017
!+        V7P2
!+  Setting of UL moved outside the test IF(LT.EQ.0)... After a post by
!+  Qilong Bi (thanks Qilong...).
!
!history M. GAMBARINI (POLIMI)
!+        02/12/2020
!+        V8P1
!+  Adaptation for wind variable in time and space, temperature constant in
!+  space and variable in time
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| AT             |-->| TIME
!| ATMFILEA       |-->| LOGICAL UNIT OF THE ASCII ATMOSPHERIC FILE
!| ATMFILEB       |-->| LOGICAL UNIT OF THE BINARY ATMOSPHERIC FILE
!| ATMOS          |-->| YES IF PRESSURE TAKEN INTO ACCOUNT
!| FILES          |-->| BIEF_FILES STRUCTURES OF ALL FILES
!| FUAIR          |<->| VELOCITY OF WIND ALONG X, IF CONSTANT
!| FVAIR          |<->| VELOCITY OF WIND ALONG Y, IF CONSTANT
!| LISTIN         |-->| IF YES, PRINTS INFORMATION
!| LT             |-->| ITERATION NUMBER
!| NPOIN          |-->| NUMBER OF POINTS IN THE MESH
!| PATMOS         |<--| ATMOSPHERIC PRESSURE
!| PATMOS_VALUE   |-->| VALUE OF ATMOSPHERIC PRESSURE IS CONSTANT
!| VENT           |-->| YES IF WIND TAKEN INTO ACCOUNT
!| WINDX          |<--| FIRST COMPONENT OF WIND VELOCITY
!| WINDY          |<--| SECOND COMPONENT OF WIND VELOCITY
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_WAQTEL,ONLY: PVAP,RAY3,NWIND,NEBU,TAIR,
     &                              TAIR_VALUE,HREL,RAINFALL,
     &                              EVAPORATION,ATMOSEXCH
      USE DECLARATIONS_TELEMAC3D, ONLY : X,Y
!
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)             :: LT,NPOIN,ATMFILEA,ATMFILEB
      LOGICAL, INTENT(IN)             :: ATMOS,VENT,LISTIN
      DOUBLE PRECISION, INTENT(INOUT) :: WINDX(NPOIN),WINDY(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: PATMOS(*)
      DOUBLE PRECISION, INTENT(IN)    :: AT,PATMOS_VALUE
      DOUBLE PRECISION, INTENT(INOUT) :: FUAIR,FVAIR
      TYPE(BIEF_FILE), INTENT(IN)     :: FILES(*)
!     OPTIONAL
      LOGICAL, INTENT(IN)          ,OPTIONAL :: AWATER_QUALITY
      TYPE(BIEF_OBJ), INTENT(INOUT),OPTIONAL :: PLUIE
      INTEGER, INTENT(IN)          ,OPTIONAL :: AOPTWIND
      DOUBLE PRECISION, INTENT(IN) ,OPTIONAL :: AWIND_SPD(2)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      LOGICAL WATER_QUALITY
      INTEGER UL,OPTWIND
      DOUBLE PRECISION COEF
      DOUBLE PRECISION UAIR,VAIR,WIND_SPD(2)
!     EXCHANGE WITH ATMOSPHERE
      DOUBLE PRECISION PATM,WW,TA
!
      DOUBLE PRECISION, PARAMETER :: EPS = 1.D-3
!
!     ######################################################################
!     IDWM WIND INTERPOLATION CUSTOM VARIABLES
!     ######################################################################
!
      INTEGER I, NUMSTA, NUMMETEOTIMES, A, B, J, K, JUNK
      DOUBLE PRECISION THETA_RAD, TMPDIR, TMPSPD, PI, DTR

!
!     COORDINATES OF THE STATIONS UTM
!
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: XX, YY, AT_WIND
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: OUT_WSPD
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: OUT_WDIRX
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: OUT_WDIRY
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DATAFROMFILE
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: POINTS
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: S_LAST
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: S_NEXT
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: D_LAST
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: D_NEXT

!     ADDED ON 2016.05.26
!     THIS IS THE X AND Y COMPONENT OF THE WIND READ FROM FILE
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WINDX_LAST
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WINDY_LAST
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WINDX_NEXT
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WINDY_NEXT
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: INPSTA_WINDX
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: INPSTA_WINDY

!
!     ######################################################################
!
!-----------------------------------------------------------------------
!
!     DATA THAT YOU DECLARE AND READ HERE ONCE IN A FILE MAY HAVE TO BE
!     KEPT BECAUSE THIS SUBROUTINE IS CALLED AT EVERY TIME STEP.
!     WITHOUT THE SAVE COMMAND, ALL LOCAL DATA ARE FORGOTTEN IN THE NEXT
!     CALL.
!
      SAVE
!
!-----------------------------------------------------------------------
!
!     DEFAULT VALUES OF PARAMETERS WHEN THEY ARE NOT GIVEN
!
      WATER_QUALITY=.FALSE.
      IF(PRESENT(AWATER_QUALITY)) WATER_QUALITY=AWATER_QUALITY
      OPTWIND=1
      IF(PRESENT(AOPTWIND)) OPTWIND=AOPTWIND
      WIND_SPD(1)=0.D0
      WIND_SPD(2)=0.D0
      IF(PRESENT(AWIND_SPD)) THEN
        WIND_SPD(1)=AWIND_SPD(1)
        WIND_SPD(2)=AWIND_SPD(2)
      ENDIF
!
!-----------------------------------------------------------------------
!
      UL = FILES(ATMFILEA)%LU

!     AT FIRST TIMESTEP
!
      IF(LT.EQ.0) THEN
!
!       ATMOSPHERIC PRESSURE AND AIR TEMPERATURE
!
        IF(ATMOS.OR.WATER_QUALITY) THEN
          CALL OV('X=C     ', X=PATMOS, C=PATMOS_VALUE, DIM1=NPOIN)
        ENDIF
        IF(WATER_QUALITY) THEN
          CALL OV('X=C     ', X=TAIR%R, C=TAIR_VALUE, DIM1=NPOIN)
        ENDIF
!
!       WIND :
!
        IF(VENT.OR.WATER_QUALITY) THEN

!
!         ######################################################################
!         IDWM WIND INTERPOLATION; THIS IS EXECUTED ONLY ONCE AT THE START
!         ######################################################################
!
          IF(OPTWIND.EQ.3) THEN
        ! READ BLANK LINE AT BEGINNING OF FILE
            READ(UL,*)
        ! READ NUMSTA AND NUMMETEOTIMES
            READ(UL,*) NUMSTA, NUMMETEOTIMES


          !ALLOCATE THE ARRAYS
            ALLOCATE(XX(NUMSTA), YY(NUMSTA), AT_WIND(NUMMETEOTIMES))
            ALLOCATE(DATAFROMFILE(NUMMETEOTIMES,NUMSTA*2+2))
            ALLOCATE(POINTS(NPOIN,2))
            ALLOCATE(S_LAST(NUMSTA), D_LAST(NUMSTA))
            ALLOCATE(S_NEXT(NUMSTA), D_NEXT(NUMSTA))
            ALLOCATE(WINDX_LAST(NUMSTA), WINDY_LAST(NUMSTA))
            ALLOCATE(WINDX_NEXT(NUMSTA), WINDY_NEXT(NUMSTA))
            ALLOCATE(INPSTA_WINDX(NUMSTA,3), INPSTA_WINDY(NUMSTA,3))
            ALLOCATE(OUT_WSPD(NPOIN),OUT_WDIRX(NPOIN),OUT_WDIRY(NPOIN))
!
          ! READ STATION COORDINATES
            DO B = 1,NUMSTA
              READ(UL,*) XX(B), YY(B)
           !WRITE(*,*) XX(B), YY(B)
            ENDDO
!
          ! READ THE WIND TIME SERIES FROM THE INPUT FILE
          ! FIRST COLUMN IS TIME IN SECONDS, REST OF COLUMNS ARE WSPD
          ! AND WDIR FOR EACH STATION READ
          ! LAST COLUMN IS TEMPERATURE (ASSUMED UNIFORM FOR NOW)
            DO A = 1,NUMMETEOTIMES
              READ(UL,*) (DATAFROMFILE(A,B), B=1,1+NUMSTA*2+1)
!             PRINT*, DATAFROMFILE(A, 1+NUMSTA*2+1)
            ENDDO
!
          ! EXTRACT AT_WIND FROM WIND(A,B); FIRST COLUMN IS TIME IN SECONDS
            DO A = 1,NUMMETEOTIMES
              AT_WIND(A) = DATAFROMFILE(A,1)
            ENDDO
!
            PRINT*, 'NPOIN',NPOIN
          ! ASSEMBLE THE POINTS ARRAY FOR IDWM FUNCTION
            DO I = 1,NPOIN
              POINTS(I,1) = X(I)
              POINTS(I,2) = Y(I)
            !  PRINT*, X(I), Y(I)
            ENDDO
!
! #######################################################################
!
          ENDIF
        ENDIF
      ENDIF
!
!-----------------------------------------------------------------------
!
!     FOR THE REMAINING TIME STEPS
!
!      IF(VENT.OR.WATER_QUALITY) THEN
       IF(WATER_QUALITY) THEN

!
!       WATER QUALITY
!
        IF(FILES(ATMFILEA)%NAME(1:1).NE.' ')THEN

          IF(WATER_QUALITY) THEN
!         TIME VARYING WATER QUALITY
            IF(ATMOSEXCH.EQ.0)THEN

              CALL INTERPMETEO2(NWIND,UAIR,VAIR,TA,PATM,NEBU,RAINFALL,
     &                          PVAP,RAY3,AT,UL)
!
              CALL OV('X=C     ', X=WINDX, C=UAIR, DIM1=NPOIN)
              CALL OV('X=C     ', X=WINDY, C=VAIR, DIM1=NPOIN)
              CALL OV('X=C     ', X=PATMOS, C=PATM, DIM1=NPOIN)
              CALL OV('X=C     ', X=TAIR%R, C=TA, DIM1=NPOIN)
              IF(PRESENT(PLUIE))THEN
                CALL OS('X=C     ',X = PLUIE, C=RAINFALL) ! MM/S
              ENDIF
!
!         TIME VARYING WATER QUALITY WITH HEAT EXCHANGE WITH ATMOSPHERE
!            ELSEIF(ATMOSEXCH.EQ.1.OR.ATMOSEXCH.EQ.2) THEN
            ELSEIF(ATMOSEXCH.EQ.2) THEN

              CALL INTERPMETEO(WW,UAIR,VAIR,TA,PATM,
     &                         HREL,NEBU,RAINFALL,EVAPORATION,AT,UL)
              CALL OV('X=C     ',X=WINDX, C=UAIR, DIM1=NPOIN)
              CALL OV('X=C     ',X=WINDY, C=VAIR, DIM1=NPOIN)
              CALL OV('X=C     ',X=PATMOS, C=PATM, DIM1=NPOIN)
              CALL OV('X=C     ',X=TAIR%R, C=TA, DIM1=NPOIN)
!            ENDIF
!
!         SIMPLIFIED HEAT EXCHANGE: TEMPERATURE AND WIND DATA ARE ENOUGH
            ELSEIF(ATMOSEXCH.EQ.1) THEN

              A=2
              DO WHILE (AT_WIND(A).LT.AT)
                A = A+1
              ENDDO
!         INTERPOLATE IN TIME THE CONSTANT TEMPERATURE VALUE
              TA = DATAFROMFILE(A-1, 1+NUMSTA*2+1) +
     &          (DATAFROMFILE(A, 1+NUMSTA*2+1) -
     &           DATAFROMFILE(A-1, 1+NUMSTA*2+1)) /
     &          (AT_WIND(A) - AT_WIND(A-1)) *
     &          (AT - AT_WIND(A-1))
            WRITE(LU, *) 'AT=',AT
            WRITE(LU, *) A
            WRITE(LU, *) 'TEMPERATURE: ', DATAFROMFILE(A, 1+NUMSTA*2+1)
              CALL OV('X=C     ',X=TAIR%R, C=TA, DIM1=NPOIN)
            ENDIF !TYPE OF EXCHANGE

          ENDIF !WATER QUALITY


!          ELSEIF (VENT) THEN
           IF (VENT) THEN
!

          IF(OPTWIND.EQ.3)THEN
!
! #######################################################################
!         IDWM WIND INTERPOLATION CODE
! #######################################################################
!
!
!       ASSEMBLE THE ARRAYS OF X,Y,WNDSPD AND X,Y,WNDDIR FOR EACH ITERATION
            PI = 4.D0*ATAN(1.D0)
            DTR = PI/180.D0

! FIND THE LAST WIND DATA. WILL INTERPOLATE FROM THERE TO THE NEXT
            A = 2
            DO WHILE (AT_WIND(A).LT.AT)
              A = A+1
            ENDDO


            WRITE(LU,*) 'METEO: WIND INTERPOLATED BETWEEN ',
     &            AT_WIND(A-1), ' AND ', AT_WIND(A)
            DO B = 1,NUMSTA
              ! ASSEMBLE THE ARRAYS FOR THIS TIME STEP
              ! DIRECTIONS FROM INPUT FILE
              D_LAST(B) = DATAFROMFILE(A-1,B*2+1)
              D_NEXT(B) = DATAFROMFILE(A,B*2+1)

              ! SPEEDS FROM INPUT FILE
              S_LAST(B) = DATAFROMFILE(A-1,B*2)
              S_NEXT(B) = DATAFROMFILE(A,B*2)

                WRITE(LU,*) 'WINDSPLAST', S_LAST(B)
                WRITE(LU,*) 'WINDDIRLAST', D_LAST(B)

              ! ADDED ON 2016.06.15
              ! CHECK IF WIND SPEED IS +VE, AND IF WIND DIRECTION
              ! IS BETWEEN 0 AND 360 DEG;
              IF (S_LAST(B) < 1.D-6) THEN
                S_LAST(B) = 1.D-6
              END IF

              IF (S_NEXT(B) < 1.D-6) THEN
                S_NEXT(B) = 1.D-6
              END IF

              IF (D_LAST(B) .LT. 1.D-6) THEN
                D_LAST(B) = 1.D-6
              END IF

              IF (D_LAST(B) .GT. 360.D0) THEN
                D_LAST(B) = 360.D0
              END IF

              IF (D_NEXT(B) .LT. 1.D-6) THEN
                D_NEXT(B) = 1.D-6
              END IF

              IF (D_NEXT(B) .GT. 360.D0) THEN
                D_NEXT(B) = 360.D0
              END IF

              ! ADDED ON 2016.05.26
              ! RATHER THAN INTERPOLATING THE DIRECTION VARIABLE
              ! CONVERT INPSTA_D TO ITS X AND Y COMPONENTS,
              ! INTERPOLATE BOTH

              ! THIS IS NEEDED BECAUSE INTERPOLATING A DIRECTION
              ! VARIABLE THAT TAKES ON VALUES BETWEEN 0 AND 360 HAS
              ! PROBLEMS WHEN INTERPOLATING NODES LOCATED CLOSE TO
              ! STATIONS WITH  DIR~0'S AND DIR~350'S


              ! CONVERT D_ AND S_ TO WINDX AND WINDY
              ! THESE ARE CARTESIAN VECTORS OF THE WIND
            IF (D_LAST(B) >= 0.D0 .AND. S_LAST(B) >= 0.D0) THEN
              IF ((D_LAST(B) >= 0.D0) .AND. (D_LAST(B) <= 90.D0)) THEN
                  THETA_RAD = D_LAST(B) * DTR
                    WINDX_LAST(B) = -SIN(THETA_RAD)*S_LAST(B)
                    WINDY_LAST(B) = -COS(THETA_RAD)*S_LAST(B)
              END IF
!
              IF ((D_LAST(B) > 90.D0) .AND. (D_LAST(B) <= 180.D0)) THEN
                  THETA_RAD = (180.D0 - D_LAST(B)) * DTR
                    WINDX_LAST(B) = -SIN(THETA_RAD)*S_LAST(B)
                    WINDY_LAST(B) =  COS(THETA_RAD)*S_LAST(B)
              END IF
!
              IF ((D_LAST(B) > 180.D0) .AND. (D_LAST(B) <= 270.D0)) THEN
                  THETA_RAD = (D_LAST(B)-180.D0) * DTR
                    WINDX_LAST(B) = SIN(THETA_RAD)*S_LAST(B)
                    WINDY_LAST(B) = COS(THETA_RAD)*S_LAST(B)
              END IF
!
              IF ((D_LAST(B) > 270.D0) .AND. (D_LAST(B) <= 360.D0)) THEN
                  THETA_RAD = (360.D0-D_LAST(B)) * DTR
                    WINDX_LAST(B) =  SIN(THETA_RAD)*S_LAST(B)
                    WINDY_LAST(B) = -COS(THETA_RAD)*S_LAST(B)
              END IF
            ELSE
                WINDX_LAST(B) = -999.D0
                WINDY_LAST(B) = -999.D0
            ENDIF


            IF (D_NEXT(B) >= 0.D0 .AND. S_NEXT(B) >= 0.D0) THEN
              IF ((D_NEXT(B) >= 0.D0) .AND. (D_NEXT(B) <= 90.D0)) THEN
                  THETA_RAD = D_NEXT(B) * DTR
                    WINDX_NEXT(B) = -SIN(THETA_RAD)*S_NEXT(B)
                    WINDY_NEXT(B) = -COS(THETA_RAD)*S_NEXT(B)
              END IF
!
              IF ((D_NEXT(B) > 90.D0) .AND. (D_NEXT(B) <= 180.D0)) THEN
                  THETA_RAD = (180.D0 - D_NEXT(B)) * DTR
                    WINDX_NEXT(B) = -SIN(THETA_RAD)*S_NEXT(B)
                    WINDY_NEXT(B) =  COS(THETA_RAD)*S_NEXT(B)
              END IF
!
              IF ((D_NEXT(B) > 180.D0) .AND. (D_NEXT(B) <= 270.D0)) THEN
                  THETA_RAD = (D_NEXT(B)-180.D0) * DTR
                    WINDX_NEXT(B) = SIN(THETA_RAD)*S_NEXT(B)
                    WINDY_NEXT(B) = COS(THETA_RAD)*S_NEXT(B)
              END IF
!
              IF ((D_NEXT(B) > 270.D0) .AND. (D_NEXT(B) <= 360.D0)) THEN
                  THETA_RAD = (360.D0-D_NEXT(B)) * DTR
                    WINDX_NEXT(B) =  SIN(THETA_RAD)*S_NEXT(B)
                    WINDY_NEXT(B) = -COS(THETA_RAD)*S_NEXT(B)
              END IF
            ELSE
                WINDX_NEXT(B) = -999.D0
                WINDY_NEXT(B) = -999.D0
            ENDIF

              IF (WINDX_LAST(B).EQ.-999.D0.OR.WINDY_LAST(B).EQ.-999.D0.
     &            OR.WINDX_NEXT(B).EQ.-999.D0.OR.
     &            WINDY_NEXT(B).EQ.-999.D0)
     &             THEN
                INPSTA_WINDX(B,1) = -999.D0
                INPSTA_WINDX(B,2) = -999.D0
              ENDIF


              ! ASSIGN X AND Y COORDINATES
              INPSTA_WINDX(B,1) = XX(B)
              INPSTA_WINDX(B,2) = YY(B)

              INPSTA_WINDY(B,1) = XX(B)
              INPSTA_WINDY(B,2) = YY(B)

              ! INTERPOLATE IN TIME
              INPSTA_WINDX(B,3) = WINDX_LAST(B) +
     &         (WINDX_NEXT(B)-WINDX_LAST(B))/(AT_WIND(A)-AT_WIND(A-1))*
     &         (AT - AT_WIND(A-1))

              INPSTA_WINDY(B,3) = WINDY_LAST(B) +
     &         (WINDY_NEXT(B)-WINDY_LAST(B))/(AT_WIND(A)-AT_WIND(A-1))*
     &         (AT - AT_WIND(A-1))

        !    WRITE(LU,*) 'WINDX', INPSTA_WINDX(B,3)
        !    WRITE(LU,*) 'WINDY', INPSTA_WINDY(B,3)

            ENDDO ! B


            WRITE(LU,*) "CALLING THE INTERPOLATOR"
            CALL RBF_2D(INPSTA_WINDX,POINTS,OUT_WDIRX,NPOIN,NUMSTA)
            CALL RBF_2D(INPSTA_WINDY,POINTS,OUT_WDIRY,NPOIN,NUMSTA)


!       FINAL WINDX AND WINDY OUTPUT
            DO K = 1,NPOIN
              WINDX(K) = OUT_WDIRX(K)
              WINDY(K) = OUT_WDIRY(K)
            END DO
!
! #######################################################################
!
          ENDIF ! OPTWIND.EQ.3
          ENDIF
        ENDIF
!
!       WIND AND/OR WATER QUALITY VARIABLES
!       VARYING IN SPACE AND TIME, FROM A BINARY FILE
!
        IF(FILES(ATMFILEB)%NAME(1:1).NE.' ') THEN
          IF(FILES(ATMFILEA)%NAME(1:1).NE.' ')THEN
            WRITE(LU,*) 'METEO: THE DATA FROM THE ASCII METEO'
            WRITE(LU,*) 'FILE WILL BE OVERWRITTEN BY THE'
            WRITE(LU,*) 'CORRESPONDING BINARY FILE DATA'
          ENDIF
          CALL METEO_FROM_BINARY_FILE(PATMOS,WINDX,WINDY,AT,NPOIN,VENT,
     &                                ATMOS,ATMFILEB,FILES,LISTIN,
     &                                OPTWIND,WIND_SPD)
        ENDIF
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
