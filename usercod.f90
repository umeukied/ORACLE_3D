!*********************************************************************
!> This routine initializes some parameters, field and other stuffs
!>
!> @author Philippe Traore
!*********************************************************************
MODULE USERCOD
IMPLICIT NONE
CONTAINS

SUBROUTINE INIT_USER
  USE DEC
  IMPLICIT NONE
  REAL(RK) :: A_CONSTANT,B_CONSTANT,C_CONSTANT
  REAL(RK), PARAMETER :: TOL = 1.0E-12_RK

  IF (LHYDROSTATIC) THEN    !!! START FROM HYDROSTATIC SOLUTIONS

     IF (ABS(CINJECTION-0.1_RK) < TOL) THEN
        A_CONSTANT = 0.4363_RK        !! C = 0.1   WEAK INJECTION
     ELSEIF (ABS(CINJECTION-1.0_RK) < TOL) THEN
        A_CONSTANT = 1.1256_RK        !! C = 1.0   MODERATE INJECTION
     ELSEIF (ABS(CINJECTION-10.0_RK) < TOL) THEN
        A_CONSTANT = 1.4882_RK        !! C = 10.0  STRONG INJECTION
     ELSE
        WRITE(*,*) 'HYDROSTATIC SOLUTION IS NOT AVAILABLE!'
        WRITE(*,*) 'CHECK C = 0.1, 1.0 OR 10.0'
        STOP
     END IF

     B_CONSTANT = (A_CONSTANT/(2.0_RK*CINJECTION))**2
     C_CONSTANT = 1.0_RK + (2.0_RK/3.0_RK)*A_CONSTANT*(B_CONSTANT**1.5_RK)

     U = 0.0_RK
     V = 0.0_RK
     W = 0.0_RK

!........WE ASSUME THAT HEATING IS FROM BELOW SIDE AND
!........IONIC MOBILITY AND PERMITTIVITY IS TEMPERATURE-INDEPENDENT

     IF (LCAL(IEN)) T = 1.0_RK - YC

!........ASSUMING THAT INJECTION IS FROM BELOW

     VP = C_CONSTANT - (2.0_RK/3.0_RK)*A_CONSTANT*((B_CONSTANT+YC)**1.5_RK)
     Q  = A_CONSTANT/(2.0_RK*CINJECTION*SQRT(YC+B_CONSTANT))

  ELSE     !!! START FROM ZERO / UNIFORM INITIAL FIELDS

     U = UIN
     V = VIN
     W = WIN

     IF (LCAL(IEN)) T = TIN

!.... FOR HYDRO-ONLY CASE: FORCE ELECTRIC VARIABLES TO ZERO
     VP = 0.0_RK
     Q  = 0.0_RK

  END IF

!.... EXTRA SAFETY FOR YOUR CASE: IF EHD IS DISABLED, KEEP EVERYTHING ELECTRIC AT ZERO
  IF ((.NOT. LEHD) .AND. (.NOT. LETHD)) THEN
     VP = 0.0_RK
     Q  = 0.0_RK
  END IF

  RETURN
END SUBROUTINE INIT_USER


!*********************************************************************
!> This routine set the boundary conditions
!>
!> @author Philippe Traore
!*********************************************************************
SUBROUTINE SETBC
  USE DEC
  USE GEOM
  IMPLICIT NONE
  INTEGER(IK) :: NBL,IDR,NBC_S,NBC_E,I,INP,ISIDE

!.....ALL BLOCKS
  DO NBL=1,NBLOCK
     CALL SETIND(NBL)

     DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)   !!! BOUNDARY REGION

        NBC_S=NUM_SPR(IDR)+1
        NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)

        SELECT CASE(NUM_TYP(IDR))

        CASE(1)     !!! INLET

           DO I=NBC_S,NBC_E
              INP   = NBL_ST(NBL)+NUM_IND(I,1)
              ISIDE = NUM_IND(I,2)

              CALL INIT_GEOM_BF(ISIDE,INP)

              U(INE) = 0.0_RK
              V(INE) = -1.0_RK
              W(INE) = 0.0_RK
              
              
              FLX_IND(I)=DEN(INP)*(U(INE)*ARX+V(INE)*ARY+W(INE)*ARZ)
              
              

!.....HYDRO ONLY
              Q(INE)  = 0.0_RK
              VP(INE) = 0.0_RK
           END DO

        CASE(4)     !!! WALL

           DO I=NBC_S,NBC_E
              INP   = NBL_ST(NBL)+NUM_IND(I,1)
              ISIDE = NUM_IND(I,2)

              CALL INIT_GEOM_BF(ISIDE,INP)

              U(INE)  = 0.0_RK
              V(INE)  = 0.0_RK
              W(INE)  = 0.0_RK
              Q(INE)  = 0.0_RK
              VP(INE) = 0.0_RK
           END DO

        CASE(10)    !!! EXPOSED ELECTRODE --> DISABLED IN HYDRO CASE

           DO I=NBC_S,NBC_E
              INP   = NBL_ST(NBL)+NUM_IND(I,1)
              ISIDE = NUM_IND(I,2)

              CALL INIT_GEOM_BF(ISIDE,INP)

              U(INE)  = 0.0_RK
              V(INE)  = 0.0_RK
              W(INE)  = 0.0_RK
              Q(INE)  = 0.0_RK
              VP(INE) = 0.0_RK
           END DO

        CASE(11)    !!! GROUNDED ELECTRODE --> DISABLED IN HYDRO CASE

           DO I=NBC_S,NBC_E
              INP   = NBL_ST(NBL)+NUM_IND(I,1)
              ISIDE = NUM_IND(I,2)

              CALL INIT_GEOM_BF(ISIDE,INP)

              Q(INE)  = 0.0_RK
              VP(INE) = 0.0_RK
           END DO

        CASE(12)    !!! SURFACE CHARGE DENSITY --> DISABLED IN HYDRO CASE

           DO I=NBC_S,NBC_E
              INP   = NBL_ST(NBL)+NUM_IND(I,1)
              ISIDE = NUM_IND(I,2)

              CALL INIT_GEOM_BF(ISIDE,INP)

              U(INE)  = 0.0_RK
              V(INE)  = 0.0_RK
              W(INE)  = 0.0_RK
              Q(INE)  = 0.0_RK
              VP(INE) = 0.0_RK
           END DO

        CASE(14)    !!! WALL FOR VELOCITIES / FAR WALL

           DO I=NBC_S,NBC_E
              INP   = NBL_ST(NBL)+NUM_IND(I,1)
              ISIDE = NUM_IND(I,2)

              CALL INIT_GEOM_BF(ISIDE,INP)

              U(INE) = 0.0_RK
              V(INE) = 0.0_RK
              W(INE) = 0.0_RK
           END DO

        END SELECT

     END DO
  END DO

  RETURN
END SUBROUTINE SETBC


!*********************************************************************
!> This routine set the outlet boundary conditions
!>
!> @author Philippe Traore
!*********************************************************************
SUBROUTINE SETOBC
  USE DEC
  USE GEOM
  IMPLICIT NONE
  INTEGER(IK) :: NBL,IDR,NBC_S,NBC_E,I,INP,ISIDE

!.....ALL BLOCKS
  DO NBL=1,NBLOCK
     CALL SETIND(NBL)

     DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)

        IF (NUM_TYP(IDR)==2) THEN

           NBC_S=NUM_SPR(IDR)+1
           NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)

           DO I=NBC_S,NBC_E
              INP   = NBL_ST(NBL)+NUM_IND(I,1)
              ISIDE = NUM_IND(I,2)

              CALL INIT_GEOM_BF(ISIDE,INP)

              DEN(INE) = DEN(INP)
              VIS(INE) = VIS(INP)

              IF (LCAL(IEN)) T(INE) = T(INP)

!.....HYDRO ONLY
              Q(INE)  = 0.0_RK
              VP(INE) = 0.0_RK
           END DO

        END IF

     END DO
  END DO

  RETURN
END SUBROUTINE SETOBC

END MODULE USERCOD
