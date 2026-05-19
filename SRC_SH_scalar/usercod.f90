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
  REAL(RK)    :: A_CONSTANT,B_CONSTANT,C_CONSTANT

  IF(LHYDROSTATIC) THEN    !!! START FROM HYSTROSATIC SOLUTONS

     IF(abs(CINJECTION-0.1)<tiny(1.)) THEN
        A_CONSTANT = 0.4363        !! C = 0.1   WEAK INJECTION
     ELSEIF(abs(CINJECTION-1.0)<tiny(1.)) THEN
        A_CONSTANT = 1.1256        !! C = 1.0   MODERATE INJECTION
     ELSEIF(abs(CINJECTION-10.0)<tiny(1.)) THEN
        A_CONSTANT = 1.4882        !! C = 10.0  STRONG INJECTION
     ELSE
        WRITE(*, *) 'HYDROSTATIC SOLUTION IS NOT AVAIABLE!'
        WRITE(*, *) 'CHECK C = 0.1, 1.0 OR 10.0 OR NOT  ! '
        STOP
     END IF

     B_CONSTANT= (A_CONSTANT/(2.0*CINJECTION))**2
     C_CONSTANT = 1.0+(2.0/3.0)*A_CONSTANT*(B_CONSTANT**1.5)

     U=0.
     V=0.
     W=0.

!........WE ASSUME THAT HEATING IS FROM BELOW SIDE AND
!........IONIC MOBILITY AND PERMITTIVITY IS TEMPERATURE-INDEPEDENT

     IF(LCAL(IEN)) T=1.0-YC

!........ASSUMING THAT INJECTION IS FROM BELOW

     VP=C_CONSTANT-(2.0/3.0)*A_CONSTANT &
          *((B_CONSTANT+YC)**1.5)
     Q=A_CONSTANT/(2.0*CINJECTION &
          *SQRT(YC+B_CONSTANT))

  ELSE     !!! START FROM ZERO-INITIAL FIELDS

     U=UIN
     V=VIN
     W=WIN
     IF(LCAL(IEN)) T=TIN
     VP=0.
     Q=0.

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

!
!#####################################

!
!.....ALL BLOCKS
  DO NBL=1,NBLOCK
     CALL SETIND(NBL)
     DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL) !!! BOUNDARY REGION
!
        NBC_S=NUM_SPR(IDR)+1
        NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)
!
        SELECT CASE(NUM_TYP(IDR))

         CASE(4)     !!! 4 WALL 

            DO I=NBC_S,NBC_E 
              INP=NBL_ST(NBL)+NUM_IND(I,1)
              ISIDE=NUM_IND(I,2)
!
              CALL INIT_GEOM_BF(ISIDE,INP)
             
              U(INE)=0.0
              V(INE)=0.0
              W(INE)=0.0
              Q(INE)=0.0

            ENDDO

          CASE(10)     !! ! REGION  WITH BC 10 -- Exposed electrode  
!    
            DO I=NBC_S,NBC_E 
              INP=NBL_ST(NBL)+NUM_IND(I,1)
              ISIDE=NUM_IND(I,2)
!
              CALL INIT_GEOM_BF(ISIDE,INP)
!       
              U(INE)=0.0
              V(INE)=0.0
              W(INE)=0.0
              VP(INE)=1.0
        
           ENDDO

        CASE(11)     !!! Grounded electrode -- Zero potential BC
!
           DO I=NBC_S,NBC_E
              INP=NBL_ST(NBL)+NUM_IND(I,1)
              ISIDE=NUM_IND(I,2)
!
              CALL INIT_GEOM_BF(ISIDE,INP)

              VP(INE) = 0.0     

           ENDDO
         
        CASE(12)     !!! Surface charge density- Gaussian profile 
!
           DO I=NBC_S,NBC_E
              INP=NBL_ST(NBL)+NUM_IND(I,1)
              ISIDE=NUM_IND(I,2)
!
              CALL INIT_GEOM_BF(ISIDE,INP)
!
              U(INE)=0.0
              V(INE)=0.0
              W(INE)=0.0
              Q(INE) = EXP(-((XC(INE)-MuLOCATION)**2)*(0.5*(sigmaR)**2)) 
    
! Q normalised with -- Qmax* f(t)  where f(t)=sin(2*PI*F_BASE*t)
!complete expression for Q on BC12 is Q = Qmax * f(t)*EXP(-((XC(INP)-MuLOCATION)**2)*(0.5*(sigmaR)**2))

           ENDDO
        
         CASE(14)     !! ! WALL for VELOCITIES
!    
            DO I=NBC_S,NBC_E 
              INP=NBL_ST(NBL)+NUM_IND(I,1)
              ISIDE=NUM_IND(I,2)
!
              CALL INIT_GEOM_BF(ISIDE,INP)
!       
              U(INE)=0.0
              V(INE)=0.0
              W(INE)=0.0
                      
           ENDDO


      END SELECT
!
!.....END REGIONS                   
      ENDDO
!
!.....END BLOCKS
  ENDDO
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
!
!..... HERE SHOULD BE GIVEN THE VALUES FOR VARIABLES
!.....   (BOUNDARY  CONDITION)  FI=FI(R,TIME) OR GRAD FI=GRAD FI(R,TIME)
!..... AT OUTLET B.C. TYP  2 (TAKE CARE THAT YOU SET THE VALUE OF VARIABLES IF "INLET"
!.....                                      OCCURS AT "OUTLET")
!
!
!.....ALL BLOCKS
  DO NBL=1,3
     CALL SETIND(NBL)
!
!.....SET FOR  DIFFERENT REGIONS
     DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)
!
        IF(NUM_TYP(IDR)==2) THEN
!*********************************************************************
!..... SET VALUES AT OUTLET BOUNDARIES ...   T Y P    2              *
!..... HERE ONLY FOR SCALAR(S), SEE ALSO  OUTBC.F                    *
!*********************************************************************
!.....................................................................
!..... SET IT INDEPENDENT OF THE DIFFERENT REGION(S)                 .
!..... IT MEANS, FOR ALL "OUTLET" REGION(S)  THE SAME CONDITION      .
!...... BUT IF YOU NEED SOMETHING ELSE, YOU CAN CHANGE IT            .
!.....................................................................
           NBC_S=NUM_SPR(IDR)+1
           NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)
!
           DO I=NBC_S,NBC_E
              INP=NBL_ST(NBL)+NUM_IND(I,1)
              ISIDE=NUM_IND(I,2)
!
              CALL INIT_GEOM_BF(ISIDE,INP)
!
              DEN(INE)  =  DEN(INP)
              VIS(INE)  =  VIS(INP)
             
            
              IF(LCAL(IEN)) T(INE)  =    T(INP)
!
           ENDDO
!.....END  T Y P  OUTLET
        ENDIF
!
!.....END REGIONS
     ENDDO
!
!.....END BLOCKS
  ENDDO
  RETURN
END SUBROUTINE SETOBC
END MODULE USERCOD

