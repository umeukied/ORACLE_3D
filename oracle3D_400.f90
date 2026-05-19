!******************************************************
!**==================================================**
!                                                    **
!                   ORACLE3D                         **
!                                                    **
!    = = =    VERSION       # 4.00   #   = = =       **
!                                                    **
!   3-DIMENSIONAL COMPUTATIONAL FLUID DYNAMICS CODE  **
!                                                    **
!                                                    **
!   MAIN FEATURES:                                   **
!                  FINITE VOLUME                     **
!                  NON ORTHOGONAL CV                 **
!                  STRUCTURED MULTI-BLOCK            **
!                  COLOCATED  (RHIE & CHOW)          **
!                                                    **
!**==================================================**
!******************************************************
PROGRAM ORACLE3D_4
!******************************************************
!     NUMERICAL SOLVER FOR 3D ELECTRO-CONVECTION &    C
!     ELECTRO-THERMO-CONVECTION DUE TO UNIPOLAR       C
!     CHARGE INJECTION AND THERMAL GRADIENT           C
!                                                     C
!     NON-DIMENSIONAL STYLES, REFER TO:               C
!         PHILIPPE - POF 2012                         C
!         PHILIPPE - JFM 2010                         C
!                                                     C
!     VERSION: 2013-FEB-7TH                           C
!                                                     C
!     UPDATE FOR IVF METHOD _3D,  4TH-JUNE-2013       C
!     UPDATE FOR SOLVING NS _3D, 15TH-JULY-2013       C
!     UPDATE FOR PARALLELIZATION,1ST-APRIL-2014       C
!**==================================================**
  USE DEC
  USE USERCOD
  IMPLICIT NONE

  INTEGER(IK)         :: I
  REAL(RK)            :: SOURCEO
  REAL(RK)            :: TIME_1,TIME_2,TMP


	
!.....INIT DATA
  CALL  INITDATA

!.....DEFINE SOME CONSTANTS
    CALL MODCON

!.....INPUT DATA
  CALL SETDATA
	
!.....OPEN FILES
  CALL OPFILE

!.....INITIATE GRID
  CALL READGRID

!..... COORDINATES ON BLOCK INTERFACES SAVED
  CALL SAVCIN

!......CALCULATE METRIC AND STUFFS
  CALL METRIC

!.....TO CONVERT TO NON-DIMENSIONAL COMPUTATIONS
  CALL NONDIMENSION_FORM

!.....INITIALISATION
  CALL INIT

!.....CONTINUE COMPUTATION FROM PREVIOUS SOLUTIONS
  IF(LREAD) CALL RREST

!================================================================
!     START  =>
!================================================================

  CALL NEWFLX
  PP=P

!---------------------------------------------------------

!.....SET INITIAL BOUNDARY CONDITIONS
  CALL SETBC         !=>   INLETS, WALLS
  CALL SETOBC        !=>   OUTLETS
!
 

!.....CALCULATE NORMALISATION FACTORS
  CALL CALNOR
!
!.....INITIAL OUTPUT
  CALL MOOUT1

!
  WRITE(*,611) NBLOCK,NCVALL
  WRITE(*,600) NBLMON,IMON,JMON,KMON

   IF((LEHD.OR.LETHD).AND.ABS(D_C)>SMALL) THEN
 !      IF(LESP) CALL ESP_UPDATE(NXYZA,ESP_P,ESP_K,RA,T,L_PERMIT,N_MOBILITY)
        CALL ES_MODULE
        CALL WRUNST2PLT
     END IF
  

!.....START TIME LOOP
!
!.....TIME DEPEND ONLY IF LTIME=.TRUE.
  IF(.NOT.LTIME) NTM=1

!++++++++++++++++++++++++++++
!     TIME LOOP --> BEGIN
!++++++++++++++++++++++++++++
  call get_time(TIME_1)
  TIME_LOOP: DO ITIM=1,NTM

IF(ITIM==2)   call get_time(TIME_1)
     TIME=TIME+DT
     IF(LTIME) THEN
        WRITE(*,617) ITIM,TIME
!..........UPDATE PREVIOUS TIME LEVELS
        CALL UPDOLD
     END IF
!
!===================================================C
!                                                   C
!      EHD & ETHD =>  ELECTROSTATICS                C
!        * CHARGE DENSITY                           C
!        * ELECTRICAL POTENTIAL                     C
!        * ELECTRIC FIELD                           C
!                                                   C
!     LEHD:  TO COMPUTE ELECTRO-CONVECTION  ?       C
!     LETHD: TO COMPUTE ELECTRO-THERMO-CONVECTION ? C
!     LESP:  TEMPERATURE DEPENDENT PERMITTIVITY?    C
!                                                   C
!===================================================C
!
    
  

     ITER=0
     SOURCEO = 0.

!==============================================
!     START => INNER ITERATIONS LOOP (SIMPLE)
!==============================================
!
     SIMPLE_LOOP: DO LS=1,NSWEEP

        ITER=ITER+1

        IF(LCAL(IU))    CALL CALUVW
        IF(LCAL(IU))    CALL OUTBC
        IF(LCAL(IP))    CALL CALCP
        IF(LCAL(IEN))   CALL CALCSC(IEN,T,TO,TOO)

        CALL SETOBC     !=>  OUTLET BOUNDARY CONDITION
!
!.....RESIDUAL NORMALIZATION, CONVERGENCE CHECK
!
        RESOR(1:NPHI)=RESOR(1:NPHI)*SNORIN(1:NPHI)

        IF(I_SC==1) THEN  !.....PERIC CRITERION BASED ON RESIDUAL
           SOURCE=MAX(RESOR(IU) ,RESOR(IV),RESOR(IW) ) ! , & RESOR(IEN),RESOR(IP))
                      
        ELSE                !.....DATE'S RM2 STOPPING CRITERIION
           CALL DIVERGENCE(PP)
           SOURCE=RMM2
        END IF

        TMP=0.
        IF(LCAL(IEN)) TMP=T(IJKMON)
        IF (LSCREEN) WRITE(*,601) MOD(ITER,1000_ik), &
             (RESOR(I),I=1,5),U(IJKMON),          &
             V(IJKMON),     W(IJKMON), P(IJKMON), &
             TMP,Q(IJKMON),VP(IJKMON), &
             SOURCE

!========================================
!     RESIDUALS ARE TOO LARGE           =
!     TERMINATE PROGRAM WITH MESSAGE    =
!========================================

        IF(SOURCE>SLARGE) THEN
           WRITE(*,602)
           CALL END
        END IF

!=============================================
!     RESIDUALS IS SMLLER THAN GIVEN         =
!     CRITERIA --> JUMP OUT THE SIMPLE LOOP  =
!     END =>  OUTER ITERATIONS  LOOP         =
!=============================================

!        IF((SOURCE<SORMAX).OR.(ABS(SOURCE-SOURCEO)<SORMAX2)) THEN

!	write(*,*) source , sormax
	
	
        IF(SOURCE<SORMAX) THEN
           WRITE(*,605)
           EXIT
        END IF

        SOURCEO = SOURCE

     END DO SIMPLE_LOOP

!
!================================================
!     =>  THE ENERGY EQUATION (T)
!     =>  OUT OF SIMPLE LOOP
!================================================
!
      IF(LCAL(IEN)) THEN
        CALL CALCSC(IEN, T, TO, TOO)
      END IF


!================================================
!      STORE UNSTEADY RESULTS FOR POST PROCESSING
!================================================

     IF(LSTORE.AND.LTIME.and.MOD(ITIM,IMA)==0) THEN
 !          IF(ITIM==10025.OR.ITIM==13025.OR.ITIM==15025.OR.ITIM==17025.OR.ITIM==20025) THEN

                WRITE(*,606)
                CALL WRUNST2PLT
  

!      ==================================
!      CALCUL DES CRITERES TOURBILLONNAIRES
!      ==================================
                IF (LVORTEX) CALL WRITE_VORTEX_FIELD
!      END IF

!==================================
!      SAVE SOLUTIONS FOR RESTART
!==================================
        IF(LWRITE) CALL WREST

!==================================
!      OUTPUT THE NUSSELT NUMBER
!==================================

        IF (LNUSSELT) THEN
           CALL NUSSELT_COMPUTATION                   ! ==> FX IS REQUIRED FOR GEOM.INC
             WRITE(103,*) 'NU_X_0.0 : ', R_NU_0       ! ==> NUSSELT OUTPUT REQUIRES TO VALIFY
             WRITE(103,*) 'NU_X_0.5 : ', R_NU_05
             WRITE(103,*) 'NU_X_1.0 : ', R_NU_10
        END IF

     END IF

!===================================================
!      OUTPUT THE INFORMATION AT MONITOR POINT
!===================================================
     TMP=0.
     IF(LCAL(IEN)) TMP=T(IJKMON)
     WRITE(101,*) TIME, U(IJKMON), V(IJKMON),W(IJKMON),  &
                        P(IJKMON), TMP,       &
                        Q(IJKMON), VP(IJKMON)

!===================================================
!      OUTPUT THE MAXIMUM VELOCITY INFORMATION
!===================================================

     IF (LUVWMAX) THEN
        CALL FIND_UVWMAX
        WRITE(102,'(5(2x,E23.15))') TIME, UMAX, VMAX, WMAX, UVW_NORM_MAX
        WRITE(103,*) 'TIME ------------>', TIME
        WRITE(103,*) 'NORM', UVW_NORM_MAX, XMAX_NORM, &
                            YMAX_NORM,     ZMAX_NORM
        WRITE(103,*) 'U', UMAX, XMAX_U, YMAX_U, ZMAX_U
        WRITE(103,*) 'V', VMAX, XMAX_V, YMAX_V, ZMAX_V
        WRITE(103,*) 'W', WMAX, XMAX_W, YMAX_W, ZMAX_W
     END IF

!.....IMPLICIT THREE TIME LEVEL SCHEME
     IF (ITIME_U==1) GAMT_U = 1.0
     IF (ITIME_T==1) GAMT_T = 1.0
     IF (ITIME_Q==1) GAMT_Q = 1.0

!++++++++++++++++++++++++++++
!     TIME LOOP --> END
!++++++++++++++++++++++++++++

  END DO TIME_LOOP
  call get_time(TIME_2)
PRINT*,'Integration time : ',TIME_2-TIME_1

  close(IUCHK)
  
  IF(LSTORE)THEN
    CLOSE( 78)
    CLOSE(101)
    CLOSE(102)
    CLOSE(103)

    CLOSE(IURESX)
    CLOSE(IURESY)
    CLOSE(IURESZ)

    CLOSE(85)
    CLOSE(86)
    CLOSE(87)

    IF (LNUSSELT) THEN
       CLOSE(104)
       CLOSE(106)
    END IF
    IF (LVORTEX) CLOSE(109)
  endif
  

  CALL END


  WRITE(*,*) 'WELL DONE !!'

!.............................................................
!.....        F O R M A T S    ...............................
!.............................................................
!.....
600 FORMAT(/,1X,' IT ',                                     &
       'I---------ABSOLUTE RESIDUAL SOURCE SUMS',           &
       '----------I',10X,                                   &
       'I----FIELD VALUES AT MONITORING LOCATION BLOCK # ', &
       I3,' (',I3,',',I3,',',I3,')-----I',/,                &
       4X,'NO ',2X,'UMOM',5X,'VMOM',5X,'WMOM',              &
       5X,'MASS',5X,'ENER',                                 &
       20X,'U',9X,'V',9X,'W',9X,'P',9X,'T',/)
601 FORMAT(1X,I5,1P5E10.2,12X,1P8E10.2)


602 FORMAT(//,10X,'*** PROGRAM TERMINATED  -', &
                  'OUTER ITERATIONS DIVERGE ***')

604 FORMAT(//,10X,'***CONVERGENCE CRITERION  IS NOT SATISFIED***')
605 FORMAT(//,10X,'***CONVERGENCE CRITERION   IS SATISFIED***')
606 FORMAT(/,10X,'WRITE RESULTS OF CALCULATIONS IN FILE ',A10)
615 FORMAT(10X,'*** #SPACE# CONVERGENCE CRITERION  IS NOT ', &
               'SATISFIED*** ITER= ',I5)
619 FORMAT(10X,'*** #SPACE# CONVERGENCE CRITERION   IS SATISFIED***', &
               ' ITER= ',I5)
611 FORMAT(/,30X,60('*'),/,30X,                                    &
       '* ****',13(' '),13(' '),'**** *',/,                        &
       30X,60('*'),/,30X,'* ',                                     &
       '  NO OF BLOCK(S)= ',I4,' ***',2X,'NO OF CV = ',I6,12X,'*', &
       /,30X,60('*'))
617 FORMAT(/,10X,'==> TIME STEP ',I5,'   TIME : ',1PE11.4,' SEC')
620 FORMAT(//,20X,'*************************************************', &
            /,20X,'******************  STEADY STATE  ***************', &
            /,20X,'*** #TIME# CONVERGENCE CRITERION  IS SATISFIED***', &
            /,20X,'*************************************************')
622 FORMAT(//,20X,'*************************************************', &
            /,20X,'********   CAN NOT FIND STEADY STATE  ***********', &
            /,20X,'* #TIME# CONVERGENCE CRITERION IS NOT SATISFIED**', &
            /,20X,'*************************************************')
!
CONTAINS

!******************************
  SUBROUTINE CALUVW      !*
!******************************
    IMPLICIT NONE

    INTEGER(IK)  :: NBL
    REAL(RK),ALLOCATABLE     :: SEHD_C(:),SETHD_D(:),APT(:)
    REAL(RK)                 :: URFRS,URFMS

    ALLOCATE(SEHD_C(NXYZA),SETHD_D(NXYZA),APT(NXYZA))

    SU=0. ;    SV=0. ;    SW=0.
    SPU=0.;    SPV=0.;    SPW=0.

    CALL GRADFI(GR1X,GR1Y,GR1Z,U)
    CALL GRADFI(GR2X,GR2Y,GR2Z,V)
    CALL GRADFI(GR3X,GR3Y,GR3Z,W)

!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
     CALL CELUVW(NI-2,NJM,NKM,NJ,1,NIJ,FX,AE,AW,F1,2)
     CALL CELUVW(NIM,NJ-2,NKM,1,NIJ,NJ,FY,AN,AS,F2,1)
     CALL CELUVW(NIM,NJM,NK-2,NIJ,NJ,1,FZ,AT,AB,F3,3)

       CALL MODUVW(NBL)

    END DO     ! END BLOCK LOOP

    !! BUOYANCY FORCE(Y DIRECTION)
    IF(LCAL(IEN).OR.LETHD) THEN                 !! PHILIPPE JFM - 2010
       SV=SV+RA*T*VOL/PRANL                     !! ATTENTION TO THE SIGN HERE
    END IF

!.....PART-1 : COULOMB FORCE
    !! IF(LEHD.OR.LETHD) THEN
!......SEHD_C=Q*VOL*CINJECTION*(R*MSTABILITY)**2
!...... Suzen-Huang model Source term.......for dimensional case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Implantation de D_C 
     !! SEHD_C=D_C*Q*VOL*Qmax*VPmax*((SIN(2.*PI*F_BASE*TIME))**2)
      !! SEHD_C=D_C*Q*VOL*((SIN(2.*PI*F_BASE*TIME))**2)
     
       !! SU=SU+SEHD_C*EX
      !!  SV=SV+SEHD_C*EY
       !! SW=SW+SEHD_C*EZ
  
 !! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$      
 !      .PART-1 : COULOMB FORCE
    IF(LEHD.OR.LETHD) THEN
      SEHD_C = 0.0_RK
    ENDIF
	 

!.....PART-2 : DIELECTRIC FORCE
!     C1: NON-ISONTHERMAL FLUID
!     C2: TEMPERATURE DEPENDENT PERMITTIVITY
!
!     EPS_P: PERMITTIVITY VARIBALE

    IF(LETHD.AND.LESP) THEN
       CALL GRADFI(GR1X,GR1Y,GR1Z,ESP_P)
       SETHD_D=-0.5*R*TSTABILITY*VOL* &
              (EX**2+EY**2+EZ**2)
       SU=SU+SETHD_D*GR1X
       SV=SV+SETHD_D*GR1Y
       SW=SW+SETHD_D*GR1Z
    END IF

!.....EXTRAPOLATE PRESSURE

    CALL BPRES(P)
    CALL GRADFI(GR1X,GR1Y,GR1Z,P)

    SU=SU - GR1X*VOL
    SV=SV - GR1Y*VOL
    SW=SW - GR1Z*VOL
    IF(LTIME) THEN
       APT=DEN*VOL*DTR
       SU=SU+APT*((1.+GAMT_U)*UO-0.5*GAMT_U*UOO)
       SV=SV+APT*((1.+GAMT_U)*VO-0.5*GAMT_U*VOO)
       SW=SW+APT*((1.+GAMT_U)*WO-0.5*GAMT_U*WOO)
       SPU = SPU+APT*(1.+0.5*GAMT_U)
       SPV=  SPV+APT*(1.+0.5*GAMT_U)
       SPW = SPW+APT*(1.+0.5*GAMT_U)
    END IF
  DEALLOCATE(SEHD_C,SETHD_D,APT)
!==================================================================
!.....FINAL COEFFICIENT AND SOURCES MATRIX FOR U MOMENTUM-EQUATION
!==================================================================

    URFRS=URFR(IU)
    URFMS=URFM(IU)

    AP=-AE-AW-AN-AS-AT-AB+SPU
    AP=AP*URFRS       ! UNDER-RELAXATION
    SU=SU+URFMS*AP*U
    APU=1./(AP+SMALL) !.....1/AP  FOR U = > APU
!===============================
!     SOLVE LINEAR SYSTEM FOR U
!===============================
    CALL SIPSOL(U,IU)

!==================================================================
!.....FINAL COEFFICIENT AND SOURCES MATRIX FOR V MOMENTUM-EQUATION
!==================================================================

    URFRS=URFR(IV)
    URFMS=URFM(IV)

    AP=-AE-AW-AN-AS-AT-AB+SPV
    AP=AP*URFRS         ! UNDER-RELAXATION
    SU=SV+URFMS*AP*V
    APV=1./(AP+SMALL)   !.....1/AP  FOR V = > APV
!===============================
!     SOLVE LINEAR SYSTEM FOR V
!===============================
    CALL SIPSOL(V,IV)

!==================================================================
!.....FINAL COEFFICIENT AND SOURCES MATRIX FOR W MOMENTUM-EQUATION
!==================================================================

    URFRS=URFR(IW)
    URFMS=URFM(IW)

    AP=-AE-AW-AN-AS-AT-AB+SPW
    AP=AP*URFRS     ! UNDER-RELAXATION
    SU=SW+URFMS*AP*W
    APW=1./(AP+SMALL) !.....1/AP  FOR W = > APW
!
!===============================
!     SOLVE LINEAR SYSTEM FOR W
!===============================
    CALL SIPSOL(W,IW)

!=====================================
!     EXCHANGE VARIABLE BETWEEN BLOCKS
!=====================================

!    CALL CHVVEC(U,V,W)
    CALL CHVVEC(APU,APV,APW)
!
  END SUBROUTINE CALUVW

!*******************************************************
  SUBROUTINE CELUVW(NIE,NJE,NKE,IDEW,IDNS,IDTB, &
                    FIF,ACFE,ACFW,FCF,DIR1)
!*******************************************************
    IMPLICIT NONE

    REAL(RK)   ,INTENT(INOUT) :: FIF(NXYZA),FCF(NXYZA)
    REAL(RK)   ,INTENT(INOUT) :: ACFE(NXYZA),ACFW(NXYZA)
    INTEGER(IK),INTENT(IN)    :: NIE,NJE,NKE,IDEW,IDNS,IDTB,DIR1

    REAL(RK)    :: GAM,GAME,DE,FLCF,CE,CW
    REAL(RK)    :: SUEH,SVEH,SWEH,SUEL,SVEL,SWEL
    REAL(RK)    :: UHE,VHE,WHE,ULE,VLE,WLE
    REAL(RK)    :: SUC,SVC,SWC,SUADD,SVADD,SWADD

    INTEGER(IK) :: I,J,K
    REAL(RK)    :: FXE,FXW
    INTEGER(IK) :: INP,INE

    GAM=GDS(IU)
!
!.....CALCULATE EAST,TOP,NORTH  CELL FACE
    DO K=2,NKE
       DO I=2,NIE
          DO J=2,NJE
!
             INP=LK(K)+LI(I)+J
             INE=INP+IDEW
!.....
!.....INTERPOLATION FACTOR
             FXE=FIF(INP) 
             FXW=1.-FXE
!
!.....DIFUSION COEFFICIENT
!
             GAME=VIS(INP)*FXW+VIS(INE)*FXE
!++++++++++++++++++++++++++
             GAME=GAME/REY
!++++++++++++++++++++++++++
             DE=GAME*SQRT(ARE2_T(INP,DIR1)/ARKSI2_T(INP,DIR1))
!
!.....CONVECTION FLUXES - UDS
!
             FLCF=FCF(INP)
             CE=MIN( FLCF,R_0_0)
             CW=MAX( FLCF,R_0_0)
!
             ACFE(INP)=-DE+CE
             ACFW(INE)=-DE-CW
!


!.....        D I F F U S I O N
!.....
             SUEH=GAME*((     2.*GR1X(INP)    *ARX_T(INP,DIR1)+ & 
                         (GR1Y(INP)+GR2X(INP))*ARY_T(INP,DIR1)+ &
                         (GR1Z(INP)+GR3X(INP))*ARZ_T(INP,DIR1))*FXW + & 
                        (     2.*GR1X(INE)    *ARX_T(INP,DIR1)+ & 
                         (GR1Y(INE)+GR2X(INE))*ARY_T(INP,DIR1)+ &
                         (GR1Z(INE)+GR3X(INE))*ARZ_T(INP,DIR1))*FXE)
             SVEH=GAME*(((GR2X(INP)+GR1Y(INP))*ARX_T(INP,DIR1)+ & 
                              2.*GR2Y(INP)    *ARY_T(INP,DIR1)+ &
                         (GR2Z(INP)+GR3Y(INP))*ARZ_T(INP,DIR1))*FXW + & 
                        ((GR2X(INE)+GR1Y(INE))*ARX_T(INP,DIR1)+ & 
                              2.*GR2Y(INE)    *ARY_T(INP,DIR1)+ &
                         (GR2Z(INE)+GR3Y(INE))*ARZ_T(INP,DIR1))*FXE)
             SWEH=GAME*(((GR3X(INP)+GR1Z(INP))*ARX_T(INP,DIR1)+ & 
                         (GR3Y(INP)+GR2Z(INP))*ARY_T(INP,DIR1)+ &
                              2.*GR3Z(INP)    *ARZ_T(INP,DIR1))*FXW + & 
                        ((GR3X(INE)+GR1Z(INE))*ARX_T(INP,DIR1)+ & 
                         (GR3Y(INE)+GR2Z(INE))*ARY_T(INP,DIR1)+ &
                              2.*GR3Z(INE)    *ARZ_T(INP,DIR1))*FXE)
!.....
             SUEL=DE*((GR1X(INP)*FXW+GR1X(INE)*FXE)*AKX_T(INP,DIR1)+ &
                      (GR1Y(INP)*FXW+GR1Y(INE)*FXE)*AKY_T(INP,DIR1)+ &
                      (GR1Z(INP)*FXW+GR1Z(INE)*FXE)*AKZ_T(INP,DIR1))
             SVEL=DE*((GR2X(INP)*FXW+GR2X(INE)*FXE)*AKX_T(INP,DIR1)+ &
                      (GR2Y(INP)*FXW+GR2Y(INE)*FXE)*AKY_T(INP,DIR1)+ &
                      (GR2Z(INP)*FXW+GR2Z(INE)*FXE)*AKZ_T(INP,DIR1))
             SWEL=DE*((GR3X(INP)*FXW+GR3X(INE)*FXE)*AKX_T(INP,DIR1)+ &
                      (GR3Y(INP)*FXW+GR3Y(INE)*FXE)*AKY_T(INP,DIR1)+ &
                      (GR3Z(INP)*FXW+GR3Z(INE)*FXE)*AKZ_T(INP,DIR1))
!.....
!.....     C O N V E C T I O N
!.....
             UHE=FLCF*(U(INE)*FXE+U(INP)*FXW)
             VHE=FLCF*(V(INE)*FXE+V(INP)*FXW)
             WHE=FLCF*(W(INE)*FXE+W(INP)*FXW)
!.....
             ULE=CW*U(INP)+CE*U(INE)
             VLE=CW*V(INP)+CE*V(INE)
             WLE=CW*W(INP)+CE*W(INE)
!..... PART TO ADD
             SUC=-GAM*(UHE-ULE)
             SVC=-GAM*(VHE-VLE)
             SWC=-GAM*(WHE-WLE)
!----------------------------------------

             SUADD=SUC+SUEH-SUEL
             SVADD=SVC+SVEH-SVEL
             SWADD=SWC+SWEH-SWEL
!
             SU(INP)=SU(INP)+SUADD
             SV(INP)=SV(INP)+SVADD
             SW(INP)=SW(INP)+SWADD
             SU(INE)=SU(INE)-SUADD
             SV(INE)=SV(INE)-SVADD
             SW(INE)=SW(INE)-SWADD
!
          END DO
       END DO
    END DO
!
    RETURN
  END SUBROUTINE CELUVW

!*****************************************
  SUBROUTINE CALCSC(IFI,FI,FIO,FIOO) !*
!*****************************************
    IMPLICIT NONE

    REAL(RK),ALLOCATABLE,INTENT(INOUT) :: FI(:)
    REAL(RK)            ,INTENT(INOUT) :: FIO(NXYZA),FIOO(NXYZA)
    INTEGER(IK)         ,INTENT(INOUT) :: IFI

    REAL(RK)                :: URFRS,URFMS
    REAL(RK),ALLOCATABLE    :: APT(:)
    INTEGER(IK)             :: NBL

    ALLOCATE(APT(NXYZA))
!.....
    PRTR=PRTINV(IFI)
!
!====================
!     CLEAR ARRAYS
!====================
    SU=0.
    SP=0.

if(IFI==7) then

    CALL GRADFI(GR1X,GR1Y,GR1Z,FI,IVP)
else
 CALL GRADFI(GR1X,GR1Y,GR1Z,FI)
end if 

!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)

     CALL CELSC(NI-2,NJM,NKM,NJ,1,NIJ,IFI,FI,FX,AE,AW,F1,2)
     CALL CELSC(NIM,NJ-2,NKM,1,NIJ,NJ,IFI,FI,FY,AN,AS,F2,1)
     CALL CELSC(NIM,NJM,NK-2,NIJ,NJ,1,IFI,FI,FZ,AT,AB,F3,3)
!

       CALL MODSC(NBL,IFI,FI)

    END DO     !   END ALL BLOCKS


    IF(LTIME) THEN
       APT=DEN*VOL*DTR
       SU=SU+APT*((1.+GAMT_T)*FIO-0.5*GAMT_T*FIOO)
       SP=SP+APT*(1.+0.5*GAMT_T)
    END IF

    DEALLOCATE(APT)
!==================================================================
!.....FINAL COEFFICIENT AND SOURCES MATRIX FOR SCALAR-EQUATION(S)
!==================================================================

    URFRS=URFR(IFI)
    URFMS=URFM(IFI)

    AP=-AE-AW-AN-AS-AT-AB+SP
    AP=AP*URFRS
    SU=SU+URFMS*AP*FI
!
!===============================
!     SOLVE LINEAR SYSTEM FOR FI
!===============================
    CALL SIPSOL(FI,IFI)

!=====================================
!     EXCHANGE VARIABLE BETWEEN BLOCKS
!=====================================
!    CALL CHVSCA(FI)

  END SUBROUTINE CALCSC

!*******************************************************
  SUBROUTINE CELSC(NIE,NJE,NKE,IDEW,IDNS,IDTB,IFI,FI, &
                   FIF,ACFE,ACFW,FCF,DIR1)
!*******************************************************
    IMPLICIT NONE

    REAL(RK)   ,INTENT(INOUT) :: FI(NXYZA),FIF(NXYZA),FCF(NXYZA)
    REAL(RK)   ,INTENT(INOUT) :: ACFW(NXYZA),ACFE(NXYZA)
    INTEGER(IK),INTENT(IN)    :: NIE,NJE,NKE,IDEW,IDNS,IDTB,IFI,DIR1

    REAL(RK)    :: GAM,GAME,DE,FLCF,CE,CW
    REAL(RK)    :: SUEL,SUEH,UHE,ULE,SUC,SUADD
    REAL(RK)    :: FXE,FXW
    INTEGER(IK) :: INP,INE
    INTEGER(IK) :: I,J,K

    GAM=GDS(IFI)
    PRTR=PRTINV(IFI)
!
!.....CALCULATE EAST,TOP,NORTH  CELL FACE
    DO K=2,NKE
       DO I=2,NIE
          DO J=2,NJE
!
             INP=LK(K)+LI(I)+J
             INE=INP+IDEW
!.....INTERPOLATION FACTOR
             FXE=FIF(INP)
             FXW=1.-FXE
!
             GAME=PRTR*(VIS(INP)*FXW+VIS(INE)*FXE)      !! CHECK HERE, CONSISTENT WITH NON-DIMENSIONAL EQUATIONS;
!                                                !!  JFM - 2010
!.....DIFUSION COEFFICIENT
             DE=GAME*SQRT(ARE2_T(INP,DIR1)/ARKSI2_T(INP,DIR1))
!.....CONVECTION FLUXES - UDS
             FLCF=FCF(INP)
             CE=MIN( FLCF,R_0_0)
             CW=MAX( FLCF,R_0_0)
!
             ACFE(INP)=-DE+CE
             ACFW(INE)=-DE-CW
!

!.....        D I F F U S I O N
!.....
             SUEH=GAME*((GR1X(INP)*FXW+GR1X(INE)*FXE)*ARX_T(INP,DIR1)+ &
                        (GR1Y(INP)*FXW+GR1Y(INE)*FXE)*ARY_T(INP,DIR1)+ &
                        (GR1Z(INP)*FXW+GR1Z(INE)*FXE)*ARZ_T(INP,DIR1))
!.....
             SUEL=DE*((GR1X(INP)*FXW+GR1X(INE)*FXE)*AKX_T(INP,DIR1)+ &
                      (GR1Y(INP)*FXW+GR1Y(INE)*FXE)*AKY_T(INP,DIR1)+ &
                      (GR1Z(INP)*FXW+GR1Z(INE)*FXE)*AKZ_T(INP,DIR1))
!.....
!.....     C O N V E C T I O N
!.....
             UHE=FLCF*(FI(INE)*FXE+FI(INP)*FXW)
!.....
             ULE=CW*FI(INP)+CE*FI(INE)
!..... PART TO ADD
             SUC=-GAM*(UHE-ULE)

             SUADD=SUC+SUEH-SUEL
!
             SU(INP)=SU(INP)+SUADD
             SU(INE)=SU(INE)-SUADD
!
          END DO
       END DO
    END DO
!
  END SUBROUTINE CELSC

!****************************************
  SUBROUTINE NEWFLX !*
!****************************************
    USE GEOM
    IMPLICIT NONE

    REAL(RK)    :: DENE,UE,VE,WE
    INTEGER(IK) :: NBL,I,INP,ISIDE,IND
!
!.....ESTIMATE NEW FLUX
!
    DO NBL=1,NBLOCK

       CALL SETIND(NBL)
!.....EAST, NORTH, TOP, CELL - FACE

    CALL FLX(NI-2,NJM,NKM,NJ,1,NIJ,FX,F1,2)
    CALL FLX(NIM,NJ-2,NKM,1,NIJ,NJ,FY,F2,1)
    CALL FLX(NIM,NJM,NK-2,NIJ,NJ,1,FZ,F3,3)

    DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
      INP=NBL_ST(NBL)+NUM_CFI(I,1)
      IND=NBL_ST(NBL)+NUM_CFI(I,3)
      ISIDE=NUM_CFI(I,2)
!
      CALL INIT_GEOM_BLOCK(ISIDE,I)
!
!.....INTERPOLATION FACTOR
      FXE=FX_CFI(I)
      FXW=1.-FXE
!
!.....DENSITY AT THE CELL FACE
      DENE=DEN(INP)*FXW+DEN(IND)*FXE
!
!.....VELOCITIES
!
      UE=U(INP)*FXW+U(IND)*FXE
      VE=V(INP)*FXW+V(IND)*FXE
      WE=W(INP)*FXW+W(IND)*FXE
!
!.....MASS FLUX
      FLX_CFI(I)=DENE*(UE*ARX+VE*ARY+WE*ARZ)
!
     END DO     !   END BLOCK INTERFACES
!
    END DO     !   END ALL BLOCKS
!
  END SUBROUTINE NEWFLX

!***********************************************
  SUBROUTINE FLX(NIE,NJE,NKE,IDEW,IDNS,IDTB, &
                 FIF,FCF,DIR1)
!***********************************************
    IMPLICIT NONE

    REAL(RK)   ,INTENT(IN)  :: FIF(NXYZA)
    REAL(RK)   ,INTENT(OUT) :: FCF(NXYZA)
    INTEGER(IK),INTENT(IN)  :: IDEW,IDNS,IDTB,NIE,NJE,NKE,DIR1

    REAL(RK)    :: DENE,UE,VE,WE
    REAL(RK)    :: FXE,FXW
    INTEGER(IK) :: INP,INE
    INTEGER(IK) :: I,J,K

!
!.....FLUX AT THE CELL FACE
!
!.....CALCULATE EAST CELL FACE
!
    DO K=2,NKE
       DO I=2,NIE
          DO J=2,NJE
!
             INP=LK(K)+LI(I)+J
             INE=INP+IDEW
!.....INTERPOLATION FACTOR
             FXE=FIF(INP)
             FXW=1.-FXE
!.....
!.....DENSITY AT THE CELL FACE
             DENE=DEN(INP)*FXW+DEN(INE)*FXE
!
!.....VELOCITIES
!
             UE=U(INP)*FXW+U(INE)*FXE
             VE=V(INP)*FXW+V(INE)*FXE
             WE=W(INP)*FXW+W(INE)*FXE
!
!.....MASS FLUX
             FCF(INP)=DENE*(UE*ARX_T(INP,DIR1) + &
                            VE*ARY_T(INP,DIR1) + &
                            WE*ARZ_T(INP,DIR1))
          END DO
       END DO
    END DO
!
  END SUBROUTINE FLX

!*******************************************************
  SUBROUTINE CALCP
!*******************************************************
    IMPLICIT NONE

    REAL(RK)    :: PPREF,VOLINP,SUM,SUMA
    INTEGER(IK) :: INP,NBL,IDEW,IDNS,IDTB,I,J,K

!
!====================
!     CLEAR ARRAYS
!====================
!
    SP=0.
!
!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!
!=============================================
!     CALCULATE FLUXES THROUGH INNER CV-FACES
!=============================================

       CALL CELP(NI-2,NJM,NKM,NJ,1,NIJ,FX,AE,AW,F1,2)
       CALL CELP(NIM,NJ-2,NKM,1,NIJ,NJ,FY,AN,AS,F2,1)
       CALL CELP(NIM,NJM,NK-2,NIJ,NJ,1,FZ,AT,AB,F3,3)
!
!==================================
!.....IMPLEMENT BOUNDARY CONDITIONS
!==================================
!
       CALL MODPR(NBL)
!.....END ALL BLOCKS
    END DO

!====================
!     CLEAR ARRAYS
!====================
    SU=0.

!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!.....THE MAIN LOOP
       IDEW=NJ
       IDNS=1
       IDTB=NIJ

!================================
!     ADD FLUXES FROM BOUNDARIES
!================================
       CALL ADDFLX(NBL)

!==========================================================================
!.....FINAL COEFFICIENT AND SOURCES MATRIX FOR PRESSURE CORRECTION-EQUATION
!==========================================================================
       SUM=0.

       DO K=2,NKM
          DO I=2,NIM
             DO J=2,NJM
                INP=LK(K)+LI(I)+J
                SU(INP)=SU(INP)-F1(INP)+F1(INP-IDEW)-F2(INP)+F2(INP-IDNS) &
                       -F3(INP)+F3(INP-IDTB)
                SUM=SUM+SU(INP)
                AP(INP)=(-AE(INP)-AW(INP)-AN(INP)-AS(INP)-AT(INP)-AB(INP)+SP(INP))
             END DO
          END DO
       END DO
!
!.....END ALL BLOCKS
       END DO
!
!====================
!     CLEAR ARRAYS
!====================
    PP=0.
!
!.....TEST   SUM=0
    IF(LTEST) WRITE(IUCHK,600) SUM
!
!============================================================
!     SOLVE LINEAR SYSTEM FOR PP PRESSURE CORRECTION EQUATION
!============================================================
!
    CALL SIPSOL(PP,IP)
!
!=====================
!     CORRECTOR STAGE
!=====================
!
!....PRESSURE REFERENCE
    PPREF=PP(IJKPR)
!....CALCULATE PRESSURE GRADIENTS
    CALL BPRES(PP)
    CALL GRADFI(GR1X,GR1Y,GR1Z,PP)
!
!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
       DO K=2,NKM
          DO I=2,NIM
             DO J=2,NJM
                INP=LK(K)+LI(I)+J
!
!==============================================
!     CORRECT VELOCITY COMPONNENTS AND PRESSURE
!==============================================
!
                VOLINP=VOL(INP)
                U(INP)=U(INP)-GR1X(INP)*VOLINP*APU(INP)
                V(INP)=V(INP)-GR1Y(INP)*VOLINP*APV(INP)
                W(INP)=W(INP)-GR1Z(INP)*VOLINP*APW(INP)
                P(INP)=P(INP)+URF(IP)*(PP(INP)-PPREF)
             END DO
          END DO
       END DO

!==============================================
!     CORRECT MASS FLUXES(ONLY INNER FLUXES)
!==============================================

       IDEW=NJ
       IDNS=1
       IDTB=NIJ

       DO K=2,NKM
          DO I=2,NIM-1
             DO J=2,NJM
                INP=LK(K)+LI(I)+J
                F1(INP)=F1(INP)+AE(INP)*(PP(INP+IDEW)-PP(INP))
             END DO
          END DO
       END DO

       DO K=2,NKM
          DO I=2,NIM
             DO J=2,NJM-1
                INP=LK(K)+LI(I)+J
                F2(INP)=F2(INP)+AN(INP)*(PP(INP+IDNS)-PP(INP))
             END DO
          END DO
       END DO

       DO K=2,NKM-1
          DO I=2,NIM
             DO J=2,NJM
                INP=LK(K)+LI(I)+J
                F3(INP)=F3(INP)+AT(INP)*(PP(INP+IDTB)-PP(INP))
             END DO
          END DO
       END DO
!====================================================
!.....CORRECT  MASS FLUXES DUE TO BOUNDARY CONDITIONS
!=====================================================
       CALL MODPRC(NBL)

!.....END ALL BLOCKS
    END DO

!=====================
!     TEST STAGE
!=====================

!.....TEST   SUM=0
    IF(LTEST) THEN

!====================
!     CLEAR ARRAYS
!====================
       SU=0.
!
       SUM=0.
       SUMA=0.
!.....ALL BLOCKS
       DO NBL=1,3
          CALL SETIND(NBL)
!
          IDEW=NJ
          IDNS=1
          IDTB=NIJ
!================================
!     ADD FLUXES FROM BOUNDARIES
!================================
          CALL ADDFLX(NBL)

          DO K=2,NKM
             DO I=2,NIM
                DO J=2,NJM
                   INP=LK(K)+LI(I)+J
!.....SOURCE TERMS
                   SU(INP)=SU(INP)-F1(INP)+F1(INP-IDEW)-F2(INP) &
                          +F2(INP-IDNS)-F3(INP)+F3(INP-IDTB)
                   SUM=SUM+SU(INP)
                   SUMA=SUMA+ABS(SU(INP))
                END DO
             END DO
          END DO

!.....END ALL BLOCKS
       END DO
       WRITE(IUCHK,610) SUM,SUMA

    END IF
!
  600 FORMAT(20X,' SUM  =',1PE10.3)
  610 FORMAT(20X,' SUM  =',1PE10.3,/,20X,'|SUM| =',1PE10.3)
  END SUBROUTINE CALCP

!***********************************************************
  SUBROUTINE CELP(NIE,NJE,NKE,IDEW,IDNS,IDTB, &
                  FIF,ACFE,ACFW,FCF,DIR1)
!*******************************************************
    IMPLICIT NONE

    REAL(RK)   ,INTENT(IN)  :: FIF(NXYZA)
    REAL(RK)   ,INTENT(OUT) :: ACFE(NXYZA),ACFW(NXYZA),FCF(NXYZA)
    INTEGER(IK),INTENT(IN)  :: IDEW,IDNS,IDTB,NIE,NJE,NKE,DIR1

    REAL(RK)    :: DUE,DVE,DWE,DENE,GRKSI,DPEKSI,UE,VE,WE
    REAL(RK)    :: FXE,FXW
    INTEGER(IK) :: INP,INE
    INTEGER(IK) :: I,J,K
!
!.....CALCULATE EAST CELL FACE
!
!.....CALCULATE EAST,TOP,NORTH  CELL FACE
    DO K=2,NKE
       DO I=2,NIE
          DO J=2,NJE
!
             INP=LK(K)+LI(I)+J
             INE=INP+IDEW
!.....INTERPOLATION FACTOR
             FXE=FIF(INP)
             FXW=1.-FXE
!
!.....EAST NODE
!      COFHLP=(VOL(INP)*FXW+VOL(INE)*FXE)*SQRT(ARE2*ARKSI2) ! PERIC STYLE
!.....CELL FACE COEFFICIENTS  1/AP(INP)

!
             DUE=(APU(INP)*FXW+APU(INE)*FXE)*ARX_T(INP,DIR1)
             DVE=(APV(INP)*FXW+APV(INE)*FXE)*ARY_T(INP,DIR1)
             DWE=(APW(INP)*FXW+APW(INE)*FXE)*ARZ_T(INP,DIR1)
!.....DENSITY AT THE CELL FACE
             DENE=DEN(INP)*FXW+DEN(INE)*FXE
!
!.....INTERPOLATED  GRAD P  -DOT-  KSI
             GRKSI=0.5*((GR1X(INP)+GR1X(INE))*AKX_T(INP,DIR1)+ &
                        (GR1Y(INP)+GR1Y(INE))*AKY_T(INP,DIR1)+ &
                        (GR1Z(INP)+GR1Z(INE))*AKZ_T(INP,DIR1))
!.....VELOCITIES
!
             DPEKSI=P(INE)-P(INP) - GRKSI
             UE=U(INP)*FXW+U(INE)*FXE-DUE*DPEKSI
             VE=V(INP)*FXW+V(INE)*FXE-DVE*DPEKSI
             WE=W(INP)*FXW+W(INE)*FXE-DWE*DPEKSI
!
!.....COEFFICIENTS OF PRESSURE-CORRECTION EQUATION
             ACFE(INP)=-DENE*(DUE*ARX_T(INP,DIR1) +&
                              DVE*ARY_T(INP,DIR1) +&
                              DWE*ARZ_T(INP,DIR1))
             ACFW(INE)=ACFE(INP)
!
!.....MASS FLUX
             FCF(INP)=DENE*(UE*ARX_T(INP,DIR1) +&
                            VE*ARY_T(INP,DIR1) +&
                            WE*ARZ_T(INP,DIR1))
          END DO
       END DO
    END DO
  END SUBROUTINE CELP


!*******************************************************
  SUBROUTINE ADDFLX(NBL)
!*******************************************************
    IMPLICIT NONE
    INTEGER(IK),INTENT(IN) :: NBL

    INTEGER(IK) :: IDR,NBCTYP,NBC_S,NBC_E,I,INP
!
!.....ADD FLUXES FROM BOUNDARIES
!.....ALL BC
    DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)
!.....TYP
       NBCTYP=NUM_TYP(IDR)
!.....
       NBC_S=NUM_SPR(IDR)+1
       NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)
!
!.....  IF  *INLET* ODER *OUTLET*
       IF(NBCTYP==1.OR.NBCTYP==2)  THEN
          DO I=NBC_S,NBC_E
             INP=NBL_ST(NBL)+NUM_IND(I,1)
             SU(INP)=SU(INP) - FLX_IND(I)
          END DO
       END IF
!..... END ALL BC
    END DO
!
!.....BLOCK INTERFACES
    DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
       INP=NBL_ST(NBL)+NUM_CFI(I,1)

    IF(NUM_CFI(I,7)==0) THEN

       SU(INP)=SU(INP) - FLX_CFI(I)

    END IF
    END DO
  END SUBROUTINE ADDFLX

!*******************************************************
  SUBROUTINE MODUVW(NBL)
!*******************************************************
    IMPLICIT NONE
    INTEGER(IK),INTENT(IN) :: NBL

    INTEGER(IK) :: IDR,NBCTYP,NBC_S,NBC_E,I,INP
!***************************************************
!.....U V W - M O M E N T U M  BOUNDARY CONDITIONS *
!***************************************************
!                                                  |
!__________________________________________________|
!
!.....ALL BC
    DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)
!.....TYP
       NBCTYP=NUM_TYP(IDR)
!.....
       NBC_S=NUM_SPR(IDR)+1
       NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)
       IF(NBCTYP==1)  THEN
!********************************
!*****     I N L E T  B.C. ******
!********************************
          DO I=NBC_S,NBC_E
             INP=NBL_ST(NBL)+NUM_IND(I,1)
             IF(FLX_IND(I)<R_0_0) THEN
                CALL INLUVW(INP,NUM_IND(I,2),FLX_IND(I))  ! IF INLET
             ELSE
                CALL OUTUVW(INP,NUM_IND(I,2),FLX_IND(I))  ! IF OUTLET
             END IF
!
             SU(INP)=SU(INP)+SUU
             SV(INP)=SV(INP)+SUV
             SW(INP)=SW(INP)+SUW

             SPU(INP) =   SPU(INP)+SUP
             SPV(INP)=    SPV(INP)+SVP
             SPW(INP) =   SPW(INP)+SWP
                
          END DO
!
       ELSE IF(NBCTYP==2) THEN
!********************************
!*****     O U T L E T   B.C. ***
!********************************

          DO I=NBC_S,NBC_E
             INP=NBL_ST(NBL)+NUM_IND(I,1)
             IF(FLX_IND(I)<R_0_0) THEN
                CALL INLUVW(INP,NUM_IND(I,2),FLX_IND(I))  ! IF INLET
             ELSE
                CALL OUTUVW(INP,NUM_IND(I,2),FLX_IND(I))  ! IF OUTLET
             END IF
!
             SU(INP)=SU(INP)+SUU
             SV(INP)=SV(INP)+SUV
             SW(INP)=SW(INP)+SUW

             SPU(INP) =   SPU(INP)+SUP
             SPV(INP)=    SPV(INP)+SVP
             SPW(INP) =   SPW(INP)+SWP

          END DO

!********************************
       ELSE IF(NBCTYP==3) THEN
!********************************
!*****  S Y M M E T R Y   B.C. **
!********************************
          DO I=NBC_S,NBC_E
             INP=NBL_ST(NBL)+NUM_IND(I,1)
             CALL SYMUVW(INP,NUM_IND(I,2))
!
             SU(INP)=SU(INP)+SUU
             SV(INP)=SV(INP)+SUV
             SW(INP)=SW(INP)+SUW

             SPU(INP) =   SPU(INP)+SUP
             SPV(INP)=    SPV(INP)+SVP
             SPW(INP) =   SPW(INP)+SWP
          END DO

!      ELSE IF(NBCTYP==4) THEN
       ELSE                          !!! HERE WE USE 'ELSE' TO CONSIDER TYPE _ 10, 12, 14,4 
!********************************
!*****     W A L L        B.C. **
!********************************
          DO I=NBC_S,NBC_E
             INP=NBL_ST(NBL)+NUM_IND(I,1)
             CALL WALUVW(INP,NUM_IND(I,2))
!
             SU(INP)=SU(INP)+SUU
             SV(INP)=SV(INP)+SUV
             SW(INP)=SW(INP)+SUW

             SPU(INP) =   SPU(INP)+SUP
             SPV(INP)=    SPV(INP)+SVP
             SPW(INP) =   SPW(INP)+SWP
!
          END DO
       END IF
!
!.....END ALL BC
    END DO

!=================================
!  **   BLOCK INTERFACES  **
!=================================

    CALL COMUVW(NBL)
!

   END SUBROUTINE MODUVW

!******************************************
  SUBROUTINE INLUVW(INP,ISIDE,FLCF)  !*
!******************************************
    USE GEOM
    IMPLICIT NONE

    REAL(RK)   ,INTENT(INOUT) :: FLCF
    INTEGER(IK),INTENT(INOUT) :: INP,ISIDE

    REAL(RK)    :: GAME,DE,CE,ADC
    REAL(RK)    :: SUEH,SVEH,SWEH
    REAL(RK)    :: SUEL,SVEL,SWEL

    CALL INIT_GEOM_BF(ISIDE,INP)
!
    GAME=VIS(INE)
!++++++++++++++++++++++++++++
    GAME=GAME/REY
!++++++++++++++++++++++++++++
!.....DIFUSION COEFFICIENT
    DE=GAME*SQRT(ARE2/ARKSI2)
!
!.....EXPLICIT PART OF DIFFUSION FLUXES (HIGH ODER)
!.....GRAD UI +(GRAD V)*I   ;I=1,2,3
    SUEH=GAME*(GR1X(INP)*ARX+GR1Y(INP)*ARY+GR1Z(INP)*ARZ+ &
               GR1X(INP)*ARX+GR2X(INP)*ARY+GR3X(INP)*ARZ)
    SVEH=GAME*(GR2X(INP)*ARX+GR2Y(INP)*ARY+GR2Z(INP)*ARZ+ &
               GR1Y(INP)*ARX+GR2Y(INP)*ARY+GR3Y(INP)*ARZ)
    SWEH=GAME*(GR3X(INP)*ARX+GR3Y(INP)*ARY+GR3Z(INP)*ARZ+ &
               GR1Z(INP)*ARX+GR2Z(INP)*ARY+GR3Z(INP)*ARZ)
!
!.....
    SUEL=DE*(GR1X(INP)*AKX+GR1Y(INP)*AKY+GR1Z(INP)*AKZ)
    SVEL=DE*(GR2X(INP)*AKX+GR2Y(INP)*AKY+GR2Z(INP)*AKZ)
    SWEL=DE*(GR3X(INP)*AKX+GR3Y(INP)*AKY+GR3Z(INP)*AKZ)
!
    CE=MIN(FLCF,R_0_0)
    ADC=-DE+CE
!
!.....EXPLICIT PART OF DIFFUSION FLUXES
!
    SUU=SUEH-SUEL-ADC*U(INE)
    SUV=SVEH-SVEL-ADC*V(INE)
    SUW=SWEH-SWEL-ADC*W(INE)
    SUP=-ADC
    SVP=-ADC
    SWP=-ADC
!
  END SUBROUTINE INLUVW


!*************************************
  SUBROUTINE OUTUVW(INP,ISIDE,FLCF)
!*************************************
!
!--------------------------------------------------------------
!      CONSTANT GRADIENT BETWEEN BOUNDARY & CV-CENTER ASSUMED
!---------------------------------------------------------------
!
    USE GEOM
    IMPLICIT NONE

    REAL(RK)   ,INTENT(INOUT) :: FLCF
    INTEGER(IK),INTENT(INOUT) :: ISIDE,INP

    REAL(RK)    :: GAME,DE,CE,CW
    REAL(RK)    :: SUEH,SVEH,SWEH
    REAL(RK)    :: SUEL,SVEL,SWEL

    CALL INIT_GEOM_BF(ISIDE,INP)
!
    GAME=VIS(INP)
!++++++++++++++++++++++++++++
    GAME=GAME/REY
!++++++++++++++++++++++++++++
!.....DIFUSION COEFFICIENT
    DE=GAME*SQRT(ARE2/ARKSI2)
!

    SUEH=GAME*(GR1X(INP)*ARX+GR1Y(INP)*ARY+GR1Z(INP)*ARZ+ &
               GR1X(INP)*ARX+GR2X(INP)*ARY+GR3X(INP)*ARZ)
    SVEH=GAME*(GR2X(INP)*ARX+GR2Y(INP)*ARY+GR2Z(INP)*ARZ+ &
               GR1Y(INP)*ARX+GR2Y(INP)*ARY+GR3Y(INP)*ARZ)
    SWEH=GAME*(GR3X(INP)*ARX+GR3Y(INP)*ARY+GR3Z(INP)*ARZ+ &
               GR1Z(INP)*ARX+GR2Z(INP)*ARY+GR3Z(INP)*ARZ)
!
!.....
    SUEL=DE*(GR1X(INP)*AKX+GR1Y(INP)*AKY+GR1Z(INP)*AKZ)
    SVEL=DE*(GR2X(INP)*AKX+GR2Y(INP)*AKY+GR2Z(INP)*AKZ)
    SWEL=DE*(GR3X(INP)*AKX+GR3Y(INP)*AKY+GR3Z(INP)*AKZ)
!
!.....EXPLICIT PART OF DIFFUSION FLUXES
!
    SUU=SUEH-SUEL
    SUV=SVEH-SVEL
    SUW=SWEH-SWEL
    SUP=0.
    SVP=0.
    SWP=0.
!
  END SUBROUTINE OUTUVW


!*************************************
  SUBROUTINE SYMUVW(INP,ISIDE)  !*
!*************************************
    USE GEOM
    IMPLICIT NONE

    INTEGER(IK),INTENT(INOUT) :: ISIDE,INP

    REAL(RK)    :: TCOEF,TAR,UVWN

    CALL  INIT_GEOM_BF(ISIDE,INP)
    ARX=ARX/(ARE)
    ARY=ARY/(ARE)
    ARZ=ARZ/(ARE)
!
!.....BOUNDARY CONDITION FOR CARTESIAN COMPONENTS
!              OF VELOCITIES
!.....VISCOSITY FROM INNER CV-POINT
    TCOEF=VIS(INP)*DELNR
!++++++++++++++++++++++++++++
    TCOEF=TCOEF/REY
!++++++++++++++++++++++++++++
    TAR=2.*TCOEF*ARE
!.....VELOCITY NORMAL TO THE BOUNDARY FACE
!
    UVWN=U(INP)*ARX+V(INP)*ARY+W(INP)*ARZ
!.....CARTESIAN COMPONENTS OF THE  PARALEL VELOCITY
!.....ON THE BOUNDARY
    U(INE)=U(INP)-UVWN*ARX
    V(INE)=V(INP)-UVWN*ARY
    W(INE)=W(INP)-UVWN*ARZ
!
!.....SOURCE TERMS
    SUU=-TAR*ARX*(ARY*V(INP)+ARZ*W(INP))
    SUV=-TAR*ARY*(ARX*U(INP)+ARZ*W(INP))
    SUW=-TAR*ARZ*(ARX*U(INP)+ARY*V(INP))
    SUP=TAR*ARX**2
    SVP=TAR*ARY**2
    SWP=TAR*ARZ**2
!
  END SUBROUTINE SYMUVW

!************************************
  SUBROUTINE WALUVW(INP,ISIDE)
!************************************
    USE GEOM
    IMPLICIT NONE

    INTEGER(IK),INTENT(INOUT) :: ISIDE,INP

    REAL(RK)    :: TCOEF,TAR

    CALL  INIT_GEOM_BF(ISIDE,INP)
    ARX=ARX/(ARE)
    ARY=ARY/(ARE)
    ARZ=ARZ/(ARE)
!
!.....BOUNDARY CONDITION FOR CARTESIAN COMPONENTS
!              OF VELOCITIES

    DELN=ABS(DELN)      ! ADDED TO OVERCOME  ROUND OFF
!.....VISCOSITY FROM INNER CV-POINT
    TCOEF=VIS(INP)*DELNR
!++++++++++++++++++++++++++++
    TCOEF=TCOEF/REY
!++++++++++++++++++++++++++++
!
    TAR=TCOEF*ARE
!.....SOURCE TERMS
    SUU=TAR*(ARX*(ARY*(V(INP)-V(INE))+ARZ*(W(INP)-W(INE))) &
         +(1.-ARX**2)*U(INE))
    SUV=TAR*(ARY*(ARX*(U(INP)-U(INE))+ARZ*(W(INP)-W(INE))) &
         +(1.-ARY**2)*V(INE))
    SUW=TAR*(ARZ*(ARX*(U(INP)-U(INE))+ARY*(V(INP)-V(INE))) &
         +(1.-ARZ**2)*W(INE))
    SUP=TAR*(1.-ARX**2)
    SVP=TAR*(1.-ARY**2)
    SWP=TAR*(1.-ARZ**2)
!
  END SUBROUTINE WALUVW

!*********************************
  SUBROUTINE COMUVW(NBL)
!*********************************
    USE GEOM
    IMPLICIT NONE

    INTEGER(IK),INTENT(IN) :: NBL
    REAL(RK)    :: GAM,GAME,DE,FLCF,CE,CW,SVEH
    REAL(RK)    :: SUEH,SWEH,SUEL,SVEL,SWEL
    REAL(RK)    :: UHE,VHE,WHE,ULE,VLE,WLE
    REAL(RK)    :: SUC,SVC,SWC,SUADD,SVADD,SWADD
    INTEGER(IK) :: I,INP,ISIDE,IND
!
    GAM=GDS(IU)

    DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
       INP=NBL_ST(NBL) + NUM_CFI(I,1)
       IND=NBL_ST(NBL) + NUM_CFI(I,3)
       ISIDE=NUM_CFI(I,2)
!
  IF(NUM_CFI(I,7)==0) THEN 
       CALL  INIT_GEOM_BLOCK(ISIDE,I)

      AKX=XC(IND)-XC(INP)
      AKY=YC(IND)-YC(INP)
      AKZ=ZC(IND)-ZC(INP)
      ARKSI2=AKX**2+AKY**2+AKZ**2
!
!.....INTERPOLATION FACTOR
       FXE=FX_CFI(I)
       FXW=1.-FXE
!
       GAME=(VIS(INP)*FXW+VIS(IND)*FXE)
!++++++++++++++++++++++++++++
       GAME=GAME/REY
!++++++++++++++++++++++++++++
!
!.....DIFUSION COEFFICIENT
       DE=GAME*SQRT(ARE2/ARKSI2)
!.....CONVECTION FLUXES - UDS
       FLCF=FLX_CFI(I)
       CE=MIN( FLCF,R_0_0)
       CW=MAX( FLCF,R_0_0)
!
       A_E(I)=-DE+CE
!
!.....
!.....        D I F F U S I O N
!.....
       SUEH=GAME*((GR1X(INP)*FXW+GR1X(IND)*FXE)*ARX+ &
                  (GR1Y(INP)*FXW+GR1Y(IND)*FXE)*ARY+ &
                  (GR1Z(INP)*FXW+GR1Z(IND)*FXE)*ARZ+ &
                  (GR1X(INP)*FXW+GR1X(IND)*FXE)*ARX+ &
                  (GR2X(INP)*FXW+GR2X(IND)*FXE)*ARY+ &
                  (GR3X(INP)*FXW+GR3X(IND)*FXE)*ARZ)
       SVEH=GAME*((GR2X(INP)*FXW+GR2X(IND)*FXE)*ARX+ &
                  (GR2Y(INP)*FXW+GR2Y(IND)*FXE)*ARY+ &
                  (GR2Z(INP)*FXW+GR2Z(IND)*FXE)*ARZ+ &
                  (GR1Y(INP)*FXW+GR1Y(IND)*FXE)*ARX+ &
                  (GR2Y(INP)*FXW+GR2Y(IND)*FXE)*ARY+ &
                  (GR3Y(INP)*FXW+GR3Y(IND)*FXE)*ARZ)
       SWEH=GAME*((GR3X(INP)*FXW+GR3X(IND)*FXE)*ARX+ &
                  (GR3Y(INP)*FXW+GR3Y(IND)*FXE)*ARY+ &
                  (GR3Z(INP)*FXW+GR3Z(IND)*FXE)*ARZ+ &
                  (GR1Z(INP)*FXW+GR1Z(IND)*FXE)*ARX+ &
                  (GR2Z(INP)*FXW+GR2Z(IND)*FXE)*ARY+ &
                  (GR3Z(INP)*FXW+GR3Z(IND)*FXE)*ARZ)
!.....
       SUEL=DE*((GR1X(INP)*FXW+GR1X(IND)*FXE)*AKX+ &
                (GR1Y(INP)*FXW+GR1Y(IND)*FXE)*AKY+ &
                (GR1Z(INP)*FXW+GR1Z(IND)*FXE)*AKZ)
       SVEL=DE*((GR2X(INP)*FXW+GR2X(IND)*FXE)*AKX+ &
                (GR2Y(INP)*FXW+GR2Y(IND)*FXE)*AKY+ &
                (GR2Z(INP)*FXW+GR2Z(IND)*FXE)*AKZ)
       SWEL=DE*((GR3X(INP)*FXW+GR3X(IND)*FXE)*AKX+ &
                (GR3Y(INP)*FXW+GR3Y(IND)*FXE)*AKY+ &
                (GR3Z(INP)*FXW+GR3Z(IND)*FXE)*AKZ)
!.....
!.....     C O N V E C T I O N
!.....
       UHE=FLCF*(U(IND)*FXE+U(INP)*FXW)
       VHE=FLCF*(V(IND)*FXE+V(INP)*FXW)
       WHE=FLCF*(W(IND)*FXE+W(INP)*FXW)
!.....
       ULE=CW*U(INP)+CE*U(IND)
       VLE=CW*V(INP)+CE*V(IND)
       WLE=CW*W(INP)+CE*W(IND)
!..... PART TO ADD
       SUC=-GAM*(UHE-ULE)
       SVC=-GAM*(VHE-VLE)
       SWC=-GAM*(WHE-WLE)
!---------------------------------------
       SUADD=SUC+SUEH-SUEL
       SVADD=SVC+SVEH-SVEL
       SWADD=SWC+SWEH-SWEL
!
       SU(INP)=SU(INP)+SUADD
       SV(INP)=SV(INP)+SVADD
       SW(INP)=SW(INP)+SWADD

       SPU(INP)=  SPU(INP)-A_E(I)
       SPV(INP)=  SPV(INP)-A_E(I)
       SPW(INP)=  SPW(INP)-A_E(I)
!
 END IF
    END DO

  END SUBROUTINE COMUVW

!*********************************************************************
  SUBROUTINE OUTBC
!*********************************************************************
    USE GEOM
    IMPLICIT NONE

    INTEGER(IK) :: NBL,IDR,NBC_S,NBC_E,I,INP,ISIDE,IDEW
    LOGICAL     :: LONEG
!
!..... EXTRA CHECK IF NEGATIV FLUX OCCURS AT OUTLET
!
!.....CORRECT OUTFLOW TO FULFILL OVERALL MASS BALANCE
!
    LONEG=.FALSE.
    FLOW=0.0
!
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!
!.....CHECK - TAKE ONLY OUTLET REGION
       DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)
!
!..... T Y P  OUTLET
          IF(NUM_TYP(IDR)==2) THEN
             NBC_S=NUM_SPR(IDR)+1
             NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)
!
             DO I=NBC_S,NBC_E
                INP=NBL_ST(NBL)+NUM_IND(I,1)
                ISIDE=NUM_IND(I,2)
!

                CALL INIT_GEOM_BF(ISIDE,INP)
!
!.....VELOCITY FROM CENTRAL POINT
                FLX_IND(I)=DEN(INE)*(U(INP)*ARX+V(INP)*ARY+W(INP)*ARZ)
!
                FLOW=FLOW+FLX_IND(I)
!..... CHECK
                LONEG=LONEG.OR.FLX_IND(I)<R_0_0
!
             END DO
          END IF
!
!.....END REGIONS
       END DO
!
!.....END BLOCKS
    END DO
!
!.....CALCULATE CORRECTION FACTOR
    FACOUT=FLOWIN/(FLOW+SMALL)
    IF(LTEST) WRITE(IUCHK,600) FACOUT,FLOWIN,FLOW
    IF(LONEG.AND.LTEST) WRITE(IUCHK,610)
!
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!
!.....CHECK - TAKE ONLY OUTLET REGION
       DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)
!
!..... T Y P  OUTLET
          IF(NUM_TYP(IDR)==2) THEN
             NBC_S=NUM_SPR(IDR)+1
             NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)
!
             DO I=NBC_S,NBC_E
                INP=NBL_ST(NBL)+NUM_IND(I,1)
                ISIDE=NUM_IND(I,2)
!
                IDEW=ID_SIDE  (ID_ODER(ISIDE,1))*ID_SIGN(ISIDE)
!
!.....CALCULATE INDEX
                INE=INP+IDEW
!
!.....VELOCITY FROM CENTRAL POINT
                FLX_IND(I)=FLX_IND(I)*FACOUT

!to manage neumann condition for the walls in SH  model (BC 4 == BC 2) 
	!				   		FACOUT=1.
                U(INE)=U(INP)*FACOUT
                V(INE)=V(INP)*FACOUT
                W(INE)=W(INP)*FACOUT
             END DO
          END IF
!
!.....END REGIONS
       END DO
!
!.....END BLOCKS
    END DO
!
  600 FORMAT(' OUTBC- ',' FACTOR= ',1PE13.6,' FLOWIN= ',1PE13.6, &
             ' FLOW= ',1PE13.6)
  610 FORMAT('   =>  NEGATIV FLUX AT OUTLET  <=')
  END SUBROUTINE OUTBC

!*****************************
  SUBROUTINE MODPR(NBL)
!*****************************
    IMPLICIT NONE
    INTEGER(IK),INTENT(IN) :: NBL
!********************************************
!..... P R E S S U R E  BOUNDARY CONDITIONS *
!********************************************
!
!.....BLOCK INTERFACES
    CALL COMPR(NBL)
!
  END SUBROUTINE MODPR

!*********************************
  SUBROUTINE COMPR(NBL)
!*********************************
    USE GEOM
    IMPLICIT NONE
    INTEGER(IK),INTENT(IN) :: NBL

    REAL(RK)    :: DENE,GRKSI,DPEKSI
    REAL(RK)    :: DUE,DVE,DWE
    REAL(RK)    :: UE,VE,WE
    INTEGER(IK) :: I,INP,ISIDE,IND

!
    DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
       INP=NBL_ST(NBL) + NUM_CFI(I,1)
       IND=NBL_ST(NBL) + NUM_CFI(I,3)
       ISIDE=NUM_CFI(I,2)
!
    IF(NUM_CFI(I,7)==0) THEN
       CALL INIT_GEOM_BLOCK(ISIDE,I)

      AKX=XC(IND)-XC(INP)
      AKY=YC(IND)-YC(INP)
      AKZ=ZC(IND)-ZC(INP)
      ARKSI2=AKX**2+AKY**2+AKZ**2

!.....INTERPOLATION FACTOR
       FXE=FX_CFI(I)
       FXW=1.-FXE
!
       DUE=(APU(INP)*FXW+APU(IND)*FXE)*ARX
       DVE=(APV(INP)*FXW+APV(IND)*FXE)*ARY
       DWE=(APW(INP)*FXW+APW(IND)*FXE)*ARZ
!.....DENSITY AT THE CELL FACE
       DENE=DEN(INP)*FXW+DEN(IND)*FXE
!
!.....INTERPOLATED  GRAD P  -DOT-  KSI
       GRKSI=(GR1X(INP)*FXW+GR1X(IND)*FXE)*AKX+ &
             (GR1Y(INP)*FXW+GR1Y(IND)*FXE)*AKY+ &
             (GR1Z(INP)*FXW+GR1Z(IND)*FXE)*AKZ
!.....VELOCITIES
!
       DPEKSI=P(IND)-P(INP) - GRKSI
       UE=U(INP)*FXW+U(IND)*FXE-DUE*DPEKSI
       VE=V(INP)*FXW+V(IND)*FXE-DVE*DPEKSI
       WE=W(INP)*FXW+W(IND)*FXE-DWE*DPEKSI
!
!.....COEFFICIENTS OF PRESSURE-CORRECTION EQUATION
       A_E(I)=-DENE*(DUE*ARX+DVE*ARY+DWE*ARZ)
!
!.....MASS FLUX
       FLX_CFI(I)=DENE*(UE*ARX+VE*ARY+WE*ARZ)
       SP(INP)=SP(INP)- A_E(I)
!
    END IF
    END DO
  END SUBROUTINE COMPR

!*****************************
  SUBROUTINE MODPRC(NBL)
!*****************************
    IMPLICIT NONE
    INTEGER(IK),INTENT(IN) :: NBL
!********************************************
!..... P R E S S U R E  BOUNDARY CONDITIONS *
!********************************************
!
!.....BLOCK INTERFACES
    CALL COMPRC(NBL)
!
  END SUBROUTINE MODPRC

!*********************************
  SUBROUTINE COMPRC(NBL)
!*********************************
    IMPLICIT NONE
    INTEGER(IK),INTENT(IN) :: NBL
    INTEGER(IK) :: I,INP,IND
!
    DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
       INP=NBL_ST(NBL) + NUM_CFI(I,1)
       IND=NBL_ST(NBL) + NUM_CFI(I,3)

    IF(NUM_CFI(I,7)==0) THEN
!..... CORRECT MASS FLUX
       FLX_CFI(I)=FLX_CFI(I)+A_E(I)*(PP(IND)-PP(INP))
     END IF
    END DO
  END SUBROUTINE COMPRC

!************************************
  SUBROUTINE MODSC(NBL,IFI,FI)
!************************************
    USE GEOM
    IMPLICIT NONE

    REAL(RK)   ,INTENT(INOUT) :: FI(NXYZA)
    INTEGER(IK),INTENT(IN)    :: NBL
    INTEGER(IK),INTENT(INOUT) :: IFI

    REAL(RK)    :: DE
    INTEGER(IK) :: IDR,NBCTYP
    INTEGER(IK) :: NBC_S,NBC_E,I,INP,ISIDE
!***************************************************
!.....     BOUNDARY CONDITIONS  FOR SCALAR(S)      *
!***************************************************
!           T Y P                                  |
!             1         I N L E T                  |
!             2         O U T L E T                |
!             3         S Y M M E T R Y            |
!             4         W A L L                    |
!       7         W A L L:ADIABTIC                 |
!       8         W A L L:ISOTHERMO, HOT           |
!       9         W A L L:ISOTHERMO, LOW           |
!__________________________________________________|
!
!.....ALL BC
    DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)
!.....TYP
       NBCTYP=NUM_TYP(IDR)
!.....
       NBC_S=NUM_SPR(IDR)+1
       NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)
!
       IF(NBCTYP==1)  THEN
!********************************
!*****     I N L E T  B.C. ******
!********************************
          DO I=NBC_S,NBC_E
             INP=NBL_ST(NBL)+NUM_IND(I,1)
             CALL INLSC(INP,NUM_IND(I,2),FLX_IND(I),IFI,FI)
!
!.....FOR FI  SP => SP
             SU(INP)=SU(INP)+SUU
             SP(INP)=SP(INP)+SUP
          END DO
!
       ELSE IF(NBCTYP==2) THEN
!********************************
!*****     O U T L E T   B.C. ***
!********************************

          DO I=NBC_S,NBC_E
             INP=NBL_ST(NBL)+NUM_IND(I,1)
             CALL OUTSC(INP,NUM_IND(I,2),FLX_IND(I),IFI,FI)
!
             SU(INP)=SU(INP)+SUU
             SP(INP)=SP(INP)+SUP
          END DO
!
!********************************
       ELSE IF(NBCTYP==3) THEN
!********************************
!*****  S Y M M E T R Y   B.C. **
!********************************
          DO I=NBC_S,NBC_E
             INP=NBL_ST(NBL)+NUM_IND(I,1)
             CALL SYMSC(INP,NUM_IND(I,2),IFI,FI)

             SU(INP)=SU(INP)+SUU
             SP(INP)=SP(INP)+SUP

          END DO
!
       ELSE IF(NBCTYP==7) THEN
!
!......ADIABTIC
!
          DO I = NBC_S,NBC_E
!
             INP = NBL_ST(NBL) + NUM_IND(I,1)
             ISIDE = NUM_IND(I,2)

             CALL INIT_GEOM_BF(ISIDE,INP)
!
             FI(INE) = FI(INP)
!
          END DO

       ELSE IF((NBCTYP==10).OR.(NBCTYP==12).OR.(NBCTYP==14).OR.(NBCTYP==4)) THEN
!
          DO I = NBC_S,NBC_E
             INP = NBL_ST(NBL) + NUM_IND(I,1)
             ISIDE = NUM_IND(I,2)
             CALL INIT_GEOM_BF(ISIDE,INP)
!
!......ISOTHERMOL
!
             DE=VIS(INP)*PRTINV(IFI)*SQRT(ARE2/ARKSI2)

             SUP = DE
             SUU = DE*FI(INE)

             SP(INP)=SP(INP)+SUP
             SU(INP)=SU(INP)+SUU

          END DO
       END IF
!.....END ALL BC
    END DO
!
!.....BLOCK INTERFACES
    CALL COMSC(NBL,IFI,FI)
  END SUBROUTINE MODSC

!************************************
  SUBROUTINE MODSCD(NBL,INP)
!************************************
    IMPLICIT NONE
    INTEGER(IK),INTENT(INOUT) :: NBL,INP
    INTEGER(IK)               :: I,INPTES

!..... IT SETS TO ZERO ALL COEFFICIENTS TO THE NEIGHBOUR BLCOKS
    DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
       INPTES=NBL_ST(NBL)+NUM_CFI(I,1)
       IF(INP==INPTES) THEN
          A_E(I)=0.0
       END IF
    END DO
!
  END SUBROUTINE MODSCD

!********************************************
  SUBROUTINE INLSC(INP,ISIDE,FLCF,IFI,FI)
!*********************************************
    USE GEOM
    IMPLICIT NONE
    REAL(RK)   ,INTENT(INOUT) :: FI(NXYZA),FLCF
    INTEGER(IK),INTENT(INOUT) :: INP,ISIDE,IFI
!
    REAL(RK)    :: GAME,DE,SUEH,SUEL,CE,ADC
!
    CALL INIT_GEOM_BF(ISIDE,INP)
!
    PRTR=PRTINV(IFI)
    GAME=PRTR*VIS(INE)
!.....DIFUSION COEFFICIENT
    DE=GAME*SQRT(ARE2/ARKSI2)
!
!.....EXPLICIT PART OF DIFFUSION FLUXES (HIGH ODER)
!.....GRAD FI
    SUEH=GAME*(GR1X(INP)*ARX+GR1Y(INP)*ARY+GR1Z(INP)*ARZ+ &
               GR1X(INP)*ARX+GR2X(INP)*ARY+GR3X(INP)*ARZ)
!
!.....
    SUEL=DE*(GR1X(INP)*AKX+GR1Y(INP)*AKY+GR1Z(INP)*AKZ)
!
    CE=MIN(FLCF,R_0_0)
    ADC=-DE+CE
!
!.....EXPLICIT PART OF DIFFUSION FLUXES
!
    SUU=SUEH-SUEL-ADC*FI(INE)
    SUP=-ADC
!
  END SUBROUTINE INLSC


!*********************************************
  SUBROUTINE OUTSC(INP,ISIDE,FLCF,IFI,FI)
!*********************************************

!---------------------------------------------------------------
!      CONSTANT GRADIENT BETWEEN BOUNDARY & CV-CENTER ASSUMED
!---------------------------------------------------------------
    USE GEOM
    IMPLICIT NONE
!
    REAL(RK)   ,INTENT(INOUT) :: FI(NXYZA),FLCF
    INTEGER(IK),INTENT(INOUT) :: INP,ISIDE,IFI

    REAL(RK)    :: GAM,GAME
    REAL(RK)    :: CE,CW,DE,SUEH,SUEL
    REAL(RK)    :: UHE,ULE,SUC
!
    CALL INIT_GEOM_BF(ISIDE,INP)
!
    PRTR=PRTINV(IFI)
    GAME=PRTR*VIS(INP)
!.....DIFUSION COEFFICIENT
    DE=GAME*SQRT(ARE2/ARKSI2)
!
!.....EXPLICIT PART OF DIFFUSION FLUXES (HIGH ODER)
!.....GRAD FI
    SUEH=GAME*(GR1X(INP)*ARX+GR1Y(INP)*ARY+GR1Z(INP)*ARZ+ &
               GR1X(INP)*ARX+GR2X(INP)*ARY+GR3X(INP)*ARZ)
!
!.....
    SUEL=DE*(GR1X(INP)*AKX+GR1Y(INP)*AKY+GR1Z(INP)*AKZ)
!
    CE=MIN( FLCF,R_0_0)
    CW=MAX( FLCF,R_0_0)
!.....
!.....     C O N V E C T I O N
!.....
    UHE=FLCF*FI(INE)
!.....
    ULE=CW*FI(INP)+CE*FI(INE)
!..... PART TO ADD
    GAM=GDS(IFI)
    SUC=-GAM*(UHE-ULE)
!
!.....EXPLICIT PART OF DIFFUSION FLUXES
!
    SUU=SUEH-SUEL+SUC
    SUP=0.
!
  END SUBROUTINE OUTSC


!***********************************
  SUBROUTINE SYMSC(INP,ISIDE,IFI,FI)
!***********************************
    USE GEOM
    IMPLICIT NONE
!
    REAL(RK)   ,INTENT(INOUT) :: FI(NXYZA)
    INTEGER(IK),INTENT(INOUT) :: INP,ISIDE,IFI
!
    CALL  INIT_GEOM_BF(ISIDE,INP)

!
!.....ON THE SYMMETRY
    FI(INE)=FI(INP)
!.....SOURCE TERMS = > SUU=SUP  = 0.
!
    SUU=0.
    SUP=0.
  END SUBROUTINE SYMSC


!************************************
  SUBROUTINE WALSC(INP,ISIDE)
!************************************
    IMPLICIT NONE
    INTEGER(IK),INTENT(INOUT) :: INP,ISIDE

    SUU=0.
    SUP=0.
  END SUBROUTINE WALSC


!*********************************
  SUBROUTINE COMSC(NBL,IFI,FI)
!*********************************
    USE GEOM
    IMPLICIT NONE
!
    REAL(RK)   ,INTENT(INOUT) :: FI(NXYZA)
    INTEGER(IK),INTENT(IN)    :: NBL
    INTEGER(IK),INTENT(INOUT) :: IFI

    REAL(RK)    :: GAM,GAME,DE,CE,FLCF,CW
    REAL(RK)    :: SUEH,SUEL,UHE,ULE,SUC,SUADD
    INTEGER(IK) :: I,INP,ISIDE,IND
!.....
    GAM=GDS(IFI)
    PRTR=PRTINV(IFI)

!
    DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
       INP=NBL_ST(NBL) + NUM_CFI(I,1)
       IND=NBL_ST(NBL) + NUM_CFI(I,3)
       ISIDE=NUM_CFI(I,2)
!
   IF(NUM_CFI(I,7)==0) THEN
       CALL INIT_GEOM_BLOCK(ISIDE,I)

      AKX=XC(IND)-XC(INP)
      AKY=YC(IND)-YC(INP)
      AKZ=ZC(IND)-ZC(INP)
      ARKSI2=AKX**2+AKY**2+AKZ**2
!
!.....INTERPOLATION FACTOR
       FXE=FX_CFI(I)
       FXW=1.-FXE
!
       GAME=PRTR*(VIS(INP)*FXW+VIS(IND)*FXE)
!
!.....DIFUSION COEFFICIENT
       DE=GAME*SQRT(ARE2/ARKSI2)
!.....CONVECTION FLUXES - UDS
       FLCF=FLX_CFI(I)
       CE=MIN( FLCF,R_0_0)
       CW=MAX( FLCF,R_0_0)
!
       A_E(I)=-DE+CE
!
!.....
!.....        D I F F U S I O N
!.....
       SUEH=GAME*((GR1X(INP)*FXW+GR1X(IND)*FXE)*ARX+ &
                  (GR1Y(INP)*FXW+GR1Y(IND)*FXE)*ARY+ &
                  (GR1Z(INP)*FXW+GR1Z(IND)*FXE)*ARZ)
!.....
       SUEL=DE*((GR1X(INP)*FXW+GR1X(IND)*FXE)*AKX+ &
                (GR1Y(INP)*FXW+GR1Y(IND)*FXE)*AKY+ &
                (GR1Z(INP)*FXW+GR1Z(IND)*FXE)*AKZ)
!.....
!.....     C O N V E C T I O N
!.....
       UHE=FLCF*(FI(IND)*FXE+FI(INP)*FXW)
!.....
       ULE=CW*FI(INP)+CE*FI(IND)
!..... PART TO ADD
       SUC=-GAM*(UHE-ULE)
!

       SUADD=SUC+SUEH-SUEL
!
       SU(INP)=SU(INP)+SUADD
!
       SP (INP)= SP(INP)-A_E(I)
!
      END IF
    END DO
  END  SUBROUTINE COMSC

!**************************************************
  SUBROUTINE METRIC
!**************************************************
    IMPLICIT NONE
!
    CALL CALCVOL(X,Y,Z,VOL)
    CALL CHVSCA(VOL)

    CALL CALCG

    CALL CHVVEC(XC,YC,ZC)
!
    CALL CALFAC
!
  END SUBROUTINE METRIC

!**************************************************
  SUBROUTINE CALCG
!**************************************************
    IMPLICIT NONE

    INTEGER(IK) :: NBL,I,J,K
    INTEGER(IK) :: INP,INS,INB,INBS,INW,INSW,INBW,INBSW
!.....
!
!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)

       DO I=2,NIM
          DO J=2,NJM
             DO K=2,NKM
                INP=LK(K)+LI(I)+J
                INS=INP-1
                INB=INP-NIJ
                INBS=INB-1
!
                INW=INP-NJ
                INSW=INW-1
                INBW=INW-NIJ
                INBSW=INBW-1
!....
                XC(INP)=.125*(X(INP)+X(INS )+X(INB )+X(INBS)+ &
                              X(INW)+X(INSW)+X(INBW)+X(INBSW))
                YC(INP)=.125*(Y(INP)+Y(INS )+Y(INB )+Y(INBS)+ &
                              Y(INW)+Y(INSW)+Y(INBW)+Y(INBSW))
                ZC(INP)=.125*(Z(INP)+Z(INS )+Z(INB )+Z(INBS)+ &
                              Z(INW)+Z(INSW)+Z(INBW)+Z(INBSW))
             END DO
          END DO
       END DO

!.....THE REST
!.....WEST
       DO K=2,NKM
          DO J=2,NJM
             INP=LK(K)+LI(1)+J
             INS=INP-1
             INB=INP-NIJ
             INBS=INB-1
!
             XC(INP)=.25*(X(INP)+X(INS )+X(INB )+X(INBS))
             YC(INP)=.25*(Y(INP)+Y(INS )+Y(INB )+Y(INBS))
             ZC(INP)=.25*(Z(INP)+Z(INS )+Z(INB )+Z(INBS))
          END DO
       END DO
!.....EAST
       DO K=2,NKM
          DO J=2,NJM
             INP=LK(K)+LI(NIM)+J
             INS=INP-1
             INB=INP-NIJ
             INBS=INB-1
!
             XC(INP+NJ)=.25*(X(INP)+X(INS )+X(INB )+X(INBS))
             YC(INP+NJ)=.25*(Y(INP)+Y(INS )+Y(INB )+Y(INBS))
             ZC(INP+NJ)=.25*(Z(INP)+Z(INS )+Z(INB )+Z(INBS))
          END DO
       END DO
!
!.....SOUTH
       DO K=2,NKM
          DO I=2,NIM
             INP=LK(K)+LI(I)+1
             INW=INP-NJ
             INB=INP-NIJ
             INBW=INB-NJ
!
             XC(INP)=.25*(X(INP)+X(INB )+X(INW )+X(INBW))
             YC(INP)=.25*(Y(INP)+Y(INB )+Y(INW )+Y(INBW))
             ZC(INP)=.25*(Z(INP)+Z(INB )+Z(INW )+Z(INBW))
          END DO
       END DO
!
!.....NORTH
       DO K=2,NKM
          DO I=2,NIM
             INP=LK(K)+LI(I)+NJM
             INW=INP-NJ
             INB=INP-NIJ
             INBW=INB-NJ
!
             XC(INP+1)=.25*(X(INP)+X(INB )+X(INW )+X(INBW))
             YC(INP+1)=.25*(Y(INP)+Y(INB )+Y(INW )+Y(INBW))
             ZC(INP+1)=.25*(Z(INP)+Z(INB )+Z(INW )+Z(INBW))
          END DO
       END DO
!
!.....BOTOM
       DO I=2,NIM
          DO J=2,NJM
             INP=LK(1)+LI(I)+J
             INW=INP-NJ
             INS=INP-1
             INSW=INS-NJ
!
             XC(INP)=.25*(X(INP)+X(INS )+X(INW )+X(INSW))
             YC(INP)=.25*(Y(INP)+Y(INS )+Y(INW )+Y(INSW))
             ZC(INP)=.25*(Z(INP)+Z(INS )+Z(INW )+Z(INSW))
          END DO
       END DO
!
!.....TOP
       DO I=2,NIM
          DO J=2,NJM
             INP=LK(NKM)+LI(I)+J
             INW=INP-NJ
             INS=INP-1
             INSW=INS-NJ
!
             XC(INP+NIJ)=.25*(X(INP)+X(INS )+X(INW )+X(INSW))
             YC(INP+NIJ)=.25*(Y(INP)+Y(INS )+Y(INW )+Y(INSW))
             ZC(INP+NIJ)=.25*(Z(INP)+Z(INS )+Z(INW )+Z(INSW))
          END DO
       END DO
!
!.....PART  1 ...CORNERS

       I = 1
       J = 1
       K = 1
       INP=LK(K)+LI(I)+J

       XC(INP)= X(INP)
       YC(INP)= Y(INP)
       ZC(INP)= Z(INP)

       I = 1
       J = NJ
       K = 1
       INP=LK(K)+LI(I)+J

       XC(INP)= X(INP-1)
       YC(INP)= Y(INP-1)
       ZC(INP)= Z(INP-1)

       I = 1
       J = NJ
       K = NK
       INP=LK(K)+LI(I)+J

       XC(INP)= X(INP-1-NIJ)
       YC(INP)= Y(INP-1-NIJ)
       ZC(INP)= Z(INP-1-NIJ)

       I = 1
       J = 1
       K = NK
       INP=LK(K)+LI(I)+J

       XC(INP)= X(INP-NIJ)
       YC(INP)= Y(INP-NIJ)
       ZC(INP)= Z(INP-NIJ)

       I = NI
       J = 1
       K = 1
       INP=LK(K)+LI(I)+J

       XC(INP)= X(INP-NJ)
       YC(INP)= Y(INP-NJ)
       ZC(INP)= Z(INP-NJ)

       I = NI
       J = NJ
       K = 1
       INP=LK(K)+LI(I)+J

       XC(INP)= X(INP-NJ-1)
       YC(INP)= Y(INP-NJ-1)
       ZC(INP)= Z(INP-NJ-1)

       I = NI
       J = NJ
       K = NK
       INP=LK(K)+LI(I)+J

       XC(INP)= X(INP-NJ-1-NIJ)
       YC(INP)= Y(INP-NJ-1-NIJ)
       ZC(INP)= Z(INP-NJ-1-NIJ)

       I = NI
       J = 1
       K = NK
       INP=LK(K)+LI(I)+J

       XC(INP)= X(INP-NJ-NIJ)
       YC(INP)= Y(INP-NJ-NIJ)
       ZC(INP)= Z(INP-NJ-NIJ)

!.....PART  2 ...LINE POINTS

!
!.....WEST LINES
!
       DO J=2,NJM

          K = 2

          INP=LK(K)+LI(1)+J
          INS=INP-1
          INB=INP-NIJ
          INBS=INB-1
!
          XC(INP-NIJ)=.5*(X(INB )+X(INBS))
          YC(INP-NIJ)=.5*(Y(INB )+Y(INBS))
          ZC(INP-NIJ)=.5*(Z(INB )+Z(INBS))

          K = NKM

          INP=LK(K)+LI(1)+J
          INS=INP-1
          INB=INP-NIJ
          INBS=INB-1
!
          XC(INP+NIJ)=.5*(X(INP)+X(INS))
          YC(INP+NIJ)=.5*(Y(INP)+Y(INS))
          ZC(INP+NIJ)=.5*(Z(INP)+Z(INS))

       END DO

       DO K=2,NKM

          J=2

          INP=LK(K)+LI(1)+J
          INS=INP-1
          INB=INP-NIJ
          INBS=INB-1
!
          XC(INP-1)=.5*(X(INS)+X(INBS))
          YC(INP-1)=.5*(Y(INS)+Y(INBS))
          ZC(INP-1)=.5*(Z(INS)+Z(INBS))

          J=NJM

          INP=LK(K)+LI(1)+J
          INS=INP-1
          INB=INP-NIJ
          INBS=INB-1
!
          XC(INP+1)=.5*(X(INP)+X(INB))
          YC(INP+1)=.5*(Y(INP)+Y(INB))
          ZC(INP+1)=.5*(Z(INP)+Z(INB))

       END DO

!.....EAST
       DO J=2,NJM

          K = 2

          INP=LK(K)+LI(NIM)+J
          INS=INP-1
          INB=INP-NIJ
          INBS=INB-1
!
          XC(INP+NJ-NIJ)=.5*(X(INB )+X(INBS))
          YC(INP+NJ-NIJ)=.5*(Y(INB )+Y(INBS))
          ZC(INP+NJ-NIJ)=.5*(Z(INB )+Z(INBS))

          K = NKM
          INP=LK(K)+LI(NIM)+J

          INS=INP-1
          INB=INP-NIJ
          INBS=INB-1
!
          XC(INP+NJ+NIJ)=.5*(X(INP)+X(INS))
          YC(INP+NJ+NIJ)=.5*(Y(INP)+Y(INS))
          ZC(INP+NJ+NIJ)=.5*(Z(INP)+Z(INS))

       END DO

       DO K=2,NKM

          J=2

          INP=LK(K)+LI(NIM)+J
          INS=INP-1
          INB=INP-NIJ
          INBS=INB-1
!
          XC(INP+NJ-1)=.5*(X(INS)+X(INBS))
          YC(INP+NJ-1)=.5*(Y(INS)+Y(INBS))
          ZC(INP+NJ-1)=.5*(Z(INS)+Z(INBS))

          J=NJM

          INP=LK(K)+LI(NIM)+J
          INS=INP-1
          INB=INP-NIJ
          INBS=INB-1
!
          XC(INP+NJ+1)=.5*(X(INP)+X(INB))
          YC(INP+NJ+1)=.5*(Y(INP)+Y(INB))
          ZC(INP+NJ+1)=.5*(Z(INP)+Z(INB))

       END DO
!
!.....SOUTH
       DO I=2,NIM
          K=2
          INP=LK(K)+LI(I)+1
          INW=INP-NJ
          INB=INP-NIJ
          INBW=INB-NJ
!
          XC(INP-NIJ)=.5*(X(INB)+X(INBW))
          YC(INP-NIJ)=.5*(Y(INB)+Y(INBW))
          ZC(INP-NIJ)=.5*(Z(INB)+Z(INBW))

          K=NKM
          INP=LK(K)+LI(I)+1
          INW=INP-NJ
          INB=INP-NIJ
          INBW=INB-NJ
!
          XC(INP+NIJ)=.5*(X(INP)+X(INW))
          YC(INP+NIJ)=.5*(Y(INP)+Y(INW))
          ZC(INP+NIJ)=.5*(Z(INP)+Z(INW))

       END DO
!
!.....NORTH
!
       DO I=2,NIM

          K =2

          INP=LK(K)+LI(I)+NJM
          INW=INP-NJ
          INB=INP-NIJ
          INBW=INB-NJ
!
          XC(INP+1-NIJ)=.5*(X(INB)+X(INBW))
          YC(INP+1-NIJ)=.5*(Y(INB)+Y(INBW))
          ZC(INP+1-NIJ)=.5*(Z(INB)+Z(INBW))

          K =NKM

          INP=LK(K)+LI(I)+NJM
          INW=INP-NJ
          INB=INP-NIJ
          INBW=INB-NJ
!
          XC(INP+1+NIJ)=.5*(X(INP)+X(INW))
          YC(INP+1+NIJ)=.5*(Y(INP)+Y(INW))
          ZC(INP+1+NIJ)=.5*(Z(INP)+Z(INW))

       END DO

!
!.....END BLOCK
!
    END DO
!
  END SUBROUTINE CALCG

!**************************************************
  SUBROUTINE CALFAC
!**************************************************
    IMPLICIT NONE
    INTEGER(IK) :: NBL
!
!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!
!
!.....KSI - DIRECTION
       CALL INTFAC(1,2,2,NJ,1,NIJ,FX)
!.....ETA - DIRECTION
       CALL INTFAC(2,1,2,1,NJ,NIJ,FY)
!.....ZETA - DIRECTION
       CALL INTFAC(2,2,1,NIJ,1,NJ,FZ)
!
!.....BLOCK INTERFACES
       CALL INTCF(NBL)
    END DO
  END SUBROUTINE CALFAC

!************************************************************
  SUBROUTINE INTFAC(NIS,NJS,NKS,IDEW,IDNS,IDTB,FINT)
!************************************************************
    IMPLICIT NONE

    REAL(RK)   ,INTENT(INOUT) :: FINT(NXYZA)
    INTEGER(IK),INTENT(IN)  :: NIS,NJS,NKS,IDEW,IDNS,IDTB

    REAL(RK)    :: DXP,DYP,DZP,DLP,DLE
    INTEGER(IK) :: I,J,K,INP,INE,INS,INB,INBS

!.....FX INTERPOLATION FACTOR
    DO K=NKS,NKM
       DO I=NIS,NIM
          DO J=NJS,NJM
             INP=LK(K)+LI(I)+J
             INE=INP+IDEW
             INS=INP-IDNS
             INB=INP-IDTB
             INBS=INB-IDNS
!.....COORDINATES OF THE POINTS AT THE CENTRE OF CF
!.....CENTRAL (EAST CV FACE)
             DXP=0.25*(X(INP)+X(INS )+X(INB )+X(INBS ))
             DYP=0.25*(Y(INP)+Y(INS )+Y(INB )+Y(INBS ))
             DZP=0.25*(Z(INP)+Z(INS )+Z(INB )+Z(INBS ))
!.....LENGHT OF THE LINES ON THE POINTS P AND E
             DLP=SQRT((DXP-XC(INP))**2+(DYP-YC(INP))**2+(DZP-ZC(INP))**2)
             DLE=SQRT((DXP-XC(INE))**2+(DYP-YC(INE))**2+(DZP-ZC(INE))**2)
             FINT(INP)=DLP/(DLP+DLE)
          END DO
       END DO
    END DO
  END SUBROUTINE INTFAC

!*****************************************
  SUBROUTINE INTCF(NBL)
!*****************************************
    IMPLICIT NONE

    INTEGER(IK),INTENT(IN) :: NBL
    INTEGER(IK) :: I,INP,IND
    REAL(RK)    :: DXP,DYP,DZP,DLP,DLE
!
    DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
       INP=NBL_ST(NBL) + NUM_CFI(I,1)
       IND=NBL_ST(NBL) + NUM_CFI(I,3)

!.....INTERPOLATION FACTOR
!.....COORDINATES OF THE POINTS AT THE CENTRE OF CF
!.....CENTRAL (EAST CV FACE)
       DXP=0.25*(X_CFI(I,1)+X_CFI(I,2)+X_CFI(I,4)+X_CFI(I,3))
       DYP=0.25*(Y_CFI(I,1)+Y_CFI(I,2)+Y_CFI(I,4)+Y_CFI(I,3))
       DZP=0.25*(Z_CFI(I,1)+Z_CFI(I,2)+Z_CFI(I,4)+Z_CFI(I,3))

!.....LENGHT OF THE LINES ON THE POINTS P AND E
       DLP=SQRT((DXP-XC(INP))**2+(DYP-YC(INP))**2+(DZP-ZC(INP))**2)
       DLE=SQRT((DXP-XC(IND))**2+(DYP-YC(IND))**2+(DZP-ZC(IND))**2)
!
       FX_CFI(I)=DLP/(DLP+DLE)
    END DO
  END SUBROUTINE INTCF

!**************************************************
  SUBROUTINE CALCVOL(XT,YT,ZT,VOLTMP)
!**************************************************

    IMPLICIT NONE
    REAL(RK)   ,INTENT(IN)  :: XT(NXYZA),YT(NXYZA),ZT(NXYZA)
    REAL(RK)   ,INTENT(OUT) :: VOLTMP(NXYZA)

        REAL(RK)    :: SIXR,XA,YA,ZA,VOLUM
    REAL(RK)    :: DXAB,DYAB,DZAB
    REAL(RK)    :: DXAC,DYAC,DZAC
    REAL(RK)    :: DXAD,DYAD,DZAD
    REAL(RK)    :: DXAE,DYAE,DZAE
    REAL(RK)    :: DXAF,DYAF,DZAF
    REAL(RK)    :: DXAG,DYAG,DZAG
    REAL(RK)    :: DXAH,DYAH,DZAH
    INTEGER(IK) :: NBL,I,J,K,INPHLP
    INTEGER(IK) :: IJK,IMJK,IMJMK,IJKM,IJKA,IJMK
    LOGICAL     :: LHELP
!
    LHELP=.FALSE.
!
!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!
!.....CALCULATION OF CELL VOLUMES
       SIXR=1./6.

       DO K=2,NKM
          DO I=2,NIM
             DO J=2,NJM
                IJK=LK(K)+LI(I)+J

                IMJK=IJK-NJ
                IMJMK=IMJK-1
                IJMK=IJK-1
                IJKM=IJK-NIJ
                IJKA=IJKM-NJ-1
!
                XA=XT(IJKA)
                YA=YT(IJKA)
                ZA=ZT(IJKA)
                DXAD=XT(IJMK)-XA
                DYAD=YT(IJMK)-YA
                DZAD=ZT(IJMK)-ZA
                DXAB=XT(IJKM-1)-XA
                DYAB=YT(IJKM-1)-YA
                DZAB=ZT(IJKM-1)-ZA
                DXAC=XT(IMJMK)-XA
                DYAC=YT(IMJMK)-YA
                DZAC=ZT(IMJMK)-ZA
                DXAE=XT(IJKM-NJ)-XA
                DYAE=YT(IJKM-NJ)-YA
                DZAE=ZT(IJKM-NJ)-ZA
                DXAF=XT(IJKM)-XA
                DYAF=YT(IJKM)-YA
                DZAF=ZT(IJKM)-ZA
                DXAG=XT(IMJK)-XA
                DYAG=YT(IMJK)-YA
                DZAG=ZT(IMJK)-ZA
                DXAH=XT(IJK)-XA
                DYAH=YT(IJK)-YA
                DZAH=ZT(IJK)-ZA
!
                VOLUM=VOL1(DXAD,DYAD,DZAD,DXAB,DYAB,DZAB,DXAH,DYAH,DZAH)+ &
                      VOL1(DXAC,DYAC,DZAC,DXAD,DYAD,DZAD,DXAH,DYAH,DZAH)+ &
                      VOL1(DXAG,DYAG,DZAG,DXAC,DYAC,DZAC,DXAH,DYAH,DZAH)+ &
                      VOL1(DXAB,DYAB,DZAB,DXAF,DYAF,DZAF,DXAH,DYAH,DZAH)+ &
                      VOL1(DXAF,DYAF,DZAF,DXAE,DYAE,DZAE,DXAH,DYAH,DZAH)+ &
                      VOL1(DXAE,DYAE,DZAE,DXAG,DYAG,DZAG,DXAH,DYAH,DZAH)
                VOLTMP(IJK)=ABS(VOLUM*SIXR)
!.....CHECK IF NEGATIV (on prend la valeur absolue pour geometries courbes)
                IF(VOLTMP(IJK)<=0.) THEN
                   WRITE(*,*) ' NBL=',NBL,I,J,K,IJK,' VOLUM NUL'
                   LHELP=.TRUE.
                END IF
             END DO
          END DO
       END DO

       IF(LHELP)  CALL END(' CHECK IT OUT NEGATIV VOLUME ')
    END DO
  600 FORMAT(2X,I3,3X,1P,(3E13.6,2X))
  END SUBROUTINE CALCVOL

  FUNCTION  VOL1(A1,A2,A3,B1,B2,B3,Q1,Q2,Q3)
    REAL(RK)    :: VOL1,A1,A2,A3,B1,B2,B3,Q1,Q2,Q3

    VOL1=(A2*B3-B2*A3)*Q1 &
        +(B1*A3-A1*B3)*Q2 &
        +(A1*B2-A2*B1)*Q3
  END FUNCTION VOL1

!*****************************************
  SUBROUTINE GRADFI(GRX,GRY,GRZ,FI,id)
!*****************************************
    USE GEOM
    IMPLICIT NONE

    REAL(RK),INTENT(INOUT) :: FI(NXYZA)
    REAL(RK),INTENT(INOUT) :: GRX(NXYZA),GRY(NXYZA),GRZ(NXYZA)
    INTEGER(IK),INTENT(IN),OPTIONAL ::id

    REAL(RK)    :: FIE,VOLR
    INTEGER(IK) :: INP,NBL,IDR,NBC_S,NBC_E,I,J,K,ISIDE,IND
    INTEGER(IK)  :: id_pass
!.....EXCHANGE VALUES ON BLOCK INTERFACES

		if(PRESENT(id)) THEN
   			CALL CHVSCA(FI,alias=id)
                            id_pass=id
		else
			CALL CHVSCA(FI)
                             id_pass=0
		end if 


!.....CLEAR
    GRX=0.
    GRY=0.
    GRZ=0.
!
!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!
     CALL GRADE(NI-2,NJM,NKM,  NJ,  1,NIJ,FI,FX,GRX,GRY,GRZ,2)
     CALL GRADE(NIM,NJ-2,NKM,   1,NIJ, NJ,FI,FY,GRX,GRY,GRZ,1)
     CALL GRADE(NIM,NJM,NK-2, NIJ, NJ,  1,FI,FZ,GRX,GRY,GRZ,3)
!
!...  COLLECT SIDES OF THE BLOCKS
!...  ADD HERE
!.....ALL BC
     DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)
!.....
        NBC_S=NUM_SPR(IDR)+1
        NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)
        DO I=NBC_S,NBC_E
!
           INP=NBL_ST(NBL)+NUM_IND(I,1)
           ISIDE=NUM_IND(I,2)
!
           CALL INIT_GEOM_BF(ISIDE,INP)
!
           FIE=FI(INE)
           GRX(INP)=GRX(INP)+FIE*ARX
           GRY(INP)=GRY(INP)+FIE*ARY
           GRZ(INP)=GRZ(INP)+FIE*ARZ
        END DO
     END DO
!
!.....BLOCK INTERFACES
     DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
        INP=NBL_ST(NBL) + NUM_CFI(I,1)
        IND=NBL_ST(NBL) + NUM_CFI(I,3)
        ISIDE=NUM_CFI(I,2)
!
 IF(NUM_CFI(I,7)==0.OR.id_pass==7) THEN 

        CALL INIT_GEOM_BLOCK(ISIDE,I)
!
!.....INTERPOLATION FACTOR
        FXE=FX_CFI(I)
        FXW=1.-FXE

        FIE=FI(INP)*FXW+FI(IND)*FXE

        GRX(INP)=GRX(INP)+FIE*ARX
        GRY(INP)=GRY(INP)+FIE*ARY
        GRZ(INP)=GRZ(INP)+FIE*ARZ
 end if
     END DO
!
     DO K=2,NKM
        DO I=2,NIM
           DO J=2,NJM
              INP=LK(K)+LI(I)+J
              VOLR=1./VOL(INP)
              GRX(INP)=GRX(INP)*VOLR
              GRY(INP)=GRY(INP)*VOLR
              GRZ(INP)=GRZ(INP)*VOLR
           END DO
        END DO
     END DO
    END DO
!
!.....EXCHANGE VALUES ON BLOCK INTERFACES
    CALL CHVVEC(GRX,GRY,GRZ)
  END SUBROUTINE GRADFI

!********************************************************************
  SUBROUTINE GRADE(NIE,NJE,NKE,IDEW,IDNS,IDTB,FI,FIF,GRX,GRY,GRZ,DIR1)
!********************************************************************
    IMPLICIT NONE

    REAL(RK)   ,INTENT(INOUT) :: FI(NXYZA),FIF(NXYZA)
    REAL(RK)   ,INTENT(INOUT) :: GRX(NXYZA),GRY(NXYZA),GRZ(NXYZA)
    INTEGER(IK),INTENT(IN)    :: NIE,NJE,NKE,IDEW,IDNS,IDTB,DIR1

    REAL(RK)    :: BN,BE,BT
    REAL(RK)    :: FIE
    INTEGER(IK) :: I,J,K

    REAL(RK)    :: FXE,FXW
    INTEGER(IK) :: INP,INE

    DO K=2,NKE
       DO I=2,NIE
          DO J=2,NJE
!
             INP=LK(K)+LI(I)+J
             INE=INP+IDEW
!
!.....INTERPOLATION FACTOR
             FXE=FIF(INP)
             FXW=1.-FXE

             FIE=FI(INP)*FXW+FI(INE)*FXE
             BE= FIE*ARX_T(INP,DIR1)
             BN= FIE*ARY_T(INP,DIR1)
             BT= FIE*ARZ_T(INP,DIR1)

             GRX(INP)=GRX(INP)+BE
             GRY(INP)=GRY(INP)+BN
             GRZ(INP)=GRZ(INP)+BT
             GRX(INE)=GRX(INE)-BE
             GRY(INE)=GRY(INE)-BN
             GRZ(INE)=GRZ(INE)-BT
          END DO
       END DO
    END DO

  END SUBROUTINE GRADE
  
  
!*********************************************************************
  SUBROUTINE BPRES(PPP,id)
!*********************************************************************
    USE GEOM
    
    IMPLICIT NONE
         
    REAL(RK),INTENT(INOUT) :: PPP(NXYZA)
    INTEGER(IK),INTENT(IN),OPTIONAL :: id
    INTEGER(IK)            :: NBL,I,J,K,INP
    INTEGER(IK)            :: IDR,NBC_E,NBC_S,ISIDE

!.....SET BOUNDARY PRESSURES
!
!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!
!.....BOTTOM BOUNDARY
       DO I=2,NIM
          DO J=2,NJM
             INP=LK(2)+LI(I)+J
             PPP(INP-NIJ)=PPP(INP)-(PPP(INP+NIJ)-PPP(INP))*FZ(INP)
          END DO
       END DO

!.....TOP BOUNDARY
       DO I=2,NIM
          DO J=2,NJM
             INP=LK(NKM)+LI(I)+J
             PPP(INP+NIJ)=PPP(INP)+(PPP(INP)-PPP(INP-NIJ))*(1.-FZ(INP-NIJ))
          END DO
       END DO

!.....WEST BOUNDARY
       DO K=1,NK
          DO J=2,NJM
             INP=LK(K)+LI(2)+J
             PPP(INP-NJ)=PPP(INP)-(PPP(INP+NJ)-PPP(INP))*FX(INP)
          END DO
       END DO

!.....EAST BOUNDARY
       DO K=1,NK
          DO J=2,NJM
             INP=LK(K)+LI(NIM)+J
             PPP(INP+NJ)=PPP(INP)+(PPP(INP)-PPP(INP-NJ))*(1.-FX(INP-NJ))
          END DO
       END DO

!.....SOUTH BOUNDARY
       DO I=1,NI
          DO K=1,NK
             INP=LK(K)+LI(I)+2
             PPP(INP-1)=PPP(INP)-(PPP(INP+1)-PPP(INP))*FY(INP)
          END DO
       END DO

!.....NORTH BOUNDARY
       DO I=1,NI
          DO K=1,NK
             INP=LK(K)+LI(I)+NJM
             PPP(INP+1)=PPP(INP)+(PPP(INP)-PPP(INP-1))*(1.-FY(INP-1))
          END DO
       END DO
!
    END DO



!*********************************************************************
!..... SET VALUES AT SYMM BOUNDARIES ...   T Y P    3               *
!*********************************************************************


      DO NBL=1,NBLOCK
      CALL SETIND(NBL) 
!
!.....SET FOR  DIFFERENT REGIONS
      DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)
!

      IF(NUM_TYP(IDR)==3) THEN


      NBC_S=NUM_SPR(IDR)+1
      NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)
! 	   
      DO I=NBC_S,NBC_E 
      INP=NBL_ST(NBL)+NUM_IND(I,1)
      ISIDE=NUM_IND(I,2)
!
      CALL INIT_GEOM_BF(ISIDE,INP)
!       
      PPP(INE)= PPP(INP)      ! first order
    

      ENDDO

     
      ENDIF    ! END BF TYPE 3
!
!.....END REGIONS                   
      ENDDO

!.....END BLOCKS
      ENDDO

!.....EXCHANGE VALUES AT BLOCK INTERFACES
    IF(PRESENT(id)) THEN
            CALL CHVSCA(PPP,alias=id)
    ELSE
            CALL CHVSCA(PPP)
    END IF

  END SUBROUTINE BPRES

   
   

!*****************************************
  SUBROUTINE SIPSOL(FI,IFI)
!*****************************************
    IMPLICIT NONE

    REAL(RK),ALLOCATABLE   ,INTENT(INOUT) :: FI(:)
    INTEGER(IK),INTENT(IN)    :: IFI

    REAL(RK),ALLOCATABLE      :: BB(:),BS(:),BW(:),BP(:)
    REAL(RK),ALLOCATABLE      :: BN(:),BE(:),BT(:),su1(:)
    REAL(RK)                  :: P1,P2,P3,RES1,RSM,RES2,VALUES(7)
    INTEGER(IK)               :: INP,NBL,NBL1,I,J,K,IJK,IMJK,IJKM,IJMK
    INTEGER(IK)               :: NS,N,IND,ROW,COLS(7)

    ALLOCATE(BB(NXYZA),BS(NXYZA),BW(NXYZA),BP(NXYZA))
    ALLOCATE(BN(NXYZA),BE(NXYZA),BT(NXYZA))
!.....CLEAR

    BW=0.
    BE=0.
    BN=0.
    BS=0.
    BB=0.
    BT=0.

!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!
!.....CALCULATE COEFFICIENTS OF  L  AND  U  MATRICES
       DO K=2, NKM
          DO I=2, NIM
             DO J=2, NJM
                IJK=LI(I)+LK(K)+J
                IMJK=IJK-NJ
                IJKM=IJK-NIJ
                IJMK=IJK-1
                BB(IJK)=+AB(IJK)/(1.+ALFA*(BN(IJKM)+BE(IJKM)))
                BW(IJK)=+AW(IJK)/(1.+ALFA*(BN(IMJK)+BT(IMJK)))
                BS(IJK)=+AS(IJK)/(1.+ALFA*(BE(IJMK)+BT(IJMK)))
                P1=ALFA*(BB(IJK)*BN(IJKM)+BW(IJK)*BN(IMJK))
                P2=ALFA*(BB(IJK)*BE(IJKM)+BS(IJK)*BE(IJMK))
                P3=ALFA*(BW(IJK)*BT(IMJK)+BS(IJK)*BT(IJMK))
                BP(IJK)=1./(AP(IJK)+P1+P2+P3-BB(IJK)*BT(IJKM) &
                                            -BW(IJK)*BE(IMJK) &
                                            -BS(IJK)*BN(IJMK))
                BN(IJK)=(+AN(IJK)-P1)*BP(IJK)
                BE(IJK)=(+AE(IJK)-P2)*BP(IJK)
                BT(IJK)=(+AT(IJK)-P3)*BP(IJK)
 

             END DO
          END DO
       END DO
!.....END ALL BLOCKS
     END DO
!
    

!..... EXCHANGE SOLUTION AT INTERFACES
!
if(IFI==7) then
    CALL CHVSCA(FI,alias=IFI)
  !CALL CHVSCA(FI)
else
CALL CHVSCA(FI)
end if
!

!.....ITERATION LOOP
    NS=NSW(IFI)
    DO N=1,NS
!.....CALCULATE RESIDUALS AND AUXILLIARY VECTOR
       RES1=0.0
!
!.....ALL BLOCKS
!
       DO NBL=1,NBLOCK
          CALL SETIND(NBL)

        DO K=2,NKM
           DO I=2,NIM
              DO J=2,NJM
                 IJK =  LK(K)+LI(I)+J

                 RES(IJK)=SU(IJK)-AP(IJK)*FI(IJK)                        &
                                 -AE(IJK)*FI(IJK+NJ) -AW(IJK)*FI(IJK-NJ) &
                                 -AN(IJK)*FI(IJK+1)  -AS(IJK)*FI(IJK-1)  &
                                 -AT(IJK)*FI(IJK+NIJ)-AB(IJK)*FI(IJK-NIJ)
              END DO
           END DO
        END DO

!-----.....  BLOCK INTERFACES  ....---
        DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
           INP=NBL_ST(NBL) + NUM_CFI(I,1)
           IND=NBL_ST(NBL) + NUM_CFI(I,3)

   IF(NUM_CFI(I,7)==0.OR.IFI==7) THEN 

           RES(INP)=RES(INP)-A_E(I)*FI(IND)
!           RES1=RES1+ABS(RES(INP))
   end if
        END DO

        DO K=2,NKM
           DO I=2,NIM
              DO J=2,NJM
                 IJK =  LK(K)+LI(I)+J
                 RES1=RES1+ABS(RES(IJK))
                 RES(IJK)=(RES(IJK)-BB(IJK)*RES(IJK-NIJ) &
                                   -BW(IJK)*RES(IJK-NJ)  &
                                   -BS(IJK)*RES(IJK-1))*BP(IJK)
              END DO
           END DO
        END DO
!.....END ALL BLOCKS
       END DO
!
!..............................................................
!.....CALCULATE INCREMENT AND UPDATE VARIABLES
!
!.....ALL BLOCKS
       DO NBL=1,NBLOCK
          CALL SETIND(NBL)
          DO K=NKM,2,-1
             DO I=NIM,2,-1
                DO J=NJM,2,-1
                   IJK =  LK(K)+LI(I)+J

                   RES(IJK)=RES(IJK)-BN(IJK)*RES(IJK+1)  &
                                    -BE(IJK)*RES(IJK+NJ) &
                                    -BT(IJK)*RES(IJK+NIJ)

                   FI(IJK)=FI(IJK)+RES(IJK)
                END DO
             END DO
          END DO

!.....END ALL BLOCKS
       END DO
!
!..... EXCHANGE SOLUTION AT INTERFACES
!
       if(IFI==7) then
           CALL CHVSCA(FI,alias=IFI)
      !CALL CHVSCA(FI)
       else
          CALL CHVSCA(FI)
       end if

       IF(N==1) RESOR(IFI)=RES1
       RSM=RES1/(RESOR(IFI)+tiny(1.))
!
!.....CHECK CONVERGENCE OF INNER ITERATIONS
!
       IF(LTEST) WRITE(IUCHK,600) CHVAR(IFI),N,RES1,RSM
       IF(RSM<SOR(IFI)) EXIT
    END DO
    RESOR(IFI)=RES1
    DEALLOCATE(BB,BS,BW,BP,BN,BE,BT)

  600 FORMAT(20X,'FI=',A3,'  SWEEP=',I4,'  RES=',1PE10.3, &
                 ' RSM=',1PE10.3)
  END SUBROUTINE SIPSOL

!*********************************************************************
  SUBROUTINE CALNOR
!*********************************************************************
    IMPLICIT NONE
    INTEGER(IK) :: NBL,IDR,NBC_E,NBC_S,I,INP,ISIDE,IDEW,INE

!.....CALCULATE NORMALISATION FACTORS
!
!.....ALL BLOCKS
!
    FLOWIN=0.
    XMONIN=0.
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!
!.....CHECK - TAKE ONLY INLET REGION
       DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)
!
!..... T Y P  INLET
          IF(NUM_TYP(IDR)==1) THEN
             NBC_S=NUM_SPR(IDR)+1
             NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)
!
             DO I=NBC_S,NBC_E
                INP=NBL_ST(NBL)+NUM_IND(I,1)
                ISIDE=NUM_IND(I,2)
                IDEW=ID_SIDE  (ID_ODER(ISIDE,1))*ID_SIGN(ISIDE)
!
!.....CALCULATE INDEX OF BOUNDARY POINT
                INE=INP+IDEW
!
                FLOWIN=FLOWIN-FLX_IND(I)
                XMONIN=XMONIN-FLX_IND(I)*SQRT(U(INE)**2+V(INE)**2+W(INE)**2)
             END DO
          END IF
!
!.....END REGIONS
       END DO
!
!.....END BLOCKS
    END DO
!
!.....NORMALIZATION FACTORS
    IF(FLOWIN<R_0_0) THEN
       WRITE(*,600)
       CALL END
    END IF
    IF(FLOWIN<SMALL) THEN
       FLOWIN=1.
       XMONIN=1.
    END IF
    RESOR(1:NPHI)=0.
    SNORIN(IU)=1./(XMONIN)
    SNORIN(IV)=SNORIN(IU)
    SNORIN(IW)=SNORIN(IU)
    SNORIN(IP)=1./(FLOWIN)
    SNORIN(IEN)=1.
!
  600 FORMAT(1X,' MASS INFLOW NEGATIV,  STOP')
  END SUBROUTINE CALNOR


!******************************
  SUBROUTINE UPDOLD      !*
!******************************
    IMPLICIT NONE
!
    UOO=UO
    VOO=VO
    WOO=WO
!..
    UO=U
    VO=V
    WO=W
!
    IF(LCAL(IEN).OR.LETHD) THEN
       TOO=TO
       TO=T
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Modif
    IF((LEHD.OR.LETHD).AND.ABS(D_C)>SMALL) THEN
       QOO = QO
       QO  = Q
       VPO = VP

       EXO = EX
       EYO = EY
       EZO = EZ
    END IF
  END SUBROUTINE UPDOLD

!**************************************************
  SUBROUTINE FINSID
!**************************************************
    IMPLICIT NONE
!
    INTEGER(IK) :: NBL1,I1,INPB1,INPB1T2,ISIB1
    INTEGER(IK) :: NBL2,I2,INPB2,INPB2T1,ISIB2,mini
    LOGICAL     :: LTOUCH(NUM_CF_ALL)
    REAL(RK)    :: x1,x2,y1,y2,z1,z2,dist,mindist

if (.false.) then
! TRY to FIX CONNECTIONS IGNORING  NUM_CFI(I1,4) THIS SHOULDN'T BE USED !
    DO NBL1=1,NBLOCK
       DO I1=NUM_SCF(NBL1)+1,NUM_SCF(NBL1)+NUM_CF(NBL1)
         INPB1  =NBL_ST(NBL1)+NUM_CFI(I1,1)
         x1=x(INPB1)
         y1=y(INPB1)
         z1=z(INPB1)
         NBL2=NUM_CFI(I1,5)
         mindist=huge(1)
         DO I2=NUM_SCF(NBL2)+1,NUM_SCF(NBL2)+NUM_CF(NBL2)
             INPB2  =NBL_ST(NBL2)+NUM_CFI(I2,1)
             x2=x(INPB2)
             y2=y(INPB2)
             z2=z(INPB2)
             dist=(x1-x2)**2 +(y1-y2)**2 + (z1-z2)**2
             if (dist<mindist) then
               mindist=dist
               mini=I2
             endif
         enddo
         NUM_CFI(I1,4)=NUM_CFI(mini,1)
      enddo
   enddo
endif

!.....DEFAULT
    LTOUCH=.FALSE.
!
!..... MAKE A SEARCH ONLY FOR FIRST GRID

    DO NBL1=1,NBLOCK
!
       INTERF1: DO I1=NUM_SCF(NBL1)+1,NUM_SCF(NBL1)+NUM_CF(NBL1)
          IF(.NOT.LTOUCH(I1)) THEN
             LTOUCH(I1)=.TRUE.
             INPB1  =NUM_CFI(I1,1)
             INPB1T2=NUM_CFI(I1,4)
             ISIB1=NUM_CFI(I1,2)
!
             NBL2=NUM_CFI(I1,5)
             INTERF2: DO I2=NUM_SCF(NBL2)+1,NUM_SCF(NBL2)+NUM_CF(NBL2)
                IF(.NOT.LTOUCH(I2).AND.NUM_CFI(I2,5)==NBL1) THEN
                   INPB2  =NUM_CFI(I2,1)
                   INPB2T1=NUM_CFI(I2,4)
                   ISIB2=NUM_CFI(I2,2)
                   IF(INPB1==INPB2T1.AND.INPB2==INPB1T2) THEN
                      LTOUCH(I2)=.TRUE.
                      NUM_CFI(I1,6)=ISIB2
                      NUM_CFI(I2,6)=ISIB1
                      CYCLE INTERF1
                   END IF
!
                END IF
             END DO INTERF2 ! LOOP-2  * BLOCK INTERFACES
             WRITE(*,*) ' I CAN NOT FIND CONECTION, WILL STOP IN FINSID'
             WRITE(*,*) '  BLOCK=',NBL1,' P1 =',INPB1,' FACE =',ISIB1
             WRITE(*,*) '  BLOCK=',NBL2,' P2 =',INPB1T2
             CALL END
          END IF
       END DO INTERF1 ! LOOP-1   * BLOCK INTERFACES
    END DO  ! LOOP-1   * BLOCKS
!
  END SUBROUTINE FINSID
  
!*****************************************
SUBROUTINE OPFILE
!*****************************************
  IMPLICIT NONE
  LOGICAL, SAVE :: FILES_INITIALIZED = .FALSE.

  IUGRD= 12
  IUCHK= 13

  IF (.NOT. FILES_INITIALIZED) THEN

     OPEN(UNIT=IUGRD, FILE=FILGRD, STATUS='OLD', FORM='FORMATTED')
     REWIND(IUGRD)

     OPEN(UNIT=IUCHK, FILE='output.chk', STATUS='UNKNOWN', FORM='FORMATTED')
     REWIND(IUCHK)

     IURESX = 75
     IURESY = 76
     IURESZ = 77

     IF (LSTORE) THEN
        OPEN(UNIT= 78, FILE='field.plt'  , STATUS='REPLACE', FORM='FORMATTED', ACTION='WRITE')
        OPEN(UNIT=101, FILE='monitor.plt', STATUS='REPLACE', FORM='FORMATTED', ACTION='WRITE')
        OPEN(UNIT=102, FILE='uvwmax.plt' , STATUS='REPLACE', FORM='FORMATTED', ACTION='WRITE')
        OPEN(UNIT=103, FILE='RESULT.dat' , STATUS='REPLACE', FORM='FORMATTED', ACTION='WRITE')

        OPEN(UNIT=IURESX, FILE='Xplane.plt', STATUS='REPLACE', FORM='FORMATTED', ACTION='WRITE')
        OPEN(UNIT=IURESY, FILE='Yplane.plt', STATUS='REPLACE', FORM='FORMATTED', ACTION='WRITE')
        OPEN(UNIT=IURESZ, FILE='Zplane.plt', STATUS='REPLACE', FORM='FORMATTED', ACTION='WRITE')

        OPEN(UNIT=85, FILE='Xline.plt', STATUS='REPLACE', FORM='FORMATTED', ACTION='WRITE')
        OPEN(UNIT=86, FILE='Yline.plt', STATUS='REPLACE', FORM='FORMATTED', ACTION='WRITE')
        OPEN(UNIT=87, FILE='Zline.plt', STATUS='REPLACE', FORM='FORMATTED', ACTION='WRITE')

        IF (LNUSSELT) THEN
           OPEN(UNIT=104, FILE='Nu_x_0.plt'  , STATUS='REPLACE', FORM='FORMATTED', ACTION='WRITE')
           OPEN(UNIT=105, FILE='Nu_x_0_5.plt', STATUS='REPLACE', FORM='FORMATTED', ACTION='WRITE')
           OPEN(UNIT=106, FILE='Nu_x_1.plt'  , STATUS='REPLACE', FORM='FORMATTED', ACTION='WRITE')
        END IF

        WRITE(78 ,'(A)') 'TITLE="Oracle3D field"'

        IF (LVORTEX) THEN
           OPEN(UNIT=109, FILE='vortex_field.plt', STATUS='REPLACE', &
                FORM='FORMATTED', ACTION='WRITE')
           WRITE(109,'(A)') 'TITLE="Oracle3D vortex field"'
           WRITE(109,'(A)') 'VARIABLES="X","Y","Z",'// &
                            '"OMEGAX","OMEGAY","OMEGAZ",'// &
                            '"QCRIT","LAMBDA2"'
        END IF
        WRITE(78 ,'(A)') 'VARIABLES="X","Y","Z","U","V","W","P"'
        WRITE(101,*) 'VARIABLES= "Time","U","V","W","P","T","Q","VP"'
        WRITE(102,*) 'VARIABLES= "Time","U","V","W","Norm"'
     END IF

     FILES_INITIALIZED = .TRUE.
  END IF

END SUBROUTINE OPFILE
!*********************************************************************
  SUBROUTINE SETDATA
!*********************************************************************
    IMPLICIT NONE
    INTEGER(IK) :: I
!
    IUIN = 10                      ! UNIT FOR     INPUT  FILE  donn.dat

    OPEN(UNIT=IUIN , FILE='donn_400.dat')
    REWIND (IUIN)

!.....READ INPUT DATA FROM IUIN
!
    READ(IUIN,5) TITLE
  5 FORMAT(A50)
    READ(IUIN,6) FILGRD
  6 FORMAT(A)
    READ(IUIN,*) LSCREEN,LREAD,LWRITE,LTEST,LRESAVE
    READ(IUIN,*) LPRI,LPRJ,LPRK
!.....SELECTION OF DEPENDENT VARIABLES
    READ(IUIN,*) (LCAL(I),I=1,NPHI)
    READ(IUIN,*) LONGEO
!.....ITERATION AND CONV. LIMITS, OUTPUT CONTROL, FLUID PROPERTIES

    READ(IUIN,*) NBLMON,IMON,JMON,KMON
    READ(IUIN,*) NBLPR,IPR,JPR,KPR
    READ(IUIN,*) ALFA
!
!....TO CHOOSE THE STOPPING CRITERIA FOR SIMPLE LOOP: 1 - RESIDUAL 2 - RM2
!....FOR I_SC = 1.0, SORMAX = 1.E-3 OR 1.E-4
!....FOR I_SC = 2.0, SORMAX = 1.E-6 OR 1.E-7
!....SORMAX2 = (0.1 OR 0.01)*SORMAX1  --> TO JUMP OUT OF SIMPLE LOOP
!....IF RESIDUAL DOES NOT CONTINUE TO DECREASE
!....SLARGE = 1.E+3 OR 1.E+4 IS OK
!
    READ(IUIN,*) I_SC, SORMAX, SORMAX2, SLARGE, NSWEEP
!
    READ(IUIN,*) LTIME,LSTORE
!
    READ(IUIN,*) NTM,IMA,DT
    READ(IUIN,*) ITIME_U, ITIME_T
    READ(IUIN,*) DENSIT,VISCOS

!.....INITIALISATION OF VARIABLES
    READ(IUIN,*) UIN,VIN,WIN,TIN
    READ(IUIN,*) (GDS(I),I=1,NPHI)
    READ(IUIN,*) (URF(I),I=1,NPHI)
    READ(IUIN,*) (SOR(I),I=1,NPHI)
    READ(IUIN,*) (NSW(I),I=1,NPHI)

!.....NEW ADDED ......LID-DRIVEN CAVITY FLOW
    READ(IUIN,*) LNONDIMENSION        ! LOCICAL PARAMETER FOR NON-DIMENSIONAL COMPUTATION
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    READ(IUIN,*) RE_INPUT                  ! THE REYNOLDS NUMBER
    READ(IUIN,*) LUVWMAX              ! LOCICAL PARAMETER FOR FINDING THE MAXIMUM VELOCITY
!.....NEW ADDED ......NATURAL CONVECTION FLOW
    READ(IUIN,*) PRANL, RA            ! THE PRANDLT NUMBER AND RAYLEIGH NUMBER
    READ(IUIN,*) LNUSSELT             ! TO OUTPUT NUSSELT NUMBER INFORMATION OR NOT
!.....NEW ADDED ......EHD CONVECTION
    READ(IUIN,*) LEHD, LHYDROSTATIC
    READ(IUIN,*) (ESP_DFC(I),I=NPHI+1,NPHI2)
    READ(IUIN,*) (ITSWEEP(I),I=NPHI+1,NPHI2)
 !.....NEW ADDED .....EHD COMPUTATION CONTROL
    READ(IUIN,*) ITIME_Q
    READ(IUIN,*) (URF(I),I=NPHI+1,NPHI2)
    READ(IUIN,*) (SOR(I),I=NPHI+1,NPHI2)
    READ(IUIN,*) (NSW(I),I=NPHI+1,NPHI2)
!.....NEW ADDED .....ETHD CONVECTION
    READ(IUIN,*) LETHD, LESP
    READ(IUIN,*) L_PERMIT, N_MOBILITY  ! TEMPERATURE-DEPENDENT PERMITTIVITY AND IONIC MOBILITY
!....NEW ADDED ... For SUZEN-HUANG Model 
    READ(IUIN,*) D_C  !DImensionless coefficient of the SUZEN-Huang Model in the momentum equation for Coulomb force 
    READ(IUIN,*) F_BASE, VPmax ! FRequency of the unstediness of force (applied voltage - AC) and amplitude of VP
    READ(IUIN,*) LAMBDA ! THe DEBYE LENGTH 
    READ(IUIN,*) MuLOCATION , sigma ! lcoation parameter and decay rate parameter GAussian charge density profile
    READ(IUIN,*) Qmax  ! maximum charge density on surface
!....NEW ADDED ... Vortex identification criteria (GORTLER)
    READ(IUIN,*) LVORTEX  ! T=calcul vorticite+Q-critere+Lambda2, F=rien
    CLOSE(IUIN)
!
!==================================================================

  END SUBROUTINE SETDATA

!*****************************************
  SUBROUTINE MODCON
!*****************************************
!---------------------------------------------------------------
!     DEFINE CONSTANTS TO AVOID POSSIBLE PROBLEMS
!     WITH SOME COMPILERS
!---------------------------------------------------------------
    IMPLICIT NONE

    R_0_0=0.0
    SMALL = TINY(0._RK)
    GREAT = HUGE(0._RK)
!      SLARGE = 1.E20

  END SUBROUTINE MODCON

!*****************************************
  SUBROUTINE DEFMEM(NICV,NJCV,NKCV)
!*****************************************
    IMPLICIT NONE
    INTEGER(IK) :: NBL
    INTEGER(IK) :: NICV(NBLOCK),NKCV(NBLOCK),NJCV(NBLOCK)

!
!..... START AND END INDICES FOR ALL BLOCKS
    ICSTALL=NBL_ST(1)+1
    ICENALL=NBL_ST(NBLOCK)+ &
            NIBL(NBLOCK)*NJBL(NBLOCK)*NKBL(NBLOCK)
!
    NCVALL=0
    DO NBL=1,NBLOCK
      NCVALL=NICV(NBL)*NJCV(NBL)*NKCV(NBL)
    END DO

    CALL ALLOCDATA

  END SUBROUTINE DEFMEM

!*********************************************************************
  SUBROUTINE INIT
!*********************************************************************
    USE GEOM
    IMPLICIT NONE
    INTEGER(IK) :: NBL

!.....COMPUTE SOME GEOMETRIC CONSTANTS
    CALL  INIT_GEOM_T(X,Y,Z,XC,YC,ZC, &
            ARX_T,ARY_T,ARZ_T,AKX_T,AKY_T,AKZ_T,ARE2_T,ARKSI2_T, &
            ARE_T,ARER_T,DELN_T,DELNR_T)

!-----
    IF(LONGEO) THEN
       WRITE(*,6000) '********************************************'
       WRITE(*,6000) '  => LOGICAL VARIABLE => LONGEO <= IS .TRUE.'
       WRITE(*,6000) '          SAVE THE GEOMETRY DATA AND'
       WRITE(*,6000) '               STOP THE CODE'
       WRITE(*,6000) '********************************************'
  6000 FORMAT(35X,A)

       CALL WRGEO

       CALL END
    END IF
!
    DTR=1./DT
!
    NBLMON=MIN(NBLMON,NBLOCK)
!
    CALL SETIND(NBLMON)
    IMON=MIN(IMON,NI-1)
    JMON=MIN(JMON,NJ-1)
    KMON=MIN(KMON,NK-1)
    IJKMO=LK(KMON)+LI(IMON)+JMON
!
    NBLPR=MIN(NBLPR,NBLOCK)
!
    CALL SETIND(NBLPR)
    IPR=MIN(IPR,NI-1)
    JPR=MIN(JPR,NJ-1)
    KPR=MIN(KPR,NK-1)
    IJKPE=LK(KPR)+LI(IPR)+JPR


!.....PRANDTL NUMBER FOR FLUID

!.....PRANL FOR WATER, PRANL=7.0 (AT 20 DEGREES CELSIUS)
!.....PRANL FOR AIR  , PRANL=0.71

!      PRANT=1.  !! ALREAD READ FROM INPUT FILE
!
!.....RECIPROCAL VALUES OF PRANDTL NUMBERS
!
    PRTINV=1.

    PRTINV(IEN)=1./PRANL

!.....RECIPROCAL VALUES OF UNDERRRELAXATION FACTORS
    URFR=1./URF
    URFM=1.-URF

!...some other reciprocals 
    sigmaR=1./sigma
	LAMBDAR=1./LAMBDA
!
    CALL INIT_USER

!.....DENSITY VISCOSITY
    VIS=VISCOS
    DEN=DENSIT
!.....PERMITTIVITY MOBILITY
    ESP_P=1.
    ESP_K=1.
!
!.....SOME OTHER CONSTANTS
!
    ITIM=0
    TIME=0.

    GAMT_U = 0. ! <---TIME SCHEME (VELOCITY)      (IE AND I3L STARTING)
    GAMT_T = 0. ! <---TIME SCHEME (TEMPERATURE)   (IE AND I3L STARTING)
    GAMT_Q = 0. ! <---TIME SCHEME (CHARGE DENSITY)(IE AND I3L STARTING)

!.....MONITORING POINT, PRESSURE REFERENCE POINT
    IJKMON=IJKMO
    IJKPR=IJKPE

    CALL CHVVEC(U,V,W)
    CALL CHVSCA(VIS)
    CALL CHVSCA(DEN)
!
!==================================================
! PERTURBATION INITIALE POUR INSTABILITES DE GORTLER
!==================================================
  IF (LVORTEX) CALL INIT_GORTLER_PERTURBATION
!

  END SUBROUTINE INIT

!**************************************************
  SUBROUTINE READGRID
!**************************************************
    IMPLICIT NONE
    REAL(RK)   ,ALLOCATABLE :: X1(:),Y1(:),Z1(:)
    INTEGER(IK),ALLOCATABLE :: LI1(:),LK1(:),NUM_IND1(:,:)
    INTEGER(IK),ALLOCATABLE :: NICV(:), NJCV(:), NKCV(:),NIJKBL(:)
    INTEGER(IK),ALLOCATABLE :: NUM_BF(:)
    INTEGER(IK)             :: NUM_BL,I,J,K,INP,IDR,idr1(1)
    INTEGER(IK)             :: ISIDE,IDEW,NBL
!RETURN

    REWIND(IUGRD)
!.....
    READ(IUGRD,*)  NBLOCK
    ALLOCATE(NUM_SDR(NBLOCK),NUM_SCF(NBLOCK),NBL_ST(NBLOCK))
    ALLOCATE(NIBL(NBLOCK),NJBL(NBLOCK),NKBL(NBLOCK),NIJKBL(NBLOCK))
    ALLOCATE(NICV(NBLOCK),NJCV(NBLOCK),NKCV(NBLOCK))
    ALLOCATE(NUM_CF(NBLOCK),NUM_DR(NBLOCK),NUM_BF(NBLOCK))

    READ(IUGRD,*)  NUM_CF_ALL
    ALLOCATE(NUM_CFI(NUM_CF_ALL,8) , &
             NUM_NUM(NUM_DR_ALL), &
             NUM_TYP(NUM_DR_ALL), &
             NUM_PR( NUM_DR_ALL), &
             NUM_SPR(NUM_DR_ALL))
NUM_CFI=0
NUM_CF_ALL=0


!.....READ GRID

    NUM_SDR(1)=0
    NUM_SPR=-1
    NUM_SCF(1)=0
    NBL_ST=0
    NXYZA=0
    NXMAX=0
    NZMAX=0
    NUM_BF_ALL=0

!  INITIALISE IN ORDER TO AVOID UNECESSARY TESTS
    NI=1
    NJ=1
    NK=1
    NIM=1
    NJM=1
    NKM=1
    NIJ=1
    NIK=1
    NJK=1
    NIJK=1
!.
    WRITE(*,'(1(/))')
    WRITE(*,640)
!
    DO NBL=1,NBLOCK
!
!.....READ NO OF BLOCK
       READ(IUGRD,*) NUM_BL

       IF(NBL/=NUM_BL) THEN
          WRITE(*,*) ' STOP BY READING THE GRID NO ', NUM_BL, &
                     ' BLOCK NUMBER  MUST BE ORDERED, IT SHOULD BE THE NUMBER ',NBL
          CALL END
       END IF

       READ(IUGRD,*) NICV(NBL),NJCV(NBL),NKCV(NBL)
       READ(IUGRD,*) NUM_CF(NBL)
       READ(IUGRD,*) NUM_DR(NBL)

!
!..... CALCULATE NUMBER OF CV
!
       NUM_BF(NBL)=2*(NICV(NBL)*NJCV(NBL) &
                  +NICV(NBL)*NKCV(NBL)+NJCV(NBL)*NKCV(NBL))
!
       NIBL(NBL)=NICV(NBL)+2
       NJBL(NBL)=NJCV(NBL)+2
       NKBL(NBL)=NKCV(NBL)+2
!
       NIJKBL(NBL)=NIBL(NBL)*NJBL(NBL)*NKBL(NBL)
!
         WRITE(*,650)  NBL,NICV(NBL),NJCV(NBL),NKCV(NBL)

       IF(NBL>1) THEN
          ALLOCATE( X1(NXYZA+NIJKBL(NBL)), &
                    Y1(NXYZA+NIJKBL(NBL)), &
                    Z1(NXYZA+NIJKBL(NBL)))
          X1(1:NXYZA)=X ; Y1(1:NXYZA)=Y ; Z1(1:NXYZA)=Z
          DEALLOCATE(X,Y,Z)
       END IF
       NXYZA=NXYZA+NIJKBL(NBL)
         ALLOCATE( X(NXYZA), Y(NXYZA), Z(NXYZA))
       IF(NBL>1) THEN
          X=X1 ; Y=Y1 ; Z=Z1
          DEALLOCATE(X1,Y1,Z1)
       END IF

       IF (NIBL(NBL)>NXMAX.OR.NBL==1) THEN
          IF(NBL>1) THEN
             ALLOCATE( LI1(NIBL(NBL)))
             LI1(1:NXMAX)=LI
             DEALLOCATE(LI)
          END IF
          NXMAX=NIBL(NBL)
         ALLOCATE( LI(NXMAX))
          IF(NBL>1) THEN
             LI=LI1
             DEALLOCATE(LI1)
          END IF
       END IF

       IF (NKBL(NBL)>NZMAX.OR.NBL==1) THEN
          IF(NBL>1) THEN
             ALLOCATE( LK1(NKBL(NBL)))
             LK1(1:NZMAX)=LK
             DEALLOCATE(LK)
          END IF
          NZMAX=NKBL(NBL)
         ALLOCATE( LK(NZMAX))
          IF(NBL>1) THEN
             LK=LK1
             DEALLOCATE(LK1)
          END IF
       END IF

       IF(NBL>1) THEN
          NBL_ST (NBL)= NBL_ST(NBL-1)+NIJKBL(NBL-1)
       END IF

       CALL SETIND(NBL)

       DO K=1,NKM
          DO I=1,NIM
             DO J=1,NJM
                INP=LK(K)+LI(I)+J
                READ(IUGRD,*) X(INP),Y(INP),Z(INP)
             END DO
          END DO
       END DO
!
       IF(NBL>1) NUM_SDR(NBL)=NUM_SDR(NBL-1)+NUM_DR(NBL-1)
       IF(NBL>1) NUM_SCF(NBL)=NUM_SCF(NBL-1)+NUM_CF(NBL-1)
       NUM_CF_ALL=NUM_CF_ALL+NUM_CF(NBL)
!

       IF(NBL>1) THEN
          ALLOCATE( NUM_IND1(NUM_BF_ALL+ NUM_BF(NBL),2))
          NUM_IND1(1:NUM_BF_ALL,:)=NUM_IND
          DEALLOCATE(NUM_IND)
       END IF
       NUM_BF_ALL=NUM_BF_ALL + NUM_BF(NBL)
       ALLOCATE( NUM_IND(NUM_BF_ALL,2))
       IF(NBL>1) THEN
          NUM_IND=NUM_IND1
          DEALLOCATE(NUM_IND1)
       END IF


       DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)
!
          READ(IUGRD,*) NUM_NUM(IDR),NUM_TYP(IDR),NUM_PR(IDR)
          idr1=maxloc(NUM_SPR)
          if (NUM_SPR(idr1(1))/=-1) then
             NUM_SPR(IDR)=NUM_SPR(idr1(1))+NUM_PR(idr1(1))
          else
            NUM_SPR(IDR)=0
          endif
          DO I=NUM_SPR(IDR)+1,NUM_SPR(IDR)+NUM_PR(IDR)
             READ(IUGRD,*) (NUM_IND(I,K), K=1,2)
          END DO
       END DO
!
       DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
          READ(IUGRD,*) NUM_CFI(I,1),NUM_CFI(I,2),NUM_CFI(I,4),NUM_CFI(I,5),NUM_CFI(I,7)
         
       END DO

       DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
          ISIDE=NUM_CFI(I,2)
          IDEW=ID_SIDE  (ID_ODER(ISIDE,1))*ID_SIGN(ISIDE)
          NUM_CFI(I,3)=NUM_CFI(I,1)+IDEW
       END DO

    END DO   ! END BLOCKS
    
     CLOSE(IUGRD)


!
    NDREG=NUM_SDR(NBLOCK)+NUM_DR(NBLOCK)
!
    WRITE(*,640)

    CALL DEFMEM(NICV,NJCV,NKCV)
    CALL FINSID

!*****************************
! WRITE GRID TO TECPLOT FORMAT
! ***************************

    OPEN (UNIT=79, FILE='grid_output.plt',STATUS='UNKNOWN')
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
       WRITE(79,*)'VARIABLES= "X","Y","Z"'
           WRITE(79,*)'ZONE I=',NIM,', J=', NJM,', K=',NKM
       DO K=1,NKM
          DO J=1,NJM
             DO I=1,NIM
                INP=LK(K)+LI(I) +J
                WRITE(79,*) X(INP),Y(INP),Z(INP)
             END DO
          END DO
       END DO
    END DO
    close(79)

  640 FORMAT(30X,53('#'))
  650 FORMAT(30X,                              &
         'READ BLOCK NBL = ',I3,' NICV = ',I3, &
         ' NJCV = ',I3,' NKCV = ',I3)

    DEALLOCATE(NICV,NJCV,NKCV,NIJKBL,NUM_BF)
  END SUBROUTINE READGRID

!**************************************************
  SUBROUTINE WRGEO
!**************************************************
    IMPLICIT NONE
    INTEGER(IK) :: IUTMP,NBL,I,J,K,INP,IREG
    INTEGER(IK) :: I1,I2,I3,I4,I5,I6,I7,I8
    CHARACTER(LEN=6) :: FVERT,FCELL

    IUTMP=51
    WRITE(*,*)
!=============================
!.....WRITE CELL VERTICES    =
!=============================
!.....
    WRITE(FVERT,'(A5,I1)') 'VERT.',0
    OPEN(IUTMP,FILE=FVERT)
    REWIND IUTMP
!.....
    WRITE(*,606) FVERT
!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!
       DO K=1,NKM
          DO I=1,NIM
             DO J=1,NJM
!.....TAKE TOP, EAST, NORTH CORNER
                INP=LK(K)+LI(I)+J
!
                WRITE(IUTMP,610) INP,X(INP),Y(INP),Z(INP)
             END DO
          END DO
       END DO
!
    END DO  !  END ALL BLOCKS
!
!.....CLOSE FILE
    CLOSE(IUTMP)
!.....
!
!=============================
!.....WRITE CELLS            =
!=============================
!
    WRITE(FCELL,'(A5,I1)') 'CELL.',0
    OPEN(IUTMP,FILE=FCELL)
    REWIND IUTMP
  !.....
    WRITE(*,606) FCELL
!
    IREG=1
!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!
       DO K=2,NKM
          DO I=2,NIM
             DO J=2,NJM
!.....CELL DEFINED WITH EIGHT VERTICES
                INP=LK(K)+LI(I)+J
!
                I1=INP-1-NJ-NIJ
                I2=I1+NJ
                I3=I2+1
                I4=I3-NJ
                I5=I1+NIJ
                I6=I5+NJ
                I7=I6+1
                I8=I7-NJ
                WRITE(IUTMP,620) INP,I1,I2,I3,I4,I5,I6,I7,I8,IREG
             END DO
          END DO
       END DO
    END DO  !  END ALL BLOCKS
!
!.....CLOSE FILE
    CLOSE(IUTMP)
!
!.....
  606 FORMAT(28X,'WRITE GEOMETRICAL DATA IN FILE = > ',A)
  610 FORMAT(I7,6X,3G16.9)
  620 FORMAT(I7,6X,9I7)


  END SUBROUTINE WRGEO

!**************************************************
  SUBROUTINE WRBFAC(ITYPCH,FFILD)
!**************************************************
    IMPLICIT NONE

    INTEGER(IK)  ,INTENT(INOUT) :: ITYPCH
    CHARACTER(*),INTENT(INOUT) :: FFILD

    INTEGER(IK)       :: IUTMP,IBCF,NBL,IDR,NBCTYP
    INTEGER(IK)       :: NBC_S,NBC_E,I,INP,ISIDE,IDEW,IDNS,IDTB
    INTEGER(IK)       :: INPG,INS,INB,INBS
    CHARACTER(LEN=6)  :: FTOWR
!
    IUTMP=51
!
    WRITE(FTOWR,'(A5,I1)') FFILD,0
!
!.....OPEN FILE
    OPEN(IUTMP,FILE=FTOWR)
    REWIND IUTMP
!
!.... BOUNDARY FACES
!
    IBCF=0
!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!.....ALL BC
       DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)
!.....TYP
          NBCTYP=NUM_TYP(IDR)
          IF(NBCTYP==ITYPCH) THEN
!.....
             NBC_S=NUM_SPR(IDR)+1
             NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)
             DO I=NBC_S,NBC_E
!
                INP=NBL_ST(NBL)+NUM_IND(I,1)
                ISIDE=NUM_IND(I,2)
!
                IDEW=ID_SIDE  (ID_ODER(ISIDE,1))*ID_SIGN(ISIDE)
                IDNS=ID_SIDE  (ID_ODER(ISIDE,2))
                IDTB=ID_SIDE  (ID_ODER(ISIDE,3))
!
                INPG=INP+MIN(IDEW,0_ik)
                INS=INPG-IDNS
                INB=INPG-IDTB
                INBS=INB-IDNS
!
                IBCF=IBCF+1
                WRITE(IUTMP,620) IBCF,INPG,INS,INB,INBS,NBCTYP,0
!
             END DO
          END IF
       END DO
    END DO  !  END ALL BLOCKS
!
!.....CLOSE FILE
    CLOSE(IUTMP)
  620 FORMAT(I7,6X,9I7)

  END SUBROUTINE WRBFAC

!**************************************************
  SUBROUTINE WRUNSTRUCT
!**************************************************
    IMPLICIT NONE
    INTEGER(IK) :: NNODES,NBL,I,J,K,NELEM,INP
    INTEGER(IK) :: I1,I2,I3,I4,I5,I6,I7,I8

!=============================
!.....WRITE CELL VERTICES    =
!=============================
!.....

! CALCULATE THE TOTAL NUMBER OF NODES
! ************************************

    NNODES = 0
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!
       DO K=1,NKM
          DO I=1,NIM
             DO J=1,NJM

                NNODES = NNODES +1
!
             END DO
          END DO
       END DO
!
    END DO  !  END ALL BLOCKS

! CALCULATE THE TOTAL NUMBER OF CELLS
! ************************************

    NELEM = 0
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!
       DO K=2,NKM
          DO I=2,NIM
             DO J=2,NJM

                NELEM = NELEM + 1

             END DO
          END DO
       END DO
    END DO  !  END ALL BLOCKS

!**********************************************

!.....
    WRITE(99,*)'VARIABLES= "X","Y","Z"'
    WRITE(99,*)'ZONE N=',NNODES,', E=', NELEM, &
               ', ZONE DATAPACKING=POINT',     &
               ', ZONETYPE=FEBRICK'
!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!
       DO K=1,NKM
          DO I=1,NIM
             DO J=1,NJM

                INP=LK(K)+LI(I)+J
!
                WRITE(99,*) X(INP),Y(INP),Z(INP)
             END DO
          END DO
       END DO
!
    END DO  !  END ALL BLOCKS
!
!=============================
!.....WRITE CELLS            =
!=============================

!.....ALL BLOCKS
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
!
       DO K=2,NKM
          DO I=2,NIM
             DO J=2,NJM
!.....CELL DEFINED WITH EIGHT VERTICES
                INP=LK(K)+LI(I)+J
!
                I1=INP-1-NJ-NIJ
                I2=I1+NJ
                I3=I2+1
                I4=I3-NJ
                I5=I1+NIJ
                I6=I5+NJ
                I7=I6+1
                I8=I7-NJ
                WRITE(99,*) I1,I2,I3,I4,I5,I6,I7,I8
             END DO
          END DO
       END DO
    END DO  !  END ALL BLOCKS

  END SUBROUTINE WRUNSTRUCT

!**************************************************
  SUBROUTINE WREST
!**************************************************
    IMPLICIT NONE
    INTEGER(IK) :: IUTMP,INPST,INPEN,INP,I
!
    WRITE(*,600)
!
!.....WRITE FIELD VALUES FOR THE NEXT RUN
!
    IUTMP=51
    WRITE(FILRES,'(A5,I1)') 'REST.',0
    OPEN (IUTMP,FILE=FILRES,FORM='UNFORMATTED')
    REWIND IUTMP
    WRITE(*,610) FILRES
    WRITE(IUTMP) TIME, AMPLITUDE

    INPST=ICSTALL
    INPEN=ICENALL
    WRITE(IUTMP) INPST,INPEN,ITIM,TIME,(F1(INP),INP=INPST,INPEN), &
             (F2(INP),INP=INPST,INPEN),(F3(INP),INP=INPST,INPEN), &
             ( U(INP),INP=INPST,INPEN),( V(INP),INP=INPST,INPEN), &
             ( W(INP),INP=INPST,INPEN),( P(INP),INP=INPST,INPEN)
!-------------------------------
!
!.....INTERFACE MASS FLUXES
!
    WRITE(IUTMP) (FLX_CFI(I),I=    &
                     NUM_SCF(1)+1, &
                     NUM_SCF(NBLOCK)+NUM_CF(NBLOCK))
!-------------------------------
!
    IF(LTIME) WRITE(IUTMP) (UO(INP),INP=INPST,INPEN), &
                           (VO(INP),INP=INPST,INPEN), &
                           (WO(INP),INP=INPST,INPEN)

    IF(LCAL(IEN)) THEN
       WRITE(IUTMP) (T(INP),INP=INPST,INPEN)
       IF(LTIME) WRITE(IUTMP) (TO(INP),INP=INPST,INPEN)
    END IF
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Modif

    IF((LEHD.OR.LETHD).AND.ABS(D_C)>SMALL) THEN
       WRITE(IUTMP) (Q(INP),INP=INPST,INPEN)
       WRITE(IUTMP) (QO(INP),INP=INPST,INPEN)
       WRITE(IUTMP) (VP(INP),INP=INPST,INPEN)
       WRITE(IUTMP) (EX(INP),INP=INPST,INPEN), &
                    (EY(INP),INP=INPST,INPEN), &
                    (EZ(INP),INP=INPST,INPEN)
       WRITE(IUTMP) (ESP_P(INP),INP=INPST,INPEN)
       WRITE(IUTMP) (ESP_K(INP),INP=INPST,INPEN)
    END IF

!
    REWIND IUTMP
      CLOSE(IUTMP)
!
  600 FORMAT(/,25X,'PREPARE RESULTS OF CALCULATION FOR RESTART', &
             /,25X,'******************************************')
  610 FORMAT(25X,' WRITE RESULTS IN FILE :',A)
  END SUBROUTINE WREST

!**************************************************
  SUBROUTINE RREST
!**************************************************
    IMPLICIT NONE
    INTEGER(IK) :: IUTMP,INPST,INPEN,INP,I
!
!
!
    IUTMP=51
    WRITE(FILRES,'(A5,I1)') 'REST.',0
    OPEN (IUTMP,FILE=FILRES,FORM='UNFORMATTED',STATUS='OLD')
    REWIND IUTMP
    WRITE(*,610) FILRES
!
    READ(IUTMP) TIME, AMPLITUDE
!
    READ(IUTMP) INPST,INPEN, ITIM,TIME,(F1(INP),INP=INPST,INPEN), &
              (F2(INP),INP=INPST,INPEN),(F3(INP),INP=INPST,INPEN), &
              ( U(INP),INP=INPST,INPEN), (V(INP),INP=INPST,INPEN), &
              ( W(INP),INP=INPST,INPEN), (P(INP),INP=INPST,INPEN)
!-------------------------------
!
!.....INTERFACE MASS FLUXES
!
    READ(IUTMP) (FLX_CFI(I),I=   &
                   NUM_SCF(1)+1, &
                   NUM_SCF(NBLOCK)+NUM_CF(NBLOCK))
!-------------------------------
!
    IF(LTIME) READ(IUTMP) (UO(INP),INP=INPST,INPEN), &
                          (VO(INP),INP=INPST,INPEN), &
                          (WO(INP),INP=INPST,INPEN)

    IF(LCAL(IEN)) THEN
       READ(IUTMP) (T(INP),INP=INPST,INPEN)
       IF(LTIME) READ(IUTMP) (TO(INP),INP=INPST,INPEN)
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MODIF
    IF((LEHD.OR.LETHD).AND.ABS(D_C)>SMALL) THEN
       READ(IUTMP) (Q(INP),INP=INPST,INPEN)
       READ(IUTMP) (QO(INP),INP=INPST,INPEN)
       READ(IUTMP) (VP(INP),INP=INPST,INPEN)
       READ(IUTMP) (EX(INP),INP=INPST,INPEN), &
                   (EY(INP),INP=INPST,INPEN), &
                   (EZ(INP),INP=INPST,INPEN)
       READ(IUTMP) (ESP_P(INP),INP=INPST,INPEN)
       READ(IUTMP) (ESP_K(INP),INP=INPST,INPEN)
    END IF
!
    REWIND IUTMP
    CLOSE(IUTMP)

!.....TIME COUNTER  => ONE LESS
    ITIM=ITIM-1
!
    PP=P

    IF (ITIME_U==1) THEN
       GAMT_U = 1.0 ! <--- VELOCITY TIME SCHEME (SWITH FOR I3L)
    END IF
    IF (ITIME_T==1) THEN
       GAMT_T = 1.0 ! <--- TEMPERATURE TIME SCHEME (SWITH FOR I3L)
    END IF
    IF(ITIME_Q==1)  THEN
       GAMT_Q = 1.0 ! <--- CHARGE DENSITY TIME SCHEME (SWITH FOR I3L)
    END IF
!
  600 FORMAT(/,25X,'READ RESULTS OF CALCULATION FROM RESTART FILE', &
             /,25X,'*********************************************')
  610 FORMAT(25X,' RESTART FILE :',A)
  END SUBROUTINE RREST

!**************************************************
  SUBROUTINE SAVCIN
!**************************************************
    IMPLICIT NONE
!
    INTEGER(IK) :: ITMP(4)
    INTEGER(IK) :: NBL,I,INP,ISIDE,IDEW,IDNS,IDTB

    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
     DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
        INP=NBL_ST(NBL)+NUM_CFI(I,1)
        ISIDE=NUM_CFI(I,2)
!
        IDEW=ID_SIDE  (ID_ODER(ISIDE,1))*ID_SIGN(ISIDE)
        IDNS=ID_SIDE  (ID_ODER(ISIDE,2))
        IDTB=ID_SIDE  (ID_ODER(ISIDE,3))
!
        ITMP(1) =INP+MIN(IDEW,0_ik)
        ITMP(2) =ITMP(1)-IDNS
        ITMP(3) =ITMP(2)-IDTB
        ITMP(4) =ITMP(1)-IDTB
!..... AND NOW SAVE IT
        X_CFI(I,:)=X(ITMP)
        Y_CFI(I,:)=Y(ITMP)
        Z_CFI(I,:)=Z(ITMP)
!
     END DO
    END DO  ! END BLOCKS
!
  END SUBROUTINE SAVCIN

!*********************************************************************
  SUBROUTINE MOOUT1
!*********************************************************************

    IMPLICIT NONE
!
    WRITE(*,661) TITLE,FLOWIN,DENSIT, &
         VISCOS,ALFA,SORMAX,SMALL,GREAT
    IF(LCAL(IU))   WRITE(*,663)' U ',URF(IU),  &
                               ' U ',GDS(IU), NSW(IU)
    IF(LCAL(IV))   WRITE(*,663)' V ',URF(IV),  &
                               ' V ',GDS(IU), NSW(IV)
    IF(LCAL(IW))   WRITE(*,663)' W ',URF(IW),  &
                               ' W ',GDS(IU), NSW(IW)
    IF(LCAL(IP))   WRITE(*,664)' P ',URF(IP),  &
                                     NSW(IP)
    IF(LCAL(IEN))  WRITE(*,663)' T ',URF(IEN), &
                               ' T ',GDS(IEN),NSW(IEN)
    IF(LEHD.OR.LETHD) THEN
       WRITE(*,599) 'EHD.OR.ETHD CONVECTION:  Q --> SOLVED'
       WRITE(*,599) 'EHD.OR.ETHD CONVECTION: VP --> SOLVED'
    END IF

    IF(LTIME)  WRITE(*,669) DT

  661 FORMAT(//,30X,A50,/,30X,50('*'),//,40X,         &
         'FLOW INLET       : FLOW = ',1PE12.4,/,40X,  &
         'FLUID DENSITY    :  RHO = ',1PE12.4,/,40X,  &
         'DYNAMIC VISCOSITY:   MU = ',1PE12.4,/,40X,  &
         'ALFA  PARAMETER  : ALFA = ',0PF12.2,/,40X, &
         'CONV. CRITERION  :  SOR = ',1PE12.4,/,40X,  &
         '                  SMALL = ',1PE12.4,/,40X,  &
         '                  GREAT = ',1PE12.4,/,40X)
  599 FORMAT(40X,A40)
  663 FORMAT(40X,'URF(',A3,')=',F5.2,'  GDS(',A3,')=',F5.2,'  NSW=',I5)
  664 FORMAT(40X,'URF(',A3,')=',F5.2,15X,'   NSW=',I5)
  665 FORMAT(40X,'URF(',A3,')=',F5.2)
  669 FORMAT(40X,50('*')/,40X,'UNSTEADY CALCULATION   ', &
            /40X,             'TIME STEP DT= ',1PE13.3)
    WRITE(*,*)
    WRITE(*,*)

  END SUBROUTINE MOOUT1

!**************************************************
  SUBROUTINE CHVSCA(FI,BUF2,N,alias)
!**************************************************
! EXCHANGE AT INTERFACES BETWEEN PROCS
    IMPLICIT NONE
    REAL(RK),INTENT(INOUT)        :: FI(NXYZA)
    INTEGER(IK),INTENT(IN),OPTIONAL          :: N,alias
    REAL(RK),OPTIONAL,INTENT(INOUT) :: BUF2(:,:,:)

    INTEGER(IK) :: NBL2,INE,J,NBL1,NADD,ISIDE,IND
		INTEGER(IK) :: id_pass

	if(PRESENT(alias)) THEN
		id_pass=alias
	else 
		id_pass=0
	end if 

    DO NBL1=1,NBLOCK
      DO NBL2=1,NBLOCK
        IF(LINK_TAB(NBL1,NBL2)/=0) THEN   ! A MESSAGE HAVE TO BE TRANSMITTED
          J=0                             ! FROM NBL1 TO NBL2
            ! WALKING ON ALL THE INTERFACES OF THE RECEIVER
            ! IN ORDER TO CONSERVE THE ORDER
            DO I=NUM_SCF(NBL2)+1,NUM_SCF(NBL2)+NUM_CF(NBL2)
             
       IF(NUM_CFI(I,7)==0.OR.id_pass==7) THEN
         
              IF(NUM_CFI(I,5)==NBL1) THEN            ! I HAVE TO SEND THIS POINT
                INE=NBL_ST(NBL1) + NUM_CFI(I,4)      ! INDEX OF THE POINT

                IF(PRESENT(N)) THEN                  ! I HAVE TO GO N POINTS INSIDE
                  ISIDE=NUM_CFI(I,2)
                  IF(MOD(ISIDE,2_ik)==0) THEN !!! ISIDE = 2,4,6

                    IF(ISIDE==2) THEN
                       NADD = NJBL(NUM_CFI(I, 5))
                    ELSE IF(ISIDE==4) THEN
                       NADD = 1
                    ELSE IF(ISIDE==6) THEN
                       NADD = NIBL(NUM_CFI(I, 5))*NJBL(NUM_CFI(I, 5))
                    END IF

                    INE = INE + (N-1)*NADD
                  ELSE                     !!! ISIDE = 1,3,5

                    IF(ISIDE==1) THEN
                       NADD = NJBL(NUM_CFI(I, 5))
                    ELSE IF(ISIDE==3) THEN
                       NADD = 1
                    ELSE IF(ISIDE==5) THEN
                       NADD = NIBL(NUM_CFI(I, 5))*NJBL(NUM_CFI(I, 5))
                    END IF

                    INE = INE -(N-1)*NADD
                  END IF
                END IF

                J=J+1 ; BIG_BUF(J,NBL1,NBL2,1)=FI(INE) ! PUT IT IN A BUFFER
              END IF
end if
            END DO
            BIG_BUF(1:J,NBL1,NBL2,2)=BIG_BUF(1:J,NBL1,NBL2,1)
        END IF
      END DO
    END DO


    DO NBL1=1,NBLOCK
      DO NBL2=1,NBLOCK
        IF(LINK_TAB(NBL1,NBL2)/=0) THEN        ! A MESSAGE IS BEING TRANSMITTED
            IF(.NOT.PRESENT(BUF2)) THEN        ! PUT THE VALUE BACK IN THE INPUT
              J=0
              ! WALKING ON ALL THE INTERFACES OF THE RECEIVER
              ! IN ORDER TO CONSERVE THE ORDER
              DO I=NUM_SCF(NBL2)+1,NUM_SCF(NBL2)+NUM_CF(NBL2)

IF(NUM_CFI(I,7)==0.OR.id_pass==7) THEN  

                IF (NUM_CFI(I,5)==NBL1) THEN
                  IND=NBL_ST(NBL2) + NUM_CFI(I,3)      ! INDEX OF THE POINT
                  J=J+1 ; FI(IND)=BIG_BUF(J,NBL1,NBL2,2) ! GET IT FROM THE BUFFER
                END IF   
end if
              END DO
            END IF
        END IF
      END DO
    END DO


  END SUBROUTINE CHVSCA

!**************************************************
  SUBROUTINE CHVVEC(FIX,FIY,FIZ)
!**************************************************
! EXCHANGE INFO ON INTERFACES FOR A VECTOR
    IMPLICIT NONE
    REAL(RK),INTENT(INOUT) :: FIX(NXYZA),FIY(NXYZA),FIZ(NXYZA)

    CALL CHVSCA(FIX)
    CALL CHVSCA(FIY)
    CALL CHVSCA(FIZ)

  END SUBROUTINE CHVVEC

!*********************************************************************
  SUBROUTINE FIND_UVWMAX
!*********************************************************************
    IMPLICIT NONE
    REAL(RK)    :: UVW_NORM
    INTEGER(IK) :: NBL,I,J,K,INP

    UMAX = -GREAT
    VMAX = -GREAT
    WMAX = -GREAT
    UVW_NORM_MAX = -GREAT

    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
       DO  K = 2, NKM
          DO  I = 2, NIM
             DO  J = 2, NJM

               INP = LK(K) + LI(I) + J
!               IF(ABS(ZC(INP)-0.5).LE.ABS(ZC(INP+NIJ)-ZC(INP-NIJ))/4.) THEN

!         KPLAN = 1 + NKM/2               !! FOR COMPARSION WITH BENCHMARK
!         K = KPLAN                       !! SOLUTIONS OF NATURAL CONVECTION

                UVW_NORM = SQRT(V(INP)**2+U(INP)**2+W(INP)**2)

                IF (UVW_NORM>UVW_NORM_MAX) THEN
                   UVW_NORM_MAX = UVW_NORM
                   IMAX_NORM = I
                   JMAX_NORM = J
                   KMAX_NORM = K
                   XMAX_NORM = XC(INP)
                   YMAX_NORM = YC(INP)
                   ZMAX_NORM = ZC(INP)
                END IF
                IF(U(INP)>UMAX) THEN
                   UMAX = U(INP)
                   IMAX_U = I
                   JMAX_U = J
                   KMAX_U = K
                   XMAX_U = XC(INP)
                   YMAX_U = YC(INP)
                   ZMAX_U = ZC(INP)
                END IF
                IF(V(INP)>VMAX) THEN
                   VMAX = V(INP)
                   IMAX_V = I
                   JMAX_V = J
                   KMAX_V = K
                   XMAX_V = XC(INP)
                   YMAX_V = YC(INP)
                   ZMAX_V = ZC(INP)
                END IF
                IF(W(INP)>WMAX) THEN
                   WMAX = W(INP)
                   IMAX_W = I
                   JMAX_W = J
                   KMAX_W = K
                   XMAX_W = XC(INP)
                   YMAX_W = YC(INP)
                   ZMAX_W = ZC(INP)
                END IF
!                END IF
             END DO
          END DO
       END DO
!.....BLOCK
    END DO

  END SUBROUTINE FIND_UVWMAX

!*********************************************************************
  SUBROUTINE NUSSELT_COMPUTATION
!*********************************************************************
    IMPLICIT NONE
    INTEGER(IK) :: NBL
!      DIMENSION FIF(NXYZA)

    R_NU_0  = 0.0
    R_NU_05 = 0.0
    R_NU_10 = 0.0

    DO NBL=1,NBLOCK
       CALL SETIND(NBL)

       WRITE(104,2) 'ZONE ,J=',NJM-1,',K=',NKM-1,',F=POINT' ! X = 0.0, COLD
       WRITE(105,2) 'ZONE ,J=',NJM-1,',K=',NKM-1,',F=POINT' ! X = 0.5, MIDPLANE
       ! 50X50X50 CVS
       WRITE(106,2) 'ZONE ,J=',NJM-1,',K=',NKM-1,',F=POINT' ! X = 1.0, HOT

       CALL NU_X_MIDDLE
       CALL NU_X_WALL(NBL)

    END DO   ! END BLOCK

  2 FORMAT(A12,I4,A9,I4,A8)

  END SUBROUTINE NUSSELT_COMPUTATION

!*********************************************************************
  SUBROUTINE NU_X_MIDDLE
!*********************************************************************
    IMPLICIT NONE
    INTEGER(IK)            :: I,J,K,INP,INE
    REAL(RK)               :: TEMP_2
!    I = NIM/2                                  !! X= 0.5

    DO K=2,NKM-1
      DO I=2,NIM-1
         DO J=2,NJM-1
          INP=LK(K)+LI(I) +J
          INE=INP+NJ
          IF(ABS(XC(INE)-0.5)+ABS(XC(INP)-0.5).LE.ABS(XC(INE)-XC(INP))) THEN
            TEMP_2 = U(INP)*T(INP)-((T(INE) - T(INP))/(XC(INE)-XC(INP)))
            R_NU_05 = R_NU_05 + TEMP_2*ARE_T(INP,2)
            WRITE(105,*)  YC(INE),' ', ZC(INE), ' ', TEMP_2
          ENDIF
       END DO
    END DO
    END DO

  END SUBROUTINE NU_X_MIDDLE


!*********************************************************************
  SUBROUTINE NU_X_WALL(NBL)
!*********************************************************************
    USE GEOM
    IMPLICIT NONE
    INTEGER(IK),INTENT(IN) :: NBL

    REAL(RK)    :: TEMP_1,TEMP_3
    INTEGER(IK) :: IDR,NBCTYP,NBC_S,NBC_E,I,INP,ISIDE

!.....ALL BC
    DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)
!.....TYP
       NBCTYP=NUM_TYP(IDR)
!.....
       NBC_S=NUM_SPR(IDR)+1
       NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)

       IF(NBCTYP==8) THEN  !!! ISOTHERMAL SIDES _ X = 1.0

          DO I=NBC_S,NBC_E
             INP=NBL_ST(NBL)+NUM_IND(I,1)
             ISIDE =  NUM_IND(I,2)
             CALL  INIT_GEOM_BF(ISIDE,INP)
             TEMP_1 = (T(INP)-T(INE))/(XC(INP)-XC(INE))
             R_NU_0 = R_NU_0 + TEMP_1*ARE
             WRITE(104,*)  YC(INE),' ', ZC(INE), ' ', TEMP_1
          END DO

       ELSE IF(NBCTYP==9) THEN   !!! ISOTHERMAL SIDES _ X = 0.0

          DO I=NBC_S,NBC_E
             INP=NBL_ST(NBL)+NUM_IND(I,1)
             ISIDE =  NUM_IND(I,2)
             CALL  INIT_GEOM_BF(ISIDE,INP)
             TEMP_3 = (T(INP)-T(INE))/(XC(INP)-XC(INE))
             R_NU_10 = R_NU_10 + TEMP_3*ARE
             WRITE(106,*)  YC(INE),' ', ZC(INE), ' ', TEMP_3
          END DO

       ELSE
!
!......NOTHING HERE
!
         WRITE(*,*)  'PROBLEM, NUSSELT NUBERS: BC TYPE!! '
       END IF

    END DO ! SCAN BF
  END SUBROUTINE NU_X_WALL


!*********************************************************************
  SUBROUTINE DIVERGENCE(PHY)
!*********************************************************************
    IMPLICIT NONE
    REAL(RK)   ,INTENT(INOUT) :: PHY(NXYZA)

    INTEGER(IK) :: NBL,I,J,K,INP

    RM2 = 0.
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
       DO K=2,NKM
          DO I=2,NIM
             DO J=2,NJM
                INP=LK(K)+LI(I)+J
                RM2(NBL)= RM2(NBL) + (AP(INP)*PHY(INP)                          &
                                    + AE(INP)*PHY(INP+NJ )+AW(INP)*PHY(INP-NJ)  &
                                    + AN(INP)*PHY(INP+1  )+AS(INP)*PHY(INP-1)   &
                                    + AT(INP)*PHY(INP+NIJ)+AB(INP)*PHY(INP-NIJ) )**2
             END DO
             END DO
          END DO
       END DO

  RMM2 = SUM(SQRT(RM2))
  END SUBROUTINE DIVERGENCE

!*********************************************************************
  SUBROUTINE NORM_STAT(NORM_UVW)
!*********************************************************************
    IMPLICIT NONE
    REAL(RK)   ,INTENT(INOUT) :: NORM_UVW

    INTEGER(IK) :: NBL,I,J,K,INP

    RM2 = 0.
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
       DO K=2,NKM
          DO I=2,NIM
             DO J=2,NJM
                INP=LK(K)+LI(I)+J
                RM2(NBL)= RM2(NBL) + (u(inp)-uo(inp))**2 &
                                     + (v(inp)-vo(inp))**2 &
                                     + (w(inp)-wo(inp))**2
             END DO
             END DO
          END DO
       END DO

  NORM_UVW = SQRT(SUM(RM2))
  END SUBROUTINE NORM_STAT

!**************************************************
  SUBROUTINE ES_MODULE
!**************************************************
    IMPLICIT NONE
     
    
      WRITE (*,'(A)',ADVANCE='YES') '------> ELECTRIC POTENTIAL '
    CALL SRC_TERM(IVP)
    CALL CALPOISSON(IVP, VP)                          !!! ELECTRICAL POTENTIAL  VP(N+1)
!
    WRITE (*,'(A)',ADVANCE='YES') '------> ELECTRIC FIELD     '
    CALL ELECTRIC_FIELD(VP)                           !!! ELECTRIC FIELD  E(N+1)

!
    WRITE (*,'(A)',ADVANCE='YES') '------> CHARGE DENSITY     '
    CALL SRC_TERM(IQ)
    CALL SUZEN_HUANG(IQ, Q)  !!! CHARGE DENSITY Q(N+1)

    WRITE (*,'(A)',ADVANCE='YES') '------> BODY FORCE      '

      ! FXsh=Q*EX
      ! FYsh=Q*EY
      ! FZsh=Q*EZ
      
    ! warning attention au passage des arguments PI, F_BASE, TIME
    ! FXsh=D_C*Q*EX*VOL*((SIN(2.*PI*F_BASE*TIME))**2)
    !  FYsh=D_C*Q*EY*VOL*((SIN(2.*PI*F_BASE*TIME))**2)
    !  FZsh=D_C*Q*EZ*VOL*((SIN(2.*PI*F_BASE*TIME))**2)
 
      ! Fmag=SQRT(FXsh*FXsh+FYsh*FYsh+FZsh*FZsh)
         Fmag=SQRT(FXsh**2+FYsh**2+FZsh**2)
           
 END SUBROUTINE ES_MODULE


!**************************************************
   SUBROUTINE FORCE_SH 
!**************************************************
    IMPLICIT NONE
   
   CALL ELECTRIC_FIELD(VP)
   CALL SUZEN_HUANG(IQ,Q)

      FXsh=Q*EX
      FYsh=Q*EY
      FZsh=Q*EZ

   END SUBROUTINE FORCE_SH


!**************************************************
  SUBROUTINE ELECTRIC_FIELD(PHI)
!**************************************************
    IMPLICIT NONE
    REAL(RK),INTENT(INOUT) :: PHI(NXYZA)

    CALL GRADFI(GR1X,GR1Y,GR1Z,PHI,IVP)
    CALL BPRES(GR1X,IVP)
    CALL BPRES(GR1Y,IVP)
    CALL BPRES(GR1Z,IVP)
!
!=============================================================C
!......NEED A SUBROUTINE TO CONTROL THE ELECTRIC FIELD ON     C
!......ELECTRODES CONSIDERING CONSTANT ELECTRICAL POTENTIAL?  C
!=============================================================C

!      CALL BC_ELECTRIC(GRX_TEMP,GRY_TEMP,GRZ_TEMP)

    EX=-GR1X
    EY=-GR1Y
    EZ=-GR1Z

 !  E = SQRT(EX*EX+EY*EY+EZ*EZ)

  END SUBROUTINE ELECTRIC_FIELD

!**************************************************
  SUBROUTINE SRC_TERM(IFI)
!**************************************************
    IMPLICIT NONE
    INTEGER(IK),INTENT(INOUT) :: IFI

    IF(IFI==6) THEN      !! CHARGE DENSITY EQUATION
       SUO= 0.
    ELSE IF(IFI==7) THEN  !! ELECTRIC POTENTIAL EQUATION
       SUO= 0.
    ELSE
       WRITE(*,*) 'SRC_TERM, PROBLEM !'
    END IF
  END SUBROUTINE SRC_TERM

!**************************************************
  SUBROUTINE SUZEN_HUANG(IFI, FI)
!**************************************************
    IMPLICIT NONE

    REAL(RK),ALLOCATABLE   ,INTENT(INOUT) :: FI(:)
    INTEGER(IK),INTENT(INOUT) :: IFI

    REAL(RK),ALLOCATABLE      :: FI_TEMP(:)
    REAL(RK)                  :: URFRS,URFMS
    INTEGER(IK)               :: NBL,IT_DFC
    LOGICAL                   :: PAS_FINI
!    INTEGER(IK)               :: INP
    URFRS=URFR(IFI)
    URFMS=URFM(IFI)

    ALLOCATE( FI_TEMP(NXYZA))
!======================================================C
!                                                      C
!  NOTES:(1)DEFERRED CORRECTION TECHNIQUE IS REQUIRED  C
!           WHEN THE MESH IS NON-ORTHOGONAL            C
!        (2)SDC SCHEME IS USED CURRENTLY; TO EXTEND TO C
!           IDC SCHEME : IDC IS DIFFICULT TO IMPLEMENT C
!                                                      C
!======================================================C

    
    IT_DFC = 0            !!! DFC ITERATION TIMES
    PAS_FINI=.TRUE.
    FI_TEMP=0.
    DO WHILE(PAS_FINI)    !!! DFC LOOP START

       PAS_FINI=.FALSE.
       IT_DFC=IT_DFC+1
!
!..CLEAR AND INILIZATION ...
       SU=SUO
       SP=0.
      
       CALL GRADFI(GR1X,GR1Y,GR1Z,FI)
    
      DO NBL=1,NBLOCK
          CALL SETIND(NBL)
!
!.....EAST CELL - FACE
          CALL CELSH(NI-2,NJM,NKM,NJ,1,NIJ,FX,AE,AW,2)
!.....NORTH CELL - FACE
          CALL CELSH(NIM,NJ-2,NKM,1,NIJ,NJ,FY,AN,AS,1)
!.....TOP   CELL - FACE
          CALL CELSH(NIM,NJM,NK-2,NIJ,NJ,1,FZ,AT,AB,3)


!
!==================================
!.....IMPLEMENT BOUNDARY CONDITIONS
!==================================
!
          CALL MODSH(NBL,FI)


    END DO     !   END ALL BLOCKS

       AP_ADD= VOL*((LAMBDAR**2))  

       AP=-AE-AW-AN-AS-AT-AB+SP+AP_ADD
       AP=AP*URFRS
       SU=SU+URFMS*AP*FI
       
            
      
!
!=====================================
!     SOLVE LINEAR SYSTEM FOR FI
!======================================
!

       CALL SIPSOL(FI,IFI)

      CALL TEST(FI,FI_TEMP,PAS_FINI,IT_DFC,IFI)
       FI_TEMP=FI

    END DO !! DFC LOOP END

    DEALLOCATE( FI_TEMP)

 
  END SUBROUTINE SUZEN_HUANG

!******************************************************************
  SUBROUTINE CELSH(NIE,NJE,NKE,IDEW,IDNS,IDTB, &
                        FIF,ACFE,ACFW,DIR1)
!*******************************************************************
 !   USE DEC
    IMPLICIT NONE

    REAL(RK)   ,INTENT(IN)    :: FIF(NXYZA)
    REAL(RK)   ,INTENT(INOUT) :: ACFE(NXYZA),ACFW(NXYZA)
    INTEGER(IK),INTENT(IN)    :: NIE,NJE,NKE,IDEW,IDNS,IDTB,DIR1

    REAL(RK)    :: DE,SUEH,SUEL,SUADD,GAME
    INTEGER(IK) :: I,J,K
    REAL(RK)    :: FXE,FXW
    INTEGER(IK) :: INP,INE



     CALL ESP_SH(NXYZA, ESP_P,NI,NJ,NK)
        
!.....CALCULATE EAST,TOP,NORTH  CELL FACE
!
    DO K=2,NKE
       DO I=2,NIE
          DO J=2,NJE
!
             INP=LK(K)+LI(I)+J
             INE=INP+IDEW
!.....
!.....INTERPOLATION FACTOR
             FXE=FIF(INP)
             FXW=1.-FXE

!!.....DIFUSION COEFFICIENT
! for the harmonic mean of the permittivity 
 ! GAME=(ESP_P(INP)*ESP_P(INE))/(ESP_P(INP)*FXE+ESP_P(INE)*FXW)

! Arithmatc mean of permittivity 
!
             GAME=ESP_P(INP)*FXW+ESP_P(INE)*FXE
             DE=GAME*SQRT(ARE2_T(INP,DIR1)/ARKSI2_T(INP,DIR1)) !!! SDC SCHEME --> IDC ????;
!
             ACFE(INP)=-DE            !!! THE SIGN CHANGED HERE COMPARE TO POISSON
             ACFW(INE)=-DE            !!! THE SIGN CHANGED HERE  COMPARE TO POISSON

              SUEH=GAME*((GR1X(INP)*FXW+GR1X(INE)*FXE)*ARX_T(INP,DIR1)+ &    !!! NO DFC _ ORTHOGONAL GRID
                         (GR1Y(INP)*FXW+GR1Y(INE)*FXE)*ARY_T(INP,DIR1)+ &
                         (GR1Z(INP)*FXW+GR1Z(INE)*FXE)*ARZ_T(INP,DIR1))
!.....
              SUEL= DE*((GR1X(INP)*FXW+GR1X(INE)*FXE)*AKX_T(INP,DIR1)+ &      !! THE SIGN CHANGED HERE COMPARE TO POISSON
                       (GR1Y(INP)*FXW+GR1Y(INE)*FXE)*AKY_T(INP,DIR1)+ &
                       (GR1Z(INP)*FXW+GR1Z(INE)*FXE)*AKZ_T(INP,DIR1))

             SUADD=+(SUEH-SUEL)                                !!! CHANGE THE SIGN

             SU(INP) = SU(INP)+SUADD
             SU(INE) = SU(INE)-SUADD

          END DO
       END DO
    END DO
!
  END SUBROUTINE CELSH

!**************************************************
  SUBROUTINE CALPOISSON(IFI, FI)
!**************************************************
    IMPLICIT NONE

    REAL(RK),ALLOCATABLE   ,INTENT(INOUT) :: FI(:)
    INTEGER(IK),INTENT(INOUT) :: IFI

    REAL(RK),ALLOCATABLE      :: FI_TEMP(:)
    REAL(RK)                  :: URFRS,URFMS
    INTEGER(IK)               :: NBL,IT_DFC
    LOGICAL                   :: PAS_FINI
    INTEGER(IK)               :: I,J,K,INP  !,NI,NJ,NK
  
    URFRS=URFR(IFI)
    URFMS=URFM(IFI)

    ALLOCATE( FI_TEMP(NXYZA))
!======================================================C
!                                                      C
!  NOTES:(1)DEFERRED CORRECTION TECHNIQUE IS REQUIRED  C
!           WHEN THE MESH IS NON-ORTHOGONAL            C
!        (2)SDC SCHEME IS USED CURRENTLY; TO EXTEND TO C
!           IDC SCHEME : IDC IS DIFFICULT TO IMPLEMENT C
!                                                      C
!======================================================C

    IT_DFC = 0            !!! DFC ITERATION TIMES
    PAS_FINI=.TRUE.
    FI_TEMP=0.
    DO WHILE(PAS_FINI)    !!! DFC LOOP START

       PAS_FINI=.FALSE.
       IT_DFC=IT_DFC+1
!
!..CLEAR AND INILIZATION ...
       SU=SUO
       SP=0.

       CALL GRADFI(GR1X,GR1Y,GR1Z,FI,IVP)

 DO NBL=1,NBLOCK
          CALL SETIND(NBL)


!.....EAST CELL - FACE
          CALL CELPOISSON(NI-2,NJM,NKM,NJ,1,NIJ,FX,AE,AW,2)
!.....NORTH CELL - FACE
          CALL CELPOISSON(NIM,NJ-2,NKM,1,NIJ,NJ,FY,AN,AS,1)
!.....TOP   CELL - FACE
          CALL CELPOISSON(NIM,NJM,NK-2,NIJ,NJ,1,FZ,AT,AB,3)
!
!==================================
!.....IMPLEMENT BOUNDARY CONDITIONS
!==================================
!
         CALL MODPOISSON(NBL,FI)

       END DO  !   END ALL BLOCKS

      
       AP=-AE-AW-AN-AS-AT-AB+SP
       AP=AP*URFRS
       SU=SU+URFMS*AP*FI     

!=====================================
!     SOLVE LINEAR SYSTEM FOR FI
!======================================
!
       CALL SIPSOL(FI,IFI)
       CALL TEST(FI,FI_TEMP,PAS_FINI,IT_DFC,IFI)
       FI_TEMP=FI
   
  
   END DO !! DFC LOOP END
    DEALLOCATE( FI_TEMP)
  END SUBROUTINE CALPOISSON

!*********************************************************************
  SUBROUTINE BPRES_POISSON(FI)
!*********************************************************************
    USE GEOM
    IMPLICIT NONE
    REAL(RK),INTENT(INOUT) :: FI(NXYZA)

    INTEGER(IK)            :: NBL,IDR,NBCTYP,NBC_S,NBC_E,I,INP,ISIDE

    DO NBL=1,NBLOCK
       CALL SETIND(NBL)

       DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)
          NBCTYP=NUM_TYP(IDR)
          NBC_S=NUM_SPR(IDR)+1
          NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)
          IF((NBCTYP==8).OR.(NBCTYP==9)) THEN    !! DO NOTHING WITH UP AND BOTTOM SIDE....
             DO I=NBC_S, NBC_E                      !!
             END DO

!
!.......DIRICHLET BOUNDARY CONDITONS.....NO EXTRAPOLATION....
!
!
!
!.......NEUMANN BOUNDARY CONDITONS
!
          ELSE

             DO I = NBC_S, NBC_E

                INP = NBL_ST(NBL) + NUM_IND(I,1)
                ISIDE = NUM_IND(I,2)

                CALL INIT_GEOM_BF(ISIDE,INP)

                FI(INE) = FI(INP)    !!! 1ST ORDER EXTRAPOLATION....OK ?

             END DO

          END IF

       END DO    !   END ALL FACES
    END DO     !   END ALL BLOCKS
  END SUBROUTINE BPRES_POISSON


!******************************************************************
  SUBROUTINE CELPOISSON(NIE,NJE,NKE,IDEW,IDNS,IDTB, &
                        FIF,ACFE,ACFW,DIR1)
!*******************************************************************
    IMPLICIT NONE

    REAL(RK)   ,INTENT(IN)    :: FIF(NXYZA)
    REAL(RK)   ,INTENT(INOUT) :: ACFE(NXYZA),ACFW(NXYZA)
    INTEGER(IK),INTENT(IN)    :: NIE,NJE,NKE,IDEW,IDNS,IDTB,DIR1

    REAL(RK)    :: DE,SUEH,SUEL,SUADD,GAME
    INTEGER(IK) :: I,J,K
    REAL(RK)    :: FXE,FXW
    INTEGER(IK) :: INP,INE

     CALL ESP_SH(NXYZA, ESP_P,NI,NJ,NK)  

!.....CALCULATE EAST,TOP,NORTH  CELL FACE
!
    DO K=2,NKE
       DO I=2,NIE
          DO J=2,NJE
!
             INP=LK(K)+LI(I)+J
             INE=INP+IDEW
!.....
!.....INTERPOLATION FACTOR
             FXE=FIF(INP)
             FXW=1.-FXE
!
!.....DIFUSION COEFFICIENT
! for the harmonic mean of the permittivity 
 ! GAME=(ESP_P(INP)*ESP_P(INE))/(ESP_P(INP)*FXE+ESP_P(INE)*FXW)

! Arithmatc mean of permittivity 
             GAME=ESP_P(INP)*FXW+ESP_P(INE)*FXE

             DE=GAME*SQRT(ARE2_T(INP,DIR1)/ARKSI2_T(INP,DIR1)) !!! SDC SCHEME --> IDC ????;
!
             ACFE(INP)=DE                                      !!! THE SIGN CHANGED HERE
             ACFW(INE)=DE

              SUEH=GAME*((GR1X(INP)*FXW+GR1X(INE)*FXE)*ARX_T(INP,DIR1)+ &    !!! NO DFC _ ORTHOGONAL GRID
                         (GR1Y(INP)*FXW+GR1Y(INE)*FXE)*ARY_T(INP,DIR1)+ &
                         (GR1Z(INP)*FXW+GR1Z(INE)*FXE)*ARZ_T(INP,DIR1))
!.....
              SUEL=DE*((GR1X(INP)*FXW+GR1X(INE)*FXE)*AKX_T(INP,DIR1)+ &
                       (GR1Y(INP)*FXW+GR1Y(INE)*FXE)*AKY_T(INP,DIR1)+ &
                       (GR1Z(INP)*FXW+GR1Z(INE)*FXE)*AKZ_T(INP,DIR1))

             SUADD=-(SUEH-SUEL)                                !!! CHANGE THE SIGN

             SU(INP) = SU(INP)+SUADD
             SU(INE) = SU(INE)-SUADD

          END DO
       END DO
    END DO
!
  END SUBROUTINE CELPOISSON

!*********************************************
  SUBROUTINE RHS_POISSON(NBL,PST,SU_DFC)
!*********************************************
    IMPLICIT NONE
    REAL(RK)   ,INTENT(IN) :: PST(NXYZA),SU_DFC(NXYZA)
    INTEGER(IK),INTENT(IN) :: NBL

    INTEGER(IK) :: I,J,K,INP
!
    DO K=2,NKM
       DO I=2,NIM
          DO J=2,NJM
             INP=LK(K)+LI(I)+J
             SU(INP)=PST(INP)+SU_DFC(INP)
          END DO
       END DO
    END DO
  END SUBROUTINE RHS_POISSON

!********************************************************************
  SUBROUTINE MODPOISSON(NBL,FI)
!********************************************************************
!     THE BOUNDARY CONDITION FOR POISSION EQUATION IS TREATED HERE  C
!     ONLY TWO TYPES OF BOUNDARY CONDITIONS ARE CONSIDERED:         C
!         1. DIRICHLET                                              C
!         2. NEUMANN                                                C
!===================================================================C
    USE GEOM
    IMPLICIT NONE
    REAL(RK)   ,INTENT(INOUT) :: FI(NXYZA)
    INTEGER(IK),INTENT(IN) :: NBL

    REAL(RK)    :: DE,PK_X,PK_Y,PK_Z
    INTEGER(IK) :: INP,IDR,NBCTYP,NBC_S,NBC_E,I,ISIDE
!
!
!.....BLOCK INTERFACES
!
    CALL COMPOISSON(NBL,FI)

    DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)
!.....TYP

       NBCTYP=NUM_TYP(IDR)
!.....
       NBC_S=NUM_SPR(IDR)+1
       NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)

       IF((NBCTYP==10).OR.(NBCTYP==11)) THEN
!
!.......DIRICHLET BOUNDARY CONDITONS
!
          DO I=NBC_S, NBC_E
             INP=NBL_ST(NBL)+NUM_IND(I,1)
             ISIDE=NUM_IND(I,2)
             CALL INIT_GEOM_BF(ISIDE,INP)

             DE=ESP_P(INE)*SQRT(ARE2/ARKSI2)
             SUP = -DE
!
!  SDC: SU <-- PART_1 + PART_2 (DFC)
!

             SUU = -DE*FI(INE)                             &
                  -((ARX-SQRT(ARE2/ARKSI2)*AKX)*GR1X(INP)+ &
                    (ARY-SQRT(ARE2/ARKSI2)*AKY)*GR1Y(INP)+ &
                    (ARZ-SQRT(ARE2/ARKSI2)*AKZ)*GR1Z(INP))

             ! ?????
             !             SUU = -DE*FI(INE)                                 &
             !                   -DE*((SQRT(ARE2/ARKSI2)*ARX-AKX)*GR1X(INP)+ &
             !                        (SQRT(ARE2/ARKSI2)*ARY-AKY)*GR1Y(INP)+ &
             !                        (SQRT(ARE2/ARKSI2)*ARZ-AKZ)*GR1Z(INP)) ?????

!=======================
! IDEA OF VIRTUAL POINT    !! TO IMPROVE? !!!
!=======================

             SP(INP)=SP(INP)+SUP
             SU(INP)=SU(INP)+SUU
          END DO

       ELSEIF((NBCTYP==4).OR.(NBCTYP==3)) THEN
!
!.......NEUMANN BOUNDARY CONDITONS (NECESSARY SINCE NO UPDATE ROUTINE)
!
          DO I=NBC_S, NBC_E
             INP=NBL_ST(NBL)+NUM_IND(I,1)
             ISIDE=NUM_IND(I,2)
             CALL INIT_GEOM_BF(ISIDE,INP)
!.....CHOICE 1: SIMPLE BUT INACCURATE
             FI(INE)=FI(INP)
!.....CHOICE 2: COMPLICATED BUT ACCURATE

!=======================
! IDEA OF VIRTUAL POINT    !! MARK !!!
!=======================

             PK_X = AKX - (AKX*ARX+AKY*ARY+AKZ*ARZ)*ARX/ARE2
             PK_Y = AKY - (AKX*ARX+AKY*ARY+AKZ*ARZ)*ARY/ARE2
             PK_Z = AKZ - (AKX*ARX+AKY*ARY+AKZ*ARZ)*ARZ/ARE2
             FI(INE)=FI(INP)+GR1X(INP)*PK_X+GR1Y(INP)*PK_Y+GR1Z(INP)*PK_Z

          END DO

       END IF

    END DO  !! END ALL BC
  END SUBROUTINE MODPOISSON

!********************************************************************
SUBROUTINE MODSH(NBL,FI)
!********************************************************************
!     THE BOUNDARY CONDITION FOR POISSON EQUATION IS TREATED HERE
!     ONLY TWO TYPES OF BOUNDARY CONDITIONS ARE CONSIDERED:
!         1. DIRICHLET
!         2. NEUMANN
!====================================================================

  USE GEOM
  IMPLICIT NONE

  INTEGER(IK),INTENT(IN)    :: NBL
  REAL(RK),INTENT(INOUT)    :: FI(NXYZA)

  REAL(RK)    :: DE,PK_X,PK_Y,PK_Z
  INTEGER(IK) :: INP,IDR,NBCTYP,NBC_S,NBC_E,I,ISIDE

!.....BLOCK INTERFACES
  CALL COMSH(NBL,FI)

!.....LOOP OVER ALL BOUNDARY REGIONS OF BLOCK NBL
  DO IDR = NUM_SDR(NBL)+1, NUM_SDR(NBL)+NUM_DR(NBL)

     NBCTYP = NUM_TYP(IDR)
     NBC_S  = NUM_SPR(IDR) + 1
     NBC_E  = NUM_SPR(IDR) + NUM_PR(IDR)

!.... TYPE 12 : DIRICHLET BOUNDARY CONDITION
     IF (NBCTYP == 12) THEN

        DO I = NBC_S, NBC_E
           INP   = NBL_ST(NBL) + NUM_IND(I,1)
           ISIDE = NUM_IND(I,2)

           CALL INIT_GEOM_BF(ISIDE,INP)

           DE  = ESP_P(INE)*SQRT(ARE2/ARKSI2)
           SUP = DE

!..... SDC: SU <-- PART_1 + PART_2 (DFC)
           SUU = DE*FI(INE) + &
                 ( (ARX - SQRT(ARE2/ARKSI2)*AKX)*GR1X(INP) + &
                   (ARY - SQRT(ARE2/ARKSI2)*AKY)*GR1Y(INP) + &
                   (ARZ - SQRT(ARE2/ARKSI2)*AKZ)*GR1Z(INP) )

           SP(INP) = SP(INP) + SUP
           SU(INP) = SU(INP) + SUU
        END DO

!.... TYPES 14, 10, 3, 4 : NEUMANN BOUNDARY CONDITION
     ELSEIF ( (NBCTYP == 14) .OR. (NBCTYP == 10) .OR. &
              (NBCTYP == 3 ) .OR. (NBCTYP == 4 ) ) THEN

        DO I = NBC_S, NBC_E
           INP   = NBL_ST(NBL) + NUM_IND(I,1)
           ISIDE = NUM_IND(I,2)

           CALL INIT_GEOM_BF(ISIDE,INP)

!..... CHOICE 1: SIMPLE BUT INACCURATE
!          FI(INE) = FI(INP)

!..... CHOICE 2: COMPLICATED BUT ACCURATE
!=======================
! IDEA OF VIRTUAL POINT
!=======================
           PK_X = AKX - (AKX*ARX + AKY*ARY + AKZ*ARZ)*ARX/ARE2
           PK_Y = AKY - (AKX*ARX + AKY*ARY + AKZ*ARZ)*ARY/ARE2
           PK_Z = AKZ - (AKX*ARX + AKY*ARY + AKZ*ARZ)*ARZ/ARE2

           FI(INE) = FI(INP) + GR1X(INP)*PK_X + &
                                GR1Y(INP)*PK_Y + &
                                GR1Z(INP)*PK_Z
        END DO

     END IF

  END DO

END SUBROUTINE MODSH
!******************************************************************
  SUBROUTINE COMPOISSON(NBL,FI)
!******************************************************************
    USE GEOM
    IMPLICIT NONE
    REAL(RK)   ,INTENT(IN)    :: FI(NXYZA)
    INTEGER(IK),INTENT(IN)    :: NBL

    REAL(RK)    :: GAME,DE,SUEH,SUEL,SUADD
    INTEGER(IK) :: I,INP,ISIDE,IND
!
    CALL SETIND(NBL)

  DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
      INP=NBL_ST(NBL)+NUM_CFI(I,1)
      IND=NBL_ST(NBL)+NUM_CFI(I,3)
      ISIDE=NUM_CFI(I,2)
     
      CALL INIT_GEOM_BLOCK(ISIDE,I)

      AKX=XC(IND)-XC(INP)
      AKY=YC(IND)-YC(INP)
      AKZ=ZC(IND)-ZC(INP)
      ARKSI2=AKX**2+AKY**2+AKZ**2

!
!.....INTERPOLATION FACTOR
      FXE=FX_CFI(I)
      FXW=1.-FXE
!
!      GAME = 1.0
      GAME=ESP_P(INP)*FXW+ESP_P(IND)*FXE
!
!.....DIFUSION COEFFICIENT
      DE=GAME*SQRT(ARE2/ARKSI2)
!
      A_E(I)= DE  !!! CHANGE THE SIGN HERE
!
!.....
      SUEH=GAME*((GR1X(INP)*FXW+GR1X(IND)*FXE)*ARX+ &
                 (GR1Y(INP)*FXW+GR1Y(IND)*FXE)*ARY+ &
                 (GR1Z(INP)*FXW+GR1Z(IND)*FXE)*ARZ)
!.....
      SUEL=DE*((GR1X(INP)*FXW+GR1X(IND)*FXE)*AKX+ &
               (GR1Y(INP)*FXW+GR1Y(IND)*FXE)*AKY+ &
               (GR1Z(INP)*FXW+GR1Z(IND)*FXE)*AKZ)
!
      SUADD=-(SUEH-SUEL) !!! CHANGE THE SIGN HERE
!
      SU(INP)=SU(INP)+SUADD
!
      SP (INP)= SP(INP)-A_E(I)
  
       END DO


  END SUBROUTINE COMPOISSON

!******************************************************************
  SUBROUTINE COMSH(NBL,FI)
!******************************************************************

    USE GEOM
    IMPLICIT NONE
    REAL(RK)   ,INTENT(IN)    :: FI(NXYZA)
    INTEGER(IK),INTENT(IN)    :: NBL

    REAL(RK)    :: GAME,DE,SUEH,SUEL,SUADD
    INTEGER(IK) :: I,INP,ISIDE,IND
!
    CALL SETIND(NBL)

    DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
      INP=NBL_ST(NBL)+NUM_CFI(I,1)
      IND=NBL_ST(NBL)+NUM_CFI(I,3)
      ISIDE=NUM_CFI(I,2)
     
!
	  	IF(NUM_CFI(I,7)==0) THEN 

      CALL INIT_GEOM_BLOCK(ISIDE,I)

      AKX=XC(IND)-XC(INP)
      AKY=YC(IND)-YC(INP)
      AKZ=ZC(IND)-ZC(INP)
      ARKSI2=AKX**2+AKY**2+AKZ**2

!
!.....INTERPOLATION FACTOR
      FXE=FX_CFI(I)
      FXW=1.-FXE
!
!      GAME = 1.0
      GAME=ESP_P(INP)*FXW+ESP_P(IND)*FXE
!
!.....DIFUSION COEFFICIENT
      DE=GAME*SQRT(ARE2/ARKSI2)
!
      A_E(I)=-DE  !!! CHANGED SIGN HERE for SH MODEL 
!
!.....
      SUEH=GAME*((GR1X(INP)*FXW+GR1X(IND)*FXE)*ARX+ &
                 (GR1Y(INP)*FXW+GR1Y(IND)*FXE)*ARY+ &
                 (GR1Z(INP)*FXW+GR1Z(IND)*FXE)*ARZ)
!.....
      SUEL=DE*((GR1X(INP)*FXW+GR1X(IND)*FXE)*AKX+ & !!! CHANGED SIGN HERE for SH MODEL 
               (GR1Y(INP)*FXW+GR1Y(IND)*FXE)*AKY+ &
               (GR1Z(INP)*FXW+GR1Z(IND)*FXE)*AKZ)
!
      SUADD=+(SUEH-SUEL) !!! CHANGE THE SIGN HERE
!
      SU(INP)=SU(INP)+SUADD
!
      SP (INP)= SP(INP)-A_E(I)
    	END IF

    END DO
  END SUBROUTINE COMSH

!*******************************************************************
  SUBROUTINE TEST(FI,FIO,PAS_FINI,IT_DFC,IFI)
!*******************************************************************
    IMPLICIT NONE

    LOGICAL,INTENT(INOUT)  :: PAS_FINI
    REAL(RK)   ,INTENT(IN) :: FI(NXYZA),FIO(NXYZA)
    INTEGER(IK),INTENT(IN) :: IT_DFC,IFI

    REAL(RK)    :: DIFF
    INTEGER(IK) :: I,J,K,INP,NBL

    DO NBL=1,NBLOCK
      CALL SETIND(NBL)
    DO K=2,NKM
       DO I=2,NIM
          DO J=2,NJM
             INP=LK(K)+LI(I)+J
             DIFF=ABS(FI(INP)-FIO(INP))

!	print*,'the DIFF (FI-FI0) till now is =', DIFF ,'IT_DFC is = ',IT_DFC
             IF(DIFF>ESP_DFC(IFI)) PAS_FINI=.TRUE.
          END DO
       END DO
    END DO
    END DO

    IF(IT_DFC>ITSWEEP(IFI)) THEN
       WRITE(*,*) IFI,'RUN OUT OF DEFERRED CORRECTION LOOP!!'
       PAS_FINI=.FALSE.
    END IF
  END SUBROUTINE TEST

!*******************************************************************
  SUBROUTINE UPDATE(FI,FIO)
!*******************************************************************
    IMPLICIT NONE
    REAL(RK)   ,INTENT(INOUT) :: FIO(NXYZA)
    REAL(RK)   ,INTENT(IN)    :: FI(NXYZA)

    INTEGER(IK) :: I,J,K,INP

!==================================C
! 2 FUNCTIONS FOR TRANS3D: [1,NI]  C
!==================================C

    DO K=1,NK
       DO I=1,NI
          DO J=1,NJ
             INP=LK(K)+LI(I)+J
             FIO(INP)=FI(INP)
          END DO
       END DO
    END DO
  END SUBROUTINE UPDATE

!*****************************************************************
  SUBROUTINE TRAN3D(IFI,FI,FIO,FIOO,DEN_TRAN,GAM_TRAN)
!*****************************************************************
    IMPLICIT NONE

    REAL(RK),ALLOCATABLE   ,INTENT(INOUT) :: FI(:)
    REAL(RK)   ,INTENT(INOUT) :: FIO(NXYZA),FIOO(NXYZA)
    REAL(RK)   ,INTENT(IN)    :: DEN_TRAN,GAM_TRAN

    REAL(RK),ALLOCATABLE    :: FI_TEMP(:),APT(:),AP_ADD(:)
    REAL(RK)    :: URFRS,URFMS
    INTEGER(IK) :: NBL,IFI,IT_DFC
    LOGICAL     :: PAS_FINI

    ALLOCATE(FI_TEMP(NXYZA),APT(NXYZA),AP_ADD(NXYZA))
!......TO COMPUTE THE EFFECTIVE VELOCITIES........
!......NON-DIMENSIONAL FORM: U + R*K*E

    U_=U+R*EX*ESP_K
    V_=V+R*EY*ESP_K
    W_=W+R*EZ*ESP_K

!.......TO COMPUTE THE EFFECTIVE FLUX..........

    CALL NEWFLX2

!.......TO DO THE SPATIAL DISCRETIZATION ......

    URFRS=URFR(IFI)
    URFMS=URFM(IFI)

    FI_TEMP=FI

    PAS_FINI=.TRUE.
    IT_DFC = 0

!###############################################
! IMPLICIT SCHEMES REQUIRE DEFERRED CORRECTION
!###############################################

    DO WHILE(PAS_FINI)        !!! DFC LOOP START

       PAS_FINI=.FALSE.
       IT_DFC=IT_DFC+1
       CALL GRADFI(GR1X,GR1Y,GR1Z,FI)

       SU=SUO
       SP=0.

       CALL GHOST_CELLS(FI,FI_TEMP)

!=============================================
!     CALCULATE FLUXES THROUGH INNER CV-FACES
!=============================================

!.....ALL BLOCKS
       DO NBL=1,NBLOCK
          CALL SETIND(NBL)
!.....EAST CELL - FACE
          CALL CELQ(NI-2,NJM,NKM,NJ,1,NIJ,IFI,FI, &
                    GAM_TRAN,FX,AE,AW,F_1,2)
!.....NORTH CELL - FACE
          CALL CELQ(NIM,NJ-2,NKM,1,NIJ,NJ,IFI,FI, &
                    GAM_TRAN,FY,AN,AS,F_2,1)
!.....TOP   CELL - FACE
          CALL CELQ(NIM,NJM,NK-2,NIJ,NJ,1,IFI,FI, &
                     GAM_TRAN,FZ,AT,AB,F_3,3)

       END DO     !   END ALL BLOCKS

       DO NBL=1,NBLOCK

          CALL SETIND(NBL)
          CALL GIVE_BACK(FI,FI_TEMP)
!
!==================================
!.....IMPLEMENT BOUNDARY CONDITIONS
!==================================
!
          CALL MODQ(NBL,FI)

       END DO     !   END ALL BLOCKS

!=========================================================C
!     UNSTEADY TERM CONTRIBUTION                          C
!                                                         C
!                                                         C
!      SOME EXPLICIT METHODS LIKE EE, RKTVD               C
!      WILL(SHOULD) BE ADDED LATER                        C
!=========================================================C

       IF(LTIME) THEN
          APT=VOL*DTR
          SU=SU+APT*((1.+GAMT_Q)*FIO -0.5*GAMT_Q*FIOO)
          SP=SP+APT*(1.+0.5*GAMT_Q)
       END IF

!==================================================================
!.....FINAL COEFFICIENT AND SOURCES MATRIX FOR SCALAR-EQUATION(S)
!==================================================================
!
!  SOME EXTRA TERM IS REQUIRED TO ADD TO AP AS THE BACKGROUND FIELD
!  FOR TRANSPORT THE CHARGE DENSITY IS NON-DIVERGENCE
!  FOR TEMPERATURE-DEPENDENT ION MOBILITY, SOME IMPROVMENT IS NEEDED
!
!
       CALL GRADFI(GR1X,GR1Y,GR1Z,ESP_K)
       CALL GRADFI(GR2X,GR2Y,GR2Z,VP,IVP)

!     AP_ADD = R*(CINJECTION*VOL*FIO  &          ! PART_1: K
!            -(GR1X*GR2X+GR1Y*GR2Y+   & ! PART_2: GRAD_K
!              GR1Z*GR2Z))

       AP_ADD = R*CINJECTION*VOL*FIO

       AP=-AE-AW-AN-AS-AT-AB+SP+AP_ADD

       AP=AP*URFRS
       SU=SU+URFMS*AP*FI

!====================================
!     SOLVE LINEAR SYSTEM FOR FI
!====================================
       CALL SIPSOL(FI,IFI)
       CALL TEST(FI,FI_TEMP,PAS_FINI,IT_DFC,IFI)
       FI_TEMP=FI

    END DO     !! DFC LOOP END !
    DEALLOCATE(FI_TEMP,APT,AP_ADD)

    CALL BPRES_CHARGE(FI)

!=====================================
!     EXCHANGE VARIABLE BETWEEN BLOCKS
!=====================================

!      CALL CHVSCA(FI)
  END SUBROUTINE TRAN3D

!****************************************
  SUBROUTINE NEWFLX2 !*
!****************************************
    USE GEOM
    IMPLICIT NONE

    INTEGER(IK) :: NBL,I,INP,ISIDE,IND
    REAL(RK)    :: UE,VE,WE
!
!.....ESTIMATE NEW FLUX
!
!===================================
!.....ONLY INNER SURFACES
!===================================

    DO NBL=1,NBLOCK

       CALL SETIND(NBL)
!.....EAST, NORTH, TOP, CELL - FACE

    CALL FLX2(NI-2,NJM,NKM,NJ,1,NIJ,FX,F_1,2)   !! F_1 : FOR NEW FLUX
    CALL FLX2(NIM,NJ-2,NKM,1,NIJ,NJ,FY,F_2,1)   !! F_2
    CALL FLX2(NIM,NJM,NK-2,NIJ,NJ,1,FZ,F_3,3)   !! F_3

!===================================
!.....BLOCK INTERFACES
!===================================

     DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
        INP=NBL_ST(NBL)+NUM_CFI(I,1)
        IND=NBL_ST(NBL)+NUM_CFI(I,3)
        ISIDE=NUM_CFI(I,2)
!
        CALL INIT_GEOM_BLOCK(ISIDE,I)
!
!.....INTERPOLATION FACTOR
        FXE=FX_CFI(I)
        FXW=1.-FXE
!
!.....VELOCITIES
!
        UE=U_(INP)*FXW+U_(IND)*FXE
        VE=V_(INP)*FXW+V_(IND)*FXE
        WE=W_(INP)*FXW+W_(IND)*FXE
!
!.....MASS FLUX
        FLX2_CFI(I)= UE*ARX+VE*ARY+WE*ARZ    !! FLX --> FLX2
!
     END DO     !   END BLOCK INTERFACES
    END DO     !   END ALL BLOCKS
!
  END SUBROUTINE NEWFLX2

!***********************************************
  SUBROUTINE FLX2(NIE,NJE,NKE,IDEW,IDNS,IDTB, &
                  FIF,FCF,DIR1)
!***********************************************
    IMPLICIT NONE
!
    REAL(RK)   ,INTENT(IN)    :: FIF(NXYZA)
    REAL(RK)   ,INTENT(INOUT) :: FCF(NXYZA)
    INTEGER(IK),INTENT(IN)    :: NIE,NJE,NKE,IDEW,IDNS,IDTB,DIR1

    REAL(RK)    :: UE_,VE_,WE_
    INTEGER(IK) :: I,J,K
    REAL(RK)    :: FXE,FXW
    INTEGER(IK) :: INE,INP
!
!.....FLUX AT THE CELL FACE
!
!.....CALCULATE EAST CELL FACE
!
    DO K=2,NKE
       DO I=2,NIE
          DO J=2,NJE
!
             INP=LK(K)+LI(I)+J
             INE=INP+IDEW
!.....INTERPOLATION FACTOR
             FXE=FIF(INP)
             FXW=1.-FXE
!
!.....VELOCITIES
!
             UE_=U_(INP)*FXW+U_(INE)*FXE
             VE_=V_(INP)*FXW+V_(INE)*FXE
             WE_=W_(INP)*FXW+W_(INE)*FXE
!
!.....MASS FLUX
             FCF(INP)=(UE_*ARX_T(INP,DIR1) +&
                       VE_*ARY_T(INP,DIR1) +&
                       WE_*ARZ_T(INP,DIR1))
          END DO
       END DO
    END DO
!
  END SUBROUTINE FLX2


!*******************************************************
  SUBROUTINE CELQ(NIE,NJE,NKE,IDEW,IDNS,IDTB, &
                  IFI,FI,GAM_TRAN,FIF,ACFE,ACFW,FCF,DIR1)
!*******************************************************
    IMPLICIT NONE

    REAL(RK)   ,INTENT(INOUT) :: ACFE(NXYZA),ACFW(NXYZA)
    REAL(RK)   ,INTENT(IN)    :: FIF(NXYZA),FCF(NXYZA),FI(NXYZA),GAM_TRAN
    INTEGER(IK),INTENT(IN)    :: NIE,NJE,NKE,IDEW,IDNS,IDTB,IFI,DIR1

    REAL(RK)    :: DE,CE,CW,FLCF,SUEH,SUEL
    REAL(RK)    :: ALPHA,RPTERM,RNTERM,SUADD!,TVD_FUNC
    INTEGER(IK) :: I,J,K
    INTEGER(IK) :: INEE,INEW
    REAL(RK)    :: FXE,FXW
    INTEGER(IK) :: INP,INE
!
!.....CALCULATE EAST,TOP,NORTH  CELL FACE
!
    DO K=2,NKE
       DO I=2,NIE
          DO J=2,NJE
!
             INP=LK(K)+LI(I)+J
             INE=INP+IDEW
!.....INTERPOLATION FACTOR
             FXE=FIF(INP)
             FXW=1.-FXE
!
!.....DIFFUSION COEFFICIENT

             DE = GAM_TRAN*SQRT(ARE2_T(INP,DIR1)/ARKSI2_T(INP,DIR1))
!
!.....CONVECTION FLUXES - UDS

             FLCF=FCF(INP)
             CE=MIN( FLCF,R_0_0)
             CW=MAX( FLCF,R_0_0)
!
             ACFE(INP)= -DE + CE
             ACFW(INE)= -DE - CW
!.....
!.....     D I F F U S I O N
!.....
             SUEH=GAM_TRAN*((GR1X(INP)*FXW+GR1X(INE)*FXE)*ARX_T(INP,DIR1)+ &
                            (GR1Y(INP)*FXW+GR1Y(INE)*FXE)*ARY_T(INP,DIR1)+ &
                            (GR1Z(INP)*FXW+GR1Z(INE)*FXE)*ARZ_T(INP,DIR1))
!.....
             SUEL=DE*((GR1X(INP)*FXW+GR1X(INE)*FXE)*AKX_T(INP,DIR1)+ &
                      (GR1Y(INP)*FXW+GR1Y(INE)*FXE)*AKY_T(INP,DIR1)+ &
                      (GR1Z(INP)*FXW+GR1Z(INE)*FXE)*AKZ_T(INP,DIR1))

!.....
!.....     C O N V E C T I O N
!.....


             INEE = INE + IDEW
             INEW = INP - IDEW

             IF (FLCF>0.) THEN
                ALPHA = 1.0
             ELSE
                ALPHA = 0.
             END IF

           IF ( abs(FI(INE)-FI(INP))<tiny(1.)) THEN
                RPTERM = 0.
                RNTERM = 0.
             ELSE
                RPTERM=( FI(INP)-FI(INEW))/(FI(INE)-FI(INP))
                RNTERM=(FI(INEE)-FI(INE))/(FI(INE)-FI(INP))
             END IF

             SUADD = 0.5*FLCF*((1.-ALPHA)*TVD_FUNC(RNTERM)-  &
                  ALPHA*TVD_FUNC(RPTERM))*(FI(INE)-FI(INP))

             SU(INP)=SU(INP) + SUADD + (SUEH-SUEL)
             SU(INE)=SU(INE)-SUADD - (SUEH-SUEL)
!
          END DO
       END DO
    END DO
!
  END SUBROUTINE CELQ

!**************************************************
  SUBROUTINE MODQ(NBL,FI)
!**************************************************
!     BOUNDARY CONDITION                          C

!=================================================C
    USE GEOM
    IMPLICIT NONE

    REAL(RK)   ,INTENT(INOUT)  :: FI(NXYZA)
    INTEGER(IK),INTENT(IN)  :: NBL

    REAL(RK)    :: FLUX_Q,CE
    INTEGER(IK) :: IDR,NBCTYP,NBC_S,NBC_E,I,INP,ISIDE
!
!
!.....BLOCK INTERFACES
!
    CALL COMQ(NBL, FI)

    DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)

       NBCTYP=NUM_TYP(IDR)
       IF((NBCTYP==12).OR.(NBCTYP==13)) THEN !!!! FOR Q, INFLOW B.C;
          NBC_S=NUM_SPR(IDR)+1
          NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)

          DO I=NBC_S, NBC_E
             INP=NBL_ST(NBL)+NUM_IND(I,1)
             ISIDE=NUM_IND(I,2)
             CALL INIT_GEOM_BF(ISIDE,INP)
             FLUX_Q=U_(INE)*ARX+V_(INE)*ARY+W_(INE)*ARZ
             CE=MIN(FLUX_Q,R_0_0)
             SP(INP)=SP(INP)-CE              !! OK
             SU(INP)=SU(INP)-CE*FI(INE)      !! OK
          END DO
       END IF
    END DO  !! END ALL BC
  END SUBROUTINE MODQ

!**************************************************
  SUBROUTINE COMQ(NBL,FI)
!**************************************************
    IMPLICIT NONE

    REAL(RK)   ,INTENT(INOUT) :: FI(NXYZA)
    INTEGER(IK),INTENT(IN) :: NBL

    REAL(RK)             :: FLCF,CE,CW,ALPHA,RPTERM,RNTERM,SUADD!,TVD_FUNC
    INTEGER(IK)          :: I,INP,ISIDE,NADD,INW,IND,J,NBL1
    REAL(RK),ALLOCATABLE :: BUF1(:,:,:)

    ALLOCATE(BUF1(MAXVAL(LINK_TAB(:,:)),NBLOCK,2)) ! MAY BE SIZE-OPTIMIZABLE
    CALL CHVSCA(FI,BUF1,2)
    DO NBL1=1,NBLOCK
      J=0
      DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
         IF (NUM_CFI(I,5)==NBL1) THEN
           J=J+1 ;
           INP=NBL_ST(NBL) + NUM_CFI(I,1)
           IND=NBL_ST(NBL) + NUM_CFI(I,3)
           ISIDE=NUM_CFI(I,2)

           FLCF = FLX2_CFI(I)
           CE=MIN( FLCF,R_0_0)
           CW=MAX( FLCF,R_0_0)

           A_E(I) = CE

           IF (FLCF>0.) THEN
              ALPHA = 1.0
           ELSE
              ALPHA = 0.
           END IF

           IF(MOD(ISIDE,2_ik)==0) THEN    !!! ISIDE = 2,4,6

              IF(ISIDE==2) THEN
                 NADD = ID_SIDE(1)
              ELSE IF(ISIDE==4) THEN
                 NADD = ID_SIDE(2)
              ELSE IF(ISIDE==6) THEN
                 NADD = ID_SIDE(3)
              END IF

              INW = INP - NADD

           ELSE                  !!! ISIDE = 1, 3, 5

              IF(ISIDE==1) THEN
                 NADD = ID_SIDE(1)
              ELSE IF(ISIDE==3) THEN
                 NADD = ID_SIDE(2)
              ELSE IF(ISIDE==5) THEN
                 NADD = ID_SIDE(3)
              END IF
              INW = INP + NADD

           END IF

           IF ( abs(FI(IND)-FI(INP))<tiny(1.)) THEN
              RPTERM = 0.
              RNTERM = 0.
           ELSE
              RPTERM = ( FI(INP)     - FI(INW) )/( FI(IND) - FI(INP) )
              RNTERM = ( BUF1(J,NBL1,2) - FI(IND) )/( FI(IND) - FI(INP) )
           END IF

           SUADD = 0.5*FLCF*((1.-ALPHA)*TVD_FUNC(RNTERM)- &
                ALPHA*TVD_FUNC(RPTERM))*(FI(IND)-FI(INP))

           SU(INP) = SU(INP) + SUADD
           SP(INP) = SP(INP) - A_E(I)

         END IF
      END DO

    END DO
   DEALLOCATE(BUF1)
  END SUBROUTINE COMQ

!*********************************************************************
  SUBROUTINE BPRES_CHARGE(FI)
!*********************************************************************
    USE GEOM
    IMPLICIT NONE
    REAL(RK),INTENT(INOUT)    :: FI(NXYZA)

    INTEGER(IK) :: NBL,IDR,NBCTYP,NBC_S,NBC_E,I,INP,ISIDE

    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
       DO IDR=NUM_SDR(NBL)+1,NUM_SDR(NBL)+NUM_DR(NBL)
          NBCTYP=NUM_TYP(IDR)
          NBC_S=NUM_SPR(IDR)+1
          NBC_E=NUM_SPR(IDR)+NUM_PR(IDR)
          IF((NBCTYP==10) .OR.(NBCTYP==11).OR.(NBCTYP==14))THEN    ! NOT DIRICHLET B.C. --> NEUMANN B.C.
             DO I=NBC_S,NBC_E
                INP=NBL_ST(NBL)+NUM_IND(I,1)
                ISIDE=NUM_IND(I,2)
                CALL INIT_GEOM_BF(ISIDE,INP)
                !==================================C
                ! NOT EASY A BETTER WAY FOR UPDATE C
                !==================================C
                FI(INE)=FI(INP)
             END DO
          END IF
       END DO      !   END ALL FACES
    END DO         !   END ALL BLOCKS
  END SUBROUTINE BPRES_CHARGE

!***********************************************
  PURE SUBROUTINE ESP_UPDATE(NXYZA,ESP_P,ESP_K,RA,T,L_PERMIT,N_MOBILITY)
!***********************************************
    IMPLICIT NONE
    INTEGER(IK),INTENT(IN) :: NXYZA
    REAL(RK),INTENT(OUT)   :: ESP_P(NXYZA),ESP_K(NXYZA)
    REAL(RK),INTENT(IN)    :: T(NXYZA),RA,L_PERMIT,N_MOBILITY
!======================================================
!   A LINEAR MODEL IS USED TO UPDATE THE TEMPERATURE  C
!   DEPENDENT PERMITTIVITY AND ION MOBILITY           C
!======================================================
    ESP_P=1.-RA*T*L_PERMIT
    ESP_K=1.+RA*T*N_MOBILITY
  END SUBROUTINE ESP_UPDATE

!***********************************************
  SUBROUTINE ESP_SH(NXYZA,ESP_P,NI,NJ,NK)
 !***********************************************
    
    IMPLICIT NONE
    INTEGER(IK),INTENT(IN) :: NXYZA
    REAL(RK),INTENT(OUT)   :: ESP_P(NXYZA)
    INTEGER(IK) :: I,J,K
    INTEGER(IK) :: INP,NK,NI,NJ
	
!............SH model permittivity values 
!............epsilon_r for air = 1.0..in our grid ( y > 0.0)
!............epsilon_r for dielectric = 2.7....(y < 0.0)

!.....Assign epsilon on CV INPs in air and dielectric
!
    DO K=1,NK
       DO I=1,NI
          DO J=1,NJ
!
             INP=LK(K)+LI(I)+J
  
    IF(YC(INP)>0.0) THEN
	    ESP_P(INP) = 1.0

	     ELSE IF(YC(INP)<0.0) THEN
		 ESP_P(INP) = 2.7	
  	END IF 
	  
          END DO
        END DO
     END DO 		

  END SUBROUTINE ESP_SH	
!*****************************************
  SUBROUTINE GHOST_CELLS(FI,FI_TEMP)
!*****************************************
    IMPLICIT NONE
    REAL(RK),INTENT(INOUT) :: FI(NXYZA)
    REAL(RK),INTENT(IN)    :: FI_TEMP(NXYZA)

    INTEGER(IK)            :: NBL,I,J,K,INP

!==================================C
!....UPDATE THE BOUNDARY VALUES    C
!....FOLLOWING MIRROR POINT METHOD C
!==================================C
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)

       DO J=2, NJM
          DO K=2, NKM
             INP=LK(K)+LI(1)+J
             FI(INP)=2.*FI_TEMP(INP)-FI_TEMP(INP+NJ)
             INP=LK(NK)+LI(NI)+J
             FI(INP)=2.*FI_TEMP(INP)-FI_TEMP(INP-NJ)
          END DO
       END DO

       DO I=2, NIM
          DO J=2, NJM
             INP=LK(1)+LI(I)+J
             FI(INP)=2.*FI_TEMP(INP)-FI_TEMP(INP+NIJ)
             INP=LK(NK)+LI(I)+J
             FI(INP)=2.*FI_TEMP(INP)-FI_TEMP(INP-NIJ)
          END DO
       END DO

       DO I=2, NIM
          DO K=2, NKM
             INP=LK(K)+LI(I)+1
             FI(INP)=2.*FI_TEMP(INP)-FI_TEMP(INP+1)
             INP=LK(K)+LI(I)+NJ
             FI(INP)=2.*FI_TEMP(INP)-FI_TEMP(INP-1)
          END DO
       END DO
    END DO  !   END ALL BLOCKS
  END SUBROUTINE GHOST_CELLS

!*******************************************
  SUBROUTINE GIVE_BACK(FI,FI_TEMP)
!*******************************************
    IMPLICIT NONE
    REAL(RK),INTENT(INOUT) :: FI(NXYZA)
    REAL(RK),INTENT(IN)    :: FI_TEMP(NXYZA)

    INTEGER(IK) :: I,J,K,INP
    DO J=2, NJM
       DO K=2, NKM
          INP=LK(K)+LI(1)+J
          FI(INP)=FI_TEMP(INP)
          INP=LK(NK)+LI(NI)+J
          FI(INP)=FI_TEMP(INP)
       END DO
    END DO

    DO I=2, NIM
       DO J=2, NJM
          INP=LK(1)+LI(I)+J
          FI(INP) = FI_TEMP(INP)
          INP=LK(NK)+LI(I)+J
          FI(INP) = FI_TEMP(INP)
       END DO
    END DO

    DO I = 2, NIM
       DO K = 2, NKM
          INP=LK(K)+LI(I)+1
          FI(INP) = FI_TEMP(INP)
          INP=LK(K)+LI(I)+NJ
          FI(INP) = FI_TEMP(INP)
       END DO
    END DO
  END SUBROUTINE GIVE_BACK

!****************************************************
  SUBROUTINE WRUNST2PLT
!****************************************************
    IMPLICIT NONE
    INTEGER(IK) :: NBL,I,J,K,I1,I2,IC,JC,KC
    INTEGER(IK) :: INP1,INP2,INP,J1,J2,K1,K2
    REAL(RK) :: TMP,XREF,YREF,ZREF,L1,L2,L3,L4,L5
!
!.....TO OBTAIN VALUES ON CORNOR AND LINE BOUNDARY POINTS
!
    CALL BC_SHOW2

    XREF=0.5
    YREF=0.5
    ZREF=0.5
    TMP=0.

! YZ PLAN ---> X=XREF
! *******************

    IF(LPRI) THEN
       WRITE(IURESX,*) 'VARIABLES = "Y","Z","U","V","W","T","Q","VP"'
       DO NBL=1,NBLOCK
          CALL SETIND(NBL)
          INP1=LK(1)+LI(1)+1
          INP2=LK(NK)+LI(NI)+NJ
          L1=abs(XREF-XC(INP1))
          L2=abs(XC(INP2)-XREF)
          L3=abs(XC(INP2)-XC(INP1))
          IF((L1+L2-L3)<tiny(1.)) THEN    ! IF THERE IS DATA INSIDE THE CURRENT BLOCK
          WRITE(IURESX,1) 'ZONE J = ', NJ,', K =', NK,',F=POINT'
          DO I1=1,NI
             I2 = min(I1+1,NI)         ! THE LAST POINT IS TREATED DIFFERENTLY
             DO K=1,NK
                DO J=1,NJ
                   INP1=LK(K)+LI(I1)+J
                   INP2=LK(K)+LI(I2)+J
                   L1=abs(XREF-XC(INP1))
                   L2=abs(XC(INP2)-XREF)
                   L3=abs(XC(INP2)-XC(INP1))
                   IF (I1==I2.and.L1<tiny(1.)) then   ! IF XREF IS ON THE LAST POINT
                      L1=0.                           ! NO INTERPOLATION
                      L2=1.
                      L3=1.
                   endif
                   IF((L1+L2-L3)<tiny(1.).and. &      ! IF XREF IS BETWEEN I1 and I2
                       L2>tiny(1.)) then              ! BUT NOT ON THE RIGHT
                   
                     IF(LCAL(IEN)) TMP=(L2* T(INP1)+L1* T(INP2))/L3
                     WRITE(IURESX,'(8(2X,E13.6))')   (L2*YC(INP1)+L1*YC(INP2))/L3, &
                                       (L2*ZC(INP1)+L1*ZC(INP2))/L3, &
                                       (L2* U(INP1)+L1* U(INP2))/L3, &
                                       (L2* V(INP1)+L1* V(INP2))/L3, &  ! SAVE DATA WITH INTERPOLATION
                                       (L2* W(INP1)+L1* W(INP2))/L3, &
                                       (       TMP       )         , &
                                       (L2* Q(INP1)+L1* Q(INP2))/L3, &
                                       (L2*VP(INP1)+L1*VP(INP2))/L3
                   ENDIF
                END DO
             END DO
          END DO
          ENDIF
       END DO   ! END BLOCKS
    END IF

! XZ PLAN---> Y=YREF
! *****************

    IF(LPRJ) THEN
       WRITE(IURESY,*) 'VARIABLES = "X","Z","U","V","W","T","Q","VP"'
       DO NBL=1,NBLOCK
          CALL SETIND(NBL)
          INP1=LK(1)+LI(1)+1
          INP2=LK(NK)+LI(NI)+NJ
          L1=abs(YREF-YC(INP1))
          L2=abs(YC(INP2)-YREF)
          L3=abs(YC(INP2)-YC(INP1))
          IF((L1+L2-L3)<tiny(1.)) THEN    ! IF THERE IS DATA INSIDE THE CURRENT BLOCK
          WRITE(IURESY,1) 'ZONE I =', NI,',K = ', NK,',F=POINT'
          DO J1=1,NJ
             J2 = min(J1+1,NJ)         ! THE LAST POINT IS TREATED DIFFERENTLY
             DO K=1,NK
                DO I=1,NI
                   INP1=LK(K)+LI(I)+J1
                   INP2=LK(K)+LI(I)+J2
                   L1=abs(YREF-YC(INP1))
                   L2=abs(YC(INP2)-YREF)
                   L3=abs(YC(INP2)-YC(INP1))
                   IF (J1==J2.and.L1<tiny(1.)) then   ! IF YREF IS ON THE LAST POINT
                      L1=0.                           ! NO INTERPOLATION
                      L2=1.
                      L3=1.
                   endif
                   IF((L1+L2-L3)<tiny(1.).and. &      ! IF YREF IS BETWEEN J1 and J2
                       L2>tiny(1.)) then              ! BUT NOT ON THE RIGHT
                   
                     IF(LCAL(IEN)) TMP=(L2* T(INP1)+L1* T(INP2))/L3
                     WRITE(IURESY,'(8(2X,E13.6))')   (L2*XC(INP1)+L1*XC(INP2))/L3, &
                                       (L2*ZC(INP1)+L1*ZC(INP2))/L3, &
                                       (L2* U(INP1)+L1* U(INP2))/L3, &
                                       (L2* V(INP1)+L1* V(INP2))/L3, &  ! SAVE DATA WITH INTERPOLATION
                                       (L2* W(INP1)+L1* W(INP2))/L3, &
                                       (       TMP       )         , &
                                       (L2* Q(INP1)+L1* Q(INP2))/L3, &
                                       (L2*VP(INP1)+L1*VP(INP2))/L3
                   ENDIF
                END DO
             END DO
          END DO
          ENDIF
       END DO   ! END BLOCKS
    END IF

! XY PLAN ---> Z=ZREF
! ******************

    IF(LPRK) THEN
       WRITE(IURESZ,*) 'VARIABLES = "X","Y","U","V","W","T","Q","VP"'
       DO NBL=1,NBLOCK
          CALL SETIND(NBL)
          INP1=LK(1)+LI(1)+1
          INP2=LK(NK)+LI(NI)+NJ
          L1=abs(ZREF-ZC(INP1))
          L2=abs(ZC(INP2)-ZREF)
          L3=abs(ZC(INP2)-ZC(INP1))
          IF((L1+L2-L3)<tiny(1.)) THEN    ! IF THERE IS DATA INSIDE THE CURRENT BLOCK
          WRITE(IURESZ,1) 'ZONE I = ', NI,',J = ', NJ,',F=POINT'
          DO K1=1,NK
             K2 = min(K1+1,NK)         ! THE LAST POINT IS TREATED DIFFERENTLY
             DO J=1,NJ
                DO I=1,NI
                   INP1=LK(K1)+LI(I)+J
                   INP2=LK(K2)+LI(I)+J
                   L1=abs(ZREF-ZC(INP1))
                   L2=abs(ZC(INP2)-ZREF)
                   L3=abs(ZC(INP2)-ZC(INP1))
                   IF (K1==K2.and.L1<tiny(1.)) then   ! IF ZREF IS ON THE LAST POINT
                      L1=0.                           ! NO INTERPOLATION
                      L2=1.
                      L3=1.
                   endif
                   IF((L1+L2-L3)<tiny(1.).and. &      ! IF ZREF IS BETWEEN K1 and K2
                       L2>tiny(1.)) then              ! BUT NOT ON THE RIGHT
                   
                     IF(LCAL(IEN)) TMP=(L2* T(INP1)+L1* T(INP2))/L3
                     WRITE(IURESZ,'(8(2X,E13.6))')   (L2*XC(INP1)+L1*XC(INP2))/L3, &
                                       (L2*YC(INP1)+L1*YC(INP2))/L3, &
                                       (L2* U(INP1)+L1* U(INP2))/L3, &
                                       (L2* V(INP1)+L1* V(INP2))/L3, &  ! SAVE DATA WITH INTERPOLATION
                                       (L2* W(INP1)+L1* W(INP2))/L3, &
                                       (       TMP       )         , &
                                       (L2* Q(INP1)+L1* Q(INP2))/L3, &
                                       (L2*VP(INP1)+L1*VP(INP2))/L3
                   ENDIF
                END DO
             END DO
          END DO
          ENDIF
       END DO   ! END BLOCKS
    END IF

!
!.....OUTPUT...LINES --> ODD FOR NICV, NJCV AND NKCV
!
    WRITE(85,*) 'VARIABLES = "X","U","V","W","T","Q","VP"'
    WRITE(86,*) 'VARIABLES = "Y","U","V","W","T","Q","VP","EY"'
    WRITE(87,*) 'VARIABLES = "Z","U","V","W","T","Q","VP"'
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
       WRITE(85,*) 'ZONE I =', NI,',F=POINT'
       WRITE(86,*) 'ZONE J =', NJ,',F=POINT'
       WRITE(87,*) 'ZONE K =', NK,',F=POINT'
       IC=(1+NI)/2     !! NUMBER OF CVS: ODD --> 49 CVS;
       JC=(1+NJ)/2     !! NUMBER OF CVS: ODD --> 49 CVS;
       KC=(1+NK)/2     !! NUMBER OF CVS: ODD --> 49 CVS;
!...X LINE...
       DO I=1,NI
          INP=LK(KC)+LI(I)+JC
          IF(LCAL(IEN)) TMP=T(INP)
          WRITE(85,'(7(2X,E13.6))')XC(INP),U(INP),V(INP),W(INP),TMP,Q(INP),VP(INP)
       END DO
!...Y LINE...
       DO J=1,NJ
          INP=LK(KC)+LI(IC)+J
          IF(LCAL(IEN)) TMP=T(INP)
          WRITE(86,'(8(2X,E13.6))')YC(INP),U(INP),V(INP),W(INP),TMP,Q(INP),VP(INP) &
                    ,EY(INP)
       END DO
!...Z LINE...
       DO K=1,NK
          INP=LK(K)+LI(IC)+JC
          IF(LCAL(IEN)) TMP=T(INP)
          WRITE(87,'(7(2X,E13.6))')ZC(INP),U(INP),V(INP),W(INP),TMP,Q(INP),VP(INP)
       END DO
    END DO   ! END BLOCKS
    
    

!*********************************
! WRITE ALL RESULTS INTO FIELD.PLT
!*********************************







     CALL WRITE_TECPLOT_FIELD 
     
     
     CALL WRITE_VTK_FIELD

    ! restore values for computation
    CALL CHVVEC(U,V,W)
    CALL CHVVEC(XC,YC,ZC)
    CALL CHVSCA(Q)
  !  CALL CHVSCA(VP)
    CALL CHVSCA(VP,alias=7)
    IF(LCAL(IEN)) CALL CHVSCA(T)
    
   
  1   FORMAT(A12,I4,A9,I4,A8)
  
  
  
  
  
  END SUBROUTINE WRUNST2PLT
  
  
  !****************************************************
SUBROUTINE WRITE_VTK_FIELD
!****************************************************
  USE DEC
  IMPLICIT NONE

  INTEGER(IK) :: NBL, I, J, K, INP
  INTEGER(IK) :: IVTK
  REAL(RK)    :: TMP
  CHARACTER(LEN=256) :: FNAME

  IVTK  = 79
  FNAME = 'field.vtk'

  IF (NBLOCK /= 1) THEN
     WRITE(*,*) 'WRITE_VTK_FIELD: only NBLOCK=1 is supported for field.vtk'
     WRITE(*,*) 'Current NBLOCK = ', NBLOCK
     RETURN
  END IF

  OPEN(UNIT=IVTK, FILE=FNAME, STATUS='REPLACE', FORM='FORMATTED', ACTION='WRITE')

  NBL = 1
  CALL SETIND(NBL)

  WRITE(IVTK,'(A)') '# vtk DataFile Version 3.0'
  WRITE(IVTK,'(A)') 'Oracle3D field export'
  WRITE(IVTK,'(A)') 'ASCII'
  WRITE(IVTK,'(A)') 'DATASET STRUCTURED_GRID'
  WRITE(IVTK,'(A,1X,I0,1X,I0,1X,I0)') 'DIMENSIONS', NI, NJ, NK
  WRITE(IVTK,'(A,1X,I0,1X,A)') 'POINTS', NI*NJ*NK, 'float'

  !----------------------------------------
  ! Grid points
  ! Order consistent with structured layout
  !----------------------------------------
  DO K = 1, NK
     DO I = 1, NI
        DO J = 1, NJ
           INP = LK(K) + LI(I) + J
           WRITE(IVTK,'(3(1X,ES20.12))') REAL(XC(INP),REAL4), REAL(YC(INP),REAL4), REAL(ZC(INP),REAL4)
        END DO
     END DO
  END DO

  WRITE(IVTK,'(A,1X,I0)') 'POINT_DATA', NI*NJ*NK

  !----------------------------------------
  ! Velocity vector
  !----------------------------------------
  WRITE(IVTK,'(A)') 'VECTORS Velocity float'
  DO K = 1, NK
     DO I = 1, NI
        DO J = 1, NJ
           INP = LK(K) + LI(I) + J
           WRITE(IVTK,'(3(1X,ES20.12))') REAL(U(INP),REAL4), REAL(V(INP),REAL4), REAL(W(INP),REAL4)
        END DO
     END DO
  END DO

  !----------------------------------------
  ! Pressure
  !----------------------------------------
  WRITE(IVTK,'(A)') 'SCALARS Pressure float 1'
  WRITE(IVTK,'(A)') 'LOOKUP_TABLE default'
  DO K = 1, NK
     DO I = 1, NI
        DO J = 1, NJ
           INP = LK(K) + LI(I) + J
           WRITE(IVTK,'(1X,ES20.12)') REAL(P(INP),REAL4)
        END DO
     END DO
  END DO

  !----------------------------------------
  ! Temperature
  !----------------------------------------
  WRITE(IVTK,'(A)') 'SCALARS Temperature float 1'
  WRITE(IVTK,'(A)') 'LOOKUP_TABLE default'
  DO K = 1, NK
     DO I = 1, NI
        DO J = 1, NJ
           INP = LK(K) + LI(I) + J
           TMP = 0.0_RK
           IF (LCAL(IEN)) TMP = T(INP)
           WRITE(IVTK,'(1X,ES20.12)') REAL(TMP,REAL4)
        END DO
     END DO
  END DO

  !----------------------------------------
  ! Charge density Q
  !----------------------------------------
  WRITE(IVTK,'(A)') 'SCALARS Q float 1'
  WRITE(IVTK,'(A)') 'LOOKUP_TABLE default'
  DO K = 1, NK
     DO I = 1, NI
        DO J = 1, NJ
           INP = LK(K) + LI(I) + J
           WRITE(IVTK,'(1X,ES20.12)') REAL(Q(INP),REAL4)
        END DO
     END DO
  END DO

  !----------------------------------------
  ! Electric potential VP
  !----------------------------------------
  WRITE(IVTK,'(A)') 'SCALARS VP float 1'
  WRITE(IVTK,'(A)') 'LOOKUP_TABLE default'
  DO K = 1, NK
     DO I = 1, NI
        DO J = 1, NJ
           INP = LK(K) + LI(I) + J
           WRITE(IVTK,'(1X,ES20.12)') REAL(VP(INP),REAL4)
        END DO
     END DO
  END DO

  CLOSE(IVTK)

END SUBROUTINE WRITE_VTK_FIELD


!****************************************************
  SUBROUTINE BC_SHOW2()
!****************************************************
    IMPLICIT NONE
    INTEGER(IK) :: NBL,I,i1,j1,k1,INP1,INP2,FACE,j,k,test_edge
    INTEGER(IK) :: inc(3),dir,INP
    REAL(RK) :: FI1,FI2,d1,d2,d3
    REAL(RK),ALLOCATABLE :: BUF1(:,:,:)

    !interpolate on interfaces
    
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
       
       DO I=NUM_SCF(NBL)+1,NUM_SCF(NBL)+NUM_CF(NBL)
         INP1  =NBL_ST(NBL)+NUM_CFI(I,1)
         INP2  =NBL_ST(NBL)+NUM_CFI(I,3)
         
         FACE=NUM_CFI(I,2)
          
          FI2=FX_CFI(I)
          FI1=1.-FI2
          
         U(INP2) = FI1*U(INP1) +FI2*U(INP2)
         V(INP2) = FI1*V(INP1) +FI2*V(INP2)
         W(INP2) = FI1*W(INP1) +FI2*W(INP2)
         Q(INP2) = FI1*Q(INP1) +FI2*Q(INP2)
         VP(INP2)= FI1*VP(INP1)+FI2*VP(INP2)
         IF(LCAL(IEN)) T(INP2) = FI1*T(INP1)+FI2*T(INP2)
         XC(INP2) = FI1*XC(INP1)+FI2*XC(INP2)
         YC(INP2) = FI1*YC(INP1)+FI2*YC(INP2)
         ZC(INP2) = FI1*ZC(INP1)+FI2*ZC(INP2)
         
      enddo
   enddo
   
   !extrapolate on edges
   
    DO NBL=1,NBLOCK
       CALL SETIND(NBL)
       do k=1,nk
       do j=1,nj
       do i=1,ni
         test_edge=0
         if (i==1) then
           test_edge=test_edge+1
           inc(test_edge)=nj
         elseif (i==ni) then
           test_edge=test_edge+1
           inc(test_edge)=-nj
         endif
         if (j==1) then
           test_edge=test_edge+1
           inc(test_edge)=1
         elseif (j==nj) then
           test_edge=test_edge+1
           inc(test_edge)=-1
         endif
         if (k==1) then
           test_edge=test_edge+1
           inc(test_edge)=nij
         elseif (k==nk) then
           test_edge=test_edge+1
           inc(test_edge)=-nij
         endif
         if (test_edge>=2) then ! corner or line
           INP=LK(K)+LI(I) +J
           U(INP) = 0.
           V(INP) = 0.
           W(INP) = 0.
           Q(INP) = 0.
           VP(INP)=0.
           IF(LCAL(IEN)) T(INP) = 0.
           do dir=1,test_edge
               d1=sqrt((xc(INP)-XC(INP+inc(dir)))**2 + &
                       (yc(INP)-yC(INP+inc(dir)))**2 + &
                       (zc(INP)-zC(INP+inc(dir)))**2 )
               d2=sqrt((xc(INP)-XC(INP+2*inc(dir)))**2 + &
                       (yc(INP)-yC(INP+2*inc(dir)))**2 + &
                       (zc(INP)-zC(INP+2*inc(dir)))**2 )
               U(INP) = U(INP) +d1*(d2* U(INP+inc(dir))/d1- U(INP+2*inc(dir)))/(d2-d1)
               V(INP) = V(INP) +d1*(d2* V(INP+inc(dir))/d1- V(INP+2*inc(dir)))/(d2-d1)
               W(INP) = W(INP) +d1*(d2* W(INP+inc(dir))/d1- W(INP+2*inc(dir)))/(d2-d1)
               Q(INP) = Q(INP) +d1*(d2* Q(INP+inc(dir))/d1- Q(INP+2*inc(dir)))/(d2-d1)
               VP(INP)=VP(INP)+d1*(d2*VP(INP+inc(dir))/d1-VP(INP+2*inc(dir)))/(d2-d1)
               IF(LCAL(IEN)) T(INP) = T(INP)+d1*(d2*T(INP+inc(dir))/d1-T(INP+2*inc(dir)))/(d2-d1)
           enddo
           U(INP) = U(INP)/test_edge
           V(INP) = V(INP)/test_edge
           W(INP) = W(INP)/test_edge
           Q(INP) = Q(INP)/test_edge
           VP(INP)=VP(INP)/test_edge
           IF(LCAL(IEN)) T(INP) = T(INP)/test_edge
         endif
       enddo
       enddo
       enddo
    enddo
   

  END SUBROUTINE BC_SHOW2

!****************************************************
  SUBROUTINE BC_SHOW
!****************************************************
    IMPLICIT NONE
    INTEGER(IK) :: NBL,I,J,K,INP,INP2
!
!.....SET CORNER AND LINE VALUES
!
    DO NBL = 1, NBLOCK
       CALL SETIND(NBL)
!
!......PART 1: 8 CORNER POINTS
!
       I = 1
       J = 1
       K = 1

       INP=LK(K)+LI(I)+J
       INP2 = LK(K+1)+LI(I+1)+J+1

       U(INP) = U(INP2)
       V(INP) = V(INP2)
       W(INP) = W(INP2)
       P(INP) = P(INP2)
IF(LCAL(IEN))       T(INP) = T(INP2)
       Q(INP) = Q(INP2)
       VP(INP) = VP(INP2)

       I = 1
       J = NJ
       K = 1

       INP=LK(K)+LI(I)+J
       INP2 = LK(K+1)+LI(I+1)+J-1

       U(INP) = U(INP2)
       V(INP) = V(INP2)
       W(INP) = W(INP2)
       P(INP) = P(INP2)
IF(LCAL(IEN))       T(INP) = T(INP2)
       Q(INP) = Q(INP2)
       VP(INP) = VP(INP2)

       I = 1
       J = NJ
       K = NK

       INP=LK(K)+LI(I)+J
       INP2 = LK(K-1)+LI(I+1)+J-1

       U(INP) = U(INP2)
       V(INP) = V(INP2)
       W(INP) = W(INP2)
       P(INP) = P(INP2)
IF(LCAL(IEN))       T(INP) = T(INP2)
       Q(INP) = Q(INP2)
       VP(INP) = VP(INP2)

       I = 1
       J = 1
       K = NK

       INP=LK(K)+LI(I)+J
       INP2 = LK(K-1)+LI(I+1)+J+1

       U(INP) = U(INP2)
       V(INP) = V(INP2)
       W(INP) = W(INP2)
       P(INP) = P(INP2)
IF(LCAL(IEN))       T(INP) = T(INP2)
       Q(INP) = Q(INP2)
       VP(INP) = VP(INP2)

       I = NI
       J = 1
       K = 1

       INP=LK(K)+LI(I)+J
       INP2 = LK(K+1)+LI(I-1)+J+1

       U(INP) = U(INP2)
       V(INP) = V(INP2)
       W(INP) = W(INP2)
       P(INP) = P(INP2)
IF(LCAL(IEN))       T(INP) = T(INP2)
       Q(INP) = Q(INP2)
       VP(INP) = VP(INP2)

       I = NI
       J = NJ
       K = 1

       INP=LK(K)+LI(I)+J
       INP2 = LK(K+1)+LI(I-1)+J-1

       U(INP) = U(INP2)
       V(INP) = V(INP2)
       W(INP) = W(INP2)
       P(INP) = P(INP2)
IF(LCAL(IEN))       T(INP) = T(INP2)
       Q(INP) = Q(INP2)
       VP(INP) = VP(INP2)

       I = NI
       J = NJ
       K = NK

       INP=LK(K)+LI(I)+J
       INP2 = LK(K-1)+LI(I-1)+J-1

       U(INP) = U(INP2)
       V(INP) = V(INP2)
       W(INP) = W(INP2)
       P(INP) = P(INP2)
IF(LCAL(IEN))       T(INP) = T(INP2)
       Q(INP) = Q(INP2)
       VP(INP) = VP(INP2)

       I = NI
       J = 1
       K = NK

       INP=LK(K)+LI(I)+J
       INP2 = LK(K-1)+LI(I-1)+J+1

       U(INP) = U(INP2)
       V(INP) = V(INP2)
       W(INP) = W(INP2)
       P(INP) = P(INP2)
IF(LCAL(IEN))       T(INP) = T(INP2)
       Q(INP) = Q(INP2)
       VP(INP) = VP(INP2)
!
!......PART 2: 8 LINE CORNER POINTS
!
!
!.....WEST LINES
!
       DO J=2,NJM

          INP=LK(1)+LI(1)+J
          INP2=LK(2)+LI(2)+J

          U(INP) = U(INP2)  !!! STRANGE BUT WORKS
          V(INP) = V(INP2)
          W(INP) = W(INP2)
          P(INP) = P(INP2)
IF(LCAL(IEN))          T(INP) = T(INP2)
          Q(INP) = Q(INP2)
          VP(INP) = VP(INP2)

          INP=LK(NK)+LI(1)+J
          INP2=LK(NK-1)+LI(2)+J

          U(INP) = U(INP2)  !!! STRANGE BUT WORKS
          V(INP) = V(INP2)
          W(INP) = W(INP2)
          P(INP) = P(INP2)
IF(LCAL(IEN))          T(INP) = T(INP2)
          Q(INP) = Q(INP2)
          VP(INP) = VP(INP2)

       END DO

       DO K=2,NKM

          INP=LK(K)+LI(1)+1
          INP2=LK(K)+LI(2)+2

          U(INP) = U(INP2)  !!! STRANGE BUT WORKS
          V(INP) = V(INP2)
          W(INP) = W(INP2)
          P(INP) = P(INP2)
IF(LCAL(IEN))          T(INP) = T(INP2)
          Q(INP) = Q(INP2)
          VP(INP) = VP(INP2)

          INP=LK(K)+LI(1)+NJ
          INP2=LK(K)+LI(2)+NJ-1

          U(INP) = U(INP2)  !!! STRANGE BUT WORKS
          V(INP) = V(INP2)
          W(INP) = W(INP2)
          P(INP) = P(INP2)
IF(LCAL(IEN))          T(INP) = T(INP2)
          Q(INP) = Q(INP2)
          VP(INP) = VP(INP2)

       END DO

!.....EAST
       DO J=2,NJM

          INP=LK(1)+LI(NI)+J
          INP2=LK(2)+LI(NI-1)+J

          U(INP) = U(INP2)  !!! STRANGE BUT WORKS
          V(INP) = V(INP2)
          W(INP) = W(INP2)
          P(INP) = P(INP2)
IF(LCAL(IEN))          T(INP) = T(INP2)
          Q(INP) = Q(INP2)
          VP(INP) = VP(INP2)

          INP=LK(NK)+LI(NI)+J
          INP2=LK(NK-1)+LI(NI-1)+J

          U(INP) = U(INP2)  !!! STRANGE BUT WORKS
          V(INP) = V(INP2)
          W(INP) = W(INP2)
          P(INP) = P(INP2)
IF(LCAL(IEN))          T(INP) = T(INP2)
          Q(INP) = Q(INP2)
          VP(INP) = VP(INP2)

       END DO

       DO K=2,NKM

          INP=LK(K)+LI(NI)+1
          INP2=LK(K)+LI(NI-1)+2

          U(INP) = U(INP2)  !!! STRANGE BUT WORKS
          V(INP) = V(INP2)
          W(INP) = W(INP2)
          P(INP) = P(INP2)
IF(LCAL(IEN))          T(INP) = T(INP2)
          Q(INP) = Q(INP2)
          VP(INP) = VP(INP2)

          INP=LK(K)+LI(NI)+NJ
          INP2=LK(K)+LI(NI-1)+NJ-1

          U(INP) = U(INP2)  !!! STRANGE BUT WORKS
          V(INP) = V(INP2)
          W(INP) = W(INP2)
          P(INP) = P(INP2)
IF(LCAL(IEN))          T(INP) = T(INP2)
          Q(INP) = Q(INP2)
          VP(INP) = VP(INP2)

       END DO

!
!.....SOUTH
       DO I=2,NIM

          INP=LK(1)+LI(I)+1
          INP2=LK(2)+LI(I)+2

          U(INP) = U(INP2)  !!! STRANGE BUT WORKS
          V(INP) = V(INP2)
          W(INP) = W(INP2)
          P(INP) = P(INP2)
IF(LCAL(IEN))          T(INP) = T(INP2)
          Q(INP) = Q(INP2)
          VP(INP) = VP(INP2)

          INP=LK(NK)+LI(I)+1
          INP2=LK(NK-1)+LI(I)+2

          U(INP) = U(INP2)  !!! STRANGE BUT WORKS
          V(INP) = V(INP2)
          W(INP) = W(INP2)
          P(INP) = P(INP2)
IF(LCAL(IEN))          T(INP) = T(INP2)
          Q(INP) = Q(INP2)
          VP(INP) = VP(INP2)

!
!.....NORTH

          INP=LK(1)+LI(I)+NJ
          INP2=LK(2)+LI(I)+NJ-1

          U(INP) = U(INP2)  !!! STRANGE BUT WORKS
          V(INP) = V(INP2)
          W(INP) = W(INP2)
          P(INP) = P(INP2)
IF(LCAL(IEN))          T(INP) = T(INP2)
          Q(INP) = Q(INP2)
          VP(INP) = VP(INP2)

          INP=LK(NK)+LI(I)+NJ
          INP2=LK(NK-1)+LI(I)+NJ-1

          U(INP) = U(INP2)  !!! STRANGE BUT WORKS
          V(INP) = V(INP2)
          W(INP) = W(INP2)
          P(INP) = P(INP2)
IF(LCAL(IEN))          T(INP) = T(INP2)
          Q(INP) = Q(INP2)
          VP(INP) = VP(INP2)

       END DO

    END DO

  END SUBROUTINE BC_SHOW

!************************************
  SUBROUTINE NONDIMENSION_FORM
    IMPLICIT NONE

    IF (LNONDIMENSION) THEN
       ! On travaille en variables adimensionnées
       DENSIT = 1.0_RK
       VISCOS = 1.0_RK

       ! Le Reynolds est désormais imposé par le fichier d'entrée
       REY = RE_INPUT

       IF (LEHD) THEN
          ! EHD isotherme
          R = 1.0_RK
       END IF

       IF (LETHD) THEN
          ! ETHD
          R = TSTABILITY/(MSTABILITY*MSTABILITY)
       END IF

    ELSE
       ! En mode dimensionnel, le code historique travaille avec REY=1
       REY = 1.0_RK
    END IF

  END SUBROUTINE NONDIMENSION_FORM


!****************************************************
SUBROUTINE WRITE_TECPLOT_FIELD
!****************************************************
  USE GEOM
  IMPLICIT NONE
  INTEGER(IK) :: NBL,I,J,K,INP

  DO NBL=1,NBLOCK
     CALL SETIND(NBL)

     WRITE(78,'(A,ES14.6,A,I0,A,I0,A,I0,A)') &
          'ZONE T="time=', TIME, '", I=', NI, ', J=', NJ, ', K=', NK, ', DATAPACKING=POINT'

     DO K=1,NK
        DO J=1,NJ
           DO I=1,NI
              INP = LK(K) + LI(I) + J
              WRITE(78,'(7(1X,E10.3))') XC(INP), YC(INP), ZC(INP), &
                                        U(INP), V(INP), W(INP), P(INP)
           END DO
        END DO
     END DO
  END DO

  FLUSH(78)

END SUBROUTINE WRITE_TECPLOT_FIELD
  

!****************************************************
  SUBROUTINE WRITE_VORTEX_FIELD
!****************************************************
!
!  Calcule et ecrit dans vortex_field.plt (unite 109) :
!    - Le rotationnel (vorticite) : OMEGAX, OMEGAY, OMEGAZ
!    - Le critere Q
!    - Le critere Lambda2 (2eme valeur propre de S^2 + Omega^2)
!
!  Les gradients de vitesse sont recalcules via GRADFI.
!  GR1X/Y/Z = grad(U),  GR2X/Y/Z = grad(V),  GR3X/Y/Z = grad(W)
!
  USE DEC
  IMPLICIT NONE

  INTEGER(IK) :: NBL, I, J, K, INP

  REAL(RK) :: DUDX, DUDY, DUDZ
  REAL(RK) :: DVDX, DVDY, DVDZ
  REAL(RK) :: DWDX, DWDY, DWDZ

  REAL(RK) :: OMEGAX, OMEGAY, OMEGAZ

  REAL(RK) :: S11, S12, S13, S22, S23, S33
  REAL(RK) :: OM12, OM13, OM23

  REAL(RK) :: QVAL

  REAL(RK) :: M11, M12, M13, M22, M23, M33
  REAL(RK) :: P_COEFF, Q_COEFF, DET_COEFF
  REAL(RK) :: PHI, EIG1, EIG2, EIG3, ETMP
  REAL(RK) :: PP2, QQ2
  REAL(RK), PARAMETER :: THIRD  = 1.0_RK / 3.0_RK
  REAL(RK), PARAMETER :: PI_L2  = 3.14159265358979323846_RK

!---- Calcul des gradients de vitesse ----
  CALL GRADFI(GR1X, GR1Y, GR1Z, U)
  CALL GRADFI(GR2X, GR2Y, GR2Z, V)
  CALL GRADFI(GR3X, GR3Y, GR3Z, W)

  DO NBL = 1, NBLOCK
     CALL SETIND(NBL)

     WRITE(109,'(A,ES14.6,A,I0,A,I0,A,I0,A)') &
          'ZONE T="time=', TIME, '", I=', NI, ', J=', NJ, &
          ', K=', NK, ', DATAPACKING=POINT'

     DO K = 2, NKM
        DO I = 2, NIM
           DO J = 2, NJM

              INP = LK(K) + LI(I) + J

              !---- Composantes du tenseur gradient de vitesse ----
              DUDX = GR1X(INP)
              DUDY = GR1Y(INP)
              DUDZ = GR1Z(INP)
              DVDX = GR2X(INP)
              DVDY = GR2Y(INP)
              DVDZ = GR2Z(INP)
              DWDX = GR3X(INP)
              DWDY = GR3Y(INP)
              DWDZ = GR3Z(INP)

              !============================================
              ! 1) ROTATIONNEL omega = rot(U)
              !    OMEGAX = dW/dy - dV/dz
              !    OMEGAY = dU/dz - dW/dx
              !    OMEGAZ = dV/dx - dU/dy
              !============================================
              OMEGAX = DWDY - DVDZ
              OMEGAY = DUDZ - DWDX
              OMEGAZ = DVDX - DUDY

              !============================================
              ! 2) Tenseurs S (symetrique) et OMEGA (antisymetrique)
              !    S_ij   = 0.5*(dUi/dxj + dUj/dxi)
              !    Om_ij  = 0.5*(dUi/dxj - dUj/dxi)
              !============================================
              S11  =  DUDX
              S22  =  DVDY
              S33  =  DWDZ
              S12  = 0.5_RK * (DUDY + DVDX)
              S13  = 0.5_RK * (DUDZ + DWDX)
              S23  = 0.5_RK * (DVDZ + DWDY)

              OM12 = 0.5_RK * (DUDY - DVDX)
              OM13 = 0.5_RK * (DUDZ - DWDX)
              OM23 = 0.5_RK * (DVDZ - DWDY)

              !============================================
              ! 3) CRITERE Q
              !    Q = 0.5*(||Omega||^2 - ||S||^2)
              !============================================
              QVAL = (OM12**2 + OM13**2 + OM23**2) &
                   - (S11**2 + S22**2 + S33**2        &
                      + 2.0_RK*(S12**2 + S13**2 + S23**2)) * 0.5_RK

              !============================================
              ! 4) CRITERE LAMBDA2
              !    M = S^2 + Omega^2  (matrice 3x3 symetrique)
              !    Lambda2 = 2eme valeur propre de M (tri croissant)
              !    Vortex si lambda2 < 0
              !============================================
              M11 = S11*S11 + S12*S12 + S13*S13  - OM12*OM12 - OM13*OM13
              M22 = S12*S12 + S22*S22 + S23*S23  - OM12*OM12 - OM23*OM23
              M33 = S13*S13 + S23*S23 + S33*S33  - OM13*OM13 - OM23*OM23
              M12 = S11*S12 + S12*S22 + S13*S23  - OM13*OM23
              M13 = S11*S13 + S12*S23 + S13*S33  + OM12*OM23
              M23 = S12*S13 + S22*S23 + S23*S33  - OM12*OM13

              !---- Valeurs propres par methode de Cardano ----
              PP2 = M11 + M22 + M33
              QQ2 = M11*M22 + M11*M33 + M22*M33 &
                  - M12**2 - M13**2 - M23**2
              DET_COEFF = M11*(M22*M33 - M23**2) &
                        - M12*(M12*M33 - M23*M13) &
                        + M13*(M12*M23 - M22*M13)

              P_COEFF = PP2**2 * THIRD - QQ2
              Q_COEFF = PP2 * THIRD * (QQ2 - PP2**2 * THIRD) * (2.0_RK/3.0_RK) &
                      + DET_COEFF

              IF (P_COEFF < 0.0_RK) P_COEFF = 0.0_RK

              IF (P_COEFF > SMALL) THEN
                 PHI  = ACOS(MAX(-1.0_RK, MIN(1.0_RK, &
                        -0.5_RK * Q_COEFF / (P_COEFF * SQRT(P_COEFF)))))
                 EIG1 = PP2*THIRD + 2.0_RK*SQRT(P_COEFF)*COS( PHI                 *THIRD)
                 EIG2 = PP2*THIRD + 2.0_RK*SQRT(P_COEFF)*COS((PHI + 2.0_RK*PI_L2)*THIRD)
                 EIG3 = PP2*THIRD + 2.0_RK*SQRT(P_COEFF)*COS((PHI + 4.0_RK*PI_L2)*THIRD)
              ELSE
                 EIG1 = PP2 * THIRD
                 EIG2 = PP2 * THIRD
                 EIG3 = PP2 * THIRD
              END IF

              !---- Tri croissant : EIG2 = Lambda2 ----
              IF (EIG1 > EIG2) THEN; ETMP=EIG1; EIG1=EIG2; EIG2=ETMP; END IF
              IF (EIG2 > EIG3) THEN; ETMP=EIG2; EIG2=EIG3; EIG3=ETMP; END IF
              IF (EIG1 > EIG2) THEN; ETMP=EIG1; EIG1=EIG2; EIG2=ETMP; END IF

              WRITE(109,'(8(1X,E13.6))') &
                   XC(INP), YC(INP), ZC(INP), &
                   OMEGAX,  OMEGAY,  OMEGAZ,  &
                   QVAL,    EIG2

           END DO
        END DO
     END DO
  END DO

  FLUSH(109)

END SUBROUTINE WRITE_VORTEX_FIELD
!****************************************************
  SUBROUTINE INIT_GORTLER_PERTURBATION
!****************************************************
!
!  Ajoute une perturbation initiale en W pour declencher
!  les instabilites de Gortler.
!
!  La perturbation est absolue, independante de ABS(U),
!  pour ne pas disparaitre si U est nul ou faible
!  au moment de l'initialisation.
!
!****************************************************

  USE DEC
  IMPLICIT NONE

  INTEGER(IK) :: NBL, I, J, K, INP
  REAL(RK)    :: RAND_VAL
  REAL(RK)    :: EPSILON_GORTLER
  REAL(RK)    :: KZ, ZC_MIN, ZC_MAX, LZ_BLK
  REAL(RK)    :: ZLOC
  REAL(RK)    :: WPERT_MIN, WPERT_MAX, WPERT_ABSMAX
  REAL(RK)    :: WADD

  EPSILON_GORTLER = 1.0E-3_RK

  WPERT_MIN    =  1.0E30_RK
  WPERT_MAX    = -1.0E30_RK
  WPERT_ABSMAX =  0.0_RK

  DO NBL = 1, NBLOCK

     CALL SETIND(NBL)

     !---------------------------------------------------------
     ! Calcul de la taille spanwise du bloc courant
     !---------------------------------------------------------
     ZC_MIN = ZC(LK(2)   + LI(2) + 2)
     ZC_MAX = ZC(LK(NKM) + LI(2) + 2)

     LZ_BLK = ZC_MAX - ZC_MIN

     IF (ABS(LZ_BLK) < SMALL) THEN
        WRITE(*,*) 'WARNING INIT_GORTLER_PERTURBATION: LZ_BLK presque nul.'
        LZ_BLK = 1.0_RK
     END IF

     !---------------------------------------------------------
     ! Mode spanwise
     !
     ! lambda_z = Lz/2
     ! donc kz = 2*pi/lambda_z = 4*pi/Lz
     !
     ! Ici on utilise PI venant du module DEC.
     ! Ne pas redeclarer PI localement.
     !---------------------------------------------------------
     KZ = 4.0_RK * PI / LZ_BLK

     DO K = 2, NKM
        DO I = 2, NIM
           DO J = 2, NJM

              INP = LK(K) + LI(I) + J

              ZLOC = ZC(INP) - ZC_MIN

              !------------------------------------------------
              ! Bruit blanc absolu
              !------------------------------------------------
              CALL RANDOM_NUMBER(RAND_VAL)

              WADD = EPSILON_GORTLER * &
                     (2.0_RK * RAND_VAL - 1.0_RK)

              W(INP) = W(INP) + WADD

              !------------------------------------------------
              ! Mode sinusoidal spanwise absolu
              !------------------------------------------------
              WADD = EPSILON_GORTLER * SIN(KZ * ZLOC)

              W(INP) = W(INP) + WADD

              !------------------------------------------------
              ! Diagnostics locaux
              !------------------------------------------------
              WPERT_MIN    = MIN(WPERT_MIN, W(INP))
              WPERT_MAX    = MAX(WPERT_MAX, W(INP))
              WPERT_ABSMAX = MAX(WPERT_ABSMAX, ABS(W(INP)))

           END DO
        END DO
     END DO

  END DO

  WRITE(*,*) ' '
  WRITE(*,*) '================================================='
  WRITE(*,'(A)')        ' ==> GORTLER PERTURBATION INITIALIZED'
  WRITE(*,'(A,ES12.4)') '     epsilon       = ', EPSILON_GORTLER
  WRITE(*,'(A,ES12.4)') '     W min         = ', WPERT_MIN
  WRITE(*,'(A,ES12.4)') '     W max         = ', WPERT_MAX
  WRITE(*,'(A,ES12.4)') '     max |W|       = ', WPERT_ABSMAX
  WRITE(*,*) '================================================='
  WRITE(*,*) ' '

END SUBROUTINE INIT_GORTLER_PERTURBATION


!==================================================
!           FLUX LIMITER FUNCTIONS
!==================================================
!
PURE FUNCTION TVD_FUNC(R)
  IMPLICIT NONE
  REAL(RK),INTENT(IN)  :: R
  REAL(RK)             :: TVD_FUNC

!       TVD_FUNC = 0.                                !UD
!       TVD_FUNC = 1.                                !CD
!       TVD_FUNC = MAX(0.,MIN(2.0, 2.0*R,0.5*(1.0+R))) !MUSCL
  TVD_FUNC = MAX(0.,MIN(4.0*R,(0.75+0.25*R),2.0))!SMART

END FUNCTION TVD_FUNC

END PROGRAM ORACLE3D_4
