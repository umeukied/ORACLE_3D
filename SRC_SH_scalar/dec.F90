!*********************************************************************
!> Define all common variables and functions
!>
!> @author Philippe Traore
!*********************************************************************
!#define REAL_KIND 4
!#define REAL_KIND 8

MODULE DEC

  IMPLICIT NONE

 ! INTEGER KINDS
  INTEGER, PARAMETER :: INT4 = SELECTED_INT_KIND(8)  ! SEEMS OK BUT NOT SURE
  INTEGER, PARAMETER :: INT8 = SELECTED_INT_KIND(10)
! FLOATING POINT KINDS
  INTEGER, PARAMETER :: REAL4 = SELECTED_REAL_KIND(6)
  INTEGER, PARAMETER :: REAL8 = SELECTED_REAL_KIND(15)
  INTEGER, PARAMETER :: REAL16 = SELECTED_REAL_KIND(32)

  INTEGER, PARAMETER :: IK = INT4
#if REAL_KIND == 4
  INTEGER, PARAMETER :: RK = REAL4
#elif REAL_KIND == 8
  INTEGER, PARAMETER :: RK = REAL8
#endif

  INTEGER(IK),PARAMETER  :: NPHI=5,NPHI2=7
  REAL(RK)               :: MSTABILITY, L_PERMIT, N_MOBILITY
  INTEGER(IK)            :: NXMAX,NZMAX,NUM_BF_ALL,NUM_CF_ALL,NXYZA,NUM_DR_ALL=35

  REAL(RK), PARAMETER :: PI = 4 *atan(1.0)
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * *      INC_IND    * * * * * * * * * * * * * * * *
!
!INDEX
  INTEGER(IK),ALLOCATABLE ::   NIBL(:),NJBL(:), NKBL(:), &
                             NBL_ST(:),  LI(:),   LK(:)
  INTEGER(IK)             :: NSW(NPHI2),ITSWEEP(NPHI2), &
                             NCVALL,NI ,NJ ,NK ,NIM,NJM,NKM, &
                             NIJ,NIK,NJK,NIJK, &
                             NBLOCK, &
                             IU,IV,IW,IP,IEN, &
                             IQ,IVP, &
                             ICST,ICEN, &
                             ICSTALL,ICENALL, &
                             ITIM,NTM,ITIMS,ITIME,IMA, &
                             ITER,NSWEEP,LS, &
                             NBLMON,IMON,JMON, &
                             KMON,IJKMON,IJKMO, &
                             NBLPR,IPR ,JPR , &
                             KPR ,IJKPR,IJKPE
!IREAL
  REAL(RK),ALLOCATABLE ::  RM1(:),RM2(:)
  REAL(RK)             ::  SLARGE,SORMAX,SORMAX2, &
                           SOURCE,SMALL,GREAT, &
                           RESOR(NPHI2),SNORIN(NPHI2), &
                           ALFA, &
                           DT,DTR,GAMT,TIME, &
                           DENSIT,VISCOS, &
                           PRTINV(NPHI),PRTR, &
                           UIN,VIN,WIN,TIN, &
                           GDS(NPHI2),URF(NPHI2),SOR(NPHI2), &
                           URFR(NPHI2),URFM(NPHI2), &
                           ESP_DFC(NPHI2), &
                           FLOWIN,XMONIN,TENOM,EDNOM, &
                           FLOW,FACOUT, &
                           R_0_0, &
                           RE, &
                           RA, &
                           RMM1, RMM2, &
                           AMPLITUDE

!SCHEME
  INTEGER(IK)       :: ITIME_U, ITIME_T, &
                       I_SC, ITIME_Q
  REAL(RK)          :: GAMT_U,  GAMT_T, GAMT_Q

!ILOGI
  LOGICAL       :: LREAD,LWRITE,LTIME,LCAL(NPHI), &
                   LSTORE,LTEST,LSCREEN,LRESAVE, &
                   LPRI,LPRJ,LPRK, &
                   LONGEO, &
                   LNONDIMENSION, &
                   LUVWMAX, &
                   LNUSSELT, &
                   LEHD, &
                   LHYDROSTATIC, &
                   LETHD, &
                   LESP

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * *      INC_CHR    * * * * * * * * * * * * * * * *
!
!CHCOMM
  CHARACTER(LEN=256) :: FILIN,FILOUT,FILRES,TITLE
  CHARACTER(LEN=20)  :: FILGRD
  CHARACTER(LEN=3)   :: CHVAR(NPHI2)
!IFUNIT
  INTEGER(IK)       :: IUIN,IUOUT,IUGRD,IURES,IURESX,IURESY,IURESZ, &
                       IUCHK
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * *      INC_TIM    * * * * * * * * * * * * * * * *
!
!TIMCOM
  REAL(RK),ALLOCATABLE ::   UO(:),  VO(:), WO(:), &
                           UOO(:), VOO(:),WOO(:), &
                            TO(:), TOO(:)
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * *      INC_UVW    * * * * * * * * * * * * * * * *
!
!UVW
  REAL(RK),ALLOCATABLE ::     U(:),   V(:),   W(:), &
                             F1(:),  F2(:),  F3(:), &
                              P(:),  PP(:),   T(:), &
                            VIS(:), DEN(:), &
                           GR1X(:),GR1Y(:),GR1Z(:), &
                           GR2X(:),GR2Y(:),GR2Z(:), &
                           GR3X(:),GR3Y(:),GR3Z(:)

!UVWMAX
  REAL(RK)          :: UMAX, VMAX, WMAX, UVW_NORM_MAX, &
                       XMAX_NORM, YMAX_NORM, ZMAX_NORM, &
                       XMAX_U, YMAX_U, ZMAX_U, &
                       XMAX_V, YMAX_V, ZMAX_V, &
                       XMAX_W, YMAX_W, ZMAX_W
  INTEGER(IK)       :: IMAX_NORM, JMAX_NORM, KMAX_NORM, &
                       IMAX_U, JMAX_U, KMAX_U, &
                       IMAX_V, JMAX_V, KMAX_V, &
                       IMAX_W, JMAX_W, KMAX_W

!NUSSELT
  REAL(RK)          :: R_NU_0, R_NU_05, R_NU_10

!EHD

  REAL(RK),ALLOCATABLE ::      U_(:),   V_(:), W_(:), &
                               F_1(:),  F_2(:),F_3(:), &
                               Q(:),   QO(:),QOO(:), &
                              VP(:),  VPO(:), &
                              EX(:),   EY(:), EZ(:), &
                    ! Body force vector for SH model FXsh=Q*Ex
                              FXsh(:),   FYsh(:), FZsh(:),Fmag(:), &
                              EXO(:),  EYO(:),EZO(:), &
                              ESP_P(:),ESP_K(:), &
                              SEHD_C(:)   ! source term for the momentum equation 

  REAL(RK)             ::  REY, D_CHARGE, R, &
                           CINJECTION,TSTABILITY, &
                    !! NEW ADDITIONS for S-H model 
                           LAMBDA, LAMBDAR, &   ! THE DEBYE LENGTH and LAMBDAR = 1/(DEBYE LENGTH)
                           D_C, &               ! D_C (!DImensionless coefficient of the SUZEN-Huang Model in the momentum equation for Coulomb force)       
                           F_BASE, VPmax, &            ! F_BASE (applied voltage frequency-- AC) and amplitude of VP
                           MuLOCATION, sigma, sigmaR, & ! Mu=location parameter, sigma= decay rate scale paramter, sigmaR= 1/sigma 
                           Qmax                  ! Maximum charge density on surface - from experiments 


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * *      INC_BND    * * * * * * * * * * * * * * * *
!
!TYPBC
  INTEGER(IK)             :: ID_ODER(6,3),ID_SIDE(3),ID_SIGN(6)
  INTEGER(IK),ALLOCATABLE ::  NUM_DR(:), &
                             NUM_SDR(:), &
                             NUM_IND(:,:), &
                              NUM_CF(:), &
                             NUM_SCF(:), &
                             NUM_CFI(:,:), &
                             NUM_NUM(:), &
                             NUM_TYP(:), &
                             NUM_PR (:), &
                             NUM_SPR(:)
!RTYPBC
  REAL(RK),ALLOCATABLE ::   FLX_IND(:), &
                            FLX_CFI(:), &
                             FX_CFI(:), &
                                A_E(:), &
                              X_CFI(:,:), &
                              Y_CFI(:,:), &
                              Z_CFI(:,:), &
                           FLX2_CFI(:)
  INTEGER(IK)       ::     NDREG

!BOUNDC
  REAL(RK)          :: PRANL,PRANT, &
                       SUU,SUV,SUW,SUP,SVP,SWP

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * *      INC_GEO    * * * * * * * * * * * * * * * *
!
!GEO
  REAL(RK),ALLOCATABLE ::    X(:), Y(:), Z(:), &
                            XC(:),YC(:),ZC(:), &
                            FX(:),FY(:),FZ(:), &
                           VOL(:)
  REAL(RK),ALLOCATABLE ::  ARX_T(:,:),   ARY_T(:,:), ARZ_T(:,:), &
                           AKX_T(:,:),   AKY_T(:,:),AKZ_T(:,:), &
                          ARE2_T(:,:),ARKSI2_T(:,:), ARE_T(:,:), &
                          DELN_T(:,:), DELNR_T(:,:),ARER_T(:,:)
!BPCHAR
CHARACTER(LEN=6) :: CSIDE(6)
CHARACTER(LEN=1) :: CSIDSE(2,6)

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * *      INC_COE    * * * * * * * * * * * * * * * *
!COEF
  REAL(RK),ALLOCATABLE ::   AE(:), AW(:), AN(:),AS(:),AT(:), &
                            AB(:), AP(:), SU(:),SV(:),SW(:), &
                            AP_ADD(:), &   ! AP_ADD = VOL*(LAMBDAR**2) for SH model
                            SP(:),RES(:),SUO(:), &
                           APU(:),APV(:),APW(:), &
                           SPU(:),SPV(:),SPW(:),BIG_BUF(:,:,:,:)
  INTEGER(IK),ALLOCATABLE :: LINK_TAB(:,:)
!
!_______   END

CONTAINS

!**********************
  SUBROUTINE INITDATA
!**********************
!..................................................................
!    IN THIS BLOCK DATA VALUES OF SOME INTEGERS
!    AND CHARACTERS ARE SET
!------------------------------------------------------------------
!
    CHVAR =(/' U ',' V ',' W ',' P ',' EN',' Q ',' VP'/)
!.....VARIABLE IDENTIFIERS
!.....IU  = 1: X VELOCITY
!.....IV  = 2: Y VELOCITY
!.....IW  = 3: Z VELOCITY
!.....IP  = 4: PRESSURE
!.....IEN = 5: ENERGY
!.....IQ  = 6: CHARGE DENSITY
!.....IVP = 7: ELECTRIC POTENTIAL
    IU=1;IV=2;IW=3;IP=4;IEN=5;IQ=6;IVP=7

!.....  WEST     1
!.....  EAST     2
!.....  SOUTH    3
!.....  NORTH    4
!.....  BOTTOM   5
!.....  TOP      6

    ID_ODER=RESHAPE( (/ 1, 1, 2, 2, 3, 3, &
                        2, 2, 3, 3, 1, 1, &
                        3, 3, 1, 1, 2, 2/)&
                   ,(/6,3/))
    ID_SIGN=(/ -1, &
                1, &
               -1, &
                1, &
               -1, &
                1/)
!
    CSIDE=(/ 'WEST  ', &
             'EAST  ', &
             'SOUTH ', &
             'NORTH ', &
             'BOTTOM', &
             'TOP   '/)
    CSIDSE=RESHAPE( (/ 'J','K', &
                       'J','K', &
                       'I','K', &
                       'I','K', &
                       'I','J', &
                       'I','J'/)&
                  ,(/2,6/))

  END SUBROUTINE INITDATA


  SUBROUTINE ALLOCDATA
    IMPLICIT NONE
    INTEGER(IK) :: I,IND1,IND2,NBL1

    ALLOCATE( RM1(NBLOCK),RM2(NBLOCK) )

    ALLOCATE(     U(NXYZA),   V(NXYZA),   W(NXYZA), &
                 F1(NXYZA),  F2(NXYZA),  F3(NXYZA), &
                  P(NXYZA),  PP(NXYZA), VIS(NXYZA),&
               GR1X(NXYZA),GR1Y(NXYZA),GR1Z(NXYZA), &
               GR2X(NXYZA),GR2Y(NXYZA),GR2Z(NXYZA), &
               GR3X(NXYZA),GR3Y(NXYZA),GR3Z(NXYZA), &
                DEN(NXYZA) )

    ALLOCATE(    U_(NXYZA),   V_(NXYZA),  W_(NXYZA), &
                F_1(NXYZA),  F_2(NXYZA), F_3(NXYZA), &
                  Q(NXYZA),   VP(NXYZA), &
                 EX(NXYZA),   EY(NXYZA),  EZ(NXYZA), &
                FXsh(NXYZA),   FYsh(NXYZA),  FZsh(NXYZA),Fmag(NXYZA), &
                EXO(NXYZA),  EYO(NXYZA), EZO(NXYZA), &
              ESP_P(NXYZA),ESP_K(NXYZA))

    ALLOCATE(  FLX_IND(NUM_BF_ALL), &
               FLX_CFI(NUM_CF_ALL), &
                FX_CFI(NUM_CF_ALL), &
                   A_E(NUM_CF_ALL), &
                 X_CFI(NUM_CF_ALL,4), &
                 Y_CFI(NUM_CF_ALL,4), &
                 Z_CFI(NUM_CF_ALL,4), &
              FLX2_CFI(NUM_CF_ALL) )

    ALLOCATE( XC(NXYZA),YC(NXYZA),ZC(NXYZA), &
              FX(NXYZA),FY(NXYZA),FZ(NXYZA), &
             VOL(NXYZA) )

    ALLOCATE(  AE(NXYZA), AW(NXYZA), AN(NXYZA),AS(NXYZA),AT(NXYZA), &
               AB(NXYZA), AP(NXYZA), SU(NXYZA),SV(NXYZA),SW(NXYZA), &
               SP(NXYZA),RES(NXYZA),SUO(NXYZA), &
               AP_ADD(NXYZA), &
              APU(NXYZA),APV(NXYZA),APW(NXYZA), &
              SPU(NXYZA),SPV(NXYZA),SPW(NXYZA) )

      ALLOCATE(ARX_T(NXYZA,3), &
               ARY_T(NXYZA,3), &
               ARZ_T(NXYZA,3), &
               AKX_T(NXYZA,3), &
               AKY_T(NXYZA,3), &
               AKZ_T(NXYZA,3), &
              ARE2_T(NXYZA,3), &
            ARKSI2_T(NXYZA,3), &
               ARE_T(NXYZA,3), &
              ARER_T(NXYZA,3), &
              DELN_T(NXYZA,3), &
             DELNR_T(NXYZA,3) )

      ALLOCATE(LINK_TAB(NBLOCK,NBLOCK))

    IF(LTIME) &
        ALLOCATE(   UO(NXYZA),  VO(NXYZA), WO(NXYZA), &
                   UOO(NXYZA), VOO(NXYZA),WOO(NXYZA), &
                    TO(NXYZA), TOO(NXYZA) , &
                   VPO(NXYZA),  QO(NXYZA), QOO(NXYZA)  )

    IF(LCAL(IEN)) ALLOCATE(T(NXYZA))

    RM1=0. ;RM2=0.

    U=0.    ;    V=0. ;   W=0.
    F1=0.   ;   F2=0. ;  F3=0.
    P=0.    ;   PP=0.
    VIS=0.  ;  DEN=0.
    GR1X=0. ;GR1Y=0.  ;GR1Z=0.
    GR2X=0. ;GR2Y=0.  ;GR2Z=0.
    GR3X=0. ;GR3Y=0.  ;GR3Z=0.

    U_=0.    ;  V_=0. ;  W_=0.
    F_1=0.   ; F_2=0. ; F_3=0.
    EX=0.    ;  EY=0. ;   EZ=0.
    FXsh=0.    ; FYsh=0.  ; FZsh=0.;Fmag=0.
    EXO=0.   ; EYO=0. ; EZO=0.
    ESP_P=0. ;ESP_K=0.

    FLX_IND=0.
    FLX_CFI=0.
    FX_CFI=0.
    A_E=0.
    X_CFI=0.
    Y_CFI=0.
    Z_CFI=0.
    FLX2_CFI=0.

    XC=0. ;YC=0. ;ZC=0.
    FX=0. ;FY=0. ;FZ=0.
    VOL=0.

    AE=0.  ;AW=0.   ;AN=0. ;AS=0. ;AT=0.
    AB=0.  ;AP=0.   ;SU=0. ;SV=0. ;SW=0.
    AP_ADD=0.
    SP=0.  ;RES=0.  ;SUO=0.
    APU=0. ; APV=0. ; APW=0.
    SPU=0. ; SPV=0. ; SPW=0.

    AKX_T=0. ; AKZ_T=0. ; AKY_T=0.
    ARX_T=0. ; ARZ_T=0. ; ARY_T=0.

    IF(LCAL(IEN)) T=0.
    IF(LTIME) THEN
      UO=0. ;  VO=0. ; WO=0. ; UOO=0. ; VOO=0. ;WOO=0. ; TO=0. ;  TOO=0.
      Q=0.     ;QO=0.   ;QOO=0.
      VP=0.    ;VPO=0.
    ENDIF

      LINK_TAB=0
      DO NBL1=1,NBLOCK
       DO I=NUM_SCF(NBL1)+1,NUM_SCF(NBL1)+NUM_CF(NBL1)
        LINK_TAB(NBL1,NUM_CFI(I,5))=LINK_TAB(NBL1,NUM_CFI(I,5))+1
      ENDDO
     ENDDO
    ALLOCATE(BIG_BUF(MAXVAL(LINK_TAB(:,:)),NBLOCK,NBLOCK,2))  ! MAY BE SIZE-OPTIMIZABLE

  END SUBROUTINE ALLOCDATA

  FUNCTION EXTR(A,DIR)
!EXTRACT THE INTERIOR OF THE ARRAY A ON BLOCK NBL
!WITH A SHIFT OF DIR(1) IN Y DIRECTION
!WITH A SHIFT OF DIR(2) IN X DIRECTION
!WITH A SHIFT OF DIR(3) IN Z DIRECTION
    IMPLICIT NONE
    REAL(RK)   ,INTENT(IN)          :: A(NJ,NI,NK)
    INTEGER(IK),INTENT(IN),OPTIONAL :: DIR(3)

    REAL(RK),ALLOCATABLE :: EXTR(:,:,:)
    INTEGER(IK) :: NIE,NJE,NKE,IE,JE,KE

    ALLOCATE(EXTR(NJ,NI,NK))
    EXTR=0.

     NIE= NI-1 ; IE=2
     NJE= NJ-1 ; JE=2 !NO-SHIFT
     NKE= NK-1 ; KE=2

    IF(PRESENT(DIR)) THEN
       NIE= NIE+DIR(2) ; IE=IE+DIR(2)
       NJE= NJE+DIR(1) ; JE=JE+DIR(1) !SHIFT
       NKE= NKE+DIR(3) ; KE=KE+DIR(3)
    ENDIF

    !INTERIOR
    EXTR(2:NJ-1,2:NI-1,2:NK-1)=A(JE:NJE,IE:NIE,KE:NKE)

  END FUNCTION EXTR

  SUBROUTINE DEALLOCDATA
    IMPLICIT NONE

    IF(ALLOCATED(RM1)) THEN ! ONLY DEALLOCATE IF ALLOCATED
      DEALLOCATE(RM1,RM2)


      DEALLOCATE(U,    V,   W, &
                 F1,   F2,  F3, &
                 P,    PP,  &
                 VIS,  DEN, &
                 GR1X,GR1Y,GR1Z, &
                 GR2X,GR2Y,GR2Z, &
                 GR3X,GR3Y,GR3Z )

      DEALLOCATE(  U_,  V_,  W_, &
                   F_1, F_2, F_3, &
                   Q  , VP      , &
                   EX,  EY,   EZ, &
                   FXsh,  FYsh, FZsh,Fmag, &
                   EXO, EYO, EZO, &
                   ESP_P,ESP_K )

      DEALLOCATE(  FLX_IND, &
                   FLX_CFI, &
                   FX_CFI, &
                   A_E, &
                   X_CFI, &
                   Y_CFI, &
                   Z_CFI, &
                   FLX2_CFI )

      DEALLOCATE(XC,YC,ZC, &
                 FX,FY,FZ, &
                 VOL)

      DEALLOCATE(  AE,AW,AN,AS,AT, &
                   AB,AP,SU,SV,SW, &
                   AP_ADD, &
                   SP,RES,SUO, &
                   APU, APV, APW, &
                   SPU, SPV, SPW )

      DEALLOCATE(NUM_SDR,NUM_SCF,NBL_ST)
      DEALLOCATE(NIBL,NJBL,NKBL,NUM_NUM,NUM_SPR)
      DEALLOCATE(NUM_CF,NUM_DR,NUM_TYP)
      DEALLOCATE(NUM_CFI,NUM_IND,NUM_PR)

      DEALLOCATE(X,Y,Z,LI,LK)
      DEALLOCATE(ARX_T,ARY_T,ARZ_T,ARE2_T )
      DEALLOCATE(AKX_T,AKY_T,AKZ_T,ARKSI2_T, &
                 ARE_T,ARER_T,DELN_T,DELNR_T)

      IF(LTIME) &
          DEALLOCATE(   UO,  VO, WO, &
                       UOO, VOO,WOO, &
                        TO, TOO , &
                       VPO,  QO, QOO  )

      IF(LCAL(IEN)) DEALLOCATE(T)

   ENDIF
  END SUBROUTINE DEALLOCDATA

  SUBROUTINE END(MESSAGE)
    IMPLICIT NONE
    CHARACTER(LEN=*),intent(in),OPTIONAL :: MESSAGE

    CALL DEALLOCDATA

    IF(PRESENT(MESSAGE)) PRINT*, MESSAGE
    STOP

  END SUBROUTINE END

SUBROUTINE REALLOCATE_1R(IN,SIZE)
  IMPLICIT NONE
  REAL(RK),ALLOCATABLE,INTENT(INOUT) :: IN(:)
  INTEGER(IK),INTENT(IN)             :: SIZE

  IF(ALLOCATED(IN)) DEALLOCATE(IN)
  ALLOCATE(IN(SIZE))
  
END SUBROUTINE REALLOCATE_1R

SUBROUTINE REALLOCATE_2R(IN,SIZE1,SIZE2)
  IMPLICIT NONE
  REAL(RK),ALLOCATABLE,INTENT(INOUT) :: IN(:,:)
  INTEGER(IK),INTENT(IN)             :: SIZE1,SIZE2

  IF(ALLOCATED(IN)) DEALLOCATE(IN)
  ALLOCATE(IN(SIZE1,SIZE2))
  
END SUBROUTINE REALLOCATE_2R

SUBROUTINE REALLOCATE_3R(IN,SIZE1,SIZE2,SIZE3)
  IMPLICIT NONE
  REAL(RK),ALLOCATABLE,INTENT(INOUT) :: IN(:,:,:)
  INTEGER(IK),INTENT(IN)             :: SIZE1,SIZE2,SIZE3

  IF(ALLOCATED(IN)) DEALLOCATE(IN)
  ALLOCATE(IN(SIZE1,SIZE2,SIZE3))
  
END SUBROUTINE REALLOCATE_3R


SUBROUTINE REALLOCATE_1I(IN,SIZE)
  IMPLICIT NONE
  INTEGER(IK),ALLOCATABLE,INTENT(INOUT) :: IN(:)
  INTEGER(IK),INTENT(IN)                :: SIZE

  IF(ALLOCATED(IN)) DEALLOCATE(IN)
  ALLOCATE(IN(SIZE))
  
END SUBROUTINE REALLOCATE_1I

SUBROUTINE REALLOCATE_2I(IN,SIZE1,SIZE2)
  IMPLICIT NONE
  INTEGER(IK),ALLOCATABLE,INTENT(INOUT) :: IN(:,:)
  INTEGER(IK),INTENT(IN)                :: SIZE1,SIZE2

  IF(ALLOCATED(IN)) DEALLOCATE(IN)
  ALLOCATE(IN(SIZE1,SIZE2))
  
END SUBROUTINE REALLOCATE_2I

SUBROUTINE REALLOCATE_3I(IN,SIZE1,SIZE2,SIZE3)
  IMPLICIT NONE
  INTEGER(IK),ALLOCATABLE,INTENT(INOUT) :: IN(:,:,:)
  INTEGER(IK),INTENT(IN)                :: SIZE1,SIZE2,SIZE3

  IF(ALLOCATED(IN)) DEALLOCATE(IN)
  ALLOCATE(IN(SIZE1,SIZE2,SIZE3))
  
END SUBROUTINE REALLOCATE_3I

  LOGICAL FUNCTION IS_NAN(R)
  IMPLICIT NONE
  REAL(rk) ::r
      IF (R .LT. HUGE(r)) THEN
        IS_NAN = .FALSE.
      ELSE
        IS_NAN = .TRUE.
      END IF
  END FUNCTION IS_NAN
    
   subroutine get_time(time)
      implicit none
      REAL(rk),intent(out) :: time
      logical,parameter :: wall_time=.true.
      integer(kind=8) :: clock_rate,clock
      if (wall_time) then
        call system_clock(clock,clock_rate)
        time=clock*1./clock_rate
      else
        call CPU_TIME(time)
      endif
  end subroutine get_time

!************************************************************
SUBROUTINE SETIND(NBLD)
!************************************************************
!...........................................................
!     THIS ROUTINE SETS THE INDICES FOR THE CURRENT BLOCK.
!----------------------------------------------------------
  IMPLICIT NONE
  INTEGER(IK),INTENT(IN)  :: NBLD
  INTEGER(IK)             :: I,K
!.....
  NI=NIBL(NBLD)
  NJ=NJBL(NBLD)
  NK=NKBL(NBLD)
!.....
  NIM=NI-1
  NJM=NJ-1
  NKM=NK-1
  NIJ=NI*NJ
  NIK=NI*NK
  NJK=NJ*NK
  NIJK=NI*NJ*NK
!
  LI(1:NI)=(/((I-1)*NJ+NBL_ST(NBLD),I=1,NI)/)
  LK(1:NK)=(/((K-1)*NIJ,K=1,NK)/)
!
  ICST= NBL_ST(NBLD)+1
  ICEN= NBL_ST(NBLD)+NIJK
!
  ID_SIDE(1)=NJ
  ID_SIDE(2)=1
  ID_SIDE(3)=NIJ

END SUBROUTINE SETIND
END MODULE DEC

