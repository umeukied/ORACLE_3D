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
  INTEGER, PARAMETER :: INT4 = SELECTED_INT_KIND(8)
  INTEGER, PARAMETER :: INT8 = SELECTED_INT_KIND(10)

 ! FLOATING POINT KINDS
  INTEGER, PARAMETER :: REAL4  = SELECTED_REAL_KIND(6)
  INTEGER, PARAMETER :: REAL8  = SELECTED_REAL_KIND(15)
  INTEGER, PARAMETER :: REAL16 = SELECTED_REAL_KIND(32)

  INTEGER, PARAMETER :: IK = INT4
  INTEGER, PARAMETER :: RK = REAL8

  INTEGER(IK),PARAMETER  :: NPHI=5,NPHI2=7
  REAL(RK)               :: MSTABILITY, L_PERMIT, N_MOBILITY
  INTEGER(IK)            :: NXMAX,NZMAX,NUM_BF_ALL,NUM_CF_ALL,NXYZA,NUM_DR_ALL=35

  REAL(RK), PARAMETER :: PI = 4.0_RK*atan(1.0_RK)

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
                           RE, RE_INPUT, &
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
                   LESP, &
                   LVORTEX        ! T=calcul vorticite+Q-critere+Lambda2, F=rien

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
                               FXsh(:), FYsh(:), FZsh(:), Fmag(:), &
                               EXO(:), EYO(:), EZO(:), &
                               ESP_P(:), ESP_K(:), &
                               SEHD_C(:)

  REAL(RK)             ::  REY, D_CHARGE, R, &
                           CINJECTION,TSTABILITY, &
                           LAMBDA, LAMBDAR, &
                           D_C, &
                           F_BASE, VPmax, &
                           MuLOCATION, sigma, sigmaR, &
                           Qmax

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

  INTEGER(IK)       :: NDREG

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
                           AKX_T(:,:),   AKY_T(:,:), AKZ_T(:,:), &
                           ARE2_T(:,:), ARKSI2_T(:,:), ARE_T(:,:), &
                           DELN_T(:,:), DELNR_T(:,:), ARER_T(:,:)

!BPCHAR
  CHARACTER(LEN=6) :: CSIDE(6)
  CHARACTER(LEN=1) :: CSIDSE(2,6)

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * *      INC_COE    * * * * * * * * * * * * * * * *
!COEF
  REAL(RK),ALLOCATABLE ::   AE(:), AW(:), AN(:),AS(:),AT(:), &
                            AB(:), AP(:), SU(:),SV(:),SW(:), &
                            AP_ADD(:), &
                            SP(:),RES(:),SUO(:), &
                            APU(:),APV(:),APW(:), &
                            SPU(:),SPV(:),SPW(:),BIG_BUF(:,:,:,:)

  INTEGER(IK),ALLOCATABLE :: LINK_TAB(:,:)

CONTAINS

!**********************
  SUBROUTINE INITDATA
!**********************
    CHVAR =(/' U ',' V ',' W ',' P ',' EN',' Q ',' VP'/)

    IU=1
    IV=2
    IW=3
    IP=4
    IEN=5
    IQ=6
    IVP=7

    ID_ODER=RESHAPE( (/ 1, 1, 2, 2, 3, 3, &
                        2, 2, 3, 3, 1, 1, &
                        3, 3, 1, 1, 2, 2/), (/6,3/) )

    ID_SIGN=(/ -1, 1, -1, 1, -1, 1 /)

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
                       'I','J'/), (/2,6/) )

  END SUBROUTINE INITDATA

  SUBROUTINE ALLOCDATA
    IMPLICIT NONE
    INTEGER(IK) :: I,IND1,IND2,NBL1

    ALLOCATE(RM1(NBLOCK),RM2(NBLOCK))

    ALLOCATE( U(NXYZA), V(NXYZA), W(NXYZA), &
              F1(NXYZA), F2(NXYZA), F3(NXYZA), &
              P(NXYZA), PP(NXYZA), VIS(NXYZA), &
              GR1X(NXYZA), GR1Y(NXYZA), GR1Z(NXYZA), &
              GR2X(NXYZA), GR2Y(NXYZA), GR2Z(NXYZA), &
              GR3X(NXYZA), GR3Y(NXYZA), GR3Z(NXYZA), &
              DEN(NXYZA) )

    ALLOCATE( U_(NXYZA), V_(NXYZA), W_(NXYZA), &
              F_1(NXYZA), F_2(NXYZA), F_3(NXYZA), &
              Q(NXYZA), VP(NXYZA), &
              EX(NXYZA), EY(NXYZA), EZ(NXYZA), &
              FXsh(NXYZA), FYsh(NXYZA), FZsh(NXYZA), Fmag(NXYZA), &
              EXO(NXYZA), EYO(NXYZA), EZO(NXYZA), &
              ESP_P(NXYZA), ESP_K(NXYZA), &
              SEHD_C(NXYZA) )

    ALLOCATE( FLX_IND(NUM_BF_ALL), &
              FLX_CFI(NUM_CF_ALL), &
              FX_CFI(NUM_CF_ALL), &
              A_E(NUM_CF_ALL), &
              X_CFI(NUM_CF_ALL,4), &
              Y_CFI(NUM_CF_ALL,4), &
              Z_CFI(NUM_CF_ALL,4), &
              FLX2_CFI(NUM_CF_ALL) )

    ALLOCATE( XC(NXYZA), YC(NXYZA), ZC(NXYZA), &
              FX(NXYZA), FY(NXYZA), FZ(NXYZA), &
              VOL(NXYZA) )

    ALLOCATE( AE(NXYZA), AW(NXYZA), AN(NXYZA), AS(NXYZA), AT(NXYZA), &
              AB(NXYZA), AP(NXYZA), SU(NXYZA), SV(NXYZA), SW(NXYZA), &
              SP(NXYZA), RES(NXYZA), SUO(NXYZA), &
              AP_ADD(NXYZA), &
              APU(NXYZA), APV(NXYZA), APW(NXYZA), &
              SPU(NXYZA), SPV(NXYZA), SPW(NXYZA) )

    ALLOCATE( ARX_T(NXYZA,3), &
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

    IF(LTIME) THEN
       ALLOCATE( UO(NXYZA), VO(NXYZA), WO(NXYZA), &
                 UOO(NXYZA), VOO(NXYZA), WOO(NXYZA), &
                 TO(NXYZA), TOO(NXYZA), &
                 VPO(NXYZA), QO(NXYZA), QOO(NXYZA) )
    END IF

    IF(LCAL(IEN)) ALLOCATE(T(NXYZA))

    RM1 = 0.0_RK
    RM2 = 0.0_RK

    U    = 0.0_RK
    V    = 0.0_RK
    W    = 0.0_RK
    F1   = 0.0_RK
    F2   = 0.0_RK
    F3   = 0.0_RK
    P    = 0.0_RK
    PP   = 0.0_RK
    VIS  = 0.0_RK
    DEN  = 0.0_RK

    GR1X = 0.0_RK
    GR1Y = 0.0_RK
    GR1Z = 0.0_RK
    GR2X = 0.0_RK
    GR2Y = 0.0_RK
    GR2Z = 0.0_RK
    GR3X = 0.0_RK
    GR3Y = 0.0_RK
    GR3Z = 0.0_RK

    U_     = 0.0_RK
    V_     = 0.0_RK
    W_     = 0.0_RK
    F_1    = 0.0_RK
    F_2    = 0.0_RK
    F_3    = 0.0_RK
    Q      = 0.0_RK
    VP     = 0.0_RK
    EX     = 0.0_RK
    EY     = 0.0_RK
    EZ     = 0.0_RK
    FXsh   = 0.0_RK
    FYsh   = 0.0_RK
    FZsh   = 0.0_RK
    Fmag   = 0.0_RK
    EXO    = 0.0_RK
    EYO    = 0.0_RK
    EZO    = 0.0_RK
    ESP_P  = 0.0_RK
    ESP_K  = 0.0_RK
    SEHD_C = 0.0_RK

    FLX_IND  = 0.0_RK
    FLX_CFI  = 0.0_RK
    FX_CFI   = 0.0_RK
    A_E      = 0.0_RK
    X_CFI    = 0.0_RK
    Y_CFI    = 0.0_RK
    Z_CFI    = 0.0_RK
    FLX2_CFI = 0.0_RK

    XC  = 0.0_RK
    YC  = 0.0_RK
    ZC  = 0.0_RK
    FX  = 0.0_RK
    FY  = 0.0_RK
    FZ  = 0.0_RK
    VOL = 0.0_RK

    AE     = 0.0_RK
    AW     = 0.0_RK
    AN     = 0.0_RK
    AS     = 0.0_RK
    AT     = 0.0_RK
    AB     = 0.0_RK
    AP     = 0.0_RK
    SU     = 0.0_RK
    SV     = 0.0_RK
    SW     = 0.0_RK
    AP_ADD = 0.0_RK
    SP     = 0.0_RK
    RES    = 0.0_RK
    SUO    = 0.0_RK
    APU    = 0.0_RK
    APV    = 0.0_RK
    APW    = 0.0_RK
    SPU    = 0.0_RK
    SPV    = 0.0_RK
    SPW    = 0.0_RK

    AKX_T   = 0.0_RK
    AKY_T   = 0.0_RK
    AKZ_T   = 0.0_RK
    ARX_T   = 0.0_RK
    ARY_T   = 0.0_RK
    ARZ_T   = 0.0_RK
    ARE2_T  = 0.0_RK
    ARKSI2_T= 0.0_RK
    ARE_T   = 0.0_RK
    ARER_T  = 0.0_RK
    DELN_T  = 0.0_RK
    DELNR_T = 0.0_RK

    IF(LCAL(IEN)) T = 0.0_RK

    IF(LTIME) THEN
      UO  = 0.0_RK
      VO  = 0.0_RK
      WO  = 0.0_RK
      UOO = 0.0_RK
      VOO = 0.0_RK
      WOO = 0.0_RK
      TO  = 0.0_RK
      TOO = 0.0_RK
      QO  = 0.0_RK
      QOO = 0.0_RK
      VPO = 0.0_RK
    END IF

    LINK_TAB = 0

    DO NBL1=1,NBLOCK
       DO I=NUM_SCF(NBL1)+1,NUM_SCF(NBL1)+NUM_CF(NBL1)
          LINK_TAB(NBL1,NUM_CFI(I,5)) = LINK_TAB(NBL1,NUM_CFI(I,5)) + 1
       END DO
    END DO

    ALLOCATE(BIG_BUF(MAXVAL(LINK_TAB(:,:)),NBLOCK,NBLOCK,2))

  END SUBROUTINE ALLOCDATA

  FUNCTION EXTR(A,DIR)
    IMPLICIT NONE
    REAL(RK),INTENT(IN)          :: A(NJ,NI,NK)
    INTEGER(IK),INTENT(IN),OPTIONAL :: DIR(3)

    REAL(RK),ALLOCATABLE :: EXTR(:,:,:)
    INTEGER(IK) :: NIE,NJE,NKE,IE,JE,KE

    ALLOCATE(EXTR(NJ,NI,NK))
    EXTR=0.0_RK

    NIE=NI-1 ; IE=2
    NJE=NJ-1 ; JE=2
    NKE=NK-1 ; KE=2

    IF(PRESENT(DIR)) THEN
       NIE=NIE+DIR(2) ; IE=IE+DIR(2)
       NJE=NJE+DIR(1) ; JE=JE+DIR(1)
       NKE=NKE+DIR(3) ; KE=KE+DIR(3)
    END IF

    EXTR(2:NJ-1,2:NI-1,2:NK-1)=A(JE:NJE,IE:NIE,KE:NKE)

  END FUNCTION EXTR

  SUBROUTINE DEALLOCDATA
    IMPLICIT NONE

    IF(ALLOCATED(RM1)) THEN

      DEALLOCATE(RM1,RM2)

      DEALLOCATE(U,V,W, &
                 F1,F2,F3, &
                 P,PP, &
                 VIS,DEN, &
                 GR1X,GR1Y,GR1Z, &
                 GR2X,GR2Y,GR2Z, &
                 GR3X,GR3Y,GR3Z)

      DEALLOCATE(U_,V_,W_, &
                 F_1,F_2,F_3, &
                 Q,VP, &
                 EX,EY,EZ, &
                 FXsh,FYsh,FZsh,Fmag, &
                 EXO,EYO,EZO, &
                 ESP_P,ESP_K, &
                 SEHD_C)

      DEALLOCATE(FLX_IND, &
                 FLX_CFI, &
                 FX_CFI, &
                 A_E, &
                 X_CFI, &
                 Y_CFI, &
                 Z_CFI, &
                 FLX2_CFI)

      DEALLOCATE(XC,YC,ZC, &
                 FX,FY,FZ, &
                 VOL)

      DEALLOCATE(AE,AW,AN,AS,AT, &
                 AB,AP,SU,SV,SW, &
                 AP_ADD, &
                 SP,RES,SUO, &
                 APU,APV,APW, &
                 SPU,SPV,SPW)

      DEALLOCATE(NUM_SDR,NUM_SCF,NBL_ST)
      DEALLOCATE(NIBL,NJBL,NKBL,NUM_NUM,NUM_SPR)
      DEALLOCATE(NUM_CF,NUM_DR,NUM_TYP)
      DEALLOCATE(NUM_CFI,NUM_IND,NUM_PR)

      DEALLOCATE(X,Y,Z,LI,LK)
      DEALLOCATE(ARX_T,ARY_T,ARZ_T,ARE2_T)
      DEALLOCATE(AKX_T,AKY_T,AKZ_T,ARKSI2_T, &
                 ARE_T,ARER_T,DELN_T,DELNR_T)

      IF(LTIME) THEN
         DEALLOCATE(UO,VO,WO, &
                    UOO,VOO,WOO, &
                    TO,TOO, &
                    VPO,QO,QOO)
      END IF

      IF(LCAL(IEN)) DEALLOCATE(T)

    END IF
  END SUBROUTINE DEALLOCDATA

  SUBROUTINE END(MESSAGE)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: MESSAGE

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
    REAL(RK),INTENT(IN) :: R
    IS_NAN = (R /= R)
  END FUNCTION IS_NAN

  SUBROUTINE GET_TIME(TIME_OUT)
    IMPLICIT NONE
    REAL(RK),INTENT(OUT) :: TIME_OUT
    LOGICAL,PARAMETER :: WALL_TIME=.TRUE.
    INTEGER(KIND=8) :: CLOCK_RATE,CLOCK

    IF (WALL_TIME) THEN
      CALL SYSTEM_CLOCK(CLOCK,CLOCK_RATE)
      TIME_OUT = REAL(CLOCK,RK)/REAL(CLOCK_RATE,RK)
    ELSE
      CALL CPU_TIME(TIME_OUT)
    END IF
  END SUBROUTINE GET_TIME

!************************************************************
  SUBROUTINE SETIND(NBLD)
!************************************************************
    IMPLICIT NONE
    INTEGER(IK),INTENT(IN)  :: NBLD
    INTEGER(IK)             :: I,K

    NI=NIBL(NBLD)
    NJ=NJBL(NBLD)
    NK=NKBL(NBLD)

    NIM=NI-1
    NJM=NJ-1
    NKM=NK-1
    NIJ=NI*NJ
    NIK=NI*NK
    NJK=NJ*NK
    NIJK=NI*NJ*NK

    LI(1:NI)=(/((I-1)*NJ+NBL_ST(NBLD),I=1,NI)/)
    LK(1:NK)=(/((K-1)*NIJ,K=1,NK)/)

    ICST= NBL_ST(NBLD)+1
    ICEN= NBL_ST(NBLD)+NIJK

    ID_SIDE(1)=NJ
    ID_SIDE(2)=1
    ID_SIDE(3)=NIJ

  END SUBROUTINE SETIND

END MODULE DEC
