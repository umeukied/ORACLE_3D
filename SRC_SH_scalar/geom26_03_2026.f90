!*************************************************************
!> Geometry metrics used at boundaries, block interfaces
!> and inside the computational domain.
!>
!> Main module variables
!> ---------------------
!> INE         : index of neighbor point/cell across the face
!> ARX,ARY,ARZ : oriented face area vector components
!> AKX,AKY,AKZ : center-to-center vector components
!> ARE         : face area
!> ARE2        : face area squared
!> ARKSI2      : squared norm of center-to-center vector
!> DELN        : projection of center-to-center vector on face normal
!> DELNR       : inverse of DELN
!> FXE,FXW     : interpolation factors (declared here, computed elsewhere)
!>
!> Notes
!> -----
!> - This version preserves the original variables and intent.
!> - The internal structure has been clarified and factorized.
!> - Numerical protections have been harmonized.
!*************************************************************
MODULE GEOM
  USE DEC
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: INIT_GEOM_BF
  PUBLIC :: INIT_GEOM_BLOCK
  PUBLIC :: INIT_GEOM_T

  !-----------------------------------------------------------------
  ! Common variables preserved from the original implementation
  !-----------------------------------------------------------------
  REAL(RK)    :: ARX, ARY, ARZ
  REAL(RK)    :: AKX, AKY, AKZ
  REAL(RK)    :: ARE, ARE2, ARKSI2
  REAL(RK)    :: FXE, FXW
  INTEGER(IK) :: INE
  REAL(RK)    :: DELN, DELNR

CONTAINS

  !===========================================================
  !> Return the local directional offsets associated with a face
  !>
  !> DIR   : main direction
  !> DIR_P : current location
  !> DIR_E : forward neighbor in face-normal direction
  !> DIR_S : first transverse direction
  !> DIR_B : second transverse direction
  !> DIR_BS: diagonal transverse direction
  !===========================================================
  PURE SUBROUTINE GET_LOCAL_DIRS(DIR, DIR_P, DIR_E, DIR_S, DIR_B, DIR_BS)
    IMPLICIT NONE
    INTEGER(IK), INTENT(IN)  :: DIR(3)
    INTEGER(IK), INTENT(OUT) :: DIR_P(3), DIR_E(3), DIR_S(3), DIR_B(3), DIR_BS(3)

    DIR_P  = (/0_IK, 0_IK, 0_IK/)
    DIR_E  = DIR
    DIR_S  = -(/DIR(2), DIR(3), DIR(1)/)
    DIR_B  = -(/DIR(3), DIR(1), DIR(2)/)
    DIR_BS = DIR_S + DIR_B
  END SUBROUTINE GET_LOCAL_DIRS


  !===========================================================
  !> Compute cross product C = A x B
  !===========================================================
  PURE SUBROUTINE CROSS_PRODUCT(AX, AY, AZ, BX, BY, BZ, CX, CY, CZ)
    IMPLICIT NONE
    REAL(RK), INTENT(IN)  :: AX, AY, AZ
    REAL(RK), INTENT(IN)  :: BX, BY, BZ
    REAL(RK), INTENT(OUT) :: CX, CY, CZ

    CX = AY*BZ - BY*AZ
    CY = AZ*BX - BZ*AX
    CZ = AX*BY - BX*AY
  END SUBROUTINE CROSS_PRODUCT


  !===========================================================
  !> Compute scalar face metrics from face area vector and
  !> center-to-center vector.
  !===========================================================
  PURE SUBROUTINE COMPUTE_SCALAR_METRICS(ARX_L, ARY_L, ARZ_L, &
                                         AKX_L, AKY_L, AKZ_L, &
                                         ARE2_L, ARE_L, ARKSI2_L, DELN_L, DELNR_L)
    IMPLICIT NONE
    REAL(RK), INTENT(IN)  :: ARX_L, ARY_L, ARZ_L
    REAL(RK), INTENT(IN)  :: AKX_L, AKY_L, AKZ_L
    REAL(RK), INTENT(OUT) :: ARE2_L, ARE_L, ARKSI2_L, DELN_L, DELNR_L

    REAL(RK) :: ARER_L

    ARE2_L   = MAX(ARX_L*ARX_L + ARY_L*ARY_L + ARZ_L*ARZ_L, 0._RK)
    ARE_L    = SQRT(ARE2_L)
    ARER_L   = 1._RK/(ARE_L + SMALL)

    ARKSI2_L = MAX(AKX_L*AKX_L + AKY_L*AKY_L + AKZ_L*AKZ_L, 0._RK)

    DELN_L   = (ARX_L*AKX_L + ARY_L*AKY_L + ARZ_L*AKZ_L)*ARER_L
    DELNR_L  = 1._RK/(DELN_L + SMALL)
  END SUBROUTINE COMPUTE_SCALAR_METRICS


  !===========================================================
  !> Extract geometric info at a boundary face
  !>
  !> Fill the common module variables:
  !> INE, ARX, ARY, ARZ, AKX, AKY, AKZ, ARE, ARE2,
  !> ARKSI2, DELN, DELNR
  !===========================================================
  SUBROUTINE INIT_GEOM_BF(ISIDE, INP)
    IMPLICIT NONE

    INTEGER(IK), INTENT(IN) :: ISIDE, INP

    REAL(RK)    :: FACS
    REAL(RK)    :: DXS, DYS, DZS
    REAL(RK)    :: DXT, DYT, DZT
    INTEGER(IK) :: INS, INB, INBS, INPG
    INTEGER(IK) :: IDEW, IDNS, IDTB

    !---------------------------------------------------------
    ! 1) Local topological directions
    !---------------------------------------------------------
    FACS = REAL(ID_SIGN(ISIDE), RK)
    IDEW = ID_SIDE(ID_ODER(ISIDE,1))*ID_SIGN(ISIDE)
    IDNS = ID_SIDE(ID_ODER(ISIDE,2))
    IDTB = ID_SIDE(ID_ODER(ISIDE,3))

    !---------------------------------------------------------
    ! 2) Neighboring indices around the boundary face
    !---------------------------------------------------------
    INE  = INP + IDEW
    INPG = INP + MIN(IDEW, 0_IK)
    INS  = INPG - IDNS
    INBS = INS  - IDTB
    INB  = INPG - IDTB

    !---------------------------------------------------------
    ! 3) Two tangent vectors spanning the face
    !---------------------------------------------------------
    DXS = 0.5_RK*(X(INPG) - X(INS) + X(INB) - X(INBS))
    DYS = 0.5_RK*(Y(INPG) - Y(INS) + Y(INB) - Y(INBS))
    DZS = 0.5_RK*(Z(INPG) - Z(INS) + Z(INB) - Z(INBS))

    DXT = 0.5_RK*(X(INPG) - X(INB) + X(INS) - X(INBS))
    DYT = 0.5_RK*(Y(INPG) - Y(INB) + Y(INS) - Y(INBS))
    DZT = 0.5_RK*(Z(INPG) - Z(INB) + Z(INS) - Z(INBS))

    !---------------------------------------------------------
    ! 4) Oriented face area vector
    !---------------------------------------------------------
    CALL CROSS_PRODUCT(DXS, DYS, DZS, DXT, DYT, DZT, ARX, ARY, ARZ)
    ARX = ARX*FACS
    ARY = ARY*FACS
    ARZ = ARZ*FACS

    !---------------------------------------------------------
    ! 5) Center-to-center vector
    !---------------------------------------------------------
    AKX = XC(INE) - XC(INP)
    AKY = YC(INE) - YC(INP)
    AKZ = ZC(INE) - ZC(INP)

    !---------------------------------------------------------
    ! 6) Derived scalar metrics
    !---------------------------------------------------------
    CALL COMPUTE_SCALAR_METRICS(ARX, ARY, ARZ, AKX, AKY, AKZ, &
                                ARE2, ARE, ARKSI2, DELN, DELNR)

  END SUBROUTINE INIT_GEOM_BF


  !===========================================================
  !> Extract geometric info at an interface between two blocks
  !>
  !> Fill the common module variables:
  !> INE, ARX, ARY, ARZ, ARE2
  !===========================================================
  SUBROUTINE INIT_GEOM_BLOCK(ISIDE, I)
    IMPLICIT NONE

    INTEGER(IK), INTENT(IN) :: ISIDE, I

    REAL(RK) :: FACS
    REAL(RK) :: DXS, DYS, DZS
    REAL(RK) :: DXT, DYT, DZT

    FACS = REAL(ID_SIGN(ISIDE), RK)

    !---------------------------------------------------------
    ! Neighbor point/cell at the other side of the block interface
    !---------------------------------------------------------
    INE = NBL_ST(NUM_CFI(I,5)) + NUM_CFI(I,4)

    !---------------------------------------------------------
    ! Two tangent vectors spanning the interface face
    !---------------------------------------------------------
    DXS = 0.5_RK*(X_CFI(I,1) - X_CFI(I,2) + X_CFI(I,4) - X_CFI(I,3))
    DYS = 0.5_RK*(Y_CFI(I,1) - Y_CFI(I,2) + Y_CFI(I,4) - Y_CFI(I,3))
    DZS = 0.5_RK*(Z_CFI(I,1) - Z_CFI(I,2) + Z_CFI(I,4) - Z_CFI(I,3))

    DXT = 0.5_RK*(X_CFI(I,1) - X_CFI(I,4) + X_CFI(I,2) - X_CFI(I,3))
    DYT = 0.5_RK*(Y_CFI(I,1) - Y_CFI(I,4) + Y_CFI(I,2) - Y_CFI(I,3))
    DZT = 0.5_RK*(Z_CFI(I,1) - Z_CFI(I,4) + Z_CFI(I,2) - Z_CFI(I,3))

    !---------------------------------------------------------
    ! Oriented face area vector
    !---------------------------------------------------------
    CALL CROSS_PRODUCT(DXS, DYS, DZS, DXT, DYT, DZT, ARX, ARY, ARZ)
    ARX = ARX*FACS
    ARY = ARY*FACS
    ARZ = ARZ*FACS

    !---------------------------------------------------------
    ! Squared face area
    !---------------------------------------------------------
    ARE2 = MAX(ARX*ARX + ARY*ARY + ARZ*ARZ, 0._RK)

  END SUBROUTINE INIT_GEOM_BLOCK


  !===========================================================
  !> Initialize geometric information in the whole domain
  !>
  !> Inputs:
  !>   X,Y,Z    : corner coordinates
  !>   XC,YC,ZC : cell-center coordinates
  !>
  !> Outputs:
  !>   ARX_T,ARY_T,ARZ_T : oriented face area vector components
  !>   AKX_T,AKY_T,AKZ_T : center-to-center vector components
  !>   ARE2_T,ARE_T      : face area squared / face area
  !>   ARER_T            : inverse of face area
  !>   ARKSI2_T          : squared norm of center-to-center vector
  !>   DELN_T            : normal projected center-to-center distance
  !>   DELNR_T           : inverse of DELN_T
  !===========================================================
  SUBROUTINE INIT_GEOM_T(X, Y, Z, XC, YC, ZC, &
                         ARX_T, ARY_T, ARZ_T, &
                         AKX_T, AKY_T, AKZ_T, &
                         ARE2_T, ARKSI2_T, &
                         ARE_T, ARER_T, DELN_T, DELNR_T)
    IMPLICIT NONE

    REAL(RK), INTENT(IN)  :: X(NXYZA), Y(NXYZA), Z(NXYZA)
    REAL(RK), INTENT(IN)  :: XC(NXYZA), YC(NXYZA), ZC(NXYZA)

    REAL(RK), INTENT(OUT) :: ARX_T(NXYZA,3), ARY_T(NXYZA,3), ARZ_T(NXYZA,3)
    REAL(RK), INTENT(OUT) :: AKX_T(NXYZA,3), AKY_T(NXYZA,3), AKZ_T(NXYZA,3)
    REAL(RK), INTENT(OUT) :: ARE2_T(NXYZA,3), ARKSI2_T(NXYZA,3)
    REAL(RK), INTENT(OUT) :: ARE_T(NXYZA,3), ARER_T(NXYZA,3)
    REAL(RK), INTENT(OUT) :: DELN_T(NXYZA,3), DELNR_T(NXYZA,3)

    REAL(RK), ALLOCATABLE :: TMPX(:,:,:), TMPY(:,:,:), TMPZ(:,:,:)
    REAL(RK), ALLOCATABLE :: DXS(:,:,:), DYS(:,:,:), DZS(:,:,:)
    REAL(RK), ALLOCATABLE :: DXT(:,:,:), DYT(:,:,:), DZT(:,:,:)

    INTEGER(IK) :: NBL, D
    INTEGER(IK) :: DIR(3)
    INTEGER(IK) :: DIR_P(3), DIR_E(3), DIR_S(3), DIR_B(3), DIR_BS(3)

    !---------------------------------------------------------
    ! Initialize outputs
    !---------------------------------------------------------
    ARX_T    = 0._RK
    ARY_T    = 0._RK
    ARZ_T    = 0._RK
    AKX_T    = 0._RK
    AKY_T    = 0._RK
    AKZ_T    = 0._RK
    ARE2_T   = 0._RK
    ARKSI2_T = 0._RK
    ARE_T    = 0._RK
    ARER_T   = 0._RK
    DELN_T   = 0._RK
    DELNR_T  = 0._RK

    !---------------------------------------------------------
    ! Loop over mesh blocks
    !---------------------------------------------------------
    DO NBL = 1, NBLOCK

      CALL SETIND(NBL)

      ALLOCATE(TMPX(NJ,NI,NK), TMPY(NJ,NI,NK), TMPZ(NJ,NI,NK))
      ALLOCATE(DXS (NJ,NI,NK), DYS (NJ,NI,NK), DZS (NJ,NI,NK))
      ALLOCATE(DXT (NJ,NI,NK), DYT (NJ,NI,NK), DZT (NJ,NI,NK))

      !=======================================================
      ! 1) Face area vectors from corner coordinates
      !=======================================================
      TMPX = RESHAPE(X(ICST:ICEN), [NJ,NI,NK])
      TMPY = RESHAPE(Y(ICST:ICEN), [NJ,NI,NK])
      TMPZ = RESHAPE(Z(ICST:ICEN), [NJ,NI,NK])

      DIR = (/0_IK, 0_IK, 1_IK/)

      DO D = 1, 3
        DIR = CSHIFT(DIR, -1)
        CALL GET_LOCAL_DIRS(DIR, DIR_P, DIR_E, DIR_S, DIR_B, DIR_BS)

        DXS = 0.5_RK*(EXTR(TMPX,DIR_P) - EXTR(TMPX,DIR_S) + EXTR(TMPX,DIR_B) - EXTR(TMPX,DIR_BS))
        DYS = 0.5_RK*(EXTR(TMPY,DIR_P) - EXTR(TMPY,DIR_S) + EXTR(TMPY,DIR_B) - EXTR(TMPY,DIR_BS))
        DZS = 0.5_RK*(EXTR(TMPZ,DIR_P) - EXTR(TMPZ,DIR_S) + EXTR(TMPZ,DIR_B) - EXTR(TMPZ,DIR_BS))

        DXT = 0.5_RK*(EXTR(TMPX,DIR_P) + EXTR(TMPX,DIR_S) - EXTR(TMPX,DIR_B) - EXTR(TMPX,DIR_BS))
        DYT = 0.5_RK*(EXTR(TMPY,DIR_P) + EXTR(TMPY,DIR_S) - EXTR(TMPY,DIR_B) - EXTR(TMPY,DIR_BS))
        DZT = 0.5_RK*(EXTR(TMPZ,DIR_P) + EXTR(TMPZ,DIR_S) - EXTR(TMPZ,DIR_B) - EXTR(TMPZ,DIR_BS))

        ARX_T(ICST:ICEN,D) = RESHAPE(DYS*DZT - DYT*DZS, [NIJK])
        ARY_T(ICST:ICEN,D) = RESHAPE(DZS*DXT - DZT*DXS, [NIJK])
        ARZ_T(ICST:ICEN,D) = RESHAPE(DXS*DYT - DXT*DYS, [NIJK])
      END DO

      !=======================================================
      ! 2) Center-to-center vectors from cell-center coordinates
      !=======================================================
      TMPX = RESHAPE(XC(ICST:ICEN), [NJ,NI,NK])
      TMPY = RESHAPE(YC(ICST:ICEN), [NJ,NI,NK])
      TMPZ = RESHAPE(ZC(ICST:ICEN), [NJ,NI,NK])

      DIR = (/0_IK, 0_IK, 1_IK/)

      DO D = 1, 3
        DIR = CSHIFT(DIR, -1)
        CALL GET_LOCAL_DIRS(DIR, DIR_P, DIR_E, DIR_S, DIR_B, DIR_BS)

        AKX_T(ICST:ICEN,D) = RESHAPE(EXTR(TMPX,DIR_E) - EXTR(TMPX,DIR_P), [NIJK])
        AKY_T(ICST:ICEN,D) = RESHAPE(EXTR(TMPY,DIR_E) - EXTR(TMPY,DIR_P), [NIJK])
        AKZ_T(ICST:ICEN,D) = RESHAPE(EXTR(TMPZ,DIR_E) - EXTR(TMPZ,DIR_P), [NIJK])
      END DO

      DEALLOCATE(TMPX, TMPY, TMPZ)
      DEALLOCATE(DXS, DYS, DZS, DXT, DYT, DZT)

    END DO

    !---------------------------------------------------------
    ! 3) Derived scalar metrics over the whole domain
    !---------------------------------------------------------
    ARE2_T   = MAX(ARX_T*ARX_T + ARY_T*ARY_T + ARZ_T*ARZ_T, 0._RK)
    ARE_T    = SQRT(ARE2_T)
    ARER_T   = 1._RK/(ARE_T + SMALL)

    ARKSI2_T = MAX(AKX_T*AKX_T + AKY_T*AKY_T + AKZ_T*AKZ_T, 0._RK)

    DELN_T   = (ARX_T*AKX_T + ARY_T*AKY_T + ARZ_T*AKZ_T)*ARER_T
    DELNR_T  = 1._RK/(DELN_T + SMALL)

  END SUBROUTINE INIT_GEOM_T

END MODULE GEOM