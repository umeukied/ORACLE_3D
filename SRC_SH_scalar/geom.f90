!*************************************************************
!> Declare and compute common variables computed at boundaries and interfaces
!>
!> Variables   | Signification
!> ------------|-----------------
!> INE         | Point number on the boundary
!> ARX,ARY,ARZ | Projected area in the X Y and Z direction
!> AKX,AKY,AKZ | Components of normal vector
!> ARE         | Area
!> ARE2        | Area²
!> ARKSI2      | ?
!> DELN,DELNR  | ?
!> FXE,FXW     | Interpolation factor
!>
!> @author Philippe Traore
!> @todo Complete documentation
!> @todo Put everything in INIT_GEOM_T
!*************************************************************
!

MODULE GEOM
  USE DEC
  IMPLICIT NONE

  REAL(RK)    :: ARX,ARY,ARZ
  REAL(RK)    :: AKX,AKY,AKZ
  REAL(RK)    :: ARE,ARE2,ARKSI2
  REAL(RK)    :: FXE,FXW
  INTEGER(IK) :: INE
  REAL(RK)    :: DELN,DELNR 

CONTAINS

!> Extract geometric info at a boundary
!>
!> Fill the common variables of the module
!>
!> @param ISIDE Boundary direction
!> @param I Index along the boundaries
!> @todo Put in INIT_GEOM_T
!> @author Philippe Traore
  SUBROUTINE INIT_GEOM_BF(ISIDE,INP)
    IMPLICIT NONE
    INTEGER(IK),INTENT(IN) :: ISIDE,INP

    REAL(RK)               :: FACS,ARER
    REAL(RK)               :: DXS,DYS,DZS
    REAL(RK)               :: DXT,DYT,DZT
    INTEGER(IK)            :: INS,INB,INBS,INPG
    INTEGER(IK)            :: IDEW,IDNS,IDTB

    FACS=REAL(ID_SIGN(ISIDE),rk)
    IDEW=ID_SIDE  (ID_ODER(ISIDE,1))*ID_SIGN(ISIDE)
    IDNS=ID_SIDE  (ID_ODER(ISIDE,2))
    IDTB=ID_SIDE  (ID_ODER(ISIDE,3))
!
!.....CALCULATE INDEX (POINT ON CF-BOUNDARY)
    INE=INP+IDEW
!
    INPG=INP +MIN(IDEW,0_ik)
    INS =INPG-IDNS
    INBS=INS -IDTB
    INB =INPG-IDTB
!
!.....SECOND   WITH    INP = > INPG
    DXS=.5*(X(INPG)-X(INS)+X(INB)-X(INBS))
    DYS=.5*(Y(INPG)-Y(INS)+Y(INB)-Y(INBS))
    DZS=.5*(Z(INPG)-Z(INS)+Z(INB)-Z(INBS))
!.....THIRD
    DXT=.5*(X(INPG)-X(INB)+X(INS)-X(INBS))
    DYT=.5*(Y(INPG)-Y(INB)+Y(INS)-Y(INBS))
    DZT=.5*(Z(INPG)-Z(INB)+Z(INS)-Z(INBS))

!.....SECOND X THIRD DEFINE   CF AREA
!.....(AREA= SECOND X THIRD)  MULTIPLIED  WITH  *FACS*
    ARX=(DYS*DZT-DYT*DZS)*FACS
    ARY=(DZS*DXT-DZT*DXS)*FACS
    ARZ=(DXS*DYT-DXT*DYS)*FACS
!.....
    AKX=XC(INE)-XC(INP)
    AKY=YC(INE)-YC(INP)
    AKZ=ZC(INE)-ZC(INP)
    ARKSI2=AKX**2+AKY**2+AKZ**2

!.....AREA
    ARE2=ARX**2+ARY**2+ARZ**2
    ARE=SQRT(ARE2)
!
!.....COMPONENTS OF THE UNIT NORMAL VECTOR
    ARER=1.0/(ARE)
!
!.....NORMAL DISTANCE FORM SCALAR PRODUCT
    DELN=(ARX*AKX+ARY*AKY+ARZ*AKZ)*ARER
    DELNR=1./(DELN)

  END SUBROUTINE INIT_GEOM_BF

!> Extract geometric info at an interface between two blocks
!>
!> Fill the common variables of the module
!>
!> @param ISIDE Boundary direction
!> @param I Index along the boundaries
!> @todo Put in INIT_GEOM_T
!> @author Philippe Traore
  SUBROUTINE INIT_GEOM_BLOCK(ISIDE,I)
    IMPLICIT NONE
    INTEGER(IK),INTENT(IN) :: ISIDE,I

    REAL(RK)               :: FACS
    REAL(RK)               :: DXS,DYS,DZS
    REAL(RK)               :: DXT,DYT,DZT

    FACS=REAL(ID_SIGN(ISIDE),rk)
!
!.....CALCULATE INDEX (POINT ON CF-BOUNDARY)
    INE=NBL_ST(NUM_CFI(I,5))+NUM_CFI(I,4)
!
!.....SECOND
    DXS=.5*(X_CFI(I,1)-X_CFI(I,2)+X_CFI(I,4)-X_CFI(I,3))
    DYS=.5*(Y_CFI(I,1)-Y_CFI(I,2)+Y_CFI(I,4)-Y_CFI(I,3))
    DZS=.5*(Z_CFI(I,1)-Z_CFI(I,2)+Z_CFI(I,4)-Z_CFI(I,3))
!.....THIRD
    DXT=.5*(X_CFI(I,1)-X_CFI(I,4)+X_CFI(I,2)-X_CFI(I,3))
    DYT=.5*(Y_CFI(I,1)-Y_CFI(I,4)+Y_CFI(I,2)-Y_CFI(I,3))
    DZT=.5*(Z_CFI(I,1)-Z_CFI(I,4)+Z_CFI(I,2)-Z_CFI(I,3))
!
!.....SECOND X THIRD DEFINE   CF AREA
!.....(AREA= SECOND X THIRD)  MULTIPLIED  WITH  *FACS*
    ARX=(DYS*DZT-DYT*DZS)*FACS
    ARY=(DZS*DXT-DZT*DXS)*FACS
    ARZ=(DXS*DYT-DXT*DYS)*FACS

!.....AREA
    ARE2=ARX**2+ARY**2+ARZ**2
  END SUBROUTINE INIT_GEOM_BLOCK


!> Initialize geometric info inside the domaine
!>
!> The necessity of the dirloop is unclear
!>
!> @param X,Y,Z position of cell's corners
!> @param XC,YC,ZC position of cell's centers
!> @return ARX_T,ARY_T,ARZ_T Projected area in the X Y and Z direction
!> @return AKX_T,AKY_T,AKZ_T Components of normal vector
!> @return ARE_T,ARE2_T,ARER_T Area, Area² and 1/Area
!> @return ARKSI2_T,DELN_T,DELNR_T ??
!>
!> @todo understand the role of the dirloop
!> @todo integrate boundaries and interfaces
!> @author Alexandre Poux
  SUBROUTINE INIT_GEOM_T(X,Y,Z,XC,YC,ZC, &
       ARX_T,ARY_T,ARZ_T,AKX_T,AKY_T,AKZ_T,ARE2_T,ARKSI2_T, &
       ARE_T,ARER_T,DELN_T,DELNR_T)
    IMPLICIT NONE

    REAL(RK)   ,INTENT(IN)  ::      X(NXYZA)
    REAL(RK)   ,INTENT(IN)  ::      Y(NXYZA)
    REAL(RK)   ,INTENT(IN)  ::      Z(NXYZA)
    REAL(RK)   ,INTENT(IN)  ::      XC(NXYZA)
    REAL(RK)   ,INTENT(IN)  ::      YC(NXYZA)
    REAL(RK)   ,INTENT(IN)  ::      ZC(NXYZA)
    REAL(RK)   ,INTENT(OUT) ::   ARX_T(NXYZA,3)
    REAL(RK)   ,INTENT(OUT) ::   ARY_T(NXYZA,3)
    REAL(RK)   ,INTENT(OUT) ::   ARZ_T(NXYZA,3)
    REAL(RK)   ,INTENT(OUT) ::   AKX_T(NXYZA,3)
    REAL(RK)   ,INTENT(OUT) ::   AKY_T(NXYZA,3)
    REAL(RK)   ,INTENT(OUT) ::   AKZ_T(NXYZA,3)
    REAL(RK)   ,INTENT(OUT) ::   ARE2_T(NXYZA,3)
    REAL(RK)   ,INTENT(OUT) :: ARKSI2_T(NXYZA,3)
    REAL(RK)   ,INTENT(OUT) ::    ARE_T(NXYZA,3)
    REAL(RK)   ,INTENT(OUT) ::   ARER_T(NXYZA,3)
    REAL(RK)   ,INTENT(OUT) ::   DELN_T(NXYZA,3)
    REAL(RK)   ,INTENT(OUT) ::  DELNR_T(NXYZA,3)

    REAL(RK),ALLOCATABLE    :: DXS(:,:,:),DYS(:,:,:),DZS(:,:,:),TMPX(:,:,:)
    REAL(RK),ALLOCATABLE    :: DXT(:,:,:),DYT(:,:,:),DZT(:,:,:),TMPY(:,:,:),TMPZ(:,:,:)
    INTEGER(IK)             :: DIR(3),D,NBL,NBLG
    INTEGER(IK)             :: DIR_P(3),DIR_E(3),DIR_S(3),DIR_B(3),DIR_BS(3)


       DO NBL=1,NBLOCK
       
         call SETIND(NBL)

         ALLOCATE(DXS(NJ,NI,NK),DYS(NJ,NI,NK),DZS(NJ,NI,NK),TMPX(NJ,NI,NK))
         ALLOCATE(DXT(NJ,NI,NK),DYT(NJ,NI,NK),DZT(NJ,NI,NK),TMPY(NJ,NI,NK),TMPZ(NJ,NI,NK))
         DXS=0.   ; DYS=0.   ;  DZS=0.
         DXT=0.   ; DYT=0.   ;  DZT=0.

          TMPX=RESHAPE(X(ICST:ICEN),[NJ,NI,NK])
          TMPY=RESHAPE(Y(ICST:ICEN),[NJ,NI,NK])
          TMPZ=RESHAPE(Z(ICST:ICEN),[NJ,NI,NK])

          DIR=(/0,0,1/)
          DIRLOOP1: DO D=1,3
        
              DIR=CSHIFT(DIR,-1)
              DIR_P =(/0,0,0/)
              DIR_E =DIR
              DIR_S =-(/DIR(2),DIR(3),DIR(1)/)
              DIR_B =-(/DIR(3),DIR(1),DIR(2)/)
              DIR_BS=DIR_S+DIR_B

              DXS=.5*(EXTR(TMPX,DIR_P)-EXTR(TMPX,DIR_S)+EXTR(TMPX,DIR_B)-EXTR(TMPX,DIR_BS) )
              DYS=.5*(EXTR(TMPY,DIR_P)-EXTR(TMPY,DIR_S)+EXTR(TMPY,DIR_B)-EXTR(TMPY,DIR_BS) )
              DZS=.5*(EXTR(TMPZ,DIR_P)-EXTR(TMPZ,DIR_S)+EXTR(TMPZ,DIR_B)-EXTR(TMPZ,DIR_BS) )

              DXT=.5*(EXTR(TMPX,DIR_P)+EXTR(TMPX,DIR_S)-EXTR(TMPX,DIR_B)-EXTR(TMPX,DIR_BS) )
              DYT=.5*(EXTR(TMPY,DIR_P)+EXTR(TMPY,DIR_S)-EXTR(TMPY,DIR_B)-EXTR(TMPY,DIR_BS) )
              DZT=.5*(EXTR(TMPZ,DIR_P)+EXTR(TMPZ,DIR_S)-EXTR(TMPZ,DIR_B)-EXTR(TMPZ,DIR_BS) )

!.....SECOND X THIRD DEFINE   CF AREA
!.....(AREA= SECOND X THIRD)
              ARX_T(ICST:ICEN,D)=RESHAPE(DYS*DZT-DYT*DZS,[NIJK])
              ARY_T(ICST:ICEN,D)=RESHAPE(DZS*DXT-DZT*DXS,[NIJK])
              ARZ_T(ICST:ICEN,D)=RESHAPE(DXS*DYT-DXT*DYS,[NIJK])

          ENDDO DIRLOOP1


            TMPX=RESHAPE(XC(ICST:ICEN),[NJ,NI,NK])
            TMPY=RESHAPE(YC(ICST:ICEN),[NJ,NI,NK])
            TMPZ=RESHAPE(ZC(ICST:ICEN),[NJ,NI,NK])

            DIR=(/0,0,1/)
            DIRLOOP2: DO D=1,3
            
                  DIR=CSHIFT(DIR,-1)
                  DIR_P =(/0,0,0/)
                  DIR_E =DIR
                  DIR_S =-(/DIR(2),DIR(3),DIR(1)/)
                  DIR_B =-(/DIR(3),DIR(1),DIR(2)/)
                  DIR_BS=DIR_S+DIR_B

                  AKX_T(ICST:ICEN,D)=RESHAPE(EXTR(TMPX,DIR_E) - EXTR(TMPX,DIR_P),[NIJK])
                  AKY_T(ICST:ICEN,D)=RESHAPE(EXTR(TMPY,DIR_E) - EXTR(TMPY,DIR_P),[NIJK])
                  AKZ_T(ICST:ICEN,D)=RESHAPE(EXTR(TMPZ,DIR_E) - EXTR(TMPZ,DIR_P),[NIJK])
            ENDDO DIRLOOP2
        
           DEALLOCATE(DXS,DYS,DZS,TMPX,TMPY)
           DEALLOCATE(DXT,DYT,DZT,TMPZ)
       ENDDO


!.....AREA
    ARE2_T  =MAX(ARX_T**2+ARY_T**2+ARZ_T**2,0._rk) ! MAX IS FOR ROUND-OFF
    ARE_T   =SQRT(ARE2_T)
    ARER_T  =1.0/(ARE_T+SMALL)
    ARKSI2_T=MAX(AKX_T**2+AKY_T**2+AKZ_T**2,0._rk)
!
!.....COMPONENTS OF THE UNIT NORMAL VECTOR
    DELN_T  =(ARX_T*AKX_T+ARY_T*AKY_T+ARZ_T*AKZ_T)*ARER_T
    DELNR_T =1./(DELN_T+SMALL)


  END SUBROUTINE INIT_GEOM_T



END MODULE GEOM
