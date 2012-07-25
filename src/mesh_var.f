C
C     ******************************************************************
C
      MODULE MESH_VAR
C
C     ******************************************************************
C     *                                                                *
C     *   MESH VARIABLES                                               *
C     *                                                                *
C     ******************************************************************
C
      INTEGER   :: N_BASE, N_COORD, N_ZONE, N_FLOW, N_FIELD, M_BASE
      INTEGER   :: ZONETYPE, DATATYPE, LOCATION
      INTEGER   :: ICELLDIM, IPHYSDIM, DTYPE
      INTEGER   :: ISIZE(3,3), IRMIN(3), IRMAX(3)
      INTEGER   :: O_ZONE,OSIZE(3,3),T_ZONE,TSIZE(3,3),H_ZONE,HSIZE(3,3)
      INTEGER   :: CDIR, FDIR, SDIR
      INTEGER   :: FLIP
      INTEGER   :: LE_IND,TE_IND

      INTEGER, DIMENSION(:,:), ALLOCATABLE   :: HOL, HOW
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: X    , Y    , Z    ,
     .                                          XT   , YT   , ZT   ,
     .                                          XH   , YH   , ZH   ,
     .                                          XMOD , YMOD , ZMOD ,
     .                                          XTMOD, YTMOD, ZTMOD,
     .                                          XHMOD, YHMOD, ZHMOD
      REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: XS, YS, ZS

      END MODULE MESH_VAR
