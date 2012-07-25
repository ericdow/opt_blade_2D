C
C     ******************************************************************
C
      SUBROUTINE READ_CGNS_2D (FNAME_IN, PREC)
C
C     ******************************************************************
C     *                                                                *
C     *   READ IN CGNS FILE                                            *
C     *                                                                *
C     ******************************************************************
C
      USE MESH_VAR
C
      IMPLICIT NONE
C
      INCLUDE 'cgnslib_f.h'
C
C     ******************************************************************
C
C     LOCAL VARIABLES
C
C     ******************************************************************
C
      INTEGER   :: IER, I_FILE, I_BASE, I_ZONE, I_COORD, I_FLOW, I_FIELD
      INTEGER   :: PREC
      INTEGER   :: I, J, K, N, IMAX, JMAX, KMAX, LOC
      INTEGER   :: MAXS, C(3), D(3)
      REAL*8    :: DX, DY, DZ, DTH, S1, S2
      REAL*8    :: DXI, DYI, DZI, DXO, DYO, DZO
      REAL*8    :: XTMP, YTMP, ZTMP, RTMP, TH1, TH2, PI
      REAL(8), POINTER :: SOL(:,:,:)
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: XC, YC, ZC
      REAL(8), DIMENSION(:,:), ALLOCATABLE   :: DIJK
      CHARACTER :: ZONENAME*32, SOLNNAME*32, BASENAME*32, COORDNAME*32
      CHARACTER :: FIELDNAME*32, FNAME_IN*32, STATE*32, ANAME*32

      PI = 4.0*ATAN(1.0)

C
C     OPEN CGNS FILE TO READ
C
      CALL CG_OPEN_F(FNAME_IN, CG_MODE_READ, I_FILE, IER)

      IF (IER.NE.CG_OK) THEN
           CALL CG_ERROR_EXIT_F
           STOP
      END IF
      
      CALL CG_NBASES_F(I_FILE, N_BASE, IER)

      DO I_BASE = 1,N_BASE      
          CALL CG_BASE_READ_F(I_FILE, I_BASE, BASENAME,
     .                        ICELLDIM, IPHYSDIM, IER)
          M_BASE = I_BASE
C         SOLUTION FORMAT          
          IF (BASENAME.EQ.'project3d') EXIT
C         GRIDFILE FORMAT          
          IF (BASENAME.EQ.'Base#1') EXIT          
      END DO
C
C     DEFINE THE MESH DIMENSION
C
      CALL CG_NZONES_F(I_FILE, M_BASE, N_ZONE, IER)
      ALLOCATE(DIJK(N_ZONE,3))

      MAXS = -1

      DO I_ZONE = 1,N_ZONE

          CALL CG_ZONE_TYPE_F(I_FILE, M_BASE, I_ZONE, ZONETYPE, IER)

          IF (ZONETYPE .NE. Structured) THEN
              PRINT *, 'ZONETYPE of ZONE ', I_ZONE, ' must be ',
     .                 Structured, ', Now ', ZONETYPE
              STOP
          END IF

          CALL CG_ZONE_READ_F(I_FILE, M_BASE, I_ZONE,
     .                        ZONENAME, ISIZE, IER)
C      
C         CHECK THE NUMBER OF FLOW SOLUTIONS AND NUMBER OF FIELDS
C
          CALL CG_NSOLS_F(I_FILE, M_BASE, I_ZONE, N_FLOW, IER)

          CALL CG_NFIELDS_F(I_FILE, M_BASE, I_ZONE, 1, N_FIELD, IER)

          CALL CG_SOL_INFO_F(I_FILE, M_BASE, I_ZONE, 1, SOLNNAME, LOC, 
     .                       IER)

C          IF ((LOC.NE.Vertex).AND.(LOC.NE.CG_Null)) THEN
C              PRINT *, 'GridLocation must be Vertex'
C              STOP
C          END IF

C
C         CHECK THE COORDINATES
C
          CALL CG_NCOORDS_F(I_FILE, M_BASE, I_ZONE, N_COORD, IER)
C
C         ALLOCATE SPACE
C
          IF (ALLOCATED(X)) DEALLOCATE(X,Y,Z)
          ALLOCATE(X(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))
          ALLOCATE(Y(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))
          ALLOCATE(Z(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))
          IRMIN(1:3) = 1
          IRMAX(1:3) = ISIZE(1:3,1)
C
C         READ THE COORDINATES
C
          IF (PREC == 0) THEN

             IF (ALLOCATED(XS)) DEALLOCATE(XS,YS,ZS)
             ALLOCATE(XS(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))
             ALLOCATE(YS(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))
             ALLOCATE(ZS(ISIZE(1,1), ISIZE(2,1), ISIZE(3,1)))

             CALL CG_COORD_INFO_F(I_FILE, M_BASE, I_ZONE, 1,
     .                            DATATYPE, COORDNAME, IER)

             CALL CG_COORD_READ_F(I_FILE, M_BASE, I_ZONE, COORDNAME,
     .            DATATYPE, IRMIN, IRMAX, XS, IER)

             CALL CG_COORD_INFO_F(I_FILE, M_BASE, I_ZONE, 2,
     .                            DATATYPE, COORDNAME, IER)

             CALL CG_COORD_READ_F(I_FILE, M_BASE, I_ZONE, COORDNAME,
     .            DATATYPE, IRMIN, IRMAX, YS, IER)

             CALL CG_COORD_INFO_F(I_FILE, M_BASE, I_ZONE, 3,
     .                            DATATYPE, COORDNAME, IER)

             CALL CG_COORD_READ_F(I_FILE, M_BASE, I_ZONE, COORDNAME,
     .            DATATYPE, IRMIN, IRMAX, ZS, IER)

             DO I = 1,ISIZE(1,1)
             DO J = 1,ISIZE(2,1)
             DO K = 1,ISIZE(3,1)
                X(I,J,K) = DBLE(XS(I,J,K))
                Y(I,J,K) = DBLE(YS(I,J,K))
                Z(I,J,K) = DBLE(ZS(I,J,K))
             END DO
             END DO
             END DO

          ELSE

             CALL CG_COORD_INFO_F(I_FILE, M_BASE, I_ZONE, 1,
     .                            DATATYPE, COORDNAME, IER)

             CALL CG_COORD_READ_F(I_FILE, M_BASE, I_ZONE, COORDNAME,
     .            DATATYPE, IRMIN, IRMAX, X, IER)

             CALL CG_COORD_INFO_F(I_FILE, M_BASE, I_ZONE, 2,
     .                            DATATYPE, COORDNAME, IER)

             CALL CG_COORD_READ_F(I_FILE, M_BASE, I_ZONE, COORDNAME,
     .            DATATYPE, IRMIN, IRMAX, Y, IER)

             CALL CG_COORD_INFO_F(I_FILE, M_BASE, I_ZONE, 3,
     .                            DATATYPE, COORDNAME, IER)

             CALL CG_COORD_READ_F(I_FILE, M_BASE, I_ZONE, COORDNAME,
     .            DATATYPE, IRMIN, IRMAX, Z, IER)

          END IF

C
C         ZONE DIMENSIONS
C
          IMAX = ISIZE(1,1)
          JMAX = ISIZE(2,1)
          KMAX = ISIZE(3,1)

          DIJK(I_ZONE,1) = (X(IMAX,1,1)-X(1,1,1))**2 + 
     .         (Y(IMAX,1,1)-Y(1,1,1))**2 + (Z(IMAX,1,1)-Z(1,1,1))**2
          DIJK(I_ZONE,2) = (X(1,JMAX,1)-X(1,1,1))**2 + 
     .         (Y(1,JMAX,1)-Y(1,1,1))**2 + (Z(1,JMAX,1)-Z(1,1,1))**2
          DIJK(I_ZONE,3) = (X(1,1,KMAX)-X(1,1,1))**2 +
     .         (Y(1,1,KMAX)-Y(1,1,1))**2 + (Z(1,1,KMAX)-Z(1,1,1))**2

          IF (MIN(DIJK(I_ZONE,1),DIJK(I_ZONE,2),DIJK(I_ZONE,3))  
     .                                               < 1.0E-06) THEN
             IF (IMAX+JMAX+KMAX > MAXS) THEN
                O_ZONE = I_ZONE
                MAXS = IMAX+JMAX+KMAX
             END IF
          END IF
          
          PRINT *, ZONENAME
          PRINT *, 'IMAX      JMAX      KMAX'
          PRINT 5, IMAX, JMAX, KMAX     

      END DO
C
C     FIND AND STORE THE H-MESH
C
      H_ZONE = 1

      CALL CG_ZONE_READ_F(I_FILE, M_BASE, H_ZONE,
     .                    ZONENAME, HSIZE, IER)
C
C     CHECK THE COORDINATES
C
      CALL CG_NCOORDS_F(I_FILE, M_BASE, H_ZONE, N_COORD, IER)
C
C     ALLOCATE SPACE
C
      IF (ALLOCATED(X)) DEALLOCATE(X,Y,Z)
      ALLOCATE(X(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
      ALLOCATE(Y(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
      ALLOCATE(Z(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
      IRMIN(1:3) = 1
      IRMAX(1:3) = HSIZE(1:3,1)
C
C     READ THE COORDINATES
C
      IF (PREC == 0) THEN

         CALL CG_COORD_INFO_F(I_FILE, M_BASE, H_ZONE, 1,
     .                        DATATYPE, COORDNAME, IER)

         CALL CG_COORD_READ_F(I_FILE, M_BASE, H_ZONE, COORDNAME,
     .        DATATYPE, IRMIN, IRMAX, XS, IER)

         CALL CG_COORD_INFO_F(I_FILE, M_BASE, H_ZONE, 2,
     .                        DATATYPE, COORDNAME, IER)

         CALL CG_COORD_READ_F(I_FILE, M_BASE, H_ZONE, COORDNAME,
     .        DATATYPE, IRMIN, IRMAX, YS, IER)

         CALL CG_COORD_INFO_F(I_FILE, M_BASE, H_ZONE, 3,
     .                        DATATYPE, COORDNAME, IER)

         CALL CG_COORD_READ_F(I_FILE, M_BASE, H_ZONE, COORDNAME,
     .        DATATYPE, IRMIN, IRMAX, ZS, IER)
         
         DO I = 1,HSIZE(1,1)
         DO J = 1,HSIZE(2,1)
         DO K = 1,HSIZE(3,1)
            X(I,J,K) = DBLE(XS(I,J,K))
            Y(I,J,K) = DBLE(YS(I,J,K))
            Z(I,J,K) = DBLE(ZS(I,J,K))
         END DO
         END DO
         END DO

      ELSE

         CALL CG_COORD_INFO_F(I_FILE, M_BASE, H_ZONE, 1,
     .                        DATATYPE, COORDNAME, IER)

         CALL CG_COORD_READ_F(I_FILE, M_BASE, H_ZONE, COORDNAME,
     .        DATATYPE, IRMIN, IRMAX, X, IER)

         CALL CG_COORD_INFO_F(I_FILE, M_BASE, H_ZONE, 2,
     .                        DATATYPE, COORDNAME, IER)

         CALL CG_COORD_READ_F(I_FILE, M_BASE, H_ZONE, COORDNAME,
     .        DATATYPE, IRMIN, IRMAX, Y, IER)

         CALL CG_COORD_INFO_F(I_FILE, M_BASE, H_ZONE, 3,
     .                        DATATYPE, COORDNAME, IER)

         CALL CG_COORD_READ_F(I_FILE, M_BASE, H_ZONE, COORDNAME,
     .        DATATYPE, IRMIN, IRMAX, Z, IER)

      END IF
C
C     DETERMINE THE ORIENTATION OF THE MESH
C
      CDIR = 1
      FDIR = 2
      SDIR = 3
C
C     DETERMINE THE TRANSLATION FROM SUCTION TO PRESSURE SIDE
C     
      LE_IND = 0
      TE_IND = 0 
      C(FDIR) = 1
      C(CDIR) = 1
      C(SDIR) = 1
      D(FDIR) = HSIZE(FDIR,1)
      D(CDIR) = 1
      D(SDIR) = 1
      DXI = X(D(1),D(2),D(3)) - X(C(1),C(2),C(3))      
      DYI = Y(D(1),D(2),D(3)) - Y(C(1),C(2),C(3))      
      DZI = Z(D(1),D(2),D(3)) - Z(C(1),C(2),C(3))     
      DO I=1,HSIZE(CDIR,1)
         C(FDIR) = 1
         C(CDIR) = I
         C(SDIR) = 1
         D(FDIR) = HSIZE(FDIR,1)
         D(CDIR) = I
         D(SDIR) = 1
         DXO = X(D(1),D(2),D(3)) - X(C(1),C(2),C(3))      
         DYO = Y(D(1),D(2),D(3)) - Y(C(1),C(2),C(3))      
         DZO = Z(D(1),D(2),D(3)) - Z(C(1),C(2),C(3))     
         IF ((DXO-DXI)**2 + (DYO-DYI)**2 + (DZO-DZI)**2.GE.1.0e-8) THEN
             IF (LE_IND.EQ.0) THEN
                 LE_IND = I-1
             END IF
             TE_IND = I+1
         END IF 
      END DO
C
C     WRITE OUT BLADE SURFACE
C
      OPEN(UNIT=304,FILE='blade_surf.dat')
      WRITE(304,'(A,I10)') 'CDIM:', 2*((TE_IND-LE_IND)+1)-1 
      WRITE(304,'(A,I10)') 'SDIM:', HSIZE(SDIR,1) 

C     PRESSURE SIDE
      DO I=LE_IND,TE_IND
      DO J=1,HSIZE(SDIR,1)
         C(FDIR) = 1
         C(CDIR) = I
         C(SDIR) = J
         WRITE(304,'(3E20.8)') X(C(1),C(2),C(3)),Y(C(1),C(2),C(3)),
     .                         Z(C(1),C(2),C(3))    
      END DO
      END DO

C     SUCTION SIDE      
      DO I=TE_IND-1,LE_IND,-1
      DO J=1,HSIZE(SDIR,1)
         C(FDIR) = HSIZE(FDIR,1)
         C(CDIR) = I
         C(SDIR) = J
         XTMP = X(C(1),C(2),C(3))-DXI
         YTMP = Y(C(1),C(2),C(3))-DYI
         ZTMP = Z(C(1),C(2),C(3))-DZI
         WRITE(304,'(3E20.8)') XTMP,YTMP,ZTMP
      END DO
      END DO

      CLOSE(304)
C
C     ******************************************************************
C
      CALL CG_CLOSE_F(I_FILE, IER)
C
C     ******************************************************************
C
C     CLEAN UP
C
      DEALLOCATE(DIJK)

5     FORMAT (1X,I3,6X,I3,6X,I3)
10    FORMAT (I9,I9,I9,I9)
15    FORMAT (E26.19,1X,E26.19)

      RETURN
C
      END
