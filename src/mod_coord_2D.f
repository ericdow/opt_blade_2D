C
C     ******************************************************************
C
      SUBROUTINE MOD_COORD_2D
C
C     ******************************************************************
C     *                                                                *
C     *   MODIFY COORDINATES                                           *
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
      INTEGER   :: I, J, K, C(3), D(3), JJ, N
      REAL      :: DX, DY, DZ, R, RR, MAG, DENOM, HTIP
      REAL      :: DXI, DYI, DZI, DXO, DYO, DZO
      REAL      :: XTMP, YTMP, ZTMP, RTMP, TH1, TH2, DTH, PI
      REAL      :: DXS, DYS, DZS, DXH, DYH, DZH, DXB, DYB, DZB
      REAL(8), DIMENSION(:), ALLOCATABLE     :: DS, S, RO, RI, XI
      REAL(8), DIMENSION(:,:), ALLOCATABLE   :: DS2, S2, SJK, TJK
      INTEGER, DIMENSION(:), ALLOCATABLE     :: IY
      CHARACTER :: STR*32

      PI = 4.0*ATAN(1.0)

      ALLOCATE(XMOD(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
      ALLOCATE(YMOD(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
      ALLOCATE(ZMOD(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))

      DO I = 1,HSIZE(1,1)
      DO J = 1,HSIZE(2,1)
      DO K = 1,HSIZE(3,1)
         XMOD(I,J,K) = X(I,J,K)
         YMOD(I,J,K) = Y(I,J,K)
         ZMOD(I,J,K) = Z(I,J,K)
      END DO
      END DO
      END DO
C
C     READ IN BLADE SURFACE
C      
      OPEN(UNIT=304,FILE='blade_surf_mod.dat')

      READ(304,'(A,I10)') STR
      READ(304,'(A,I10)') STR

C     PRESSURE SIDE
      DO I=LE_IND,TE_IND
      DO J=1,HSIZE(SDIR,1)
         C(FDIR) = 1
         C(CDIR) = I
         C(SDIR) = J
         READ(304,'(3E20.8)') XMOD(C(1),C(2),C(3)),
     .                        YMOD(C(1),C(2),C(3)),
     .                        ZMOD(C(1),C(2),C(3))   
      END DO
      END DO
      
      C(FDIR) = 1
      C(CDIR) = 1
      C(SDIR) = 1
      D(FDIR) = HSIZE(FDIR,1)
      D(CDIR) = 1
      D(SDIR) = 1
      DXI = X(D(1),D(2),D(3)) - X(C(1),C(2),C(3))      
      DYI = Y(D(1),D(2),D(3)) - Y(C(1),C(2),C(3))      
      DZI = Z(D(1),D(2),D(3)) - Z(C(1),C(2),C(3))    

C     SUCTION SIDE      
      DO I=TE_IND-1,LE_IND,-1
      DO J=1,HSIZE(SDIR,1)
         READ(304,'(3E20.8)') XTMP,YTMP,ZTMP
         C(FDIR) = HSIZE(FDIR,1)
         C(CDIR) = I
         C(SDIR) = J
         XMOD(C(1),C(2),C(3)) = XTMP+DXI
         YMOD(C(1),C(2),C(3)) = YTMP+DYI
         ZMOD(C(1),C(2),C(3)) = ZTMP+DZI
      END DO
      END DO

      CLOSE(304)
C
C     INTERPOLATE TRANSLATION TO PERIODIC REGION
C
      ALLOCATE(S(HSIZE(CDIR,1)))
      ALLOCATE(DS(HSIZE(CDIR,1)))
      
      DO J = 1,HSIZE(SDIR,1)
         C(CDIR) = LE_IND
         C(FDIR) = 1
         C(SDIR) = J
         DX = XMOD(C(1),C(2),C(3))-X(C(1),C(2),C(3))
         DY = YMOD(C(1),C(2),C(3))-Y(C(1),C(2),C(3))
         DZ = ZMOD(C(1),C(2),C(3))-Z(C(1),C(2),C(3))
         S(1) = 0.0
         DS(1) = 0.0
         DO K = 2,LE_IND
            C(SDIR) = J
            C(CDIR) = K
            D(FDIR) = 1
            D(SDIR) = J
            D(CDIR) = K-1
            DS(K) = SQRT((X(C(1),C(2),C(3))-X(D(1),D(2),D(3)))**2 + 
     .                   (Y(C(1),C(2),C(3))-Y(D(1),D(2),D(3)))**2 +
     .                   (Z(C(1),C(2),C(3))-Z(D(1),D(2),D(3)))**2)
            S(K) = S(K-1) + DS(K)
         END DO
         DO K = 2,LE_IND-1
            C(SDIR) = J
            C(CDIR) = K
            XMOD(C(1),C(2),C(3)) = X(C(1),C(2),C(3)) + DX*S(K)/S(LE_IND)
            YMOD(C(1),C(2),C(3)) = Y(C(1),C(2),C(3)) + DY*S(K)/S(LE_IND)
            ZMOD(C(1),C(2),C(3)) = Z(C(1),C(2),C(3)) + DZ*S(K)/S(LE_IND)
         END DO
      END DO

      DO J = 1,HSIZE(SDIR,1)
         C(CDIR) = TE_IND
         C(FDIR) = 1
         C(SDIR) = J
         DX = XMOD(C(1),C(2),C(3))-X(C(1),C(2),C(3))
         DY = YMOD(C(1),C(2),C(3))-Y(C(1),C(2),C(3))
         DZ = ZMOD(C(1),C(2),C(3))-Z(C(1),C(2),C(3))
         S(TE_IND) = 0.0
         DS(TE_IND) = 0.0
         DO K = TE_IND+1,HSIZE(CDIR,1)
            C(SDIR) = J
            C(CDIR) = K
            D(FDIR) = 1
            D(SDIR) = J
            D(CDIR) = K-1
            DS(K) = SQRT((X(C(1),C(2),C(3))-X(D(1),D(2),D(3)))**2 + 
     .                   (Y(C(1),C(2),C(3))-Y(D(1),D(2),D(3)))**2 +
     .                   (Z(C(1),C(2),C(3))-Z(D(1),D(2),D(3)))**2)
            S(K) = S(K-1) + DS(K)
         END DO
         DO K = TE_IND+1,HSIZE(CDIR,1)-1
            C(SDIR) = J
            C(CDIR) = K
            XMOD(C(1),C(2),C(3)) = X(C(1),C(2),C(3)) 
     .                     + DX*(S(HSIZE(CDIR,1))-S(K))/S(HSIZE(CDIR,1))
            YMOD(C(1),C(2),C(3)) = Y(C(1),C(2),C(3)) 
     .                     + DY*(S(HSIZE(CDIR,1))-S(K))/S(HSIZE(CDIR,1))
            ZMOD(C(1),C(2),C(3)) = Z(C(1),C(2),C(3)) 
     .                     + DZ*(S(HSIZE(CDIR,1))-S(K))/S(HSIZE(CDIR,1))
         END DO
      END DO
C      
C     MAKE PERIODIC REGIONS MATCH
C 
      C(FDIR) = 1
      C(CDIR) = 1
      C(SDIR) = 1
      D(FDIR) = HSIZE(FDIR,1)
      D(CDIR) = 1
      D(SDIR) = 1
      DXI = X(D(1),D(2),D(3)) - X(C(1),C(2),C(3))      
      DYI = Y(D(1),D(2),D(3)) - Y(C(1),C(2),C(3))      
      DZI = Z(D(1),D(2),D(3)) - Z(C(1),C(2),C(3))
      
      DO I = 1,HSIZE(SDIR,1)
      DO J = 1,LE_IND
      C(SDIR) = I
      C(CDIR) = J
      C(FDIR) = HSIZE(FDIR,1) 
      D(SDIR) = I
      D(CDIR) = J
      D(FDIR) = 1
         XMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3)) + DXI
         YMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3)) + DYI
         ZMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3)) + DZI
      END DO
      DO J = TE_IND,HSIZE(CDIR,1)
      C(SDIR) = I
      C(CDIR) = J
      C(FDIR) = HSIZE(FDIR,1) 
      D(SDIR) = I
      D(CDIR) = J
      D(FDIR) = 1
         XMOD(C(1),C(2),C(3)) = XMOD(D(1),D(2),D(3)) + DXI
         YMOD(C(1),C(2),C(3)) = YMOD(D(1),D(2),D(3)) + DYI
         ZMOD(C(1),C(2),C(3)) = ZMOD(D(1),D(2),D(3)) + DZI
      END DO
      END DO
C      
C     INTERPOLATE H-MESH COORDINATES ON HUB AND SHROUD SURFACES
C
      DEALLOCATE(S)
      DEALLOCATE(DS)
      ALLOCATE(S(HSIZE(FDIR,1)))
      ALLOCATE(DS(HSIZE(FDIR,1)))

      DO I = 1,HSIZE(CDIR,1)
      DO JJ = 1,2 
         IF (JJ.EQ.1) J = 1 
         IF (JJ.EQ.2) J = HSIZE(SDIR,1)
         C(CDIR) = I
         C(SDIR) = J
         C(FDIR) = 1

         DXI = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
         DYI = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
         DZI = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))

         C(FDIR) = HSIZE(FDIR,1)

         DXO = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
         DYO = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
         DZO = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
         
         S(1) = 0.0
         DS(1) = 0.0
         DO K = 2,HSIZE(FDIR,1)
            C(CDIR) = I
            C(SDIR) = J
            C(FDIR) = K
            D(CDIR) = I
            D(SDIR) = J
            D(FDIR) = K-1
            DS(K) = SQRT((X(C(1),C(2),C(3))-X(D(1),D(2),D(3)))**2 + 
     .                   (Y(C(1),C(2),C(3))-Y(D(1),D(2),D(3)))**2 +
     .                   (Z(C(1),C(2),C(3))-Z(D(1),D(2),D(3)))**2)
            S(K) = S(K-1) + DS(K)
         END DO
         DO K = 2,HSIZE(FDIR,1)-1
            C(CDIR) = I
            C(SDIR) = J
            C(FDIR) = K
            DX = DXI*(S(HSIZE(FDIR,1)) - S(K))/S(HSIZE(FDIR,1)) + 
     .           DXO*S(K)/S(HSIZE(FDIR,1))
            DY = DYI*(S(HSIZE(FDIR,1)) - S(K))/S(HSIZE(FDIR,1)) + 
     .           DYO*S(K)/S(HSIZE(FDIR,1))
            DZ = DZI*(S(HSIZE(FDIR,1)) - S(K))/S(HSIZE(FDIR,1)) + 
     .           DZO*S(K)/S(HSIZE(FDIR,1))
            XMOD(C(1),C(2),C(3)) = X(C(1),C(2),C(3)) + DX
            YMOD(C(1),C(2),C(3)) = Y(C(1),C(2),C(3)) + DY
            ZMOD(C(1),C(2),C(3)) = Z(C(1),C(2),C(3)) + DZ
         END DO
      END DO
      END DO
C
C     BILINEAR INTERPOLATION OF H-MESH INTERIOR
C
      ALLOCATE(SJK(HSIZE(FDIR,1),HSIZE(SDIR,1)))
      ALLOCATE(TJK(HSIZE(FDIR,1),HSIZE(SDIR,1)))
      DO I = 2,HSIZE(CDIR,1)-1
         C(FDIR) = 1
         C(CDIR) = I
         C(SDIR) = HSIZE(SDIR,1)

         SJK(1,1:HSIZE(SDIR,1)) = 0.0
         TJK(1:HSIZE(FDIR,1),1) = 0.0
C        FILL SJK
         DO J = 2,HSIZE(FDIR,1)
         DO K = 1,HSIZE(SDIR,1)
            C(FDIR) = J
            C(SDIR) = K
            D(CDIR) = I
            D(FDIR) = J-1
            D(SDIR) = K
            SJK(J,K) = SJK(J-1,K) + 
     .           SQRT((X(C(1),C(2),C(3))-X(D(1),D(2),D(3)))**2 + 
     .                (Y(C(1),C(2),C(3))-Y(D(1),D(2),D(3)))**2 +
     .                (Z(C(1),C(2),C(3))-Z(D(1),D(2),D(3)))**2)
         END DO
         END DO
C        FILL TJK
         DO J = 1,HSIZE(FDIR,1)
         DO K = 2,HSIZE(SDIR,1)
            C(FDIR) = J
            C(SDIR) = K
            D(CDIR) = I
            D(FDIR) = J
            D(SDIR) = K-1
            TJK(J,K) = TJK(J,K-1) + 
     .           SQRT((X(C(1),C(2),C(3))-X(D(1),D(2),D(3)))**2 + 
     .                (Y(C(1),C(2),C(3))-Y(D(1),D(2),D(3)))**2 +
     .                (Z(C(1),C(2),C(3))-Z(D(1),D(2),D(3)))**2)
         END DO
         END DO
C        SHIFT INTERIOR POINTS USING BILINEAR INTERPOLATION
         DO J = 2,HSIZE(FDIR,1)-1
C           SHIFT VECTOR AT HUB            
            C(FDIR) = J
            C(CDIR) = I
            C(SDIR) = 1
            DXH = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
            DYH = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
            DZH = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
C           SHIFT VECTOR AT SHROUD            
            C(SDIR) = HSIZE(SDIR,1)
            DXS = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
            DYS = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
            DZS = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
            DO K = 2,HSIZE(SDIR,1)-1
C              SHIFT VECTOR AT FIRST BLADE SURFACE
               C(FDIR) = 1
               C(SDIR) = K
               DXI = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
               DYI = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
               DZI = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
C              SHIFT VECTOR AT SECOND BLADE SURFACE
               C(FDIR) = HSIZE(FDIR,1)
               C(SDIR) = K
               DXO = XMOD(C(1),C(2),C(3)) - X(C(1),C(2),C(3)) 
               DYO = YMOD(C(1),C(2),C(3)) - Y(C(1),C(2),C(3)) 
               DZO = ZMOD(C(1),C(2),C(3)) - Z(C(1),C(2),C(3))
C              MOVE THE INTERIOR POINTS
               DENOM = 1.0/SJK(J,K)+1.0/(SJK(HSIZE(FDIR,1),K)-SJK(J,K))
     .               + 1.0/TJK(J,K)+1.0/(TJK(J,HSIZE(SDIR,1))-TJK(J,K))
               DX = (DXH/TJK(J,K) + DXS/(TJK(J,HSIZE(SDIR,1))-TJK(J,K)) 
     .            +  DXI/SJK(J,K) + DXO/(SJK(HSIZE(FDIR,1),K)-SJK(J,K)))
     .            /  DENOM
               DY = (DYH/TJK(J,K) + DYS/(TJK(J,HSIZE(SDIR,1))-TJK(J,K)) 
     .            +  DYI/SJK(J,K) + DYO/(SJK(HSIZE(FDIR,1),K)-SJK(J,K)))
     .            /  DENOM
               DZ = (DZH/TJK(J,K) + DZS/(TJK(J,HSIZE(SDIR,1))-TJK(J,K)) 
     .            +  DZI/SJK(J,K) + DZO/(SJK(HSIZE(FDIR,1),K)-SJK(J,K)))
     .            /  DENOM
               C(FDIR) = J
               C(CDIR) = I
               C(SDIR) = K
               XMOD(C(1),C(2),C(3)) = X(C(1),C(2),C(3)) + DX
               YMOD(C(1),C(2),C(3)) = Y(C(1),C(2),C(3)) + DY
               ZMOD(C(1),C(2),C(3)) = Z(C(1),C(2),C(3)) + DZ
            END DO
         END DO
      END DO
C     
      RETURN
C
      END
