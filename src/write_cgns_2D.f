C
C     ******************************************************************
C
      SUBROUTINE WRITE_CGNS_2D (FNAME_OUT, PREC)
C
C     ******************************************************************
C     *                                                                *
C     *   WRITE CGNS FILE                                              *
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
      INTEGER   :: I, J, K, C(3), D(3), PREC
      REAL(8), POINTER :: XC(:,:,:), YC(:,:,:), ZC(:,:,:)

      CHARACTER :: FNAME_OUT*32

      CALL CG_OPEN_F(FNAME_OUT, CG_MODE_MODIFY, I_FILE, IER)
C
C     WRITE OUT THE MODIFIED COORDINATES
C
      IF (PREC == 0) THEN

         IF (ALLOCATED(XS)) DEALLOCATE(XS,YS,ZS)
         ALLOCATE(XS(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
         ALLOCATE(YS(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
         ALLOCATE(ZS(HSIZE(1,1), HSIZE(2,1), HSIZE(3,1)))
         DO I = 1,HSIZE(1,1)
         DO J = 1,HSIZE(2,1)
         DO K = 1,HSIZE(3,1)
            XS(I,J,K) = REAL(XMOD(I,J,K))
            YS(I,J,K) = REAL(YMOD(I,J,K))
            ZS(I,J,K) = REAL(ZMOD(I,J,K))
         END DO
         END DO
         END DO

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, H_ZONE, DATATYPE,
     .        'CoordinateX', XS, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, H_ZONE, DATATYPE,
     .        'CoordinateY', YS, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, H_ZONE, DATATYPE,
     .        'CoordinateZ', ZS, I_COORD, IER)


      ELSE

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, H_ZONE, DATATYPE,
     .        'CoordinateX', XMOD, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, H_ZONE, DATATYPE,
     .        'CoordinateY', YMOD, I_COORD, IER)

         CALL CG_COORD_WRITE_F(I_FILE, M_BASE, H_ZONE, DATATYPE,
     .        'CoordinateZ', ZMOD, I_COORD, IER)

      END IF

      RETURN

      END
