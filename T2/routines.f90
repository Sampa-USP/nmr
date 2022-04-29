      MODULE ROUTINES
      IMPLICIT NONE
      CONTAINS


      SUBROUTINE READ_INPUT(TSIM,DT,TWIN,DUMPFILE,OUTFILE)
      IMPLICIT NONE
      INTEGER,        INTENT(OUT) :: TSIM
      REAL,           INTENT(OUT) :: DT
      INTEGER,        INTENT(OUT) :: TWIN
      CHARACTER(256), INTENT(OUT) :: DUMPFILE
      CHARACTER(256), INTENT(OUT) :: OUTFILE
      INTEGER                     :: IOSTAT
      CHARACTER(256)              :: LINE,KEYWORD,KEYWORD1
      INTEGER                     :: I
      LOGICAL                     :: FOUNDSHARP

      DUMPFILE=""
      OUTFILE=""


      PRINT*, "READING INPUT FILE"
      PRINT*, "==================="
      PRINT*, " "
      DO
        READ(*,"(A)",IOSTAT=IOSTAT) LINE
        IF(IOSTAT/=0) EXIT
        FOUNDSHARP=.FALSE.
        DO I=1,LEN(LINE)
           IF(LINE(I:I)=="#") FOUNDSHARP=.TRUE.
           IF(FOUNDSHARP) LINE(I:I)=" "
        ENDDO
        IF(LEN_TRIM(LINE)==0) CYCLE
        READ(LINE,*) KEYWORD

        SELECT CASE(KEYWORD)
          CASE("tsim")
              READ(LINE,*) KEYWORD1,TSIM
              PRINT*, "TSIM:     ", TSIM
          CASE("dt")
              READ(LINE,*) KEYWORD1,DT
              PRINT*, "DT:       ", DT
          CASE("twin")
              READ(LINE,*) KEYWORD1,TWIN
              PRINT*, "TWIN:     ", TWIN
          CASE("dumpfile")
              READ(LINE,*) KEYWORD1,DUMPFILE
              PRINT*, "DUMPFILE: ", TRIM(DUMPFILE)
          CASE("outfile")
              READ(LINE,*) KEYWORD1,OUTFILE
              PRINT*, "OUTFILE: ", TRIM(OUTFILE)
          CASE DEFAULT
              WRITE(0,*) "UNKNOWN KEYWORD :",TRIM(KEYWORD)
              STOP
        END SELECT
      END DO


      IF(DUMPFILE=="") THEN
         WRITE(0,*) "SPECIFY DUMP FILE"
         STOP
      ENDIF

      IF(OUTFILE=="") THEN
         WRITE(0,*) "SPECIFY OUT FILE"
         STOP
      END IF

      RETURN
      END SUBROUTINE READ_INPUT







      SUBROUTINE READ_DUMP(INPUTFILE,TSIM,NATOMS,NTYPES,&
      & MOLECULE1,MOL1,NMOL1,BOX)
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN)       :: INPUTFILE
      INTEGER,      INTENT(IN)       :: TSIM
      INTEGER,      INTENT(IN)       :: NATOMS
      INTEGER,      INTENT(IN)       :: NTYPES
      INTEGER,      INTENT(IN)       :: MOL1
      INTEGER,      INTENT(IN)       :: NMOL1
      REAL,DIMENSION(3,2)            :: BOX
      REAL,DIMENSION(3,3,NMOL1,TSIM) :: MOLECULE1
      INTEGER                        :: I,J,K,M
      INTEGER                        :: ATID,ATYPE
      CHARACTER(LEN=50)              :: DUMMY
      REAL,ALLOCATABLE               :: RCM_MOLECULE1(:)
      REAL,ALLOCATABLE               :: ATOMS1(:,:)
      INTEGER                        :: TYPE1
      INTEGER                        :: MOL
      REAL                           :: LX,LY,LZ
      REAL, ALLOCATABLE              :: LBOX(:)
      TYPE1=3  ! WATER

      ALLOCATE(LBOX(3))
      LBOX=0.000000
      ALLOCATE(ATOMS1(MOL1,3))
      ATOMS1=0.00000

      PRINT*, " "
      PRINT*, "READING DUMP FILE"
      PRINT*, " "
      PRINT*, " "
      OPEN(8,FILE=INPUTFILE,STATUS='OLD',ACTION='READ')
      DO I=1,TSIM
         print*,i,tsim
         READ(8,*)
         READ(8,*)
         READ(8,*)
         READ(8,*)
         READ(8,*)
         READ(8,*)BOX(1,1),BOX(1,2)
         READ(8,*)BOX(2,1),BOX(2,2)
         READ(8,*)BOX(3,1),BOX(3,2)
         READ(8,*)

         DO K=1,3
            LBOX(K) = BOX(K,2) - BOX(K,1)
         END DO


         DO J=1,NMOL1  ! NUMBERS WATER MOLECULES
            DO K=1,MOL1  ! NUMBER OF ATOMS IN A MOLECULE
               READ(8,*)ATID,  ATYPE, MOL, (ATOMS1(K,M),M=1,3)
               ATOMS1(K,1)=ATOMS1(K,1)-BOX(1,1)
               ATOMS1(K,1)=MODULO(ATOMS1(K,1),LBOX(1))
               ATOMS1(K,2)=ATOMS1(K,2)-BOX(2,1)
               ATOMS1(K,2)=MODULO(ATOMS1(K,2),LBOX(2))
               ATOMS1(K,3)=ATOMS1(K,3)-BOX(3,1)
               ATOMS1(K,3)=MODULO(ATOMS1(K,3),LBOX(3))
            END DO

            DO M=2,3
               DO K=1,3
                  ATOMS1(M,K) = ATOMS1(M,K) - ATOMS1(1,K)
                  ATOMS1(M,K) = ATOMS1(M,K) - &
                  & NINT(ATOMS1(M,K)/(LBOX(K)))*(LBOX(K))
                  ATOMS1(M,K) = ATOMS1(M,K) + ATOMS1(1,K)
                  MOLECULE1(K,M,J,I)=ATOMS1(M,K)
               END DO
            ENDDO
         ENDDO


      ENDDO
      PRINT*, " "
      PRINT*, " "
      DEALLOCATE(LBOX,ATOMS1)
      CLOSE(8)

      RETURN

      END SUBROUTINE READ_DUMP


      SUBROUTINE DISTANCE(DINTRA,DINTER,THETA1,&
      &          THETA2,POSITIONS_T0,NATOMS,BOX)
      IMPLICIT NONE
      INTEGER,                          INTENT (IN)   :: NATOMS
      REAL,  DIMENSION(3,3,NATOMS),      INTENT (IN)  :: POSITIONS_T0
      REAL,  DIMENSION(3,2),             INTENT (IN)  :: BOX
      REAL,  DIMENSION(NATOMS),          INTENT (OUT) :: DINTRA  ! INTRAMOLECULAR DISTANCE
      REAL,  DIMENSION(3,NATOMS,NATOMS), INTENT (OUT) :: DINTER  ! INTERMOLECULAR DISTANCE
      REAL,  DIMENSION(NATOMS),          INTENT (OUT) :: THETA1
      REAL,  DIMENSION(3,NATOMS,NATOMS), INTENT (OUT) :: THETA2
      REAL                                            :: DX,DY,DZ,CUTOFF
      REAL,  ALLOCATABLE                              :: LBOX(:)
      REAL,  ALLOCATABLE                              :: INTRA(:,:)
      REAL,  ALLOCATABLE                              :: INTER(:,:,:,:)
      INTEGER                                         :: ISLAB
      INTEGER                                         :: I,J,K,L





      ALLOCATE(LBOX(3))
      LBOX=0.00
      DO K=1,3
         LBOX(K) = BOX(K,2) - BOX(K,1)
      ENDDO

      ALLOCATE(INTRA(3,NATOMS))

      !$OMP PARALLEL DO PRIVATE(I,J)
      DO I = 1,NATOMS
         DO J = 1,3
            INTRA(J,I) = POSITIONS_T0(J,2,I)-POSITIONS_T0(J,3,I)
            INTRA(J,I) = INTRA(J,I)-NINT(INTRA(J,I)/LBOX(J))*LBOX(J)
         ENDDO

         DINTRA(I) = SQRT(INTRA(1,I)**2+INTRA(2,I)**2+INTRA(3,I)**2)
         THETA1(I) = INTRA(3,I)/DINTRA(I)
      END DO
      !$OMP END PARALLEL DO

      DEALLOCATE(INTRA)



      ALLOCATE(INTER(3,3,NATOMS,NATOMS))
!      !$OMP PARALLEL DO PRIVATE(L,I,J)
     !     DO I = 1, NATOMS
      !   DO L =I+1, NATOMS
       !     DO J = 1,3
        !       INTER(J,2,I,L) =&
         !      &POSITIONS_T0(J,2,L) - POSITIONS_T0(J,3,I)
         !      INTER(J,2,I,L) =&
         !1      &INTER(J,2,I,L)-&
        !       &NINT(INTER(J,2,I,L)/(LBOX(J)))*(LBOX(J))
        !    ENDDO

         !   DINTER(2,I,L) = SQRT(INTER(1,2,I,L)**2+&
         !1   &INTER(2,2,I,L)**2+INTER(3,2,I,L)**2)

         !   THETA2(2,I,L)=INTER(3,2,I,L)/DINTER(2,I,L)
       !  ENDDO
      !ENDDO
     ! !$OMP END PARALLEL DO


      !$OMP PARALLEL DO PRIVATE(L,I,J,K)
      DO L=1, NATOMS
         DO I=L+1, NATOMS
             DO K = 2,3
                DO J = 1,3
                   INTER(J,K,I,L) = &
                   &POSITIONS_T0(J,2,L)-POSITIONS_T0(J,K,I)
                   INTER(J,K,I,L) = &
                   &INTER(J,K,I,L)-&
                   &NINT(INTER(J,K,I,L)/(LBOX(J)))*(LBOX(J))
                ENDDO
                DINTER(K,I,L) =&
                &SQRT(INTER(1,K,I,L)**2+INTER(2,K,I,L)**2+&
                &INTER(3,K,I,L)**2)
                THETA2 (K,I,L) = INTER(3,K,I,L)/DINTER(K,I,L)
             ENDDO
         ENDDO
      END DO
      !$OMP END PARALLEL DO

      DEALLOCATE(LBOX)
      DEALLOCATE(INTER)
      RETURN



      END SUBROUTINE DISTANCE













      SUBROUTINE CORRELATION(G_INTRA, G_INTER, DINTRA_T0,&
      & DINTRA_T, DINTER_T0, DINTER_T, THETA1_T0, THETA2_T0, THETA1_T, &
      & THETA2_T,NATOMS)
      IMPLICIT NONE
      INTEGER,                          INTENT (IN)  :: NATOMS
      REAL, DIMENSION(NATOMS),          INTENT (IN)  :: DINTRA_T0
      REAL, DIMENSION(3,NATOMS,NATOMS), INTENT (IN)  :: DINTER_T
      REAL, DIMENSION(NATOMS),          INTENT (IN)  :: DINTRA_T
      REAL, DIMENSION(3,NATOMS,NATOMS), INTENT (IN)  :: DINTER_T0
      REAL, DIMENSION(NATOMS),          INTENT (IN)  :: THETA1_T0
      REAL, DIMENSION(NATOMS),          INTENT (IN)  :: THETA1_T
      REAL, DIMENSION(3,NATOMS,NATOMS), INTENT (IN)  :: THETA2_T0
      REAL, DIMENSION(3,NATOMS,NATOMS), INTENT (IN)  :: THETA2_T
      REAL,                             INTENT (OUT) :: G_INTRA
      REAL,                             INTENT (OUT) :: G_INTER
      INTEGER                                        :: ISLAB,ID,NN
      INTEGER                                        :: IATOMS
      INTEGER                                        :: I,J,K,L, NT_AVE

      G_INTRA = 0.0
      G_INTER = 0.0


      !$OMP PARALLEL DO REDUCTION(+:G_INTRA) PRIVATE(J)
      DO J = 1, NATOMS
         G_INTRA = G_INTRA+&
         &((3.D0*THETA1_T(J)**2-1)*(3.D0*THETA1_T0(J)**2-1)/&
         &((DINTRA_T(J)**3)*(DINTRA_T0(J)**3)))
      END DO
      !$OMP END PARALLEL DO

      !!$OMP PARALLEL DO REDUCTION(+:G_INTER) PRIVATE(L,J)
     ! DO J = 1, NATOMS
      !   DO L =J+1, NATOMS
       !     G_INTER=G_INTER+((3*THETA2_T(2,J,L)**2-1)*&
       !1     &(3*THETA2_T0(2,J,L)**2-1)/((DINTER_T(2,J,L)**3)*&
       !     &(DINTER_T0(2,J,L)**3)))
       !   ENDDO
       !ENDDO
       !!$OMP END PARALLEL DO

       !$OMP PARALLEL DO REDUCTION(+:G_INTER) PRIVATE(L,J,K)
       DO L=1, NATOMS
          DO J=L+1, NATOMS
             DO K =2,3
                G_INTER=G_INTER+((3*THETA2_T(K,J,L)**2-1)*&
                &(3*THETA2_T0(K,J,L)**2-1)/((DINTER_T(K,J,L)**3)*&
                &(DINTER_T0(K,J,L)**3)))
             ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      RETURN

      END SUBROUTINE CORRELATION









      SUBROUTINE WRITE_T2(OUTPUTFILE,AVE_INTRA_HH, AVE_INTER_HH,&
      &TWIN,NWIN,DT, NATOMS)
      CHARACTER(*),          INTENT(IN) :: OUTPUTFILE
      INTEGER,               INTENT(IN) :: NATOMS
      INTEGER,               INTENT(IN) :: TWIN
      INTEGER,               INTENT(IN) :: NWIN
      REAL, DIMENSION(TWIN), INTENT(IN) :: AVE_INTRA_HH
      REAL, DIMENSION(TWIN), INTENT(IN) :: AVE_INTER_HH
      REAL,                  INTENT(IN) :: DT
      INTEGER                           :: ISLAB,I,J
      REAL                              :: DZ
      REAL                              :: TCORR,P2_INTRA, P2_INTER
      REAL                              :: GT_INTRA, GT_INTER


      PRINT*, " "
      PRINT*, "WRITING RESULTS"
      PRINT*, "================="
      PRINT*, " "

  20  FORMAT(F12.1,3X,4F15.5)
  30  FORMAT(A13)
      OPEN(10,FILE=OUTPUTFILE,ACTION='WRITE',STATUS='UNKNOWN')
      WRITE(10,30)"#TCORR P2_INTRA, P2_INTER"
      WRITE(10,*) AVE_INTRA_HH(1), AVE_INTER_HH(1)
      DO J=1,TWIN
         TCORR=DT*(J-1)
         P2_INTRA = AVE_INTRA_HH(J)/(AVE_INTRA_HH(1))
         GT_INTRA = AVE_INTRA_HH(J)
         P2_INTER = AVE_INTER_HH(J)/(AVE_INTER_HH(1))
         GT_INTER = AVE_INTER_HH(J)
         WRITE(10,20)TCORR, P2_INTRA,&
         & P2_INTER, GT_INTRA,GT_INTER
      END DO

      CLOSE(10)

      RETURN

      END SUBROUTINE WRITE_T2


      END MODULE ROUTINES
