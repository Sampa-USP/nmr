      MODULE ROUTINES
      IMPLICIT NONE
      CONTAINS
      

      SUBROUTINE READ_INPUT(TSIM,DT,TWIN,DUMPFILE,OUTFILE,CENTERFILE)
      IMPLICIT NONE
      INTEGER,        INTENT(OUT) :: TSIM
      REAL,           INTENT(OUT) :: DT
      INTEGER,        INTENT(OUT) :: TWIN
      CHARACTER(256), INTENT(OUT) :: DUMPFILE 
      CHARACTER(256), INTENT(OUT) :: OUTFILE
      CHARACTER(256), INTENT(OUT) :: CENTERFILE
      INTEGER                     :: IOSTAT
      CHARACTER(256)              :: LINE,KEYWORD,KEYWORD1
      INTEGER                     :: I
      LOGICAL                     :: FOUNDSHARP

      DUMPFILE=""
      OUTFILE=""
      CENTERFILE=""
     
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
      & MOLECULE,MOL1,NMOL1,HNO,BOX) 
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN)       :: INPUTFILE
      INTEGER,      INTENT(IN)       :: TSIM
      INTEGER,      INTENT(IN)       :: NATOMS
      INTEGER,      INTENT(IN)       :: NTYPES
      INTEGER,      INTENT(IN)       :: MOL1
      INTEGER,      INTENT(IN)       :: NMOL1
      INTEGER,      INTENT(IN)       :: HNO
      REAL,DIMENSION(3,2)            :: BOX  
      REAL,DIMENSION(3,HNO,NMOL1,TSIM) :: MOLECULE  
      INTEGER                        :: I,J,K,M,N
      INTEGER                        :: ATID,ATYPE
      CHARACTER(LEN=50)              :: DUMMY
      REAL,ALLOCATABLE               :: RCM_MOLECULE1(:)
      REAL,ALLOCATABLE               :: ATOMS1(:,:)
      REAL,ALLOCATABLE               :: ATOMS_H(:,:)
      REAL,ALLOCATABLE               :: ATOMS_C(:,:)
      REAL,ALLOCATABLE               :: ATOMS_O(:,:)
      REAL,ALLOCATABLE               :: R1(:)
      INTEGER                        :: TYPE1
      INTEGER                        :: MOL,ID,L
      REAL                           :: LX,LY,LZ
      REAL, ALLOCATABLE              :: LBOX(:)
      TYPE1=3  ! WATER
     
      ALLOCATE(LBOX(3))
      LBOX=0.000000
      ALLOCATE(ATOMS1(MOL1,3))
      ALLOCATE(ATOMS_H(MOL1,3))
      ALLOCATE(ATOMS_C(MOL1,3))
      ALLOCATE(ATOMS_O(MOL1,3))
      ALLOCATE(R1(3))
      ATOMS1=0.00000
      
      PRINT*, " "
      PRINT*, "READING DUMP FILE"
      PRINT*, " "
      PRINT*, " "
      OPEN(8,FILE=INPUTFILE,STATUS='OLD',ACTION='READ')
      DO I=1,TSIM
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
               READ(8,*)ATID, ATYPE,MOL,(ATOMS1(K,M),M=1,3)  
               ATOMS1(K,1)=ATOMS1(K,1)-BOX(1,1)
               ATOMS1(K,1)=MODULO(ATOMS1(K,1),LBOX(1))
               ATOMS1(K,2)=ATOMS1(K,2)-BOX(2,1)
               ATOMS1(K,2)=MODULO(ATOMS1(K,2),LBOX(2))           
               ATOMS1(K,3)=ATOMS1(K,3)-BOX(3,1)
               ATOMS1(K,3)=MODULO(ATOMS1(K,3),LBOX(3))
             END DO   
             
            
         DO N=1,8
            DO K=1,3
               ATOMS_C(N,K)=ATOMS1(N,K)
            END DO
         END DO

         DO N=9,10
            DO K=1,3
               ATOMS_O(N,K)=ATOMS1(N,K)
            END DO
         END DO

         DO N=11,26
           DO K=1,3
              ATOMS_H(N,K)=ATOMS1(N,K)
           END DO
         END DO
         
         DO K=1,3
           ATOMS1(1,K)=ATOMS_C(1,K)
           ATOMS1(2,K)=ATOMS_H(11,K)
           ATOMS1(3,K)=ATOMS_H(12,K)
           ATOMS1(4,K)=ATOMS_H(13,K)
           ATOMS1(5,K)=ATOMS_C(2,K)
           ATOMS1(6,K)=ATOMS_H(14,K)
           ATOMS1(7,K)=ATOMS_H(15,K)
           ATOMS1(8,K)=ATOMS_C(3,K)
           ATOMS1(9,K)=ATOMS_H(16,K)
           ATOMS1(10,K)=ATOMS_H(17,K)
           ATOMS1(11,K)=ATOMS_C(4,K)
           ATOMS1(12,K)=ATOMS_H(18,K)
           ATOMS1(13,K)=ATOMS_H(19,K)
           ATOMS1(14,K)=ATOMS_C(5,K)
           ATOMS1(15,K)=ATOMS_H(20,K)
           ATOMS1(16,K)=ATOMS_H(21,K)
           ATOMS1(17,K)=ATOMS_C(6,K)
           ATOMS1(18,K)=ATOMS_H(22,K)
           ATOMS1(19,K)=ATOMS_H(23,K)
           ATOMS1(20,K)=ATOMS_C(7,K)
           ATOMS1(21,K)=ATOMS_H(24,K)
           ATOMS1(22,K)=ATOMS_H(25,K)
           ATOMS1(23,K)=ATOMS_C(8,K)
           ATOMS1(24,K)=ATOMS_O(9,K)
           ATOMS1(25,K)=ATOMS_O(10,K)
           ATOMS1(26,K)=ATOMS_H(26,K)
         END DO
        
      
         ID = 1

        DO K=1,3
           R1(K) = ATOMS1(ID,K)
        END DO
        DO L=1,4
           ID = ID + 1
           DO K=1,3
           ATOMS1(ID,K) = ATOMS1(ID,K) - R1(K)
           ATOMS1(ID,K) = ATOMS1(ID,K) - NINT(ATOMS1(ID,K)&
                          &/LBOX(K))*LBOX(K)
           ATOMS1(ID,K) = ATOMS1(ID,K) + R1(K)
           END DO
        END DO

        DO L=1,7

           DO K=1,3
               R1(K) = ATOMS1(ID,K)
           END DO

        DO N=1,3
           ID = ID + 1
           DO K=1,3
              ATOMS1(ID,K) = ATOMS1(ID,K) - R1(K)
              ATOMS1(ID,K) = ATOMS1(ID,K) - NINT(ATOMS1(ID,K)/&
                            LBOX(K))*LBOX(K)
              ATOMS1(ID,K) = ATOMS1(ID,K) + R1(K)
          END DO
        END DO

      END DO
      



       
         DO K = 1,3
            MOLECULE(K,1,J,I) = ATOMS1(2,K)
            MOLECULE(K,2,J,I) = ATOMS1(3,K)
            MOLECULE(K,3,J,I) = ATOMS1(4,K)
            MOLECULE(K,4,J,I) = ATOMS1(6,K)
            MOLECULE(K,5,J,I) = ATOMS1(7,K)
            MOLECULE(K,6,J,I) = ATOMS1(9,K)
            MOLECULE(K,7,J,I) = ATOMS1(10,K)
            MOLECULE(K,8,J,I) = ATOMS1(12,K)
            MOLECULE(K,9,J,I) = ATOMS1(13,K)
            MOLECULE(K,10,J,I) = ATOMS1(15,K)
            MOLECULE(K,11,J,I) = ATOMS1(16,K)
            MOLECULE(K,12,J,I) = ATOMS1(18,K)
            MOLECULE(K,13,J,I) = ATOMS1(19,K)
            MOLECULE(K,14,J,I) = ATOMS1(21,K)
            MOLECULE(K,15,J,I) = ATOMS1(22,K)
            MOLECULE(K,16,J,I) = ATOMS1(26,K)
        END DO

      END DO
      ENDDO
      PRINT*, " "
      PRINT*, " " 
      DEALLOCATE(LBOX,ATOMS1,ATOMS_C,ATOMS_H,ATOMS_O,R1)
      CLOSE(8) 

      RETURN

      END SUBROUTINE READ_DUMP
      

      SUBROUTINE DISTANCE(DINTRA,DINTER,THETA1,&
      &          THETA2,POSITIONS_T0,NATOMS,MOL,BOX)
      IMPLICIT NONE
      INTEGER,                           INTENT (IN)  :: NATOMS
      INTEGER,                           INTENT (IN)  :: MOL
      REAL,  DIMENSION(3,MOL,NATOMS),      INTENT (IN)  :: POSITIONS_T0
      REAL,  DIMENSION(3,2),             INTENT (IN)  :: BOX           
      REAL,  DIMENSION(MOL,MOL,NATOMS),          INTENT (OUT) :: DINTRA  ! INTRAMOLECULAR DISTANCE
      REAL,  DIMENSION(MOL,MOL,NATOMS,NATOMS), INTENT (OUT) :: DINTER  ! INTERMOLECULAR DISTANCE
      REAL,  DIMENSION(MOL,MOL,NATOMS),          INTENT (OUT) :: THETA1
      REAL,  DIMENSION(MOL,MOL,NATOMS,NATOMS), INTENT (OUT) :: THETA2
      REAL                                            :: DX,DY,DZ
      REAL,  ALLOCATABLE                              :: LBOX(:)
      REAL,  ALLOCATABLE                              :: INTRA(:,:,:,:)
     REAL,  ALLOCATABLE                              :: INTER(:,:,:,:,:)
      INTEGER                                         :: ISLAB
      INTEGER                                         :: I,J,K,L,M


      ALLOCATE(LBOX(3))
      LBOX=0.00
      DO K=1,3
         LBOX(K) = BOX(K,2) - BOX(K,1)
      ENDDO
           

      ALLOCATE(INTRA(3,MOL,MOL,NATOMS))

      !$OMP PARALLEL DO PRIVATE(I,J)
      DO I = 1,NATOMS
         DO K = 1, MOL
         DO L = K+1, MOL
            DO J = 1,3
               INTRA(J,L,K,I) = POSITIONS_T0(J,K,I)-POSITIONS_T0(J,L,I)
               INTRA(J,L,K,I) = INTRA(J,L,K,I)-NINT(INTRA(J,L,K,I)&
                                /LBOX(J))*LBOX(J)    
            ENDDO
            
         DINTRA(K,L,I) = SQRT(INTRA(1,L,K,I)**2+INTRA(2,L,K,I)&
                      &**2+INTRA(3,L,K,I)**2)
         THETA1(K,L,I) = INTRA(3,L,K,I)/DINTRA(K,L,I)
        END DO
       END DO 
      END DO 
      !$OMP END PARALLEL DO

      DEALLOCATE(INTRA)



      ALLOCATE(INTER(3,MOL,MOL,NATOMS,NATOMS))
      !$OMP PARALLEL DO PRIVATE(L,I,J)   
      DO I = 1, NATOMS
         DO L =I+1, NATOMS
          DO K = 1, MOL
           DO M = K+1, MOL
            DO J = 1,3
               INTER(J,M,K,L,I) =&
               &POSITIONS_T0(J,K,I) - POSITIONS_T0(J,M,L)
               INTER(J,M,K,L,I) =&
               &INTER(J,M,K,L,I)-&
               &NINT(INTER(J,M,K,L,I)/(LBOX(J)))*(LBOX(J))
            ENDDO
            
            DINTER(K,M,I,L) = SQRT(INTER(1,M,K,L,I)**2+&
            &INTER(2,M,K,L,I)**2+INTER(3,M,K,L,I)**2)
            
            THETA2(K,M,I,L)=INTER(3,M,K,L,I)/DINTER(K,M,I,L)
          END DO
         END DO
       ENDDO
      ENDDO
      !$OMP END PARALLEL DO
      

      DEALLOCATE(LBOX)
      DEALLOCATE(INTER)
      RETURN 

    

      END SUBROUTINE DISTANCE
      

      SUBROUTINE CORRELATION(G_INTRA, G_INTER, DINTRA_T0,&
      & DINTRA_T, DINTER_T0, DINTER_T, THETA1_T0, THETA2_T0, THETA1_T, &
      & THETA2_T,NATOMS,HNO)
      IMPLICIT NONE
      INTEGER,                               INTENT (IN)  :: NATOMS
      INTEGER,                               INTENT (IN)  :: HNO
      REAL, DIMENSION(HNO,HNO,NATOMS),      INTENT (IN)  :: DINTRA_T0 
      REAL, DIMENSION(HNO,HNO,NATOMS,NATOMS),INTENT (IN)  :: DINTER_T 
      REAL, DIMENSION(HNO,HNO,NATOMS),       INTENT (IN)  :: DINTRA_T
      REAL, DIMENSION(HNO,HNO,NATOMS,NATOMS), INTENT (IN)  :: DINTER_T0 
      REAL, DIMENSION(HNO,HNO,NATOMS),   INTENT (IN)  :: THETA1_T0
      REAL, DIMENSION(HNO,HNO,NATOMS),   INTENT (IN)  :: THETA1_T    
      REAL, DIMENSION(HNO,HNO,NATOMS,NATOMS), INTENT (IN)  :: THETA2_T0 
      REAL, DIMENSION(HNO,HNO,NATOMS,NATOMS), INTENT (IN)  :: THETA2_T 
      REAL,                             INTENT (OUT) :: G_INTRA         
      REAL,                             INTENT (OUT) :: G_INTER
      INTEGER                                        :: ISLAB,ID,NN,M
      INTEGER                                        :: IATOMS
      INTEGER                                        :: I,J,K,L, NT_AVE

      G_INTRA = 0.0
      G_INTER = 0.0

 
      !$OMP PARALLEL DO REDUCTION(+:G_INTRA) PRIVATE(J) 
      DO J = 1, NATOMS 
        DO   K = 1, HNO
          DO   L = K+1, HNO
             G_INTRA = G_INTRA+&
            &((3.D0*THETA1_T(K,L,J)**2-1)*(3.D0*THETA1_T0(K,L,J)**2-1)/&
            &((DINTRA_T(K,L,J)**3)*(DINTRA_T0(K,L,J)**3)))
          END DO
        END DO
      END DO
      !$OMP END PARALLEL DO
      
!         print*, G_INTRA
      !$OMP PARALLEL DO REDUCTION(+:G_INTER) PRIVATE(L,J)   
      DO J = 1, NATOMS
         DO L =J+1, NATOMS
          DO  K = 1, HNO
            DO  M = K+1, HNO
            G_INTER=G_INTER+((3*THETA2_T(K,M,J,L)**2-1)*&
            &(3*THETA2_T0(K,M,J,L)**2-1)/((DINTER_T(K,M,J,L)**3)*&
            &(DINTER_T0(K,M,J,L)**3)))
            ENDDO
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

