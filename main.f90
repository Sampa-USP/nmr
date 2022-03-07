       PROGRAM T2_OCTANOIC
       USE ROUTINES
       IMPLICIT NONE
       ! DEFINE PARAMETERS
       REAL,ALLOCATABLE        :: BOX(:,:)  
       CHARACTER(256)          :: DUMPFILE
       CHARACTER(256)          :: OUTFILE
       CHARACTER(256)          :: CENTERFILE
       INTEGER                 :: NATOMS                ! NUMBERS OF ATOMS 
       INTEGER                 :: NTYPES                ! NUMBERS OF TYPES OF ATOMS 
       INTEGER                 :: MOL1                  ! NUMBERS OF ATOMS TO CONFORM THE MOLECULE TYPE 1
       INTEGER                 :: NMOL1                 ! NUMBERS OF MOLECULE OF TYPE 1
       INTEGER                 :: TSIM                  ! NUMBERS OF SNAPSHOT CONFIGURATION 
       REAL, ALLOCATABLE       :: MOLECULE1(:,:,:,:)      ! POSITIONS OF CENTER MASS OF ALL MOLECULES THE TYPE 1
       REAL                    :: DT                    ! TIME STEP
       INTEGER                 :: TCORR                 ! TIME CORRELATION
       INTEGER                 :: TWIN                  ! TIME WINDOW
       REAL, ALLOCATABLE       :: AVE_INTRA_HH(:)         ! MEAN SQUARE DISPLACEMENT FOR THE MOLECULE TYPE 2
       REAL, ALLOCATABLE       :: AVE_INTER_HH(:)
       REAL, ALLOCATABLE       :: POSITIONS1_T0(:,:,:)     ! MOLECULES TYPE 1 AT THE TIME T0 
       REAL, ALLOCATABLE       :: POSITIONS1_T(:,:,:)      ! MOLECULES TYPE 1 SURVIVOR AT THE TIME T 
       REAL, ALLOCATABLE       :: DINTRA_T0(:,:,:)
       REAL, ALLOCATABLE       :: DINTER_T0(:,:,:,:)
       REAL, ALLOCATABLE       :: DINTRA_T(:,:,:)
       REAL, ALLOCATABLE       :: DINTER_T(:,:,:,:)
       REAL, ALLOCATABLE       :: THETA1_T0(:,:,:)
       REAL, ALLOCATABLE       :: THETA1_T(:,:,:)
       REAL, ALLOCATABLE       :: THETA2_T0(:,:,:,:)
       REAL, ALLOCATABLE       :: THETA2_T(:,:,:,:)
       REAL, ALLOCATABLE       :: THETA1(:,:,:)
       REAL, ALLOCATABLE       :: THETA2(:,:,:,:)
       INTEGER                 :: NT, ID0, NT_AVE
       REAL, ALLOCATABLE       :: P2_HH(:)
       REAL                    :: G_INTRA
       REAL                    :: G_INTER
       REAL, ALLOCATABLE       :: DIST_HH(:)
       REAL, ALLOCATABLE       :: DINTRA(:,:,:)
       REAL, ALLOCATABLE       :: DINTER(:,:,:,:)
       INTEGER, ALLOCATABLE    :: ID(:)
       INTEGER                 :: TT,T0,NWIN,T, NMOL2
       REAL                    :: XLO,XHI
       REAL                    :: YLO,YHI
       REAL                    :: ZLO,ZHI
       REAL                    :: LX,LY, LZ
       INTEGER                 :: ISLAB,I,J,K,L,HNO
       REAL                    :: START_DUMP,FINISH_DUMP
       REAL, ALLOCATABLE       :: HH_T0(:), HH_T(:)
       
       CALL READ_INPUT(TSIM,DT,TWIN,DUMPFILE,OUTFILE,CENTERFILE)
       

       MOL1 = 26
       NTYPES= 1
       HNO = 16
!       NMOL2=85 


       OPEN(12, FILE=DUMPFILE,ACTION="READ", STATUS='OLD')
       DO I=1,3
          READ(12,*)
       ENDDO
       READ(12,*) NATOMS
       CLOSE(12)

       NMOL2=NATOMS/MOL1
       !IF(NMOL2.gt.50)NMOL2=50
       
       NMOL1=NATOMS/26
       PRINT*, "NATOMS: ", NATOMS
       PRINT*, "NMOL1:  ", NMOL1
       
       ! ALLOCATION OF DYNAMICAL ARRAYS
       
       ALLOCATE(DINTER_T0(HNO,HNO,NMOL1,NMOL1))
       ALLOCATE(DINTER_T(HNO,HNO,NMOL1,NMOL1))
       ALLOCATE(DINTRA_T0(HNO,HNO,NMOL1))
       ALLOCATE(DINTRA_T(HNO,HNO,NMOL1))
       ALLOCATE(THETA2_T0(HNO,HNO,NMOL1,NMOL1))
       ALLOCATE(THETA2_T(HNO,HNO,NMOL1,NMOL1))
       ALLOCATE(THETA1_T0(HNO,HNO,NMOL1))
       ALLOCATE(THETA1_T(HNO,HNO,NMOL1))
       ALLOCATE(MOLECULE1(3,HNO,NMOL1,TSIM))
       ALLOCATE(DINTER(HNO,HNO,NMOL1,NMOL1))
       ALLOCATE(DINTRA(HNO,HNO,NMOL1))
       ALLOCATE(AVE_INTRA_HH(TWIN))         
       ALLOCATE(AVE_INTER_HH(TWIN))
       ALLOCATE(THETA2(HNO,HNO,NMOL1,NMOL1))
       ALLOCATE(THETA1(HNO,HNO,NMOL1))
       ALLOCATE(POSITIONS1_T0(3,HNO,NMOL1))
       ALLOCATE(POSITIONS1_T(3,HNO,NMOL1))
      
       ALLOCATE(BOX(3,2))
        
       AVE_INTRA_HH = 0.0
       AVE_INTER_HH = 0.0
       
       
       CALL CPU_TIME(START_DUMP)
       
       CALL READ_DUMP(DUMPFILE,TSIM,NATOMS,NTYPES,&
            &MOLECULE1,MOL1,NMOL1,HNO,BOX)
       
       
       
       ! print*, 'AFTER', MOLECULE1(1,16,1,1), MOLECULE1(2,16,1,1), &
        !                MOLECULE1(3,16,1,1)

       NWIN=INT(TSIM/TWIN)
       PRINT*, "NWIN: " , NWIN       

       LX=BOX(1,2)-BOX(1,1)
       LY=BOX(2,2)-BOX(2,1)
       LZ=BOX(3,2)-BOX(3,1)
       PRINT*, " "
       PRINT*, "LX,LY,LZ: ", LX, LY,LZ
       PRINT*, " "


 
       PRINT*, " "
       PRINT*, " PROPERTIES CALCULATION"
       PRINT*, "========================"

       PRINT*, " "
       PRINT*, " "


       !===============================================================
       DO TT=0,NWIN-1
       
          T0=1+TWIN*TT     ! TRANSLATION TIME ORIGEN  
          !$OMP PARALLEL DO PRIVATE(I,J,K) 
          DO I=1,NMOL2
             DO J=1,HNO
             DO K=1,3
                   POSITIONS1_T0(K,J,I)=MOLECULE1(K,J,I,T0)
             ENDDO
             END DO
          END DO
          !$OMP END PARALLEL DO
        
          CALL  DISTANCE(DINTRA,DINTER,THETA1,THETA2,&
                & POSITIONS1_T0,NMOL2,HNO,BOX)
            
          
          DINTRA_T0 = DINTRA
          THETA1_T0 = THETA1
          DINTER_T0 = DINTER
          THETA2_T0 = THETA2

       
          DO TCORR=1,TWIN
             T=TCORR+TWIN*TT
             !$OMP PARALLEL DO PRIVATE(I,J,K)
             DO I=1,NMOL2
                DO J=1,HNO
                   DO K=1,3
                      POSITIONS1_T(K,J,I)=MOLECULE1(K,J,I,T)
                   ENDDO
                ENDDO
             ENDDO
             !$OMP END PARALLEL DO
        
             CALL  DISTANCE(DINTRA,DINTER,THETA1,THETA2,&
                   & POSITIONS1_T,NMOL2,HNO,BOX)
                
              DINTRA_T = DINTRA                  
              THETA1_T = THETA1  
              DINTER_T = DINTER
              THETA2_T = THETA2
         
             
             CALL CORRELATION(G_INTRA,G_INTER,DINTRA_T0,&
                  & DINTRA_T, DINTER_T0, DINTER_T,THETA1_T0,&
                  & THETA2_T0, THETA1_T, THETA2_T,NMOL2,HNO)
              
             AVE_INTRA_HH(TCORR)=AVE_INTRA_HH(TCORR)+2*G_INTRA
             AVE_INTER_HH(TCORR)=AVE_INTER_HH(TCORR)+2*G_INTER
             
          END DO ! END LOOP TCORR
       END DO ! END LOOP TIME ORIGIN 
       !===============================================================         
       PRINT*, " "
       PRINT*, " " 


       DO TCORR = 1,TWIN
          AVE_INTRA_HH(TCORR)=AVE_INTRA_HH(TCORR)/((HNO-1)*NMOL2*NWIN)
          AVE_INTER_HH(TCORR)=AVE_INTER_HH(TCORR)/&
                             (HNO*(NMOL2-1)*NWIN)
       END DO   
       
       CALL WRITE_T2(OUTFILE,AVE_INTRA_HH, AVE_INTER_HH,&
       &    TWIN,NWIN,DT,NMOL2) 
       
       
       DEALLOCATE(DINTER_T0)
       DEALLOCATE(DINTER_T)
       DEALLOCATE(DINTRA_T0)
       DEALLOCATE(DINTRA_T)
       DEALLOCATE(THETA2_T0)
       DEALLOCATE(THETA2_T)
       DEALLOCATE(THETA1_T0)
       DEALLOCATE(THETA1_T)
       DEALLOCATE(MOLECULE1)
       DEALLOCATE(DINTER)
       DEALLOCATE(DINTRA)
       DEALLOCATE(AVE_INTRA_HH)
       DEALLOCATE(AVE_INTER_HH)
       
       END PROGRAM T2_OCTANOIC  
 
