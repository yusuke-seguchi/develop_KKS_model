! ===================================================
! PHASE-FIELD CALCULATION PROGRAM FOR BINARY ALLOY
! FOR FREE DNDRITE ANALYSIS
! 大出真知子、鈴木俊夫、”フェーズフィールドシミュレーション入門”、鋳造工学　第73巻（2001）第5号
! ===================================================

MODULE COMMON_V 
    INTEGER::L=0,I=1,J=1,M,N,MC,NC 
    INTEGER::TRIANGLE=15,MODSAVE=2000,LSAVE=0 
    REAL(8)::XM,EP,EP2,W 
    REAL(8)::DX,DY,DT,DT1,DS,DL,DL2,DS2 
    REAL(8)::CLE,CSE,TEMPMELT,TEMP,VM,C0 
    REAL(8)::V,YK,SIGMA,XME,KE,BETA,COMNOISE 

    REAL(8),ALLOCATABLE:: PHI(:,:,:),COM(:,:,:) 
    REAL(8),ALLOCATABLE:: FCL(:,:),FCC(:,:),CL(:,:),CS(:,:) 
    REAL,PARAMETER::R=8.314
END MODULE COMMON_V

! ================= 
! MAIN PROGRAM
! =================
PROGRAM MAIN 

    USE COMMON_V
    IMPLICIT NONE 
    REAL(8)::PD,PHIX,PHIY,PHIXX,PHIYY,PHIXY 
    REAL(8)::E1,E2,E3,E4,E5,TH,ETA,DETA 
    REAL(8)::P,PP,PG,GD,GG,FP,DC,C,DPHI 
    REAL(8)::FCCW,FCCE,FCCN,FCCS,FCCL 
    REAL(8)::XJ1,XJ2,XJ3,XJ4 
    REAL(8)::D1,D2,D3,D4,D5 
    REAL(8)::A,AA,BB,CC 

! --------------------
!  READ CALCULATION CONDITION 
! --------------------
    CALL CAL_COND 

! ------------------
!  SET THE INITIAL PHASE-CONDITION 
! ------------------
    CALL INIT_COND
    CALL OUTSAVE 

    ALLOCATE( CS(0:M+1,0:N+1) )
    ALLOCATE( CL(0:M+1,0:N+1) ) 
    ALLOCATE( FCL(0:M+1,0:N+1) )
    ALLOCATE( FCC(0:M+1,0:N+1) ) 

! -------------------
!  SET THE PHASE-FIELD PARAMETERS (& TIMESTEP) 
! ------------------
    EP=SQRT(18./2.2*DX*SIGMA)
    EP2=EP*EP 
    W=2.2*SIGMA/DX 
    CALL MOBILITY

    A=CLE/CSE*(1-CSE)/(1-CLE) 

    WRITE(6,*)'SET ALL THE CALCULATION CONDITIONS' 
    WRITE(6,*)'NOW CALCULATING ....... ' 

! =================
! CALCULATE THE GOVERNING EQUATIONS
! ================= 

500 L=L+1

! ----------------
! CALCULATE CS &CL 
! -----------------
    DO I=0, MC+1 
    DO J=0, NC+1 
         P=PHI(0, I, J) 
         PP=P**3*(10.-15.*P+6.*P*P) 
    
        IF (PHI(0, I, J).LT.0.001) THEN 
         CL(I, J) = COM(0, I, J) 
         CS(I, J) = CL(I, J) / (A + (1. - A)*CL(I, J)) 
        
        ELSE IF(PHI(0, I, J).GT.0.999) THEN 
         CS(I, J) = COM(0, I, J) 
         CL(I, J) = A*CS(I, J) / (1. + (A - 1.)*CS(I, J)) 
        
        ELSE 
         AA = PP*(1. -A) 
         BB = PP + (1.-PP)*A+COM(0, I, J)*(1-A)
         CC = COM(0, I, J)

         CS(I, J) = (BB - SQRT(BB*BB - 4.*AA*CC))/(2.*AA) 
         CL(I, J) = (CC - PP*CS(I, J)) / (1. -PP) 
        END IF

        FCL(I, J) = R*TEMP/VM * LOG( CL(I, J)/(1. -Cl(I, J)) ) 
        FCCL = R*TEMP/VM/(CL(I, J) * (1. - CL(I, J))) 
        FCCS = R*TEMP/VM/(CS(I, J) * (1. - CS(I, J))) 
        FCC(I, J) = FCCL*FCCS/((1.-PP)*FCCS + PP*FCCL) 
        
    END DO 
    END DO 

! ------------------
!  GOVERNING EQUATIONS 
! ------------------

    DO I=1, MC 
    DO J=1, NC 
        P = PHI(0, I, J) ; C = COM(0, I, J)

! ------------------
!  TIME SAVING 
! ------------------
        PD = (PHI(0, I+1, J) + PHI(0, I-1, J) + PHI(0, I, J+1) + PHI(0, I, J-1))/4. 
        IF(PD.LE.1.E-5) THEN 
            DPHI = 0. 
            DC = DL*( (COM(0, I, J+1) + COM(0, I, J-1) - 2.*C)/(DY*DY) + ( (COM(0, I+1, J) + COM(0, I-1, J) - 2.*C))/(DX*DX)) 
            
        ELSE IF(PD.GE.1.-1.E-5)THEN
            DPHI=0. 
            DC = DS*( (COM(0, I, J+1) + COM(0, I, J-1) - 2.*C)/(DY*DY) + ( (COM(0, I+1, J) + COM(0, I-1, J) - 2.*C))/(DX*DX) ) 

! --------------------
!  NON-TIME SAVING 
! --------------------
        ELSE 
            PG = 30*P*P*(1-P)*(1-P) 
            GD = 2.*P*(1.-P)*(1.-2.*P)

            GG = PG*LOG( (1. - CSE)/ (1. - CLE) * (1. - CL(I, J)) / (1. - CS(I, J))) 

            FP = R*TEMP/VM*GG - W*GD

        PHIX = (PHI(0, I-1, J) - PHI(0, I+1, J))/(2.*DX) 
        PHIY = (PHI(0, I, J-1) - PHI(0, I, J+1))/(2.*DY) 
        PHIXX = (PHI(0, I-1, J) + PHI(0, I+1, J) -2.*P) / (DX*DX) 
        PHIYY = (PHI(0, I, J-1) + PHI(0, I, J+1) -2.*P) / (DY*DY) 
        PHIXY = (PHI(0, I+1, J+1) + PHI(0, I-1, J-1) - PHI(0, I-1, J+1) - PHI(0, I+1, J-1))/(2.*DX*2.*DY) 
        TH = ATAN(PHIY/(PHIX + 1.E-20)) 
        ETA = 1. + V*COS(YK* TH) 
        DETA = V*YK*SIN(YK*TH) 
        E1 = EP2*ETA*ETA*(PHIXX + PHIYY) 
        E2 = EP2*ETA*(-DETA)*(SIN(2.*TH)*(PHIYY-PHIXX)+2.*COS(2.*TH)*PHIXY) 
        E3 = 0.5*EP2 
        E4 = DETA*DETA + ETA*(-V*YK*YK*COS(YK*TH)) 
        E5 = 2.*SIN(2.*TH)*PHIXY - PHIXX - PHIYY - COS(2.*TH)*(PHIYY-PHIXX) 
        DPHI = XM*(E1 + E2 - E3*E4*E5 + FP)

        D1 = DS; D2 = DS; D3 = DS; D4 = DS; D5 = DS

        IF(P.LE.0.9) D1 = DL 
        IF(PHI(0, I-1, J).LE.0.9) D2 = DL; IF(PHI(0, I+1, J).LE.0.9) D3 = DL 
        IF(PHI(0, I, J+1).LE.0.9) D4 = DL; IF(PHI(0, I, J-1).LE.0.9) D5 = DL 

        FCCW = 2.*D1/FCC(I,J)*D2/FCC(I-1,J)/(D1/FCC(I,J)+D2/FCC(I-1,J)) 
        FCCE = 2.*D1/FCC(I,J)*D3/FCC(I+1,J)/(D1/FCC(I,J)+D3/FCC(I+1,J)) 
        FCCS = 2.*D1/FCC(I,J)*D4/FCC(I,J+1)/(D1/FCC(I,J)+D4/FCC(I,J+1)) 
        FCCN = 2.*D1/FCC(I,J)*D5/FCC(I,J-1)/(D1/FCC(I,J)+D5/FCC(I,J-1))

        XJ1 = (-FCL(I,J)+FCL(I-1,J))/DX*FCCW 
        XJ2 = (-FCL(I,J)+FCL(I+1,J))/DX*FCCE 
        XJ3 = (-FCL(I,J)+FCL(I,J+1))/DY*FCCS 
        XJ4 = (-FCL(I,J)+FCL(I,J-1))/DY*FCCN 

        DC = (XJ1+XJ2)/DX+(XJ3+XJ4)/DY

        END IF 

            PHI(1,I,J) = P + DPHI*DT; COM(1,I,J) = C + DC*DT 

    END DO 
    END DO 

!--------------------
!  END GOVERNING EQUATION CALCULATIONS 
!--------------------
!  BOUDARY CONDITONS 
!--------------------

    DO I = 0, M+1 
        PHI(1, I, 0) = PHI(1, I, 1); PHI(1, I, N+1) = PHI(1, I, N) 
        COM(1, I, 0) = COM(1, I, 1); COM(1, I, N+1) = COM(1, I, N) 
    END DO

    DO J = 0, N+1 
        PHI(1, 0, J) = PHI(1, 1, J); PHI(1, M+1, J) = PHI(1, M, J) 
        COM(1, 0, J) = COM(1, 1, J); COM(1, M+1, J) = COM(1, M, J) 
    END DO

!----------------
!  RENEWAL OF PHASE &CONCENTRATION FIELDS 
!----------------
    DO I=0,M+1 
    DO J=0,N+1
        PHI(0, I, J)=PHI(1, I, J); COM(0, I, J)=COM(1, I, J) 
    END DO 
    END DO 

!------------------
!  NOISE 
!------------------
    CALL NOISE 
!------------------
!  AREA SET FOR TIME SAVING 
!------------------
    IF(L.GE.100) CALL AREASET 
!--------------------- 
!  OUT PUT 
!------------------
    IF(MOD(L,MODSAVE).EQ.0) CALL OUTSAVE 
!-----------------
!  END CONDITION 
!------------------
    IF(PHI(0, 1, N-10).LE.0.5) GOTO 500 
    WRITE(6,*)'CALCULATION HAS FINISHED!' 
    
END PROGRAM MAIN 

!=============== 
!  SUBRUTINE 
!================ 
!  READ CALCULATION CONDITION 
!--------------------------
SUBROUTINE CAL_COND 

    USE COMMON_V
    IMPLICIT NONE 
    M = 750 ! (X-0 IRECT I ON MESH NUIIBER) 
    N = 750 !CY-DIRECTION MESH NUMBER) 
    DX = 1.e-8 !MESH SIZE 
    DY = 1.e-8 !MESH SIZE 
    DL = 3.e-9 !DL 
    DS = 3.e-13 !DS 
    C0 = 0.0196 !INITIAL SOLUTE CONTENT 
    TEMP = 900.0 !INITIAL TEMPERATURE 
    TEMPMELT = 933.3 !MELTING POINT 
    VM = 10.547e-6 !MOLER VOLUME 
    XME = 640.0 !LIOUIDUS SLOPE 
    KE = 0.14 !PARTITION COEFFICIENT 
    BETA = 0.0 !KINETIC COEFFICIENT 
    V=0.03 !ANISOTROPY EP=EP(l+V~OS(YK*THETA) 
    YK=4.0 !ANISOTROPY 
    SIGMA = 0.093 !INTERFACE ENERGY 
    COMNOISE = 0.0 !NOISE 

    NC = N; MC = M; DL2 = 2.0*DL/DX; DS2 = 2.0*DS/DX 
    CLE = (TEMPMELT - TEMP)/XME; CSE = CLE*KE 

END SUBROUTINE CAL_COND

! ---------------------
!  INIT_COND 
!-------------------------
SUBROUTINE INIT_COND 

    USE COMMON_V 
    IMPLICIT NONE

    ALLOCATE( PHI(0:1,0:M+1,0:N+1) ) 
    ALLOCATE( COM(0:1,0:M+1,0:N+1) ) 
    
    DO I=1,M 
        DO J=1,N 
            PHI(0, I, J) = 0.; PHI(1, I, J) = 0.
            COM(0, I, J) = C0; COM(1, I, J) = C0
        END DO 
    END DO 

    DO I=1, TRIANGLE+5 
        DO J=1, TRIANGLE+5 
            IF(J.LT.(-I+TRIANGLE))THEN 
                PHI(0, I, J) = 1.; PHI(1, I, J) = 1. 
                COM(0, I, J) = CSE; COM(1, I, J) = CSE 
            END IF 

            IF(J.EQ.(-I+TRIANGLE))THEN 
                PHI(0, I, J) = 0.5; PHI(1, I, J) = 0.5 
            END IF

        END DO 
    END DO 
     
    DO I = 0, M+1 
        PHI(0, I, 0) = PHI(0, I, 1); PHI(0, I, N+1) = PHI(0, I, N) 
        COM(0, I, 0) = COM(0, I, 1); COM(0, I, N+1) = COM(0, I, N) 
    END DO 
     
    DO J = 0, N+1 
        PHI(0, 0, J) = PHI(0, 1, J); PHI(0, M+1, J)= PHI(0, M, J) 
        COM(0, 0, J) = COM(0, 1, J); COM(0, M+1,J) = COM(0, M, J) 
    END DO

    RETURN 

END SUBROUTINE INIT_COND

!-------------------
!  PF MOBILITY 
! ------------------
SUBROUTINE MOBILITY

    USE COMMON_V 
    IMPLICIT NONE

    REAL(8)::P1,P2,PP1,PP2,FUN1,FUN2,ZETA,FCCLE,FCCSE,ALPHA 
    
    EP2 = EP*EP; ZETA = 0.

    FCCLE = R*TEMP/VM/(CLE*(1.-CLE)) 
    FCCSE = R*TEMP/VM/(CSE*(1.-CSE))

    DO P1 = 0.001, 0.998, 0.001 
        P2 = P1 + 0.001 
        PP1 = P1**3*(10. - 15.*P1 + 6.*P1*P1) 
        PP2 = P2**3*(10. - 15.*P2 + 6.*P2*P2) 
        FUN1 = PP1*(1-PP1)/( (1-PP1)*fCCSE+PP1*fCCLE )/(P1*(1.-P1)) 
        FUN2 = PP2*(1-PP2)/( (1-PP2)*FCCSE+PP1*FCCLE )/(P2*(1.-P2)) 
        ZETA = ZETA+(FUN1 + FUN2)*0.001/2. 
    END DO 

    ALPHA = BETA*R*TEMP*(1-KE) / (VM*XME) 
    XM = 1./(EP*EP/SIGMA*(ALPHA + EP/(DL*SQRT(2.*W))*ZETA*FCCSE*FCCLE*(CLE-CSE)**2 )) 
    DT = DX**2/(5.*XM*EP**2) 
    DT1 = DX**2/(5.*DL) 
    DT = DMIN1(DT,DT1)

    RETURN

END SUBROUTINE MOBILITY

! ------------------
!  COIINOISE 
! ------------------
SUBROUTINE NOISE

    USE COMMON_V 
    IMPLICIT NONE

    REAL(8)::COUNTNOISE,COMTOT,COMDAM,CNOISE
    INTEGER::LAM=12869,C=6925,MMU=32768,X=19724 

    COUNTNOISE=0. ; COMTOT=0. 

    DO I = 1, MC 
        DO J = 1, NC 

            IF(PHI(0, I, J).GT.0.01.AND.PHI(0, I, J).LE.0.5)THEN 
                COMDAM = COM(0, I, J)

                X = MOD( (X*LAM + C), MMU) 
                CNOISE = (REAL(X)/MMU-0.5)*COMNOISE

                COM(0, I, J) = COMDAM*(1. + CNOISE) 
                COMTOT = COMTOT + COMDAM*CNOISE 
                COUNTNOISE = COUNTNOISE + 1. 
            END IF 

        END DO 
    END DO

    DO I=1, MC 
        DO J=1, NC 

            IF (PHI(0, I, J).GT.0.01 .AND. PHI(0, I, J).LE.0.5) THEN 
                COM(0, I, J) = COM(0, I, J) - (COMTOT/COUNTNOISE) 
            END IF 

        END DO 
    END DO 

    DO I = 0, M+1 
        COM(0, I, 0) = COM(0, I, 1); COM(0, I, N+1) = COM(0, I, N) 
    END DO

    DO J = 0, N+1 
        COM(0, 0, J) = COM(0, 1, J); COM(0, M+1, J) = COM(0, M, J) 
    END DO

    RETURN 

END SUBROUTINE NOISE 

!--------------------
!  AREA SET 
!--------------------
SUBROUTINE AREASET

    USE COMMON_V
    IMPLICIT NONE

    DO J = 1, N 
        IF (ABS (COM(0, 1, J)/C0 - 1.).GT.1.E-5) NC = J 
    END DO

    DO I = 1, M
        IF (ABS (COM(0, I, 1)/C0 - 1.).GT.1.E-5) MC = I 
    END DO

    NC = NC + 10; MC = MC + 10
    IF (NC.GT.N) NC = N; IF (MC.GT.M) MC = M

    RETURN 

END SUBROUTINE AREASET
 

!---------
!  OUTSAVE 
!----------------
SUBROUTINE OUTSAVE 

    USE COMMON_V 
    IMPLICIT NONE 

    CHARACTER*3::OUT_NUM 
    CHARACTER*30::FPOUT 
    INTEGER::ONE,TEN,HAND

    LSAVE = LSAVE + 1 
    ONE = MOD(LSAVE, 10) 
    TEN = MOD(INT(REAL(LSAVE)/10.), 10) 
    HAND = MOD(INT(REAL(LSAVE)/100.), 10)

    ONE = 48 + ONE; TEN = 48 + TEN; HAND = 48 + HAND; 
    OUT_NUM = CHAR(HAND)//CHAR(TEN)//CHAR(ONE)

    FPOUT = "./calc_result/output"//OUT_NUM//".vtk" 
    OPEN(101, FILE=FPOUT, ERR=1000)
    write(101,'(a)') '# vtk DataFile Version 3.0'
    write(101,'(a)') 'output.vtk'
    write(101,'(a)') 'ASCII'
    write(101,'(a)') 'DATASET STRUCTURED_POINTS'
    write(101,'(a,3i5)') 'DIMENSIONS',M,N,1
    write(101,'(a,3f4.1)')'ORIGIN' ,0.0,0.0,0.0
    write(101,'(a,3i2)')'ASPECT_RATIO',1,1,1
    write(101,'(a,1i11)')'POINT_DATA',M*N*1
    write(101,'(a)')'SCALARS concentration double'
    write(101,'(a)')'LOOKUP_TABLE default' 

    DO J = 1, N
        DO I = 1, M 
            WRITE (101,*) COM(0, I, J) 
        END DO 
    END DO

    write(101,'(a)')'SCALARS phase_field double'
    write(101,'(a)')'LOOKUP_TABLE default'

    do J = 1, N
        do I = 1, M
            WRITE(101,*) PHI(0, I, J)
        END DO
    END DO

    CLOSE(101)
    WRITE(6,*) "SAVING! "//OUT_NUM//"STEPS"

    RETURN

410 FORMAT(I4, I4, I10) 
300 FORMAT(I5, I5, E12.5 ,E12.5) 
1000 WRITE(6,*)'ERROR IN FILE OPEN'

END SUBROUTINE OUTSAVE




!---------
!  OUTSAVE 
!----------------
! SUBROUTINE OUTSAVE 

!     USE COMMON_V 
!     IMPLICIT NONE 

!     CHARACTER*3::OUT_NUM 
!     CHARACTER*15::FPOUT 
!     INTEGER::ONE,TEN,HAND

!     LSAVE = LSAVE + 1 
!     ONE = MOD(LSAVE, 10) 
!     TEN = MOD(INT(REAL(LSAVE)/10.), 10) 
!     HAND = MOD(INT(REAL(LSAVE)/100.), 10)

!     ONE = 48 + ONE; TEN = 48 + TEN; HAND = 48 + HAND; 
!     OUT_NUM = CHAR(HAND)//CHAR(TEN)//CHAR(ONE)

!     FPOUT = "output"//OUT_NUM//".dat" 
!     OPEN(14, FILE=FPOUT, ERR=1000) 

!     WRITE(14, 410) MC,NC,L

!     DO I = 1, MC 
!         DO J = 1, NC 
!             WRITE (14, 300) I, J, PHI(0, I, J), COM(0, I, J) 
!         END DO 
!     END DO

!     CLOSE(14)

!     RETURN

! 410 FORMAT(I4, I4, I10) 
! 300 FORMAT(I5, I5, E12.5 ,E12.5) 
! 1000 WRITE(6,*)'ERROR IN FILE OPEN'

! END SUBROUTINE OUTSAVE


! OUTSAVE
! subroutine outsave
!     use COMMON_V
!     implicit none

!     character*3 :: out_num
!     character*15 :: fpout
!     integer :: one, ten, hand

    
!     one = mod(lsave, 10)
!     ten = mod(int(real(lsave)/10.), 10)
!     hand = mod(int(real(lsave)/100.), 10)

!     one = 48 + one; ten = 48 + ten; hand = 48 + hand;
!     out_num = char(hand)//char(ten)//char(one)
    
!     fpout = "output"//out_num//".vtk"
!     open(101,file=fpout, err=1000)
!     write(101,'(a)') '# vtk DataFile Version 3.0'
!     write(101,'(a)') 'output.vtk'
!     write(101,'(a)') 'ASCII'
!     write(101,'(a)') 'DATASET STRUCTURED_POINTS'
!     write(101,'(a,3i5)') 'DIMENSIONS',m,n,1
!     write(101,'(a,3f4.1)')'ORIGIN' ,0.0,0.0,0.0
!     write(101,'(a,3i2)')'ASPECT_RATIO',1,1,1
!     write(101,'(a,1i11)')'POINT_DATA',m*n*1
!     write(101,'(a)')'SCALARS concentration double'
!     write(101,'(a)')'LOOKUP_TABLE default'

!     do j=1,n
!         do i=1,m
!             write(101,*) com(0,i,j)
!         end do
!     end do

!     write(101,'(a)')'SCALARS phase_field double'
!     write(101,'(a)')'LOOKUP_TABLE default'

!     do j=1,n
!         do i=1,m
!             write(101,*) phi(0,i,j)
!         end do
!     end do

!     close(101)
!     return

! ! 300     format(e12.5)
! 1000    write(6,*)'ERROR IN FILE OPEN'

! end subroutine outsave