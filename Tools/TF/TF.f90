PROGRAM TRANSFOFOURIER

IMPLICIT NONE

INTEGER*8 :: i,Norb,Orbital,Ngrid_Omega,Ngrid_time,Ngrid_r
REAL*8 :: step_Omega,step_time,t,Real_Part,Imaginary_Part,Omega,blank1,blank2
COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: ARRAY
COMPLEX*16, DIMENSION(500000) :: TF

OPEN(10,file='Dipole_Norb=001.dat')
OPEN(20,file='TF.dat')

Norb=1
Orbital=1

Ngrid_Omega=100000
step_Omega=1.D-4
Ngrid_time=4000
step_time=0.1D0

ALLOCATE(ARRAY(Ngrid_time,1))

!Imaginary_Part=0.D0

DO i=1,Ngrid_time

   READ(10,*)t,Real_Part,Imaginary_Part
   ARRAY(i,1)=DCMPLX(Real_Part,Imaginary_Part)

END DO

CALL DFT_Time(Ngrid_Omega,step_Omega,Ngrid_time,step_time,ARRAY,Norb,Orbital,&
     &TF)

DO i=1,Ngrid_Omega

   Omega=REAL(i)*step_Omega

   WRITE(20,*)Omega,ABS(TF(i))**2 ! Fourier transform & HHG

END DO

END PROGRAM TRANSFOFOURIER

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$                                                                            $
!$                                                                            $
!$                           SUBROUTINE DFT_Time                              $
!$                                                                            $
!$                                                                            $
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!==============================================================================
!               
! This subroutines calculates the Discrete Fourier Transform of an array whose
! values are expressed with respect to time.
!
!==============================================================================

SUBROUTINE DFT_Time(N_Omega,step_Omega,Ngrid_time,step_time,ARRAY,Norb,Orbital,&
     &TF)

IMPLICIT NONE

! Input/Output variables

INTEGER*8 :: N_Omega,Ngrid_time,Norb,Orbital
REAL*8 :: step_Omega,step_time
COMPLEX*16, DIMENSION(N_Omega) :: TF
COMPLEX*16, DIMENSION(Ngrid_time,Norb) :: ARRAY

! Local variables

INTEGER :: i,j
REAL*8 :: Omega

REAL*8, PARAMETER :: PI=3.14159265358979D0

TF=DCMPLX(0.D0,0.D0)

DO i=1,N_Omega

   Omega=REAL(i)*step_Omega

   DO j=1,Ngrid_time

      ! Compute DFT with basic integratin scheme

      TF(i)=TF(i)+ARRAY(j,Orbital)*EXP(DCMPLX(0.D0,-Omega*REAL(j)*step_time))

   END DO

   TF(i)=TF(i)/(SQRT(2.D0*PI)*Ngrid_time) ! Normalizing the Fourier Transform 

END DO

RETURN

END SUBROUTINE DFT_Time
