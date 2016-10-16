MODULE MATHLIB

USE EXPOKIT

IMPLICIT NONE

CONTAINS

!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!%                                                                            $
!%                                                                            $
!%                           SUBROUTINE DFT_Time2                             $
!%                                                                            $
!%                                                                            $
!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!==============================================================================
!               
! This subroutines calculates the Discrete Fourier Transform of an array whose
! values are expressed with respect to time. This version is for the general
! dipole.
!
!==============================================================================

SUBROUTINE DFT_Time2(N_Omega,step_Omega,Ngrid_time,step_time,ARRAY,TF2)

IMPLICIT NONE

! Input/Output variables

INTEGER*8 :: N_Omega,Ngrid_time
REAL*8 :: step_Omega,step_time
COMPLEX*16, DIMENSION(N_Omega) :: TF2
COMPLEX*16, DIMENSION(Ngrid_time+1) :: ARRAY

! Local variables

INTEGER :: i,j
REAL*8 :: Omega

REAL*8, PARAMETER :: PI=3.14159265358979D0

TF2=DCMPLX(0.D0,0.D0)

DO i=1,N_Omega

   Omega=REAL(i)*step_Omega

   DO j=1,Ngrid_time

      ! Compute DFT with basic integratin scheme

      TF2(i)=TF2(i)+ARRAY(j)*EXP(DCMPLX(0.D0,-Omega*REAL(j)*step_time))

   END DO

   TF2(i)=TF2(i)/(SQRT(2.D0*PI)*Ngrid_time) ! Normalizing the Fourier Transform 

END DO

RETURN

END SUBROUTINE DFT_Time2

!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!%                                                                            $
!%                                                                            $
!%                           SUBROUTINE DFT_Time                              $
!%                                                                            $
!%                                                                            $
!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

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
COMPLEX*16, DIMENSION(Ngrid_time+1,Norb) :: ARRAY

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

  subroutine propexact(nr,wf,nord,dt,h0matrix,h0matrix1,h0matrix2,h1matrix,&
       &h2matrix,h2matrix1,h2matrix2,lmin,lmax,lmaxcoup)
    
    implicit none
    integer*8 :: nr,lmin,lmax,iord,nord,lmaxcoup,i
    double precision,intent(in) :: dt
    double complex,intent(inout) :: wf(1:nr,0:lmax)
    double complex,dimension(1:size(wf,1),0:size(wf,2)-1) :: wf1,wf2,wf3
    double complex :: xfac,xdt
    

    double precision :: h0matrix(1:nr,0:lmax)
    double precision :: h1matrix(1:nr,0:lmax,0:lmax)
    double precision :: h0matrix1,h0matrix2

    !double complex :: h2matrix(1:nr,0:lmax)
    !double complex :: h2matrix1(0:lmax)
    !double complex :: h2matrix2(0:lmax)

    double precision :: h2matrix(1:nr,0:lmax)
    double precision :: h2matrix1(0:lmax)
    double precision :: h2matrix2(0:lmax)

    xdt=-dcmplx(0.d0,1.d0)*dt

    wf1=wf
    xfac=1
    iord=0

    do while(iord<nord)

      iord=iord+1
      call velmatmul(nr,wf1,wf2,h0matrix,h0matrix1,h0matrix2,h1matrix,&
       &h2matrix,h2matrix1,h2matrix2,lmin,lmax,lmaxcoup)
      wf1=wf2
      xfac=xfac*xdt/iord
      wf3=wf+xfac*wf1

      iord=iord+1
      call velmatmul(nr,wf1,wf2,h0matrix,h0matrix1,h0matrix2,h1matrix,&
       &h2matrix,h2matrix1,h2matrix2,lmin,lmax,lmaxcoup)
      wf1=wf2
      xfac=xfac*xdt/iord
      wf=wf3+xfac*wf1
    end do

    !do i=1,nr
    !print*,wf(i,0)
    !end do
    !stop

    return
    !write(*,'("wf",f6.2,2es15.6)')(wf(i,0),i=1,nr)
    !stop

  end subroutine propexact

  subroutine velmatmul(nr,wf1,wf2,h0matrix,h0matrix1,h0matrix2,h1matrix,&
       &h2matrix,h2matrix1,h2matrix2,lmin,lmax,lmaxcoup)
    
    implicit none

    integer*8 :: nr,lmin,lmax,i,j,l,l1,l2

    double precision :: h0matrix(1:nr,0:lmax)
    double precision :: h1matrix(1:nr,0:lmax,0:lmax)
    double precision :: h0matrix1,h0matrix2

    !double complex :: h2matrix(1:nr,0:lmax)
    !double complex :: h2matrix1(0:lmax)
    !double complex :: h2matrix2(0:lmax)

    double precision :: h2matrix(1:nr,0:lmax)
    double precision :: h2matrix1(0:lmax)
    double precision :: h2matrix2(0:lmax)

  
    double complex :: wf1(1:nr,0:lmax)
    double complex :: wf2(1:nr,0:lmax)

    integer*8 :: n1,n2,n3,lmaxcoup

    n1=nr-1
    n2=nr-2
    n3=nr-3

    !do i=1,nr
    !print*,'wf1',wf1(i,3)
    !end do

    !read(*,*)

    do l=lmin,lmax
      wf2(1:nr,l)=h0matrix(1:nr,l)*wf1(1:nr,l)
      wf2(1:n2,l)=wf2(1:n2,l)+h0matrix2*wf1(3:nr,l)
      wf2(1:n1,l)=wf2(1:n1,l)+h0matrix1*wf1(2:nr,l)
      wf2(2:nr,l)=wf2(2:nr,l)+h0matrix1*wf1(1:n1,l)
      wf2(3:nr,l)=wf2(3:nr,l)+h0matrix2*wf1(1:n2,l)
    end do

    do l=lmin,lmax
      j=l-1
      if(j>=0)then
        wf2(1:nr,j)=wf2(1:nr,j)+h2matrix(1:nr,j)*wf1(1:nr,l)
        !wf2(1:n2,j)=wf2(1:n2,j)+h2matrix2(j)*wf1(3:nr,l)
        !wf2(1:n1,j)=wf2(1:n1,j)+h2matrix1(j)*wf1(2:nr,l)
        !wf2(2:nr,j)=wf2(2:nr,j)-h2matrix1(j)*wf1(1:n1,l)
        !wf2(3:nr,j)=wf2(3:nr,j)-h2matrix2(j)*wf1(1:n2,l)
      end if
      j=l+1
      if(j<=lmax)then
        wf2(1:nr,j)=wf2(1:nr,j)+h2matrix(1:nr,l)*wf1(1:nr,l)
        !wf2(1:nr,j)=wf2(1:nr,j)-h2matrix(1:nr,l)*wf1(1:nr,l)
        !wf2(1:n2,j)=wf2(1:n2,j)+h2matrix2(l)*wf1(3:nr,l)
        !wf2(1:n1,j)=wf2(1:n1,j)+h2matrix1(l)*wf1(2:nr,l)
        !wf2(2:nr,j)=wf2(2:nr,j)-h2matrix1(l)*wf1(1:n1,l)
        !wf2(3:nr,j)=wf2(3:nr,j)-h2matrix2(l)*wf1(1:n2,l)
      end if
    end do

    !do i=1,nr
    !print*,'wf2',wf2(i,1)
    !end do

    !read(*,*)


    do l1=lmin,lmax
      do l2=lmin,lmax
        j=abs(l1-l2)
        if(0<j.AND.j<=lmaxcoup)then
          do i=1,nr
            wf2(i,l1)=wf2(i,l1)+h1matrix(i,l1,l2)*wf1(i,l2)
            !print*,'h1,wf1',h1matrix(i,l1,l2),wf1(i,l2)
          end do
          !read(*,*)
        end if
      end do
    end do

    !do i=1,nr
    !   print*,'wf2',wf2(i,4)
    !end do

    !read(*,*)

  end subroutine velmatmul


!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!%                                                                            $
!%                                                                            $
!%                            SUBROUTINE TAYLOR                               $
!%                                                                            $
!%                                                                            $
!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!==============================================================================
!               
! This subroutine calculates the expansion of the exponential of a matrix, up
! to the required order. It is suited for the block diagonal structure of the
! Hamiltonian.
!
!==============================================================================

SUBROUTINE TAYLOR(Ngrid_r,Norb,Orbital,lmax,step_time,m_Taylor,&
     &H_DIAG,H_DIAG1,H_DIAG2,H_NONDIAG,H_LASER,H_LASER1,&
     &H_LASER2,RADIAL,W)

IMPLICIT NONE

! Input/Output variables declaration

INTEGER*8 :: Ngrid_r,Norb,Orbital,lmax,m_Taylor
REAL*8 :: step_time,H_DIAG1,H_DIAG2
REAL*8, DIMENSION(0:lmax) :: H_LASER1,H_LASER2
REAL*8, DIMENSION(Ngrid_r,0:lmax) :: H_DIAG,H_LASER
REAL*8, DIMENSION(Ngrid_r,0:lmax,0:lmax) :: H_NONDIAG
COMPLEX*16, DIMENSION(Ngrid_r,0:lmax,Norb) :: RADIAL
COMPLEX*16, DIMENSION((lmax+1)*Ngrid_r) :: W,V

! Subroutine variables declaration

INTEGER*8 :: i,j,k,l,lprime,Column
!REAL*8 :: FACTORIAL
REAL*8, DIMENSION((lmax+1)*Ngrid_r,(lmax+1)*Ngrid_r) :: HAMILTONIAN
COMPLEX*16 :: ALPHA,BETA
COMPLEX*16, DIMENSION((lmax+1)*Ngrid_r,(lmax+1)*Ngrid_r) :: A
COMPLEX*16, DIMENSION((lmax+1)*Ngrid_r) :: Y

! Compute the upper part of the Hamiltonian. Then the whole matrix
! is filled by symmetry

HAMILTONIAN=0.D0

! Fill all diagonal blocks

DO l=0,lmax

   ! Fill the diagonals

   DO i=1,Ngrid_r

      j=i+l*Ngrid_r

      HAMILTONIAN(j,j)=H_DIAG(i,l)

      !PRINT*,'HAMILTONIAN',HAMILTONIAN(j,j)

   END DO

   !READ(*,*)

   ! Fill the first superdiagonal

   DO i=1+l*Ngrid_r,Ngrid_r-1+l*Ngrid_r

      j=i+1

      HAMILTONIAN(i,j)=H_DIAG1

   END DO

   ! Fill the second superdiagonal

   DO i=1+l*Ngrid_r,Ngrid_r-2+l*Ngrid_r

      j=i+2

      HAMILTONIAN(i,j)=H_DIAG2

   END DO

   ! Fill the non-diagonal blocks from laser and potential couplings

   DO lprime=l+1,lmax

      IF (lprime .EQ. l+1) THEN

         DO k=1,Ngrid_r

            i=k+l*Ngrid_r
            j=k+lprime*Ngrid_r

            HAMILTONIAN(i,j)=H_LASER(k,l)+H_NONDIAG(k,l,lprime)

         END DO

      ELSE

         DO k=1,Ngrid_r

            i=k+l*Ngrid_r
            j=k+lprime*Ngrid_r

            HAMILTONIAN(i,j)=H_NONDIAG(k,l,lprime)

         END DO

      END IF

   END DO

END DO

! Filling the lower part by symmetry

DO i=1,(lmax+1)*Ngrid_r

   DO j=i,(lmax+1)*Ngrid_r

      IF (j /= i) THEN

      HAMILTONIAN(j,i)=HAMILTONIAN(i,j)  

      END IF

   END DO

END DO

! Initializing quantities

A=-DCMPLX(0.D0,1.D0)*step_time*HAMILTONIAN

! Define V and W vectors suited for Taylor expansion

V=DCMPLX(0.D0,0.D0)
W=DCMPLX(0.D0,0.D0)

DO l=0,lmax

   DO i=1,Ngrid_r

      j=i+l*Ngrid_r

      V(j)=RADIAL(i,l,Orbital)
      W(j)=RADIAL(i,l,Orbital)
      !PRINT*,'V',V(j)
      
      !IF(j .EQ. 500 .OR. j .EQ. 1000 .OR. j .EQ. 1500) READ(*,*)

   END DO

   !PRINT*,'V',V(i)
   !PRINT*,'fin'
   !READ(*,*)

END DO

! m-order Taylor expansion with BLAS subroutine

!ALPHA=DCMPLX(1.D0,0.D0)
!BETA=DCMPLX(0.D0,0.D0)

!DO k=1,m_Taylor

!   Y=DCMPLX(0.D0,0.D0)

!   CALL ZGEMV('N',(lmax+1)*Ngrid_r,(lmax+1)*Ngrid_r,ALPHA,A,&
!        &(lmax+1)*Ngrid_r,V,1,BETA,Y,1)

!   DO i=1,(lmax+1)*Ngrid_r

      !PRINT*,'Y',Y(i)

!      W(i)=W(i)+Y(i)/FACTORIAL(k) ! Taylor summation
!      V(i)=Y(i)

      !PRINT*,'W',W(i)

!   END DO

   !READ(*,*)

!END DO

!PRINT*,'Here'

!RETURN

! m-order Taylor expansion with self-made subroutine

DO k=1,m_Taylor

   Y=DCMPLX(0.D0,0.D0)

   DO l=0,lmax

      IF (l == 0) THEN

         DO i=1,Ngrid_r

            IF (i == 1) THEN
               
               Y(i)=A(i,i)*V(i)+A(i,i+1)*V(i+1)+A(i,i+2)*V(i+2)

               DO Column=1,lmax

                  Y(i)=Y(i)+A(i,i+Column*Ngrid_r)*V(i+Column*Ngrid_r)

               END DO

               !PRINT*,'Y(i)',Y(i)

            ELSE IF (i == 2) THEN
               
               Y(i)=A(i,i-1)*V(i-1)+A(i,i)*V(i)+A(i,i+1)*V(i+1)+A(i,i+2)*V(i+2)

               DO Column=1,lmax

                  Y(i)=Y(i)+A(i,i+Column*Ngrid_r)*V(i+Column*Ngrid_r)

               END DO

               !PRINT*,'Y(i)',Y(i)

            ELSE IF (i /= Ngrid_r .AND. i /= Ngrid_r-1) THEN

               Y(i)=A(i,i-2)*V(i-2)+A(i,i-1)*V(i-1)+A(i,i)*V(i)+A(i,i+1)*V(i+1)&
                   &+A(i,i+2)*V(i+2)

               DO Column=1,lmax

                  Y(i)=Y(i)+A(i,i+Column*Ngrid_r)*V(i+Column*Ngrid_r)

               END DO

               !PRINT*,'Y(i)',Y(i)

               !PRINT*,'A(i,i-2),V(i-2)',A(i,i-2),V(i-2)
               !PRINT*,'A(i,i-1),V(i-1)',A(i,i-1),V(i-1)
               !PRINT*,'A(i,i),V(i)',A(i,i),V(i)
               !PRINT*,'A(i,i+1),V(i+1)',A(i,i+1),V(i+1)
               !PRINT*,'A(i,i+2),V(i+2)',A(i,i+2),V(i+2)
               !PRINT*,'A(i,i+Ngrid_r),V(i+Ngrid_r)',A(i,i+Ngrid_r),V(i+Ngrid_r)
               

               !PRINT*,'Y',Y(i)
               !READ(*,*)

            ELSE IF (i == Ngrid_r-1) THEN

               Y(i)=A(i,i-2)*V(i-2)+A(i,i-1)*V(i-1)+A(i,i)*V(i)+A(i,i+1)*V(i+1)

               DO Column=1,lmax

                  Y(i)=Y(i)+A(i,i+Column*Ngrid_r)*V(i+Column*Ngrid_r)

               END DO

               !PRINT*,'Y(i)',Y(i)

               !PRINT*,'A(i,i-2),V(i-2)',A(i,i-2),V(i-2)
               !PRINT*,'A(i,i-1),V(i-1)',A(i,i-1),V(i-1)
               !PRINT*,'A(i,i),V(i)',A(i,i),V(i)
               !PRINT*,'A(i,i+1),V(i+1)',A(i,i+1),V(i+1)
               !PRINT*,'A(i,i+Ngrid_r),V(i+Ngrid_r)',A(i,i+Ngrid_r),V(i+Ngrid_r)
               

               !PRINT*,'Y',Y(i)
               !READ(*,*)

            ELSE IF (i == Ngrid_r) THEN

               Y(i)=A(i,i-2)*V(i-2)+A(i,i-1)*V(i-1)+A(i,i)*V(i)

               DO Column=1,lmax

                  Y(i)=Y(i)+A(i,i+Column*Ngrid_r)*V(i+Column*Ngrid_r)

               END DO

               !PRINT*,'Y(i)',Y(i)

            END IF

         END DO

      ELSE IF (l /= lmax) THEN

         DO i=1,Ngrid_r
            
            j=i+l*Ngrid_r

            IF (i == 1) THEN

               Y(j)=A(j,j+1)*V(j+1)+A(j,j+2)*V(j+2)

               DO Column=-l,lmax-l

                  Y(j)=Y(j)+A(j,j+Column*Ngrid_r)*V(j+Column*Ngrid_r)

               END DO

               !PRINT*,'Y(j)',Y(j)

            ELSE IF (i == 2) THEN

               Y(j)=A(j,j-1)*V(j-1)+A(j,j+1)*V(j+1)+A(j,j+2)*V(j+2)

               DO Column=-l,lmax-l

                  Y(j)=Y(j)+A(j,j+Column*Ngrid_r)*V(j+Column*Ngrid_r)

               END DO

               !PRINT*,'Y(j)',Y(j)

            ELSE IF (i /= Ngrid_r .AND. i /= Ngrid_r-1) THEN

               Y(j)=A(j,j-2)*V(j-2)+A(j,j-1)*V(j-1)+A(j,j+1)*V(j+1)+&
                    &A(j,j+2)*V(j+2)

               DO Column=-l,lmax-l

                  Y(j)=Y(j)+A(j,j+Column*Ngrid_r)*V(j+Column*Ngrid_r)

               END DO

               !PRINT*,'Y(j)',Y(j)

            ELSE IF (i == Ngrid_r-1) THEN

               Y(j)=A(j,j-2)*V(j-2)+A(j,j-1)*V(j-1)+A(j,j+1)*V(j+1)

               DO Column=-l,lmax-l

                  Y(j)=Y(j)+A(j,j+Column*Ngrid_r)*V(j+Column*Ngrid_r)

               END DO

               !PRINT*,'Y(j)',Y(j)

            ELSE IF (i == Ngrid_r) THEN

               Y(j)=A(j,j-2)*V(j-2)+A(j,j-1)*V(j-1)

               DO Column=-l,lmax-l

                  Y(j)=Y(j)+A(j,j+Column*Ngrid_r)*V(j+Column*Ngrid_r)

               END DO

               !PRINT*,'Y(j)',Y(j)

            END IF

         END DO

      ELSE IF (l == lmax) THEN

         DO i=1,Ngrid_r

            j=i+l*Ngrid_r

            IF (i == 1) THEN

               Y(j)=A(j,j)*V(j)+A(j,j+1)*V(j+1)+A(j,j+2)*V(j+2)

               DO Column=-l,-1

                  Y(j)=Y(j)+A(j,j+Column*Ngrid_r)*V(j+Column*Ngrid_r)

               END DO
       
               !PRINT*,'Y(j)',Y(j)

            ELSE IF (i == 2) THEN

               Y(j)=A(j,j-1)*V(j-1)+A(j,j)*V(j)+A(j,j+1)*V(j+1)+A(j,j+2)*V(j+2)

               DO Column=-l,-1

                  Y(j)=Y(j)+A(j,j+Column*Ngrid_r)*V(j+Column*Ngrid_r)

               END DO

               !PRINT*,'Y(j)',Y(j)

            ELSE IF (i /= Ngrid_r .AND. i /= Ngrid_r-1) THEN

                Y(j)=A(j,j-2)*V(j-2)+A(j,j-1)*V(j-1)+A(j,j)*V(j)+&
                     &A(j,j+1)*V(j+1)+A(j,j+2)*V(j+2)

               DO Column=-l,-1

                  Y(j)=Y(j)+A(j,j+Column*Ngrid_r)*V(j+Column*Ngrid_r)

               END DO

               !PRINT*,'Y(j)',Y(j)

            ELSE IF (i == Ngrid_r-1) THEN

               Y(j)=A(j,j-2)*V(j-2)+A(j,j-1)*V(j-1)+A(j,j)*V(j)+A(j,j+1)*V(j+1)

               DO Column=-l,-1

                  Y(j)=Y(j)+A(j,j+Column*Ngrid_r)*V(j+Column*Ngrid_r)

               END DO

               !PRINT*,'Y(j)',Y(j)

            ELSE IF (i == Ngrid_r) THEN

               Y(j)=A(j,j-2)*V(j-2)+A(j,j-1)*V(j-1)+A(j,j)*V(j)

               DO Column=-l,-1

                  Y(j)=Y(j)+A(j,j+Column*Ngrid_r)*V(j+Column*Ngrid_r)

               END DO

               !PRINT*,'Y(j)',Y(j)

            END IF
           
         END DO

      END IF
   
   END DO

   DO i=1,(lmax+1)*Ngrid_r

      !PRINT*,'Y',Y(i)

      W(i)=W(i)+Y(i)/FACTORIAL(k) ! Taylor summation
      V(i)=Y(i)

      !PRINT*,'W',W(i)

   END DO

   !READ(*,*)

END DO

RETURN

END SUBROUTINE TAYLOR

!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!%                                                                            $
!%                                                                            $
!%                           SUBROUTINE EXPOMAT                               $
!%                                                                            $
!%                                                                            $
!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!==============================================================================
!
! This subroutine sets all quantities required for Expokit subroutine ZGEXPV
!
!==============================================================================

SUBROUTINE EXPOMAT(Ngrid_r,m_Krylov,step_time,tolerance,RADIAL,Norb,&
     &Orbital,lmin,lmax,L_Coupling,H_DIAG,H_DIAG1,H_DIAG2,&
     &H_NONDIAG,H_LASER,H_LASER1,H_LASER2,W)

IMPLICIT NONE

! Input/Output variables

INTEGER*8 :: Ngrid_r,m_Krylov,Norb,Orbital,lmin,lmax,L_Coupling
REAL*8 :: step_time,tolerance,H_DIAG1,H_DIAG2
REAL*8, DIMENSION(0:lmax) :: H_LASER1,H_LASER2
REAL*8, DIMENSION(Ngrid_r,0:lmax) :: H_DIAG,H_LASER
REAL*8, DIMENSION(Ngrid_r,0:lmax,0:lmax) :: H_NONDIAG
COMPLEX*16, DIMENSION((lmax+1)*Ngrid_r) :: W
COMPLEX*16, DIMENSION(Ngrid_r,0:lmax,Norb) :: RADIAL

! Local variables

INTEGER*8 :: i,j,l,lwsp,liwsp,iflag,itrace
REAL*8 :: Fnorm
INTEGER*8, DIMENSION(:), ALLOCATABLE :: iwsp
!REAL*8, DIMENSION((lmax+1)*Ngrid_r,(lmax+1)*Ngrid_r) :: HAMILTONIAN
COMPLEX*16, DIMENSION((lmax+1)*Ngrid_r) :: V
!COMPLEX*16, DIMENSION((lmax+1)*Ngrid_r,(lmax+1)*Ngrid_r) :: A
COMPLEX*16, DIMENSION(:), ALLOCATABLE :: wsp

! Compute the Frobenius norm of the Hamiltonian

Fnorm=1.D0

! Define V vector

DO l=0,lmax

   DO i=1,Ngrid_r

      j=i+l*Ngrid_r

      V(j)=RADIAL(i,l,Orbital)
      !PRINT*,'V',V(j)
      
      !IF(j .EQ. 500 .OR. j .EQ. 1000 .OR. j .EQ. 1500) READ(*,*)

   END DO

   !PRINT*,'V',V(i)
   !PRINT*,'fin'
   !READ(*,*)

END DO

!PRINT*,'V',V

! Compute some preliminary quantities involved in calculations

lwsp=(lmax+1)*Ngrid_r*(m_Krylov+1)+(lmax+1)*Ngrid_r+(m_Krylov+2)**2+&
&4*(m_Krylov+2)**2+10+1
liwsp=m_Krylov+2

!PRINT*,'Here100'

! Allocate arrays

ALLOCATE(wsp(lwsp),iwsp(liwsp))

! Call for the ZGEXPV subroutine that computes the action of the
! Hamiltonian matrix on the radial part of the wavefunction

!PRINT*,'m_Krylov,step_time,tolerance',m_Krylov,step_time,tolerance

itrace=0 ! Silent mode: 0 or stpe-by-step info: 1

CALL ZGEXPV((lmax+1)*Ngrid_r,Ngrid_r,lmin,lmax,L_Coupling,m_Krylov,step_time,&
     &V,W,tolerance,Fnorm,wsp,lwsp,iwsp,liwsp,matvec,itrace,iflag,H_DIAG,&
     &H_DIAG1,H_DIAG2,H_NONDIAG,H_LASER,H_LASER1,H_LASER2)

!PRINT*,'iflag,wsp(9)',iflag,wsp(9)

!READ(*,*)

RETURN

END SUBROUTINE EXPOMAT

!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!%                                                                            $
!%                                                                            $
!%                            SUBROUTINE MATVEC                               $
!%                                                                            $
!%                                                                            $
!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!==============================================================================
!
! This subroutine calculates the multiplication of a matrix A by a vector X. 
! It is suited for the use of Expokit ZGEXPV subroutine. The calculation is the
! following
!
!                                    A*X=Y                  
! 
! and it takes into account the sparsity of the matrix A (the same scheme is
! used in the Taylor expansion) 
!
!==============================================================================

SUBROUTINE MATVEC(X,Y,n,Ngrid_r,lmin,lmax,L_Coupling,H_DIAG,H_DIAG1,H_DIAG2,&
     &H_NONDIAG,H_LASER,H_LASER1,H_LASER2)

IMPLICIT NONE

! Input/Output variables declaration

INTEGER*8 :: n,Ngrid_r,lmin,lmax,L_Coupling
REAL*8 :: H_DIAG1,H_DIAG2
REAL*8, DIMENSION(Ngrid_r,0:lmax) :: H_DIAG,H_LASER
REAL*8, DIMENSION(0:lmax) :: H_LASER1,H_LASER2
REAL*8, DIMENSION(Ngrid_r,0:lmax,0:lmax) :: H_NONDIAG
COMPLEX*16, DIMENSION(n) :: X,Y

double complex :: wf(1:Ngrid_r,0:lmax)
double complex,dimension(1:size(wf,1),0:size(wf,2)-1) :: wf1,wf2

! Subroutine variables declaration

INTEGER :: i,j,l

Y=DCMPLX(0.D0,0.D0)

wf=DCMPLX(0.D0,0.D0)
wf1=DCMPLX(0.D0,0.D0)
wf2=DCMPLX(0.D0,0.D0)

DO i=1,Ngrid_r

   DO l=0,lmax

      j=i+l*Ngrid_r

      wf(i,l)=X(j)

   END DO

END DO

wf1=wf

call velmatmul(Ngrid_r,wf1,wf2,H_DIAG,H_DIAG1,H_DIAG2,H_NONDIAG,&
       &H_LASER,H_LASER1,H_LASER2,lmin,lmax,L_Coupling)

wf=wf2

DO i=1,Ngrid_r

   DO l=0,lmax

      j=i+l*Ngrid_r

      Y(j)=-DCMPLX(0.D0,1.D0)*wf(i,l)

   END DO

END DO

RETURN

END SUBROUTINE MATVEC

!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!%                                                                            $
!%                                                                            $
!%                          SUBROUTINE CONVERGENCE                            $
!%                                                                            $
!%                                                                            $
!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!==============================================================================
!
! This subroutine calculates the convergence parameter defined as
!
!            Conv=int_0^infty dr r^2[n_i(r)-n_i-1(r)]^2
!
!==============================================================================

SUBROUTINE CONVERGENCE(r,Ngrid_r,step_r,EFFPOT_OLD,EFFPOT,lmax,mmax,Conv)

IMPLICIT NONE

! Input/Output varables declaration

INTEGER*8 :: Ngrid_r,lmax,mmax
REAL*8 :: step_r,Conv
REAL*8, DIMENSION(Ngrid_r) :: r
REAL*8, DIMENSION(Ngrid_r,0:lmax,0:lmax,0:mmax) :: EFFPOT_OLD,EFFPOT

! Subroutine variables declaration

INTEGER :: i
REAL*8, DIMENSION(Ngrid_r) :: INTEGD

DO i=1,Ngrid_r

   INTEGD(i)=(r(i)*(EFFPOT(i,0,0,0)-EFFPOT_OLD(i,0,0,0)))**2
   !INTEGD(i)=r(i)

END DO

!PRINT*,'N_NEW(1),N_OLD(1)',N_NEW(1),N_OLD(1)
!PRINT*,'INTEGD',INTEGD

!READ(*,*)

CALL COMPOSITE_SIMP(Ngrid_r,step_r,INTEGD,Conv)

RETURN

END SUBROUTINE CONVERGENCE

!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!%                                                                            $
!%                                                                            $
!%                        SUBROUTINE COMPOSITE_SIMP                           $
!%                                                                            $
!%                                                                            $
!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!==============================================================================
!
! This subroutine calculates numerical integration using the following
! composite Simpson's rule:
!
! int_a^b f(x) dx = (h/3)[f(x1)+4*Sum_j=1^N/2 f(x_2j)
!                              +2*Sum_j=1^{N/2-1} f(x_{2j+1})+f(x_{N+1})]
!
! where the range [a,b] is shared into N equal intervals (N must be even), with
! N+1 points and f(x1)=a, f(x_{N+1})=b and h=(b-a)/N.
!
!==============================================================================

SUBROUTINE COMPOSITE_SIMP(Ngrid,step,INTEGD,Integrale)

IMPLICIT NONE

! Input/Output declaration

INTEGER*8 :: Ngrid
REAL*8 :: step,Integrale
REAL*8, DIMENSION(Ngrid) :: INTEGD

! Subroutine variable declaration

INTEGER*8 :: i
REAL*8 :: Sumb,Sume,Sumo

REAL*8, PARAMETER :: PI=3.14159265358979D0

!PRINT*,'Ngrid',Ngrid

Sumb=INTEGD(1)+INTEGD(Ngrid) ! Boundaries summation

!   PRINT*,'INTEGD(1)',INTEGD(1)
!   PRINT*,'Sumb',Sumb

Sume=0.D0

   DO i=1,(Ngrid-1)/2

      Sume=Sume+INTEGD(2*i) ! Even summation
      !PRINT*,'Sume',Sume

   END DO

   Sume=4.D0*Sume

!   PRINT*,'Sume',Sume

   Sumo=0.D0

   DO i=1,(Ngrid-1)/2-1

      Sumo=Sumo+INTEGD(2*i+1) ! Odd summation

   END DO

   Sumo=2.D0*Sumo

!   PRINT*,'Sumo',Sumo

   Integrale=step*(Sumb+Sume+Sumo)/3.D0 ! Simpson's rule

!   PRINT*,'Integrale',Integrale
!   READ(*,*)

RETURN

END SUBROUTINE COMPOSITE_SIMP

!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!%                                                                            $
!%                                                                            $
!%                       SUBROUTINE INTEG_GAUSS_LEG                           $
!%                                                                            $
!%                                                                            $
!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!==============================================================================
!
! This subroutine calculates numerical integration using the following
! Gauss-Legendre formula
!
! int_a^b dx f(x) = (b-a)/2 * Sum_i=1^n W(Xi)*f((b-a)/2 * Xi + (b+a)/2)
!
!==============================================================================

SUBROUTINE INTEG_GAUSS_LEG(N_GL,ib,sb,COEF,GL_ARG,GL_WEIGHTS)

IMPLICIT NONE

! Input/Output declaration

INTEGER*8 :: N_GL
REAL*8 :: ib,sb,COEF
REAL*8, DIMENSION(N_GL) :: GL_ARG,GL_WEIGHTS

! Subroutine variables declaration

INTEGER*8 :: i,j
REAL*8, DIMENSION(N_GL) :: ROOTS,WEIGHTS

! Loading roots and weights for N_GL order Gauss-Legendre integration

CALL GAUSS_LEG_ROOTS_WEIGHTS(N_GL,ROOTS,WEIGHTS)

DO i=1,N_GL/2
   
   j=i+N_GL/2

   GL_ARG(j)=(sb-ib)*ROOTS(i)/2.D0+(sb+ib)/2.D0
   GL_WEIGHTS(j)=WEIGHTS(i)

END DO

DO i=N_GL/2+1,N_GL
   
   j=i-N_GL/2

   GL_ARG(j)=(sb-ib)*ROOTS(i)/2.D0+(sb+ib)/2.D0
   GL_WEIGHTS(j)=WEIGHTS(i)

END DO

COEF=(sb-ib)/2.D0

RETURN

END SUBROUTINE INTEG_GAUSS_LEG

!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!%                                                                            $
!%                                                                            $
!%                   SUBROUTINE GAUSS_LEG_ROOTS_WEIGHTS                       $
!%                                                                            $
!%                                                                            $
!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!==============================================================================
!
! This subroutine sends back the Gauss-Legendre roots and weights.
!
!==============================================================================

SUBROUTINE GAUSS_LEG_ROOTS_WEIGHTS(N_GL,GL_ROOTS,GL_WEIGHTS)

IMPLICIT NONE

! Input/Output declaration

INTEGER*8 :: N_GL
REAL*8, DIMENSION(N_GL) :: GL_ROOTS,GL_WEIGHTS

IF (N_GL .EQ. 2) THEN

GL_WEIGHTS(1)=1.00000000000000000000000000000000d0
GL_WEIGHTS(2)=GL_WEIGHTS(1)

GL_ROOTS(1)=0.57735026918962576450914878050196d0
GL_ROOTS(2)=-GL_ROOTS(1)

END IF

IF (N_GL .EQ. 4) THEN

GL_WEIGHTS(1)=0.65214515486254614262693605077800d0
GL_WEIGHTS(2)=0.34785484513745385737306394922200d0
GL_WEIGHTS(3)=GL_WEIGHTS(2)
GL_WEIGHTS(4)=GL_WEIGHTS(1)

GL_ROOTS(1)=0.33998104358485626480266575910324d0
GL_ROOTS(2)=0.86113631159405257522394648889281d0
GL_ROOTS(3)=-GL_ROOTS(2)
GL_ROOTS(4)=-GL_ROOTS(1)

END IF

IF (N_GL .EQ. 6) THEN

GL_WEIGHTS(1)=0.46791393457269104738987034398955d0
GL_WEIGHTS(2)=0.36076157304813860756983351383772d0
GL_WEIGHTS(3)=0.17132449237917034504029614217273d0
GL_WEIGHTS(4)=GL_WEIGHTS(3)
GL_WEIGHTS(5)=GL_WEIGHTS(2)
GL_WEIGHTS(6)=GL_WEIGHTS(1)

GL_ROOTS(1)=0.23861918608319690863050172168071d0
GL_ROOTS(2)=0.66120938646626451366139959501991d0
GL_ROOTS(3)=0.93246951420315202781230155449399d0
GL_ROOTS(4)=-GL_ROOTS(3)
GL_ROOTS(5)=-GL_ROOTS(2)
GL_ROOTS(6)=-GL_ROOTS(1)

END IF

IF (N_GL .EQ. 8) THEN

GL_WEIGHTS(1)=0.36268378337836198296515044927720d0
GL_WEIGHTS(2)=0.31370664587788728733796220198660d0
GL_WEIGHTS(3)=0.22238103445337447054435599442624d0
GL_WEIGHTS(4)=0.10122853629037625915253135430996d0
GL_WEIGHTS(5)=GL_WEIGHTS(4)
GL_WEIGHTS(6)=GL_WEIGHTS(3)
GL_WEIGHTS(7)=GL_WEIGHTS(2)
GL_WEIGHTS(8)=GL_WEIGHTS(1)

GL_ROOTS(1)=0.18343464249564980493947614236018d0
GL_ROOTS(2)=0.52553240991632898581773904918925d0
GL_ROOTS(3)=0.79666647741362673959155393647583d0
GL_ROOTS(4)=0.96028985649753623168356086856947d0
GL_ROOTS(5)=-GL_ROOTS(4)
GL_ROOTS(6)=-GL_ROOTS(3)
GL_ROOTS(7)=-GL_ROOTS(2)
GL_ROOTS(8)=-GL_ROOTS(1)

END IF

IF (N_GL .EQ. 10) THEN

GL_WEIGHTS(1)=0.29552422471475287017389299465134d0
GL_WEIGHTS(2)=0.26926671930999635509122692156947d0
GL_WEIGHTS(3)=0.21908636251598204399553493422816d0
GL_WEIGHTS(4)=0.14945134915058059314577633965770d0
GL_WEIGHTS(5)=0.06667134430868813759356880989333d0
GL_WEIGHTS(6)=GL_WEIGHTS(5)
GL_WEIGHTS(7)=GL_WEIGHTS(4)
GL_WEIGHTS(8)=GL_WEIGHTS(3)
GL_WEIGHTS(9)=GL_WEIGHTS(2)
GL_WEIGHTS(10)=GL_WEIGHTS(1)

GL_ROOTS(1)=0.14887433898163121088482600112972d0
GL_ROOTS(2)=0.43339539412924719079926594316578d0
GL_ROOTS(3)=0.67940956829902440623432736511487d0
GL_ROOTS(4)=0.86506336668898451073209668842349d0
GL_ROOTS(5)=0.97390652851717172007796401208445d0
GL_ROOTS(6)=-GL_ROOTS(5)
GL_ROOTS(7)=-GL_ROOTS(4)
GL_ROOTS(8)=-GL_ROOTS(3)
GL_ROOTS(9)=-GL_ROOTS(2)
GL_ROOTS(10)=-GL_ROOTS(1)

END IF

IF (N_GL .EQ. 12) THEN

GL_WEIGHTS(1)=0.249147045813402785001d0
GL_WEIGHTS(2)=0.233492536538354808761d0
GL_WEIGHTS(3)=0.203167426723065921749d0
GL_WEIGHTS(4)=0.160078328543346226335d0
GL_WEIGHTS(5)=0.10693932599531843096d0
GL_WEIGHTS(6)=0.047175336386511827195d0
GL_WEIGHTS(7)=GL_WEIGHTS(6)
GL_WEIGHTS(8)=GL_WEIGHTS(5)
GL_WEIGHTS(9)=GL_WEIGHTS(4)
GL_WEIGHTS(10)=GL_WEIGHTS(3)
GL_WEIGHTS(11)=GL_WEIGHTS(2)
GL_WEIGHTS(12)=GL_WEIGHTS(1)

GL_ROOTS(1)=0.125233408511468915472d0
GL_ROOTS(2)=0.36783149899818019375d0
GL_ROOTS(3)=0.5873179542866174473d0
GL_ROOTS(4)=0.76990267419430468704d0
GL_ROOTS(5)=0.90411725637047485668d0
GL_ROOTS(6)=0.98156063424671925069d0
GL_ROOTS(7)=-GL_ROOTS(6)
GL_ROOTS(8)=-GL_ROOTS(5)
GL_ROOTS(9)=-GL_ROOTS(4)
GL_ROOTS(10)=-GL_ROOTS(3)
GL_ROOTS(11)=-GL_ROOTS(2)
GL_ROOTS(12)=-GL_ROOTS(1)

END IF

IF (N_GL .EQ. 16) THEN

GL_WEIGHTS(1)=0.18945061045506849629d0
GL_WEIGHTS(2)=0.18260341504492358887d0
GL_WEIGHTS(3)=0.16915651939500253819d0
GL_WEIGHTS(4)=0.14959598881657673208d0
GL_WEIGHTS(5)=0.12462897125553387205d0
GL_WEIGHTS(6)=0.09515851168249278481d0
GL_WEIGHTS(7)=0.06225352393864789286d0
GL_WEIGHTS(8)=0.027152459411754094852d0
GL_WEIGHTS(9)=GL_WEIGHTS(8)
GL_WEIGHTS(10)=GL_WEIGHTS(7)
GL_WEIGHTS(11)=GL_WEIGHTS(6)
GL_WEIGHTS(12)=GL_WEIGHTS(5)
GL_WEIGHTS(13)=GL_WEIGHTS(4)
GL_WEIGHTS(14)=GL_WEIGHTS(3)
GL_WEIGHTS(15)=GL_WEIGHTS(2)
GL_WEIGHTS(16)=GL_WEIGHTS(1)

GL_ROOTS(1)=0.09501250983763744019d0
GL_ROOTS(2)=0.28160355077925891323d0
GL_ROOTS(3)=0.4580167776572273863d0
GL_ROOTS(4)=0.6178762444026437484d0
GL_ROOTS(5)=0.7554044083550030339d0
GL_ROOTS(6)=0.8656312023878317439d0
GL_ROOTS(7)=0.9445750230732325761d0
GL_ROOTS(8)=0.9894009349916499326d0
GL_ROOTS(9)=-GL_ROOTS(8)
GL_ROOTS(10)=-GL_ROOTS(7)
GL_ROOTS(11)=-GL_ROOTS(6)
GL_ROOTS(12)=-GL_ROOTS(5)
GL_ROOTS(13)=-GL_ROOTS(4)
GL_ROOTS(14)=-GL_ROOTS(3)
GL_ROOTS(15)=-GL_ROOTS(2)
GL_ROOTS(16)=-GL_ROOTS(1)

END IF

IF (N_GL .EQ. 32) THEN

GL_WEIGHTS(1)=0.09654008851472780057d0
GL_WEIGHTS(2)=0.09563872007927485942d0
GL_WEIGHTS(3)=0.09384439908080456564d0
GL_WEIGHTS(4)=0.09117387869576388471d0
GL_WEIGHTS(5)=0.08765209300440381114d0
GL_WEIGHTS(6)=0.08331192422694675522d0
GL_WEIGHTS(7)=0.07819389578707030647d0
GL_WEIGHTS(8)=0.07234579410884850623d0
GL_WEIGHTS(9)=0.06582222277636184684d0
GL_WEIGHTS(10)=0.05868409347853554715d0
GL_WEIGHTS(11)=0.05099805926237617620d0
GL_WEIGHTS(12)=0.04283589802222668066d0
GL_WEIGHTS(13)=0.03427386291302143310d0
GL_WEIGHTS(14)=0.025392065309262059456d0
GL_WEIGHTS(15)=0.016274394730905670605d0
GL_WEIGHTS(16)=0.007018610009470096600d0
GL_WEIGHTS(17)=GL_WEIGHTS(16)
GL_WEIGHTS(18)=GL_WEIGHTS(15)
GL_WEIGHTS(19)=GL_WEIGHTS(14)
GL_WEIGHTS(20)=GL_WEIGHTS(13)
GL_WEIGHTS(21)=GL_WEIGHTS(12)
GL_WEIGHTS(22)=GL_WEIGHTS(11)
GL_WEIGHTS(23)=GL_WEIGHTS(10)
GL_WEIGHTS(24)=GL_WEIGHTS(9)
GL_WEIGHTS(25)=GL_WEIGHTS(8)
GL_WEIGHTS(26)=GL_WEIGHTS(7)
GL_WEIGHTS(27)=GL_WEIGHTS(6)
GL_WEIGHTS(28)=GL_WEIGHTS(5)
GL_WEIGHTS(29)=GL_WEIGHTS(4)
GL_WEIGHTS(30)=GL_WEIGHTS(3)
GL_WEIGHTS(31)=GL_WEIGHTS(2)
GL_WEIGHTS(32)=GL_WEIGHTS(1)

GL_ROOTS(1)=0.04830766568773831623d0
GL_ROOTS(2)=0.14447196158279649349d0
GL_ROOTS(3)=0.23928736225213707454d0
GL_ROOTS(4)=0.3318686022821276498d0
GL_ROOTS(5)=0.4213512761306353454d0
GL_ROOTS(6)=0.5068999089322293900d0
GL_ROOTS(7)=0.5877157572407623290d0
GL_ROOTS(8)=0.6630442669302152010d0
GL_ROOTS(9)=0.7321821187402896804d0
GL_ROOTS(10)=0.7944837959679424070d0
GL_ROOTS(11)=0.8493676137325699701d0
GL_ROOTS(12)=0.8963211557660521240d0
GL_ROOTS(13)=0.9349060759377396892d0
GL_ROOTS(14)=0.9647622555875064308d0
GL_ROOTS(15)=0.9856115115452683354d0
GL_ROOTS(16)=0.9972638618494815635d0
GL_ROOTS(17)=-GL_ROOTS(16)
GL_ROOTS(18)=-GL_ROOTS(15)
GL_ROOTS(19)=-GL_ROOTS(14)
GL_ROOTS(20)=-GL_ROOTS(13)
GL_ROOTS(21)=-GL_ROOTS(12)
GL_ROOTS(22)=-GL_ROOTS(11)
GL_ROOTS(23)=-GL_ROOTS(10)
GL_ROOTS(24)=-GL_ROOTS(9)
GL_ROOTS(25)=-GL_ROOTS(8)
GL_ROOTS(26)=-GL_ROOTS(7)
GL_ROOTS(27)=-GL_ROOTS(6)
GL_ROOTS(28)=-GL_ROOTS(5)
GL_ROOTS(29)=-GL_ROOTS(4)
GL_ROOTS(30)=-GL_ROOTS(3)
GL_ROOTS(31)=-GL_ROOTS(2)
GL_ROOTS(32)=-GL_ROOTS(1)

END IF

IF (N_GL .EQ. 64) THEN

GL_WEIGHTS(1)=0.04869095700913972038d0
GL_WEIGHTS(2)=0.04857546744150342693d0
GL_WEIGHTS(3)=0.04834476223480295717d0
GL_WEIGHTS(4)=0.04799938859645830773d0
GL_WEIGHTS(5)=0.04754016571483030866d0
GL_WEIGHTS(6)=0.04696818281621001733d0
GL_WEIGHTS(7)=0.04628479658131441730d0
GL_WEIGHTS(8)=0.04549162792741814448d0
GL_WEIGHTS(9)=0.04459055816375656306d0
GL_WEIGHTS(10)=0.04358372452932345338d0
GL_WEIGHTS(11)=0.04247351512365358901d0
GL_WEIGHTS(12)=0.04126256324262352861d0
GL_WEIGHTS(13)=0.03995374113272034139d0
GL_WEIGHTS(14)=0.03855015317861562913d0
GL_WEIGHTS(15)=0.03705512854024004604d0
GL_WEIGHTS(16)=0.03547221325688238381d0
GL_WEIGHTS(17)=0.03380516183714160939d0
GL_WEIGHTS(18)=0.03205792835485155359d0
GL_WEIGHTS(19)=0.030234657072402478868d0
GL_WEIGHTS(20)=0.028339672614259483228d0
GL_WEIGHTS(21)=0.026377469715054658672d0
GL_WEIGHTS(22)=0.024352702568710873338d0
GL_WEIGHTS(23)=0.022270173808383254159d0
GL_WEIGHTS(24)=0.020134823153530209372d0
GL_WEIGHTS(25)=0.017951715775697343085d0
GL_WEIGHTS(26)=0.015726030476024719322d0
GL_WEIGHTS(27)=0.013463047896718642598d0
GL_WEIGHTS(28)=0.011168139460131128819d0
GL_WEIGHTS(29)=0.008846759826363947723d0
GL_WEIGHTS(30)=0.006504457968978362856d0
GL_WEIGHTS(31)=0.004147033260562467635d0
GL_WEIGHTS(32)=0.0017832807216964329473d0
GL_WEIGHTS(33)=GL_WEIGHTS(32)
GL_WEIGHTS(34)=GL_WEIGHTS(31)
GL_WEIGHTS(35)=GL_WEIGHTS(30)
GL_WEIGHTS(36)=GL_WEIGHTS(29)
GL_WEIGHTS(37)=GL_WEIGHTS(28)
GL_WEIGHTS(38)=GL_WEIGHTS(27)
GL_WEIGHTS(39)=GL_WEIGHTS(26)
GL_WEIGHTS(40)=GL_WEIGHTS(25)
GL_WEIGHTS(41)=GL_WEIGHTS(24)
GL_WEIGHTS(42)=GL_WEIGHTS(23)
GL_WEIGHTS(43)=GL_WEIGHTS(22)
GL_WEIGHTS(44)=GL_WEIGHTS(21)
GL_WEIGHTS(45)=GL_WEIGHTS(20)
GL_WEIGHTS(46)=GL_WEIGHTS(19)
GL_WEIGHTS(47)=GL_WEIGHTS(18)
GL_WEIGHTS(48)=GL_WEIGHTS(17)
GL_WEIGHTS(49)=GL_WEIGHTS(16)
GL_WEIGHTS(50)=GL_WEIGHTS(15)
GL_WEIGHTS(51)=GL_WEIGHTS(14)
GL_WEIGHTS(52)=GL_WEIGHTS(13)
GL_WEIGHTS(53)=GL_WEIGHTS(12)
GL_WEIGHTS(54)=GL_WEIGHTS(11)
GL_WEIGHTS(55)=GL_WEIGHTS(10)
GL_WEIGHTS(56)=GL_WEIGHTS(9)
GL_WEIGHTS(57)=GL_WEIGHTS(8)
GL_WEIGHTS(58)=GL_WEIGHTS(7)
GL_WEIGHTS(59)=GL_WEIGHTS(6)
GL_WEIGHTS(60)=GL_WEIGHTS(5)
GL_WEIGHTS(61)=GL_WEIGHTS(4)
GL_WEIGHTS(62)=GL_WEIGHTS(3)
GL_WEIGHTS(63)=GL_WEIGHTS(2)
GL_WEIGHTS(64)=GL_WEIGHTS(1)

GL_ROOTS(1)=0.024350292663424432509d0
GL_ROOTS(2)=0.07299312178779903945d0
GL_ROOTS(3)=0.12146281929612055447d0
GL_ROOTS(4)=0.16964442042399281804d0
GL_ROOTS(5)=0.21742364374000708415d0
GL_ROOTS(6)=0.26468716220876741637d0
GL_ROOTS(7)=0.31132287199021095616d0
GL_ROOTS(8)=0.3572201583376681160d0
GL_ROOTS(9)=0.4022701579639916037d0
GL_ROOTS(10)=0.4463660172534640880d0
GL_ROOTS(11)=0.4894031457070529575d0
GL_ROOTS(12)=0.5312794640198945457d0
GL_ROOTS(13)=0.5718956462026340343d0
GL_ROOTS(14)=0.6111553551723932502d0
GL_ROOTS(15)=0.6489654712546573399d0
GL_ROOTS(16)=0.6852363130542332426d0
GL_ROOTS(17)=0.7198818501716108268d0
GL_ROOTS(18)=0.7528199072605318966d0
GL_ROOTS(19)=0.7839723589433414076d0
GL_ROOTS(20)=0.8132653151227975597d0
GL_ROOTS(21)=0.8406292962525803628d0
GL_ROOTS(22)=0.8659993981540928198d0
GL_ROOTS(23)=0.8893154459951141059d0
GL_ROOTS(24)=0.9105221370785028058d0
GL_ROOTS(25)=0.9295691721319395758d0
GL_ROOTS(26)=0.9464113748584028161d0
GL_ROOTS(27)=0.9610087996520537189d0
GL_ROOTS(28)=0.9733268277899109637d0
GL_ROOTS(29)=0.9833362538846259569d0
GL_ROOTS(30)=0.9910133714767443207d0
GL_ROOTS(31)=0.9963401167719552793d0
GL_ROOTS(32)=0.9993050417357721395d0
GL_ROOTS(33)=-GL_ROOTS(32)
GL_ROOTS(34)=-GL_ROOTS(31)
GL_ROOTS(35)=-GL_ROOTS(30)
GL_ROOTS(36)=-GL_ROOTS(29)
GL_ROOTS(37)=-GL_ROOTS(28)
GL_ROOTS(38)=-GL_ROOTS(27)
GL_ROOTS(39)=-GL_ROOTS(26)
GL_ROOTS(40)=-GL_ROOTS(25)
GL_ROOTS(41)=-GL_ROOTS(24)
GL_ROOTS(42)=-GL_ROOTS(23)
GL_ROOTS(43)=-GL_ROOTS(22)
GL_ROOTS(44)=-GL_ROOTS(21)
GL_ROOTS(45)=-GL_ROOTS(20)
GL_ROOTS(46)=-GL_ROOTS(19)
GL_ROOTS(47)=-GL_ROOTS(18)
GL_ROOTS(48)=-GL_ROOTS(17)
GL_ROOTS(49)=-GL_ROOTS(16)
GL_ROOTS(50)=-GL_ROOTS(15)
GL_ROOTS(51)=-GL_ROOTS(14)
GL_ROOTS(52)=-GL_ROOTS(13)
GL_ROOTS(53)=-GL_ROOTS(12)
GL_ROOTS(54)=-GL_ROOTS(11)
GL_ROOTS(55)=-GL_ROOTS(10)
GL_ROOTS(56)=-GL_ROOTS(9)
GL_ROOTS(57)=-GL_ROOTS(8)
GL_ROOTS(58)=-GL_ROOTS(7)
GL_ROOTS(59)=-GL_ROOTS(6)
GL_ROOTS(60)=-GL_ROOTS(5)
GL_ROOTS(61)=-GL_ROOTS(4)
GL_ROOTS(62)=-GL_ROOTS(3)
GL_ROOTS(63)=-GL_ROOTS(2)
GL_ROOTS(64)=-GL_ROOTS(1)

END IF

IF (N_GL .EQ. 128) THEN

GL_WEIGHTS(1)=0.024446180196262518211d0
GL_WEIGHTS(2)=0.024431569097850045055d0
GL_WEIGHTS(3)=0.024402355633849582093d0
GL_WEIGHTS(4)=0.024358557264690625853d0
GL_WEIGHTS(5)=0.024300200167971865323d0
GL_WEIGHTS(6)=0.024227319222815248120d0
GL_WEIGHTS(7)=0.024139957989019284998d0
GL_WEIGHTS(8)=0.024038168681024052638d0
GL_WEIGHTS(9)=0.023922012136703455672d0
GL_WEIGHTS(10)=0.023791557781003400639d0
GL_WEIGHTS(11)=0.023646883584447615144d0
GL_WEIGHTS(12)=0.023488076016535913153d0
GL_WEIGHTS(13)=0.023315229994062760122d0
GL_WEIGHTS(14)=0.023128448824387027879d0
GL_WEIGHTS(15)=0.022927844143686846920d0
GL_WEIGHTS(16)=0.022713535850236461310d0
GL_WEIGHTS(17)=0.022485652032744966872d0
GL_WEIGHTS(18)=0.022244328893799765105d0
GL_WEIGHTS(19)=0.021989710668460491434d0
GL_WEIGHTS(20)=0.021721949538052075375d0
GL_WEIGHTS(21)=0.021441205539208460137d0
GL_WEIGHTS(22)=0.021147646468221348537d0
GL_WEIGHTS(23)=0.020841447780751149114d0
GL_WEIGHTS(24)=0.020522792486960069432d0
GL_WEIGHTS(25)=0.020191871042130041181d0
GL_WEIGHTS(26)=0.019848881232830862220d0
GL_WEIGHTS(27)=0.019494028058706602823d0
GL_WEIGHTS(28)=0.019127523609950945487d0
GL_WEIGHTS(29)=0.018749586940544708651d0
GL_WEIGHTS(30)=0.018360443937331343221d0
GL_WEIGHTS(31)=0.017960327185008685940d0
GL_WEIGHTS(32)=0.017549475827117704649d0
GL_WEIGHTS(33)=0.017128135423111376831d0
GL_WEIGHTS(34)=0.016696557801589204589d0
GL_WEIGHTS(35)=0.016255000909785187052d0
GL_WEIGHTS(36)=0.015803728659399346859d0
GL_WEIGHTS(37)=0.015343010768865144086d0
GL_WEIGHTS(38)=0.014873122602147314252d0
GL_WEIGHTS(39)=0.014394345004166846177d0
GL_WEIGHTS(40)=0.013906964132951985244d0
GL_WEIGHTS(41)=0.013411271288616332314d0
GL_WEIGHTS(42)=0.012907562739267347220d0
GL_WEIGHTS(43)=0.012396139543950922969d0
GL_WEIGHTS(44)=0.011877307372740279576d0
GL_WEIGHTS(45)=0.011351376324080416693d0
GL_WEIGHTS(46)=0.010818660739503076248d0
GL_WEIGHTS(47)=0.010279479015832157133d0
GL_WEIGHTS(48)=0.009734153415006805864d0
GL_WEIGHTS(49)=0.009183009871660874334d0
GL_WEIGHTS(50)=0.008626377798616749705d0
GL_WEIGHTS(51)=0.008064589890486057973d0
GL_WEIGHTS(52)=0.007497981925634728688d0
GL_WEIGHTS(53)=0.006926892566898813563d0
GL_WEIGHTS(54)=0.006351663161707188787d0
GL_WEIGHTS(55)=0.005772637542865698589d0
GL_WEIGHTS(56)=0.005190161832676330205d0
GL_WEIGHTS(57)=0.004604584256702955118d0
GL_WEIGHTS(58)=0.004016254983738642313d0
GL_WEIGHTS(59)=0.003425526040910215774d0
GL_WEIGHTS(60)=0.0028327514714579910953d0
GL_WEIGHTS(61)=0.0022382884309626187436d0
GL_WEIGHTS(62)=0.0016425030186690295388d0
GL_WEIGHTS(63)=0.0010458126793403487793d0
GL_WEIGHTS(64)=0.0004493809602920903764d0
GL_WEIGHTS(65)=GL_WEIGHTS(64)
GL_WEIGHTS(66)=GL_WEIGHTS(63)
GL_WEIGHTS(67)=GL_WEIGHTS(62)
GL_WEIGHTS(68)=GL_WEIGHTS(61)
GL_WEIGHTS(69)=GL_WEIGHTS(60)
GL_WEIGHTS(70)=GL_WEIGHTS(59)
GL_WEIGHTS(71)=GL_WEIGHTS(58)
GL_WEIGHTS(72)=GL_WEIGHTS(57)
GL_WEIGHTS(73)=GL_WEIGHTS(56)
GL_WEIGHTS(74)=GL_WEIGHTS(55)
GL_WEIGHTS(75)=GL_WEIGHTS(54)
GL_WEIGHTS(76)=GL_WEIGHTS(53)
GL_WEIGHTS(77)=GL_WEIGHTS(52)
GL_WEIGHTS(78)=GL_WEIGHTS(51)
GL_WEIGHTS(79)=GL_WEIGHTS(50)
GL_WEIGHTS(80)=GL_WEIGHTS(49)
GL_WEIGHTS(81)=GL_WEIGHTS(48)
GL_WEIGHTS(82)=GL_WEIGHTS(47)
GL_WEIGHTS(83)=GL_WEIGHTS(46)
GL_WEIGHTS(84)=GL_WEIGHTS(45)
GL_WEIGHTS(85)=GL_WEIGHTS(44)
GL_WEIGHTS(86)=GL_WEIGHTS(43)
GL_WEIGHTS(87)=GL_WEIGHTS(42)
GL_WEIGHTS(88)=GL_WEIGHTS(41)
GL_WEIGHTS(89)=GL_WEIGHTS(40)
GL_WEIGHTS(90)=GL_WEIGHTS(39)
GL_WEIGHTS(91)=GL_WEIGHTS(38)
GL_WEIGHTS(92)=GL_WEIGHTS(37)
GL_WEIGHTS(93)=GL_WEIGHTS(36)
GL_WEIGHTS(94)=GL_WEIGHTS(35)
GL_WEIGHTS(95)=GL_WEIGHTS(34)
GL_WEIGHTS(96)=GL_WEIGHTS(33)
GL_WEIGHTS(97)=GL_WEIGHTS(32)
GL_WEIGHTS(98)=GL_WEIGHTS(31)
GL_WEIGHTS(99)=GL_WEIGHTS(30)
GL_WEIGHTS(100)=GL_WEIGHTS(29)
GL_WEIGHTS(101)=GL_WEIGHTS(28)
GL_WEIGHTS(102)=GL_WEIGHTS(27)
GL_WEIGHTS(103)=GL_WEIGHTS(26)
GL_WEIGHTS(104)=GL_WEIGHTS(25)
GL_WEIGHTS(105)=GL_WEIGHTS(24)
GL_WEIGHTS(106)=GL_WEIGHTS(23)
GL_WEIGHTS(107)=GL_WEIGHTS(22)
GL_WEIGHTS(108)=GL_WEIGHTS(21)
GL_WEIGHTS(109)=GL_WEIGHTS(20)
GL_WEIGHTS(110)=GL_WEIGHTS(19)
GL_WEIGHTS(111)=GL_WEIGHTS(18)
GL_WEIGHTS(112)=GL_WEIGHTS(17)
GL_WEIGHTS(113)=GL_WEIGHTS(16)
GL_WEIGHTS(114)=GL_WEIGHTS(15)
GL_WEIGHTS(115)=GL_WEIGHTS(14)
GL_WEIGHTS(116)=GL_WEIGHTS(13)
GL_WEIGHTS(117)=GL_WEIGHTS(12)
GL_WEIGHTS(118)=GL_WEIGHTS(11)
GL_WEIGHTS(119)=GL_WEIGHTS(10)
GL_WEIGHTS(120)=GL_WEIGHTS(9)
GL_WEIGHTS(121)=GL_WEIGHTS(8)
GL_WEIGHTS(122)=GL_WEIGHTS(7)
GL_WEIGHTS(123)=GL_WEIGHTS(6)
GL_WEIGHTS(124)=GL_WEIGHTS(5)
GL_WEIGHTS(125)=GL_WEIGHTS(4)
GL_WEIGHTS(126)=GL_WEIGHTS(3)
GL_WEIGHTS(127)=GL_WEIGHTS(2)
GL_WEIGHTS(128)=GL_WEIGHTS(1)

GL_ROOTS(1)=0.012223698960615764198d0
GL_ROOTS(2)=0.03666379096873349333d0
GL_ROOTS(3)=0.06108196960413956810d0
GL_ROOTS(4)=0.08546364050451549864d0
GL_ROOTS(5)=0.10979423112764374667d0
GL_ROOTS(6)=0.13405919946118778512d0
GL_ROOTS(7)=0.15824404271422493400d0
GL_ROOTS(8)=0.18233430598533718241d0
GL_ROOTS(9)=0.20631559090207921715d0
GL_ROOTS(10)=0.23017356422665998641d0
GL_ROOTS(11)=0.25389396642269432086d0
GL_ROOTS(12)=0.27746262017790440281d0
GL_ROOTS(13)=0.30086543887767720267d0
GL_ROOTS(14)=0.3240884350244133752d0
GL_ROOTS(15)=0.3471177285976355084d0
GL_ROOTS(16)=0.3699395553498590266d0
GL_ROOTS(17)=0.3925402750332674427d0
GL_ROOTS(18)=0.4149063795522750155d0
GL_ROOTS(19)=0.4370245010371041629d0
GL_ROOTS(20)=0.4588814198335521954d0
GL_ROOTS(21)=0.4804640724041720259d0
GL_ROOTS(22)=0.5017595591361444643d0
GL_ROOTS(23)=0.5227551520511754785d0
GL_ROOTS(24)=0.5434383024128103634d0
GL_ROOTS(25)=0.5637966482266180839d0
GL_ROOTS(26)=0.5838180216287630896d0
GL_ROOTS(27)=0.6034904561585486242d0
GL_ROOTS(28)=0.6228021939105849108d0
GL_ROOTS(29)=0.6417416925623075572d0
GL_ROOTS(30)=0.6602976322726460521d0
GL_ROOTS(31)=0.6784589224477192594d0
GL_ROOTS(32)=0.6962147083695143324d0
GL_ROOTS(33)=0.7135543776835874133d0
GL_ROOTS(34)=0.7304675667419088065d0
GL_ROOTS(35)=0.7469441667970619812d0
GL_ROOTS(36)=0.7629743300440947228d0
GL_ROOTS(37)=0.7785484755064119669d0
GL_ROOTS(38)=0.7936572947621932902d0
GL_ROOTS(39)=0.8082917575079136601d0
GL_ROOTS(40)=0.8224431169556438425d0
GL_ROOTS(41)=0.8361029150609068471d0
GL_ROOTS(42)=0.8492629875779689692d0
GL_ROOTS(43)=0.8619154689395484606d0
GL_ROOTS(44)=0.8740527969580317987d0
GL_ROOTS(45)=0.8856677173453972174d0
GL_ROOTS(46)=0.8967532880491581844d0
GL_ROOTS(47)=0.9073028834017568139d0
GL_ROOTS(48)=0.9173101980809605370d0
GL_ROOTS(49)=0.9267692508789478433d0
GL_ROOTS(50)=0.9356743882779163758d0
GL_ROOTS(51)=0.9440202878302201821d0
GL_ROOTS(52)=0.9518019613412643862d0
GL_ROOTS(53)=0.9590147578536999281d0
GL_ROOTS(54)=0.9656543664319652686d0
GL_ROOTS(55)=0.9717168187471365809d0
GL_ROOTS(56)=0.9771984914639073872d0
GL_ROOTS(57)=0.9820961084357185360d0
GL_ROOTS(58)=0.9864067427245862089d0
GL_ROOTS(59)=0.9901278184917343833d0
GL_ROOTS(60)=0.9932571129002129353d0
GL_ROOTS(61)=0.9957927585349811869d0
GL_ROOTS(62)=0.9977332486255140199d0
GL_ROOTS(63)=0.9990774599773758950d0
GL_ROOTS(64)=0.9998248879471319145d0
GL_ROOTS(65)=-GL_ROOTS(64)
GL_ROOTS(66)=-GL_ROOTS(63)
GL_ROOTS(67)=-GL_ROOTS(62)
GL_ROOTS(68)=-GL_ROOTS(61)
GL_ROOTS(69)=-GL_ROOTS(60)
GL_ROOTS(70)=-GL_ROOTS(59)
GL_ROOTS(71)=-GL_ROOTS(58)
GL_ROOTS(72)=-GL_ROOTS(57)
GL_ROOTS(73)=-GL_ROOTS(56)
GL_ROOTS(74)=-GL_ROOTS(55)
GL_ROOTS(75)=-GL_ROOTS(54)
GL_ROOTS(76)=-GL_ROOTS(53)
GL_ROOTS(77)=-GL_ROOTS(52)
GL_ROOTS(78)=-GL_ROOTS(51)
GL_ROOTS(79)=-GL_ROOTS(50)
GL_ROOTS(80)=-GL_ROOTS(49)
GL_ROOTS(81)=-GL_ROOTS(48)
GL_ROOTS(82)=-GL_ROOTS(47)
GL_ROOTS(83)=-GL_ROOTS(46)
GL_ROOTS(84)=-GL_ROOTS(45)
GL_ROOTS(85)=-GL_ROOTS(44)
GL_ROOTS(86)=-GL_ROOTS(43)
GL_ROOTS(87)=-GL_ROOTS(42)
GL_ROOTS(88)=-GL_ROOTS(41)
GL_ROOTS(89)=-GL_ROOTS(40)
GL_ROOTS(90)=-GL_ROOTS(39)
GL_ROOTS(91)=-GL_ROOTS(38)
GL_ROOTS(92)=-GL_ROOTS(37)
GL_ROOTS(93)=-GL_ROOTS(36)
GL_ROOTS(94)=-GL_ROOTS(35)
GL_ROOTS(95)=-GL_ROOTS(34)
GL_ROOTS(96)=-GL_ROOTS(33)
GL_ROOTS(97)=-GL_ROOTS(32)
GL_ROOTS(98)=-GL_ROOTS(31)
GL_ROOTS(99)=-GL_ROOTS(30)
GL_ROOTS(100)=-GL_ROOTS(29)
GL_ROOTS(101)=-GL_ROOTS(28)
GL_ROOTS(102)=-GL_ROOTS(27)
GL_ROOTS(103)=-GL_ROOTS(26)
GL_ROOTS(104)=-GL_ROOTS(25)
GL_ROOTS(105)=-GL_ROOTS(24)
GL_ROOTS(106)=-GL_ROOTS(23)
GL_ROOTS(107)=-GL_ROOTS(22)
GL_ROOTS(108)=-GL_ROOTS(21)
GL_ROOTS(109)=-GL_ROOTS(20)
GL_ROOTS(110)=-GL_ROOTS(19)
GL_ROOTS(111)=-GL_ROOTS(18)
GL_ROOTS(112)=-GL_ROOTS(17)
GL_ROOTS(113)=-GL_ROOTS(16)
GL_ROOTS(114)=-GL_ROOTS(15)
GL_ROOTS(115)=-GL_ROOTS(14)
GL_ROOTS(116)=-GL_ROOTS(13)
GL_ROOTS(117)=-GL_ROOTS(12)
GL_ROOTS(118)=-GL_ROOTS(11)
GL_ROOTS(119)=-GL_ROOTS(10)
GL_ROOTS(120)=-GL_ROOTS(9)
GL_ROOTS(121)=-GL_ROOTS(8)
GL_ROOTS(122)=-GL_ROOTS(7)
GL_ROOTS(123)=-GL_ROOTS(6)
GL_ROOTS(124)=-GL_ROOTS(5)
GL_ROOTS(125)=-GL_ROOTS(4)
GL_ROOTS(126)=-GL_ROOTS(3)
GL_ROOTS(127)=-GL_ROOTS(2)
GL_ROOTS(128)=-GL_ROOTS(1)

END IF

IF (N_GL .EQ. 256) THEN

GL_WEIGHTS(1)=0.012247671640289755904d0
GL_WEIGHTS(2)=0.012245834369747920142d0
GL_WEIGHTS(3)=0.012242160104272800770d0
GL_WEIGHTS(4)=0.012236649395040158109d0
GL_WEIGHTS(5)=0.012229303068710278904d0
GL_WEIGHTS(6)=0.012220122227303969192d0
GL_WEIGHTS(7)=0.012209108248037240408d0
GL_WEIGHTS(8)=0.012196262783114713518d0
GL_WEIGHTS(9)=0.012181587759481772174d0
GL_WEIGHTS(10)=0.012165085378535502061d0
GL_WEIGHTS(11)=0.012146758115794459816d0
GL_WEIGHTS(12)=0.012126608720527321035d0
GL_WEIGHTS(13)=0.012104640215340463098d0
GL_WEIGHTS(14)=0.012080855895724544656d0
GL_WEIGHTS(15)=0.012055259329560149814d0
GL_WEIGHTS(16)=0.012027854356582571161d0
GL_WEIGHTS(17)=0.011998645087805811935d0
GL_WEIGHTS(18)=0.011967635904905893729d0
GL_WEIGHTS(19)=0.011934831459563562256d0
GL_WEIGHTS(20)=0.011900236672766489754d0
GL_WEIGHTS(21)=0.011863856734071078732d0
GL_WEIGHTS(22)=0.011825697100823977771d0
GL_WEIGHTS(23)=0.011785763497343426182d0
GL_WEIGHTS(24)=0.011744061914060550305d0
GL_WEIGHTS(25)=0.011700598606620740288d0
GL_WEIGHTS(26)=0.011655380094945242121d0
GL_WEIGHTS(27)=0.011608413162253105722d0
GL_WEIGHTS(28)=0.011559704854043635773d0
GL_WEIGHTS(29)=0.011509262477039497959d0
GL_WEIGHTS(30)=0.011457093598090639152d0
GL_WEIGHTS(31)=0.011403206043039185965d0
GL_WEIGHTS(32)=0.011347607895545491942d0
GL_WEIGHTS(33)=0.011290307495875509508d0
GL_WEIGHTS(34)=0.011231313439649668573d0
GL_WEIGHTS(35)=0.011170634576553449463d0
GL_WEIGHTS(36)=0.011108280009009843630d0
GL_WEIGHTS(37)=0.011044259090813901264d0
GL_WEIGHTS(38)=0.010978581425729570638d0
GL_WEIGHTS(39)=0.010911256866049039701d0
GL_WEIGHTS(40)=0.010842295511114795995d0
GL_WEIGHTS(41)=0.010771707705804626637d0
GL_WEIGHTS(42)=0.010699504038979785603d0
GL_WEIGHTS(43)=0.010625695341896561134d0
GL_WEIGHTS(44)=0.010550292686581481518d0
GL_WEIGHTS(45)=0.010473307384170403004d0
GL_WEIGHTS(46)=0.010394750983211728997d0
GL_WEIGHTS(47)=0.010314635267934015068d0
GL_WEIGHTS(48)=0.010232972256478219657d0
GL_WEIGHTS(49)=0.010149774199094865655d0
GL_WEIGHTS(50)=0.010065053576306383309d0
GL_WEIGHTS(51)=0.009978823097034910125d0
GL_WEIGHTS(52)=0.009891095696695828603d0
GL_WEIGHTS(53)=0.009801884535257327825d0
GL_WEIGHTS(54)=0.009711202995266279964d0
GL_WEIGHTS(55)=0.009619064679840727857d0
GL_WEIGHTS(56)=0.009525483410629284812d0
GL_WEIGHTS(57)=0.009430473225737752747d0
GL_WEIGHTS(58)=0.009334048377623269712d0
GL_WEIGHTS(59)=0.009236223330956302687d0
GL_WEIGHTS(60)=0.009137012760450806402d0
GL_WEIGHTS(61)=0.009036431548662873680d0
GL_WEIGHTS(62)=0.008934494783758207548d0
GL_WEIGHTS(63)=0.008831217757248750025d0
GL_WEIGHTS(64)=0.008726615961698807140d0
GL_WEIGHTS(65)=0.008620705088401014305d0
GL_WEIGHTS(66)=0.008513501025022490694d0
GL_WEIGHTS(67)=0.008405019853221535756d0
GL_WEIGHTS(68)=0.008295277846235225425d0
GL_WEIGHTS(69)=0.008184291466438269936d0
GL_WEIGHTS(70)=0.008072077362873499501d0
GL_WEIGHTS(71)=0.007958652368754348354d0
GL_WEIGHTS(72)=0.007844033498939711867d0
GL_WEIGHTS(73)=0.007728237947381555631d0
GL_WEIGHTS(74)=0.007611283084545659462d0
GL_WEIGHTS(75)=0.007493186454805883359d0
GL_WEIGHTS(76)=0.007373965773812346438d0
GL_WEIGHTS(77)=0.007253638925833913784d0
GL_WEIGHTS(78)=0.007132223961075390072d0
GL_WEIGHTS(79)=0.007009739092969822621d0
GL_WEIGHTS(80)=0.006886202695446320347d0
GL_WEIGHTS(81)=0.006761633300173798781d0
GL_WEIGHTS(82)=0.006636049593781065045d0
GL_WEIGHTS(83)=0.006509470415053660268d0
GL_WEIGHTS(84)=0.006381914752107880570d0
GL_WEIGHTS(85)=0.006253401739542401272d0
GL_WEIGHTS(86)=0.006123950655567932542d0
GL_WEIGHTS(87)=0.005993580919115338221d0
GL_WEIGHTS(88)=0.005862312086922653061d0
GL_WEIGHTS(89)=0.005730163850601437177d0
GL_WEIGHTS(90)=0.005597156033682910078d0
GL_WEIGHTS(91)=0.005463308588644310278d0
GL_WEIGHTS(92)=0.005328641593915930317d0
GL_WEIGHTS(93)=0.005193175250869280930d0
GL_WEIGHTS(94)=0.005056929880786842388d0
GL_WEIGHTS(95)=0.004919925921813865670d0
GL_WEIGHTS(96)=0.004782183925892691373d0
GL_WEIGHTS(97)=0.004643724555680060314d0
GL_WEIGHTS(98)=0.004504568581447897069d0
GL_WEIGHTS(99)=0.004364736877968056682d0
GL_WEIGHTS(100)=0.004224250421381536272d0
GL_WEIGHTS(101)=0.004083130286052668409d0
GL_WEIGHTS(102)=0.003941397641408833628d0
GL_WEIGHTS(103)=0.003799073748766257998d0
GL_WEIGHTS(104)=0.003656179958142502169d0
GL_WEIGHTS(105)=0.003512737705056307331d0
GL_WEIGHTS(106)=0.003368768507315551012d0
GL_WEIGHTS(107)=0.003224293961794198157d0
GL_WEIGHTS(108)=0.0030793357411993375832d0
GL_WEIGHTS(109)=0.0029339155908297166460d0
GL_WEIGHTS(110)=0.0027880553253277068806d0
GL_WEIGHTS(111)=0.0026417768254274905641d0
GL_WEIGHTS(112)=0.0024951020347037068508d0
GL_WEIGHTS(113)=0.0023480529563273120170d0
GL_WEIGHTS(114)=0.0022006516498399104997d0
GL_WEIGHTS(115)=0.0020529202279661431745d0
GL_WEIGHTS(116)=0.0019048808534997184044d0
GL_WEIGHTS(117)=0.0017565557363307299936d0
GL_WEIGHTS(118)=0.0016079671307493272424d0
GL_WEIGHTS(119)=0.0014591373333107332011d0
GL_WEIGHTS(120)=0.0013100886819025044578d0
GL_WEIGHTS(121)=0.0011608435575677247240d0
GL_WEIGHTS(122)=0.0010114243932084404526d0
GL_WEIGHTS(123)=0.0008618537014200890378d0
GL_WEIGHTS(124)=0.0007121541634733206669d0
GL_WEIGHTS(125)=0.0005623489540314098028d0
GL_WEIGHTS(126)=0.0004124632544261763284d0
GL_WEIGHTS(127)=0.00026253494429644590629d0
GL_WEIGHTS(128)=0.00011278901782227217551d0
GL_WEIGHTS(129)=GL_WEIGHTS(128)
GL_WEIGHTS(130)=GL_WEIGHTS(127)
GL_WEIGHTS(131)=GL_WEIGHTS(126)
GL_WEIGHTS(132)=GL_WEIGHTS(125)
GL_WEIGHTS(133)=GL_WEIGHTS(124)
GL_WEIGHTS(134)=GL_WEIGHTS(123)
GL_WEIGHTS(135)=GL_WEIGHTS(122)
GL_WEIGHTS(136)=GL_WEIGHTS(121)
GL_WEIGHTS(137)=GL_WEIGHTS(120)
GL_WEIGHTS(138)=GL_WEIGHTS(119)
GL_WEIGHTS(139)=GL_WEIGHTS(118)
GL_WEIGHTS(140)=GL_WEIGHTS(117)
GL_WEIGHTS(141)=GL_WEIGHTS(116)
GL_WEIGHTS(142)=GL_WEIGHTS(115)
GL_WEIGHTS(143)=GL_WEIGHTS(114)
GL_WEIGHTS(144)=GL_WEIGHTS(113)
GL_WEIGHTS(145)=GL_WEIGHTS(112)
GL_WEIGHTS(146)=GL_WEIGHTS(111)
GL_WEIGHTS(147)=GL_WEIGHTS(110)
GL_WEIGHTS(148)=GL_WEIGHTS(109)
GL_WEIGHTS(149)=GL_WEIGHTS(108)
GL_WEIGHTS(150)=GL_WEIGHTS(107)
GL_WEIGHTS(151)=GL_WEIGHTS(106)
GL_WEIGHTS(152)=GL_WEIGHTS(105)
GL_WEIGHTS(153)=GL_WEIGHTS(104)
GL_WEIGHTS(154)=GL_WEIGHTS(103)
GL_WEIGHTS(155)=GL_WEIGHTS(102)
GL_WEIGHTS(156)=GL_WEIGHTS(101)
GL_WEIGHTS(157)=GL_WEIGHTS(100)
GL_WEIGHTS(158)=GL_WEIGHTS(99)
GL_WEIGHTS(159)=GL_WEIGHTS(98)
GL_WEIGHTS(160)=GL_WEIGHTS(97)
GL_WEIGHTS(161)=GL_WEIGHTS(96)
GL_WEIGHTS(162)=GL_WEIGHTS(95)
GL_WEIGHTS(163)=GL_WEIGHTS(94)
GL_WEIGHTS(164)=GL_WEIGHTS(93)
GL_WEIGHTS(165)=GL_WEIGHTS(92)
GL_WEIGHTS(166)=GL_WEIGHTS(91)
GL_WEIGHTS(167)=GL_WEIGHTS(90)
GL_WEIGHTS(168)=GL_WEIGHTS(89)
GL_WEIGHTS(169)=GL_WEIGHTS(88)
GL_WEIGHTS(170)=GL_WEIGHTS(87)
GL_WEIGHTS(171)=GL_WEIGHTS(86)
GL_WEIGHTS(172)=GL_WEIGHTS(85)
GL_WEIGHTS(173)=GL_WEIGHTS(84)
GL_WEIGHTS(174)=GL_WEIGHTS(83)
GL_WEIGHTS(175)=GL_WEIGHTS(82)
GL_WEIGHTS(176)=GL_WEIGHTS(81)
GL_WEIGHTS(177)=GL_WEIGHTS(80)
GL_WEIGHTS(178)=GL_WEIGHTS(79)
GL_WEIGHTS(179)=GL_WEIGHTS(78)
GL_WEIGHTS(180)=GL_WEIGHTS(77)
GL_WEIGHTS(181)=GL_WEIGHTS(76)
GL_WEIGHTS(182)=GL_WEIGHTS(75)
GL_WEIGHTS(183)=GL_WEIGHTS(74)
GL_WEIGHTS(184)=GL_WEIGHTS(73)
GL_WEIGHTS(185)=GL_WEIGHTS(72)
GL_WEIGHTS(186)=GL_WEIGHTS(71)
GL_WEIGHTS(187)=GL_WEIGHTS(70)
GL_WEIGHTS(188)=GL_WEIGHTS(69)
GL_WEIGHTS(189)=GL_WEIGHTS(68)
GL_WEIGHTS(190)=GL_WEIGHTS(67)
GL_WEIGHTS(191)=GL_WEIGHTS(66)
GL_WEIGHTS(192)=GL_WEIGHTS(65)
GL_WEIGHTS(193)=GL_WEIGHTS(64)
GL_WEIGHTS(194)=GL_WEIGHTS(63)
GL_WEIGHTS(195)=GL_WEIGHTS(62)
GL_WEIGHTS(196)=GL_WEIGHTS(61)
GL_WEIGHTS(197)=GL_WEIGHTS(60)
GL_WEIGHTS(198)=GL_WEIGHTS(59)
GL_WEIGHTS(199)=GL_WEIGHTS(58)
GL_WEIGHTS(200)=GL_WEIGHTS(57)
GL_WEIGHTS(201)=GL_WEIGHTS(56)
GL_WEIGHTS(202)=GL_WEIGHTS(55)
GL_WEIGHTS(203)=GL_WEIGHTS(54)
GL_WEIGHTS(204)=GL_WEIGHTS(53)
GL_WEIGHTS(205)=GL_WEIGHTS(52)
GL_WEIGHTS(206)=GL_WEIGHTS(51)
GL_WEIGHTS(207)=GL_WEIGHTS(50)
GL_WEIGHTS(208)=GL_WEIGHTS(49)
GL_WEIGHTS(209)=GL_WEIGHTS(48)
GL_WEIGHTS(210)=GL_WEIGHTS(47)
GL_WEIGHTS(211)=GL_WEIGHTS(46)
GL_WEIGHTS(212)=GL_WEIGHTS(45)
GL_WEIGHTS(213)=GL_WEIGHTS(44)
GL_WEIGHTS(214)=GL_WEIGHTS(43)
GL_WEIGHTS(215)=GL_WEIGHTS(42)
GL_WEIGHTS(216)=GL_WEIGHTS(41)
GL_WEIGHTS(217)=GL_WEIGHTS(40)
GL_WEIGHTS(218)=GL_WEIGHTS(39)
GL_WEIGHTS(219)=GL_WEIGHTS(38)
GL_WEIGHTS(220)=GL_WEIGHTS(37)
GL_WEIGHTS(221)=GL_WEIGHTS(36)
GL_WEIGHTS(222)=GL_WEIGHTS(35)
GL_WEIGHTS(223)=GL_WEIGHTS(34)
GL_WEIGHTS(224)=GL_WEIGHTS(33)
GL_WEIGHTS(225)=GL_WEIGHTS(32)
GL_WEIGHTS(226)=GL_WEIGHTS(31)
GL_WEIGHTS(227)=GL_WEIGHTS(30)
GL_WEIGHTS(228)=GL_WEIGHTS(29)
GL_WEIGHTS(229)=GL_WEIGHTS(28)
GL_WEIGHTS(230)=GL_WEIGHTS(27)
GL_WEIGHTS(231)=GL_WEIGHTS(26)
GL_WEIGHTS(232)=GL_WEIGHTS(25)
GL_WEIGHTS(233)=GL_WEIGHTS(24)
GL_WEIGHTS(234)=GL_WEIGHTS(23)
GL_WEIGHTS(235)=GL_WEIGHTS(22)
GL_WEIGHTS(236)=GL_WEIGHTS(21)
GL_WEIGHTS(237)=GL_WEIGHTS(20)
GL_WEIGHTS(238)=GL_WEIGHTS(19)
GL_WEIGHTS(239)=GL_WEIGHTS(18)
GL_WEIGHTS(240)=GL_WEIGHTS(17)
GL_WEIGHTS(241)=GL_WEIGHTS(16)
GL_WEIGHTS(242)=GL_WEIGHTS(15)
GL_WEIGHTS(243)=GL_WEIGHTS(14)
GL_WEIGHTS(244)=GL_WEIGHTS(13)
GL_WEIGHTS(245)=GL_WEIGHTS(12)
GL_WEIGHTS(246)=GL_WEIGHTS(11)
GL_WEIGHTS(247)=GL_WEIGHTS(10)
GL_WEIGHTS(248)=GL_WEIGHTS(9)
GL_WEIGHTS(249)=GL_WEIGHTS(8)
GL_WEIGHTS(250)=GL_WEIGHTS(7)
GL_WEIGHTS(251)=GL_WEIGHTS(6)
GL_WEIGHTS(252)=GL_WEIGHTS(5)
GL_WEIGHTS(253)=GL_WEIGHTS(4)
GL_WEIGHTS(254)=GL_WEIGHTS(3)
GL_WEIGHTS(255)=GL_WEIGHTS(2)
GL_WEIGHTS(256)=GL_WEIGHTS(1)

GL_ROOTS(1)=0.006123912375189529501d0
GL_ROOTS(2)=0.018370818478813665118d0
GL_ROOTS(3)=0.030614968779979029366d0
GL_ROOTS(4)=0.04285452653637909838d0
GL_ROOTS(5)=0.05508765569463398410d0
GL_ROOTS(6)=0.06731252116571640024d0
GL_ROOTS(7)=0.07952728910023296590d0
GL_ROOTS(8)=0.09173012716351955203d0
GL_ROOTS(9)=0.10391920481050940364d0
GL_ROOTS(10)=0.11609269356033280494d0
GL_ROOTS(11)=0.12824876727060709474d0
GL_ROOTS(12)=0.14038560241137588591d0
GL_ROOTS(13)=0.15250137833865639537d0
GL_ROOTS(14)=0.16459427756755384983d0
GL_ROOTS(15)=0.17666248604490199740d0
GL_ROOTS(16)=0.18870419342138882646d0
GL_ROOTS(17)=0.20071759332312667007d0
GL_ROOTS(18)=0.21270088362262595794d0
GL_ROOTS(19)=0.22465226670913196715d0
GL_ROOTS(20)=0.23656994975828401848d0
GL_ROOTS(21)=0.24845214500105666683d0
GL_ROOTS(22)=0.26029706999194254198d0
GL_ROOTS(23)=0.27210294787633660951d0
GL_ROOTS(24)=0.28386800765708174180d0
GL_ROOTS(25)=0.29559048446013561456d0
GL_ROOTS(26)=0.30726861979931907626d0
GL_ROOTS(27)=0.3189006618401062756d0
GL_ROOTS(28)=0.3304848656624169762d0
GL_ROOTS(29)=0.3420194935223716365d0
GL_ROOTS(30)=0.3535028151129699895d0
GL_ROOTS(31)=0.3649331078236540185d0
GL_ROOTS(32)=0.3763086569987163903d0
GL_ROOTS(33)=0.3876277561945155836d0
GL_ROOTS(34)=0.3988887074354591277d0
GL_ROOTS(35)=0.4100898214687165500d0
GL_ROOTS(36)=0.4212294180176238250d0
GL_ROOTS(37)=0.4323058260337413100d0
GL_ROOTS(38)=0.4433173839475273572d0
GL_ROOTS(39)=0.4542624399175899988d0
GL_ROOTS(40)=0.4651393520784793136d0
GL_ROOTS(41)=0.4759464887869833064d0
GL_ROOTS(42)=0.4866822288668903501d0
GL_ROOTS(43)=0.4973449618521814771d0
GL_ROOTS(44)=0.5079330882286160362d0
GL_ROOTS(45)=0.5184450196736744762d0
GL_ROOTS(46)=0.5288791792948222620d0
GL_ROOTS(47)=0.5392340018660591811d0
GL_ROOTS(48)=0.5495079340627185570d0
GL_ROOTS(49)=0.5596994346944811451d0
GL_ROOTS(50)=0.5698069749365687591d0
GL_ROOTS(51)=0.5798290385590829449d0
GL_ROOTS(52)=0.5897641221544543008d0
GL_ROOTS(53)=0.5996107353629683217d0
GL_ROOTS(54)=0.6093674010963339395d0
GL_ROOTS(55)=0.6190326557592612194d0
GL_ROOTS(56)=0.6286050494690149754d0
GL_ROOTS(57)=0.6380831462729113687d0
GL_ROOTS(58)=0.6474655243637248626d0
GL_ROOTS(59)=0.6567507762929732219d0
GL_ROOTS(60)=0.6659375091820485599d0
GL_ROOTS(61)=0.6750243449311627639d0
GL_ROOTS(62)=0.6840099204260759531d0
GL_ROOTS(63)=0.6928928877425769601d0
GL_ROOTS(64)=0.7016719143486851594d0
GL_ROOTS(65)=0.7103456833045433134d0
GL_ROOTS(66)=0.7189128934599714484d0
GL_ROOTS(67)=0.7273722596496521266d0
GL_ROOTS(68)=0.7357225128859178346d0
GL_ROOTS(69)=0.7439624005491115685d0
GL_ROOTS(70)=0.7520906865754920596d0
GL_ROOTS(71)=0.7601061516426554549d0
GL_ROOTS(72)=0.7680075933524456360d0
GL_ROOTS(73)=0.7757938264113257391d0
GL_ROOTS(74)=0.7834636828081838208d0
GL_ROOTS(75)=0.7910160119895459945d0
GL_ROOTS(76)=0.7984496810321707588d0
GL_ROOTS(77)=0.8057635748129986233d0
GL_ROOTS(78)=0.8129565961764315431d0
GL_ROOTS(79)=0.8200276660989170674d0
GL_ROOTS(80)=0.8269757238508125143d0
GL_ROOTS(81)=0.8337997271555048943d0
GL_ROOTS(82)=0.8404986523457627139d0
GL_ROOTS(83)=0.8470714945172962072d0
GL_ROOTS(84)=0.8535172676795029651d0
GL_ROOTS(85)=0.8598350049033763507d0
GL_ROOTS(86)=0.8660237584665545193d0
GL_ROOTS(87)=0.8720825999954882891d0
GL_ROOTS(88)=0.8780106206047065440d0
GL_ROOTS(89)=0.8838069310331582849d0
GL_ROOTS(90)=0.8894706617776108888d0
GL_ROOTS(91)=0.8950009632230845774d0
GL_ROOTS(92)=0.9003970057703035448d0
GL_ROOTS(93)=0.9056579799601446471d0
GL_ROOTS(94)=0.9107830965950650119d0
GL_ROOTS(95)=0.9157715868574903845d0
GL_ROOTS(96)=0.9206227024251464955d0
GL_ROOTS(97)=0.9253357155833162029d0
GL_ROOTS(98)=0.9299099193340056412d0
GL_ROOTS(99)=0.9343446275020030943d0
GL_ROOTS(100)=0.9386391748378148050d0
GL_ROOTS(101)=0.9427929171174624432d0
GL_ROOTS(102)=0.9468052312391274814d0
GL_ROOTS(103)=0.9506755153166282764d0
GL_ROOTS(104)=0.9544031887697162418d0
GL_ROOTS(105)=0.9579876924111781294d0
GL_ROOTS(106)=0.9614284885307321440d0
GL_ROOTS(107)=0.9647250609757064309d0
GL_ROOTS(108)=0.9678769152284894549d0
GL_ROOTS(109)=0.9708835784807430293d0
GL_ROOTS(110)=0.9737445997043704053d0
GL_ROOTS(111)=0.9764595497192341556d0
GL_ROOTS(112)=0.9790280212576220388d0
GL_ROOTS(113)=0.9814496290254644058d0
GL_ROOTS(114)=0.9837240097603154962d0
GL_ROOTS(115)=0.9858508222861259565d0
GL_ROOTS(116)=0.9878297475648606089d0
GL_ROOTS(117)=0.9896604887450652183d0
GL_ROOTS(118)=0.9913427712075830869d0
GL_ROOTS(119)=0.9928763426088221171d0
GL_ROOTS(120)=0.9942609729224096650d0
GL_ROOTS(121)=0.9954964544810963566d0
GL_ROOTS(122)=0.9965826020233815404d0
GL_ROOTS(123)=0.9975192527567208276d0
GL_ROOTS(124)=0.9983062664730064441d0
GL_ROOTS(125)=0.9989435258434088566d0
GL_ROOTS(126)=0.9994309374662614082d0
GL_ROOTS(127)=0.9997684374092631861d0
GL_ROOTS(128)=0.9999560500189922307d0
GL_ROOTS(129)=-GL_ROOTS(128)
GL_ROOTS(130)=-GL_ROOTS(127)
GL_ROOTS(131)=-GL_ROOTS(126)
GL_ROOTS(132)=-GL_ROOTS(125)
GL_ROOTS(133)=-GL_ROOTS(124)
GL_ROOTS(134)=-GL_ROOTS(123)
GL_ROOTS(135)=-GL_ROOTS(122)
GL_ROOTS(136)=-GL_ROOTS(121)
GL_ROOTS(137)=-GL_ROOTS(120)
GL_ROOTS(138)=-GL_ROOTS(119)
GL_ROOTS(139)=-GL_ROOTS(118)
GL_ROOTS(140)=-GL_ROOTS(117)
GL_ROOTS(141)=-GL_ROOTS(116)
GL_ROOTS(142)=-GL_ROOTS(115)
GL_ROOTS(143)=-GL_ROOTS(114)
GL_ROOTS(144)=-GL_ROOTS(113)
GL_ROOTS(145)=-GL_ROOTS(112)
GL_ROOTS(146)=-GL_ROOTS(111)
GL_ROOTS(147)=-GL_ROOTS(110)
GL_ROOTS(148)=-GL_ROOTS(109)
GL_ROOTS(149)=-GL_ROOTS(108)
GL_ROOTS(150)=-GL_ROOTS(107)
GL_ROOTS(151)=-GL_ROOTS(106)
GL_ROOTS(152)=-GL_ROOTS(105)
GL_ROOTS(153)=-GL_ROOTS(104)
GL_ROOTS(154)=-GL_ROOTS(103)
GL_ROOTS(155)=-GL_ROOTS(102)
GL_ROOTS(156)=-GL_ROOTS(101)
GL_ROOTS(157)=-GL_ROOTS(100)
GL_ROOTS(158)=-GL_ROOTS(99)
GL_ROOTS(159)=-GL_ROOTS(98)
GL_ROOTS(160)=-GL_ROOTS(97)
GL_ROOTS(161)=-GL_ROOTS(96)
GL_ROOTS(162)=-GL_ROOTS(95)
GL_ROOTS(163)=-GL_ROOTS(94)
GL_ROOTS(164)=-GL_ROOTS(93)
GL_ROOTS(165)=-GL_ROOTS(92)
GL_ROOTS(166)=-GL_ROOTS(91)
GL_ROOTS(167)=-GL_ROOTS(90)
GL_ROOTS(168)=-GL_ROOTS(89)
GL_ROOTS(169)=-GL_ROOTS(88)
GL_ROOTS(170)=-GL_ROOTS(87)
GL_ROOTS(171)=-GL_ROOTS(86)
GL_ROOTS(172)=-GL_ROOTS(85)
GL_ROOTS(173)=-GL_ROOTS(84)
GL_ROOTS(174)=-GL_ROOTS(83)
GL_ROOTS(175)=-GL_ROOTS(82)
GL_ROOTS(176)=-GL_ROOTS(81)
GL_ROOTS(177)=-GL_ROOTS(80)
GL_ROOTS(178)=-GL_ROOTS(79)
GL_ROOTS(179)=-GL_ROOTS(78)
GL_ROOTS(180)=-GL_ROOTS(77)
GL_ROOTS(181)=-GL_ROOTS(76)
GL_ROOTS(182)=-GL_ROOTS(75)
GL_ROOTS(183)=-GL_ROOTS(74)
GL_ROOTS(184)=-GL_ROOTS(73)
GL_ROOTS(185)=-GL_ROOTS(72)
GL_ROOTS(186)=-GL_ROOTS(71)
GL_ROOTS(187)=-GL_ROOTS(70)
GL_ROOTS(188)=-GL_ROOTS(69)
GL_ROOTS(189)=-GL_ROOTS(68)
GL_ROOTS(190)=-GL_ROOTS(67)
GL_ROOTS(191)=-GL_ROOTS(66)
GL_ROOTS(192)=-GL_ROOTS(65)
GL_ROOTS(193)=-GL_ROOTS(64)
GL_ROOTS(194)=-GL_ROOTS(63)
GL_ROOTS(195)=-GL_ROOTS(62)
GL_ROOTS(196)=-GL_ROOTS(61)
GL_ROOTS(197)=-GL_ROOTS(60)
GL_ROOTS(198)=-GL_ROOTS(59)
GL_ROOTS(199)=-GL_ROOTS(58)
GL_ROOTS(200)=-GL_ROOTS(57)
GL_ROOTS(201)=-GL_ROOTS(56)
GL_ROOTS(202)=-GL_ROOTS(55)
GL_ROOTS(203)=-GL_ROOTS(54)
GL_ROOTS(204)=-GL_ROOTS(53)
GL_ROOTS(205)=-GL_ROOTS(52)
GL_ROOTS(206)=-GL_ROOTS(51)
GL_ROOTS(207)=-GL_ROOTS(50)
GL_ROOTS(208)=-GL_ROOTS(49)
GL_ROOTS(209)=-GL_ROOTS(48)
GL_ROOTS(210)=-GL_ROOTS(47)
GL_ROOTS(211)=-GL_ROOTS(46)
GL_ROOTS(212)=-GL_ROOTS(45)
GL_ROOTS(213)=-GL_ROOTS(44)
GL_ROOTS(214)=-GL_ROOTS(43)
GL_ROOTS(215)=-GL_ROOTS(42)
GL_ROOTS(216)=-GL_ROOTS(41)
GL_ROOTS(217)=-GL_ROOTS(40)
GL_ROOTS(218)=-GL_ROOTS(39)
GL_ROOTS(219)=-GL_ROOTS(38)
GL_ROOTS(220)=-GL_ROOTS(37)
GL_ROOTS(221)=-GL_ROOTS(36)
GL_ROOTS(222)=-GL_ROOTS(35)
GL_ROOTS(223)=-GL_ROOTS(34)
GL_ROOTS(224)=-GL_ROOTS(33)
GL_ROOTS(225)=-GL_ROOTS(32)
GL_ROOTS(226)=-GL_ROOTS(31)
GL_ROOTS(227)=-GL_ROOTS(30)
GL_ROOTS(228)=-GL_ROOTS(29)
GL_ROOTS(229)=-GL_ROOTS(28)
GL_ROOTS(230)=-GL_ROOTS(27)
GL_ROOTS(231)=-GL_ROOTS(26)
GL_ROOTS(232)=-GL_ROOTS(25)
GL_ROOTS(233)=-GL_ROOTS(24)
GL_ROOTS(234)=-GL_ROOTS(23)
GL_ROOTS(235)=-GL_ROOTS(22)
GL_ROOTS(236)=-GL_ROOTS(21)
GL_ROOTS(237)=-GL_ROOTS(20)
GL_ROOTS(238)=-GL_ROOTS(19)
GL_ROOTS(239)=-GL_ROOTS(18)
GL_ROOTS(240)=-GL_ROOTS(17)
GL_ROOTS(241)=-GL_ROOTS(16)
GL_ROOTS(242)=-GL_ROOTS(15)
GL_ROOTS(243)=-GL_ROOTS(14)
GL_ROOTS(244)=-GL_ROOTS(13)
GL_ROOTS(245)=-GL_ROOTS(12)
GL_ROOTS(246)=-GL_ROOTS(11)
GL_ROOTS(247)=-GL_ROOTS(10)
GL_ROOTS(248)=-GL_ROOTS(9)
GL_ROOTS(249)=-GL_ROOTS(8)
GL_ROOTS(250)=-GL_ROOTS(7)
GL_ROOTS(251)=-GL_ROOTS(6)
GL_ROOTS(252)=-GL_ROOTS(5)
GL_ROOTS(253)=-GL_ROOTS(4)
GL_ROOTS(254)=-GL_ROOTS(3)
GL_ROOTS(255)=-GL_ROOTS(2)
GL_ROOTS(256)=-GL_ROOTS(1)

END IF

RETURN

END SUBROUTINE GAUSS_LEG_ROOTS_WEIGHTS

!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!%                                                                            $
!%                                                                            $
!%                            FUNCTION TILDE_P_l                              $
!%                                                                            $
!%                                                                            $
!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!==============================================================================
!
! This function calculates ~P_l^0(COS(theta)) defined by
!
!             ~P_l^0(cos(theta))=COEFF_SH(l,0)*P_l^0(COS(theta))
!
!==============================================================================

FUNCTION TILDE_P_l(l,theta)

IMPLICIT NONE

! Input/Output declaration

INTEGER*8 :: l,m
REAL*8 :: theta,TILDE_P_l

m=0
TILDE_P_l=COEFF_SH(l,m)*POLY_LEG(l,m,DCOS(theta))

RETURN

END FUNCTION TILDE_P_l

!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!%                                                                            $
!%                                                                            $
!%                            FUNCTION TILDE_P_lm                             $
!%                                                                            $
!%                                                                            $
!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!==============================================================================
!
! This function calculates ~P_l^m(COS(theta)) defined by
!
!             ~P_l^m(cos(theta))=COEFF_SH(l,m)*P_l^m(COS(theta))
!
!==============================================================================

FUNCTION TILDE_P_lm(l,m,theta)

IMPLICIT NONE

! Input/Output declaration

INTEGER*8 :: l,m
REAL*8 :: theta,TILDE_P_lm

TILDE_P_lm=COEFF_SH(l,m)*POLY_LEG(l,ABS(m),DCOS(theta))

!IF (l==3) THEN

!PRINT*,'Coeff',COEFF_SH(l,m)
!PRINT*,'theta,l',theta,l
!PRINT*,'Poly_Leg',POLY_LEG(l,ABS(m),DCOS(theta))

!END IF

RETURN

END FUNCTION TILDE_P_lm

!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!%                                                                            $
!%                                                                            $
!%                              FUNCTION COEFF_SH                             $
!%                                                                            $
!%                                                                            $
!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!==============================================================================
!
! This function calculates the Spherical Harmonic coefficient
!
!==============================================================================

FUNCTION COEFF_SH(l,m)

IMPLICIT NONE

! Input/Ouput declaration

INTEGER*8 :: l,m
REAL*8 :: COEFF_SH

! Function Declaration

REAL*8 :: frac1,frac2,factor1

! Parameter PI

REAL*8, PARAMETER :: PI=3.14159265358979D0

! Calculating preliminary quantities involved in coefficient.

frac1=(2.D0*REAL(l)+1.D0)/(4.D0*PI)
!PRINT*,'frac1',frac1

frac2=FACTORIAL(l-ABS(m))/FACTORIAL(l+ABS(m))

factor1=frac1*frac2

! Calculating coefficient.

COEFF_SH=SQRT(factor1)

!PRINT*,'COEFF_SH',COEFF_SH
!IF(l==1)STOP

RETURN

END FUNCTION COEFF_SH

!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!%                                                                            $
!%                                                                            $
!%                              FUNCTION POLY_LEG                             $
!%                                                                            $
!%                                                                            $
!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!==============================================================================
!
! This function calculates the Associated Legendre Polynomial from Numerical 
! Recipes in Fortran 
!
!==============================================================================

FUNCTION POLY_LEG2(l,m,x)

IMPLICIT NONE

! Input/Ouput declarations

INTEGER*8 :: l,m
REAL*8 :: POLY_LEG2,x

! Function declarations

INTEGER*8 :: i,j
REAL*8 :: Pmm,Factx2,Fact,Pmmp1,Plm

! Preliminary variable test

IF (m .LT. 0 .OR. m .GT. l .OR. ABS(x) .GT. 1.D0) THEN

PRINT*, 'Bad arguments in function POLY_LEG2, module POLY_LEG2ENDRE'

END IF

Pmm=1.D0 ! Compute P_0^0

IF (m .GT. 0) THEN ! Compute P_m^m

! Factors involved in relation P_m+1^m+1=-1*(m+2)*SQRT(1-x^2)*Pmm

   Factx2=SQRT(1.D0-x**2)
   Fact=1.D0

   DO i=1,m

      Pmm=-1.D0*Fact*Factx2*Pmm
      Fact=Fact+2.D0

   END DO

END IF

IF (l .EQ. m) THEN

   POLY_LEG2=Pmm

ELSE

   Pmmp1=x*(2.D0*REAL(m)+1.D0)*Pmm ! Compute P_m+1^m=x*(2*m+1)*Pmm

   IF (l .EQ. m+1) THEN

      POLY_LEG2=Pmmp1

   ELSE

      DO j=m+2,l

         Plm=(x*(2.D0*REAL(j)-1.D0)*Pmmp1-(REAL(j)+REAL(m)-1.D0)*Pmm)/(REAL(l)-REAL(m))
         Pmm=Pmmp1
         Pmmp1=Plm

      END DO

POLY_LEG2=Plm

    END IF

END IF

RETURN

END FUNCTION POLY_LEG2

!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!%                                                                            $
!%                                                                            $
!%                              FUNCTION POLY_LEG                             $
!%                                                                            $
!%                                                                            $
!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!==============================================================================
!
! This function calculates the Associated Legendre Polynomial from Numerical 
! Recipes in Fortran 
!
!==============================================================================

FUNCTION POLY_LEG(l,m,x)

IMPLICIT NONE

    integer*8 :: l,m
    double precision :: x,POLY_LEG

    integer*8 :: i
    double precision :: p1,p2,p3,raux

    if(m<0.OR.m>l)then
      write(0,'("POLY_LEG not defined for m<0 or m>l")'&
          &)
      stop
    end if

    if(m==0.AND.l==0)then
      POLY_LEG=1
      return
    end if

    raux=sqrt(1-x**2)**m
    i=3
    do while(i<=2*m-1)
      raux=i*raux
      i=i+2
    end do

    if(l==m)then
      POLY_LEG=raux
    else
      if(l==m+1)then
        POLY_LEG=x*(2*m+1)*raux
      else
        p1=raux
        p2=x*(2*m+1)*raux
        i=m+2
        do while(i<=l)
          p3=(x*(2*i-1)*p2-(i+m-1)*p1)/(i-m)
          p1=p2
          p2=p3
          i=i+1
        end do
        POLY_LEG=p3
      end if
    end if

RETURN

END FUNCTION POLY_LEG

!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!%                                                                            $
!%                                                                            $
!%                              FUNCTION FACTORIAL                            $
!%                                                                            $
!%                                                                            $
!%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

!==============================================================================
!
! This function calculates factorial m (m!)
!
!==============================================================================

FUNCTION FACTORIAL(m)

IMPLICIT NONE

INTEGER*8 :: i,m
REAL*8 :: FACTORIAL

IF (m .EQ. 0) THEN

FACTORIAL=1.D0

ELSE

FACTORIAL=1.D0

DO i=1,m

   FACTORIAL=FACTORIAL*i

   !PRINT*,'FACTORIAL', FACTORIAL

END DO

END IF

RETURN

END FUNCTION FACTORIAL

  function spharm(l,m,t,p)
    
    implicit none
    integer*8,intent(in) :: l,m
    double precision,intent(in) :: t,p
    double precision :: spharm
    integer*8 :: i
    double precision :: raux

    raux=(2*l+1)/12.56637061435917295385d0
    do i=l-abs(m)+1,l+abs(m)
      raux=raux/i
    end do
    if(m>0)&
        spharm=sqrt(2*raux)*cos(m*p)*AssociatedLegendreFunction(l,abs(m),co&
        &s(t))
    if(m==0)&
        spharm=sqrt(raux)*AssociatedLegendreFunction(l,m,cos(t))
    if(m<0)&
        spharm=sqrt(2*raux)*sin(m*p)*AssociatedLegendreFunction(l,abs(m),co&
        &s(t))
  end function spharm


  double precision function AssociatedLegendreFunction(l,m,x)
    
    implicit none
    integer*8,intent(in) :: l,m
    double precision,intent(in) :: x

    integer*8 :: i
    double precision :: p1,p2,p3,raux

    if(m<0.OR.m>l)then
      write(0,'("AssociatedLegendreFunction not defined for m<0 or m>l")'&
          &)
      stop
    end if

    if(m==0.AND.l==0)then
      AssociatedLegendreFunction=1
      return
    end if

    raux=(-1.d0)**m*sqrt(1-x**2)**m
    i=3
    do while(i<=2*m-1)
      raux=i*raux
      i=i+2
    end do

    if(l==m)then
      AssociatedLegendreFunction=raux
    else
      if(l==m+1)then
        AssociatedLegendreFunction=x*(2*m+1)*raux
      else
        p1=raux
        p2=x*(2*m+1)*raux
        i=m+2
        do while(i<=l)
          p3=(x*(2*i-1)*p2-(i+m-1)*p1)/(i-m)
          p1=p2
          p2=p3
          i=i+1
        end do
        AssociatedLegendreFunction=p3
      end if
    end if

  end function AssociatedLegendreFunction


END MODULE MATHLIB
