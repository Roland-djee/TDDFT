PROGRAM BOX_SIZE

IMPLICIT NONE

INTEGER*8 :: i,j,N
REAL*8 :: t,ti,dt,dti,E0,Omega,T0,TR
REAL*8, DIMENSION(:), ALLOCATABLE :: KIN
REAL*8, DIMENSION(:,:), ALLOCATABLE :: X

CHARACTER(LEN=*), PARAMETER :: FMT1="(ES20.10E3,ES20.10E3,ES20.10E3,&
&ES20.10E3,ES20.10E3)"

LOGICAL :: KINETIC

REAL*8, PARAMETER :: PI=3.14159265358979D0 

OPEN(10,file='Remoteness.dat')
OPEN(20,file='Kinetic_Energy.dat')
OPEN(30,file='Field.dat')

E0=0.1D0               ! Pulse Amplitude [a.u.]
Omega=0.2D0        ! Photon Energy   [a.u.]
T0=2.D0*PI/Omega       ! Cycle duration  [a.u.]

N=100                  ! Number of time steps

dt=T0/REAL(N)          ! Time step
dti=dt/2.D0            ! Ionization Time step repeated every half-cycle
                       ! Study restrain on half-cycle.

!PRINT*,'T0,dt,dti',T0,dt,dti

! Classical equation of motion of an electron in a cosine-field

t=-dt
ti=-dti

ALLOCATE(X(N,N),KIN(N))

DO i=1,N

   t=t+dt

   DO j=1,N

      ti=ti+dti

      ! A. Zaïr PhD p.67 Remoteness

      X(i,j)=(E0/(Omega**2))*(DCOS(Omega*t)-DCOS(Omega*ti))+&
           &(E0/Omega)*DSIN(Omega*ti)*(t-ti)

      IF (ti .GT. t) THEN ! Select correct values for ti

         WRITE(10,FMT1)ti,t,-1.D0

      ELSE

         WRITE(10,FMT1)ti,t,-X(i,j)

      END IF

   END DO

   ti=-dti

END DO

ti=-dti
t=-dt

DO i=1,N

   ti=ti+dti
   t=t+dt

   ! Calculation of return time using dychotomy

   CALL RETURN_TIME(Omega,ti,T0,TR,KINETIC)

   ! Kinetic energy

   IF (KINETIC) THEN
         
      KIN(i)=((E0/Omega)*(DSIN(Omega*TR)-DSIN(Omega*ti)))**2/2.D0

   ELSE

      KIN(i)=0.D0

   END IF

   ! Kinetic energy normalized to Up
      
   IF (ti .LE. 30.D0) THEN

      WRITE(20,FMT1)ti,110.D0,KIN(i)/(E0**2/(4.D0*Omega**2))

   END IF

   ! Field
      
   WRITE(30,FMT1)30.D0,t,E0*COS(Omega*t)

END DO

! If we need a maximum value

PRINT*,'Maximum Classical Length [a.u.]',MAXVAL(X)
PRINT*,'Minimum Classical Length [a.u.]',MINVAL(X)

DEALLOCATE(X)

END PROGRAM BOX_SIZE

SUBROUTINE RETURN_TIME(Omega,ti,T0,T,KINETIC)

IMPLICIT NONE

REAL*8 :: Omega,ti,T0,T,solution

LOGICAL :: KINETIC

CALL DICHOTOMIE(Omega,ti,solution)

IF (ISNAN(solution)) THEN

   T=0.D0
   KINETIC=.FALSE.

ELSE

   T=solution
   KINETIC=.TRUE.

END IF

END SUBROUTINE RETURN_TIME

SUBROUTINE DICHOTOMIE(Omega,ti,solution)

IMPLICIT NONE

REAL*8 :: ti,x_r,x_l,x_m,f_x_l,f_x_m,solution
REAL*8 :: EVAL,Omega
REAL*8, PARAMETER :: PI=3.14159265358979D0

x_r=0.D0
x_l=2.D0*PI/Omega
!x_l=2.D0

DO WHILE(ABS(x_r-x_l) .GT. 1.D-5)

   x_m=(x_r+x_l)/2.D0

   f_x_l=EVAL(Omega,ti,x_l)
   f_x_m=EVAL(Omega,ti,x_m)

   IF (f_x_l*f_x_m .GT. 0.D0) THEN
      
      x_l=x_m

   ELSE

      x_r=x_m

   END IF

END DO

solution=x_r

RETURN

END SUBROUTINE DICHOTOMIE

FUNCTION EVAL(Omega,ti,tr)

IMPLICIT NONE

REAL*8 :: Omega,A,B,ti,tr
REAL*8 :: EVAL

A=Omega*DSIN(Omega*ti)
B=-DCOS(Omega*ti)-Omega*DSIN(Omega*ti)*ti

EVAL=DCOS(Omega*tr)+A*tr+B

!EVAL=-tr**2+1.D0

RETURN

END FUNCTION

