! ==========================================================
! satellites_linear.f95
! 
! This program calculates the control and state of a
! fully actuated satellite using optimal control with
! quadratic cost
!	I(u) = 1/2 * integ( x*Qx + u*Ru ) dt
!
! It uses a linearised model of the form
!	x' = Ax + Bu
! Presented by Yang in [Yang, Y. Analytic LQR design for
! spacecraft control system based on quaternion model.
! Journal of Aerospace Engineering, 25(3):448-453, jul 2012]
! 
! INPUT:
! - "input_file": This file is created by the user.
! It contains the simulation parameters as follows:
!  T_0: Initial time
!  T_1: Final time
!  j1 j2 j3: Moments of inertia
!  q1 q2 q3 q4 q5 q6: Elemets of diagonal matrix q
!  r1 r2 r3: Elements of diagonal matrix R
!  omega_0(1) omega_0(2), omega_0(3): Initial ang. vel.
!  omega_1(1) omega_1(2), omega_1(3): Final ang. vel.
!  eul_0(1) eul_0(2) eul_0(3): Initial Euler angles
!  eul_1(1) eul_1(2) eul_1(3): Final Euler angles
!  m: Number of intervals in the discretization
!
! OUTPUT:
! - "results_control": Control
! - "results_angveloc": Angular velocities
! - "results_eulerang": Euler angles
! - "results_quaternion": Quaternion
! - "results_time": Time vector
! - "results_misc": Miscelaneous results (Cost, CPU time...)
! 
! The convention for quaternion and
! Euler angles is as follows:
! q = (epsilon,eta)
! eul = (phi,theta,rho) = (yaw,pitch,roll)
!
! This program depends on LAPACK library: It must be
! compiled with -llapack and -lblas.
! ==========================================================

! ======== INITIALIZATION ========
! Variables declared here are passed to all subroutines
module initialization
implicit none 
real:: q1, q2, q3, q4, q5, q6 
real:: r1, r2, r3
real:: j1, j2, j3
end module initialization

! ======== MAIN PROGRAM ========
program main
use initialization
implicit none
interface
	subroutine func(n,x,t,y)
	integer, intent(in):: n
	real, intent(in):: x(n), t
	real, intent(out):: y(n)
	end subroutine func
end interface

integer, allocatable:: ipiv(:)
real, allocatable:: state_0(:), state_1(:)
real, allocatable:: state_aux_0(:), aux_solution_0_T(:)
real, allocatable:: aux_solutions(:,:), costate_0(:)
real, allocatable:: lu_fact(:,:)
real, allocatable:: state(:,:), state_costate(:,:), time(:)
integer:: n, nm, m, i, j, k
integer:: info
real:: omega_0(3), omega_1(3), eul_0(3), eul_1(3)
real:: quat_0(4), quat_1(4), eulang(3)
real:: t_0, t_1, pi, cpu_clock_ini, cpu_clock_fin
real:: cost, h
n = 6 ! Dimension of the problem
nm = 3 ! Dimension of control
pi = 4.*atan(1.)

! Open input file
open (11,file='input_file',status='old')
! Load input parameters from file
read(11,*) T_0
read(11,*) T_1
read(11,*) j1, j2, j3
read(11,*) q1, q2, q3, q4, q5, q6
read(11,*) r1, r2, r3
read(11,*) omega_0(1), omega_0(2), omega_0(3)
read(11,*) omega_1(1), omega_1(2), omega_1(3)
read(11,*) eul_0(1), eul_0(2), eul_0(3)
read(11,*) eul_1(1), eul_1(2), eul_1(3)
read(11,*) m
close(11)

allocate (state_0(n), state_1(n))
allocate (lu_fact(n,n), ipiv(n))
allocate (state_aux_0(2*n), aux_solution_0_T(n))
allocate (aux_solutions(n,n), costate_0(n))
allocate (state(m+1,n), state_costate(m+1,2*n), time(m+1))

! Here the state is (w,epsilon).
! Eta is calculated afterwards imposing |q|=1
call eul2quat(eul_0,quat_0)
call eul2quat(eul_1,quat_1)
state_0 = (/omega_0,quat_0(1:3)/)
state_1 = (/omega_1,quat_1(1:3)/)

call CPU_TIME(cpu_clock_ini)

! Solve n+1 auxiliary problems
! Problem i=0
state_aux_0 = 0. ! Construct initial condition
do k = 1,n
	state_aux_0(k) = state_0(k)
end do

! Call the solver
call rk4(func,2*n,state_aux_0,t_0,t_1,m,state_costate,time) 
do k = 1,n ! Save solution at final time
	aux_solution_0_T(k) = state_costate(m+1,k)
end do

! Problems i=1,...,N
aux_solutions = 0.
do i = 1,n
	! Construct init. condition (Canonical basis of Rn)
	state_aux_0 = 0. 
	state_aux_0(n+i) = 1.
	! Call the solver
	call rk4(func,2*n,state_aux_0,t_0,t_1, & 
	m,state_costate,time)
	do k = 1,n ! Save solution at final time
		aux_solutions(k,i) = state_costate(m+1,k)
	end do
end do

! Solve p_0 = A^-1 (x_T-x_0(T)) 
! costate_0 = inv(aux_solutions)*(state_1-aux_solution_0_T)
! LU factorization
lu_fact = aux_solutions
call sgetrf(n,n,lu_fact,n,ipiv,info) 
costate_0 =  state_1-aux_solution_0_T
! Solve linear problem
call sgetrs('n',n,1,lu_fact,n,ipiv,costate_0,n,info) 

! Final problem
do k = 1,n
	state_aux_0(k) = state_0(k)
	state_aux_0(k+n) = costate_0(k)
end do

call rk4(func,2*n,state_aux_0,t_0,t_1,m,state_costate,time)	 

call CPU_TIME(cpu_clock_fin)

cost = 0. ! Calculate cost
h = (t_1-t_0)/m
do i = 1,m+1
	cost = cost + h*(q1*state_costate(i,1)**2 & 
	+ q2*state_costate(i,2)**2 & 
	+ q3*state_costate(i,3)**2 &
	+ q4*state_costate(i,4)**2 & 
	+ q5*state_costate(i,5)**2 & 
	+ q6*state_costate(i,6)**2 &
	+ r1*(-state_costate(i,7)/(j1*r1))**2 & 
	+ r2*(-state_costate(i,8)/(j2*r2))**2 &
	+ r3*(-state_costate(i,9)/(j3*r3))**2)
end do

! Save results
open(11,file='results_misc', &
status='replace',action='write')
write(11,*) 'LINEAR MODEL;'
write(11,*) 'Initial time T_0;', t_0
write(11,*) 'Final time T_1;', t_1
write(11,*) 'Moments of inertia;', j1, j2, j3
write(11,*) 'Q (diagonal matix);', q1, q2, q3, q4, q5, q6
write(11,*) 'R (diagonal matix);', r1, r2, r3
write(11,*) 'Initial angular velocity w_0;', omega_0
write(11,*) 'Initial angular velocity w_1;', omega_1
write(11,*) 'Initial orientation (yaw,pitch,roll);', eul_0
write(11,*) 'Final orientation (yaw,pitch,roll);', eul_1
write(11,*) '# Intervals discretization;', m
write(11,*) 'Cost I(u) = ', cost
write(11,*) 'CPU time (s);', cpu_clock_fin-cpu_clock_ini
close(11)

open(11,file='results_control_LQR', &
status='replace',action='write')
do i = 1,m+1
	write(11,*) -state_costate(i,7)/(j1*r1), & 
		    -state_costate(i,8)/(j2*r2), &
		    -state_costate(i,9)/(j3*r3)
end do
close(11)

open(11,file='results_control', &
status='replace',action='write')
do i = 1,m 
	write(11,*) &
	(state_costate(i+1,1) - & 
	state_costate(i,1))*m/(T_1-T_0), &
	(state_costate(i+1,2) - & 
	state_costate(i,2))*m/(T_1-T_0), &
	(state_costate(i+1,3) - & 
	state_costate(i,3))*m/(T_1-T_0)
end do
close(11)

open(11,file='results_quaternion', &
status='replace',action='write')
do i = 1,m+1
	write(11,*) & 
	state_costate(i,4:6), & 
	(1-state_costate(i,4)**2 & 
	-state_costate(i,5)**2 &
	-state_costate(i,6)**2)**.5
end do
close(11)

open(11,file='results_eulerang', &
status='replace',action='write')
do i = 1,m+1
	call quat2eul((/state_costate(i,4:6), &
	(1-state_costate(i,4)**2 &
	-state_costate(i,5)**2 &
	-state_costate(i,6)**2)**.5/) &
	,eulang)
	write (11,*) eulang
end do
close(11)

open(11,file='results_angveloc', &
status='replace',action='write')
do i = 1,m+1
	write (11,*) state_costate(i,1:3)
end do
close(11)

open(11, file='results_time', &
status='replace',action='write')
do i = 1, m+1
        write (11,*) T_0+(T_1-T_0)*(i-1)/m
end do
close (11)

deallocate (ipiv)
deallocate (state_0, state_1)
deallocate (state_aux_0, aux_solution_0_T)
deallocate (aux_solutions, costate_0)
deallocate (lu_fact)
deallocate (state, state_costate, time)

end program main

! ======== 4TH ORDER RUNGE-KUTTA METHOD ========
subroutine rk4(func,n,x0,t0,t1,m,x,time)
implicit none
integer, intent(in):: n, m
real, intent(in):: x0(n), t0, t1
real, intent(out):: x(m+1,n), time(m+1)
real:: h, a1(N), a2(N), a3(N), a4(N)
integer:: i, k
h = (t1-t0)/m

do i = 1,n
	x(1,i) = x0(i)
end do

do k = 1,m
	time(k) = t0 + (k-1)*h
	call func(n,x(k,:),time(k),a1)
	a1 = h*a1
	call func(n,x(k,:)+.5*a1,time(k)+h/2.,a2)
	a2 = h*a2
	call func(n,x(k,:)+.5*a2,time(k)+h/2.,a3)
	a3 = h*a3 
	call func(n,x(k,:)+a3,time(k)+h,a4)
	a4 = h*a4
	do i = 1,n
		x(k+1,i) = x(k,i) &
		+ (a1(i) + 2*a2(i)  +2*a3(i) + a4(i))/6.
	end do
end do
time(m+1) = t1
end subroutine rk4

! ======== DYNAMICS ======== 
subroutine func(n,x,t,y)
use initialization
implicit none
integer, intent(in):: n
real, intent(in):: x(n), t
real, intent(out):: y(n)
y(1) =  -x(7)/(r1*j1**2)
y(2) =  -x(8)/(r2*j2**2)
y(3) =  -x(9)/(r3*j3**2)
y(4) =   x(1)*.5
y(5) =   x(2)*.5
y(6) =   x(3)*.5
y(7) =  -x(10)*.5 - q1*x(1)
y(8) =  -x(11)*.5 - q2*x(2)
y(9) =  -x(12)*.5 - q3*x(3)
y(10) = -q4*x(4)
y(11) = -q5*x(5)
y(12) = -q6*x(6)
end subroutine func

! ======== EULER ANGLES TO QUATERNION ========
subroutine eul2quat(eul,q)
real, intent(in):: eul(3)
real, intent(out):: q(4)
real:: cy, sy, cp, sp, cr, sr
cy = cos(eul(1)*.5);
sy = sin(eul(1)*.5);
cp = cos(eul(2)*.5);
sp = sin(eul(2)*.5);
cr = cos(eul(3)*.5);
sr = sin(eul(3)*.5);
q(1) = cy * cp * sr - sy * sp * cr
q(2) = sy * cp * sr + cy * sp * cr
q(3) = sy * cp * cr - cy * sp * sr
q(4) = cy * cp * cr + sy * sp * sr
end subroutine eul2quat

! ======== QUATERNION TO EULER ANGLES ========
subroutine quat2eul(q,eul)
real, intent(in):: q(4)
real, intent(out):: eul(3)
eul(1) = atan2(2.*(q(4)*q(3)+q(1)*q(2)), & 
1.-2.*(q(2)**2.+q(3)**2.))
eul(2) = asin(2.*(q(4)*q(2)-q(3)*q(1)))
eul(3) = atan2(2.*(q(4)*q(1)+q(2)*q(3)), &
1.-2.*(q(1)**2.+q(2)**2.))
end subroutine quat2eul
