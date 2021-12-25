! ==========================================================
! satellites_nonlinear.f95
! 
! This program calculates the control and state of a
! fully actuated satellite using optimal control with
! quadratic cost
!	I(u) = 1/2 * integ( x*Qx + u*Ru ) dt
!
! And constraints on control
! 	ui_min <= u_i <= ui_max,  i=1,2,3.
!
! It uses the fully non-linear model
!	Jw' = -w x Jw + u
!	q' = 1/2 * Omega(w) * q
! 
! INPUT:
! - "input_file": This file is created by the user.
! It contains the simulation parameters as follows:
!  T_0: Initial time
!  T_1: Final time
!  j1 j2 j3: Moments of inertia
!  q1 q2 q3 q4 q5 q6 q7: Elemets of diagonal matrix q
!  r1 r2 r3: Elements of diagonal matrix R
!  omega_0(1) omega_0(2), omega_0(3): Initial ang. vel.
!  omega_1(1) omega_1(2), omega_1(3): Final ang. vel.
!  eul_0(1) eul_0(2) eul_0(3): Initial Euler angles
!  eul_1(1) eul_1(2) eul_1(3): Final Euler angles
!  m: Number of intervals in the discretization
!  u1min u2min u3min u1max u2max u3max: Control constraints
!  mask(1) ... mask(7): Natural/Essential bound. constraints
!  prec: Precission of the algorithm
!
! OUTPUT:
! - "results_control": Control
! - "results_angveloc": Angular velocities
! - "results_eulerang": Euler angles
! - "results_quaternion": Quaternion
! - "results_quatnorm": Norm of quaternion
! - "results_time": Time vector
! - "results_converg": Convergence
! - "results_misc": Miscelaneous results (Cost, CPU time...)
!
! The convention for quaternion and
! Euler angles is as follows:
! q = (epsilon,eta)
! eul = (phi,theta,rho) = (yaw,pitch,roll)
! ==========================================================

! ======== INITIALIZATION ========
! Variables declares here are passed to all rubroutines
module initialization
implicit none
integer:: N, m, nn
real:: T_0, T_1
real:: j1, j2, j3
real:: q1, q2, q3, q4, q5, q6, q7
real:: r1, r2, r3
real:: u1min, u1max, u2min, u2max, u3min, u3max
end module initialization

! ======== MAIN ========
program main
use initialization
implicit none

integer, allocatable:: mask(:)
real, allocatable:: state_0(:), state_1(:), state(:, :)
integer:: i
real:: cpu_clock_ini, cpu_clock_fin
real:: prec, cost, eulang(3), pi
real:: omega_0(3), omega_1(3)
real:: eul_0(3), eul_1(3), quat_0(4), quat_1(4)
N=7; ! Dimension of state
nn=14;  ! Number of constraints
pi = 4*atan(1.)
allocate (state_0(N), state_1(N), mask(N))

! Open input file
open (11,file='input_file',status='old')
! Read input parameters
read(11,*) T_0
read(11,*) T_1
read(11,*) j1, j2, j3
read(11,*) q1, q2, q3, q4, q5, q6, q7
read(11,*) r1, r2, r3
read(11,*) omega_0(1), omega_0(2), omega_0(3)
read(11,*) omega_1(1), omega_1(2), omega_1(3)
read(11,*) eul_0(1), eul_0(2), eul_0(3)
read(11,*) eul_1(1), eul_1(2), eul_1(3)
read(11,*) m
read(11,*) u1min, u2min, u3min, u1max, u2max, u3max
read(11,*) mask(1), mask(2), mask(3), mask(4),&
	 mask(5), mask(6), mask(7)
read(11,*) prec
close(11)

allocate (state(m+1,N))

call eul2quat(eul_0,quat_0)
call eul2quat(eul_1,quat_1)
state_0 = (/omega_0,quat_0/)
state_1 = (/omega_1,quat_1/)

do i=1, m+1 
	state(i, :)=((m+1.-i)/m)*state_0+((i-1.)/m)*state_1
end do

call CPU_TIME(cpu_clock_ini)
call optimalcontrol(prec, mask, state, cost)
call CPU_TIME(cpu_clock_fin)

open(11,file='results_misc', &
status='replace',action='write')
write(11,*) 'NONLINEAR MODEL;'
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
write(11,*) 'Precision;', prec
write(11,*) 'Mask (1=bound.cond./0=natural);', mask
write(11,*) 'Cost I(u);', cost
write(11,*) 'CPU time (s);', cpu_clock_fin-cpu_clock_ini
close(11)

open(11,file='results_control', &
status='replace',action='write')
do i = 1,m 
	write(11,*) &
	(state(i+1,1) - state(i,1))*m/(T_1-T_0), &
	(state(i+1,2) - state(i,2))*m/(T_1-T_0), &
	(state(i+1,3) - state(i,3))*m/(T_1-T_0)
end do
close(11)

open(11,file='results_quaternion', &
status='replace',action='write')
do i = 1,m+1
	write(11,*) state(i,4:7)
end do
close(11)
open(11,file='results_eulerang', &
status='replace',action='write')

do i = 1,m+1
	call quat2eul(state(i,4:7),eulang)
	write (11,*) eulang
end do
close(11)

open(11,file='results_angveloc', &
status='replace',action='write')
do i = 1,m+1
	write (11,*) state(i,1:3)
end do
close(11)

open(11,file='results_quatnorm', &
status='replace',action='write')
do i = 1, m+1
        write (11,*) & 
	state(i,4)**2+state(i,5)**2+ &
	state(i,6)**2+state(i,7)**2
end do
close (11)

open(11,file='results_time', &
status='replace',action='write')
do i = 1, m+1
        write (11,*) T_0+(T_1-T_0)*(i-1)/m
end do
close (11)

deallocate (state_0, state_1, state, mask)
end program main

! ======== CONSTRAINTS ========
subroutine constraints(s, x, p, R)
use initialization
implicit none
real:: s, x(N), p(N), R(nn)
R (1)  =  j1*p(1) + (j3-j2)*(x(2)*x(3)) - u1max
R (2)  = -j1*p(1) - (j3-j2)*(x(2)*x(3)) + u1min
R (3)  =  j2*p(2) + (j1-j3)*(x(1)*x(3)) - u2max
R (4)  = -j2*p(2) - (j1-j3)*(x(1)*x(3)) + u2min
R (5)  =  j3*p(3) + (j2-j1)*(x(1)*x(2)) - u3max
R (6)  = -j3*p(3) - (j2-j1)*(x(1)*x(2)) + u3min
R (7)  =  p(4) - .5*(+ x(1)*x(7) - x(2)*x(6) + x(3)*x(5))
R (8)  = -p(4) + .5*(+ x(1)*x(7) - x(2)*x(6) + x(3)*x(5))
R (9)  =  p(5) - .5*(+ x(1)*x(6) - x(3)*x(4) + x(2)*x(7))
R (10) = -p(5) + .5*(+ x(1)*x(6) - x(3)*x(4) + x(2)*x(7))
R (11) =  p(6) - .5*(+ x(2)*x(4) - x(1)*x(5) + x(3)*x(7))
R (12) = -p(6) + .5*(+ x(2)*x(4) - x(1)*x(5) + x(3)*x(7))
R (13) =  p(7) - .5*(- x(1)*x(4) - x(2)*x(5) - x(3)*x(6))
R (14) = -p(7) + .5*(- x(1)*x(4) - x(2)*x(5) - x(3)*x(6))
end subroutine constraints

! ======== INTEGRAND ========
subroutine integrand(s, x, p, y, F)
use initialization
implicit none
real, intent(in):: s, x(N), p(N), y(m, nn)
real, intent(out):: F
real:: l(nn), R(nn) 
integer:: i
call constraints(s, x, p, R)
call indexfunction(s, y, l)

F =     (r1*(j1*p(1) - j2*x(2)*x(3) + j3*x(2)*x(3))**2)/2  
F = F + (r2*(j2*p(2) + j1*x(1)*x(3) - j3*x(1)*x(3))**2)/2
F = F + (r3*(j3*p(3) - j1*x(1)*x(2) + j2*x(1)*x(2))**2)/2

F = F + (q1*x(1)**2)/2  
F = F + (q2*x(2)**2)/2
F = F + (q3*x(3)**2)/2
F = F + (q4*x(4)**2)/2
F = F + (q5*x(5)**2)/2
F = F + (q6*x(6)**2)/2
F = F + (q7*x(7)**2)/2
F = F - nn
do i = 1,nn
        F = F + exp(l(i)*R(i))
end do
end subroutine integrand

! ======== AUX INTEGRAND 1 ========
subroutine auxintegrand1(s, x, p, y, H)
use initialization
implicit none
real, intent(in)::  s, x(N), p(N), y(m, nn)
real, intent(out):: H(N)
real:: l(nn)!, R(nn)
!call constraints(s, x, p, R)
call indexfunction(s, y, l)
H(1)=q1*x(1)-(l(7)*x(7)*exp(l(7)*(p(4)-&
	(x(1)*x(7))/2+(x(2)*x(6))/2-(x(3)*x(5))/2)))/2
H(1)=H(1)+(l(8)*x(7)*exp(-l(8)*(p(4)-&
	(x(1)*x(7))/2+(x(2)*x(6))/2-(x(3)*x(5))/2)))/2
H(1)=H(1)-(l(9)*x(6)*exp(l(9)*(p(5)-&
	(x(1)*x(6))/2+(x(3)*x(4))/2-(x(2)*x(7))/2)))/2
H(1)=H(1)+(l(10)*x(6)*exp(-l(10)*(p(5)-&
	(x(1)*x(6))/2+(x(3)*x(4))/2-(x(2)*x(7))/2)))/2
H(1)=H(1)+(l(11)*x(5)*exp(l(11)*(p(6)+&
	(x(1)*x(5))/2-(x(2)*x(4))/2-(x(3)*x(7))/2)))/2
H(1)=H(1)-(l(12)*x(5)*exp(-l(12)*(p(6)+&
	(x(1)*x(5))/2-(x(2)*x(4))/2-(x(3)*x(7))/2)))/2
H(1)=H(1)+(l(13)*x(4)*exp(l(13)*(p(7)+&
	(x(1)*x(4))/2+(x(2)*x(5))/2+(x(3)*x(6))/2)))/2
H(1)=H(1)-(l(14)*x(4)*exp(-l(14)*(p(7)+&
	(x(1)*x(4))/2+(x(2)*x(5))/2+(x(3)*x(6))/2)))/2
H(1)=H(1)-r3*(j1*x(2)-j2*x(2))*&
	(j3*p(3)-j1*x(1)*x(2)+j2*x(1)*x(2))
H(1)=H(1)+r2*(j1*x(3)-j3*x(3))*&
	(j2*p(2)+j1*x(1)*x(3)-j3*x(1)*x(3))
H(1)=H(1)+l(3)*x(3)*&
	exp(l(3)*(j2*p(2)+x(1)*x(3)*(j1-j3)-1))*(j1-j3)
H(1)=H(1)-l(4)*x(3)*&
	exp(-l(4)*(j2*p(2)+x(1)*x(3)*(j1-j3)+1))*(j1-j3)
H(1)=H(1)-l(5)*x(2)*&
	exp(-l(5)*(x(1)*x(2)*(j1-j2)-j3*p(3)+1))*(j1-j2)
H(1)=H(1)+l(6)*x(2)*&
	exp(-l(6)*(j3*p(3)-x(1)*x(2)*(j1-j2)+1))*(j1-j2)

H(2)=q2*x(2)+(l(7)*x(6)*exp(l(7)*(p(4)-&
	(x(1)*x(7))/2+(x(2)*x(6))/2-(x(3)*x(5))/2)))/2
H(2)=H(2)-(l(8)*x(6)*exp(-l(8)*(p(4)-&
	(x(1)*x(7))/2+(x(2)*x(6))/2-(x(3)*x(5))/2)))/2
H(2)=H(2)-(l(9)*x(7)*exp(l(9)*(p(5)-&
	(x(1)*x(6))/2+(x(3)*x(4))/2-(x(2)*x(7))/2)))/2
H(2)=H(2)-(l(11)*x(4)*exp(l(11)*(p(6)+&
	(x(1)*x(5))/2-(x(2)*x(4))/2-(x(3)*x(7))/2)))/2
H(2)=H(2)+(l(10)*x(7)*exp(-l(10)*(p(5)-&
	(x(1)*x(6))/2+(x(3)*x(4))/2-(x(2)*x(7))/2)))/2
H(2)=H(2)+(l(12)*x(4)*exp(-l(12)*(p(6)+&
	(x(1)*x(5))/2-(x(2)*x(4))/2-(x(3)*x(7))/2)))/2
H(2)=H(2)+(l(13)*x(5)*exp(l(13)*(p(7)+&
	(x(1)*x(4))/2+(x(2)*x(5))/2+(x(3)*x(6))/2)))/2
H(2)=H(2)-(l(14)*x(5)*exp(-l(14)*(p(7)+&
	(x(1)*x(4))/2+(x(2)*x(5))/2+(x(3)*x(6))/2)))/2
H(2)=H(2)-r3*(j1*x(1)-j2*x(1))*&
	(j3*p(3)-j1*x(1)*x(2)+j2*x(1)*x(2))
H(2)=H(2)-r1*(j2*x(3)-j3*x(3))*&
	(j1*p(1)-j2*x(2)*x(3)+j3*x(2)*x(3))
H(2)=H(2)-l(1)*x(3)*&
	exp(-l(1)*(x(2)*x(3)*(j2-j3)-j1*p(1)+1))*(j2-j3)
H(2)=H(2)+l(2)*x(3)*&
	exp(-l(2)*(j1*p(1)-x(2)*x(3)*(j2-j3)+1))*(j2-j3)
H(2)=H(2)-l(5)*x(1)*&
	exp(-l(5)*(x(1)*x(2)*(j1-j2)-j3*p(3)+1))*(j1-j2)
H(2)=H(2)+l(6)*x(1)*&
	exp(-l(6)*(j3*p(3)-x(1)*x(2)*(j1-j2)+1))*(j1-j2)

H(3)=q3*x(3)-(l(7)*x(5)*exp(l(7)*(p(4)-&
	(x(1)*x(7))/2+(x(2)*x(6))/2-(x(3)*x(5))/2)))/2
H(3)=H(3)+(l(8)*x(5)*exp(-l(8)*(p(4)-&
	(x(1)*x(7))/2+(x(2)*x(6))/2-(x(3)*x(5))/2)))/2
H(3)=H(3)+(l(9)*x(4)*exp(l(9)*(p(5)-&
	(x(1)*x(6))/2+(x(3)*x(4))/2-(x(2)*x(7))/2)))/2
H(3)=H(3)-(l(10)*x(4)*exp(-l(10)*(p(5)-&
	(x(1)*x(6))/2+(x(3)*x(4))/2-(x(2)*x(7))/2)))/2
H(3)=H(3)-(l(11)*x(7)*exp(l(11)*(p(6)+&
	(x(1)*x(5))/2-(x(2)*x(4))/2-(x(3)*x(7))/2)))/2
H(3)=H(3)+(l(12)*x(7)*exp(-l(12)*(p(6)+&
	(x(1)*x(5))/2-(x(2)*x(4))/2-(x(3)*x(7))/2)))/2
H(3)=H(3)+(l(13)*x(6)*exp(l(13)*(p(7)+&
	(x(1)*x(4))/2+(x(2)*x(5))/2+(x(3)*x(6))/2)))/2
H(3)=H(3)-(l(14)*x(6)*exp(-l(14)*(p(7)+&
	(x(1)*x(4))/2+(x(2)*x(5))/2+(x(3)*x(6))/2)))/2
H(3)=H(3)+r2*(j1*x(1)-j3*x(1))*&
	(j2*p(2)+j1*x(1)*x(3)-j3*x(1)*x(3))
H(3)=H(3)-r1*(j2*x(2)-j3*x(2))*&
	(j1*p(1)-j2*x(2)*x(3)+j3*x(2)*x(3))
H(3)=H(3)-l(1)*x(2)*&
	exp(-l(1)*(x(2)*x(3)*(j2-j3)-j1*p(1)+1))*(j2-j3)
H(3)=H(3)+l(3)*x(1)*&
	exp(l(3)*(j2*p(2)+x(1)*x(3)*(j1-j3)-1))*(j1-j3)
H(3)=H(3)+l(2)*x(2)*&
	exp(-l(2)*(j1*p(1)-x(2)*x(3)*(j2-j3)+1))*(j2-j3)
H(3)=H(3)-l(4)*x(1)*&
	exp(-l(4)*(j2*p(2)+x(1)*x(3)*(j1-j3)+1))*(j1-j3)

H(4)=q4*x(4)+(l(9)*x(3)*exp(l(9)*(p(5)-&
	(x(1)*x(6))/2+(x(3)*x(4))/2-(x(2)*x(7))/2)))/2
H(4)=H(4)-(l(10)*x(3)*exp(-l(10)*(p(5)-&
	(x(1)*x(6))/2+(x(3)*x(4))/2-(x(2)*x(7))/2)))/2
H(4)=H(4)-(l(11)*x(2)*exp(l(11)*(p(6)+&
	(x(1)*x(5))/2-(x(2)*x(4))/2-(x(3)*x(7))/2)))/2
H(4)=H(4)+(l(12)*x(2)*exp(-l(12)*(p(6)+&
	(x(1)*x(5))/2-(x(2)*x(4))/2-(x(3)*x(7))/2)))/2
H(4)=H(4)+(l(13)*x(1)*exp(l(13)*(p(7)+&
	(x(1)*x(4))/2+(x(2)*x(5))/2+(x(3)*x(6))/2)))/2
H(4)=H(4)-(l(14)*x(1)*exp(-l(14)*(p(7)+&
	(x(1)*x(4))/2+(x(2)*x(5))/2+(x(3)*x(6))/2)))/2

H(5)=q5*x(5)-(l(7)*x(3)*exp(l(7)*(p(4)-&
	(x(1)*x(7))/2+(x(2)*x(6))/2-(x(3)*x(5))/2)))/2
H(5)=H(5)+(l(8)*x(3)*exp(-l(8)*(p(4)-&
	(x(1)*x(7))/2+(x(2)*x(6))/2-(x(3)*x(5))/2)))/2
H(5)=H(5)+(l(11)*x(1)*exp(l(11)*(p(6)+&
	(x(1)*x(5))/2-(x(2)*x(4))/2-(x(3)*x(7))/2)))/2
H(5)=H(5)-(l(12)*x(1)*exp(-l(12)*(p(6)+&
	(x(1)*x(5))/2-(x(2)*x(4))/2-(x(3)*x(7))/2)))/2
H(5)=H(5)+(l(13)*x(2)*exp(l(13)*(p(7)+&
	(x(1)*x(4))/2+(x(2)*x(5))/2+(x(3)*x(6))/2)))/2
H(5)=H(5)-(l(14)*x(2)*exp(-l(14)*(p(7)+&
	(x(1)*x(4))/2+(x(2)*x(5))/2+(x(3)*x(6))/2)))/2

H(6)=q6*x(6)+(l(7)*x(2)*exp(l(7)*(p(4)-&
	(x(1)*x(7))/2+(x(2)*x(6))/2-(x(3)*x(5))/2)))/2
H(6)=H(6)-(l(8)*x(2)*exp(-l(8)*(p(4)-&
	(x(1)*x(7))/2+(x(2)*x(6))/2-(x(3)*x(5))/2)))/2
H(6)=H(6)-(l(9)*x(1)*exp(l(9)*(p(5)-&
	(x(1)*x(6))/2+(x(3)*x(4))/2-(x(2)*x(7))/2)))/2
H(6)=H(6)+(l(10)*x(1)*exp(-l(10)*(p(5)-&
	(x(1)*x(6))/2+(x(3)*x(4))/2-(x(2)*x(7))/2)))/2
H(6)=H(6)+(l(13)*x(3)*exp(l(13)*(p(7)+&
	(x(1)*x(4))/2+(x(2)*x(5))/2+(x(3)*x(6))/2)))/2
H(6)=H(6)-(l(14)*x(3)*exp(-l(14)*(p(7)+&
	(x(1)*x(4))/2+(x(2)*x(5))/2+(x(3)*x(6))/2)))/2

H(7)=q7*x(7)-(l(7)*x(1)*exp(l(7)*(p(4)-&
	(x(1)*x(7))/2+(x(2)*x(6))/2-(x(3)*x(5))/2)))/2
H(7)=H(7)+(l(8)*x(1)*exp(-l(8)*(p(4)-&
	(x(1)*x(7))/2+(x(2)*x(6))/2-(x(3)*x(5))/2)))/2
H(7)=H(7)-(l(9)*x(2)*exp(l(9)*(p(5)-&
	(x(1)*x(6))/2+(x(3)*x(4))/2-(x(2)*x(7))/2)))/2
H(7)=H(7)+(l(10)*x(2)*exp(-l(10)*(p(5)-&
	(x(1)*x(6))/2+(x(3)*x(4))/2-(x(2)*x(7))/2)))/2
H(7)=H(7)-(l(11)*x(3)*exp(l(11)*(p(6)+&
	(x(1)*x(5))/2-(x(2)*x(4))/2-(x(3)*x(7))/2)))/2
H(7)=H(7)+(l(12)*x(3)*exp(-l(12)*(p(6)+&
	(x(1)*x(5))/2-(x(2)*x(4))/2-(x(3)*x(7))/2)))/2
end subroutine auxintegrand1

! ======== AUX INTEGRAND 2 ========
subroutine auxintegrand2(s, x, p, y, T, G)
use initialization
implicit none
real, intent(in):: s, x(N), p(N), y(m, nn), T
real, intent(out):: G(N)
integer:: i
real:: l(nn), H(N)!, R(nn)
!call constraints(s, x, p, R)
call indexfunction(s, y, l)
call auxintegrand1(s, x, p, y, H)
G(1)=j1*r1*(j1*p(1)-j2*x(2)*x(3)+j3*x(2)*x(3))
G(1)=G(1)+j1*l(1)*exp(-l(1)*(x(2)*x(3)*(j2-j3)-j1*p(1)+1))
G(1)=G(1)-j1*l(2)*exp(-l(2)*(j1*p(1)-x(2)*x(3)*(j2-j3)+1))
G(1)=G(1)-(T-s)*H(1)

G(2)=j2*r2*(j2*p(2)+j1*x(1)*x(3)-j3*x(1)*x(3))
G(2)=G(2)+j2*l(3)*exp(l(3)*(j2*p(2)+x(1)*x(3)*(j1-j3)-1))
G(2)=G(2)-j2*l(4)*exp(-l(4)*(j2*p(2)+x(1)*x(3)*(j1-j3)+1))
G(2)=G(2)-(T-s)*H(2)

G(3)=j3*r3*(j3*p(3)-j1*x(1)*x(2)+j2*x(1)*x(2))
G(3)=G(3)+j3*l(5)*exp(-l(5)*(x(1)*x(2)*(j1-j2)-j3*p(3)+1))
G(3)=G(3)-j3*l(6)*exp(-l(6)*(j3*p(3)-x(1)*x(2)*(j1-j2)+1))
G(3)=G(3)-(T-s)*H(3)

G(4)=l(7)*exp(l(7)*(p(4)-(x(1)*x(7))/2+&
	(x(2)*x(6))/2-(x(3)*x(5))/2))
G(4)=G(4)-l(8)*exp(-l(8)*(p(4)-(x(1)*x(7))/2+&
	(x(2)*x(6))/2-(x(3)*x(5))/2))
G(4)=G(4)-(T-s)*H(4)

G(5)=l(9)*exp(l(9)*(p(5)-(x(1)*x(6))/2+&
	(x(3)*x(4))/2-(x(2)*x(7))/2))
G(5)=G(5)-l(10)*exp(-l(10)*(p(5)-(x(1)*x(6))/2+&
	(x(3)*x(4))/2-(x(2)*x(7))/2))
G(5)=G(5)-(T-s)*H(5)

G(6)=l(11)*exp(l(11)*(p(6)+(x(1)*x(5))/2-&
	(x(2)*x(4))/2-(x(3)*x(7))/2))
G(6)=G(6)-l(12)*exp(-l(12)*(p(6)+(x(1)*x(5))/2-&
	(x(2)*x(4))/2-(x(3)*x(7))/2))
G(6)=G(6)-(T-s)*H(6)

G(7)=l(13)*exp(l(13)*(p(7)+(x(1)*x(4))/2+&
	(x(2)*x(5))/2+(x(3)*x(6))/2))
G(7)=G(7)-l(14)*exp(-l(14)*(p(7)+(x(1)*x(4))/2+&
	(x(2)*x(5))/2+(x(3)*x(6))/2))
G(7)=G(7)-(T-s)*H(7)
end subroutine auxintegrand2

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

! ======== INDEX FUNCTION ========
subroutine indexfunction(s, z, L)
use initialization
implicit none
real:: s, z(m, nn), L(nn)
integer:: k
k=int((s-T_0)*m/(T_1-T_0))+1
L=z(k, :)
end subroutine indexfunction

! ======== OPTIMAL CONTROL ========
subroutine optimalcontrol(prec, mask, state, cost)
use initialization
implicit none
external integrand
integer:: i, j, k, mask(N)
real:: deriv(N), midpoint(N), r(nn), y(m, nn)
real:: convergence, prec, state(m+1, N), hh
real:: cost_prev, cost
real:: H(m, N), y_0, G_aux(m, N), G_auxx(m, N)
real:: G(m, N), integ(N),norm, aux(m, N), tam, eta, F

open(11,file='results_converg', &
status='replace',action='write') 

hh = (T_1-T_0)/m
y_0 = 1.
y = 1.
eta = 1.e-01; 

! // Optimal control algorithm ========
do 

convergence=0.
do i=1, m
	deriv=(state(i+1, :)-state(i, :))/hh
	midpoint=0.5*(state(i, :)+state(i+1, :))
	call &
	constraints(T_0+hh*(i-0.5), midpoint, deriv, r)
	do k=1, nn
		if (convergence<abs(r(k)*y(i, k))) then
			convergence=abs(r(k)*y(i, k))
		endif
	end do
end do
print*, 'convergence ', convergence
write(11,*) convergence
if (convergence<prec) exit

cost_prev=1.e+10

do
	cost=0.
	integ=0.
	do i=1, m
		midpoint=0.5*(state(i,:)+state(i+1,:))
		deriv=(state(i+1,:)-state(i,:))/hh
		call auxintegrand1(T_0+hh*(i-0.5), &
		midpoint,deriv,y,H(i,:))
		call integrand(T_0+hh*(i-0.5), &
		midpoint,deriv,y,F)
		cost=cost+hh*F
		call auxintegrand2(T_0+hh*(i-0.5), &
		midpoint,deriv,y,T_1,G_aux(i,:))
		call auxintegrand2(T_0+hh*(i-0.5), &
		midpoint,deriv,y,T_0,G_auxx(i,:))
		do j=1, N
			if (mask(j).eq.1) then
				G(i,j)=G_aux(i,j)
				integ(j)=integ(j)+hh*G(i,j)
			else
				G(i,j)=G_auxx(i,j)
			end if
		end do
	end do
	cost_prev=cost
	norm=0.
	do i=1, m
		aux(i,:)=0.
		do k=1, i
			do j=1, N
				if (mask(j).eq.1) then
					aux(i,j)=&
					aux(i,j)+H(k,j)
				else
					aux(i,:)=&
					aux(i,:)-hh*G(k,:)
				end if
			end do
		end do
		do j=1, N
			if (mask(j).eq.1) then
				aux(i,j)=&
				aux(i,j)*hh*(hh*i-(T_1-T_0))
			end if
		end do
		do k=i+1, m
			do j=1, N
				if (mask(j).eq.1) then
					aux(i,j)=&
					aux(i,j)+hh*G(k,j)
				else
					aux(i,:)=aux(i,:)-&
					hh**2*i*H(k,:)
				end if
			end do
		end do
		aux(i,:)=aux(i,:)+&
		(hh*i-(T_1-T_0))*integ/(T_1-T_0)
		norm=norm+dot_product(aux(i,:), aux(i,:))
	end do
	tam=eta
	do 
		if (abs(tam)<1.e-05) exit
		state(2:m+1,:)=state(2:m+1,:)+tam*aux
		cost=0.
		do i=1, m
			midpoint=&
			0.5*(state(i,:)+state(i+1,:))
			deriv=(state(i+1,:)-state(i,:))/hh
			call integrand(T_0+hh*(i-0.5), &
			midpoint,deriv,y,F)
			cost=cost+hh*F
		end do
		if (cost-cost_prev<0.) exit
		state(2:m+1,:)=state(2:m+1,:)-tam*aux
		tam=(0.5*norm*tam**2)/&
		(cost-cost_prev+norm*tam)
	end do
	if (abs(tam)<prec) exit
end do
s
do i=1, m
	deriv=(state(i+1,:)-state(i,:))/hh
	midpoint=0.5*(state(i,:)+state(i+1,:))
	call constraints(T_0+hh*(i-0.5),midpoint,deriv,r)
	do k=1, nn
		y(i,k)=y(i,k)*exp(y(i, k)*r(k))
	end do
end do
if (y_0/maxval(y)<1.e-05) exit
end do
! ======== Optimal control algorithm //
close(12)
end subroutine optimalcontrol

