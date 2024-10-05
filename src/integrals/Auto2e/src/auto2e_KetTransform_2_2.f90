module auto2e_KetTransform_2_2
use arithmetic
use math_constants
implicit none
contains

subroutine auto2e_KetTransform_2_2_0_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=2
! Class=(LS|KS), L=2, K=(2,0,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:3) = Ez(1)*Ey(1)*(Ex(4)*R(1:3) - Ex(5)*R(2:4) + Ex(6)*R(3:5))
S(4:5) = Ez(1)*Ey(1)*(Ex(4)*R(6:7) - Ex(5)*R(7:8) + Ex(6)*R(8:9))
S(6) = Ez(1)*Ey(1)*(Ex(4)*R(10) - Ex(5)*R(11) + Ex(6)*R(12))
S(7:8) = Ez(1)*Ey(1)*(Ex(4)*R(16:17) - Ex(5)*R(17:18) + Ex(6)*R(18:19))
S(9) = Ez(1)*Ey(1)*(Ex(4)*R(20) - Ex(5)*R(21) + Ex(6)*R(22))
S(10) = Ez(1)*Ey(1)*(Ex(4)*R(26) - Ex(5)*R(27) + Ex(6)*R(28))
end subroutine auto2e_KetTransform_2_2_0_0

subroutine auto2e_KetTransform_2_1_1_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=2
! Class=(LS|KS), L=2, K=(1,1,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:3) = Ez(1)*Ey(2)*(Ex(2)*R(1:3) - Ex(3)*R(2:4)) + Ez(1)*Ey(3)*(-Ex(2)*R(6:8)  &
   + Ex(3)*R(7:9))
S(4:5) = Ez(1)*Ey(2)*(Ex(2)*R(6:7) - Ex(3)*R(7:8)) + Ez(1)*Ey(3)*(-Ex(2)*R(10:11)  &
   + Ex(3)*R(11:12))
S(6) = Ez(1)*Ey(2)*(Ex(2)*R(10) - Ex(3)*R(11)) + Ez(1)*Ey(3)*(-Ex(2)*R(13) + Ex(3) &
   *R(14))
S(7:8) = Ez(1)*Ey(2)*(Ex(2)*R(16:17) - Ex(3)*R(17:18)) + Ez(1)*Ey(3)*(-Ex(2)*R(20:21)  &
   + Ex(3)*R(21:22))
S(9) = Ez(1)*Ey(2)*(Ex(2)*R(20) - Ex(3)*R(21)) + Ez(1)*Ey(3)*(-Ex(2)*R(23) + Ex(3) &
   *R(24))
S(10) = Ez(1)*Ey(2)*(Ex(2)*R(26) - Ex(3)*R(27)) + Ez(1)*Ey(3)*(-Ex(2)*R(29) + Ex(3) &
   *R(30))
end subroutine auto2e_KetTransform_2_1_1_0

subroutine auto2e_KetTransform_2_1_0_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=2
! Class=(LS|KS), L=2, K=(1,0,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:3) = Ez(2)*Ey(1)*(Ex(2)*R(1:3) - Ex(3)*R(2:4)) + Ez(3)*Ey(1)*(-Ex(2)*R(16:18)  &
   + Ex(3)*R(17:19))
S(4:5) = Ez(2)*Ey(1)*(Ex(2)*R(6:7) - Ex(3)*R(7:8)) + Ez(3)*Ey(1)*(-Ex(2)*R(20:21)  &
   + Ex(3)*R(21:22))
S(6) = Ez(2)*Ey(1)*(Ex(2)*R(10) - Ex(3)*R(11)) + Ez(3)*Ey(1)*(-Ex(2)*R(23) + Ex(3) &
   *R(24))
S(7:8) = Ez(2)*Ey(1)*(Ex(2)*R(16:17) - Ex(3)*R(17:18)) + Ez(3)*Ey(1)*(-Ex(2)*R(26:27)  &
   + Ex(3)*R(27:28))
S(9) = Ez(2)*Ey(1)*(Ex(2)*R(20) - Ex(3)*R(21)) + Ez(3)*Ey(1)*(-Ex(2)*R(29) + Ex(3) &
   *R(30))
S(10) = Ez(2)*Ey(1)*(Ex(2)*R(26) - Ex(3)*R(27)) + Ez(3)*Ey(1)*(-Ex(2)*R(32) + Ex(3) &
   *R(33))
end subroutine auto2e_KetTransform_2_1_0_1

subroutine auto2e_KetTransform_2_0_2_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=2
! Class=(LS|KS), L=2, K=(0,2,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:3) = Ez(1)*Ex(1)*(Ey(4)*R(1:3) - Ey(5)*R(6:8) + Ey(6)*R(10:12))
S(4:5) = Ez(1)*Ex(1)*(Ey(4)*R(6:7) - Ey(5)*R(10:11) + Ey(6)*R(13:14))
S(6) = Ez(1)*Ex(1)*(Ey(4)*R(10) - Ey(5)*R(13) + Ey(6)*R(15))
S(7:8) = Ez(1)*Ex(1)*(Ey(4)*R(16:17) - Ey(5)*R(20:21) + Ey(6)*R(23:24))
S(9) = Ez(1)*Ex(1)*(Ey(4)*R(20) - Ey(5)*R(23) + Ey(6)*R(25))
S(10) = Ez(1)*Ex(1)*(Ey(4)*R(26) - Ey(5)*R(29) + Ey(6)*R(31))
end subroutine auto2e_KetTransform_2_0_2_0

subroutine auto2e_KetTransform_2_0_1_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=2
! Class=(LS|KS), L=2, K=(0,1,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:3) = Ez(2)*Ex(1)*(Ey(2)*R(1:3) - Ey(3)*R(6:8)) + Ez(3)*Ex(1)*(-Ey(2)*R(16:18)  &
   + Ey(3)*R(20:22))
S(4:5) = Ez(2)*Ex(1)*(Ey(2)*R(6:7) - Ey(3)*R(10:11)) + Ez(3)*Ex(1)*(-Ey(2)*R(20:21)  &
   + Ey(3)*R(23:24))
S(6) = Ez(2)*Ex(1)*(Ey(2)*R(10) - Ey(3)*R(13)) + Ez(3)*Ex(1)*(-Ey(2)*R(23) + Ey(3) &
   *R(25))
S(7:8) = Ez(2)*Ex(1)*(Ey(2)*R(16:17) - Ey(3)*R(20:21)) + Ez(3)*Ex(1)*(-Ey(2)*R(26:27)  &
   + Ey(3)*R(29:30))
S(9) = Ez(2)*Ex(1)*(Ey(2)*R(20) - Ey(3)*R(23)) + Ez(3)*Ex(1)*(-Ey(2)*R(29) + Ey(3) &
   *R(31))
S(10) = Ez(2)*Ex(1)*(Ey(2)*R(26) - Ey(3)*R(29)) + Ez(3)*Ex(1)*(-Ey(2)*R(32) + Ey(3) &
   *R(34))
end subroutine auto2e_KetTransform_2_0_1_1

subroutine auto2e_KetTransform_2_0_0_2(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=2
! Class=(LS|KS), L=2, K=(0,0,2)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:3) = Ex(1)*Ey(1)*(Ez(4)*R(1:3) - Ez(5)*R(16:18) + Ez(6)*R(26:28))
S(4:5) = Ex(1)*Ey(1)*(Ez(4)*R(6:7) - Ez(5)*R(20:21) + Ez(6)*R(29:30))
S(6) = Ex(1)*Ey(1)*(Ez(4)*R(10) - Ez(5)*R(23) + Ez(6)*R(31))
S(7:8) = Ex(1)*Ey(1)*(Ez(4)*R(16:17) - Ez(5)*R(26:27) + Ez(6)*R(32:33))
S(9) = Ex(1)*Ey(1)*(Ez(4)*R(20) - Ez(5)*R(29) + Ez(6)*R(34))
S(10) = Ex(1)*Ey(1)*(Ez(4)*R(26) - Ez(5)*R(32) + Ez(6)*R(35))
end subroutine auto2e_KetTransform_2_0_0_2
end module auto2e_KetTransform_2_2