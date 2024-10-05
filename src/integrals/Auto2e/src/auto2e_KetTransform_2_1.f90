module auto2e_KetTransform_2_1
use arithmetic
use math_constants
implicit none
contains

subroutine auto2e_KetTransform_2_1_0_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=2
! Class=(LS|KS), L=2, K=(1,0,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:3) = Ez(1)*Ey(1)*(Ex(2)*R(1:3) - Ex(3)*R(2:4))
S(4:5) = Ez(1)*Ey(1)*(Ex(2)*R(5:6) - Ex(3)*R(6:7))
S(6) = Ez(1)*Ey(1)*(Ex(2)*R(8) - Ex(3)*R(9))
S(7:8) = Ez(1)*Ey(1)*(Ex(2)*R(11:12) - Ex(3)*R(12:13))
S(9) = Ez(1)*Ey(1)*(Ex(2)*R(14) - Ex(3)*R(15))
S(10) = Ez(1)*Ey(1)*(Ex(2)*R(17) - Ex(3)*R(18))
end subroutine auto2e_KetTransform_2_1_0_0

subroutine auto2e_KetTransform_2_0_1_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=2
! Class=(LS|KS), L=2, K=(0,1,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:3) = Ez(1)*Ex(1)*(Ey(2)*R(1:3) - Ey(3)*R(5:7))
S(4:5) = Ez(1)*Ex(1)*(Ey(2)*R(5:6) - Ey(3)*R(8:9))
S(6) = Ez(1)*Ex(1)*(Ey(2)*R(8) - Ey(3)*R(10))
S(7:8) = Ez(1)*Ex(1)*(Ey(2)*R(11:12) - Ey(3)*R(14:15))
S(9) = Ez(1)*Ex(1)*(Ey(2)*R(14) - Ey(3)*R(16))
S(10) = Ez(1)*Ex(1)*(Ey(2)*R(17) - Ey(3)*R(19))
end subroutine auto2e_KetTransform_2_0_1_0

subroutine auto2e_KetTransform_2_0_0_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=2
! Class=(LS|KS), L=2, K=(0,0,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:3) = Ex(1)*Ey(1)*(Ez(2)*R(1:3) - Ez(3)*R(11:13))
S(4:5) = Ex(1)*Ey(1)*(Ez(2)*R(5:6) - Ez(3)*R(14:15))
S(6) = Ex(1)*Ey(1)*(Ez(2)*R(8) - Ez(3)*R(16))
S(7:8) = Ex(1)*Ey(1)*(Ez(2)*R(11:12) - Ez(3)*R(17:18))
S(9) = Ex(1)*Ey(1)*(Ez(2)*R(14) - Ez(3)*R(19))
S(10) = Ex(1)*Ey(1)*(Ez(2)*R(17) - Ez(3)*R(20))
end subroutine auto2e_KetTransform_2_0_0_1
end module auto2e_KetTransform_2_1