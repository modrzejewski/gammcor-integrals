module auto2e_KetTransform_3_1
use arithmetic
use math_constants
implicit none
contains

subroutine auto2e_KetTransform_3_1_0_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(1,0,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(1)*Ey(1)*(Ex(2)*R(1:4) - Ex(3)*R(2:5))
S(5:7) = Ez(1)*Ey(1)*(Ex(2)*R(6:8) - Ex(3)*R(7:9))
S(8:9) = Ez(1)*Ey(1)*(Ex(2)*R(10:11) - Ex(3)*R(11:12))
S(10) = Ez(1)*Ey(1)*(Ex(2)*R(13) - Ex(3)*R(14))
S(11:13) = Ez(1)*Ey(1)*(Ex(2)*R(16:18) - Ex(3)*R(17:19))
S(14:15) = Ez(1)*Ey(1)*(Ex(2)*R(20:21) - Ex(3)*R(21:22))
S(16) = Ez(1)*Ey(1)*(Ex(2)*R(23) - Ex(3)*R(24))
S(17:18) = Ez(1)*Ey(1)*(Ex(2)*R(26:27) - Ex(3)*R(27:28))
S(19) = Ez(1)*Ey(1)*(Ex(2)*R(29) - Ex(3)*R(30))
S(20) = Ez(1)*Ey(1)*(Ex(2)*R(32) - Ex(3)*R(33))
end subroutine auto2e_KetTransform_3_1_0_0

subroutine auto2e_KetTransform_3_0_1_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(0,1,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(1)*Ex(1)*(Ey(2)*R(1:4) - Ey(3)*R(6:9))
S(5:7) = Ez(1)*Ex(1)*(Ey(2)*R(6:8) - Ey(3)*R(10:12))
S(8:9) = Ez(1)*Ex(1)*(Ey(2)*R(10:11) - Ey(3)*R(13:14))
S(10) = Ez(1)*Ex(1)*(Ey(2)*R(13) - Ey(3)*R(15))
S(11:13) = Ez(1)*Ex(1)*(Ey(2)*R(16:18) - Ey(3)*R(20:22))
S(14:15) = Ez(1)*Ex(1)*(Ey(2)*R(20:21) - Ey(3)*R(23:24))
S(16) = Ez(1)*Ex(1)*(Ey(2)*R(23) - Ey(3)*R(25))
S(17:18) = Ez(1)*Ex(1)*(Ey(2)*R(26:27) - Ey(3)*R(29:30))
S(19) = Ez(1)*Ex(1)*(Ey(2)*R(29) - Ey(3)*R(31))
S(20) = Ez(1)*Ex(1)*(Ey(2)*R(32) - Ey(3)*R(34))
end subroutine auto2e_KetTransform_3_0_1_0

subroutine auto2e_KetTransform_3_0_0_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(0,0,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ex(1)*Ey(1)*(Ez(2)*R(1:4) - Ez(3)*R(16:19))
S(5:7) = Ex(1)*Ey(1)*(Ez(2)*R(6:8) - Ez(3)*R(20:22))
S(8:9) = Ex(1)*Ey(1)*(Ez(2)*R(10:11) - Ez(3)*R(23:24))
S(10) = Ex(1)*Ey(1)*(Ez(2)*R(13) - Ez(3)*R(25))
S(11:13) = Ex(1)*Ey(1)*(Ez(2)*R(16:18) - Ez(3)*R(26:28))
S(14:15) = Ex(1)*Ey(1)*(Ez(2)*R(20:21) - Ez(3)*R(29:30))
S(16) = Ex(1)*Ey(1)*(Ez(2)*R(23) - Ez(3)*R(31))
S(17:18) = Ex(1)*Ey(1)*(Ez(2)*R(26:27) - Ez(3)*R(32:33))
S(19) = Ex(1)*Ey(1)*(Ez(2)*R(29) - Ez(3)*R(34))
S(20) = Ex(1)*Ey(1)*(Ez(2)*R(32) - Ez(3)*R(35))
end subroutine auto2e_KetTransform_3_0_0_1
end module auto2e_KetTransform_3_1