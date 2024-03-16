module auto2e_KetTransform_4_1
use arithmetic
use math_constants
implicit none
contains

subroutine auto2e_KetTransform_4_1_0_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(1,0,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(1)*Ey(1)*(Ex(2)*R(1:5) - Ex(3)*R(2:6))
S(6:9) = Ez(1)*Ey(1)*(Ex(2)*R(7:10) - Ex(3)*R(8:11))
S(10:12) = Ez(1)*Ey(1)*(Ex(2)*R(12:14) - Ex(3)*R(13:15))
S(13:14) = Ez(1)*Ey(1)*(Ex(2)*R(16:17) - Ex(3)*R(17:18))
S(15) = Ez(1)*Ey(1)*(Ex(2)*R(19) - Ex(3)*R(20))
S(16:19) = Ez(1)*Ey(1)*(Ex(2)*R(22:25) - Ex(3)*R(23:26))
S(20:22) = Ez(1)*Ey(1)*(Ex(2)*R(27:29) - Ex(3)*R(28:30))
S(23:24) = Ez(1)*Ey(1)*(Ex(2)*R(31:32) - Ex(3)*R(32:33))
S(25) = Ez(1)*Ey(1)*(Ex(2)*R(34) - Ex(3)*R(35))
S(26:28) = Ez(1)*Ey(1)*(Ex(2)*R(37:39) - Ex(3)*R(38:40))
S(29:30) = Ez(1)*Ey(1)*(Ex(2)*R(41:42) - Ex(3)*R(42:43))
S(31) = Ez(1)*Ey(1)*(Ex(2)*R(44) - Ex(3)*R(45))
S(32:33) = Ez(1)*Ey(1)*(Ex(2)*R(47:48) - Ex(3)*R(48:49))
S(34) = Ez(1)*Ey(1)*(Ex(2)*R(50) - Ex(3)*R(51))
S(35) = Ez(1)*Ey(1)*(Ex(2)*R(53) - Ex(3)*R(54))
end subroutine auto2e_KetTransform_4_1_0_0

subroutine auto2e_KetTransform_4_0_1_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(0,1,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(1)*Ex(1)*(Ey(2)*R(1:5) - Ey(3)*R(7:11))
S(6:9) = Ez(1)*Ex(1)*(Ey(2)*R(7:10) - Ey(3)*R(12:15))
S(10:12) = Ez(1)*Ex(1)*(Ey(2)*R(12:14) - Ey(3)*R(16:18))
S(13:14) = Ez(1)*Ex(1)*(Ey(2)*R(16:17) - Ey(3)*R(19:20))
S(15) = Ez(1)*Ex(1)*(Ey(2)*R(19) - Ey(3)*R(21))
S(16:19) = Ez(1)*Ex(1)*(Ey(2)*R(22:25) - Ey(3)*R(27:30))
S(20:22) = Ez(1)*Ex(1)*(Ey(2)*R(27:29) - Ey(3)*R(31:33))
S(23:24) = Ez(1)*Ex(1)*(Ey(2)*R(31:32) - Ey(3)*R(34:35))
S(25) = Ez(1)*Ex(1)*(Ey(2)*R(34) - Ey(3)*R(36))
S(26:28) = Ez(1)*Ex(1)*(Ey(2)*R(37:39) - Ey(3)*R(41:43))
S(29:30) = Ez(1)*Ex(1)*(Ey(2)*R(41:42) - Ey(3)*R(44:45))
S(31) = Ez(1)*Ex(1)*(Ey(2)*R(44) - Ey(3)*R(46))
S(32:33) = Ez(1)*Ex(1)*(Ey(2)*R(47:48) - Ey(3)*R(50:51))
S(34) = Ez(1)*Ex(1)*(Ey(2)*R(50) - Ey(3)*R(52))
S(35) = Ez(1)*Ex(1)*(Ey(2)*R(53) - Ey(3)*R(55))
end subroutine auto2e_KetTransform_4_0_1_0

subroutine auto2e_KetTransform_4_0_0_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(0,0,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ex(1)*Ey(1)*(Ez(2)*R(1:5) - Ez(3)*R(22:26))
S(6:9) = Ex(1)*Ey(1)*(Ez(2)*R(7:10) - Ez(3)*R(27:30))
S(10:12) = Ex(1)*Ey(1)*(Ez(2)*R(12:14) - Ez(3)*R(31:33))
S(13:14) = Ex(1)*Ey(1)*(Ez(2)*R(16:17) - Ez(3)*R(34:35))
S(15) = Ex(1)*Ey(1)*(Ez(2)*R(19) - Ez(3)*R(36))
S(16:19) = Ex(1)*Ey(1)*(Ez(2)*R(22:25) - Ez(3)*R(37:40))
S(20:22) = Ex(1)*Ey(1)*(Ez(2)*R(27:29) - Ez(3)*R(41:43))
S(23:24) = Ex(1)*Ey(1)*(Ez(2)*R(31:32) - Ez(3)*R(44:45))
S(25) = Ex(1)*Ey(1)*(Ez(2)*R(34) - Ez(3)*R(46))
S(26:28) = Ex(1)*Ey(1)*(Ez(2)*R(37:39) - Ez(3)*R(47:49))
S(29:30) = Ex(1)*Ey(1)*(Ez(2)*R(41:42) - Ez(3)*R(50:51))
S(31) = Ex(1)*Ey(1)*(Ez(2)*R(44) - Ez(3)*R(52))
S(32:33) = Ex(1)*Ey(1)*(Ez(2)*R(47:48) - Ez(3)*R(53:54))
S(34) = Ex(1)*Ey(1)*(Ez(2)*R(50) - Ez(3)*R(55))
S(35) = Ex(1)*Ey(1)*(Ez(2)*R(53) - Ez(3)*R(56))
end subroutine auto2e_KetTransform_4_0_0_1
end module auto2e_KetTransform_4_1