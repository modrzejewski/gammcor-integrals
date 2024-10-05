module auto2e_KetTransform_5_1
use arithmetic
use math_constants
implicit none
contains

subroutine auto2e_KetTransform_5_1_0_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(1,0,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(1)*Ey(1)*(Ex(2)*R(1:6) - Ex(3)*R(2:7))
S(7:11) = Ez(1)*Ey(1)*(Ex(2)*R(8:12) - Ex(3)*R(9:13))
S(12:15) = Ez(1)*Ey(1)*(Ex(2)*R(14:17) - Ex(3)*R(15:18))
S(16:18) = Ez(1)*Ey(1)*(Ex(2)*R(19:21) - Ex(3)*R(20:22))
S(19:20) = Ez(1)*Ey(1)*(Ex(2)*R(23:24) - Ex(3)*R(24:25))
S(21) = Ez(1)*Ey(1)*(Ex(2)*R(26) - Ex(3)*R(27))
S(22:26) = Ez(1)*Ey(1)*(Ex(2)*R(29:33) - Ex(3)*R(30:34))
S(27:30) = Ez(1)*Ey(1)*(Ex(2)*R(35:38) - Ex(3)*R(36:39))
S(31:33) = Ez(1)*Ey(1)*(Ex(2)*R(40:42) - Ex(3)*R(41:43))
S(34:35) = Ez(1)*Ey(1)*(Ex(2)*R(44:45) - Ex(3)*R(45:46))
S(36) = Ez(1)*Ey(1)*(Ex(2)*R(47) - Ex(3)*R(48))
S(37:40) = Ez(1)*Ey(1)*(Ex(2)*R(50:53) - Ex(3)*R(51:54))
S(41:43) = Ez(1)*Ey(1)*(Ex(2)*R(55:57) - Ex(3)*R(56:58))
S(44:45) = Ez(1)*Ey(1)*(Ex(2)*R(59:60) - Ex(3)*R(60:61))
S(46) = Ez(1)*Ey(1)*(Ex(2)*R(62) - Ex(3)*R(63))
S(47:49) = Ez(1)*Ey(1)*(Ex(2)*R(65:67) - Ex(3)*R(66:68))
S(50:51) = Ez(1)*Ey(1)*(Ex(2)*R(69:70) - Ex(3)*R(70:71))
S(52) = Ez(1)*Ey(1)*(Ex(2)*R(72) - Ex(3)*R(73))
S(53:54) = Ez(1)*Ey(1)*(Ex(2)*R(75:76) - Ex(3)*R(76:77))
S(55) = Ez(1)*Ey(1)*(Ex(2)*R(78) - Ex(3)*R(79))
S(56) = Ez(1)*Ey(1)*(Ex(2)*R(81) - Ex(3)*R(82))
end subroutine auto2e_KetTransform_5_1_0_0

subroutine auto2e_KetTransform_5_0_1_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(0,1,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(1)*Ex(1)*(Ey(2)*R(1:6) - Ey(3)*R(8:13))
S(7:11) = Ez(1)*Ex(1)*(Ey(2)*R(8:12) - Ey(3)*R(14:18))
S(12:15) = Ez(1)*Ex(1)*(Ey(2)*R(14:17) - Ey(3)*R(19:22))
S(16:18) = Ez(1)*Ex(1)*(Ey(2)*R(19:21) - Ey(3)*R(23:25))
S(19:20) = Ez(1)*Ex(1)*(Ey(2)*R(23:24) - Ey(3)*R(26:27))
S(21) = Ez(1)*Ex(1)*(Ey(2)*R(26) - Ey(3)*R(28))
S(22:26) = Ez(1)*Ex(1)*(Ey(2)*R(29:33) - Ey(3)*R(35:39))
S(27:30) = Ez(1)*Ex(1)*(Ey(2)*R(35:38) - Ey(3)*R(40:43))
S(31:33) = Ez(1)*Ex(1)*(Ey(2)*R(40:42) - Ey(3)*R(44:46))
S(34:35) = Ez(1)*Ex(1)*(Ey(2)*R(44:45) - Ey(3)*R(47:48))
S(36) = Ez(1)*Ex(1)*(Ey(2)*R(47) - Ey(3)*R(49))
S(37:40) = Ez(1)*Ex(1)*(Ey(2)*R(50:53) - Ey(3)*R(55:58))
S(41:43) = Ez(1)*Ex(1)*(Ey(2)*R(55:57) - Ey(3)*R(59:61))
S(44:45) = Ez(1)*Ex(1)*(Ey(2)*R(59:60) - Ey(3)*R(62:63))
S(46) = Ez(1)*Ex(1)*(Ey(2)*R(62) - Ey(3)*R(64))
S(47:49) = Ez(1)*Ex(1)*(Ey(2)*R(65:67) - Ey(3)*R(69:71))
S(50:51) = Ez(1)*Ex(1)*(Ey(2)*R(69:70) - Ey(3)*R(72:73))
S(52) = Ez(1)*Ex(1)*(Ey(2)*R(72) - Ey(3)*R(74))
S(53:54) = Ez(1)*Ex(1)*(Ey(2)*R(75:76) - Ey(3)*R(78:79))
S(55) = Ez(1)*Ex(1)*(Ey(2)*R(78) - Ey(3)*R(80))
S(56) = Ez(1)*Ex(1)*(Ey(2)*R(81) - Ey(3)*R(83))
end subroutine auto2e_KetTransform_5_0_1_0

subroutine auto2e_KetTransform_5_0_0_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(0,0,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ex(1)*Ey(1)*(Ez(2)*R(1:6) - Ez(3)*R(29:34))
S(7:11) = Ex(1)*Ey(1)*(Ez(2)*R(8:12) - Ez(3)*R(35:39))
S(12:15) = Ex(1)*Ey(1)*(Ez(2)*R(14:17) - Ez(3)*R(40:43))
S(16:18) = Ex(1)*Ey(1)*(Ez(2)*R(19:21) - Ez(3)*R(44:46))
S(19:20) = Ex(1)*Ey(1)*(Ez(2)*R(23:24) - Ez(3)*R(47:48))
S(21) = Ex(1)*Ey(1)*(Ez(2)*R(26) - Ez(3)*R(49))
S(22:26) = Ex(1)*Ey(1)*(Ez(2)*R(29:33) - Ez(3)*R(50:54))
S(27:30) = Ex(1)*Ey(1)*(Ez(2)*R(35:38) - Ez(3)*R(55:58))
S(31:33) = Ex(1)*Ey(1)*(Ez(2)*R(40:42) - Ez(3)*R(59:61))
S(34:35) = Ex(1)*Ey(1)*(Ez(2)*R(44:45) - Ez(3)*R(62:63))
S(36) = Ex(1)*Ey(1)*(Ez(2)*R(47) - Ez(3)*R(64))
S(37:40) = Ex(1)*Ey(1)*(Ez(2)*R(50:53) - Ez(3)*R(65:68))
S(41:43) = Ex(1)*Ey(1)*(Ez(2)*R(55:57) - Ez(3)*R(69:71))
S(44:45) = Ex(1)*Ey(1)*(Ez(2)*R(59:60) - Ez(3)*R(72:73))
S(46) = Ex(1)*Ey(1)*(Ez(2)*R(62) - Ez(3)*R(74))
S(47:49) = Ex(1)*Ey(1)*(Ez(2)*R(65:67) - Ez(3)*R(75:77))
S(50:51) = Ex(1)*Ey(1)*(Ez(2)*R(69:70) - Ez(3)*R(78:79))
S(52) = Ex(1)*Ey(1)*(Ez(2)*R(72) - Ez(3)*R(80))
S(53:54) = Ex(1)*Ey(1)*(Ez(2)*R(75:76) - Ez(3)*R(81:82))
S(55) = Ex(1)*Ey(1)*(Ez(2)*R(78) - Ez(3)*R(83))
S(56) = Ex(1)*Ey(1)*(Ez(2)*R(81) - Ez(3)*R(84))
end subroutine auto2e_KetTransform_5_0_0_1
end module auto2e_KetTransform_5_1