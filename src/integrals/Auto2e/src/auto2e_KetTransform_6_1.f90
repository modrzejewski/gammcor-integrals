module auto2e_KetTransform_6_1
use arithmetic
use math_constants
implicit none
contains

subroutine auto2e_KetTransform_6_1_0_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=6
! Class=(LS|KS), L=6, K=(1,0,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:7) = Ez(1)*Ey(1)*(Ex(2)*R(1:7) - Ex(3)*R(2:8))
S(8:13) = Ez(1)*Ey(1)*(Ex(2)*R(9:14) - Ex(3)*R(10:15))
S(14:18) = Ez(1)*Ey(1)*(Ex(2)*R(16:20) - Ex(3)*R(17:21))
S(19:22) = Ez(1)*Ey(1)*(Ex(2)*R(22:25) - Ex(3)*R(23:26))
S(23:25) = Ez(1)*Ey(1)*(Ex(2)*R(27:29) - Ex(3)*R(28:30))
S(26:27) = Ez(1)*Ey(1)*(Ex(2)*R(31:32) - Ex(3)*R(32:33))
S(28) = Ez(1)*Ey(1)*(Ex(2)*R(34) - Ex(3)*R(35))
S(29:34) = Ez(1)*Ey(1)*(Ex(2)*R(37:42) - Ex(3)*R(38:43))
S(35:39) = Ez(1)*Ey(1)*(Ex(2)*R(44:48) - Ex(3)*R(45:49))
S(40:43) = Ez(1)*Ey(1)*(Ex(2)*R(50:53) - Ex(3)*R(51:54))
S(44:46) = Ez(1)*Ey(1)*(Ex(2)*R(55:57) - Ex(3)*R(56:58))
S(47:48) = Ez(1)*Ey(1)*(Ex(2)*R(59:60) - Ex(3)*R(60:61))
S(49) = Ez(1)*Ey(1)*(Ex(2)*R(62) - Ex(3)*R(63))
S(50:54) = Ez(1)*Ey(1)*(Ex(2)*R(65:69) - Ex(3)*R(66:70))
S(55:58) = Ez(1)*Ey(1)*(Ex(2)*R(71:74) - Ex(3)*R(72:75))
S(59:61) = Ez(1)*Ey(1)*(Ex(2)*R(76:78) - Ex(3)*R(77:79))
S(62:63) = Ez(1)*Ey(1)*(Ex(2)*R(80:81) - Ex(3)*R(81:82))
S(64) = Ez(1)*Ey(1)*(Ex(2)*R(83) - Ex(3)*R(84))
S(65:68) = Ez(1)*Ey(1)*(Ex(2)*R(86:89) - Ex(3)*R(87:90))
S(69:71) = Ez(1)*Ey(1)*(Ex(2)*R(91:93) - Ex(3)*R(92:94))
S(72:73) = Ez(1)*Ey(1)*(Ex(2)*R(95:96) - Ex(3)*R(96:97))
S(74) = Ez(1)*Ey(1)*(Ex(2)*R(98) - Ex(3)*R(99))
S(75:77) = Ez(1)*Ey(1)*(Ex(2)*R(101:103) - Ex(3)*R(102:104))
S(78:79) = Ez(1)*Ey(1)*(Ex(2)*R(105:106) - Ex(3)*R(106:107))
S(80) = Ez(1)*Ey(1)*(Ex(2)*R(108) - Ex(3)*R(109))
S(81:82) = Ez(1)*Ey(1)*(Ex(2)*R(111:112) - Ex(3)*R(112:113))
S(83) = Ez(1)*Ey(1)*(Ex(2)*R(114) - Ex(3)*R(115))
S(84) = Ez(1)*Ey(1)*(Ex(2)*R(117) - Ex(3)*R(118))
end subroutine auto2e_KetTransform_6_1_0_0

subroutine auto2e_KetTransform_6_0_1_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=6
! Class=(LS|KS), L=6, K=(0,1,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:7) = Ez(1)*Ex(1)*(Ey(2)*R(1:7) - Ey(3)*R(9:15))
S(8:13) = Ez(1)*Ex(1)*(Ey(2)*R(9:14) - Ey(3)*R(16:21))
S(14:18) = Ez(1)*Ex(1)*(Ey(2)*R(16:20) - Ey(3)*R(22:26))
S(19:22) = Ez(1)*Ex(1)*(Ey(2)*R(22:25) - Ey(3)*R(27:30))
S(23:25) = Ez(1)*Ex(1)*(Ey(2)*R(27:29) - Ey(3)*R(31:33))
S(26:27) = Ez(1)*Ex(1)*(Ey(2)*R(31:32) - Ey(3)*R(34:35))
S(28) = Ez(1)*Ex(1)*(Ey(2)*R(34) - Ey(3)*R(36))
S(29:34) = Ez(1)*Ex(1)*(Ey(2)*R(37:42) - Ey(3)*R(44:49))
S(35:39) = Ez(1)*Ex(1)*(Ey(2)*R(44:48) - Ey(3)*R(50:54))
S(40:43) = Ez(1)*Ex(1)*(Ey(2)*R(50:53) - Ey(3)*R(55:58))
S(44:46) = Ez(1)*Ex(1)*(Ey(2)*R(55:57) - Ey(3)*R(59:61))
S(47:48) = Ez(1)*Ex(1)*(Ey(2)*R(59:60) - Ey(3)*R(62:63))
S(49) = Ez(1)*Ex(1)*(Ey(2)*R(62) - Ey(3)*R(64))
S(50:54) = Ez(1)*Ex(1)*(Ey(2)*R(65:69) - Ey(3)*R(71:75))
S(55:58) = Ez(1)*Ex(1)*(Ey(2)*R(71:74) - Ey(3)*R(76:79))
S(59:61) = Ez(1)*Ex(1)*(Ey(2)*R(76:78) - Ey(3)*R(80:82))
S(62:63) = Ez(1)*Ex(1)*(Ey(2)*R(80:81) - Ey(3)*R(83:84))
S(64) = Ez(1)*Ex(1)*(Ey(2)*R(83) - Ey(3)*R(85))
S(65:68) = Ez(1)*Ex(1)*(Ey(2)*R(86:89) - Ey(3)*R(91:94))
S(69:71) = Ez(1)*Ex(1)*(Ey(2)*R(91:93) - Ey(3)*R(95:97))
S(72:73) = Ez(1)*Ex(1)*(Ey(2)*R(95:96) - Ey(3)*R(98:99))
S(74) = Ez(1)*Ex(1)*(Ey(2)*R(98) - Ey(3)*R(100))
S(75:77) = Ez(1)*Ex(1)*(Ey(2)*R(101:103) - Ey(3)*R(105:107))
S(78:79) = Ez(1)*Ex(1)*(Ey(2)*R(105:106) - Ey(3)*R(108:109))
S(80) = Ez(1)*Ex(1)*(Ey(2)*R(108) - Ey(3)*R(110))
S(81:82) = Ez(1)*Ex(1)*(Ey(2)*R(111:112) - Ey(3)*R(114:115))
S(83) = Ez(1)*Ex(1)*(Ey(2)*R(114) - Ey(3)*R(116))
S(84) = Ez(1)*Ex(1)*(Ey(2)*R(117) - Ey(3)*R(119))
end subroutine auto2e_KetTransform_6_0_1_0

subroutine auto2e_KetTransform_6_0_0_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=6
! Class=(LS|KS), L=6, K=(0,0,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:7) = Ex(1)*Ey(1)*(Ez(2)*R(1:7) - Ez(3)*R(37:43))
S(8:13) = Ex(1)*Ey(1)*(Ez(2)*R(9:14) - Ez(3)*R(44:49))
S(14:18) = Ex(1)*Ey(1)*(Ez(2)*R(16:20) - Ez(3)*R(50:54))
S(19:22) = Ex(1)*Ey(1)*(Ez(2)*R(22:25) - Ez(3)*R(55:58))
S(23:25) = Ex(1)*Ey(1)*(Ez(2)*R(27:29) - Ez(3)*R(59:61))
S(26:27) = Ex(1)*Ey(1)*(Ez(2)*R(31:32) - Ez(3)*R(62:63))
S(28) = Ex(1)*Ey(1)*(Ez(2)*R(34) - Ez(3)*R(64))
S(29:34) = Ex(1)*Ey(1)*(Ez(2)*R(37:42) - Ez(3)*R(65:70))
S(35:39) = Ex(1)*Ey(1)*(Ez(2)*R(44:48) - Ez(3)*R(71:75))
S(40:43) = Ex(1)*Ey(1)*(Ez(2)*R(50:53) - Ez(3)*R(76:79))
S(44:46) = Ex(1)*Ey(1)*(Ez(2)*R(55:57) - Ez(3)*R(80:82))
S(47:48) = Ex(1)*Ey(1)*(Ez(2)*R(59:60) - Ez(3)*R(83:84))
S(49) = Ex(1)*Ey(1)*(Ez(2)*R(62) - Ez(3)*R(85))
S(50:54) = Ex(1)*Ey(1)*(Ez(2)*R(65:69) - Ez(3)*R(86:90))
S(55:58) = Ex(1)*Ey(1)*(Ez(2)*R(71:74) - Ez(3)*R(91:94))
S(59:61) = Ex(1)*Ey(1)*(Ez(2)*R(76:78) - Ez(3)*R(95:97))
S(62:63) = Ex(1)*Ey(1)*(Ez(2)*R(80:81) - Ez(3)*R(98:99))
S(64) = Ex(1)*Ey(1)*(Ez(2)*R(83) - Ez(3)*R(100))
S(65:68) = Ex(1)*Ey(1)*(Ez(2)*R(86:89) - Ez(3)*R(101:104))
S(69:71) = Ex(1)*Ey(1)*(Ez(2)*R(91:93) - Ez(3)*R(105:107))
S(72:73) = Ex(1)*Ey(1)*(Ez(2)*R(95:96) - Ez(3)*R(108:109))
S(74) = Ex(1)*Ey(1)*(Ez(2)*R(98) - Ez(3)*R(110))
S(75:77) = Ex(1)*Ey(1)*(Ez(2)*R(101:103) - Ez(3)*R(111:113))
S(78:79) = Ex(1)*Ey(1)*(Ez(2)*R(105:106) - Ez(3)*R(114:115))
S(80) = Ex(1)*Ey(1)*(Ez(2)*R(108) - Ez(3)*R(116))
S(81:82) = Ex(1)*Ey(1)*(Ez(2)*R(111:112) - Ez(3)*R(117:118))
S(83) = Ex(1)*Ey(1)*(Ez(2)*R(114) - Ez(3)*R(119))
S(84) = Ex(1)*Ey(1)*(Ez(2)*R(117) - Ez(3)*R(120))
end subroutine auto2e_KetTransform_6_0_0_1
end module auto2e_KetTransform_6_1