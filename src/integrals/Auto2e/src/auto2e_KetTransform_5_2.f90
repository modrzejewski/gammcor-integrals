module auto2e_KetTransform_5_2
use arithmetic
use math_constants
implicit none
contains

subroutine auto2e_KetTransform_5_2_0_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(2,0,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(1)*Ey(1)*(Ex(4)*R(1:6) - Ex(5)*R(2:7) + Ex(6)*R(3:8))
S(7:11) = Ez(1)*Ey(1)*(Ex(4)*R(9:13) - Ex(5)*R(10:14) + Ex(6)*R(11:15))
S(12:15) = Ez(1)*Ey(1)*(Ex(4)*R(16:19) - Ex(5)*R(17:20) + Ex(6)*R(18:21))
S(16:18) = Ez(1)*Ey(1)*(Ex(4)*R(22:24) - Ex(5)*R(23:25) + Ex(6)*R(24:26))
S(19:20) = Ez(1)*Ey(1)*(Ex(4)*R(27:28) - Ex(5)*R(28:29) + Ex(6)*R(29:30))
S(21) = Ez(1)*Ey(1)*(Ex(4)*R(31) - Ex(5)*R(32) + Ex(6)*R(33))
S(22:26) = Ez(1)*Ey(1)*(Ex(4)*R(37:41) - Ex(5)*R(38:42) + Ex(6)*R(39:43))
S(27:30) = Ez(1)*Ey(1)*(Ex(4)*R(44:47) - Ex(5)*R(45:48) + Ex(6)*R(46:49))
S(31:33) = Ez(1)*Ey(1)*(Ex(4)*R(50:52) - Ex(5)*R(51:53) + Ex(6)*R(52:54))
S(34:35) = Ez(1)*Ey(1)*(Ex(4)*R(55:56) - Ex(5)*R(56:57) + Ex(6)*R(57:58))
S(36) = Ez(1)*Ey(1)*(Ex(4)*R(59) - Ex(5)*R(60) + Ex(6)*R(61))
S(37:40) = Ez(1)*Ey(1)*(Ex(4)*R(65:68) - Ex(5)*R(66:69) + Ex(6)*R(67:70))
S(41:43) = Ez(1)*Ey(1)*(Ex(4)*R(71:73) - Ex(5)*R(72:74) + Ex(6)*R(73:75))
S(44:45) = Ez(1)*Ey(1)*(Ex(4)*R(76:77) - Ex(5)*R(77:78) + Ex(6)*R(78:79))
S(46) = Ez(1)*Ey(1)*(Ex(4)*R(80) - Ex(5)*R(81) + Ex(6)*R(82))
S(47:49) = Ez(1)*Ey(1)*(Ex(4)*R(86:88) - Ex(5)*R(87:89) + Ex(6)*R(88:90))
S(50:51) = Ez(1)*Ey(1)*(Ex(4)*R(91:92) - Ex(5)*R(92:93) + Ex(6)*R(93:94))
S(52) = Ez(1)*Ey(1)*(Ex(4)*R(95) - Ex(5)*R(96) + Ex(6)*R(97))
S(53:54) = Ez(1)*Ey(1)*(Ex(4)*R(101:102) - Ex(5)*R(102:103) + Ex(6)*R(103:104))
S(55) = Ez(1)*Ey(1)*(Ex(4)*R(105) - Ex(5)*R(106) + Ex(6)*R(107))
S(56) = Ez(1)*Ey(1)*(Ex(4)*R(111) - Ex(5)*R(112) + Ex(6)*R(113))
end subroutine auto2e_KetTransform_5_2_0_0

subroutine auto2e_KetTransform_5_1_1_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(1,1,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(1)*Ey(2)*(Ex(2)*R(1:6) - Ex(3)*R(2:7)) + Ez(1)*Ey(3)*(-Ex(2)*R(9:14)  &
   + Ex(3)*R(10:15))
S(7:11) = Ez(1)*Ey(2)*(Ex(2)*R(9:13) - Ex(3)*R(10:14)) + Ez(1)*Ey(3)*(-Ex(2)*R(16:20)  &
   + Ex(3)*R(17:21))
S(12:15) = Ez(1)*Ey(2)*(Ex(2)*R(16:19) - Ex(3)*R(17:20)) + Ez(1)*Ey(3)*(-Ex(2)*R(22:25)  &
   + Ex(3)*R(23:26))
S(16:18) = Ez(1)*Ey(2)*(Ex(2)*R(22:24) - Ex(3)*R(23:25)) + Ez(1)*Ey(3)*(-Ex(2)*R(27:29)  &
   + Ex(3)*R(28:30))
S(19:20) = Ez(1)*Ey(2)*(Ex(2)*R(27:28) - Ex(3)*R(28:29)) + Ez(1)*Ey(3)*(-Ex(2)*R(31:32)  &
   + Ex(3)*R(32:33))
S(21) = Ez(1)*Ey(2)*(Ex(2)*R(31) - Ex(3)*R(32)) + Ez(1)*Ey(3)*(-Ex(2)*R(34) + Ex(3) &
   *R(35))
S(22:26) = Ez(1)*Ey(2)*(Ex(2)*R(37:41) - Ex(3)*R(38:42)) + Ez(1)*Ey(3)*(-Ex(2)*R(44:48)  &
   + Ex(3)*R(45:49))
S(27:30) = Ez(1)*Ey(2)*(Ex(2)*R(44:47) - Ex(3)*R(45:48)) + Ez(1)*Ey(3)*(-Ex(2)*R(50:53)  &
   + Ex(3)*R(51:54))
S(31:33) = Ez(1)*Ey(2)*(Ex(2)*R(50:52) - Ex(3)*R(51:53)) + Ez(1)*Ey(3)*(-Ex(2)*R(55:57)  &
   + Ex(3)*R(56:58))
S(34:35) = Ez(1)*Ey(2)*(Ex(2)*R(55:56) - Ex(3)*R(56:57)) + Ez(1)*Ey(3)*(-Ex(2)*R(59:60)  &
   + Ex(3)*R(60:61))
S(36) = Ez(1)*Ey(2)*(Ex(2)*R(59) - Ex(3)*R(60)) + Ez(1)*Ey(3)*(-Ex(2)*R(62) + Ex(3) &
   *R(63))
S(37:40) = Ez(1)*Ey(2)*(Ex(2)*R(65:68) - Ex(3)*R(66:69)) + Ez(1)*Ey(3)*(-Ex(2)*R(71:74)  &
   + Ex(3)*R(72:75))
S(41:43) = Ez(1)*Ey(2)*(Ex(2)*R(71:73) - Ex(3)*R(72:74)) + Ez(1)*Ey(3)*(-Ex(2)*R(76:78)  &
   + Ex(3)*R(77:79))
S(44:45) = Ez(1)*Ey(2)*(Ex(2)*R(76:77) - Ex(3)*R(77:78)) + Ez(1)*Ey(3)*(-Ex(2)*R(80:81)  &
   + Ex(3)*R(81:82))
S(46) = Ez(1)*Ey(2)*(Ex(2)*R(80) - Ex(3)*R(81)) + Ez(1)*Ey(3)*(-Ex(2)*R(83) + Ex(3) &
   *R(84))
S(47:49) = Ez(1)*Ey(2)*(Ex(2)*R(86:88) - Ex(3)*R(87:89)) + Ez(1)*Ey(3)*(-Ex(2)*R(91:93)  &
   + Ex(3)*R(92:94))
S(50:51) = Ez(1)*Ey(2)*(Ex(2)*R(91:92) - Ex(3)*R(92:93)) + Ez(1)*Ey(3)*(-Ex(2)*R(95:96)  &
   + Ex(3)*R(96:97))
S(52) = Ez(1)*Ey(2)*(Ex(2)*R(95) - Ex(3)*R(96)) + Ez(1)*Ey(3)*(-Ex(2)*R(98) + Ex(3) &
   *R(99))
S(53:54) = Ez(1)*Ey(2)*(Ex(2)*R(101:102) - Ex(3)*R(102:103)) + Ez(1)*Ey(3)*(-Ex(2) &
   *R(105:106) + Ex(3)*R(106:107))
S(55) = Ez(1)*Ey(2)*(Ex(2)*R(105) - Ex(3)*R(106)) + Ez(1)*Ey(3)*(-Ex(2)*R(108) + Ex(3) &
   *R(109))
S(56) = Ez(1)*Ey(2)*(Ex(2)*R(111) - Ex(3)*R(112)) + Ez(1)*Ey(3)*(-Ex(2)*R(114) + Ex(3) &
   *R(115))
end subroutine auto2e_KetTransform_5_1_1_0

subroutine auto2e_KetTransform_5_1_0_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(1,0,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(2)*Ey(1)*(Ex(2)*R(1:6) - Ex(3)*R(2:7)) + Ez(3)*Ey(1)*(-Ex(2)*R(37:42)  &
   + Ex(3)*R(38:43))
S(7:11) = Ez(2)*Ey(1)*(Ex(2)*R(9:13) - Ex(3)*R(10:14)) + Ez(3)*Ey(1)*(-Ex(2)*R(44:48)  &
   + Ex(3)*R(45:49))
S(12:15) = Ez(2)*Ey(1)*(Ex(2)*R(16:19) - Ex(3)*R(17:20)) + Ez(3)*Ey(1)*(-Ex(2)*R(50:53)  &
   + Ex(3)*R(51:54))
S(16:18) = Ez(2)*Ey(1)*(Ex(2)*R(22:24) - Ex(3)*R(23:25)) + Ez(3)*Ey(1)*(-Ex(2)*R(55:57)  &
   + Ex(3)*R(56:58))
S(19:20) = Ez(2)*Ey(1)*(Ex(2)*R(27:28) - Ex(3)*R(28:29)) + Ez(3)*Ey(1)*(-Ex(2)*R(59:60)  &
   + Ex(3)*R(60:61))
S(21) = Ez(2)*Ey(1)*(Ex(2)*R(31) - Ex(3)*R(32)) + Ez(3)*Ey(1)*(-Ex(2)*R(62) + Ex(3) &
   *R(63))
S(22:26) = Ez(2)*Ey(1)*(Ex(2)*R(37:41) - Ex(3)*R(38:42)) + Ez(3)*Ey(1)*(-Ex(2)*R(65:69)  &
   + Ex(3)*R(66:70))
S(27:30) = Ez(2)*Ey(1)*(Ex(2)*R(44:47) - Ex(3)*R(45:48)) + Ez(3)*Ey(1)*(-Ex(2)*R(71:74)  &
   + Ex(3)*R(72:75))
S(31:33) = Ez(2)*Ey(1)*(Ex(2)*R(50:52) - Ex(3)*R(51:53)) + Ez(3)*Ey(1)*(-Ex(2)*R(76:78)  &
   + Ex(3)*R(77:79))
S(34:35) = Ez(2)*Ey(1)*(Ex(2)*R(55:56) - Ex(3)*R(56:57)) + Ez(3)*Ey(1)*(-Ex(2)*R(80:81)  &
   + Ex(3)*R(81:82))
S(36) = Ez(2)*Ey(1)*(Ex(2)*R(59) - Ex(3)*R(60)) + Ez(3)*Ey(1)*(-Ex(2)*R(83) + Ex(3) &
   *R(84))
S(37:40) = Ez(2)*Ey(1)*(Ex(2)*R(65:68) - Ex(3)*R(66:69)) + Ez(3)*Ey(1)*(-Ex(2)*R(86:89)  &
   + Ex(3)*R(87:90))
S(41:43) = Ez(2)*Ey(1)*(Ex(2)*R(71:73) - Ex(3)*R(72:74)) + Ez(3)*Ey(1)*(-Ex(2)*R(91:93)  &
   + Ex(3)*R(92:94))
S(44:45) = Ez(2)*Ey(1)*(Ex(2)*R(76:77) - Ex(3)*R(77:78)) + Ez(3)*Ey(1)*(-Ex(2)*R(95:96)  &
   + Ex(3)*R(96:97))
S(46) = Ez(2)*Ey(1)*(Ex(2)*R(80) - Ex(3)*R(81)) + Ez(3)*Ey(1)*(-Ex(2)*R(98) + Ex(3) &
   *R(99))
S(47:49) = Ez(2)*Ey(1)*(Ex(2)*R(86:88) - Ex(3)*R(87:89)) + Ez(3)*Ey(1)*(-Ex(2)*R(101:103)  &
   + Ex(3)*R(102:104))
S(50:51) = Ez(2)*Ey(1)*(Ex(2)*R(91:92) - Ex(3)*R(92:93)) + Ez(3)*Ey(1)*(-Ex(2)*R(105:106)  &
   + Ex(3)*R(106:107))
S(52) = Ez(2)*Ey(1)*(Ex(2)*R(95) - Ex(3)*R(96)) + Ez(3)*Ey(1)*(-Ex(2)*R(108) + Ex(3) &
   *R(109))
S(53:54) = Ez(2)*Ey(1)*(Ex(2)*R(101:102) - Ex(3)*R(102:103)) + Ez(3)*Ey(1)*(-Ex(2) &
   *R(111:112) + Ex(3)*R(112:113))
S(55) = Ez(2)*Ey(1)*(Ex(2)*R(105) - Ex(3)*R(106)) + Ez(3)*Ey(1)*(-Ex(2)*R(114) + Ex(3) &
   *R(115))
S(56) = Ez(2)*Ey(1)*(Ex(2)*R(111) - Ex(3)*R(112)) + Ez(3)*Ey(1)*(-Ex(2)*R(117) + Ex(3) &
   *R(118))
end subroutine auto2e_KetTransform_5_1_0_1

subroutine auto2e_KetTransform_5_0_2_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(0,2,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(1)*Ex(1)*(Ey(4)*R(1:6) - Ey(5)*R(9:14) + Ey(6)*R(16:21))
S(7:11) = Ez(1)*Ex(1)*(Ey(4)*R(9:13) - Ey(5)*R(16:20) + Ey(6)*R(22:26))
S(12:15) = Ez(1)*Ex(1)*(Ey(4)*R(16:19) - Ey(5)*R(22:25) + Ey(6)*R(27:30))
S(16:18) = Ez(1)*Ex(1)*(Ey(4)*R(22:24) - Ey(5)*R(27:29) + Ey(6)*R(31:33))
S(19:20) = Ez(1)*Ex(1)*(Ey(4)*R(27:28) - Ey(5)*R(31:32) + Ey(6)*R(34:35))
S(21) = Ez(1)*Ex(1)*(Ey(4)*R(31) - Ey(5)*R(34) + Ey(6)*R(36))
S(22:26) = Ez(1)*Ex(1)*(Ey(4)*R(37:41) - Ey(5)*R(44:48) + Ey(6)*R(50:54))
S(27:30) = Ez(1)*Ex(1)*(Ey(4)*R(44:47) - Ey(5)*R(50:53) + Ey(6)*R(55:58))
S(31:33) = Ez(1)*Ex(1)*(Ey(4)*R(50:52) - Ey(5)*R(55:57) + Ey(6)*R(59:61))
S(34:35) = Ez(1)*Ex(1)*(Ey(4)*R(55:56) - Ey(5)*R(59:60) + Ey(6)*R(62:63))
S(36) = Ez(1)*Ex(1)*(Ey(4)*R(59) - Ey(5)*R(62) + Ey(6)*R(64))
S(37:40) = Ez(1)*Ex(1)*(Ey(4)*R(65:68) - Ey(5)*R(71:74) + Ey(6)*R(76:79))
S(41:43) = Ez(1)*Ex(1)*(Ey(4)*R(71:73) - Ey(5)*R(76:78) + Ey(6)*R(80:82))
S(44:45) = Ez(1)*Ex(1)*(Ey(4)*R(76:77) - Ey(5)*R(80:81) + Ey(6)*R(83:84))
S(46) = Ez(1)*Ex(1)*(Ey(4)*R(80) - Ey(5)*R(83) + Ey(6)*R(85))
S(47:49) = Ez(1)*Ex(1)*(Ey(4)*R(86:88) - Ey(5)*R(91:93) + Ey(6)*R(95:97))
S(50:51) = Ez(1)*Ex(1)*(Ey(4)*R(91:92) - Ey(5)*R(95:96) + Ey(6)*R(98:99))
S(52) = Ez(1)*Ex(1)*(Ey(4)*R(95) - Ey(5)*R(98) + Ey(6)*R(100))
S(53:54) = Ez(1)*Ex(1)*(Ey(4)*R(101:102) - Ey(5)*R(105:106) + Ey(6)*R(108:109))
S(55) = Ez(1)*Ex(1)*(Ey(4)*R(105) - Ey(5)*R(108) + Ey(6)*R(110))
S(56) = Ez(1)*Ex(1)*(Ey(4)*R(111) - Ey(5)*R(114) + Ey(6)*R(116))
end subroutine auto2e_KetTransform_5_0_2_0

subroutine auto2e_KetTransform_5_0_1_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(0,1,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(2)*Ex(1)*(Ey(2)*R(1:6) - Ey(3)*R(9:14)) + Ez(3)*Ex(1)*(-Ey(2)*R(37:42)  &
   + Ey(3)*R(44:49))
S(7:11) = Ez(2)*Ex(1)*(Ey(2)*R(9:13) - Ey(3)*R(16:20)) + Ez(3)*Ex(1)*(-Ey(2)*R(44:48)  &
   + Ey(3)*R(50:54))
S(12:15) = Ez(2)*Ex(1)*(Ey(2)*R(16:19) - Ey(3)*R(22:25)) + Ez(3)*Ex(1)*(-Ey(2)*R(50:53)  &
   + Ey(3)*R(55:58))
S(16:18) = Ez(2)*Ex(1)*(Ey(2)*R(22:24) - Ey(3)*R(27:29)) + Ez(3)*Ex(1)*(-Ey(2)*R(55:57)  &
   + Ey(3)*R(59:61))
S(19:20) = Ez(2)*Ex(1)*(Ey(2)*R(27:28) - Ey(3)*R(31:32)) + Ez(3)*Ex(1)*(-Ey(2)*R(59:60)  &
   + Ey(3)*R(62:63))
S(21) = Ez(2)*Ex(1)*(Ey(2)*R(31) - Ey(3)*R(34)) + Ez(3)*Ex(1)*(-Ey(2)*R(62) + Ey(3) &
   *R(64))
S(22:26) = Ez(2)*Ex(1)*(Ey(2)*R(37:41) - Ey(3)*R(44:48)) + Ez(3)*Ex(1)*(-Ey(2)*R(65:69)  &
   + Ey(3)*R(71:75))
S(27:30) = Ez(2)*Ex(1)*(Ey(2)*R(44:47) - Ey(3)*R(50:53)) + Ez(3)*Ex(1)*(-Ey(2)*R(71:74)  &
   + Ey(3)*R(76:79))
S(31:33) = Ez(2)*Ex(1)*(Ey(2)*R(50:52) - Ey(3)*R(55:57)) + Ez(3)*Ex(1)*(-Ey(2)*R(76:78)  &
   + Ey(3)*R(80:82))
S(34:35) = Ez(2)*Ex(1)*(Ey(2)*R(55:56) - Ey(3)*R(59:60)) + Ez(3)*Ex(1)*(-Ey(2)*R(80:81)  &
   + Ey(3)*R(83:84))
S(36) = Ez(2)*Ex(1)*(Ey(2)*R(59) - Ey(3)*R(62)) + Ez(3)*Ex(1)*(-Ey(2)*R(83) + Ey(3) &
   *R(85))
S(37:40) = Ez(2)*Ex(1)*(Ey(2)*R(65:68) - Ey(3)*R(71:74)) + Ez(3)*Ex(1)*(-Ey(2)*R(86:89)  &
   + Ey(3)*R(91:94))
S(41:43) = Ez(2)*Ex(1)*(Ey(2)*R(71:73) - Ey(3)*R(76:78)) + Ez(3)*Ex(1)*(-Ey(2)*R(91:93)  &
   + Ey(3)*R(95:97))
S(44:45) = Ez(2)*Ex(1)*(Ey(2)*R(76:77) - Ey(3)*R(80:81)) + Ez(3)*Ex(1)*(-Ey(2)*R(95:96)  &
   + Ey(3)*R(98:99))
S(46) = Ez(2)*Ex(1)*(Ey(2)*R(80) - Ey(3)*R(83)) + Ez(3)*Ex(1)*(-Ey(2)*R(98) + Ey(3) &
   *R(100))
S(47:49) = Ez(2)*Ex(1)*(Ey(2)*R(86:88) - Ey(3)*R(91:93)) + Ez(3)*Ex(1)*(-Ey(2)*R(101:103)  &
   + Ey(3)*R(105:107))
S(50:51) = Ez(2)*Ex(1)*(Ey(2)*R(91:92) - Ey(3)*R(95:96)) + Ez(3)*Ex(1)*(-Ey(2)*R(105:106)  &
   + Ey(3)*R(108:109))
S(52) = Ez(2)*Ex(1)*(Ey(2)*R(95) - Ey(3)*R(98)) + Ez(3)*Ex(1)*(-Ey(2)*R(108) + Ey(3) &
   *R(110))
S(53:54) = Ez(2)*Ex(1)*(Ey(2)*R(101:102) - Ey(3)*R(105:106)) + Ez(3)*Ex(1)*(-Ey(2) &
   *R(111:112) + Ey(3)*R(114:115))
S(55) = Ez(2)*Ex(1)*(Ey(2)*R(105) - Ey(3)*R(108)) + Ez(3)*Ex(1)*(-Ey(2)*R(114) + Ey(3) &
   *R(116))
S(56) = Ez(2)*Ex(1)*(Ey(2)*R(111) - Ey(3)*R(114)) + Ez(3)*Ex(1)*(-Ey(2)*R(117) + Ey(3) &
   *R(119))
end subroutine auto2e_KetTransform_5_0_1_1

subroutine auto2e_KetTransform_5_0_0_2(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(0,0,2)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ex(1)*Ey(1)*(Ez(4)*R(1:6) - Ez(5)*R(37:42) + Ez(6)*R(65:70))
S(7:11) = Ex(1)*Ey(1)*(Ez(4)*R(9:13) - Ez(5)*R(44:48) + Ez(6)*R(71:75))
S(12:15) = Ex(1)*Ey(1)*(Ez(4)*R(16:19) - Ez(5)*R(50:53) + Ez(6)*R(76:79))
S(16:18) = Ex(1)*Ey(1)*(Ez(4)*R(22:24) - Ez(5)*R(55:57) + Ez(6)*R(80:82))
S(19:20) = Ex(1)*Ey(1)*(Ez(4)*R(27:28) - Ez(5)*R(59:60) + Ez(6)*R(83:84))
S(21) = Ex(1)*Ey(1)*(Ez(4)*R(31) - Ez(5)*R(62) + Ez(6)*R(85))
S(22:26) = Ex(1)*Ey(1)*(Ez(4)*R(37:41) - Ez(5)*R(65:69) + Ez(6)*R(86:90))
S(27:30) = Ex(1)*Ey(1)*(Ez(4)*R(44:47) - Ez(5)*R(71:74) + Ez(6)*R(91:94))
S(31:33) = Ex(1)*Ey(1)*(Ez(4)*R(50:52) - Ez(5)*R(76:78) + Ez(6)*R(95:97))
S(34:35) = Ex(1)*Ey(1)*(Ez(4)*R(55:56) - Ez(5)*R(80:81) + Ez(6)*R(98:99))
S(36) = Ex(1)*Ey(1)*(Ez(4)*R(59) - Ez(5)*R(83) + Ez(6)*R(100))
S(37:40) = Ex(1)*Ey(1)*(Ez(4)*R(65:68) - Ez(5)*R(86:89) + Ez(6)*R(101:104))
S(41:43) = Ex(1)*Ey(1)*(Ez(4)*R(71:73) - Ez(5)*R(91:93) + Ez(6)*R(105:107))
S(44:45) = Ex(1)*Ey(1)*(Ez(4)*R(76:77) - Ez(5)*R(95:96) + Ez(6)*R(108:109))
S(46) = Ex(1)*Ey(1)*(Ez(4)*R(80) - Ez(5)*R(98) + Ez(6)*R(110))
S(47:49) = Ex(1)*Ey(1)*(Ez(4)*R(86:88) - Ez(5)*R(101:103) + Ez(6)*R(111:113))
S(50:51) = Ex(1)*Ey(1)*(Ez(4)*R(91:92) - Ez(5)*R(105:106) + Ez(6)*R(114:115))
S(52) = Ex(1)*Ey(1)*(Ez(4)*R(95) - Ez(5)*R(108) + Ez(6)*R(116))
S(53:54) = Ex(1)*Ey(1)*(Ez(4)*R(101:102) - Ez(5)*R(111:112) + Ez(6)*R(117:118))
S(55) = Ex(1)*Ey(1)*(Ez(4)*R(105) - Ez(5)*R(114) + Ez(6)*R(119))
S(56) = Ex(1)*Ey(1)*(Ez(4)*R(111) - Ez(5)*R(117) + Ez(6)*R(120))
end subroutine auto2e_KetTransform_5_0_0_2
end module auto2e_KetTransform_5_2