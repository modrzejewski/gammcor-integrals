module auto2e_KetTransform_5_3
use arithmetic
use math_constants
implicit none
contains

subroutine auto2e_KetTransform_5_3_0_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(3,0,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(1)*Ey(1)*(Ex(7)*R(1:6) - Ex(8)*R(2:7) + Ex(9)*R(3:8) - Ex(10)*R(4:9))
S(7:11) = Ez(1)*Ey(1)*(Ex(7)*R(10:14) - Ex(8)*R(11:15) + Ex(9)*R(12:16) - Ex(10) &
   *R(13:17))
S(12:15) = Ez(1)*Ey(1)*(Ex(7)*R(18:21) - Ex(8)*R(19:22) + Ex(9)*R(20:23) - Ex(10) &
   *R(21:24))
S(16:18) = Ez(1)*Ey(1)*(Ex(7)*R(25:27) - Ex(8)*R(26:28) + Ex(9)*R(27:29) - Ex(10) &
   *R(28:30))
S(19:20) = Ez(1)*Ey(1)*(Ex(7)*R(31:32) - Ex(8)*R(32:33) + Ex(9)*R(33:34) - Ex(10) &
   *R(34:35))
S(21) = Ez(1)*Ey(1)*(Ex(7)*R(36) - Ex(8)*R(37) + Ex(9)*R(38) - Ex(10)*R(39))
S(22:26) = Ez(1)*Ey(1)*(Ex(7)*R(46:50) - Ex(8)*R(47:51) + Ex(9)*R(48:52) - Ex(10) &
   *R(49:53))
S(27:30) = Ez(1)*Ey(1)*(Ex(7)*R(54:57) - Ex(8)*R(55:58) + Ex(9)*R(56:59) - Ex(10) &
   *R(57:60))
S(31:33) = Ez(1)*Ey(1)*(Ex(7)*R(61:63) - Ex(8)*R(62:64) + Ex(9)*R(63:65) - Ex(10) &
   *R(64:66))
S(34:35) = Ez(1)*Ey(1)*(Ex(7)*R(67:68) - Ex(8)*R(68:69) + Ex(9)*R(69:70) - Ex(10) &
   *R(70:71))
S(36) = Ez(1)*Ey(1)*(Ex(7)*R(72) - Ex(8)*R(73) + Ex(9)*R(74) - Ex(10)*R(75))
S(37:40) = Ez(1)*Ey(1)*(Ex(7)*R(82:85) - Ex(8)*R(83:86) + Ex(9)*R(84:87) - Ex(10) &
   *R(85:88))
S(41:43) = Ez(1)*Ey(1)*(Ex(7)*R(89:91) - Ex(8)*R(90:92) + Ex(9)*R(91:93) - Ex(10) &
   *R(92:94))
S(44:45) = Ez(1)*Ey(1)*(Ex(7)*R(95:96) - Ex(8)*R(96:97) + Ex(9)*R(97:98) - Ex(10) &
   *R(98:99))
S(46) = Ez(1)*Ey(1)*(Ex(7)*R(100) - Ex(8)*R(101) + Ex(9)*R(102) - Ex(10)*R(103))
S(47:49) = Ez(1)*Ey(1)*(Ex(7)*R(110:112) - Ex(8)*R(111:113) + Ex(9)*R(112:114) - Ex(10) &
   *R(113:115))
S(50:51) = Ez(1)*Ey(1)*(Ex(7)*R(116:117) - Ex(8)*R(117:118) + Ex(9)*R(118:119) - Ex(10) &
   *R(119:120))
S(52) = Ez(1)*Ey(1)*(Ex(7)*R(121) - Ex(8)*R(122) + Ex(9)*R(123) - Ex(10)*R(124))
S(53:54) = Ez(1)*Ey(1)*(Ex(7)*R(131:132) - Ex(8)*R(132:133) + Ex(9)*R(133:134) - Ex(10) &
   *R(134:135))
S(55) = Ez(1)*Ey(1)*(Ex(7)*R(136) - Ex(8)*R(137) + Ex(9)*R(138) - Ex(10)*R(139))
S(56) = Ez(1)*Ey(1)*(Ex(7)*R(146) - Ex(8)*R(147) + Ex(9)*R(148) - Ex(10)*R(149))
end subroutine auto2e_KetTransform_5_3_0_0

subroutine auto2e_KetTransform_5_2_1_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(2,1,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(1)*Ey(2)*(Ex(4)*R(1:6) - Ex(5)*R(2:7) + Ex(6)*R(3:8)) + Ez(1)*Ey(3)*( &
   -Ex(4)*R(10:15) + Ex(5)*R(11:16) - Ex(6)*R(12:17))
S(7:11) = Ez(1)*Ey(2)*(Ex(4)*R(10:14) - Ex(5)*R(11:15) + Ex(6)*R(12:16)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(18:22) + Ex(5)*R(19:23) - Ex(6)*R(20:24))
S(12:15) = Ez(1)*Ey(2)*(Ex(4)*R(18:21) - Ex(5)*R(19:22) + Ex(6)*R(20:23)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(25:28) + Ex(5)*R(26:29) - Ex(6)*R(27:30))
S(16:18) = Ez(1)*Ey(2)*(Ex(4)*R(25:27) - Ex(5)*R(26:28) + Ex(6)*R(27:29)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(31:33) + Ex(5)*R(32:34) - Ex(6)*R(33:35))
S(19:20) = Ez(1)*Ey(2)*(Ex(4)*R(31:32) - Ex(5)*R(32:33) + Ex(6)*R(33:34)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(36:37) + Ex(5)*R(37:38) - Ex(6)*R(38:39))
S(21) = Ez(1)*Ey(2)*(Ex(4)*R(36) - Ex(5)*R(37) + Ex(6)*R(38)) + Ez(1)*Ey(3)*(-Ex(4) &
   *R(40) + Ex(5)*R(41) - Ex(6)*R(42))
S(22:26) = Ez(1)*Ey(2)*(Ex(4)*R(46:50) - Ex(5)*R(47:51) + Ex(6)*R(48:52)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(54:58) + Ex(5)*R(55:59) - Ex(6)*R(56:60))
S(27:30) = Ez(1)*Ey(2)*(Ex(4)*R(54:57) - Ex(5)*R(55:58) + Ex(6)*R(56:59)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(61:64) + Ex(5)*R(62:65) - Ex(6)*R(63:66))
S(31:33) = Ez(1)*Ey(2)*(Ex(4)*R(61:63) - Ex(5)*R(62:64) + Ex(6)*R(63:65)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(67:69) + Ex(5)*R(68:70) - Ex(6)*R(69:71))
S(34:35) = Ez(1)*Ey(2)*(Ex(4)*R(67:68) - Ex(5)*R(68:69) + Ex(6)*R(69:70)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(72:73) + Ex(5)*R(73:74) - Ex(6)*R(74:75))
S(36) = Ez(1)*Ey(2)*(Ex(4)*R(72) - Ex(5)*R(73) + Ex(6)*R(74)) + Ez(1)*Ey(3)*(-Ex(4) &
   *R(76) + Ex(5)*R(77) - Ex(6)*R(78))
S(37:40) = Ez(1)*Ey(2)*(Ex(4)*R(82:85) - Ex(5)*R(83:86) + Ex(6)*R(84:87)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(89:92) + Ex(5)*R(90:93) - Ex(6)*R(91:94))
S(41:43) = Ez(1)*Ey(2)*(Ex(4)*R(89:91) - Ex(5)*R(90:92) + Ex(6)*R(91:93)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(95:97) + Ex(5)*R(96:98) - Ex(6)*R(97:99))
S(44:45) = Ez(1)*Ey(2)*(Ex(4)*R(95:96) - Ex(5)*R(96:97) + Ex(6)*R(97:98)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(100:101) + Ex(5)*R(101:102) - Ex(6)*R(102:103))
S(46) = Ez(1)*Ey(2)*(Ex(4)*R(100) - Ex(5)*R(101) + Ex(6)*R(102)) + Ez(1)*Ey(3)*( &
   -Ex(4)*R(104) + Ex(5)*R(105) - Ex(6)*R(106))
S(47:49) = Ez(1)*Ey(2)*(Ex(4)*R(110:112) - Ex(5)*R(111:113) + Ex(6)*R(112:114))  &
   + Ez(1)*Ey(3)*(-Ex(4)*R(116:118) + Ex(5)*R(117:119) - Ex(6)*R(118:120))
S(50:51) = Ez(1)*Ey(2)*(Ex(4)*R(116:117) - Ex(5)*R(117:118) + Ex(6)*R(118:119))  &
   + Ez(1)*Ey(3)*(-Ex(4)*R(121:122) + Ex(5)*R(122:123) - Ex(6)*R(123:124))
S(52) = Ez(1)*Ey(2)*(Ex(4)*R(121) - Ex(5)*R(122) + Ex(6)*R(123)) + Ez(1)*Ey(3)*( &
   -Ex(4)*R(125) + Ex(5)*R(126) - Ex(6)*R(127))
S(53:54) = Ez(1)*Ey(2)*(Ex(4)*R(131:132) - Ex(5)*R(132:133) + Ex(6)*R(133:134))  &
   + Ez(1)*Ey(3)*(-Ex(4)*R(136:137) + Ex(5)*R(137:138) - Ex(6)*R(138:139))
S(55) = Ez(1)*Ey(2)*(Ex(4)*R(136) - Ex(5)*R(137) + Ex(6)*R(138)) + Ez(1)*Ey(3)*( &
   -Ex(4)*R(140) + Ex(5)*R(141) - Ex(6)*R(142))
S(56) = Ez(1)*Ey(2)*(Ex(4)*R(146) - Ex(5)*R(147) + Ex(6)*R(148)) + Ez(1)*Ey(3)*( &
   -Ex(4)*R(150) + Ex(5)*R(151) - Ex(6)*R(152))
end subroutine auto2e_KetTransform_5_2_1_0

subroutine auto2e_KetTransform_5_2_0_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(2,0,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(2)*Ey(1)*(Ex(4)*R(1:6) - Ex(5)*R(2:7) + Ex(6)*R(3:8)) + Ez(3)*Ey(1)*( &
   -Ex(4)*R(46:51) + Ex(5)*R(47:52) - Ex(6)*R(48:53))
S(7:11) = Ez(2)*Ey(1)*(Ex(4)*R(10:14) - Ex(5)*R(11:15) + Ex(6)*R(12:16)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(54:58) + Ex(5)*R(55:59) - Ex(6)*R(56:60))
S(12:15) = Ez(2)*Ey(1)*(Ex(4)*R(18:21) - Ex(5)*R(19:22) + Ex(6)*R(20:23)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(61:64) + Ex(5)*R(62:65) - Ex(6)*R(63:66))
S(16:18) = Ez(2)*Ey(1)*(Ex(4)*R(25:27) - Ex(5)*R(26:28) + Ex(6)*R(27:29)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(67:69) + Ex(5)*R(68:70) - Ex(6)*R(69:71))
S(19:20) = Ez(2)*Ey(1)*(Ex(4)*R(31:32) - Ex(5)*R(32:33) + Ex(6)*R(33:34)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(72:73) + Ex(5)*R(73:74) - Ex(6)*R(74:75))
S(21) = Ez(2)*Ey(1)*(Ex(4)*R(36) - Ex(5)*R(37) + Ex(6)*R(38)) + Ez(3)*Ey(1)*(-Ex(4) &
   *R(76) + Ex(5)*R(77) - Ex(6)*R(78))
S(22:26) = Ez(2)*Ey(1)*(Ex(4)*R(46:50) - Ex(5)*R(47:51) + Ex(6)*R(48:52)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(82:86) + Ex(5)*R(83:87) - Ex(6)*R(84:88))
S(27:30) = Ez(2)*Ey(1)*(Ex(4)*R(54:57) - Ex(5)*R(55:58) + Ex(6)*R(56:59)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(89:92) + Ex(5)*R(90:93) - Ex(6)*R(91:94))
S(31:33) = Ez(2)*Ey(1)*(Ex(4)*R(61:63) - Ex(5)*R(62:64) + Ex(6)*R(63:65)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(95:97) + Ex(5)*R(96:98) - Ex(6)*R(97:99))
S(34:35) = Ez(2)*Ey(1)*(Ex(4)*R(67:68) - Ex(5)*R(68:69) + Ex(6)*R(69:70)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(100:101) + Ex(5)*R(101:102) - Ex(6)*R(102:103))
S(36) = Ez(2)*Ey(1)*(Ex(4)*R(72) - Ex(5)*R(73) + Ex(6)*R(74)) + Ez(3)*Ey(1)*(-Ex(4) &
   *R(104) + Ex(5)*R(105) - Ex(6)*R(106))
S(37:40) = Ez(2)*Ey(1)*(Ex(4)*R(82:85) - Ex(5)*R(83:86) + Ex(6)*R(84:87)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(110:113) + Ex(5)*R(111:114) - Ex(6)*R(112:115))
S(41:43) = Ez(2)*Ey(1)*(Ex(4)*R(89:91) - Ex(5)*R(90:92) + Ex(6)*R(91:93)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(116:118) + Ex(5)*R(117:119) - Ex(6)*R(118:120))
S(44:45) = Ez(2)*Ey(1)*(Ex(4)*R(95:96) - Ex(5)*R(96:97) + Ex(6)*R(97:98)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(121:122) + Ex(5)*R(122:123) - Ex(6)*R(123:124))
S(46) = Ez(2)*Ey(1)*(Ex(4)*R(100) - Ex(5)*R(101) + Ex(6)*R(102)) + Ez(3)*Ey(1)*( &
   -Ex(4)*R(125) + Ex(5)*R(126) - Ex(6)*R(127))
S(47:49) = Ez(2)*Ey(1)*(Ex(4)*R(110:112) - Ex(5)*R(111:113) + Ex(6)*R(112:114))  &
   + Ez(3)*Ey(1)*(-Ex(4)*R(131:133) + Ex(5)*R(132:134) - Ex(6)*R(133:135))
S(50:51) = Ez(2)*Ey(1)*(Ex(4)*R(116:117) - Ex(5)*R(117:118) + Ex(6)*R(118:119))  &
   + Ez(3)*Ey(1)*(-Ex(4)*R(136:137) + Ex(5)*R(137:138) - Ex(6)*R(138:139))
S(52) = Ez(2)*Ey(1)*(Ex(4)*R(121) - Ex(5)*R(122) + Ex(6)*R(123)) + Ez(3)*Ey(1)*( &
   -Ex(4)*R(140) + Ex(5)*R(141) - Ex(6)*R(142))
S(53:54) = Ez(2)*Ey(1)*(Ex(4)*R(131:132) - Ex(5)*R(132:133) + Ex(6)*R(133:134))  &
   + Ez(3)*Ey(1)*(-Ex(4)*R(146:147) + Ex(5)*R(147:148) - Ex(6)*R(148:149))
S(55) = Ez(2)*Ey(1)*(Ex(4)*R(136) - Ex(5)*R(137) + Ex(6)*R(138)) + Ez(3)*Ey(1)*( &
   -Ex(4)*R(150) + Ex(5)*R(151) - Ex(6)*R(152))
S(56) = Ez(2)*Ey(1)*(Ex(4)*R(146) - Ex(5)*R(147) + Ex(6)*R(148)) + Ez(3)*Ey(1)*( &
   -Ex(4)*R(156) + Ex(5)*R(157) - Ex(6)*R(158))
end subroutine auto2e_KetTransform_5_2_0_1

subroutine auto2e_KetTransform_5_1_2_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(1,2,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(1)*Ey(4)*(Ex(2)*R(1:6) - Ex(3)*R(2:7)) + Ey(5)*Ez(1)*(-Ex(2)*R(10:15)  &
   + Ex(3)*R(11:16)) + Ez(1)*Ey(6)*(Ex(2)*R(18:23) - Ex(3)*R(19:24))
S(7:11) = Ez(1)*Ey(4)*(Ex(2)*R(10:14) - Ex(3)*R(11:15)) + Ey(5)*Ez(1)*(-Ex(2)*R(18:22)  &
   + Ex(3)*R(19:23)) + Ez(1)*Ey(6)*(Ex(2)*R(25:29) - Ex(3)*R(26:30))
S(12:15) = Ez(1)*Ey(4)*(Ex(2)*R(18:21) - Ex(3)*R(19:22)) + Ey(5)*Ez(1)*(-Ex(2)*R(25:28)  &
   + Ex(3)*R(26:29)) + Ez(1)*Ey(6)*(Ex(2)*R(31:34) - Ex(3)*R(32:35))
S(16:18) = Ez(1)*Ey(4)*(Ex(2)*R(25:27) - Ex(3)*R(26:28)) + Ey(5)*Ez(1)*(-Ex(2)*R(31:33)  &
   + Ex(3)*R(32:34)) + Ez(1)*Ey(6)*(Ex(2)*R(36:38) - Ex(3)*R(37:39))
S(19:20) = Ez(1)*Ey(4)*(Ex(2)*R(31:32) - Ex(3)*R(32:33)) + Ey(5)*Ez(1)*(-Ex(2)*R(36:37)  &
   + Ex(3)*R(37:38)) + Ez(1)*Ey(6)*(Ex(2)*R(40:41) - Ex(3)*R(41:42))
S(21) = Ez(1)*Ey(4)*(Ex(2)*R(36) - Ex(3)*R(37)) + Ey(5)*Ez(1)*(-Ex(2)*R(40) + Ex(3) &
   *R(41)) + Ez(1)*Ey(6)*(Ex(2)*R(43) - Ex(3)*R(44))
S(22:26) = Ez(1)*Ey(4)*(Ex(2)*R(46:50) - Ex(3)*R(47:51)) + Ey(5)*Ez(1)*(-Ex(2)*R(54:58)  &
   + Ex(3)*R(55:59)) + Ez(1)*Ey(6)*(Ex(2)*R(61:65) - Ex(3)*R(62:66))
S(27:30) = Ez(1)*Ey(4)*(Ex(2)*R(54:57) - Ex(3)*R(55:58)) + Ey(5)*Ez(1)*(-Ex(2)*R(61:64)  &
   + Ex(3)*R(62:65)) + Ez(1)*Ey(6)*(Ex(2)*R(67:70) - Ex(3)*R(68:71))
S(31:33) = Ez(1)*Ey(4)*(Ex(2)*R(61:63) - Ex(3)*R(62:64)) + Ey(5)*Ez(1)*(-Ex(2)*R(67:69)  &
   + Ex(3)*R(68:70)) + Ez(1)*Ey(6)*(Ex(2)*R(72:74) - Ex(3)*R(73:75))
S(34:35) = Ez(1)*Ey(4)*(Ex(2)*R(67:68) - Ex(3)*R(68:69)) + Ey(5)*Ez(1)*(-Ex(2)*R(72:73)  &
   + Ex(3)*R(73:74)) + Ez(1)*Ey(6)*(Ex(2)*R(76:77) - Ex(3)*R(77:78))
S(36) = Ez(1)*Ey(4)*(Ex(2)*R(72) - Ex(3)*R(73)) + Ey(5)*Ez(1)*(-Ex(2)*R(76) + Ex(3) &
   *R(77)) + Ez(1)*Ey(6)*(Ex(2)*R(79) - Ex(3)*R(80))
S(37:40) = Ez(1)*Ey(4)*(Ex(2)*R(82:85) - Ex(3)*R(83:86)) + Ey(5)*Ez(1)*(-Ex(2)*R(89:92)  &
   + Ex(3)*R(90:93)) + Ez(1)*Ey(6)*(Ex(2)*R(95:98) - Ex(3)*R(96:99))
S(41:43) = Ez(1)*Ey(4)*(Ex(2)*R(89:91) - Ex(3)*R(90:92)) + Ey(5)*Ez(1)*(-Ex(2)*R(95:97)  &
   + Ex(3)*R(96:98)) + Ez(1)*Ey(6)*(Ex(2)*R(100:102) - Ex(3)*R(101:103))
S(44:45) = Ez(1)*Ey(4)*(Ex(2)*R(95:96) - Ex(3)*R(96:97)) + Ey(5)*Ez(1)*(-Ex(2)*R(100:101)  &
   + Ex(3)*R(101:102)) + Ez(1)*Ey(6)*(Ex(2)*R(104:105) - Ex(3)*R(105:106))
S(46) = Ez(1)*Ey(4)*(Ex(2)*R(100) - Ex(3)*R(101)) + Ey(5)*Ez(1)*(-Ex(2)*R(104) + Ex(3) &
   *R(105)) + Ez(1)*Ey(6)*(Ex(2)*R(107) - Ex(3)*R(108))
S(47:49) = Ez(1)*Ey(4)*(Ex(2)*R(110:112) - Ex(3)*R(111:113)) + Ey(5)*Ez(1)*(-Ex(2) &
   *R(116:118) + Ex(3)*R(117:119)) + Ez(1)*Ey(6)*(Ex(2)*R(121:123) - Ex(3)*R(122:124))
S(50:51) = Ez(1)*Ey(4)*(Ex(2)*R(116:117) - Ex(3)*R(117:118)) + Ey(5)*Ez(1)*(-Ex(2) &
   *R(121:122) + Ex(3)*R(122:123)) + Ez(1)*Ey(6)*(Ex(2)*R(125:126) - Ex(3)*R(126:127))
S(52) = Ez(1)*Ey(4)*(Ex(2)*R(121) - Ex(3)*R(122)) + Ey(5)*Ez(1)*(-Ex(2)*R(125) + Ex(3) &
   *R(126)) + Ez(1)*Ey(6)*(Ex(2)*R(128) - Ex(3)*R(129))
S(53:54) = Ez(1)*Ey(4)*(Ex(2)*R(131:132) - Ex(3)*R(132:133)) + Ey(5)*Ez(1)*(-Ex(2) &
   *R(136:137) + Ex(3)*R(137:138)) + Ez(1)*Ey(6)*(Ex(2)*R(140:141) - Ex(3)*R(141:142))
S(55) = Ez(1)*Ey(4)*(Ex(2)*R(136) - Ex(3)*R(137)) + Ey(5)*Ez(1)*(-Ex(2)*R(140) + Ex(3) &
   *R(141)) + Ez(1)*Ey(6)*(Ex(2)*R(143) - Ex(3)*R(144))
S(56) = Ez(1)*Ey(4)*(Ex(2)*R(146) - Ex(3)*R(147)) + Ey(5)*Ez(1)*(-Ex(2)*R(150) + Ex(3) &
   *R(151)) + Ez(1)*Ey(6)*(Ex(2)*R(153) - Ex(3)*R(154))
end subroutine auto2e_KetTransform_5_1_2_0

subroutine auto2e_KetTransform_5_1_1_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(1,1,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(2)*Ey(2)*(Ex(2)*R(1:6) - Ex(3)*R(2:7)) + Ez(2)*Ey(3)*(-Ex(2)*R(10:15)  &
   + Ex(3)*R(11:16)) + Ez(3)*Ey(2)*(-Ex(2)*R(46:51) + Ex(3)*R(47:52)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(54:59) - Ex(3)*R(55:60))
S(7:11) = Ez(2)*Ey(2)*(Ex(2)*R(10:14) - Ex(3)*R(11:15)) + Ez(2)*Ey(3)*(-Ex(2)*R(18:22)  &
   + Ex(3)*R(19:23)) + Ez(3)*Ey(2)*(-Ex(2)*R(54:58) + Ex(3)*R(55:59)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(61:65) - Ex(3)*R(62:66))
S(12:15) = Ez(2)*Ey(2)*(Ex(2)*R(18:21) - Ex(3)*R(19:22)) + Ez(2)*Ey(3)*(-Ex(2)*R(25:28)  &
   + Ex(3)*R(26:29)) + Ez(3)*Ey(2)*(-Ex(2)*R(61:64) + Ex(3)*R(62:65)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(67:70) - Ex(3)*R(68:71))
S(16:18) = Ez(2)*Ey(2)*(Ex(2)*R(25:27) - Ex(3)*R(26:28)) + Ez(2)*Ey(3)*(-Ex(2)*R(31:33)  &
   + Ex(3)*R(32:34)) + Ez(3)*Ey(2)*(-Ex(2)*R(67:69) + Ex(3)*R(68:70)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(72:74) - Ex(3)*R(73:75))
S(19:20) = Ez(2)*Ey(2)*(Ex(2)*R(31:32) - Ex(3)*R(32:33)) + Ez(2)*Ey(3)*(-Ex(2)*R(36:37)  &
   + Ex(3)*R(37:38)) + Ez(3)*Ey(2)*(-Ex(2)*R(72:73) + Ex(3)*R(73:74)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(76:77) - Ex(3)*R(77:78))
S(21) = Ez(2)*Ey(2)*(Ex(2)*R(36) - Ex(3)*R(37)) + Ez(2)*Ey(3)*(-Ex(2)*R(40) + Ex(3) &
   *R(41)) + Ez(3)*Ey(2)*(-Ex(2)*R(76) + Ex(3)*R(77)) + Ey(3)*Ez(3)*(Ex(2)*R(79)  &
   - Ex(3)*R(80))
S(22:26) = Ez(2)*Ey(2)*(Ex(2)*R(46:50) - Ex(3)*R(47:51)) + Ez(2)*Ey(3)*(-Ex(2)*R(54:58)  &
   + Ex(3)*R(55:59)) + Ez(3)*Ey(2)*(-Ex(2)*R(82:86) + Ex(3)*R(83:87)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(89:93) - Ex(3)*R(90:94))
S(27:30) = Ez(2)*Ey(2)*(Ex(2)*R(54:57) - Ex(3)*R(55:58)) + Ez(2)*Ey(3)*(-Ex(2)*R(61:64)  &
   + Ex(3)*R(62:65)) + Ez(3)*Ey(2)*(-Ex(2)*R(89:92) + Ex(3)*R(90:93)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(95:98) - Ex(3)*R(96:99))
S(31:33) = Ez(2)*Ey(2)*(Ex(2)*R(61:63) - Ex(3)*R(62:64)) + Ez(2)*Ey(3)*(-Ex(2)*R(67:69)  &
   + Ex(3)*R(68:70)) + Ez(3)*Ey(2)*(-Ex(2)*R(95:97) + Ex(3)*R(96:98)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(100:102) - Ex(3)*R(101:103))
S(34:35) = Ez(2)*Ey(2)*(Ex(2)*R(67:68) - Ex(3)*R(68:69)) + Ez(2)*Ey(3)*(-Ex(2)*R(72:73)  &
   + Ex(3)*R(73:74)) + Ez(3)*Ey(2)*(-Ex(2)*R(100:101) + Ex(3)*R(101:102)) + Ey(3) &
   *Ez(3)*(Ex(2)*R(104:105) - Ex(3)*R(105:106))
S(36) = Ez(2)*Ey(2)*(Ex(2)*R(72) - Ex(3)*R(73)) + Ez(2)*Ey(3)*(-Ex(2)*R(76) + Ex(3) &
   *R(77)) + Ez(3)*Ey(2)*(-Ex(2)*R(104) + Ex(3)*R(105)) + Ey(3)*Ez(3)*(Ex(2)*R(107)  &
   - Ex(3)*R(108))
S(37:40) = Ez(2)*Ey(2)*(Ex(2)*R(82:85) - Ex(3)*R(83:86)) + Ez(2)*Ey(3)*(-Ex(2)*R(89:92)  &
   + Ex(3)*R(90:93)) + Ez(3)*Ey(2)*(-Ex(2)*R(110:113) + Ex(3)*R(111:114)) + Ey(3) &
   *Ez(3)*(Ex(2)*R(116:119) - Ex(3)*R(117:120))
S(41:43) = Ez(2)*Ey(2)*(Ex(2)*R(89:91) - Ex(3)*R(90:92)) + Ez(2)*Ey(3)*(-Ex(2)*R(95:97)  &
   + Ex(3)*R(96:98)) + Ez(3)*Ey(2)*(-Ex(2)*R(116:118) + Ex(3)*R(117:119)) + Ey(3) &
   *Ez(3)*(Ex(2)*R(121:123) - Ex(3)*R(122:124))
S(44:45) = Ez(2)*Ey(2)*(Ex(2)*R(95:96) - Ex(3)*R(96:97)) + Ez(2)*Ey(3)*(-Ex(2)*R(100:101)  &
   + Ex(3)*R(101:102)) + Ez(3)*Ey(2)*(-Ex(2)*R(121:122) + Ex(3)*R(122:123)) + Ey(3) &
   *Ez(3)*(Ex(2)*R(125:126) - Ex(3)*R(126:127))
S(46) = Ez(2)*Ey(2)*(Ex(2)*R(100) - Ex(3)*R(101)) + Ez(2)*Ey(3)*(-Ex(2)*R(104) + Ex(3) &
   *R(105)) + Ez(3)*Ey(2)*(-Ex(2)*R(125) + Ex(3)*R(126)) + Ey(3)*Ez(3)*(Ex(2)*R(128)  &
   - Ex(3)*R(129))
S(47:49) = Ez(2)*Ey(2)*(Ex(2)*R(110:112) - Ex(3)*R(111:113)) + Ez(2)*Ey(3)*(-Ex(2) &
   *R(116:118) + Ex(3)*R(117:119)) + Ez(3)*Ey(2)*(-Ex(2)*R(131:133) + Ex(3)*R(132:134))  &
   + Ey(3)*Ez(3)*(Ex(2)*R(136:138) - Ex(3)*R(137:139))
S(50:51) = Ez(2)*Ey(2)*(Ex(2)*R(116:117) - Ex(3)*R(117:118)) + Ez(2)*Ey(3)*(-Ex(2) &
   *R(121:122) + Ex(3)*R(122:123)) + Ez(3)*Ey(2)*(-Ex(2)*R(136:137) + Ex(3)*R(137:138))  &
   + Ey(3)*Ez(3)*(Ex(2)*R(140:141) - Ex(3)*R(141:142))
S(52) = Ez(2)*Ey(2)*(Ex(2)*R(121) - Ex(3)*R(122)) + Ez(2)*Ey(3)*(-Ex(2)*R(125) + Ex(3) &
   *R(126)) + Ez(3)*Ey(2)*(-Ex(2)*R(140) + Ex(3)*R(141)) + Ey(3)*Ez(3)*(Ex(2)*R(143)  &
   - Ex(3)*R(144))
S(53:54) = Ez(2)*Ey(2)*(Ex(2)*R(131:132) - Ex(3)*R(132:133)) + Ez(2)*Ey(3)*(-Ex(2) &
   *R(136:137) + Ex(3)*R(137:138)) + Ez(3)*Ey(2)*(-Ex(2)*R(146:147) + Ex(3)*R(147:148))  &
   + Ey(3)*Ez(3)*(Ex(2)*R(150:151) - Ex(3)*R(151:152))
S(55) = Ez(2)*Ey(2)*(Ex(2)*R(136) - Ex(3)*R(137)) + Ez(2)*Ey(3)*(-Ex(2)*R(140) + Ex(3) &
   *R(141)) + Ez(3)*Ey(2)*(-Ex(2)*R(150) + Ex(3)*R(151)) + Ey(3)*Ez(3)*(Ex(2)*R(153)  &
   - Ex(3)*R(154))
S(56) = Ez(2)*Ey(2)*(Ex(2)*R(146) - Ex(3)*R(147)) + Ez(2)*Ey(3)*(-Ex(2)*R(150) + Ex(3) &
   *R(151)) + Ez(3)*Ey(2)*(-Ex(2)*R(156) + Ex(3)*R(157)) + Ey(3)*Ez(3)*(Ex(2)*R(159)  &
   - Ex(3)*R(160))
end subroutine auto2e_KetTransform_5_1_1_1

subroutine auto2e_KetTransform_5_1_0_2(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(1,0,2)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(4)*Ey(1)*(Ex(2)*R(1:6) - Ex(3)*R(2:7)) + Ez(5)*Ey(1)*(-Ex(2)*R(46:51)  &
   + Ex(3)*R(47:52)) + Ez(6)*Ey(1)*(Ex(2)*R(82:87) - Ex(3)*R(83:88))
S(7:11) = Ez(4)*Ey(1)*(Ex(2)*R(10:14) - Ex(3)*R(11:15)) + Ez(5)*Ey(1)*(-Ex(2)*R(54:58)  &
   + Ex(3)*R(55:59)) + Ez(6)*Ey(1)*(Ex(2)*R(89:93) - Ex(3)*R(90:94))
S(12:15) = Ez(4)*Ey(1)*(Ex(2)*R(18:21) - Ex(3)*R(19:22)) + Ez(5)*Ey(1)*(-Ex(2)*R(61:64)  &
   + Ex(3)*R(62:65)) + Ez(6)*Ey(1)*(Ex(2)*R(95:98) - Ex(3)*R(96:99))
S(16:18) = Ez(4)*Ey(1)*(Ex(2)*R(25:27) - Ex(3)*R(26:28)) + Ez(5)*Ey(1)*(-Ex(2)*R(67:69)  &
   + Ex(3)*R(68:70)) + Ez(6)*Ey(1)*(Ex(2)*R(100:102) - Ex(3)*R(101:103))
S(19:20) = Ez(4)*Ey(1)*(Ex(2)*R(31:32) - Ex(3)*R(32:33)) + Ez(5)*Ey(1)*(-Ex(2)*R(72:73)  &
   + Ex(3)*R(73:74)) + Ez(6)*Ey(1)*(Ex(2)*R(104:105) - Ex(3)*R(105:106))
S(21) = Ez(4)*Ey(1)*(Ex(2)*R(36) - Ex(3)*R(37)) + Ez(5)*Ey(1)*(-Ex(2)*R(76) + Ex(3) &
   *R(77)) + Ez(6)*Ey(1)*(Ex(2)*R(107) - Ex(3)*R(108))
S(22:26) = Ez(4)*Ey(1)*(Ex(2)*R(46:50) - Ex(3)*R(47:51)) + Ez(5)*Ey(1)*(-Ex(2)*R(82:86)  &
   + Ex(3)*R(83:87)) + Ez(6)*Ey(1)*(Ex(2)*R(110:114) - Ex(3)*R(111:115))
S(27:30) = Ez(4)*Ey(1)*(Ex(2)*R(54:57) - Ex(3)*R(55:58)) + Ez(5)*Ey(1)*(-Ex(2)*R(89:92)  &
   + Ex(3)*R(90:93)) + Ez(6)*Ey(1)*(Ex(2)*R(116:119) - Ex(3)*R(117:120))
S(31:33) = Ez(4)*Ey(1)*(Ex(2)*R(61:63) - Ex(3)*R(62:64)) + Ez(5)*Ey(1)*(-Ex(2)*R(95:97)  &
   + Ex(3)*R(96:98)) + Ez(6)*Ey(1)*(Ex(2)*R(121:123) - Ex(3)*R(122:124))
S(34:35) = Ez(4)*Ey(1)*(Ex(2)*R(67:68) - Ex(3)*R(68:69)) + Ez(5)*Ey(1)*(-Ex(2)*R(100:101)  &
   + Ex(3)*R(101:102)) + Ez(6)*Ey(1)*(Ex(2)*R(125:126) - Ex(3)*R(126:127))
S(36) = Ez(4)*Ey(1)*(Ex(2)*R(72) - Ex(3)*R(73)) + Ez(5)*Ey(1)*(-Ex(2)*R(104) + Ex(3) &
   *R(105)) + Ez(6)*Ey(1)*(Ex(2)*R(128) - Ex(3)*R(129))
S(37:40) = Ez(4)*Ey(1)*(Ex(2)*R(82:85) - Ex(3)*R(83:86)) + Ez(5)*Ey(1)*(-Ex(2)*R(110:113)  &
   + Ex(3)*R(111:114)) + Ez(6)*Ey(1)*(Ex(2)*R(131:134) - Ex(3)*R(132:135))
S(41:43) = Ez(4)*Ey(1)*(Ex(2)*R(89:91) - Ex(3)*R(90:92)) + Ez(5)*Ey(1)*(-Ex(2)*R(116:118)  &
   + Ex(3)*R(117:119)) + Ez(6)*Ey(1)*(Ex(2)*R(136:138) - Ex(3)*R(137:139))
S(44:45) = Ez(4)*Ey(1)*(Ex(2)*R(95:96) - Ex(3)*R(96:97)) + Ez(5)*Ey(1)*(-Ex(2)*R(121:122)  &
   + Ex(3)*R(122:123)) + Ez(6)*Ey(1)*(Ex(2)*R(140:141) - Ex(3)*R(141:142))
S(46) = Ez(4)*Ey(1)*(Ex(2)*R(100) - Ex(3)*R(101)) + Ez(5)*Ey(1)*(-Ex(2)*R(125) + Ex(3) &
   *R(126)) + Ez(6)*Ey(1)*(Ex(2)*R(143) - Ex(3)*R(144))
S(47:49) = Ez(4)*Ey(1)*(Ex(2)*R(110:112) - Ex(3)*R(111:113)) + Ez(5)*Ey(1)*(-Ex(2) &
   *R(131:133) + Ex(3)*R(132:134)) + Ez(6)*Ey(1)*(Ex(2)*R(146:148) - Ex(3)*R(147:149))
S(50:51) = Ez(4)*Ey(1)*(Ex(2)*R(116:117) - Ex(3)*R(117:118)) + Ez(5)*Ey(1)*(-Ex(2) &
   *R(136:137) + Ex(3)*R(137:138)) + Ez(6)*Ey(1)*(Ex(2)*R(150:151) - Ex(3)*R(151:152))
S(52) = Ez(4)*Ey(1)*(Ex(2)*R(121) - Ex(3)*R(122)) + Ez(5)*Ey(1)*(-Ex(2)*R(140) + Ex(3) &
   *R(141)) + Ez(6)*Ey(1)*(Ex(2)*R(153) - Ex(3)*R(154))
S(53:54) = Ez(4)*Ey(1)*(Ex(2)*R(131:132) - Ex(3)*R(132:133)) + Ez(5)*Ey(1)*(-Ex(2) &
   *R(146:147) + Ex(3)*R(147:148)) + Ez(6)*Ey(1)*(Ex(2)*R(156:157) - Ex(3)*R(157:158))
S(55) = Ez(4)*Ey(1)*(Ex(2)*R(136) - Ex(3)*R(137)) + Ez(5)*Ey(1)*(-Ex(2)*R(150) + Ex(3) &
   *R(151)) + Ez(6)*Ey(1)*(Ex(2)*R(159) - Ex(3)*R(160))
S(56) = Ez(4)*Ey(1)*(Ex(2)*R(146) - Ex(3)*R(147)) + Ez(5)*Ey(1)*(-Ex(2)*R(156) + Ex(3) &
   *R(157)) + Ez(6)*Ey(1)*(Ex(2)*R(162) - Ex(3)*R(163))
end subroutine auto2e_KetTransform_5_1_0_2

subroutine auto2e_KetTransform_5_0_3_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(0,3,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(1)*Ex(1)*(Ey(7)*R(1:6) - Ey(8)*R(10:15) + Ey(9)*R(18:23) - Ey(10)*R(25:30))
S(7:11) = Ez(1)*Ex(1)*(Ey(7)*R(10:14) - Ey(8)*R(18:22) + Ey(9)*R(25:29) - Ey(10) &
   *R(31:35))
S(12:15) = Ez(1)*Ex(1)*(Ey(7)*R(18:21) - Ey(8)*R(25:28) + Ey(9)*R(31:34) - Ey(10) &
   *R(36:39))
S(16:18) = Ez(1)*Ex(1)*(Ey(7)*R(25:27) - Ey(8)*R(31:33) + Ey(9)*R(36:38) - Ey(10) &
   *R(40:42))
S(19:20) = Ez(1)*Ex(1)*(Ey(7)*R(31:32) - Ey(8)*R(36:37) + Ey(9)*R(40:41) - Ey(10) &
   *R(43:44))
S(21) = Ez(1)*Ex(1)*(Ey(7)*R(36) - Ey(8)*R(40) + Ey(9)*R(43) - Ey(10)*R(45))
S(22:26) = Ez(1)*Ex(1)*(Ey(7)*R(46:50) - Ey(8)*R(54:58) + Ey(9)*R(61:65) - Ey(10) &
   *R(67:71))
S(27:30) = Ez(1)*Ex(1)*(Ey(7)*R(54:57) - Ey(8)*R(61:64) + Ey(9)*R(67:70) - Ey(10) &
   *R(72:75))
S(31:33) = Ez(1)*Ex(1)*(Ey(7)*R(61:63) - Ey(8)*R(67:69) + Ey(9)*R(72:74) - Ey(10) &
   *R(76:78))
S(34:35) = Ez(1)*Ex(1)*(Ey(7)*R(67:68) - Ey(8)*R(72:73) + Ey(9)*R(76:77) - Ey(10) &
   *R(79:80))
S(36) = Ez(1)*Ex(1)*(Ey(7)*R(72) - Ey(8)*R(76) + Ey(9)*R(79) - Ey(10)*R(81))
S(37:40) = Ez(1)*Ex(1)*(Ey(7)*R(82:85) - Ey(8)*R(89:92) + Ey(9)*R(95:98) - Ey(10) &
   *R(100:103))
S(41:43) = Ez(1)*Ex(1)*(Ey(7)*R(89:91) - Ey(8)*R(95:97) + Ey(9)*R(100:102) - Ey(10) &
   *R(104:106))
S(44:45) = Ez(1)*Ex(1)*(Ey(7)*R(95:96) - Ey(8)*R(100:101) + Ey(9)*R(104:105) - Ey(10) &
   *R(107:108))
S(46) = Ez(1)*Ex(1)*(Ey(7)*R(100) - Ey(8)*R(104) + Ey(9)*R(107) - Ey(10)*R(109))
S(47:49) = Ez(1)*Ex(1)*(Ey(7)*R(110:112) - Ey(8)*R(116:118) + Ey(9)*R(121:123) - Ey(10) &
   *R(125:127))
S(50:51) = Ez(1)*Ex(1)*(Ey(7)*R(116:117) - Ey(8)*R(121:122) + Ey(9)*R(125:126) - Ey(10) &
   *R(128:129))
S(52) = Ez(1)*Ex(1)*(Ey(7)*R(121) - Ey(8)*R(125) + Ey(9)*R(128) - Ey(10)*R(130))
S(53:54) = Ez(1)*Ex(1)*(Ey(7)*R(131:132) - Ey(8)*R(136:137) + Ey(9)*R(140:141) - Ey(10) &
   *R(143:144))
S(55) = Ez(1)*Ex(1)*(Ey(7)*R(136) - Ey(8)*R(140) + Ey(9)*R(143) - Ey(10)*R(145))
S(56) = Ez(1)*Ex(1)*(Ey(7)*R(146) - Ey(8)*R(150) + Ey(9)*R(153) - Ey(10)*R(155))
end subroutine auto2e_KetTransform_5_0_3_0

subroutine auto2e_KetTransform_5_0_2_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(0,2,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(2)*Ex(1)*(Ey(4)*R(1:6) - Ey(5)*R(10:15) + Ey(6)*R(18:23)) + Ez(3)*Ex(1) &
   *(-Ey(4)*R(46:51) + Ey(5)*R(54:59) - Ey(6)*R(61:66))
S(7:11) = Ez(2)*Ex(1)*(Ey(4)*R(10:14) - Ey(5)*R(18:22) + Ey(6)*R(25:29)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(54:58) + Ey(5)*R(61:65) - Ey(6)*R(67:71))
S(12:15) = Ez(2)*Ex(1)*(Ey(4)*R(18:21) - Ey(5)*R(25:28) + Ey(6)*R(31:34)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(61:64) + Ey(5)*R(67:70) - Ey(6)*R(72:75))
S(16:18) = Ez(2)*Ex(1)*(Ey(4)*R(25:27) - Ey(5)*R(31:33) + Ey(6)*R(36:38)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(67:69) + Ey(5)*R(72:74) - Ey(6)*R(76:78))
S(19:20) = Ez(2)*Ex(1)*(Ey(4)*R(31:32) - Ey(5)*R(36:37) + Ey(6)*R(40:41)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(72:73) + Ey(5)*R(76:77) - Ey(6)*R(79:80))
S(21) = Ez(2)*Ex(1)*(Ey(4)*R(36) - Ey(5)*R(40) + Ey(6)*R(43)) + Ez(3)*Ex(1)*(-Ey(4) &
   *R(76) + Ey(5)*R(79) - Ey(6)*R(81))
S(22:26) = Ez(2)*Ex(1)*(Ey(4)*R(46:50) - Ey(5)*R(54:58) + Ey(6)*R(61:65)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(82:86) + Ey(5)*R(89:93) - Ey(6)*R(95:99))
S(27:30) = Ez(2)*Ex(1)*(Ey(4)*R(54:57) - Ey(5)*R(61:64) + Ey(6)*R(67:70)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(89:92) + Ey(5)*R(95:98) - Ey(6)*R(100:103))
S(31:33) = Ez(2)*Ex(1)*(Ey(4)*R(61:63) - Ey(5)*R(67:69) + Ey(6)*R(72:74)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(95:97) + Ey(5)*R(100:102) - Ey(6)*R(104:106))
S(34:35) = Ez(2)*Ex(1)*(Ey(4)*R(67:68) - Ey(5)*R(72:73) + Ey(6)*R(76:77)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(100:101) + Ey(5)*R(104:105) - Ey(6)*R(107:108))
S(36) = Ez(2)*Ex(1)*(Ey(4)*R(72) - Ey(5)*R(76) + Ey(6)*R(79)) + Ez(3)*Ex(1)*(-Ey(4) &
   *R(104) + Ey(5)*R(107) - Ey(6)*R(109))
S(37:40) = Ez(2)*Ex(1)*(Ey(4)*R(82:85) - Ey(5)*R(89:92) + Ey(6)*R(95:98)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(110:113) + Ey(5)*R(116:119) - Ey(6)*R(121:124))
S(41:43) = Ez(2)*Ex(1)*(Ey(4)*R(89:91) - Ey(5)*R(95:97) + Ey(6)*R(100:102)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(116:118) + Ey(5)*R(121:123) - Ey(6)*R(125:127))
S(44:45) = Ez(2)*Ex(1)*(Ey(4)*R(95:96) - Ey(5)*R(100:101) + Ey(6)*R(104:105)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(121:122) + Ey(5)*R(125:126) - Ey(6)*R(128:129))
S(46) = Ez(2)*Ex(1)*(Ey(4)*R(100) - Ey(5)*R(104) + Ey(6)*R(107)) + Ez(3)*Ex(1)*( &
   -Ey(4)*R(125) + Ey(5)*R(128) - Ey(6)*R(130))
S(47:49) = Ez(2)*Ex(1)*(Ey(4)*R(110:112) - Ey(5)*R(116:118) + Ey(6)*R(121:123))  &
   + Ez(3)*Ex(1)*(-Ey(4)*R(131:133) + Ey(5)*R(136:138) - Ey(6)*R(140:142))
S(50:51) = Ez(2)*Ex(1)*(Ey(4)*R(116:117) - Ey(5)*R(121:122) + Ey(6)*R(125:126))  &
   + Ez(3)*Ex(1)*(-Ey(4)*R(136:137) + Ey(5)*R(140:141) - Ey(6)*R(143:144))
S(52) = Ez(2)*Ex(1)*(Ey(4)*R(121) - Ey(5)*R(125) + Ey(6)*R(128)) + Ez(3)*Ex(1)*( &
   -Ey(4)*R(140) + Ey(5)*R(143) - Ey(6)*R(145))
S(53:54) = Ez(2)*Ex(1)*(Ey(4)*R(131:132) - Ey(5)*R(136:137) + Ey(6)*R(140:141))  &
   + Ez(3)*Ex(1)*(-Ey(4)*R(146:147) + Ey(5)*R(150:151) - Ey(6)*R(153:154))
S(55) = Ez(2)*Ex(1)*(Ey(4)*R(136) - Ey(5)*R(140) + Ey(6)*R(143)) + Ez(3)*Ex(1)*( &
   -Ey(4)*R(150) + Ey(5)*R(153) - Ey(6)*R(155))
S(56) = Ez(2)*Ex(1)*(Ey(4)*R(146) - Ey(5)*R(150) + Ey(6)*R(153)) + Ez(3)*Ex(1)*( &
   -Ey(4)*R(156) + Ey(5)*R(159) - Ey(6)*R(161))
end subroutine auto2e_KetTransform_5_0_2_1

subroutine auto2e_KetTransform_5_0_1_2(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(0,1,2)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ez(4)*Ex(1)*(Ey(2)*R(1:6) - Ey(3)*R(10:15)) + Ez(5)*Ex(1)*(-Ey(2)*R(46:51)  &
   + Ey(3)*R(54:59)) + Ex(1)*Ez(6)*(Ey(2)*R(82:87) - Ey(3)*R(89:94))
S(7:11) = Ez(4)*Ex(1)*(Ey(2)*R(10:14) - Ey(3)*R(18:22)) + Ez(5)*Ex(1)*(-Ey(2)*R(54:58)  &
   + Ey(3)*R(61:65)) + Ex(1)*Ez(6)*(Ey(2)*R(89:93) - Ey(3)*R(95:99))
S(12:15) = Ez(4)*Ex(1)*(Ey(2)*R(18:21) - Ey(3)*R(25:28)) + Ez(5)*Ex(1)*(-Ey(2)*R(61:64)  &
   + Ey(3)*R(67:70)) + Ex(1)*Ez(6)*(Ey(2)*R(95:98) - Ey(3)*R(100:103))
S(16:18) = Ez(4)*Ex(1)*(Ey(2)*R(25:27) - Ey(3)*R(31:33)) + Ez(5)*Ex(1)*(-Ey(2)*R(67:69)  &
   + Ey(3)*R(72:74)) + Ex(1)*Ez(6)*(Ey(2)*R(100:102) - Ey(3)*R(104:106))
S(19:20) = Ez(4)*Ex(1)*(Ey(2)*R(31:32) - Ey(3)*R(36:37)) + Ez(5)*Ex(1)*(-Ey(2)*R(72:73)  &
   + Ey(3)*R(76:77)) + Ex(1)*Ez(6)*(Ey(2)*R(104:105) - Ey(3)*R(107:108))
S(21) = Ez(4)*Ex(1)*(Ey(2)*R(36) - Ey(3)*R(40)) + Ez(5)*Ex(1)*(-Ey(2)*R(76) + Ey(3) &
   *R(79)) + Ex(1)*Ez(6)*(Ey(2)*R(107) - Ey(3)*R(109))
S(22:26) = Ez(4)*Ex(1)*(Ey(2)*R(46:50) - Ey(3)*R(54:58)) + Ez(5)*Ex(1)*(-Ey(2)*R(82:86)  &
   + Ey(3)*R(89:93)) + Ex(1)*Ez(6)*(Ey(2)*R(110:114) - Ey(3)*R(116:120))
S(27:30) = Ez(4)*Ex(1)*(Ey(2)*R(54:57) - Ey(3)*R(61:64)) + Ez(5)*Ex(1)*(-Ey(2)*R(89:92)  &
   + Ey(3)*R(95:98)) + Ex(1)*Ez(6)*(Ey(2)*R(116:119) - Ey(3)*R(121:124))
S(31:33) = Ez(4)*Ex(1)*(Ey(2)*R(61:63) - Ey(3)*R(67:69)) + Ez(5)*Ex(1)*(-Ey(2)*R(95:97)  &
   + Ey(3)*R(100:102)) + Ex(1)*Ez(6)*(Ey(2)*R(121:123) - Ey(3)*R(125:127))
S(34:35) = Ez(4)*Ex(1)*(Ey(2)*R(67:68) - Ey(3)*R(72:73)) + Ez(5)*Ex(1)*(-Ey(2)*R(100:101)  &
   + Ey(3)*R(104:105)) + Ex(1)*Ez(6)*(Ey(2)*R(125:126) - Ey(3)*R(128:129))
S(36) = Ez(4)*Ex(1)*(Ey(2)*R(72) - Ey(3)*R(76)) + Ez(5)*Ex(1)*(-Ey(2)*R(104) + Ey(3) &
   *R(107)) + Ex(1)*Ez(6)*(Ey(2)*R(128) - Ey(3)*R(130))
S(37:40) = Ez(4)*Ex(1)*(Ey(2)*R(82:85) - Ey(3)*R(89:92)) + Ez(5)*Ex(1)*(-Ey(2)*R(110:113)  &
   + Ey(3)*R(116:119)) + Ex(1)*Ez(6)*(Ey(2)*R(131:134) - Ey(3)*R(136:139))
S(41:43) = Ez(4)*Ex(1)*(Ey(2)*R(89:91) - Ey(3)*R(95:97)) + Ez(5)*Ex(1)*(-Ey(2)*R(116:118)  &
   + Ey(3)*R(121:123)) + Ex(1)*Ez(6)*(Ey(2)*R(136:138) - Ey(3)*R(140:142))
S(44:45) = Ez(4)*Ex(1)*(Ey(2)*R(95:96) - Ey(3)*R(100:101)) + Ez(5)*Ex(1)*(-Ey(2) &
   *R(121:122) + Ey(3)*R(125:126)) + Ex(1)*Ez(6)*(Ey(2)*R(140:141) - Ey(3)*R(143:144))
S(46) = Ez(4)*Ex(1)*(Ey(2)*R(100) - Ey(3)*R(104)) + Ez(5)*Ex(1)*(-Ey(2)*R(125) + Ey(3) &
   *R(128)) + Ex(1)*Ez(6)*(Ey(2)*R(143) - Ey(3)*R(145))
S(47:49) = Ez(4)*Ex(1)*(Ey(2)*R(110:112) - Ey(3)*R(116:118)) + Ez(5)*Ex(1)*(-Ey(2) &
   *R(131:133) + Ey(3)*R(136:138)) + Ex(1)*Ez(6)*(Ey(2)*R(146:148) - Ey(3)*R(150:152))
S(50:51) = Ez(4)*Ex(1)*(Ey(2)*R(116:117) - Ey(3)*R(121:122)) + Ez(5)*Ex(1)*(-Ey(2) &
   *R(136:137) + Ey(3)*R(140:141)) + Ex(1)*Ez(6)*(Ey(2)*R(150:151) - Ey(3)*R(153:154))
S(52) = Ez(4)*Ex(1)*(Ey(2)*R(121) - Ey(3)*R(125)) + Ez(5)*Ex(1)*(-Ey(2)*R(140) + Ey(3) &
   *R(143)) + Ex(1)*Ez(6)*(Ey(2)*R(153) - Ey(3)*R(155))
S(53:54) = Ez(4)*Ex(1)*(Ey(2)*R(131:132) - Ey(3)*R(136:137)) + Ez(5)*Ex(1)*(-Ey(2) &
   *R(146:147) + Ey(3)*R(150:151)) + Ex(1)*Ez(6)*(Ey(2)*R(156:157) - Ey(3)*R(159:160))
S(55) = Ez(4)*Ex(1)*(Ey(2)*R(136) - Ey(3)*R(140)) + Ez(5)*Ex(1)*(-Ey(2)*R(150) + Ey(3) &
   *R(153)) + Ex(1)*Ez(6)*(Ey(2)*R(159) - Ey(3)*R(161))
S(56) = Ez(4)*Ex(1)*(Ey(2)*R(146) - Ey(3)*R(150)) + Ez(5)*Ex(1)*(-Ey(2)*R(156) + Ey(3) &
   *R(159)) + Ex(1)*Ez(6)*(Ey(2)*R(162) - Ey(3)*R(164))
end subroutine auto2e_KetTransform_5_0_1_2

subroutine auto2e_KetTransform_5_0_0_3(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=5
! Class=(LS|KS), L=5, K=(0,0,3)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:6) = Ex(1)*Ey(1)*(Ez(7)*R(1:6) - Ez(8)*R(46:51) + Ez(9)*R(82:87) - Ez(10)*R(110:115))
S(7:11) = Ex(1)*Ey(1)*(Ez(7)*R(10:14) - Ez(8)*R(54:58) + Ez(9)*R(89:93) - Ez(10) &
   *R(116:120))
S(12:15) = Ex(1)*Ey(1)*(Ez(7)*R(18:21) - Ez(8)*R(61:64) + Ez(9)*R(95:98) - Ez(10) &
   *R(121:124))
S(16:18) = Ex(1)*Ey(1)*(Ez(7)*R(25:27) - Ez(8)*R(67:69) + Ez(9)*R(100:102) - Ez(10) &
   *R(125:127))
S(19:20) = Ex(1)*Ey(1)*(Ez(7)*R(31:32) - Ez(8)*R(72:73) + Ez(9)*R(104:105) - Ez(10) &
   *R(128:129))
S(21) = Ex(1)*Ey(1)*(Ez(7)*R(36) - Ez(8)*R(76) + Ez(9)*R(107) - Ez(10)*R(130))
S(22:26) = Ex(1)*Ey(1)*(Ez(7)*R(46:50) - Ez(8)*R(82:86) + Ez(9)*R(110:114) - Ez(10) &
   *R(131:135))
S(27:30) = Ex(1)*Ey(1)*(Ez(7)*R(54:57) - Ez(8)*R(89:92) + Ez(9)*R(116:119) - Ez(10) &
   *R(136:139))
S(31:33) = Ex(1)*Ey(1)*(Ez(7)*R(61:63) - Ez(8)*R(95:97) + Ez(9)*R(121:123) - Ez(10) &
   *R(140:142))
S(34:35) = Ex(1)*Ey(1)*(Ez(7)*R(67:68) - Ez(8)*R(100:101) + Ez(9)*R(125:126) - Ez(10) &
   *R(143:144))
S(36) = Ex(1)*Ey(1)*(Ez(7)*R(72) - Ez(8)*R(104) + Ez(9)*R(128) - Ez(10)*R(145))
S(37:40) = Ex(1)*Ey(1)*(Ez(7)*R(82:85) - Ez(8)*R(110:113) + Ez(9)*R(131:134) - Ez(10) &
   *R(146:149))
S(41:43) = Ex(1)*Ey(1)*(Ez(7)*R(89:91) - Ez(8)*R(116:118) + Ez(9)*R(136:138) - Ez(10) &
   *R(150:152))
S(44:45) = Ex(1)*Ey(1)*(Ez(7)*R(95:96) - Ez(8)*R(121:122) + Ez(9)*R(140:141) - Ez(10) &
   *R(153:154))
S(46) = Ex(1)*Ey(1)*(Ez(7)*R(100) - Ez(8)*R(125) + Ez(9)*R(143) - Ez(10)*R(155))
S(47:49) = Ex(1)*Ey(1)*(Ez(7)*R(110:112) - Ez(8)*R(131:133) + Ez(9)*R(146:148) - Ez(10) &
   *R(156:158))
S(50:51) = Ex(1)*Ey(1)*(Ez(7)*R(116:117) - Ez(8)*R(136:137) + Ez(9)*R(150:151) - Ez(10) &
   *R(159:160))
S(52) = Ex(1)*Ey(1)*(Ez(7)*R(121) - Ez(8)*R(140) + Ez(9)*R(153) - Ez(10)*R(161))
S(53:54) = Ex(1)*Ey(1)*(Ez(7)*R(131:132) - Ez(8)*R(146:147) + Ez(9)*R(156:157) - Ez(10) &
   *R(162:163))
S(55) = Ex(1)*Ey(1)*(Ez(7)*R(136) - Ez(8)*R(150) + Ez(9)*R(159) - Ez(10)*R(164))
S(56) = Ex(1)*Ey(1)*(Ez(7)*R(146) - Ez(8)*R(156) + Ez(9)*R(162) - Ez(10)*R(165))
end subroutine auto2e_KetTransform_5_0_0_3
end module auto2e_KetTransform_5_3