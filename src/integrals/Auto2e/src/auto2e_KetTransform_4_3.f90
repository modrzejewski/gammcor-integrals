module auto2e_KetTransform_4_3
use arithmetic
use math_constants
implicit none
contains

subroutine auto2e_KetTransform_4_3_0_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(3,0,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(1)*Ey(1)*(Ex(7)*R(1:5) - Ex(8)*R(2:6) + Ex(9)*R(3:7) - Ex(10)*R(4:8))
S(6:9) = Ez(1)*Ey(1)*(Ex(7)*R(9:12) - Ex(8)*R(10:13) + Ex(9)*R(11:14) - Ex(10)*R(12:15))
S(10:12) = Ez(1)*Ey(1)*(Ex(7)*R(16:18) - Ex(8)*R(17:19) + Ex(9)*R(18:20) - Ex(10) &
   *R(19:21))
S(13:14) = Ez(1)*Ey(1)*(Ex(7)*R(22:23) - Ex(8)*R(23:24) + Ex(9)*R(24:25) - Ex(10) &
   *R(25:26))
S(15) = Ez(1)*Ey(1)*(Ex(7)*R(27) - Ex(8)*R(28) + Ex(9)*R(29) - Ex(10)*R(30))
S(16:19) = Ez(1)*Ey(1)*(Ex(7)*R(37:40) - Ex(8)*R(38:41) + Ex(9)*R(39:42) - Ex(10) &
   *R(40:43))
S(20:22) = Ez(1)*Ey(1)*(Ex(7)*R(44:46) - Ex(8)*R(45:47) + Ex(9)*R(46:48) - Ex(10) &
   *R(47:49))
S(23:24) = Ez(1)*Ey(1)*(Ex(7)*R(50:51) - Ex(8)*R(51:52) + Ex(9)*R(52:53) - Ex(10) &
   *R(53:54))
S(25) = Ez(1)*Ey(1)*(Ex(7)*R(55) - Ex(8)*R(56) + Ex(9)*R(57) - Ex(10)*R(58))
S(26:28) = Ez(1)*Ey(1)*(Ex(7)*R(65:67) - Ex(8)*R(66:68) + Ex(9)*R(67:69) - Ex(10) &
   *R(68:70))
S(29:30) = Ez(1)*Ey(1)*(Ex(7)*R(71:72) - Ex(8)*R(72:73) + Ex(9)*R(73:74) - Ex(10) &
   *R(74:75))
S(31) = Ez(1)*Ey(1)*(Ex(7)*R(76) - Ex(8)*R(77) + Ex(9)*R(78) - Ex(10)*R(79))
S(32:33) = Ez(1)*Ey(1)*(Ex(7)*R(86:87) - Ex(8)*R(87:88) + Ex(9)*R(88:89) - Ex(10) &
   *R(89:90))
S(34) = Ez(1)*Ey(1)*(Ex(7)*R(91) - Ex(8)*R(92) + Ex(9)*R(93) - Ex(10)*R(94))
S(35) = Ez(1)*Ey(1)*(Ex(7)*R(101) - Ex(8)*R(102) + Ex(9)*R(103) - Ex(10)*R(104))
end subroutine auto2e_KetTransform_4_3_0_0

subroutine auto2e_KetTransform_4_2_1_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(2,1,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(1)*Ey(2)*(Ex(4)*R(1:5) - Ex(5)*R(2:6) + Ex(6)*R(3:7)) + Ez(1)*Ey(3)*( &
   -Ex(4)*R(9:13) + Ex(5)*R(10:14) - Ex(6)*R(11:15))
S(6:9) = Ez(1)*Ey(2)*(Ex(4)*R(9:12) - Ex(5)*R(10:13) + Ex(6)*R(11:14)) + Ez(1)*Ey(3) &
   *(-Ex(4)*R(16:19) + Ex(5)*R(17:20) - Ex(6)*R(18:21))
S(10:12) = Ez(1)*Ey(2)*(Ex(4)*R(16:18) - Ex(5)*R(17:19) + Ex(6)*R(18:20)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(22:24) + Ex(5)*R(23:25) - Ex(6)*R(24:26))
S(13:14) = Ez(1)*Ey(2)*(Ex(4)*R(22:23) - Ex(5)*R(23:24) + Ex(6)*R(24:25)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(27:28) + Ex(5)*R(28:29) - Ex(6)*R(29:30))
S(15) = Ez(1)*Ey(2)*(Ex(4)*R(27) - Ex(5)*R(28) + Ex(6)*R(29)) + Ez(1)*Ey(3)*(-Ex(4) &
   *R(31) + Ex(5)*R(32) - Ex(6)*R(33))
S(16:19) = Ez(1)*Ey(2)*(Ex(4)*R(37:40) - Ex(5)*R(38:41) + Ex(6)*R(39:42)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(44:47) + Ex(5)*R(45:48) - Ex(6)*R(46:49))
S(20:22) = Ez(1)*Ey(2)*(Ex(4)*R(44:46) - Ex(5)*R(45:47) + Ex(6)*R(46:48)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(50:52) + Ex(5)*R(51:53) - Ex(6)*R(52:54))
S(23:24) = Ez(1)*Ey(2)*(Ex(4)*R(50:51) - Ex(5)*R(51:52) + Ex(6)*R(52:53)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(55:56) + Ex(5)*R(56:57) - Ex(6)*R(57:58))
S(25) = Ez(1)*Ey(2)*(Ex(4)*R(55) - Ex(5)*R(56) + Ex(6)*R(57)) + Ez(1)*Ey(3)*(-Ex(4) &
   *R(59) + Ex(5)*R(60) - Ex(6)*R(61))
S(26:28) = Ez(1)*Ey(2)*(Ex(4)*R(65:67) - Ex(5)*R(66:68) + Ex(6)*R(67:69)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(71:73) + Ex(5)*R(72:74) - Ex(6)*R(73:75))
S(29:30) = Ez(1)*Ey(2)*(Ex(4)*R(71:72) - Ex(5)*R(72:73) + Ex(6)*R(73:74)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(76:77) + Ex(5)*R(77:78) - Ex(6)*R(78:79))
S(31) = Ez(1)*Ey(2)*(Ex(4)*R(76) - Ex(5)*R(77) + Ex(6)*R(78)) + Ez(1)*Ey(3)*(-Ex(4) &
   *R(80) + Ex(5)*R(81) - Ex(6)*R(82))
S(32:33) = Ez(1)*Ey(2)*(Ex(4)*R(86:87) - Ex(5)*R(87:88) + Ex(6)*R(88:89)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(91:92) + Ex(5)*R(92:93) - Ex(6)*R(93:94))
S(34) = Ez(1)*Ey(2)*(Ex(4)*R(91) - Ex(5)*R(92) + Ex(6)*R(93)) + Ez(1)*Ey(3)*(-Ex(4) &
   *R(95) + Ex(5)*R(96) - Ex(6)*R(97))
S(35) = Ez(1)*Ey(2)*(Ex(4)*R(101) - Ex(5)*R(102) + Ex(6)*R(103)) + Ez(1)*Ey(3)*( &
   -Ex(4)*R(105) + Ex(5)*R(106) - Ex(6)*R(107))
end subroutine auto2e_KetTransform_4_2_1_0

subroutine auto2e_KetTransform_4_2_0_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(2,0,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(2)*Ey(1)*(Ex(4)*R(1:5) - Ex(5)*R(2:6) + Ex(6)*R(3:7)) + Ez(3)*Ey(1)*( &
   -Ex(4)*R(37:41) + Ex(5)*R(38:42) - Ex(6)*R(39:43))
S(6:9) = Ez(2)*Ey(1)*(Ex(4)*R(9:12) - Ex(5)*R(10:13) + Ex(6)*R(11:14)) + Ez(3)*Ey(1) &
   *(-Ex(4)*R(44:47) + Ex(5)*R(45:48) - Ex(6)*R(46:49))
S(10:12) = Ez(2)*Ey(1)*(Ex(4)*R(16:18) - Ex(5)*R(17:19) + Ex(6)*R(18:20)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(50:52) + Ex(5)*R(51:53) - Ex(6)*R(52:54))
S(13:14) = Ez(2)*Ey(1)*(Ex(4)*R(22:23) - Ex(5)*R(23:24) + Ex(6)*R(24:25)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(55:56) + Ex(5)*R(56:57) - Ex(6)*R(57:58))
S(15) = Ez(2)*Ey(1)*(Ex(4)*R(27) - Ex(5)*R(28) + Ex(6)*R(29)) + Ez(3)*Ey(1)*(-Ex(4) &
   *R(59) + Ex(5)*R(60) - Ex(6)*R(61))
S(16:19) = Ez(2)*Ey(1)*(Ex(4)*R(37:40) - Ex(5)*R(38:41) + Ex(6)*R(39:42)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(65:68) + Ex(5)*R(66:69) - Ex(6)*R(67:70))
S(20:22) = Ez(2)*Ey(1)*(Ex(4)*R(44:46) - Ex(5)*R(45:47) + Ex(6)*R(46:48)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(71:73) + Ex(5)*R(72:74) - Ex(6)*R(73:75))
S(23:24) = Ez(2)*Ey(1)*(Ex(4)*R(50:51) - Ex(5)*R(51:52) + Ex(6)*R(52:53)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(76:77) + Ex(5)*R(77:78) - Ex(6)*R(78:79))
S(25) = Ez(2)*Ey(1)*(Ex(4)*R(55) - Ex(5)*R(56) + Ex(6)*R(57)) + Ez(3)*Ey(1)*(-Ex(4) &
   *R(80) + Ex(5)*R(81) - Ex(6)*R(82))
S(26:28) = Ez(2)*Ey(1)*(Ex(4)*R(65:67) - Ex(5)*R(66:68) + Ex(6)*R(67:69)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(86:88) + Ex(5)*R(87:89) - Ex(6)*R(88:90))
S(29:30) = Ez(2)*Ey(1)*(Ex(4)*R(71:72) - Ex(5)*R(72:73) + Ex(6)*R(73:74)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(91:92) + Ex(5)*R(92:93) - Ex(6)*R(93:94))
S(31) = Ez(2)*Ey(1)*(Ex(4)*R(76) - Ex(5)*R(77) + Ex(6)*R(78)) + Ez(3)*Ey(1)*(-Ex(4) &
   *R(95) + Ex(5)*R(96) - Ex(6)*R(97))
S(32:33) = Ez(2)*Ey(1)*(Ex(4)*R(86:87) - Ex(5)*R(87:88) + Ex(6)*R(88:89)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(101:102) + Ex(5)*R(102:103) - Ex(6)*R(103:104))
S(34) = Ez(2)*Ey(1)*(Ex(4)*R(91) - Ex(5)*R(92) + Ex(6)*R(93)) + Ez(3)*Ey(1)*(-Ex(4) &
   *R(105) + Ex(5)*R(106) - Ex(6)*R(107))
S(35) = Ez(2)*Ey(1)*(Ex(4)*R(101) - Ex(5)*R(102) + Ex(6)*R(103)) + Ez(3)*Ey(1)*( &
   -Ex(4)*R(111) + Ex(5)*R(112) - Ex(6)*R(113))
end subroutine auto2e_KetTransform_4_2_0_1

subroutine auto2e_KetTransform_4_1_2_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(1,2,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(1)*Ey(4)*(Ex(2)*R(1:5) - Ex(3)*R(2:6)) + Ey(5)*Ez(1)*(-Ex(2)*R(9:13)  &
   + Ex(3)*R(10:14)) + Ez(1)*Ey(6)*(Ex(2)*R(16:20) - Ex(3)*R(17:21))
S(6:9) = Ez(1)*Ey(4)*(Ex(2)*R(9:12) - Ex(3)*R(10:13)) + Ey(5)*Ez(1)*(-Ex(2)*R(16:19)  &
   + Ex(3)*R(17:20)) + Ez(1)*Ey(6)*(Ex(2)*R(22:25) - Ex(3)*R(23:26))
S(10:12) = Ez(1)*Ey(4)*(Ex(2)*R(16:18) - Ex(3)*R(17:19)) + Ey(5)*Ez(1)*(-Ex(2)*R(22:24)  &
   + Ex(3)*R(23:25)) + Ez(1)*Ey(6)*(Ex(2)*R(27:29) - Ex(3)*R(28:30))
S(13:14) = Ez(1)*Ey(4)*(Ex(2)*R(22:23) - Ex(3)*R(23:24)) + Ey(5)*Ez(1)*(-Ex(2)*R(27:28)  &
   + Ex(3)*R(28:29)) + Ez(1)*Ey(6)*(Ex(2)*R(31:32) - Ex(3)*R(32:33))
S(15) = Ez(1)*Ey(4)*(Ex(2)*R(27) - Ex(3)*R(28)) + Ey(5)*Ez(1)*(-Ex(2)*R(31) + Ex(3) &
   *R(32)) + Ez(1)*Ey(6)*(Ex(2)*R(34) - Ex(3)*R(35))
S(16:19) = Ez(1)*Ey(4)*(Ex(2)*R(37:40) - Ex(3)*R(38:41)) + Ey(5)*Ez(1)*(-Ex(2)*R(44:47)  &
   + Ex(3)*R(45:48)) + Ez(1)*Ey(6)*(Ex(2)*R(50:53) - Ex(3)*R(51:54))
S(20:22) = Ez(1)*Ey(4)*(Ex(2)*R(44:46) - Ex(3)*R(45:47)) + Ey(5)*Ez(1)*(-Ex(2)*R(50:52)  &
   + Ex(3)*R(51:53)) + Ez(1)*Ey(6)*(Ex(2)*R(55:57) - Ex(3)*R(56:58))
S(23:24) = Ez(1)*Ey(4)*(Ex(2)*R(50:51) - Ex(3)*R(51:52)) + Ey(5)*Ez(1)*(-Ex(2)*R(55:56)  &
   + Ex(3)*R(56:57)) + Ez(1)*Ey(6)*(Ex(2)*R(59:60) - Ex(3)*R(60:61))
S(25) = Ez(1)*Ey(4)*(Ex(2)*R(55) - Ex(3)*R(56)) + Ey(5)*Ez(1)*(-Ex(2)*R(59) + Ex(3) &
   *R(60)) + Ez(1)*Ey(6)*(Ex(2)*R(62) - Ex(3)*R(63))
S(26:28) = Ez(1)*Ey(4)*(Ex(2)*R(65:67) - Ex(3)*R(66:68)) + Ey(5)*Ez(1)*(-Ex(2)*R(71:73)  &
   + Ex(3)*R(72:74)) + Ez(1)*Ey(6)*(Ex(2)*R(76:78) - Ex(3)*R(77:79))
S(29:30) = Ez(1)*Ey(4)*(Ex(2)*R(71:72) - Ex(3)*R(72:73)) + Ey(5)*Ez(1)*(-Ex(2)*R(76:77)  &
   + Ex(3)*R(77:78)) + Ez(1)*Ey(6)*(Ex(2)*R(80:81) - Ex(3)*R(81:82))
S(31) = Ez(1)*Ey(4)*(Ex(2)*R(76) - Ex(3)*R(77)) + Ey(5)*Ez(1)*(-Ex(2)*R(80) + Ex(3) &
   *R(81)) + Ez(1)*Ey(6)*(Ex(2)*R(83) - Ex(3)*R(84))
S(32:33) = Ez(1)*Ey(4)*(Ex(2)*R(86:87) - Ex(3)*R(87:88)) + Ey(5)*Ez(1)*(-Ex(2)*R(91:92)  &
   + Ex(3)*R(92:93)) + Ez(1)*Ey(6)*(Ex(2)*R(95:96) - Ex(3)*R(96:97))
S(34) = Ez(1)*Ey(4)*(Ex(2)*R(91) - Ex(3)*R(92)) + Ey(5)*Ez(1)*(-Ex(2)*R(95) + Ex(3) &
   *R(96)) + Ez(1)*Ey(6)*(Ex(2)*R(98) - Ex(3)*R(99))
S(35) = Ez(1)*Ey(4)*(Ex(2)*R(101) - Ex(3)*R(102)) + Ey(5)*Ez(1)*(-Ex(2)*R(105) + Ex(3) &
   *R(106)) + Ez(1)*Ey(6)*(Ex(2)*R(108) - Ex(3)*R(109))
end subroutine auto2e_KetTransform_4_1_2_0

subroutine auto2e_KetTransform_4_1_1_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(1,1,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(2)*Ey(2)*(Ex(2)*R(1:5) - Ex(3)*R(2:6)) + Ez(2)*Ey(3)*(-Ex(2)*R(9:13)  &
   + Ex(3)*R(10:14)) + Ez(3)*Ey(2)*(-Ex(2)*R(37:41) + Ex(3)*R(38:42)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(44:48) - Ex(3)*R(45:49))
S(6:9) = Ez(2)*Ey(2)*(Ex(2)*R(9:12) - Ex(3)*R(10:13)) + Ez(2)*Ey(3)*(-Ex(2)*R(16:19)  &
   + Ex(3)*R(17:20)) + Ez(3)*Ey(2)*(-Ex(2)*R(44:47) + Ex(3)*R(45:48)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(50:53) - Ex(3)*R(51:54))
S(10:12) = Ez(2)*Ey(2)*(Ex(2)*R(16:18) - Ex(3)*R(17:19)) + Ez(2)*Ey(3)*(-Ex(2)*R(22:24)  &
   + Ex(3)*R(23:25)) + Ez(3)*Ey(2)*(-Ex(2)*R(50:52) + Ex(3)*R(51:53)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(55:57) - Ex(3)*R(56:58))
S(13:14) = Ez(2)*Ey(2)*(Ex(2)*R(22:23) - Ex(3)*R(23:24)) + Ez(2)*Ey(3)*(-Ex(2)*R(27:28)  &
   + Ex(3)*R(28:29)) + Ez(3)*Ey(2)*(-Ex(2)*R(55:56) + Ex(3)*R(56:57)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(59:60) - Ex(3)*R(60:61))
S(15) = Ez(2)*Ey(2)*(Ex(2)*R(27) - Ex(3)*R(28)) + Ez(2)*Ey(3)*(-Ex(2)*R(31) + Ex(3) &
   *R(32)) + Ez(3)*Ey(2)*(-Ex(2)*R(59) + Ex(3)*R(60)) + Ey(3)*Ez(3)*(Ex(2)*R(62)  &
   - Ex(3)*R(63))
S(16:19) = Ez(2)*Ey(2)*(Ex(2)*R(37:40) - Ex(3)*R(38:41)) + Ez(2)*Ey(3)*(-Ex(2)*R(44:47)  &
   + Ex(3)*R(45:48)) + Ez(3)*Ey(2)*(-Ex(2)*R(65:68) + Ex(3)*R(66:69)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(71:74) - Ex(3)*R(72:75))
S(20:22) = Ez(2)*Ey(2)*(Ex(2)*R(44:46) - Ex(3)*R(45:47)) + Ez(2)*Ey(3)*(-Ex(2)*R(50:52)  &
   + Ex(3)*R(51:53)) + Ez(3)*Ey(2)*(-Ex(2)*R(71:73) + Ex(3)*R(72:74)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(76:78) - Ex(3)*R(77:79))
S(23:24) = Ez(2)*Ey(2)*(Ex(2)*R(50:51) - Ex(3)*R(51:52)) + Ez(2)*Ey(3)*(-Ex(2)*R(55:56)  &
   + Ex(3)*R(56:57)) + Ez(3)*Ey(2)*(-Ex(2)*R(76:77) + Ex(3)*R(77:78)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(80:81) - Ex(3)*R(81:82))
S(25) = Ez(2)*Ey(2)*(Ex(2)*R(55) - Ex(3)*R(56)) + Ez(2)*Ey(3)*(-Ex(2)*R(59) + Ex(3) &
   *R(60)) + Ez(3)*Ey(2)*(-Ex(2)*R(80) + Ex(3)*R(81)) + Ey(3)*Ez(3)*(Ex(2)*R(83)  &
   - Ex(3)*R(84))
S(26:28) = Ez(2)*Ey(2)*(Ex(2)*R(65:67) - Ex(3)*R(66:68)) + Ez(2)*Ey(3)*(-Ex(2)*R(71:73)  &
   + Ex(3)*R(72:74)) + Ez(3)*Ey(2)*(-Ex(2)*R(86:88) + Ex(3)*R(87:89)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(91:93) - Ex(3)*R(92:94))
S(29:30) = Ez(2)*Ey(2)*(Ex(2)*R(71:72) - Ex(3)*R(72:73)) + Ez(2)*Ey(3)*(-Ex(2)*R(76:77)  &
   + Ex(3)*R(77:78)) + Ez(3)*Ey(2)*(-Ex(2)*R(91:92) + Ex(3)*R(92:93)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(95:96) - Ex(3)*R(96:97))
S(31) = Ez(2)*Ey(2)*(Ex(2)*R(76) - Ex(3)*R(77)) + Ez(2)*Ey(3)*(-Ex(2)*R(80) + Ex(3) &
   *R(81)) + Ez(3)*Ey(2)*(-Ex(2)*R(95) + Ex(3)*R(96)) + Ey(3)*Ez(3)*(Ex(2)*R(98)  &
   - Ex(3)*R(99))
S(32:33) = Ez(2)*Ey(2)*(Ex(2)*R(86:87) - Ex(3)*R(87:88)) + Ez(2)*Ey(3)*(-Ex(2)*R(91:92)  &
   + Ex(3)*R(92:93)) + Ez(3)*Ey(2)*(-Ex(2)*R(101:102) + Ex(3)*R(102:103)) + Ey(3) &
   *Ez(3)*(Ex(2)*R(105:106) - Ex(3)*R(106:107))
S(34) = Ez(2)*Ey(2)*(Ex(2)*R(91) - Ex(3)*R(92)) + Ez(2)*Ey(3)*(-Ex(2)*R(95) + Ex(3) &
   *R(96)) + Ez(3)*Ey(2)*(-Ex(2)*R(105) + Ex(3)*R(106)) + Ey(3)*Ez(3)*(Ex(2)*R(108)  &
   - Ex(3)*R(109))
S(35) = Ez(2)*Ey(2)*(Ex(2)*R(101) - Ex(3)*R(102)) + Ez(2)*Ey(3)*(-Ex(2)*R(105) + Ex(3) &
   *R(106)) + Ez(3)*Ey(2)*(-Ex(2)*R(111) + Ex(3)*R(112)) + Ey(3)*Ez(3)*(Ex(2)*R(114)  &
   - Ex(3)*R(115))
end subroutine auto2e_KetTransform_4_1_1_1

subroutine auto2e_KetTransform_4_1_0_2(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(1,0,2)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(4)*Ey(1)*(Ex(2)*R(1:5) - Ex(3)*R(2:6)) + Ez(5)*Ey(1)*(-Ex(2)*R(37:41)  &
   + Ex(3)*R(38:42)) + Ez(6)*Ey(1)*(Ex(2)*R(65:69) - Ex(3)*R(66:70))
S(6:9) = Ez(4)*Ey(1)*(Ex(2)*R(9:12) - Ex(3)*R(10:13)) + Ez(5)*Ey(1)*(-Ex(2)*R(44:47)  &
   + Ex(3)*R(45:48)) + Ez(6)*Ey(1)*(Ex(2)*R(71:74) - Ex(3)*R(72:75))
S(10:12) = Ez(4)*Ey(1)*(Ex(2)*R(16:18) - Ex(3)*R(17:19)) + Ez(5)*Ey(1)*(-Ex(2)*R(50:52)  &
   + Ex(3)*R(51:53)) + Ez(6)*Ey(1)*(Ex(2)*R(76:78) - Ex(3)*R(77:79))
S(13:14) = Ez(4)*Ey(1)*(Ex(2)*R(22:23) - Ex(3)*R(23:24)) + Ez(5)*Ey(1)*(-Ex(2)*R(55:56)  &
   + Ex(3)*R(56:57)) + Ez(6)*Ey(1)*(Ex(2)*R(80:81) - Ex(3)*R(81:82))
S(15) = Ez(4)*Ey(1)*(Ex(2)*R(27) - Ex(3)*R(28)) + Ez(5)*Ey(1)*(-Ex(2)*R(59) + Ex(3) &
   *R(60)) + Ez(6)*Ey(1)*(Ex(2)*R(83) - Ex(3)*R(84))
S(16:19) = Ez(4)*Ey(1)*(Ex(2)*R(37:40) - Ex(3)*R(38:41)) + Ez(5)*Ey(1)*(-Ex(2)*R(65:68)  &
   + Ex(3)*R(66:69)) + Ez(6)*Ey(1)*(Ex(2)*R(86:89) - Ex(3)*R(87:90))
S(20:22) = Ez(4)*Ey(1)*(Ex(2)*R(44:46) - Ex(3)*R(45:47)) + Ez(5)*Ey(1)*(-Ex(2)*R(71:73)  &
   + Ex(3)*R(72:74)) + Ez(6)*Ey(1)*(Ex(2)*R(91:93) - Ex(3)*R(92:94))
S(23:24) = Ez(4)*Ey(1)*(Ex(2)*R(50:51) - Ex(3)*R(51:52)) + Ez(5)*Ey(1)*(-Ex(2)*R(76:77)  &
   + Ex(3)*R(77:78)) + Ez(6)*Ey(1)*(Ex(2)*R(95:96) - Ex(3)*R(96:97))
S(25) = Ez(4)*Ey(1)*(Ex(2)*R(55) - Ex(3)*R(56)) + Ez(5)*Ey(1)*(-Ex(2)*R(80) + Ex(3) &
   *R(81)) + Ez(6)*Ey(1)*(Ex(2)*R(98) - Ex(3)*R(99))
S(26:28) = Ez(4)*Ey(1)*(Ex(2)*R(65:67) - Ex(3)*R(66:68)) + Ez(5)*Ey(1)*(-Ex(2)*R(86:88)  &
   + Ex(3)*R(87:89)) + Ez(6)*Ey(1)*(Ex(2)*R(101:103) - Ex(3)*R(102:104))
S(29:30) = Ez(4)*Ey(1)*(Ex(2)*R(71:72) - Ex(3)*R(72:73)) + Ez(5)*Ey(1)*(-Ex(2)*R(91:92)  &
   + Ex(3)*R(92:93)) + Ez(6)*Ey(1)*(Ex(2)*R(105:106) - Ex(3)*R(106:107))
S(31) = Ez(4)*Ey(1)*(Ex(2)*R(76) - Ex(3)*R(77)) + Ez(5)*Ey(1)*(-Ex(2)*R(95) + Ex(3) &
   *R(96)) + Ez(6)*Ey(1)*(Ex(2)*R(108) - Ex(3)*R(109))
S(32:33) = Ez(4)*Ey(1)*(Ex(2)*R(86:87) - Ex(3)*R(87:88)) + Ez(5)*Ey(1)*(-Ex(2)*R(101:102)  &
   + Ex(3)*R(102:103)) + Ez(6)*Ey(1)*(Ex(2)*R(111:112) - Ex(3)*R(112:113))
S(34) = Ez(4)*Ey(1)*(Ex(2)*R(91) - Ex(3)*R(92)) + Ez(5)*Ey(1)*(-Ex(2)*R(105) + Ex(3) &
   *R(106)) + Ez(6)*Ey(1)*(Ex(2)*R(114) - Ex(3)*R(115))
S(35) = Ez(4)*Ey(1)*(Ex(2)*R(101) - Ex(3)*R(102)) + Ez(5)*Ey(1)*(-Ex(2)*R(111) + Ex(3) &
   *R(112)) + Ez(6)*Ey(1)*(Ex(2)*R(117) - Ex(3)*R(118))
end subroutine auto2e_KetTransform_4_1_0_2

subroutine auto2e_KetTransform_4_0_3_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(0,3,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(1)*Ex(1)*(Ey(7)*R(1:5) - Ey(8)*R(9:13) + Ey(9)*R(16:20) - Ey(10)*R(22:26))
S(6:9) = Ez(1)*Ex(1)*(Ey(7)*R(9:12) - Ey(8)*R(16:19) + Ey(9)*R(22:25) - Ey(10)*R(27:30))
S(10:12) = Ez(1)*Ex(1)*(Ey(7)*R(16:18) - Ey(8)*R(22:24) + Ey(9)*R(27:29) - Ey(10) &
   *R(31:33))
S(13:14) = Ez(1)*Ex(1)*(Ey(7)*R(22:23) - Ey(8)*R(27:28) + Ey(9)*R(31:32) - Ey(10) &
   *R(34:35))
S(15) = Ez(1)*Ex(1)*(Ey(7)*R(27) - Ey(8)*R(31) + Ey(9)*R(34) - Ey(10)*R(36))
S(16:19) = Ez(1)*Ex(1)*(Ey(7)*R(37:40) - Ey(8)*R(44:47) + Ey(9)*R(50:53) - Ey(10) &
   *R(55:58))
S(20:22) = Ez(1)*Ex(1)*(Ey(7)*R(44:46) - Ey(8)*R(50:52) + Ey(9)*R(55:57) - Ey(10) &
   *R(59:61))
S(23:24) = Ez(1)*Ex(1)*(Ey(7)*R(50:51) - Ey(8)*R(55:56) + Ey(9)*R(59:60) - Ey(10) &
   *R(62:63))
S(25) = Ez(1)*Ex(1)*(Ey(7)*R(55) - Ey(8)*R(59) + Ey(9)*R(62) - Ey(10)*R(64))
S(26:28) = Ez(1)*Ex(1)*(Ey(7)*R(65:67) - Ey(8)*R(71:73) + Ey(9)*R(76:78) - Ey(10) &
   *R(80:82))
S(29:30) = Ez(1)*Ex(1)*(Ey(7)*R(71:72) - Ey(8)*R(76:77) + Ey(9)*R(80:81) - Ey(10) &
   *R(83:84))
S(31) = Ez(1)*Ex(1)*(Ey(7)*R(76) - Ey(8)*R(80) + Ey(9)*R(83) - Ey(10)*R(85))
S(32:33) = Ez(1)*Ex(1)*(Ey(7)*R(86:87) - Ey(8)*R(91:92) + Ey(9)*R(95:96) - Ey(10) &
   *R(98:99))
S(34) = Ez(1)*Ex(1)*(Ey(7)*R(91) - Ey(8)*R(95) + Ey(9)*R(98) - Ey(10)*R(100))
S(35) = Ez(1)*Ex(1)*(Ey(7)*R(101) - Ey(8)*R(105) + Ey(9)*R(108) - Ey(10)*R(110))
end subroutine auto2e_KetTransform_4_0_3_0

subroutine auto2e_KetTransform_4_0_2_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(0,2,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(2)*Ex(1)*(Ey(4)*R(1:5) - Ey(5)*R(9:13) + Ey(6)*R(16:20)) + Ez(3)*Ex(1) &
   *(-Ey(4)*R(37:41) + Ey(5)*R(44:48) - Ey(6)*R(50:54))
S(6:9) = Ez(2)*Ex(1)*(Ey(4)*R(9:12) - Ey(5)*R(16:19) + Ey(6)*R(22:25)) + Ez(3)*Ex(1) &
   *(-Ey(4)*R(44:47) + Ey(5)*R(50:53) - Ey(6)*R(55:58))
S(10:12) = Ez(2)*Ex(1)*(Ey(4)*R(16:18) - Ey(5)*R(22:24) + Ey(6)*R(27:29)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(50:52) + Ey(5)*R(55:57) - Ey(6)*R(59:61))
S(13:14) = Ez(2)*Ex(1)*(Ey(4)*R(22:23) - Ey(5)*R(27:28) + Ey(6)*R(31:32)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(55:56) + Ey(5)*R(59:60) - Ey(6)*R(62:63))
S(15) = Ez(2)*Ex(1)*(Ey(4)*R(27) - Ey(5)*R(31) + Ey(6)*R(34)) + Ez(3)*Ex(1)*(-Ey(4) &
   *R(59) + Ey(5)*R(62) - Ey(6)*R(64))
S(16:19) = Ez(2)*Ex(1)*(Ey(4)*R(37:40) - Ey(5)*R(44:47) + Ey(6)*R(50:53)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(65:68) + Ey(5)*R(71:74) - Ey(6)*R(76:79))
S(20:22) = Ez(2)*Ex(1)*(Ey(4)*R(44:46) - Ey(5)*R(50:52) + Ey(6)*R(55:57)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(71:73) + Ey(5)*R(76:78) - Ey(6)*R(80:82))
S(23:24) = Ez(2)*Ex(1)*(Ey(4)*R(50:51) - Ey(5)*R(55:56) + Ey(6)*R(59:60)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(76:77) + Ey(5)*R(80:81) - Ey(6)*R(83:84))
S(25) = Ez(2)*Ex(1)*(Ey(4)*R(55) - Ey(5)*R(59) + Ey(6)*R(62)) + Ez(3)*Ex(1)*(-Ey(4) &
   *R(80) + Ey(5)*R(83) - Ey(6)*R(85))
S(26:28) = Ez(2)*Ex(1)*(Ey(4)*R(65:67) - Ey(5)*R(71:73) + Ey(6)*R(76:78)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(86:88) + Ey(5)*R(91:93) - Ey(6)*R(95:97))
S(29:30) = Ez(2)*Ex(1)*(Ey(4)*R(71:72) - Ey(5)*R(76:77) + Ey(6)*R(80:81)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(91:92) + Ey(5)*R(95:96) - Ey(6)*R(98:99))
S(31) = Ez(2)*Ex(1)*(Ey(4)*R(76) - Ey(5)*R(80) + Ey(6)*R(83)) + Ez(3)*Ex(1)*(-Ey(4) &
   *R(95) + Ey(5)*R(98) - Ey(6)*R(100))
S(32:33) = Ez(2)*Ex(1)*(Ey(4)*R(86:87) - Ey(5)*R(91:92) + Ey(6)*R(95:96)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(101:102) + Ey(5)*R(105:106) - Ey(6)*R(108:109))
S(34) = Ez(2)*Ex(1)*(Ey(4)*R(91) - Ey(5)*R(95) + Ey(6)*R(98)) + Ez(3)*Ex(1)*(-Ey(4) &
   *R(105) + Ey(5)*R(108) - Ey(6)*R(110))
S(35) = Ez(2)*Ex(1)*(Ey(4)*R(101) - Ey(5)*R(105) + Ey(6)*R(108)) + Ez(3)*Ex(1)*( &
   -Ey(4)*R(111) + Ey(5)*R(114) - Ey(6)*R(116))
end subroutine auto2e_KetTransform_4_0_2_1

subroutine auto2e_KetTransform_4_0_1_2(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(0,1,2)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(4)*Ex(1)*(Ey(2)*R(1:5) - Ey(3)*R(9:13)) + Ez(5)*Ex(1)*(-Ey(2)*R(37:41)  &
   + Ey(3)*R(44:48)) + Ex(1)*Ez(6)*(Ey(2)*R(65:69) - Ey(3)*R(71:75))
S(6:9) = Ez(4)*Ex(1)*(Ey(2)*R(9:12) - Ey(3)*R(16:19)) + Ez(5)*Ex(1)*(-Ey(2)*R(44:47)  &
   + Ey(3)*R(50:53)) + Ex(1)*Ez(6)*(Ey(2)*R(71:74) - Ey(3)*R(76:79))
S(10:12) = Ez(4)*Ex(1)*(Ey(2)*R(16:18) - Ey(3)*R(22:24)) + Ez(5)*Ex(1)*(-Ey(2)*R(50:52)  &
   + Ey(3)*R(55:57)) + Ex(1)*Ez(6)*(Ey(2)*R(76:78) - Ey(3)*R(80:82))
S(13:14) = Ez(4)*Ex(1)*(Ey(2)*R(22:23) - Ey(3)*R(27:28)) + Ez(5)*Ex(1)*(-Ey(2)*R(55:56)  &
   + Ey(3)*R(59:60)) + Ex(1)*Ez(6)*(Ey(2)*R(80:81) - Ey(3)*R(83:84))
S(15) = Ez(4)*Ex(1)*(Ey(2)*R(27) - Ey(3)*R(31)) + Ez(5)*Ex(1)*(-Ey(2)*R(59) + Ey(3) &
   *R(62)) + Ex(1)*Ez(6)*(Ey(2)*R(83) - Ey(3)*R(85))
S(16:19) = Ez(4)*Ex(1)*(Ey(2)*R(37:40) - Ey(3)*R(44:47)) + Ez(5)*Ex(1)*(-Ey(2)*R(65:68)  &
   + Ey(3)*R(71:74)) + Ex(1)*Ez(6)*(Ey(2)*R(86:89) - Ey(3)*R(91:94))
S(20:22) = Ez(4)*Ex(1)*(Ey(2)*R(44:46) - Ey(3)*R(50:52)) + Ez(5)*Ex(1)*(-Ey(2)*R(71:73)  &
   + Ey(3)*R(76:78)) + Ex(1)*Ez(6)*(Ey(2)*R(91:93) - Ey(3)*R(95:97))
S(23:24) = Ez(4)*Ex(1)*(Ey(2)*R(50:51) - Ey(3)*R(55:56)) + Ez(5)*Ex(1)*(-Ey(2)*R(76:77)  &
   + Ey(3)*R(80:81)) + Ex(1)*Ez(6)*(Ey(2)*R(95:96) - Ey(3)*R(98:99))
S(25) = Ez(4)*Ex(1)*(Ey(2)*R(55) - Ey(3)*R(59)) + Ez(5)*Ex(1)*(-Ey(2)*R(80) + Ey(3) &
   *R(83)) + Ex(1)*Ez(6)*(Ey(2)*R(98) - Ey(3)*R(100))
S(26:28) = Ez(4)*Ex(1)*(Ey(2)*R(65:67) - Ey(3)*R(71:73)) + Ez(5)*Ex(1)*(-Ey(2)*R(86:88)  &
   + Ey(3)*R(91:93)) + Ex(1)*Ez(6)*(Ey(2)*R(101:103) - Ey(3)*R(105:107))
S(29:30) = Ez(4)*Ex(1)*(Ey(2)*R(71:72) - Ey(3)*R(76:77)) + Ez(5)*Ex(1)*(-Ey(2)*R(91:92)  &
   + Ey(3)*R(95:96)) + Ex(1)*Ez(6)*(Ey(2)*R(105:106) - Ey(3)*R(108:109))
S(31) = Ez(4)*Ex(1)*(Ey(2)*R(76) - Ey(3)*R(80)) + Ez(5)*Ex(1)*(-Ey(2)*R(95) + Ey(3) &
   *R(98)) + Ex(1)*Ez(6)*(Ey(2)*R(108) - Ey(3)*R(110))
S(32:33) = Ez(4)*Ex(1)*(Ey(2)*R(86:87) - Ey(3)*R(91:92)) + Ez(5)*Ex(1)*(-Ey(2)*R(101:102)  &
   + Ey(3)*R(105:106)) + Ex(1)*Ez(6)*(Ey(2)*R(111:112) - Ey(3)*R(114:115))
S(34) = Ez(4)*Ex(1)*(Ey(2)*R(91) - Ey(3)*R(95)) + Ez(5)*Ex(1)*(-Ey(2)*R(105) + Ey(3) &
   *R(108)) + Ex(1)*Ez(6)*(Ey(2)*R(114) - Ey(3)*R(116))
S(35) = Ez(4)*Ex(1)*(Ey(2)*R(101) - Ey(3)*R(105)) + Ez(5)*Ex(1)*(-Ey(2)*R(111) + Ey(3) &
   *R(114)) + Ex(1)*Ez(6)*(Ey(2)*R(117) - Ey(3)*R(119))
end subroutine auto2e_KetTransform_4_0_1_2

subroutine auto2e_KetTransform_4_0_0_3(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(0,0,3)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ex(1)*Ey(1)*(Ez(7)*R(1:5) - Ez(8)*R(37:41) + Ez(9)*R(65:69) - Ez(10)*R(86:90))
S(6:9) = Ex(1)*Ey(1)*(Ez(7)*R(9:12) - Ez(8)*R(44:47) + Ez(9)*R(71:74) - Ez(10)*R(91:94))
S(10:12) = Ex(1)*Ey(1)*(Ez(7)*R(16:18) - Ez(8)*R(50:52) + Ez(9)*R(76:78) - Ez(10) &
   *R(95:97))
S(13:14) = Ex(1)*Ey(1)*(Ez(7)*R(22:23) - Ez(8)*R(55:56) + Ez(9)*R(80:81) - Ez(10) &
   *R(98:99))
S(15) = Ex(1)*Ey(1)*(Ez(7)*R(27) - Ez(8)*R(59) + Ez(9)*R(83) - Ez(10)*R(100))
S(16:19) = Ex(1)*Ey(1)*(Ez(7)*R(37:40) - Ez(8)*R(65:68) + Ez(9)*R(86:89) - Ez(10) &
   *R(101:104))
S(20:22) = Ex(1)*Ey(1)*(Ez(7)*R(44:46) - Ez(8)*R(71:73) + Ez(9)*R(91:93) - Ez(10) &
   *R(105:107))
S(23:24) = Ex(1)*Ey(1)*(Ez(7)*R(50:51) - Ez(8)*R(76:77) + Ez(9)*R(95:96) - Ez(10) &
   *R(108:109))
S(25) = Ex(1)*Ey(1)*(Ez(7)*R(55) - Ez(8)*R(80) + Ez(9)*R(98) - Ez(10)*R(110))
S(26:28) = Ex(1)*Ey(1)*(Ez(7)*R(65:67) - Ez(8)*R(86:88) + Ez(9)*R(101:103) - Ez(10) &
   *R(111:113))
S(29:30) = Ex(1)*Ey(1)*(Ez(7)*R(71:72) - Ez(8)*R(91:92) + Ez(9)*R(105:106) - Ez(10) &
   *R(114:115))
S(31) = Ex(1)*Ey(1)*(Ez(7)*R(76) - Ez(8)*R(95) + Ez(9)*R(108) - Ez(10)*R(116))
S(32:33) = Ex(1)*Ey(1)*(Ez(7)*R(86:87) - Ez(8)*R(101:102) + Ez(9)*R(111:112) - Ez(10) &
   *R(117:118))
S(34) = Ex(1)*Ey(1)*(Ez(7)*R(91) - Ez(8)*R(105) + Ez(9)*R(114) - Ez(10)*R(119))
S(35) = Ex(1)*Ey(1)*(Ez(7)*R(101) - Ez(8)*R(111) + Ez(9)*R(117) - Ez(10)*R(120))
end subroutine auto2e_KetTransform_4_0_0_3
end module auto2e_KetTransform_4_3