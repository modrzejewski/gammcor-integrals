module auto2e_KetTransform_4_2
use arithmetic
use math_constants
implicit none
contains

subroutine auto2e_KetTransform_4_2_0_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(2,0,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(1)*Ey(1)*(Ex(4)*R(1:5) - Ex(5)*R(2:6) + Ex(6)*R(3:7))
S(6:9) = Ez(1)*Ey(1)*(Ex(4)*R(8:11) - Ex(5)*R(9:12) + Ex(6)*R(10:13))
S(10:12) = Ez(1)*Ey(1)*(Ex(4)*R(14:16) - Ex(5)*R(15:17) + Ex(6)*R(16:18))
S(13:14) = Ez(1)*Ey(1)*(Ex(4)*R(19:20) - Ex(5)*R(20:21) + Ex(6)*R(21:22))
S(15) = Ez(1)*Ey(1)*(Ex(4)*R(23) - Ex(5)*R(24) + Ex(6)*R(25))
S(16:19) = Ez(1)*Ey(1)*(Ex(4)*R(29:32) - Ex(5)*R(30:33) + Ex(6)*R(31:34))
S(20:22) = Ez(1)*Ey(1)*(Ex(4)*R(35:37) - Ex(5)*R(36:38) + Ex(6)*R(37:39))
S(23:24) = Ez(1)*Ey(1)*(Ex(4)*R(40:41) - Ex(5)*R(41:42) + Ex(6)*R(42:43))
S(25) = Ez(1)*Ey(1)*(Ex(4)*R(44) - Ex(5)*R(45) + Ex(6)*R(46))
S(26:28) = Ez(1)*Ey(1)*(Ex(4)*R(50:52) - Ex(5)*R(51:53) + Ex(6)*R(52:54))
S(29:30) = Ez(1)*Ey(1)*(Ex(4)*R(55:56) - Ex(5)*R(56:57) + Ex(6)*R(57:58))
S(31) = Ez(1)*Ey(1)*(Ex(4)*R(59) - Ex(5)*R(60) + Ex(6)*R(61))
S(32:33) = Ez(1)*Ey(1)*(Ex(4)*R(65:66) - Ex(5)*R(66:67) + Ex(6)*R(67:68))
S(34) = Ez(1)*Ey(1)*(Ex(4)*R(69) - Ex(5)*R(70) + Ex(6)*R(71))
S(35) = Ez(1)*Ey(1)*(Ex(4)*R(75) - Ex(5)*R(76) + Ex(6)*R(77))
end subroutine auto2e_KetTransform_4_2_0_0

subroutine auto2e_KetTransform_4_1_1_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(1,1,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(1)*Ey(2)*(Ex(2)*R(1:5) - Ex(3)*R(2:6)) + Ez(1)*Ey(3)*(-Ex(2)*R(8:12)  &
   + Ex(3)*R(9:13))
S(6:9) = Ez(1)*Ey(2)*(Ex(2)*R(8:11) - Ex(3)*R(9:12)) + Ez(1)*Ey(3)*(-Ex(2)*R(14:17)  &
   + Ex(3)*R(15:18))
S(10:12) = Ez(1)*Ey(2)*(Ex(2)*R(14:16) - Ex(3)*R(15:17)) + Ez(1)*Ey(3)*(-Ex(2)*R(19:21)  &
   + Ex(3)*R(20:22))
S(13:14) = Ez(1)*Ey(2)*(Ex(2)*R(19:20) - Ex(3)*R(20:21)) + Ez(1)*Ey(3)*(-Ex(2)*R(23:24)  &
   + Ex(3)*R(24:25))
S(15) = Ez(1)*Ey(2)*(Ex(2)*R(23) - Ex(3)*R(24)) + Ez(1)*Ey(3)*(-Ex(2)*R(26) + Ex(3) &
   *R(27))
S(16:19) = Ez(1)*Ey(2)*(Ex(2)*R(29:32) - Ex(3)*R(30:33)) + Ez(1)*Ey(3)*(-Ex(2)*R(35:38)  &
   + Ex(3)*R(36:39))
S(20:22) = Ez(1)*Ey(2)*(Ex(2)*R(35:37) - Ex(3)*R(36:38)) + Ez(1)*Ey(3)*(-Ex(2)*R(40:42)  &
   + Ex(3)*R(41:43))
S(23:24) = Ez(1)*Ey(2)*(Ex(2)*R(40:41) - Ex(3)*R(41:42)) + Ez(1)*Ey(3)*(-Ex(2)*R(44:45)  &
   + Ex(3)*R(45:46))
S(25) = Ez(1)*Ey(2)*(Ex(2)*R(44) - Ex(3)*R(45)) + Ez(1)*Ey(3)*(-Ex(2)*R(47) + Ex(3) &
   *R(48))
S(26:28) = Ez(1)*Ey(2)*(Ex(2)*R(50:52) - Ex(3)*R(51:53)) + Ez(1)*Ey(3)*(-Ex(2)*R(55:57)  &
   + Ex(3)*R(56:58))
S(29:30) = Ez(1)*Ey(2)*(Ex(2)*R(55:56) - Ex(3)*R(56:57)) + Ez(1)*Ey(3)*(-Ex(2)*R(59:60)  &
   + Ex(3)*R(60:61))
S(31) = Ez(1)*Ey(2)*(Ex(2)*R(59) - Ex(3)*R(60)) + Ez(1)*Ey(3)*(-Ex(2)*R(62) + Ex(3) &
   *R(63))
S(32:33) = Ez(1)*Ey(2)*(Ex(2)*R(65:66) - Ex(3)*R(66:67)) + Ez(1)*Ey(3)*(-Ex(2)*R(69:70)  &
   + Ex(3)*R(70:71))
S(34) = Ez(1)*Ey(2)*(Ex(2)*R(69) - Ex(3)*R(70)) + Ez(1)*Ey(3)*(-Ex(2)*R(72) + Ex(3) &
   *R(73))
S(35) = Ez(1)*Ey(2)*(Ex(2)*R(75) - Ex(3)*R(76)) + Ez(1)*Ey(3)*(-Ex(2)*R(78) + Ex(3) &
   *R(79))
end subroutine auto2e_KetTransform_4_1_1_0

subroutine auto2e_KetTransform_4_1_0_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(1,0,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(2)*Ey(1)*(Ex(2)*R(1:5) - Ex(3)*R(2:6)) + Ez(3)*Ey(1)*(-Ex(2)*R(29:33)  &
   + Ex(3)*R(30:34))
S(6:9) = Ez(2)*Ey(1)*(Ex(2)*R(8:11) - Ex(3)*R(9:12)) + Ez(3)*Ey(1)*(-Ex(2)*R(35:38)  &
   + Ex(3)*R(36:39))
S(10:12) = Ez(2)*Ey(1)*(Ex(2)*R(14:16) - Ex(3)*R(15:17)) + Ez(3)*Ey(1)*(-Ex(2)*R(40:42)  &
   + Ex(3)*R(41:43))
S(13:14) = Ez(2)*Ey(1)*(Ex(2)*R(19:20) - Ex(3)*R(20:21)) + Ez(3)*Ey(1)*(-Ex(2)*R(44:45)  &
   + Ex(3)*R(45:46))
S(15) = Ez(2)*Ey(1)*(Ex(2)*R(23) - Ex(3)*R(24)) + Ez(3)*Ey(1)*(-Ex(2)*R(47) + Ex(3) &
   *R(48))
S(16:19) = Ez(2)*Ey(1)*(Ex(2)*R(29:32) - Ex(3)*R(30:33)) + Ez(3)*Ey(1)*(-Ex(2)*R(50:53)  &
   + Ex(3)*R(51:54))
S(20:22) = Ez(2)*Ey(1)*(Ex(2)*R(35:37) - Ex(3)*R(36:38)) + Ez(3)*Ey(1)*(-Ex(2)*R(55:57)  &
   + Ex(3)*R(56:58))
S(23:24) = Ez(2)*Ey(1)*(Ex(2)*R(40:41) - Ex(3)*R(41:42)) + Ez(3)*Ey(1)*(-Ex(2)*R(59:60)  &
   + Ex(3)*R(60:61))
S(25) = Ez(2)*Ey(1)*(Ex(2)*R(44) - Ex(3)*R(45)) + Ez(3)*Ey(1)*(-Ex(2)*R(62) + Ex(3) &
   *R(63))
S(26:28) = Ez(2)*Ey(1)*(Ex(2)*R(50:52) - Ex(3)*R(51:53)) + Ez(3)*Ey(1)*(-Ex(2)*R(65:67)  &
   + Ex(3)*R(66:68))
S(29:30) = Ez(2)*Ey(1)*(Ex(2)*R(55:56) - Ex(3)*R(56:57)) + Ez(3)*Ey(1)*(-Ex(2)*R(69:70)  &
   + Ex(3)*R(70:71))
S(31) = Ez(2)*Ey(1)*(Ex(2)*R(59) - Ex(3)*R(60)) + Ez(3)*Ey(1)*(-Ex(2)*R(72) + Ex(3) &
   *R(73))
S(32:33) = Ez(2)*Ey(1)*(Ex(2)*R(65:66) - Ex(3)*R(66:67)) + Ez(3)*Ey(1)*(-Ex(2)*R(75:76)  &
   + Ex(3)*R(76:77))
S(34) = Ez(2)*Ey(1)*(Ex(2)*R(69) - Ex(3)*R(70)) + Ez(3)*Ey(1)*(-Ex(2)*R(78) + Ex(3) &
   *R(79))
S(35) = Ez(2)*Ey(1)*(Ex(2)*R(75) - Ex(3)*R(76)) + Ez(3)*Ey(1)*(-Ex(2)*R(81) + Ex(3) &
   *R(82))
end subroutine auto2e_KetTransform_4_1_0_1

subroutine auto2e_KetTransform_4_0_2_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(0,2,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(1)*Ex(1)*(Ey(4)*R(1:5) - Ey(5)*R(8:12) + Ey(6)*R(14:18))
S(6:9) = Ez(1)*Ex(1)*(Ey(4)*R(8:11) - Ey(5)*R(14:17) + Ey(6)*R(19:22))
S(10:12) = Ez(1)*Ex(1)*(Ey(4)*R(14:16) - Ey(5)*R(19:21) + Ey(6)*R(23:25))
S(13:14) = Ez(1)*Ex(1)*(Ey(4)*R(19:20) - Ey(5)*R(23:24) + Ey(6)*R(26:27))
S(15) = Ez(1)*Ex(1)*(Ey(4)*R(23) - Ey(5)*R(26) + Ey(6)*R(28))
S(16:19) = Ez(1)*Ex(1)*(Ey(4)*R(29:32) - Ey(5)*R(35:38) + Ey(6)*R(40:43))
S(20:22) = Ez(1)*Ex(1)*(Ey(4)*R(35:37) - Ey(5)*R(40:42) + Ey(6)*R(44:46))
S(23:24) = Ez(1)*Ex(1)*(Ey(4)*R(40:41) - Ey(5)*R(44:45) + Ey(6)*R(47:48))
S(25) = Ez(1)*Ex(1)*(Ey(4)*R(44) - Ey(5)*R(47) + Ey(6)*R(49))
S(26:28) = Ez(1)*Ex(1)*(Ey(4)*R(50:52) - Ey(5)*R(55:57) + Ey(6)*R(59:61))
S(29:30) = Ez(1)*Ex(1)*(Ey(4)*R(55:56) - Ey(5)*R(59:60) + Ey(6)*R(62:63))
S(31) = Ez(1)*Ex(1)*(Ey(4)*R(59) - Ey(5)*R(62) + Ey(6)*R(64))
S(32:33) = Ez(1)*Ex(1)*(Ey(4)*R(65:66) - Ey(5)*R(69:70) + Ey(6)*R(72:73))
S(34) = Ez(1)*Ex(1)*(Ey(4)*R(69) - Ey(5)*R(72) + Ey(6)*R(74))
S(35) = Ez(1)*Ex(1)*(Ey(4)*R(75) - Ey(5)*R(78) + Ey(6)*R(80))
end subroutine auto2e_KetTransform_4_0_2_0

subroutine auto2e_KetTransform_4_0_1_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(0,1,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ez(2)*Ex(1)*(Ey(2)*R(1:5) - Ey(3)*R(8:12)) + Ez(3)*Ex(1)*(-Ey(2)*R(29:33)  &
   + Ey(3)*R(35:39))
S(6:9) = Ez(2)*Ex(1)*(Ey(2)*R(8:11) - Ey(3)*R(14:17)) + Ez(3)*Ex(1)*(-Ey(2)*R(35:38)  &
   + Ey(3)*R(40:43))
S(10:12) = Ez(2)*Ex(1)*(Ey(2)*R(14:16) - Ey(3)*R(19:21)) + Ez(3)*Ex(1)*(-Ey(2)*R(40:42)  &
   + Ey(3)*R(44:46))
S(13:14) = Ez(2)*Ex(1)*(Ey(2)*R(19:20) - Ey(3)*R(23:24)) + Ez(3)*Ex(1)*(-Ey(2)*R(44:45)  &
   + Ey(3)*R(47:48))
S(15) = Ez(2)*Ex(1)*(Ey(2)*R(23) - Ey(3)*R(26)) + Ez(3)*Ex(1)*(-Ey(2)*R(47) + Ey(3) &
   *R(49))
S(16:19) = Ez(2)*Ex(1)*(Ey(2)*R(29:32) - Ey(3)*R(35:38)) + Ez(3)*Ex(1)*(-Ey(2)*R(50:53)  &
   + Ey(3)*R(55:58))
S(20:22) = Ez(2)*Ex(1)*(Ey(2)*R(35:37) - Ey(3)*R(40:42)) + Ez(3)*Ex(1)*(-Ey(2)*R(55:57)  &
   + Ey(3)*R(59:61))
S(23:24) = Ez(2)*Ex(1)*(Ey(2)*R(40:41) - Ey(3)*R(44:45)) + Ez(3)*Ex(1)*(-Ey(2)*R(59:60)  &
   + Ey(3)*R(62:63))
S(25) = Ez(2)*Ex(1)*(Ey(2)*R(44) - Ey(3)*R(47)) + Ez(3)*Ex(1)*(-Ey(2)*R(62) + Ey(3) &
   *R(64))
S(26:28) = Ez(2)*Ex(1)*(Ey(2)*R(50:52) - Ey(3)*R(55:57)) + Ez(3)*Ex(1)*(-Ey(2)*R(65:67)  &
   + Ey(3)*R(69:71))
S(29:30) = Ez(2)*Ex(1)*(Ey(2)*R(55:56) - Ey(3)*R(59:60)) + Ez(3)*Ex(1)*(-Ey(2)*R(69:70)  &
   + Ey(3)*R(72:73))
S(31) = Ez(2)*Ex(1)*(Ey(2)*R(59) - Ey(3)*R(62)) + Ez(3)*Ex(1)*(-Ey(2)*R(72) + Ey(3) &
   *R(74))
S(32:33) = Ez(2)*Ex(1)*(Ey(2)*R(65:66) - Ey(3)*R(69:70)) + Ez(3)*Ex(1)*(-Ey(2)*R(75:76)  &
   + Ey(3)*R(78:79))
S(34) = Ez(2)*Ex(1)*(Ey(2)*R(69) - Ey(3)*R(72)) + Ez(3)*Ex(1)*(-Ey(2)*R(78) + Ey(3) &
   *R(80))
S(35) = Ez(2)*Ex(1)*(Ey(2)*R(75) - Ey(3)*R(78)) + Ez(3)*Ex(1)*(-Ey(2)*R(81) + Ey(3) &
   *R(83))
end subroutine auto2e_KetTransform_4_0_1_1

subroutine auto2e_KetTransform_4_0_0_2(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=4
! Class=(LS|KS), L=4, K=(0,0,2)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:5) = Ex(1)*Ey(1)*(Ez(4)*R(1:5) - Ez(5)*R(29:33) + Ez(6)*R(50:54))
S(6:9) = Ex(1)*Ey(1)*(Ez(4)*R(8:11) - Ez(5)*R(35:38) + Ez(6)*R(55:58))
S(10:12) = Ex(1)*Ey(1)*(Ez(4)*R(14:16) - Ez(5)*R(40:42) + Ez(6)*R(59:61))
S(13:14) = Ex(1)*Ey(1)*(Ez(4)*R(19:20) - Ez(5)*R(44:45) + Ez(6)*R(62:63))
S(15) = Ex(1)*Ey(1)*(Ez(4)*R(23) - Ez(5)*R(47) + Ez(6)*R(64))
S(16:19) = Ex(1)*Ey(1)*(Ez(4)*R(29:32) - Ez(5)*R(50:53) + Ez(6)*R(65:68))
S(20:22) = Ex(1)*Ey(1)*(Ez(4)*R(35:37) - Ez(5)*R(55:57) + Ez(6)*R(69:71))
S(23:24) = Ex(1)*Ey(1)*(Ez(4)*R(40:41) - Ez(5)*R(59:60) + Ez(6)*R(72:73))
S(25) = Ex(1)*Ey(1)*(Ez(4)*R(44) - Ez(5)*R(62) + Ez(6)*R(74))
S(26:28) = Ex(1)*Ey(1)*(Ez(4)*R(50:52) - Ez(5)*R(65:67) + Ez(6)*R(75:77))
S(29:30) = Ex(1)*Ey(1)*(Ez(4)*R(55:56) - Ez(5)*R(69:70) + Ez(6)*R(78:79))
S(31) = Ex(1)*Ey(1)*(Ez(4)*R(59) - Ez(5)*R(72) + Ez(6)*R(80))
S(32:33) = Ex(1)*Ey(1)*(Ez(4)*R(65:66) - Ez(5)*R(75:76) + Ez(6)*R(81:82))
S(34) = Ex(1)*Ey(1)*(Ez(4)*R(69) - Ez(5)*R(78) + Ez(6)*R(83))
S(35) = Ex(1)*Ey(1)*(Ez(4)*R(75) - Ez(5)*R(81) + Ez(6)*R(84))
end subroutine auto2e_KetTransform_4_0_0_2
end module auto2e_KetTransform_4_2