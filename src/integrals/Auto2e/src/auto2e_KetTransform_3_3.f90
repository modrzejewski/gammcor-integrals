module auto2e_KetTransform_3_3
use arithmetic
use math_constants
implicit none
contains

subroutine auto2e_KetTransform_3_3_0_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(3,0,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(1)*Ey(1)*(Ex(7)*R(1:4) - Ex(8)*R(2:5) + Ex(9)*R(3:6) - Ex(10)*R(4:7))
S(5:7) = Ez(1)*Ey(1)*(Ex(7)*R(8:10) - Ex(8)*R(9:11) + Ex(9)*R(10:12) - Ex(10)*R(11:13))
S(8:9) = Ez(1)*Ey(1)*(Ex(7)*R(14:15) - Ex(8)*R(15:16) + Ex(9)*R(16:17) - Ex(10)*R(17:18))
S(10) = Ez(1)*Ey(1)*(Ex(7)*R(19) - Ex(8)*R(20) + Ex(9)*R(21) - Ex(10)*R(22))
S(11:13) = Ez(1)*Ey(1)*(Ex(7)*R(29:31) - Ex(8)*R(30:32) + Ex(9)*R(31:33) - Ex(10) &
   *R(32:34))
S(14:15) = Ez(1)*Ey(1)*(Ex(7)*R(35:36) - Ex(8)*R(36:37) + Ex(9)*R(37:38) - Ex(10) &
   *R(38:39))
S(16) = Ez(1)*Ey(1)*(Ex(7)*R(40) - Ex(8)*R(41) + Ex(9)*R(42) - Ex(10)*R(43))
S(17:18) = Ez(1)*Ey(1)*(Ex(7)*R(50:51) - Ex(8)*R(51:52) + Ex(9)*R(52:53) - Ex(10) &
   *R(53:54))
S(19) = Ez(1)*Ey(1)*(Ex(7)*R(55) - Ex(8)*R(56) + Ex(9)*R(57) - Ex(10)*R(58))
S(20) = Ez(1)*Ey(1)*(Ex(7)*R(65) - Ex(8)*R(66) + Ex(9)*R(67) - Ex(10)*R(68))
end subroutine auto2e_KetTransform_3_3_0_0

subroutine auto2e_KetTransform_3_2_1_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(2,1,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(1)*Ey(2)*(Ex(4)*R(1:4) - Ex(5)*R(2:5) + Ex(6)*R(3:6)) + Ez(1)*Ey(3)*( &
   -Ex(4)*R(8:11) + Ex(5)*R(9:12) - Ex(6)*R(10:13))
S(5:7) = Ez(1)*Ey(2)*(Ex(4)*R(8:10) - Ex(5)*R(9:11) + Ex(6)*R(10:12)) + Ez(1)*Ey(3) &
   *(-Ex(4)*R(14:16) + Ex(5)*R(15:17) - Ex(6)*R(16:18))
S(8:9) = Ez(1)*Ey(2)*(Ex(4)*R(14:15) - Ex(5)*R(15:16) + Ex(6)*R(16:17)) + Ez(1)*Ey(3) &
   *(-Ex(4)*R(19:20) + Ex(5)*R(20:21) - Ex(6)*R(21:22))
S(10) = Ez(1)*Ey(2)*(Ex(4)*R(19) - Ex(5)*R(20) + Ex(6)*R(21)) + Ez(1)*Ey(3)*(-Ex(4) &
   *R(23) + Ex(5)*R(24) - Ex(6)*R(25))
S(11:13) = Ez(1)*Ey(2)*(Ex(4)*R(29:31) - Ex(5)*R(30:32) + Ex(6)*R(31:33)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(35:37) + Ex(5)*R(36:38) - Ex(6)*R(37:39))
S(14:15) = Ez(1)*Ey(2)*(Ex(4)*R(35:36) - Ex(5)*R(36:37) + Ex(6)*R(37:38)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(40:41) + Ex(5)*R(41:42) - Ex(6)*R(42:43))
S(16) = Ez(1)*Ey(2)*(Ex(4)*R(40) - Ex(5)*R(41) + Ex(6)*R(42)) + Ez(1)*Ey(3)*(-Ex(4) &
   *R(44) + Ex(5)*R(45) - Ex(6)*R(46))
S(17:18) = Ez(1)*Ey(2)*(Ex(4)*R(50:51) - Ex(5)*R(51:52) + Ex(6)*R(52:53)) + Ez(1) &
   *Ey(3)*(-Ex(4)*R(55:56) + Ex(5)*R(56:57) - Ex(6)*R(57:58))
S(19) = Ez(1)*Ey(2)*(Ex(4)*R(55) - Ex(5)*R(56) + Ex(6)*R(57)) + Ez(1)*Ey(3)*(-Ex(4) &
   *R(59) + Ex(5)*R(60) - Ex(6)*R(61))
S(20) = Ez(1)*Ey(2)*(Ex(4)*R(65) - Ex(5)*R(66) + Ex(6)*R(67)) + Ez(1)*Ey(3)*(-Ex(4) &
   *R(69) + Ex(5)*R(70) - Ex(6)*R(71))
end subroutine auto2e_KetTransform_3_2_1_0

subroutine auto2e_KetTransform_3_2_0_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(2,0,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(2)*Ey(1)*(Ex(4)*R(1:4) - Ex(5)*R(2:5) + Ex(6)*R(3:6)) + Ez(3)*Ey(1)*( &
   -Ex(4)*R(29:32) + Ex(5)*R(30:33) - Ex(6)*R(31:34))
S(5:7) = Ez(2)*Ey(1)*(Ex(4)*R(8:10) - Ex(5)*R(9:11) + Ex(6)*R(10:12)) + Ez(3)*Ey(1) &
   *(-Ex(4)*R(35:37) + Ex(5)*R(36:38) - Ex(6)*R(37:39))
S(8:9) = Ez(2)*Ey(1)*(Ex(4)*R(14:15) - Ex(5)*R(15:16) + Ex(6)*R(16:17)) + Ez(3)*Ey(1) &
   *(-Ex(4)*R(40:41) + Ex(5)*R(41:42) - Ex(6)*R(42:43))
S(10) = Ez(2)*Ey(1)*(Ex(4)*R(19) - Ex(5)*R(20) + Ex(6)*R(21)) + Ez(3)*Ey(1)*(-Ex(4) &
   *R(44) + Ex(5)*R(45) - Ex(6)*R(46))
S(11:13) = Ez(2)*Ey(1)*(Ex(4)*R(29:31) - Ex(5)*R(30:32) + Ex(6)*R(31:33)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(50:52) + Ex(5)*R(51:53) - Ex(6)*R(52:54))
S(14:15) = Ez(2)*Ey(1)*(Ex(4)*R(35:36) - Ex(5)*R(36:37) + Ex(6)*R(37:38)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(55:56) + Ex(5)*R(56:57) - Ex(6)*R(57:58))
S(16) = Ez(2)*Ey(1)*(Ex(4)*R(40) - Ex(5)*R(41) + Ex(6)*R(42)) + Ez(3)*Ey(1)*(-Ex(4) &
   *R(59) + Ex(5)*R(60) - Ex(6)*R(61))
S(17:18) = Ez(2)*Ey(1)*(Ex(4)*R(50:51) - Ex(5)*R(51:52) + Ex(6)*R(52:53)) + Ez(3) &
   *Ey(1)*(-Ex(4)*R(65:66) + Ex(5)*R(66:67) - Ex(6)*R(67:68))
S(19) = Ez(2)*Ey(1)*(Ex(4)*R(55) - Ex(5)*R(56) + Ex(6)*R(57)) + Ez(3)*Ey(1)*(-Ex(4) &
   *R(69) + Ex(5)*R(70) - Ex(6)*R(71))
S(20) = Ez(2)*Ey(1)*(Ex(4)*R(65) - Ex(5)*R(66) + Ex(6)*R(67)) + Ez(3)*Ey(1)*(-Ex(4) &
   *R(75) + Ex(5)*R(76) - Ex(6)*R(77))
end subroutine auto2e_KetTransform_3_2_0_1

subroutine auto2e_KetTransform_3_1_2_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(1,2,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(1)*Ey(4)*(Ex(2)*R(1:4) - Ex(3)*R(2:5)) + Ey(5)*Ez(1)*(-Ex(2)*R(8:11)  &
   + Ex(3)*R(9:12)) + Ez(1)*Ey(6)*(Ex(2)*R(14:17) - Ex(3)*R(15:18))
S(5:7) = Ez(1)*Ey(4)*(Ex(2)*R(8:10) - Ex(3)*R(9:11)) + Ey(5)*Ez(1)*(-Ex(2)*R(14:16)  &
   + Ex(3)*R(15:17)) + Ez(1)*Ey(6)*(Ex(2)*R(19:21) - Ex(3)*R(20:22))
S(8:9) = Ez(1)*Ey(4)*(Ex(2)*R(14:15) - Ex(3)*R(15:16)) + Ey(5)*Ez(1)*(-Ex(2)*R(19:20)  &
   + Ex(3)*R(20:21)) + Ez(1)*Ey(6)*(Ex(2)*R(23:24) - Ex(3)*R(24:25))
S(10) = Ez(1)*Ey(4)*(Ex(2)*R(19) - Ex(3)*R(20)) + Ey(5)*Ez(1)*(-Ex(2)*R(23) + Ex(3) &
   *R(24)) + Ez(1)*Ey(6)*(Ex(2)*R(26) - Ex(3)*R(27))
S(11:13) = Ez(1)*Ey(4)*(Ex(2)*R(29:31) - Ex(3)*R(30:32)) + Ey(5)*Ez(1)*(-Ex(2)*R(35:37)  &
   + Ex(3)*R(36:38)) + Ez(1)*Ey(6)*(Ex(2)*R(40:42) - Ex(3)*R(41:43))
S(14:15) = Ez(1)*Ey(4)*(Ex(2)*R(35:36) - Ex(3)*R(36:37)) + Ey(5)*Ez(1)*(-Ex(2)*R(40:41)  &
   + Ex(3)*R(41:42)) + Ez(1)*Ey(6)*(Ex(2)*R(44:45) - Ex(3)*R(45:46))
S(16) = Ez(1)*Ey(4)*(Ex(2)*R(40) - Ex(3)*R(41)) + Ey(5)*Ez(1)*(-Ex(2)*R(44) + Ex(3) &
   *R(45)) + Ez(1)*Ey(6)*(Ex(2)*R(47) - Ex(3)*R(48))
S(17:18) = Ez(1)*Ey(4)*(Ex(2)*R(50:51) - Ex(3)*R(51:52)) + Ey(5)*Ez(1)*(-Ex(2)*R(55:56)  &
   + Ex(3)*R(56:57)) + Ez(1)*Ey(6)*(Ex(2)*R(59:60) - Ex(3)*R(60:61))
S(19) = Ez(1)*Ey(4)*(Ex(2)*R(55) - Ex(3)*R(56)) + Ey(5)*Ez(1)*(-Ex(2)*R(59) + Ex(3) &
   *R(60)) + Ez(1)*Ey(6)*(Ex(2)*R(62) - Ex(3)*R(63))
S(20) = Ez(1)*Ey(4)*(Ex(2)*R(65) - Ex(3)*R(66)) + Ey(5)*Ez(1)*(-Ex(2)*R(69) + Ex(3) &
   *R(70)) + Ez(1)*Ey(6)*(Ex(2)*R(72) - Ex(3)*R(73))
end subroutine auto2e_KetTransform_3_1_2_0

subroutine auto2e_KetTransform_3_1_1_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(1,1,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(2)*Ey(2)*(Ex(2)*R(1:4) - Ex(3)*R(2:5)) + Ez(2)*Ey(3)*(-Ex(2)*R(8:11)  &
   + Ex(3)*R(9:12)) + Ez(3)*Ey(2)*(-Ex(2)*R(29:32) + Ex(3)*R(30:33)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(35:38) - Ex(3)*R(36:39))
S(5:7) = Ez(2)*Ey(2)*(Ex(2)*R(8:10) - Ex(3)*R(9:11)) + Ez(2)*Ey(3)*(-Ex(2)*R(14:16)  &
   + Ex(3)*R(15:17)) + Ez(3)*Ey(2)*(-Ex(2)*R(35:37) + Ex(3)*R(36:38)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(40:42) - Ex(3)*R(41:43))
S(8:9) = Ez(2)*Ey(2)*(Ex(2)*R(14:15) - Ex(3)*R(15:16)) + Ez(2)*Ey(3)*(-Ex(2)*R(19:20)  &
   + Ex(3)*R(20:21)) + Ez(3)*Ey(2)*(-Ex(2)*R(40:41) + Ex(3)*R(41:42)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(44:45) - Ex(3)*R(45:46))
S(10) = Ez(2)*Ey(2)*(Ex(2)*R(19) - Ex(3)*R(20)) + Ez(2)*Ey(3)*(-Ex(2)*R(23) + Ex(3) &
   *R(24)) + Ez(3)*Ey(2)*(-Ex(2)*R(44) + Ex(3)*R(45)) + Ey(3)*Ez(3)*(Ex(2)*R(47)  &
   - Ex(3)*R(48))
S(11:13) = Ez(2)*Ey(2)*(Ex(2)*R(29:31) - Ex(3)*R(30:32)) + Ez(2)*Ey(3)*(-Ex(2)*R(35:37)  &
   + Ex(3)*R(36:38)) + Ez(3)*Ey(2)*(-Ex(2)*R(50:52) + Ex(3)*R(51:53)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(55:57) - Ex(3)*R(56:58))
S(14:15) = Ez(2)*Ey(2)*(Ex(2)*R(35:36) - Ex(3)*R(36:37)) + Ez(2)*Ey(3)*(-Ex(2)*R(40:41)  &
   + Ex(3)*R(41:42)) + Ez(3)*Ey(2)*(-Ex(2)*R(55:56) + Ex(3)*R(56:57)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(59:60) - Ex(3)*R(60:61))
S(16) = Ez(2)*Ey(2)*(Ex(2)*R(40) - Ex(3)*R(41)) + Ez(2)*Ey(3)*(-Ex(2)*R(44) + Ex(3) &
   *R(45)) + Ez(3)*Ey(2)*(-Ex(2)*R(59) + Ex(3)*R(60)) + Ey(3)*Ez(3)*(Ex(2)*R(62)  &
   - Ex(3)*R(63))
S(17:18) = Ez(2)*Ey(2)*(Ex(2)*R(50:51) - Ex(3)*R(51:52)) + Ez(2)*Ey(3)*(-Ex(2)*R(55:56)  &
   + Ex(3)*R(56:57)) + Ez(3)*Ey(2)*(-Ex(2)*R(65:66) + Ex(3)*R(66:67)) + Ey(3)*Ez(3) &
   *(Ex(2)*R(69:70) - Ex(3)*R(70:71))
S(19) = Ez(2)*Ey(2)*(Ex(2)*R(55) - Ex(3)*R(56)) + Ez(2)*Ey(3)*(-Ex(2)*R(59) + Ex(3) &
   *R(60)) + Ez(3)*Ey(2)*(-Ex(2)*R(69) + Ex(3)*R(70)) + Ey(3)*Ez(3)*(Ex(2)*R(72)  &
   - Ex(3)*R(73))
S(20) = Ez(2)*Ey(2)*(Ex(2)*R(65) - Ex(3)*R(66)) + Ez(2)*Ey(3)*(-Ex(2)*R(69) + Ex(3) &
   *R(70)) + Ez(3)*Ey(2)*(-Ex(2)*R(75) + Ex(3)*R(76)) + Ey(3)*Ez(3)*(Ex(2)*R(78)  &
   - Ex(3)*R(79))
end subroutine auto2e_KetTransform_3_1_1_1

subroutine auto2e_KetTransform_3_1_0_2(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(1,0,2)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(4)*Ey(1)*(Ex(2)*R(1:4) - Ex(3)*R(2:5)) + Ez(5)*Ey(1)*(-Ex(2)*R(29:32)  &
   + Ex(3)*R(30:33)) + Ez(6)*Ey(1)*(Ex(2)*R(50:53) - Ex(3)*R(51:54))
S(5:7) = Ez(4)*Ey(1)*(Ex(2)*R(8:10) - Ex(3)*R(9:11)) + Ez(5)*Ey(1)*(-Ex(2)*R(35:37)  &
   + Ex(3)*R(36:38)) + Ez(6)*Ey(1)*(Ex(2)*R(55:57) - Ex(3)*R(56:58))
S(8:9) = Ez(4)*Ey(1)*(Ex(2)*R(14:15) - Ex(3)*R(15:16)) + Ez(5)*Ey(1)*(-Ex(2)*R(40:41)  &
   + Ex(3)*R(41:42)) + Ez(6)*Ey(1)*(Ex(2)*R(59:60) - Ex(3)*R(60:61))
S(10) = Ez(4)*Ey(1)*(Ex(2)*R(19) - Ex(3)*R(20)) + Ez(5)*Ey(1)*(-Ex(2)*R(44) + Ex(3) &
   *R(45)) + Ez(6)*Ey(1)*(Ex(2)*R(62) - Ex(3)*R(63))
S(11:13) = Ez(4)*Ey(1)*(Ex(2)*R(29:31) - Ex(3)*R(30:32)) + Ez(5)*Ey(1)*(-Ex(2)*R(50:52)  &
   + Ex(3)*R(51:53)) + Ez(6)*Ey(1)*(Ex(2)*R(65:67) - Ex(3)*R(66:68))
S(14:15) = Ez(4)*Ey(1)*(Ex(2)*R(35:36) - Ex(3)*R(36:37)) + Ez(5)*Ey(1)*(-Ex(2)*R(55:56)  &
   + Ex(3)*R(56:57)) + Ez(6)*Ey(1)*(Ex(2)*R(69:70) - Ex(3)*R(70:71))
S(16) = Ez(4)*Ey(1)*(Ex(2)*R(40) - Ex(3)*R(41)) + Ez(5)*Ey(1)*(-Ex(2)*R(59) + Ex(3) &
   *R(60)) + Ez(6)*Ey(1)*(Ex(2)*R(72) - Ex(3)*R(73))
S(17:18) = Ez(4)*Ey(1)*(Ex(2)*R(50:51) - Ex(3)*R(51:52)) + Ez(5)*Ey(1)*(-Ex(2)*R(65:66)  &
   + Ex(3)*R(66:67)) + Ez(6)*Ey(1)*(Ex(2)*R(75:76) - Ex(3)*R(76:77))
S(19) = Ez(4)*Ey(1)*(Ex(2)*R(55) - Ex(3)*R(56)) + Ez(5)*Ey(1)*(-Ex(2)*R(69) + Ex(3) &
   *R(70)) + Ez(6)*Ey(1)*(Ex(2)*R(78) - Ex(3)*R(79))
S(20) = Ez(4)*Ey(1)*(Ex(2)*R(65) - Ex(3)*R(66)) + Ez(5)*Ey(1)*(-Ex(2)*R(75) + Ex(3) &
   *R(76)) + Ez(6)*Ey(1)*(Ex(2)*R(81) - Ex(3)*R(82))
end subroutine auto2e_KetTransform_3_1_0_2

subroutine auto2e_KetTransform_3_0_3_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(0,3,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(1)*Ex(1)*(Ey(7)*R(1:4) - Ey(8)*R(8:11) + Ey(9)*R(14:17) - Ey(10)*R(19:22))
S(5:7) = Ez(1)*Ex(1)*(Ey(7)*R(8:10) - Ey(8)*R(14:16) + Ey(9)*R(19:21) - Ey(10)*R(23:25))
S(8:9) = Ez(1)*Ex(1)*(Ey(7)*R(14:15) - Ey(8)*R(19:20) + Ey(9)*R(23:24) - Ey(10)*R(26:27))
S(10) = Ez(1)*Ex(1)*(Ey(7)*R(19) - Ey(8)*R(23) + Ey(9)*R(26) - Ey(10)*R(28))
S(11:13) = Ez(1)*Ex(1)*(Ey(7)*R(29:31) - Ey(8)*R(35:37) + Ey(9)*R(40:42) - Ey(10) &
   *R(44:46))
S(14:15) = Ez(1)*Ex(1)*(Ey(7)*R(35:36) - Ey(8)*R(40:41) + Ey(9)*R(44:45) - Ey(10) &
   *R(47:48))
S(16) = Ez(1)*Ex(1)*(Ey(7)*R(40) - Ey(8)*R(44) + Ey(9)*R(47) - Ey(10)*R(49))
S(17:18) = Ez(1)*Ex(1)*(Ey(7)*R(50:51) - Ey(8)*R(55:56) + Ey(9)*R(59:60) - Ey(10) &
   *R(62:63))
S(19) = Ez(1)*Ex(1)*(Ey(7)*R(55) - Ey(8)*R(59) + Ey(9)*R(62) - Ey(10)*R(64))
S(20) = Ez(1)*Ex(1)*(Ey(7)*R(65) - Ey(8)*R(69) + Ey(9)*R(72) - Ey(10)*R(74))
end subroutine auto2e_KetTransform_3_0_3_0

subroutine auto2e_KetTransform_3_0_2_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(0,2,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(2)*Ex(1)*(Ey(4)*R(1:4) - Ey(5)*R(8:11) + Ey(6)*R(14:17)) + Ez(3)*Ex(1) &
   *(-Ey(4)*R(29:32) + Ey(5)*R(35:38) - Ey(6)*R(40:43))
S(5:7) = Ez(2)*Ex(1)*(Ey(4)*R(8:10) - Ey(5)*R(14:16) + Ey(6)*R(19:21)) + Ez(3)*Ex(1) &
   *(-Ey(4)*R(35:37) + Ey(5)*R(40:42) - Ey(6)*R(44:46))
S(8:9) = Ez(2)*Ex(1)*(Ey(4)*R(14:15) - Ey(5)*R(19:20) + Ey(6)*R(23:24)) + Ez(3)*Ex(1) &
   *(-Ey(4)*R(40:41) + Ey(5)*R(44:45) - Ey(6)*R(47:48))
S(10) = Ez(2)*Ex(1)*(Ey(4)*R(19) - Ey(5)*R(23) + Ey(6)*R(26)) + Ez(3)*Ex(1)*(-Ey(4) &
   *R(44) + Ey(5)*R(47) - Ey(6)*R(49))
S(11:13) = Ez(2)*Ex(1)*(Ey(4)*R(29:31) - Ey(5)*R(35:37) + Ey(6)*R(40:42)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(50:52) + Ey(5)*R(55:57) - Ey(6)*R(59:61))
S(14:15) = Ez(2)*Ex(1)*(Ey(4)*R(35:36) - Ey(5)*R(40:41) + Ey(6)*R(44:45)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(55:56) + Ey(5)*R(59:60) - Ey(6)*R(62:63))
S(16) = Ez(2)*Ex(1)*(Ey(4)*R(40) - Ey(5)*R(44) + Ey(6)*R(47)) + Ez(3)*Ex(1)*(-Ey(4) &
   *R(59) + Ey(5)*R(62) - Ey(6)*R(64))
S(17:18) = Ez(2)*Ex(1)*(Ey(4)*R(50:51) - Ey(5)*R(55:56) + Ey(6)*R(59:60)) + Ez(3) &
   *Ex(1)*(-Ey(4)*R(65:66) + Ey(5)*R(69:70) - Ey(6)*R(72:73))
S(19) = Ez(2)*Ex(1)*(Ey(4)*R(55) - Ey(5)*R(59) + Ey(6)*R(62)) + Ez(3)*Ex(1)*(-Ey(4) &
   *R(69) + Ey(5)*R(72) - Ey(6)*R(74))
S(20) = Ez(2)*Ex(1)*(Ey(4)*R(65) - Ey(5)*R(69) + Ey(6)*R(72)) + Ez(3)*Ex(1)*(-Ey(4) &
   *R(75) + Ey(5)*R(78) - Ey(6)*R(80))
end subroutine auto2e_KetTransform_3_0_2_1

subroutine auto2e_KetTransform_3_0_1_2(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(0,1,2)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(4)*Ex(1)*(Ey(2)*R(1:4) - Ey(3)*R(8:11)) + Ez(5)*Ex(1)*(-Ey(2)*R(29:32)  &
   + Ey(3)*R(35:38)) + Ex(1)*Ez(6)*(Ey(2)*R(50:53) - Ey(3)*R(55:58))
S(5:7) = Ez(4)*Ex(1)*(Ey(2)*R(8:10) - Ey(3)*R(14:16)) + Ez(5)*Ex(1)*(-Ey(2)*R(35:37)  &
   + Ey(3)*R(40:42)) + Ex(1)*Ez(6)*(Ey(2)*R(55:57) - Ey(3)*R(59:61))
S(8:9) = Ez(4)*Ex(1)*(Ey(2)*R(14:15) - Ey(3)*R(19:20)) + Ez(5)*Ex(1)*(-Ey(2)*R(40:41)  &
   + Ey(3)*R(44:45)) + Ex(1)*Ez(6)*(Ey(2)*R(59:60) - Ey(3)*R(62:63))
S(10) = Ez(4)*Ex(1)*(Ey(2)*R(19) - Ey(3)*R(23)) + Ez(5)*Ex(1)*(-Ey(2)*R(44) + Ey(3) &
   *R(47)) + Ex(1)*Ez(6)*(Ey(2)*R(62) - Ey(3)*R(64))
S(11:13) = Ez(4)*Ex(1)*(Ey(2)*R(29:31) - Ey(3)*R(35:37)) + Ez(5)*Ex(1)*(-Ey(2)*R(50:52)  &
   + Ey(3)*R(55:57)) + Ex(1)*Ez(6)*(Ey(2)*R(65:67) - Ey(3)*R(69:71))
S(14:15) = Ez(4)*Ex(1)*(Ey(2)*R(35:36) - Ey(3)*R(40:41)) + Ez(5)*Ex(1)*(-Ey(2)*R(55:56)  &
   + Ey(3)*R(59:60)) + Ex(1)*Ez(6)*(Ey(2)*R(69:70) - Ey(3)*R(72:73))
S(16) = Ez(4)*Ex(1)*(Ey(2)*R(40) - Ey(3)*R(44)) + Ez(5)*Ex(1)*(-Ey(2)*R(59) + Ey(3) &
   *R(62)) + Ex(1)*Ez(6)*(Ey(2)*R(72) - Ey(3)*R(74))
S(17:18) = Ez(4)*Ex(1)*(Ey(2)*R(50:51) - Ey(3)*R(55:56)) + Ez(5)*Ex(1)*(-Ey(2)*R(65:66)  &
   + Ey(3)*R(69:70)) + Ex(1)*Ez(6)*(Ey(2)*R(75:76) - Ey(3)*R(78:79))
S(19) = Ez(4)*Ex(1)*(Ey(2)*R(55) - Ey(3)*R(59)) + Ez(5)*Ex(1)*(-Ey(2)*R(69) + Ey(3) &
   *R(72)) + Ex(1)*Ez(6)*(Ey(2)*R(78) - Ey(3)*R(80))
S(20) = Ez(4)*Ex(1)*(Ey(2)*R(65) - Ey(3)*R(69)) + Ez(5)*Ex(1)*(-Ey(2)*R(75) + Ey(3) &
   *R(78)) + Ex(1)*Ez(6)*(Ey(2)*R(81) - Ey(3)*R(83))
end subroutine auto2e_KetTransform_3_0_1_2

subroutine auto2e_KetTransform_3_0_0_3(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(0,0,3)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ex(1)*Ey(1)*(Ez(7)*R(1:4) - Ez(8)*R(29:32) + Ez(9)*R(50:53) - Ez(10)*R(65:68))
S(5:7) = Ex(1)*Ey(1)*(Ez(7)*R(8:10) - Ez(8)*R(35:37) + Ez(9)*R(55:57) - Ez(10)*R(69:71))
S(8:9) = Ex(1)*Ey(1)*(Ez(7)*R(14:15) - Ez(8)*R(40:41) + Ez(9)*R(59:60) - Ez(10)*R(72:73))
S(10) = Ex(1)*Ey(1)*(Ez(7)*R(19) - Ez(8)*R(44) + Ez(9)*R(62) - Ez(10)*R(74))
S(11:13) = Ex(1)*Ey(1)*(Ez(7)*R(29:31) - Ez(8)*R(50:52) + Ez(9)*R(65:67) - Ez(10) &
   *R(75:77))
S(14:15) = Ex(1)*Ey(1)*(Ez(7)*R(35:36) - Ez(8)*R(55:56) + Ez(9)*R(69:70) - Ez(10) &
   *R(78:79))
S(16) = Ex(1)*Ey(1)*(Ez(7)*R(40) - Ez(8)*R(59) + Ez(9)*R(72) - Ez(10)*R(80))
S(17:18) = Ex(1)*Ey(1)*(Ez(7)*R(50:51) - Ez(8)*R(65:66) + Ez(9)*R(75:76) - Ez(10) &
   *R(81:82))
S(19) = Ex(1)*Ey(1)*(Ez(7)*R(55) - Ez(8)*R(69) + Ez(9)*R(78) - Ez(10)*R(83))
S(20) = Ex(1)*Ey(1)*(Ez(7)*R(65) - Ez(8)*R(75) + Ez(9)*R(81) - Ez(10)*R(84))
end subroutine auto2e_KetTransform_3_0_0_3
end module auto2e_KetTransform_3_3