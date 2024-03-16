module auto2e_KetTransform_3_2
use arithmetic
use math_constants
implicit none
contains

subroutine auto2e_KetTransform_3_2_0_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(2,0,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(1)*Ey(1)*(Ex(4)*R(1:4) - Ex(5)*R(2:5) + Ex(6)*R(3:6))
S(5:7) = Ez(1)*Ey(1)*(Ex(4)*R(7:9) - Ex(5)*R(8:10) + Ex(6)*R(9:11))
S(8:9) = Ez(1)*Ey(1)*(Ex(4)*R(12:13) - Ex(5)*R(13:14) + Ex(6)*R(14:15))
S(10) = Ez(1)*Ey(1)*(Ex(4)*R(16) - Ex(5)*R(17) + Ex(6)*R(18))
S(11:13) = Ez(1)*Ey(1)*(Ex(4)*R(22:24) - Ex(5)*R(23:25) + Ex(6)*R(24:26))
S(14:15) = Ez(1)*Ey(1)*(Ex(4)*R(27:28) - Ex(5)*R(28:29) + Ex(6)*R(29:30))
S(16) = Ez(1)*Ey(1)*(Ex(4)*R(31) - Ex(5)*R(32) + Ex(6)*R(33))
S(17:18) = Ez(1)*Ey(1)*(Ex(4)*R(37:38) - Ex(5)*R(38:39) + Ex(6)*R(39:40))
S(19) = Ez(1)*Ey(1)*(Ex(4)*R(41) - Ex(5)*R(42) + Ex(6)*R(43))
S(20) = Ez(1)*Ey(1)*(Ex(4)*R(47) - Ex(5)*R(48) + Ex(6)*R(49))
end subroutine auto2e_KetTransform_3_2_0_0

subroutine auto2e_KetTransform_3_1_1_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(1,1,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(1)*Ey(2)*(Ex(2)*R(1:4) - Ex(3)*R(2:5)) + Ez(1)*Ey(3)*(-Ex(2)*R(7:10)  &
   + Ex(3)*R(8:11))
S(5:7) = Ez(1)*Ey(2)*(Ex(2)*R(7:9) - Ex(3)*R(8:10)) + Ez(1)*Ey(3)*(-Ex(2)*R(12:14)  &
   + Ex(3)*R(13:15))
S(8:9) = Ez(1)*Ey(2)*(Ex(2)*R(12:13) - Ex(3)*R(13:14)) + Ez(1)*Ey(3)*(-Ex(2)*R(16:17)  &
   + Ex(3)*R(17:18))
S(10) = Ez(1)*Ey(2)*(Ex(2)*R(16) - Ex(3)*R(17)) + Ez(1)*Ey(3)*(-Ex(2)*R(19) + Ex(3) &
   *R(20))
S(11:13) = Ez(1)*Ey(2)*(Ex(2)*R(22:24) - Ex(3)*R(23:25)) + Ez(1)*Ey(3)*(-Ex(2)*R(27:29)  &
   + Ex(3)*R(28:30))
S(14:15) = Ez(1)*Ey(2)*(Ex(2)*R(27:28) - Ex(3)*R(28:29)) + Ez(1)*Ey(3)*(-Ex(2)*R(31:32)  &
   + Ex(3)*R(32:33))
S(16) = Ez(1)*Ey(2)*(Ex(2)*R(31) - Ex(3)*R(32)) + Ez(1)*Ey(3)*(-Ex(2)*R(34) + Ex(3) &
   *R(35))
S(17:18) = Ez(1)*Ey(2)*(Ex(2)*R(37:38) - Ex(3)*R(38:39)) + Ez(1)*Ey(3)*(-Ex(2)*R(41:42)  &
   + Ex(3)*R(42:43))
S(19) = Ez(1)*Ey(2)*(Ex(2)*R(41) - Ex(3)*R(42)) + Ez(1)*Ey(3)*(-Ex(2)*R(44) + Ex(3) &
   *R(45))
S(20) = Ez(1)*Ey(2)*(Ex(2)*R(47) - Ex(3)*R(48)) + Ez(1)*Ey(3)*(-Ex(2)*R(50) + Ex(3) &
   *R(51))
end subroutine auto2e_KetTransform_3_1_1_0

subroutine auto2e_KetTransform_3_1_0_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(1,0,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(2)*Ey(1)*(Ex(2)*R(1:4) - Ex(3)*R(2:5)) + Ez(3)*Ey(1)*(-Ex(2)*R(22:25)  &
   + Ex(3)*R(23:26))
S(5:7) = Ez(2)*Ey(1)*(Ex(2)*R(7:9) - Ex(3)*R(8:10)) + Ez(3)*Ey(1)*(-Ex(2)*R(27:29)  &
   + Ex(3)*R(28:30))
S(8:9) = Ez(2)*Ey(1)*(Ex(2)*R(12:13) - Ex(3)*R(13:14)) + Ez(3)*Ey(1)*(-Ex(2)*R(31:32)  &
   + Ex(3)*R(32:33))
S(10) = Ez(2)*Ey(1)*(Ex(2)*R(16) - Ex(3)*R(17)) + Ez(3)*Ey(1)*(-Ex(2)*R(34) + Ex(3) &
   *R(35))
S(11:13) = Ez(2)*Ey(1)*(Ex(2)*R(22:24) - Ex(3)*R(23:25)) + Ez(3)*Ey(1)*(-Ex(2)*R(37:39)  &
   + Ex(3)*R(38:40))
S(14:15) = Ez(2)*Ey(1)*(Ex(2)*R(27:28) - Ex(3)*R(28:29)) + Ez(3)*Ey(1)*(-Ex(2)*R(41:42)  &
   + Ex(3)*R(42:43))
S(16) = Ez(2)*Ey(1)*(Ex(2)*R(31) - Ex(3)*R(32)) + Ez(3)*Ey(1)*(-Ex(2)*R(44) + Ex(3) &
   *R(45))
S(17:18) = Ez(2)*Ey(1)*(Ex(2)*R(37:38) - Ex(3)*R(38:39)) + Ez(3)*Ey(1)*(-Ex(2)*R(47:48)  &
   + Ex(3)*R(48:49))
S(19) = Ez(2)*Ey(1)*(Ex(2)*R(41) - Ex(3)*R(42)) + Ez(3)*Ey(1)*(-Ex(2)*R(50) + Ex(3) &
   *R(51))
S(20) = Ez(2)*Ey(1)*(Ex(2)*R(47) - Ex(3)*R(48)) + Ez(3)*Ey(1)*(-Ex(2)*R(53) + Ex(3) &
   *R(54))
end subroutine auto2e_KetTransform_3_1_0_1

subroutine auto2e_KetTransform_3_0_2_0(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(0,2,0)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(1)*Ex(1)*(Ey(4)*R(1:4) - Ey(5)*R(7:10) + Ey(6)*R(12:15))
S(5:7) = Ez(1)*Ex(1)*(Ey(4)*R(7:9) - Ey(5)*R(12:14) + Ey(6)*R(16:18))
S(8:9) = Ez(1)*Ex(1)*(Ey(4)*R(12:13) - Ey(5)*R(16:17) + Ey(6)*R(19:20))
S(10) = Ez(1)*Ex(1)*(Ey(4)*R(16) - Ey(5)*R(19) + Ey(6)*R(21))
S(11:13) = Ez(1)*Ex(1)*(Ey(4)*R(22:24) - Ey(5)*R(27:29) + Ey(6)*R(31:33))
S(14:15) = Ez(1)*Ex(1)*(Ey(4)*R(27:28) - Ey(5)*R(31:32) + Ey(6)*R(34:35))
S(16) = Ez(1)*Ex(1)*(Ey(4)*R(31) - Ey(5)*R(34) + Ey(6)*R(36))
S(17:18) = Ez(1)*Ex(1)*(Ey(4)*R(37:38) - Ey(5)*R(41:42) + Ey(6)*R(44:45))
S(19) = Ez(1)*Ex(1)*(Ey(4)*R(41) - Ey(5)*R(44) + Ey(6)*R(46))
S(20) = Ez(1)*Ex(1)*(Ey(4)*R(47) - Ey(5)*R(50) + Ey(6)*R(52))
end subroutine auto2e_KetTransform_3_0_2_0

subroutine auto2e_KetTransform_3_0_1_1(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(0,1,1)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ez(2)*Ex(1)*(Ey(2)*R(1:4) - Ey(3)*R(7:10)) + Ez(3)*Ex(1)*(-Ey(2)*R(22:25)  &
   + Ey(3)*R(27:30))
S(5:7) = Ez(2)*Ex(1)*(Ey(2)*R(7:9) - Ey(3)*R(12:14)) + Ez(3)*Ex(1)*(-Ey(2)*R(27:29)  &
   + Ey(3)*R(31:33))
S(8:9) = Ez(2)*Ex(1)*(Ey(2)*R(12:13) - Ey(3)*R(16:17)) + Ez(3)*Ex(1)*(-Ey(2)*R(31:32)  &
   + Ey(3)*R(34:35))
S(10) = Ez(2)*Ex(1)*(Ey(2)*R(16) - Ey(3)*R(19)) + Ez(3)*Ex(1)*(-Ey(2)*R(34) + Ey(3) &
   *R(36))
S(11:13) = Ez(2)*Ex(1)*(Ey(2)*R(22:24) - Ey(3)*R(27:29)) + Ez(3)*Ex(1)*(-Ey(2)*R(37:39)  &
   + Ey(3)*R(41:43))
S(14:15) = Ez(2)*Ex(1)*(Ey(2)*R(27:28) - Ey(3)*R(31:32)) + Ez(3)*Ex(1)*(-Ey(2)*R(41:42)  &
   + Ey(3)*R(44:45))
S(16) = Ez(2)*Ex(1)*(Ey(2)*R(31) - Ey(3)*R(34)) + Ez(3)*Ex(1)*(-Ey(2)*R(44) + Ey(3) &
   *R(46))
S(17:18) = Ez(2)*Ex(1)*(Ey(2)*R(37:38) - Ey(3)*R(41:42)) + Ez(3)*Ex(1)*(-Ey(2)*R(47:48)  &
   + Ey(3)*R(50:51))
S(19) = Ez(2)*Ex(1)*(Ey(2)*R(41) - Ey(3)*R(44)) + Ez(3)*Ex(1)*(-Ey(2)*R(50) + Ey(3) &
   *R(52))
S(20) = Ez(2)*Ex(1)*(Ey(2)*R(47) - Ey(3)*R(50)) + Ez(3)*Ex(1)*(-Ey(2)*R(53) + Ey(3) &
   *R(55))
end subroutine auto2e_KetTransform_3_0_1_1

subroutine auto2e_KetTransform_3_0_0_2(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<=3
! Class=(LS|KS), L=3, K=(0,0,2)
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
S(1:4) = Ex(1)*Ey(1)*(Ez(4)*R(1:4) - Ez(5)*R(22:25) + Ez(6)*R(37:40))
S(5:7) = Ex(1)*Ey(1)*(Ez(4)*R(7:9) - Ez(5)*R(27:29) + Ez(6)*R(41:43))
S(8:9) = Ex(1)*Ey(1)*(Ez(4)*R(12:13) - Ez(5)*R(31:32) + Ez(6)*R(44:45))
S(10) = Ex(1)*Ey(1)*(Ez(4)*R(16) - Ez(5)*R(34) + Ez(6)*R(46))
S(11:13) = Ex(1)*Ey(1)*(Ez(4)*R(22:24) - Ez(5)*R(37:39) + Ez(6)*R(47:49))
S(14:15) = Ex(1)*Ey(1)*(Ez(4)*R(27:28) - Ez(5)*R(41:42) + Ez(6)*R(50:51))
S(16) = Ex(1)*Ey(1)*(Ez(4)*R(31) - Ez(5)*R(44) + Ez(6)*R(52))
S(17:18) = Ex(1)*Ey(1)*(Ez(4)*R(37:38) - Ez(5)*R(47:48) + Ez(6)*R(53:54))
S(19) = Ex(1)*Ey(1)*(Ez(4)*R(41) - Ez(5)*R(50) + Ez(6)*R(55))
S(20) = Ex(1)*Ey(1)*(Ez(4)*R(47) - Ez(5)*R(53) + Ez(6)*R(56))
end subroutine auto2e_KetTransform_3_0_0_2
end module auto2e_KetTransform_3_2