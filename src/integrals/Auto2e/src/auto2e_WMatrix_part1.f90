!
! Compute the electron repulsion integrals (LS|KS), 0<=L,K<=10, and store them in W(Row0:, Col0:).
! Because this subroutine is used before the horizontal recurrence step, the contraction coefficients
! (CntrA, CntrB, ...) depend on the primitive Gaussian index, but do not depend on the angular function.
!
! The contracted Gaussian functions are defined as
!
! PhiA(lx,ly,lz) = (x-Ra(1))**lx * (y-Ra(2))**ly * (z-Ra(3))**lz * 
!     Sum(Ga) CntrA(Ga) Exp(-ExpA(Ga) * ((x-Ra(1))**2 + (y-Ra(2))**2 + (z-Ra(3))**2))
!
! The order of the angular-dependent prefactors is defined by the following loop:
!
! I = 1
! for lx in L..0
!     for ly in L-lx..0
!         lz = L-lx-ly
!         Index[(lx,ly,lz)] = I
!         I = I + 1
!
! Set Kappa=0 for conventional integrals with the Coulomb operator 1/r12.
! Set Kappa=1/Omega**2 for integrals with the screened interaction potential=Erf(Omega*r12)/r12.
!
! Integrals with the Erf(Omega*r12)/r12 operator
! ----------------------------------------------
! The formulas for Alpha(omega) and N(omega) (see the code) are drived from Eq. 12 in
! Gill, P.M.W, Adamson, R.D. Chem. Phys. Lett. 261, 105 (1996); doi: 10.1016/0009-2614(96)00931-1,
! where Alpha=1/P+1/Q+1/Omega**2. N(omega) has to cancel the 1/Sqrt(Alpha) factor in the Boys function,
! and that's why it depends on Omega.
!
! Code generated automatically.
!
module auto2e_WMatrix_part1
use arithmetic
use math_constants
use boys
use auto2e_Hermite
use auto2e_BraTransform
use auto2e_KetTransform_10_1
use auto2e_KetTransform_10_10
use auto2e_KetTransform_10_2
use auto2e_KetTransform_10_3
use auto2e_KetTransform_10_4
use auto2e_KetTransform_10_5
use auto2e_KetTransform_10_6
use auto2e_KetTransform_10_7
use auto2e_KetTransform_10_8
use auto2e_KetTransform_10_9
use auto2e_KetTransform_2_1
use auto2e_KetTransform_2_2
use auto2e_KetTransform_3_1
use auto2e_KetTransform_3_2
use auto2e_KetTransform_3_3
use auto2e_KetTransform_4_1
use auto2e_KetTransform_4_2
use auto2e_KetTransform_4_3
use auto2e_KetTransform_4_4
use auto2e_KetTransform_5_1
use auto2e_KetTransform_5_2
use auto2e_KetTransform_5_3
use auto2e_KetTransform_5_4
use auto2e_KetTransform_5_5
use auto2e_KetTransform_6_1
use auto2e_KetTransform_6_2
use auto2e_KetTransform_6_3
use auto2e_KetTransform_6_4
use auto2e_KetTransform_6_5_part1
use auto2e_KetTransform_6_5_part2
use auto2e_KetTransform_6_6_part1
use auto2e_KetTransform_6_6_part2
use auto2e_KetTransform_6_6_part3
use auto2e_KetTransform_7_1
use auto2e_KetTransform_7_2
use auto2e_KetTransform_7_3
use auto2e_KetTransform_7_4
use auto2e_KetTransform_7_5
use auto2e_KetTransform_7_6
use auto2e_KetTransform_7_7
use auto2e_KetTransform_8_1
use auto2e_KetTransform_8_2
use auto2e_KetTransform_8_3
use auto2e_KetTransform_8_4
use auto2e_KetTransform_8_5
use auto2e_KetTransform_8_6
use auto2e_KetTransform_8_7
use auto2e_KetTransform_8_8
use auto2e_KetTransform_9_1
use auto2e_KetTransform_9_2
use auto2e_KetTransform_9_3
use auto2e_KetTransform_9_4
use auto2e_KetTransform_9_5
use auto2e_KetTransform_9_6
use auto2e_KetTransform_9_7
use auto2e_KetTransform_9_8
use auto2e_KetTransform_9_9

implicit none
contains

subroutine auto2e_W_0_0(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(1) :: R
real(F64), dimension(1) :: S
real(F64), dimension(1) :: F
real(F64), dimension(1) :: ExAB, EyAB, EzAB
real(F64), dimension(1) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
ExCD(1) = ONE
EyCD(1) = ONE
EzCD(1) = ONE
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
ExAB(1) = E000
EyAB(1) = ONE
EzAB(1) = ONE
call f0(Alpha * dot_product(Rpq, Rpq), F)
R(1) = F(1)
S = ExCD(1)*EyCD(1)*EzCD(1)*R
W(Row0, Col0) = W(Row0, Col0) + ExAB(1)*EyAB(1)*EzAB(1)*S(1)
end do
end do
end do
end do
end subroutine auto2e_W_0_0

subroutine auto2e_W_1_0(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(4) :: R
real(F64), dimension(4) :: S
real(F64), dimension(2) :: F
real(F64), dimension(3) :: ExAB, EyAB, EzAB
real(F64), dimension(1) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
ExCD(1) = ONE
EyCD(1) = ONE
EzCD(1) = ONE
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_1(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_1(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_1(EzAB, P, Rpa(3), ONE)
call f1(Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_1(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
S = ExCD(1)*EyCD(1)*EzCD(1)*R
W(Row0+0, Col0) = W(Row0+0, Col0) + EyAB(1)*EzAB(1)*(ExAB(2)*S(1) + ExAB(3)*S(2))
W(Row0+1, Col0) = W(Row0+1, Col0) + ExAB(1)*EzAB(1)*(EyAB(2)*S(1) + EyAB(3)*S(3))
W(Row0+2, Col0) = W(Row0+2, Col0) + EyAB(1)*ExAB(1)*(EzAB(2)*S(1) + EzAB(3)*S(4))
end do
end do
end do
end do
end subroutine auto2e_W_1_0

subroutine auto2e_W_0_1(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(3, 1) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_1_0(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 3
do i = 1, 1
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_0_1

subroutine auto2e_W_1_1(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(10) :: R
real(F64), dimension(4) :: S
real(F64), dimension(3) :: F
real(F64), dimension(3) :: ExAB, EyAB, EzAB
real(F64), dimension(3) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_1(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_1(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_1(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_1(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_1(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_1(EzAB, P, Rpa(3), ONE)
call f2(Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_2(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
S(1:2) = EyCD(1)*EzCD(1)*(ExCD(2)*R(1:2) - ExCD(3)*R(2:3))
S(3) = EyCD(1)*EzCD(1)*(ExCD(2)*R(4) - ExCD(3)*R(5))
S(4) = EyCD(1)*EzCD(1)*(ExCD(2)*R(7) - ExCD(3)*R(8))
W(Row0+0, Col0+0) = W(Row0+0, Col0+0) + EyAB(1)*EzAB(1)*(ExAB(2)*S(1) + ExAB(3)*S(2))
W(Row0+1, Col0+0) = W(Row0+1, Col0+0) + ExAB(1)*EzAB(1)*(EyAB(2)*S(1) + EyAB(3)*S(3))
W(Row0+2, Col0+0) = W(Row0+2, Col0+0) + EyAB(1)*ExAB(1)*(EzAB(2)*S(1) + EzAB(3)*S(4))
S(1:2) = ExCD(1)*EzCD(1)*(EyCD(2)*R(1:2) - EyCD(3)*R(4:5))
S(3) = ExCD(1)*EzCD(1)*(EyCD(2)*R(4) - EyCD(3)*R(6))
S(4) = ExCD(1)*EzCD(1)*(EyCD(2)*R(7) - EyCD(3)*R(9))
W(Row0+0, Col0+1) = W(Row0+0, Col0+1) + EyAB(1)*EzAB(1)*(ExAB(2)*S(1) + ExAB(3)*S(2))
W(Row0+1, Col0+1) = W(Row0+1, Col0+1) + ExAB(1)*EzAB(1)*(EyAB(2)*S(1) + EyAB(3)*S(3))
W(Row0+2, Col0+1) = W(Row0+2, Col0+1) + EyAB(1)*ExAB(1)*(EzAB(2)*S(1) + EzAB(3)*S(4))
S(1:2) = EyCD(1)*ExCD(1)*(EzCD(2)*R(1:2) - EzCD(3)*R(7:8))
S(3) = EyCD(1)*ExCD(1)*(EzCD(2)*R(4) - EzCD(3)*R(9))
S(4) = EyCD(1)*ExCD(1)*(EzCD(2)*R(7) - EzCD(3)*R(10))
W(Row0+0, Col0+2) = W(Row0+0, Col0+2) + EyAB(1)*EzAB(1)*(ExAB(2)*S(1) + ExAB(3)*S(2))
W(Row0+1, Col0+2) = W(Row0+1, Col0+2) + ExAB(1)*EzAB(1)*(EyAB(2)*S(1) + EyAB(3)*S(3))
W(Row0+2, Col0+2) = W(Row0+2, Col0+2) + EyAB(1)*ExAB(1)*(EzAB(2)*S(1) + EzAB(3)*S(4))
end do
end do
end do
end do
end subroutine auto2e_W_1_1

subroutine auto2e_W_2_0(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(10) :: R
real(F64), dimension(10) :: S
real(F64), dimension(3) :: F
real(F64), dimension(6) :: ExAB, EyAB, EzAB
real(F64), dimension(1) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
ExCD(1) = ONE
EyCD(1) = ONE
EzCD(1) = ONE
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_2(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_2(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_2(EzAB, P, Rpa(3), ONE)
call f2(Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_2(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
S = ExCD(1)*EyCD(1)*EzCD(1)*R
call auto2e_BraTransform_2(W(Row0:, Col0), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_2_0

subroutine auto2e_W_0_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(6, 1) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_2_0(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 6
do i = 1, 1
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_0_2

subroutine auto2e_W_2_1(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(20) :: R
real(F64), dimension(10) :: S
real(F64), dimension(4) :: F
real(F64), dimension(6) :: ExAB, EyAB, EzAB
real(F64), dimension(3) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_1(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_1(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_1(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_2(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_2(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_2(EzAB, P, Rpa(3), ONE)
call f3(Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_3(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_2_1_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_2(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_2_0_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_2(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_2_0_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_2(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_2_1

subroutine auto2e_W_1_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(6, 3) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_2_1(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 6
do i = 1, 3
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_1_2

subroutine auto2e_W_2_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(35) :: R
real(F64), dimension(10) :: S
real(F64), dimension(5) :: F
real(F64), dimension(6) :: ExAB, EyAB, EzAB
real(F64), dimension(6) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_2(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_2(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_2(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_2(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_2(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_2(EzAB, P, Rpa(3), ONE)
call f4(Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_4(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_2_2_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_2(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_2_1_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_2(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_2_1_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_2(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_2_0_2_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_2(W(Row0:, Col0+3), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_2_0_1_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_2(W(Row0:, Col0+4), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_2_0_0_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_2(W(Row0:, Col0+5), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_2_2

subroutine auto2e_W_3_0(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(20) :: R
real(F64), dimension(20) :: S
real(F64), dimension(4) :: F
real(F64), dimension(10) :: ExAB, EyAB, EzAB
real(F64), dimension(1) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
ExCD(1) = ONE
EyCD(1) = ONE
EzCD(1) = ONE
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_3(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_3(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_3(EzAB, P, Rpa(3), ONE)
call f3(Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_3(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
S = ExCD(1)*EyCD(1)*EzCD(1)*R
call auto2e_BraTransform_3(W(Row0:, Col0), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_3_0

subroutine auto2e_W_0_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(10, 1) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_3_0(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 10
do i = 1, 1
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_0_3

subroutine auto2e_W_3_1(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(35) :: R
real(F64), dimension(20) :: S
real(F64), dimension(5) :: F
real(F64), dimension(10) :: ExAB, EyAB, EzAB
real(F64), dimension(3) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_1(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_1(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_1(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_3(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_3(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_3(EzAB, P, Rpa(3), ONE)
call f4(Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_4(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_3_1_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_0_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_0_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_3_1

subroutine auto2e_W_1_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(10, 3) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_3_1(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 10
do i = 1, 3
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_1_3

subroutine auto2e_W_3_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(56) :: R
real(F64), dimension(20) :: S
real(F64), dimension(6) :: F
real(F64), dimension(10) :: ExAB, EyAB, EzAB
real(F64), dimension(6) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_2(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_2(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_2(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_3(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_3(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_3(EzAB, P, Rpa(3), ONE)
call f5(Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_5(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_3_2_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_1_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_1_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_0_2_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+3), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_0_1_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+4), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_0_0_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+5), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_3_2

subroutine auto2e_W_2_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(10, 6) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_3_2(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 10
do i = 1, 6
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_2_3

subroutine auto2e_W_3_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(84) :: R
real(F64), dimension(20) :: S
real(F64), dimension(7) :: F
real(F64), dimension(10) :: ExAB, EyAB, EzAB
real(F64), dimension(10) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_3(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_3(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_3(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_3(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_3(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_3(EzAB, P, Rpa(3), ONE)
call f6(Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_6(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_3_3_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_2_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_2_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_1_2_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+3), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_1_1_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+4), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_1_0_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+5), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_0_3_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+6), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_0_2_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+7), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_0_1_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+8), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_3_0_0_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_3(W(Row0:, Col0+9), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_3_3

subroutine auto2e_W_4_0(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(35) :: R
real(F64), dimension(35) :: S
real(F64), dimension(5) :: F
real(F64), dimension(15) :: ExAB, EyAB, EzAB
real(F64), dimension(1) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
ExCD(1) = ONE
EyCD(1) = ONE
EzCD(1) = ONE
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_4(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_4(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_4(EzAB, P, Rpa(3), ONE)
call f4(Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_4(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
S = ExCD(1)*EyCD(1)*EzCD(1)*R
call auto2e_BraTransform_4(W(Row0:, Col0), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_4_0

subroutine auto2e_W_0_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(15, 1) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_4_0(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 15
do i = 1, 1
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_0_4

subroutine auto2e_W_4_1(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(56) :: R
real(F64), dimension(35) :: S
real(F64), dimension(6) :: F
real(F64), dimension(15) :: ExAB, EyAB, EzAB
real(F64), dimension(3) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_1(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_1(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_1(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_4(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_4(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_4(EzAB, P, Rpa(3), ONE)
call f5(Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_5(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_4_1_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_0_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_0_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_4_1

subroutine auto2e_W_1_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(15, 3) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_4_1(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 15
do i = 1, 3
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_1_4

subroutine auto2e_W_4_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(84) :: R
real(F64), dimension(35) :: S
real(F64), dimension(7) :: F
real(F64), dimension(15) :: ExAB, EyAB, EzAB
real(F64), dimension(6) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_2(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_2(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_2(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_4(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_4(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_4(EzAB, P, Rpa(3), ONE)
call f6(Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_6(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_4_2_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_1_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_1_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_0_2_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+3), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_0_1_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+4), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_0_0_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+5), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_4_2

subroutine auto2e_W_2_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(15, 6) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_4_2(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 15
do i = 1, 6
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_2_4

subroutine auto2e_W_4_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(120) :: R
real(F64), dimension(35) :: S
real(F64), dimension(8) :: F
real(F64), dimension(15) :: ExAB, EyAB, EzAB
real(F64), dimension(10) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_3(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_3(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_3(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_4(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_4(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_4(EzAB, P, Rpa(3), ONE)
call fm(7, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_7(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_4_3_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_2_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_2_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_1_2_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+3), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_1_1_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+4), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_1_0_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+5), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_0_3_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+6), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_0_2_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+7), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_0_1_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+8), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_0_0_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+9), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_4_3

subroutine auto2e_W_3_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(15, 10) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_4_3(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 15
do i = 1, 10
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_3_4

subroutine auto2e_W_4_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(165) :: R
real(F64), dimension(35) :: S
real(F64), dimension(9) :: F
real(F64), dimension(15) :: ExAB, EyAB, EzAB
real(F64), dimension(15) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_4(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_4(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_4(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_4(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_4(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_4(EzAB, P, Rpa(3), ONE)
call fm(8, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_8(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_4_4_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_3_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_3_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_2_2_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+3), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_2_1_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+4), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_2_0_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+5), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_1_3_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+6), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_1_2_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+7), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_1_1_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+8), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_1_0_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+9), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_0_4_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+10), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_0_3_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+11), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_0_2_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+12), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_0_1_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+13), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_4_0_0_4(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_4(W(Row0:, Col0+14), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_4_4

subroutine auto2e_W_5_0(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(56) :: R
real(F64), dimension(56) :: S
real(F64), dimension(6) :: F
real(F64), dimension(21) :: ExAB, EyAB, EzAB
real(F64), dimension(1) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
ExCD(1) = ONE
EyCD(1) = ONE
EzCD(1) = ONE
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_5(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_5(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_5(EzAB, P, Rpa(3), ONE)
call f5(Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_5(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
S = ExCD(1)*EyCD(1)*EzCD(1)*R
call auto2e_BraTransform_5(W(Row0:, Col0), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_5_0

subroutine auto2e_W_0_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(21, 1) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_5_0(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 21
do i = 1, 1
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_0_5

subroutine auto2e_W_5_1(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(84) :: R
real(F64), dimension(56) :: S
real(F64), dimension(7) :: F
real(F64), dimension(21) :: ExAB, EyAB, EzAB
real(F64), dimension(3) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_1(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_1(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_1(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_5(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_5(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_5(EzAB, P, Rpa(3), ONE)
call f6(Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_6(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_5_1_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_5_1

subroutine auto2e_W_1_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(21, 3) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_5_1(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 21
do i = 1, 3
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_1_5

subroutine auto2e_W_5_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(120) :: R
real(F64), dimension(56) :: S
real(F64), dimension(8) :: F
real(F64), dimension(21) :: ExAB, EyAB, EzAB
real(F64), dimension(6) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_2(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_2(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_2(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_5(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_5(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_5(EzAB, P, Rpa(3), ONE)
call fm(7, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_7(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_5_2_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_1_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_1_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_2_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+3), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_1_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+4), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_0_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+5), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_5_2

subroutine auto2e_W_2_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(21, 6) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_5_2(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 21
do i = 1, 6
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_2_5

subroutine auto2e_W_5_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(165) :: R
real(F64), dimension(56) :: S
real(F64), dimension(9) :: F
real(F64), dimension(21) :: ExAB, EyAB, EzAB
real(F64), dimension(10) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_3(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_3(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_3(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_5(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_5(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_5(EzAB, P, Rpa(3), ONE)
call fm(8, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_8(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_5_3_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_2_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_2_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_1_2_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+3), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_1_1_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+4), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_1_0_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+5), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_3_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+6), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_2_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+7), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_1_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+8), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_0_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+9), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_5_3

subroutine auto2e_W_3_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(21, 10) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_5_3(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 21
do i = 1, 10
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_3_5

subroutine auto2e_W_5_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(220) :: R
real(F64), dimension(56) :: S
real(F64), dimension(10) :: F
real(F64), dimension(21) :: ExAB, EyAB, EzAB
real(F64), dimension(15) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_4(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_4(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_4(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_5(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_5(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_5(EzAB, P, Rpa(3), ONE)
call fm(9, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_9(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_5_4_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_3_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_3_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_2_2_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+3), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_2_1_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+4), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_2_0_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+5), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_1_3_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+6), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_1_2_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+7), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_1_1_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+8), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_1_0_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+9), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_4_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+10), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_3_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+11), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_2_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+12), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_1_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+13), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_0_4(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+14), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_5_4

subroutine auto2e_W_4_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(21, 15) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_5_4(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 21
do i = 1, 15
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_4_5

subroutine auto2e_W_5_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(286) :: R
real(F64), dimension(56) :: S
real(F64), dimension(11) :: F
real(F64), dimension(21) :: ExAB, EyAB, EzAB
real(F64), dimension(21) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_5(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_5(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_5(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_5(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_5(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_5(EzAB, P, Rpa(3), ONE)
call fm(10, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_10(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_5_5_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_4_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_4_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_3_2_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+3), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_3_1_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+4), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_3_0_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+5), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_2_3_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+6), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_2_2_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+7), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_2_1_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+8), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_2_0_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+9), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_1_4_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+10), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_1_3_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+11), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_1_2_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+12), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_1_1_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+13), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_1_0_4(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+14), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_5_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+15), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_4_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+16), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_3_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+17), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_2_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+18), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_1_4(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+19), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_5_0_0_5(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_5(W(Row0:, Col0+20), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_5_5

subroutine auto2e_W_6_0(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(84) :: R
real(F64), dimension(84) :: S
real(F64), dimension(7) :: F
real(F64), dimension(28) :: ExAB, EyAB, EzAB
real(F64), dimension(1) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
ExCD(1) = ONE
EyCD(1) = ONE
EzCD(1) = ONE
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_6(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_6(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_6(EzAB, P, Rpa(3), ONE)
call f6(Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_6(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
S = ExCD(1)*EyCD(1)*EzCD(1)*R
call auto2e_BraTransform_6(W(Row0:, Col0), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_6_0

subroutine auto2e_W_0_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(28, 1) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_6_0(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 28
do i = 1, 1
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_0_6

subroutine auto2e_W_6_1(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(120) :: R
real(F64), dimension(84) :: S
real(F64), dimension(8) :: F
real(F64), dimension(28) :: ExAB, EyAB, EzAB
real(F64), dimension(3) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_1(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_1(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_1(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_6(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_6(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_6(EzAB, P, Rpa(3), ONE)
call fm(7, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_7(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_6_1_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_6_1

subroutine auto2e_W_1_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(28, 3) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_6_1(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 28
do i = 1, 3
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_1_6

subroutine auto2e_W_6_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(165) :: R
real(F64), dimension(84) :: S
real(F64), dimension(9) :: F
real(F64), dimension(28) :: ExAB, EyAB, EzAB
real(F64), dimension(6) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_2(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_2(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_2(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_6(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_6(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_6(EzAB, P, Rpa(3), ONE)
call fm(8, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_8(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_6_2_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_2_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+3), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_1_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+4), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_0_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+5), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_6_2

subroutine auto2e_W_2_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(28, 6) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_6_2(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 28
do i = 1, 6
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_2_6

subroutine auto2e_W_6_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(220) :: R
real(F64), dimension(84) :: S
real(F64), dimension(10) :: F
real(F64), dimension(28) :: ExAB, EyAB, EzAB
real(F64), dimension(10) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_3(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_3(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_3(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_6(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_6(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_6(EzAB, P, Rpa(3), ONE)
call fm(9, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_9(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_6_3_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_2_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_2_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_2_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+3), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_1_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+4), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_0_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+5), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_3_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+6), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_2_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+7), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_1_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+8), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_0_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+9), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_6_3

subroutine auto2e_W_3_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(28, 10) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_6_3(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 28
do i = 1, 10
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_3_6

subroutine auto2e_W_6_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, Rqc, Rpa, Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension(286) :: R
real(F64), dimension(84) :: S
real(F64), dimension(11) :: F
real(F64), dimension(28) :: ExAB, EyAB, EzAB
real(F64), dimension(15) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_4(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_4(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_4(EzCD, Q, Rqc(3), ONE)
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
Rpa = Rp - Ra
call auto2e_EtCoeffs_6(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_6(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_6(EzAB, P, Rpa(3), ONE)
call fm(10, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_10(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_6_4_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_3_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_3_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_2_2_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+3), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_2_1_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+4), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_2_0_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+5), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_3_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+6), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_2_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+7), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_1_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+8), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_0_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+9), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_4_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+10), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_3_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+11), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_2_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+12), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_1_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+13), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_0_4(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+14), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_6_4

subroutine auto2e_W_4_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(28, 15) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_6_4(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 28
do i = 1, 15
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_4_6
end module auto2e_WMatrix_part1
