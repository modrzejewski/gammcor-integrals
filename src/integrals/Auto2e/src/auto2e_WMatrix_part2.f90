module auto2e_WMatrix_part2
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

subroutine auto2e_W_6_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(364) :: R
real(F64), dimension(84) :: S
real(F64), dimension(12) :: F
real(F64), dimension(28) :: ExAB, EyAB, EzAB
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
call auto2e_EtCoeffs_6(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_6(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_6(EzAB, P, Rpa(3), ONE)
call fm(11, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_11(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_6_5_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_4_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_4_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_3_2_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+3), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_3_1_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+4), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_3_0_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+5), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_2_3_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+6), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_2_2_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+7), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_2_1_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+8), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_2_0_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+9), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_4_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+10), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_3_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+11), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_2_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+12), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_1_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+13), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_0_4(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+14), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_5_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+15), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_4_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+16), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_3_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+17), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_2_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+18), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_1_4(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+19), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_0_5(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+20), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_6_5

subroutine auto2e_W_5_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(28, 21) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_6_5(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 28
do i = 1, 21
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_5_6

subroutine auto2e_W_6_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(455) :: R
real(F64), dimension(84) :: S
real(F64), dimension(13) :: F
real(F64), dimension(28) :: ExAB, EyAB, EzAB
real(F64), dimension(28) :: ExCD, EyCD, EzCD
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_6(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_6(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_6(EzCD, Q, Rqc(3), ONE)
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
call fm(12, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_12(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
call auto2e_KetTransform_6_6_0_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+0), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_5_1_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+1), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_5_0_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+2), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_4_2_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+3), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_4_1_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+4), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_4_0_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+5), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_3_3_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+6), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_3_2_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+7), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_3_1_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+8), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_3_0_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+9), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_2_4_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+10), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_2_3_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+11), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_2_2_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+12), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_2_1_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+13), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_2_0_4(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+14), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_5_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+15), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_4_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+16), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_3_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+17), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_2_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+18), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_1_4(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+19), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_1_0_5(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+20), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_6_0(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+21), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_5_1(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+22), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_4_2(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+23), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_3_3(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+24), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_2_4(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+25), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_1_5(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+26), S, ExAB, EyAB, EzAB)
call auto2e_KetTransform_6_0_0_6(S, ExCD, EyCD, EzCD, R)
call auto2e_BraTransform_6(W(Row0:, Col0+27), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_6_6

subroutine auto2e_W_7_0(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(120) :: R
real(F64), dimension(120) :: S
real(F64), dimension(8) :: F
real(F64), dimension(36) :: ExAB, EyAB, EzAB
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
call auto2e_EtCoeffs_7(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_7(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_7(EzAB, P, Rpa(3), ONE)
call fm(7, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_7(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
S = ExCD(1)*EyCD(1)*EzCD(1)*R
call auto2e_BraTransform_7(W(Row0:, Col0), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_7_0

subroutine auto2e_W_0_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(36, 1) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_7_0(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 36
do i = 1, 1
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_0_7

subroutine auto2e_W_7_1(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(120) :: S
real(F64), dimension(9) :: F
real(F64), dimension(36) :: ExAB, EyAB, EzAB
real(F64), dimension(3) :: ExCD, EyCD, EzCD
real(F64), dimension(120) :: T
integer :: lx, ly, lz, k
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
call auto2e_EtCoeffs_7(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_7(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_7(EzAB, P, Rpa(3), ONE)
call fm(8, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_8(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 1, 0, -1
do ly = 1-lx, 0, -1
lz = 1 - lx - ly
call auto2e_RCopy_KetTransform_7_1(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_7(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_7_1

subroutine auto2e_W_1_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(36, 3) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_7_1(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 36
do i = 1, 3
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_1_7

subroutine auto2e_W_7_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(120) :: S
real(F64), dimension(10) :: F
real(F64), dimension(36) :: ExAB, EyAB, EzAB
real(F64), dimension(6) :: ExCD, EyCD, EzCD
real(F64), dimension(120) :: T
integer :: lx, ly, lz, k
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
call auto2e_EtCoeffs_7(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_7(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_7(EzAB, P, Rpa(3), ONE)
call fm(9, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_9(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 2, 0, -1
do ly = 2-lx, 0, -1
lz = 2 - lx - ly
call auto2e_RCopy_KetTransform_7_2(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_7(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_7_2

subroutine auto2e_W_2_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(36, 6) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_7_2(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 36
do i = 1, 6
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_2_7

subroutine auto2e_W_7_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(120) :: S
real(F64), dimension(11) :: F
real(F64), dimension(36) :: ExAB, EyAB, EzAB
real(F64), dimension(10) :: ExCD, EyCD, EzCD
real(F64), dimension(120) :: T
integer :: lx, ly, lz, k
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
call auto2e_EtCoeffs_7(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_7(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_7(EzAB, P, Rpa(3), ONE)
call fm(10, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_10(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 3, 0, -1
do ly = 3-lx, 0, -1
lz = 3 - lx - ly
call auto2e_RCopy_KetTransform_7_3(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_7(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_7_3

subroutine auto2e_W_3_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(36, 10) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_7_3(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 36
do i = 1, 10
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_3_7

subroutine auto2e_W_7_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(364) :: R
real(F64), dimension(120) :: S
real(F64), dimension(12) :: F
real(F64), dimension(36) :: ExAB, EyAB, EzAB
real(F64), dimension(15) :: ExCD, EyCD, EzCD
real(F64), dimension(120) :: T
integer :: lx, ly, lz, k
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
call auto2e_EtCoeffs_7(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_7(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_7(EzAB, P, Rpa(3), ONE)
call fm(11, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_11(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 4, 0, -1
do ly = 4-lx, 0, -1
lz = 4 - lx - ly
call auto2e_RCopy_KetTransform_7_4(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_7(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_7_4

subroutine auto2e_W_4_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(36, 15) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_7_4(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 36
do i = 1, 15
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_4_7

subroutine auto2e_W_7_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(455) :: R
real(F64), dimension(120) :: S
real(F64), dimension(13) :: F
real(F64), dimension(36) :: ExAB, EyAB, EzAB
real(F64), dimension(21) :: ExCD, EyCD, EzCD
real(F64), dimension(120) :: T
integer :: lx, ly, lz, k
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
call auto2e_EtCoeffs_7(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_7(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_7(EzAB, P, Rpa(3), ONE)
call fm(12, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_12(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 5, 0, -1
do ly = 5-lx, 0, -1
lz = 5 - lx - ly
call auto2e_RCopy_KetTransform_7_5(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_7(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_7_5

subroutine auto2e_W_5_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(36, 21) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_7_5(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 36
do i = 1, 21
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_5_7

subroutine auto2e_W_7_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(560) :: R
real(F64), dimension(120) :: S
real(F64), dimension(14) :: F
real(F64), dimension(36) :: ExAB, EyAB, EzAB
real(F64), dimension(28) :: ExCD, EyCD, EzCD
real(F64), dimension(120) :: T
integer :: lx, ly, lz, k
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_6(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_6(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_6(EzCD, Q, Rqc(3), ONE)
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
call auto2e_EtCoeffs_7(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_7(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_7(EzAB, P, Rpa(3), ONE)
call fm(13, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_13(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 6, 0, -1
do ly = 6-lx, 0, -1
lz = 6 - lx - ly
call auto2e_RCopy_KetTransform_7_6(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_7(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_7_6

subroutine auto2e_W_6_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(36, 28) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_7_6(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 36
do i = 1, 28
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_6_7

subroutine auto2e_W_7_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(680) :: R
real(F64), dimension(120) :: S
real(F64), dimension(15) :: F
real(F64), dimension(36) :: ExAB, EyAB, EzAB
real(F64), dimension(36) :: ExCD, EyCD, EzCD
real(F64), dimension(120) :: T
integer :: lx, ly, lz, k
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_7(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_7(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_7(EzCD, Q, Rqc(3), ONE)
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
call auto2e_EtCoeffs_7(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_7(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_7(EzAB, P, Rpa(3), ONE)
call fm(14, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_14(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 7, 0, -1
do ly = 7-lx, 0, -1
lz = 7 - lx - ly
call auto2e_RCopy_KetTransform_7_7(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_7(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_7_7

subroutine auto2e_W_8_0(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(165) :: R
real(F64), dimension(165) :: S
real(F64), dimension(9) :: F
real(F64), dimension(45) :: ExAB, EyAB, EzAB
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
call auto2e_EtCoeffs_8(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_8(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_8(EzAB, P, Rpa(3), ONE)
call fm(8, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_8(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
S = ExCD(1)*EyCD(1)*EzCD(1)*R
call auto2e_BraTransform_8(W(Row0:, Col0), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_8_0

subroutine auto2e_W_0_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(45, 1) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_8_0(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 45
do i = 1, 1
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_0_8

subroutine auto2e_W_8_1(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(165) :: S
real(F64), dimension(10) :: F
real(F64), dimension(45) :: ExAB, EyAB, EzAB
real(F64), dimension(3) :: ExCD, EyCD, EzCD
real(F64), dimension(165) :: T
integer :: lx, ly, lz, k
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
call auto2e_EtCoeffs_8(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_8(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_8(EzAB, P, Rpa(3), ONE)
call fm(9, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_9(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 1, 0, -1
do ly = 1-lx, 0, -1
lz = 1 - lx - ly
call auto2e_RCopy_KetTransform_8_1(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_8(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_8_1

subroutine auto2e_W_1_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(45, 3) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_8_1(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 45
do i = 1, 3
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_1_8

subroutine auto2e_W_8_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(165) :: S
real(F64), dimension(11) :: F
real(F64), dimension(45) :: ExAB, EyAB, EzAB
real(F64), dimension(6) :: ExCD, EyCD, EzCD
real(F64), dimension(165) :: T
integer :: lx, ly, lz, k
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
call auto2e_EtCoeffs_8(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_8(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_8(EzAB, P, Rpa(3), ONE)
call fm(10, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_10(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 2, 0, -1
do ly = 2-lx, 0, -1
lz = 2 - lx - ly
call auto2e_RCopy_KetTransform_8_2(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_8(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_8_2

subroutine auto2e_W_2_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(45, 6) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_8_2(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 45
do i = 1, 6
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_2_8

subroutine auto2e_W_8_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(364) :: R
real(F64), dimension(165) :: S
real(F64), dimension(12) :: F
real(F64), dimension(45) :: ExAB, EyAB, EzAB
real(F64), dimension(10) :: ExCD, EyCD, EzCD
real(F64), dimension(165) :: T
integer :: lx, ly, lz, k
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
call auto2e_EtCoeffs_8(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_8(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_8(EzAB, P, Rpa(3), ONE)
call fm(11, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_11(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 3, 0, -1
do ly = 3-lx, 0, -1
lz = 3 - lx - ly
call auto2e_RCopy_KetTransform_8_3(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_8(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_8_3

subroutine auto2e_W_3_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(45, 10) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_8_3(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 45
do i = 1, 10
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_3_8

subroutine auto2e_W_8_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(455) :: R
real(F64), dimension(165) :: S
real(F64), dimension(13) :: F
real(F64), dimension(45) :: ExAB, EyAB, EzAB
real(F64), dimension(15) :: ExCD, EyCD, EzCD
real(F64), dimension(165) :: T
integer :: lx, ly, lz, k
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
call auto2e_EtCoeffs_8(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_8(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_8(EzAB, P, Rpa(3), ONE)
call fm(12, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_12(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 4, 0, -1
do ly = 4-lx, 0, -1
lz = 4 - lx - ly
call auto2e_RCopy_KetTransform_8_4(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_8(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_8_4

subroutine auto2e_W_4_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(45, 15) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_8_4(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 45
do i = 1, 15
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_4_8

subroutine auto2e_W_8_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(560) :: R
real(F64), dimension(165) :: S
real(F64), dimension(14) :: F
real(F64), dimension(45) :: ExAB, EyAB, EzAB
real(F64), dimension(21) :: ExCD, EyCD, EzCD
real(F64), dimension(165) :: T
integer :: lx, ly, lz, k
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
call auto2e_EtCoeffs_8(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_8(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_8(EzAB, P, Rpa(3), ONE)
call fm(13, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_13(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 5, 0, -1
do ly = 5-lx, 0, -1
lz = 5 - lx - ly
call auto2e_RCopy_KetTransform_8_5(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_8(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_8_5

subroutine auto2e_W_5_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(45, 21) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_8_5(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 45
do i = 1, 21
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_5_8

subroutine auto2e_W_8_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(680) :: R
real(F64), dimension(165) :: S
real(F64), dimension(15) :: F
real(F64), dimension(45) :: ExAB, EyAB, EzAB
real(F64), dimension(28) :: ExCD, EyCD, EzCD
real(F64), dimension(165) :: T
integer :: lx, ly, lz, k
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_6(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_6(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_6(EzCD, Q, Rqc(3), ONE)
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
call auto2e_EtCoeffs_8(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_8(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_8(EzAB, P, Rpa(3), ONE)
call fm(14, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_14(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 6, 0, -1
do ly = 6-lx, 0, -1
lz = 6 - lx - ly
call auto2e_RCopy_KetTransform_8_6(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_8(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_8_6

subroutine auto2e_W_6_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(45, 28) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_8_6(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 45
do i = 1, 28
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_6_8

subroutine auto2e_W_8_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(816) :: R
real(F64), dimension(165) :: S
real(F64), dimension(16) :: F
real(F64), dimension(45) :: ExAB, EyAB, EzAB
real(F64), dimension(36) :: ExCD, EyCD, EzCD
real(F64), dimension(165) :: T
integer :: lx, ly, lz, k
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_7(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_7(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_7(EzCD, Q, Rqc(3), ONE)
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
call auto2e_EtCoeffs_8(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_8(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_8(EzAB, P, Rpa(3), ONE)
call fm(15, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_15(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 7, 0, -1
do ly = 7-lx, 0, -1
lz = 7 - lx - ly
call auto2e_RCopy_KetTransform_8_7(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_8(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_8_7

subroutine auto2e_W_7_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(45, 36) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_8_7(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 45
do i = 1, 36
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_7_8

subroutine auto2e_W_8_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(969) :: R
real(F64), dimension(165) :: S
real(F64), dimension(17) :: F
real(F64), dimension(45) :: ExAB, EyAB, EzAB
real(F64), dimension(45) :: ExCD, EyCD, EzCD
real(F64), dimension(165) :: T
integer :: lx, ly, lz, k
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_8(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_8(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_8(EzCD, Q, Rqc(3), ONE)
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
call auto2e_EtCoeffs_8(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_8(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_8(EzAB, P, Rpa(3), ONE)
call fm(16, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_16(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 8, 0, -1
do ly = 8-lx, 0, -1
lz = 8 - lx - ly
call auto2e_RCopy_KetTransform_8_8(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_8(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_8_8

subroutine auto2e_W_9_0(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(220) :: R
real(F64), dimension(220) :: S
real(F64), dimension(10) :: F
real(F64), dimension(55) :: ExAB, EyAB, EzAB
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
call auto2e_EtCoeffs_9(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_9(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_9(EzAB, P, Rpa(3), ONE)
call fm(9, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_9(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
S = ExCD(1)*EyCD(1)*EzCD(1)*R
call auto2e_BraTransform_9(W(Row0:, Col0), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_9_0

subroutine auto2e_W_0_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(55, 1) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_9_0(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 55
do i = 1, 1
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_0_9

subroutine auto2e_W_9_1(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(220) :: S
real(F64), dimension(11) :: F
real(F64), dimension(55) :: ExAB, EyAB, EzAB
real(F64), dimension(3) :: ExCD, EyCD, EzCD
real(F64), dimension(220) :: T
integer :: lx, ly, lz, k
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
call auto2e_EtCoeffs_9(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_9(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_9(EzAB, P, Rpa(3), ONE)
call fm(10, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_10(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 1, 0, -1
do ly = 1-lx, 0, -1
lz = 1 - lx - ly
call auto2e_RCopy_KetTransform_9_1(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_9(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_9_1

subroutine auto2e_W_1_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(55, 3) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_9_1(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 55
do i = 1, 3
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_1_9

subroutine auto2e_W_9_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(364) :: R
real(F64), dimension(220) :: S
real(F64), dimension(12) :: F
real(F64), dimension(55) :: ExAB, EyAB, EzAB
real(F64), dimension(6) :: ExCD, EyCD, EzCD
real(F64), dimension(220) :: T
integer :: lx, ly, lz, k
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
call auto2e_EtCoeffs_9(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_9(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_9(EzAB, P, Rpa(3), ONE)
call fm(11, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_11(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 2, 0, -1
do ly = 2-lx, 0, -1
lz = 2 - lx - ly
call auto2e_RCopy_KetTransform_9_2(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_9(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_9_2

subroutine auto2e_W_2_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(55, 6) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_9_2(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 55
do i = 1, 6
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_2_9

subroutine auto2e_W_9_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(455) :: R
real(F64), dimension(220) :: S
real(F64), dimension(13) :: F
real(F64), dimension(55) :: ExAB, EyAB, EzAB
real(F64), dimension(10) :: ExCD, EyCD, EzCD
real(F64), dimension(220) :: T
integer :: lx, ly, lz, k
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
call auto2e_EtCoeffs_9(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_9(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_9(EzAB, P, Rpa(3), ONE)
call fm(12, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_12(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 3, 0, -1
do ly = 3-lx, 0, -1
lz = 3 - lx - ly
call auto2e_RCopy_KetTransform_9_3(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_9(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_9_3

subroutine auto2e_W_3_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(55, 10) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_9_3(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 55
do i = 1, 10
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_3_9

subroutine auto2e_W_9_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(560) :: R
real(F64), dimension(220) :: S
real(F64), dimension(14) :: F
real(F64), dimension(55) :: ExAB, EyAB, EzAB
real(F64), dimension(15) :: ExCD, EyCD, EzCD
real(F64), dimension(220) :: T
integer :: lx, ly, lz, k
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
call auto2e_EtCoeffs_9(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_9(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_9(EzAB, P, Rpa(3), ONE)
call fm(13, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_13(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 4, 0, -1
do ly = 4-lx, 0, -1
lz = 4 - lx - ly
call auto2e_RCopy_KetTransform_9_4(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_9(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_9_4

subroutine auto2e_W_4_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(55, 15) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_9_4(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 55
do i = 1, 15
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_4_9

subroutine auto2e_W_9_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(680) :: R
real(F64), dimension(220) :: S
real(F64), dimension(15) :: F
real(F64), dimension(55) :: ExAB, EyAB, EzAB
real(F64), dimension(21) :: ExCD, EyCD, EzCD
real(F64), dimension(220) :: T
integer :: lx, ly, lz, k
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
call auto2e_EtCoeffs_9(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_9(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_9(EzAB, P, Rpa(3), ONE)
call fm(14, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_14(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 5, 0, -1
do ly = 5-lx, 0, -1
lz = 5 - lx - ly
call auto2e_RCopy_KetTransform_9_5(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_9(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_9_5

subroutine auto2e_W_5_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(55, 21) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_9_5(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 55
do i = 1, 21
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_5_9
end module auto2e_WMatrix_part2
