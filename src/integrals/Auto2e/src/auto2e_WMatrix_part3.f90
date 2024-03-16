module auto2e_WMatrix_part3
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

subroutine auto2e_W_9_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(220) :: S
real(F64), dimension(16) :: F
real(F64), dimension(55) :: ExAB, EyAB, EzAB
real(F64), dimension(28) :: ExCD, EyCD, EzCD
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
call auto2e_EtCoeffs_9(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_9(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_9(EzAB, P, Rpa(3), ONE)
call fm(15, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_15(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 6, 0, -1
do ly = 6-lx, 0, -1
lz = 6 - lx - ly
call auto2e_RCopy_KetTransform_9_6(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_9(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_9_6

subroutine auto2e_W_6_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(55, 28) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_9_6(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 55
do i = 1, 28
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_6_9

subroutine auto2e_W_9_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(220) :: S
real(F64), dimension(17) :: F
real(F64), dimension(55) :: ExAB, EyAB, EzAB
real(F64), dimension(36) :: ExCD, EyCD, EzCD
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
call auto2e_EtCoeffs_9(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_9(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_9(EzAB, P, Rpa(3), ONE)
call fm(16, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_16(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 7, 0, -1
do ly = 7-lx, 0, -1
lz = 7 - lx - ly
call auto2e_RCopy_KetTransform_9_7(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_9(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_9_7

subroutine auto2e_W_7_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(55, 36) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_9_7(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 55
do i = 1, 36
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_7_9

subroutine auto2e_W_9_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(1140) :: R
real(F64), dimension(220) :: S
real(F64), dimension(18) :: F
real(F64), dimension(55) :: ExAB, EyAB, EzAB
real(F64), dimension(45) :: ExCD, EyCD, EzCD
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
call auto2e_EtCoeffs_9(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_9(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_9(EzAB, P, Rpa(3), ONE)
call fm(17, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_17(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 8, 0, -1
do ly = 8-lx, 0, -1
lz = 8 - lx - ly
call auto2e_RCopy_KetTransform_9_8(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_9(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_9_8

subroutine auto2e_W_8_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(55, 45) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_9_8(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 55
do i = 1, 45
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_8_9

subroutine auto2e_W_9_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(1330) :: R
real(F64), dimension(220) :: S
real(F64), dimension(19) :: F
real(F64), dimension(55) :: ExAB, EyAB, EzAB
real(F64), dimension(55) :: ExCD, EyCD, EzCD
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
call auto2e_EtCoeffs_9(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_9(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_9(EzCD, Q, Rqc(3), ONE)
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
call fm(18, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_18(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 9, 0, -1
do ly = 9-lx, 0, -1
lz = 9 - lx - ly
call auto2e_RCopy_KetTransform_9_9(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_9(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_9_9

subroutine auto2e_W_10_0(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(286) :: R
real(F64), dimension(286) :: S
real(F64), dimension(11) :: F
real(F64), dimension(66) :: ExAB, EyAB, EzAB
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
call auto2e_EtCoeffs_10(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_10(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_10(EzAB, P, Rpa(3), ONE)
call fm(10, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_10(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
S = ExCD(1)*EyCD(1)*EzCD(1)*R
call auto2e_BraTransform_10(W(Row0:, Col0), S, ExAB, EyAB, EzAB)
end do
end do
end do
end do
end subroutine auto2e_W_10_0

subroutine auto2e_W_0_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(66, 1) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_10_0(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 66
do i = 1, 1
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_0_10

subroutine auto2e_W_10_1(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(286) :: S
real(F64), dimension(12) :: F
real(F64), dimension(66) :: ExAB, EyAB, EzAB
real(F64), dimension(3) :: ExCD, EyCD, EzCD
real(F64), dimension(286) :: T
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
call auto2e_EtCoeffs_10(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_10(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_10(EzAB, P, Rpa(3), ONE)
call fm(11, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_11(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 1, 0, -1
do ly = 1-lx, 0, -1
lz = 1 - lx - ly
call auto2e_RCopy_KetTransform_10_1(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_10(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_10_1

subroutine auto2e_W_1_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(66, 3) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_10_1(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 66
do i = 1, 3
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_1_10

subroutine auto2e_W_10_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(286) :: S
real(F64), dimension(13) :: F
real(F64), dimension(66) :: ExAB, EyAB, EzAB
real(F64), dimension(6) :: ExCD, EyCD, EzCD
real(F64), dimension(286) :: T
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
call auto2e_EtCoeffs_10(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_10(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_10(EzAB, P, Rpa(3), ONE)
call fm(12, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_12(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 2, 0, -1
do ly = 2-lx, 0, -1
lz = 2 - lx - ly
call auto2e_RCopy_KetTransform_10_2(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_10(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_10_2

subroutine auto2e_W_2_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(66, 6) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_10_2(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 66
do i = 1, 6
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_2_10

subroutine auto2e_W_10_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(286) :: S
real(F64), dimension(14) :: F
real(F64), dimension(66) :: ExAB, EyAB, EzAB
real(F64), dimension(10) :: ExCD, EyCD, EzCD
real(F64), dimension(286) :: T
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
call auto2e_EtCoeffs_10(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_10(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_10(EzAB, P, Rpa(3), ONE)
call fm(13, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_13(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 3, 0, -1
do ly = 3-lx, 0, -1
lz = 3 - lx - ly
call auto2e_RCopy_KetTransform_10_3(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_10(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_10_3

subroutine auto2e_W_3_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(66, 10) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_10_3(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 66
do i = 1, 10
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_3_10

subroutine auto2e_W_10_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(286) :: S
real(F64), dimension(15) :: F
real(F64), dimension(66) :: ExAB, EyAB, EzAB
real(F64), dimension(15) :: ExCD, EyCD, EzCD
real(F64), dimension(286) :: T
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
call auto2e_EtCoeffs_10(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_10(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_10(EzAB, P, Rpa(3), ONE)
call fm(14, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_14(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 4, 0, -1
do ly = 4-lx, 0, -1
lz = 4 - lx - ly
call auto2e_RCopy_KetTransform_10_4(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_10(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_10_4

subroutine auto2e_W_4_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(66, 15) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_10_4(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 66
do i = 1, 15
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_4_10

subroutine auto2e_W_10_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(286) :: S
real(F64), dimension(16) :: F
real(F64), dimension(66) :: ExAB, EyAB, EzAB
real(F64), dimension(21) :: ExCD, EyCD, EzCD
real(F64), dimension(286) :: T
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
call auto2e_EtCoeffs_10(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_10(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_10(EzAB, P, Rpa(3), ONE)
call fm(15, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_15(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 5, 0, -1
do ly = 5-lx, 0, -1
lz = 5 - lx - ly
call auto2e_RCopy_KetTransform_10_5(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_10(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_10_5

subroutine auto2e_W_5_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(66, 21) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_10_5(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 66
do i = 1, 21
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_5_10

subroutine auto2e_W_10_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(286) :: S
real(F64), dimension(17) :: F
real(F64), dimension(66) :: ExAB, EyAB, EzAB
real(F64), dimension(28) :: ExCD, EyCD, EzCD
real(F64), dimension(286) :: T
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
call auto2e_EtCoeffs_10(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_10(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_10(EzAB, P, Rpa(3), ONE)
call fm(16, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_16(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 6, 0, -1
do ly = 6-lx, 0, -1
lz = 6 - lx - ly
call auto2e_RCopy_KetTransform_10_6(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_10(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_10_6

subroutine auto2e_W_6_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(66, 28) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_10_6(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 66
do i = 1, 28
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_6_10

subroutine auto2e_W_10_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(1140) :: R
real(F64), dimension(286) :: S
real(F64), dimension(18) :: F
real(F64), dimension(66) :: ExAB, EyAB, EzAB
real(F64), dimension(36) :: ExCD, EyCD, EzCD
real(F64), dimension(286) :: T
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
call auto2e_EtCoeffs_10(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_10(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_10(EzAB, P, Rpa(3), ONE)
call fm(17, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_17(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 7, 0, -1
do ly = 7-lx, 0, -1
lz = 7 - lx - ly
call auto2e_RCopy_KetTransform_10_7(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_10(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_10_7

subroutine auto2e_W_7_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(66, 36) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_10_7(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 66
do i = 1, 36
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_7_10

subroutine auto2e_W_10_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(1330) :: R
real(F64), dimension(286) :: S
real(F64), dimension(19) :: F
real(F64), dimension(66) :: ExAB, EyAB, EzAB
real(F64), dimension(45) :: ExCD, EyCD, EzCD
real(F64), dimension(286) :: T
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
call auto2e_EtCoeffs_10(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_10(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_10(EzAB, P, Rpa(3), ONE)
call fm(18, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_18(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 8, 0, -1
do ly = 8-lx, 0, -1
lz = 8 - lx - ly
call auto2e_RCopy_KetTransform_10_8(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_10(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_10_8

subroutine auto2e_W_8_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(66, 45) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_10_8(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 66
do i = 1, 45
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_8_10

subroutine auto2e_W_10_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(1540) :: R
real(F64), dimension(286) :: S
real(F64), dimension(20) :: F
real(F64), dimension(66) :: ExAB, EyAB, EzAB
real(F64), dimension(55) :: ExCD, EyCD, EzCD
real(F64), dimension(286) :: T
integer :: lx, ly, lz, k
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_9(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_9(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_9(EzCD, Q, Rqc(3), ONE)
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
call auto2e_EtCoeffs_10(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_10(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_10(EzAB, P, Rpa(3), ONE)
call fm(19, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_19(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 9, 0, -1
do ly = 9-lx, 0, -1
lz = 9 - lx - ly
call auto2e_RCopy_KetTransform_10_9(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_10(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_10_9

subroutine auto2e_W_9_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension(66, 55) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_10_9(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, 66
do i = 1, 55
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_9_10

subroutine auto2e_W_10_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
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
real(F64), dimension(1771) :: R
real(F64), dimension(286) :: S
real(F64), dimension(21) :: F
real(F64), dimension(66) :: ExAB, EyAB, EzAB
real(F64), dimension(66) :: ExCD, EyCD, EzCD
real(F64), dimension(286) :: T
integer :: lx, ly, lz, k
Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
Rqc = Rq - Rc
call auto2e_EtCoeffs_10(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_10(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_10(EzCD, Q, Rqc(3), ONE)
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
call auto2e_EtCoeffs_10(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_10(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_10(EzAB, P, Rpa(3), ONE)
call fm(20, Alpha * dot_product(Rpq, Rpq), F)
call auto2e_Rtuv_20(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))
k = 0
do lx = 10, 0, -1
do ly = 10-lx, 0, -1
lz = 10 - lx - ly
call auto2e_RCopy_KetTransform_10_10(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_10(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
end do
end do
end do
end do
end subroutine auto2e_W_10_10
end module auto2e_WMatrix_part3
