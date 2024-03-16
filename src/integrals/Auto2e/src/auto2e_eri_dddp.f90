!
! Compute a batch of two-electron repulsion integrals, (AB|CD), in the Cartesian Gaussian
! and real spherical harmonic basis.
! Angular momenta: A=2, B=2, C=2, and D=1
!
! The order of the angular-dependent prefactors (x-Ra(1))**lx*(y-Ra(2))**ly*(z-Ra(3))**lz
! is defined by the following loop:
! I = 1
! for lx in A..0
!     for ly in A-lx..0
!         lz = A-lx-ly
!         Index[(lx,ly,lz)] = I
!         I = I + 1
!
! The contracted Gaussian functions are defined as
!
! PhiA(lx,ly,lz) = NormA(Index[(lx,ly,lz)]) * (x-Ra(1))**lx * (y-Ra(2))**ly * (z-Ra(3))**lz * 
!     Sum(Ga=1..NprimA) CntrA(Ga) Exp(-ExpA(Ga) * ((x-Ra(1))**2 + (y-Ra(2))**2 + (z-Ra(3))**2))
!
! Kappa specifies the form of the interaction potential:
! Kappa=0 for the Coulomb interaction operator 1/r12
! Kappa=1/Omega**2 for the screened interaction potential Erf(Omega*r12)/r12
!
! Code generated automatically.
!
module auto2e_eri_dddp
use arithmetic
use math_constants
use auto2e_KetTransfer
use auto2e_SpherTransf
use auto2e_WMatrix_part1
use auto2e_WMatrix_part2
use auto2e_WMatrix_part3
implicit none

contains

subroutine auto2e_eri_2_2_2_1(G, Ra, CntrA, ExpA, NprimA, &
Rb, CntrB, ExpB, NprimB, Rc, CntrC, ExpC, NprimC, &
Rd, CntrD, ExpD, NprimD, Kappa)
!
! Two-electron integrals class: (dd|dp)
! L-dependent normalization constants are not applied.
! Output memory layout: CDAB
!
real(F64), dimension(36, *), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(31, 16) :: W
!
! Compute (LS|KS) integrals, L=2..4 and K=2..3
!
W = ZERO
call auto2e_W_2_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 1)
call auto2e_W_3_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 7, 1)
call auto2e_W_4_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 17, 1)
call auto2e_W_2_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 7)
call auto2e_W_3_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 7, 7)
call auto2e_W_4_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 17, 7)
call auto2e_BraKetTransfer_2_2_2_1(G, W, Ra, Rb, Rc, Rd)
end subroutine auto2e_eri_2_2_2_1

subroutine auto2e_eri_Spher_2_2_2_1(G, Ra, CntrA, ExpA, NprimA, &
Rb, CntrB, ExpB, NprimB, Rc, CntrC, ExpC, NprimC, &
Rd, CntrD, ExpD, NprimD, Kappa)
!
! Two-electron integrals class: (dd|dp)
! L-dependent normalization constants are not applied.
! Output memory layout: CDAB
!
real(F64), dimension(25, *), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(31, 16) :: W
!
! Compute (LS|KS) integrals, L=2..4 and K=2..3
!
W = ZERO
call auto2e_W_2_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 1)
call auto2e_W_3_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 7, 1)
call auto2e_W_4_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 17, 1)
call auto2e_W_2_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 7)
call auto2e_W_3_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 7, 7)
call auto2e_W_4_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 17, 7)
call auto2e_BraKetTransfer_Spher_2_2_2_1(G, W, Ra, Rb, Rc, Rd)
end subroutine auto2e_eri_Spher_2_2_2_1

subroutine auto2e_BraKetTransfer_2_2_2_1(G, W, Ra, Rb, Rc, Rd)
!
! Transform (LS|KS), L=2..4, K=2..3 using the horizontal recurrence relation
! of Eq. 18 in Head-Gordon, M. and Pople, J.A. J. Chem. Phys. 89, 5777 (1988);
! doi: 10.1063/1.455553
!
! The output integrals are expressed in the Cartesian GTOs basis.
!
real(F64), dimension(18, *), intent(out) :: G
real(F64), dimension(31, 16), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(3) :: AB, CD
real(F64), dimension(31, 18) :: T
real(F64), dimension(18, 31) :: U
AB = Ra - Rb
CD = Rc - Rd
call auto2e_KetTransfer_2_1(T, W, CD(1), CD(2), CD(3))
U = transpose(T)
call auto2e_KetTransfer_2_2(G(:, 1:36), U, AB(1), AB(2), AB(3))
end subroutine auto2e_BraKetTransfer_2_2_2_1

subroutine auto2e_BraKetTransfer_Spher_2_2_2_1(G, W, Ra, Rb, Rc, Rd)
!
! Transform (LS|KS), L=2..4, K=2..3 using the horizontal recurrence relation
! of Eq. 18 in Head-Gordon, M. and Pople, J.A. J. Chem. Phys. 89, 5777 (1988);
! doi: 10.1063/1.455553
!
! The output integrals are expressed in the real spherical harmonics basis.
!
real(F64), dimension(15, *), intent(out) :: G
real(F64), dimension(31, 16), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(3) :: AB, CD
real(F64), dimension(31, 18) :: Txyz
real(F64), dimension(31, 15) :: T
real(F64), dimension(15, 31) :: U
real(F64), dimension(15, 36) :: Gxyz
AB = Ra - Rb
CD = Rc - Rd
call auto2e_KetTransfer_2_1(Txyz, W, CD(1), CD(2), CD(3))
call auto2e_SpherTransf_U_Vector_2_1(T, Txyz, 31)
U = transpose(T)
call auto2e_KetTransfer_2_2(Gxyz, U, AB(1), AB(2), AB(3))
call auto2e_SpherTransf_U_Vector_2_2(G(:, 1:25), Gxyz, 15)
end subroutine auto2e_BraKetTransfer_Spher_2_2_2_1

subroutine auto2e_Normalize_2_2_2_1_BACD(Gbacd, Gcdab, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: CDAB
! Output memory layout: BACD
!
real(F64), dimension(*), intent(out) :: Gbacd
real(F64), dimension(*), intent(in) :: Gcdab
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do b = 1, 6
do a = 1, 6
do d = 1, 3
do c = 1, 6
w = b + (a-1)*6 + (c-1)*36 + (d-1)*216
Gbacd(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_2_2_2_1_BACD

subroutine auto2e_Normalize_2_2_2_1_BADC(Gbadc, Gcdab, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: CDAB
! Output memory layout: BADC
!
real(F64), dimension(*), intent(out) :: Gbadc
real(F64), dimension(*), intent(in) :: Gcdab
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do b = 1, 6
do a = 1, 6
do d = 1, 3
do c = 1, 6
w = b + (a-1)*6 + (d-1)*36 + (c-1)*108
Gbadc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_2_2_2_1_BADC

subroutine auto2e_Normalize_2_2_2_1_CDBA(Gcdba, Gcdab, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: CDAB
! Output memory layout: CDBA
!
real(F64), dimension(*), intent(out) :: Gcdba
real(F64), dimension(*), intent(in) :: Gcdab
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do b = 1, 6
do a = 1, 6
do d = 1, 3
do c = 1, 6
w = c + (d-1)*6 + (b-1)*18 + (a-1)*108
Gcdba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_2_2_2_1_CDBA

subroutine auto2e_Normalize_2_2_2_1_DCBA(Gdcba, Gcdab, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: CDAB
! Output memory layout: DCBA
!
real(F64), dimension(*), intent(out) :: Gdcba
real(F64), dimension(*), intent(in) :: Gcdab
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do b = 1, 6
do a = 1, 6
do d = 1, 3
do c = 1, 6
w = d + (c-1)*3 + (b-1)*18 + (a-1)*108
Gdcba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_2_2_2_1_DCBA

subroutine auto2e_Normalize_Spher_2_2_2_1_BACD(Gbacd, Gcdab, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: CDAB
! Output memory layout: BACD
!
real(F64), dimension(*), intent(out) :: Gbacd
real(F64), dimension(*), intent(in) :: Gcdab
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do b = 1, 5
do a = 1, 5
do d = 1, 3
do c = 1, 5
w = b + (a-1)*5 + (c-1)*25 + (d-1)*125
Gbacd(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_2_2_2_1_BACD

subroutine auto2e_Normalize_Spher_2_2_2_1_BADC(Gbadc, Gcdab, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: CDAB
! Output memory layout: BADC
!
real(F64), dimension(*), intent(out) :: Gbadc
real(F64), dimension(*), intent(in) :: Gcdab
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do b = 1, 5
do a = 1, 5
do d = 1, 3
do c = 1, 5
w = b + (a-1)*5 + (d-1)*25 + (c-1)*75
Gbadc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_2_2_2_1_BADC

subroutine auto2e_Normalize_Spher_2_2_2_1_CDBA(Gcdba, Gcdab, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: CDAB
! Output memory layout: CDBA
!
real(F64), dimension(*), intent(out) :: Gcdba
real(F64), dimension(*), intent(in) :: Gcdab
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do b = 1, 5
do a = 1, 5
do d = 1, 3
do c = 1, 5
w = c + (d-1)*5 + (b-1)*15 + (a-1)*75
Gcdba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_2_2_2_1_CDBA

subroutine auto2e_Normalize_Spher_2_2_2_1_DCBA(Gdcba, Gcdab, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: CDAB
! Output memory layout: DCBA
!
real(F64), dimension(*), intent(out) :: Gdcba
real(F64), dimension(*), intent(in) :: Gcdab
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do b = 1, 5
do a = 1, 5
do d = 1, 3
do c = 1, 5
w = d + (c-1)*3 + (b-1)*15 + (a-1)*75
Gdcba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_2_2_2_1_DCBA

subroutine auto2e_frontend_2_1_2_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(648) :: H
call auto2e_eri_2_2_2_1(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_2_2_2_1_BADC(G, H, NormC, NormD, NormA, NormB)
end subroutine auto2e_frontend_2_1_2_2

subroutine auto2e_frontend_Spher_2_1_2_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(375) :: H
call auto2e_eri_Spher_2_2_2_1(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_Spher_2_2_2_1_BADC(G, H, NormC, NormD, NormA, NormB)
end subroutine auto2e_frontend_Spher_2_1_2_2

subroutine auto2e_frontend_2_2_2_1(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(648) :: H
call auto2e_eri_2_2_2_1(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_2_2_2_1_DCBA(G, H, NormA, NormB, NormC, NormD)
end subroutine auto2e_frontend_2_2_2_1

subroutine auto2e_frontend_Spher_2_2_2_1(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(375) :: H
call auto2e_eri_Spher_2_2_2_1(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_Spher_2_2_2_1_DCBA(G, H, NormA, NormB, NormC, NormD)
end subroutine auto2e_frontend_Spher_2_2_2_1

subroutine auto2e_frontend_2_2_1_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(648) :: H
call auto2e_eri_2_2_2_1(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, Kappa)
call auto2e_Normalize_2_2_2_1_CDBA(G, H, NormA, NormB, NormD, NormC)
end subroutine auto2e_frontend_2_2_1_2

subroutine auto2e_frontend_Spher_2_2_1_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(375) :: H
call auto2e_eri_Spher_2_2_2_1(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, Kappa)
call auto2e_Normalize_Spher_2_2_2_1_CDBA(G, H, NormA, NormB, NormD, NormC)
end subroutine auto2e_frontend_Spher_2_2_1_2

subroutine auto2e_frontend_1_2_2_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(648) :: H
call auto2e_eri_2_2_2_1(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, Kappa)
call auto2e_Normalize_2_2_2_1_BACD(G, H, NormC, NormD, NormB, NormA)
end subroutine auto2e_frontend_1_2_2_2

subroutine auto2e_frontend_Spher_1_2_2_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(375) :: H
call auto2e_eri_Spher_2_2_2_1(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, Kappa)
call auto2e_Normalize_Spher_2_2_2_1_BACD(G, H, NormC, NormD, NormB, NormA)
end subroutine auto2e_frontend_Spher_1_2_2_2
end module auto2e_eri_dddp