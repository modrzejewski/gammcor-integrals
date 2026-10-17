!
! Compute a batch of two-electron repulsion integrals, (AB|CD), in the Cartesian Gaussian
! and real spherical harmonic basis.
! Angular momenta: A=6, B=6, C=2, and D=1
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
module auto2e_eri_iidp
use arithmetic
use math_constants    
use auto2e_KetTransfer
use auto2e_SpherTransf
use auto2e_WMatrix_part1
use auto2e_WMatrix_part2
use auto2e_WMatrix_part3
use auto2e_WMatrix_part4
implicit none

contains

subroutine auto2e_eri_6_6_2_1(G, Ra, CntrA, ExpA, NprimA, &
Rb, CntrB, ExpB, NprimB, Rc, CntrC, ExpC, NprimC, &
Rd, CntrD, ExpD, NprimD, Kappa)
!
! Two-electron integrals class: (ii|dp)
! L-dependent normalization constants are not applied.
! Output memory layout: CDAB
!
real(F64), dimension(784, *), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(399, 16) :: W
!
! Compute (LS|KS) integrals, L=6..12 and K=2..3
!
W = ZERO
call auto2e_W_6_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 1)
call auto2e_W_7_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 1)
call auto2e_W_8_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 1)
call auto2e_W_9_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 1)
call auto2e_W_10_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 165, 1)
call auto2e_W_11_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 231, 1)
call auto2e_W_12_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 309, 1)
call auto2e_W_6_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 7)
call auto2e_W_7_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 7)
call auto2e_W_8_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 7)
call auto2e_W_9_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 7)
call auto2e_W_10_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 165, 7)
call auto2e_W_11_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 231, 7)
call auto2e_W_12_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 309, 7)
call auto2e_BraKetTransfer_6_6_2_1(G, W, Ra, Rb, Rc, Rd)
end subroutine auto2e_eri_6_6_2_1

subroutine auto2e_eri_Spher_6_6_2_1(G, Ra, CntrA, ExpA, NprimA, &
Rb, CntrB, ExpB, NprimB, Rc, CntrC, ExpC, NprimC, &
Rd, CntrD, ExpD, NprimD, Kappa)
!
! Two-electron integrals class: (ii|dp)
! L-dependent normalization constants are not applied.
! Output memory layout: CDAB
!
real(F64), dimension(169, *), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(399, 16) :: W
!
! Compute (LS|KS) integrals, L=6..12 and K=2..3
!
W = ZERO
call auto2e_W_6_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 1)
call auto2e_W_7_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 1)
call auto2e_W_8_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 1)
call auto2e_W_9_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 1)
call auto2e_W_10_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 165, 1)
call auto2e_W_11_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 231, 1)
call auto2e_W_12_2(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 309, 1)
call auto2e_W_6_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 7)
call auto2e_W_7_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 7)
call auto2e_W_8_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 7)
call auto2e_W_9_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 7)
call auto2e_W_10_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 165, 7)
call auto2e_W_11_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 231, 7)
call auto2e_W_12_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 309, 7)
call auto2e_BraKetTransfer_Spher_6_6_2_1(G, W, Ra, Rb, Rc, Rd)
end subroutine auto2e_eri_Spher_6_6_2_1

subroutine auto2e_BraKetTransfer_6_6_2_1(G, W, Ra, Rb, Rc, Rd)
!
! Transform (LS|KS), L=6..12, K=2..3 using the horizontal recurrence relation
! of Eq. 18 in Head-Gordon, M. and Pople, J.A. J. Chem. Phys. 89, 5777 (1988);
! doi: 10.1063/1.455553
!
! The output integrals are expressed in the Cartesian GTOs basis.
!
real(F64), dimension(18, *), intent(out) :: G
real(F64), dimension(399, 16), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(3) :: AB, CD
real(F64), dimension(399, 18) :: T
real(F64), dimension(18, 399) :: U
AB = Ra - Rb
CD = Rc - Rd
call auto2e_KetTransfer_2_1(T, W, CD(1), CD(2), CD(3))
U = transpose(T)
call auto2e_KetTransfer_6_6(G(:, 1:784), U, AB(1), AB(2), AB(3))
end subroutine auto2e_BraKetTransfer_6_6_2_1

subroutine auto2e_BraKetTransfer_Spher_6_6_2_1(G, W, Ra, Rb, Rc, Rd)
!
! Transform (LS|KS), L=6..12, K=2..3 using the horizontal recurrence relation
! of Eq. 18 in Head-Gordon, M. and Pople, J.A. J. Chem. Phys. 89, 5777 (1988);
! doi: 10.1063/1.455553
!
! The output integrals are expressed in the real spherical harmonics basis.
!
real(F64), dimension(15, *), intent(out) :: G
real(F64), dimension(399, 16), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(3) :: AB, CD
real(F64), dimension(399, 18) :: Txyz
real(F64), dimension(399, 15) :: T
real(F64), dimension(15, 399) :: U
real(F64), dimension(15, 784) :: Gxyz
AB = Ra - Rb
CD = Rc - Rd
call auto2e_KetTransfer_2_1(Txyz, W, CD(1), CD(2), CD(3))
call auto2e_SpherTransf_U_Vector_2_1(T, Txyz, 399)
U = transpose(T)
call auto2e_KetTransfer_6_6(Gxyz, U, AB(1), AB(2), AB(3))
call auto2e_SpherTransf_U_Vector_6_6(G(:, 1:169), Gxyz, 15)
end subroutine auto2e_BraKetTransfer_Spher_6_6_2_1

subroutine auto2e_Normalize_6_6_2_1_BACD(Gbacd, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 28
do a = 1, 28
do d = 1, 3
do c = 1, 6
w = b + (a-1)*28 + (c-1)*784 + (d-1)*4704
Gbacd(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_6_6_2_1_BACD

subroutine auto2e_Normalize_6_6_2_1_BADC(Gbadc, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 28
do a = 1, 28
do d = 1, 3
do c = 1, 6
w = b + (a-1)*28 + (d-1)*784 + (c-1)*2352
Gbadc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_6_6_2_1_BADC

subroutine auto2e_Normalize_6_6_2_1_CDBA(Gcdba, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 28
do a = 1, 28
do d = 1, 3
do c = 1, 6
w = c + (d-1)*6 + (b-1)*18 + (a-1)*504
Gcdba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_6_6_2_1_CDBA

subroutine auto2e_Normalize_6_6_2_1_DCBA(Gdcba, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 28
do a = 1, 28
do d = 1, 3
do c = 1, 6
w = d + (c-1)*3 + (b-1)*18 + (a-1)*504
Gdcba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_6_6_2_1_DCBA

subroutine auto2e_Normalize_Spher_6_6_2_1_BACD(Gbacd, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 13
do a = 1, 13
do d = 1, 3
do c = 1, 5
w = b + (a-1)*13 + (c-1)*169 + (d-1)*845
Gbacd(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_6_6_2_1_BACD

subroutine auto2e_Normalize_Spher_6_6_2_1_BADC(Gbadc, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 13
do a = 1, 13
do d = 1, 3
do c = 1, 5
w = b + (a-1)*13 + (d-1)*169 + (c-1)*507
Gbadc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_6_6_2_1_BADC

subroutine auto2e_Normalize_Spher_6_6_2_1_CDBA(Gcdba, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 13
do a = 1, 13
do d = 1, 3
do c = 1, 5
w = c + (d-1)*5 + (b-1)*15 + (a-1)*195
Gcdba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_6_6_2_1_CDBA

subroutine auto2e_Normalize_Spher_6_6_2_1_DCBA(Gdcba, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 13
do a = 1, 13
do d = 1, 3
do c = 1, 5
w = d + (c-1)*3 + (b-1)*15 + (a-1)*195
Gdcba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_6_6_2_1_DCBA

subroutine auto2e_frontend_1_2_6_6(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(14112) :: H
call auto2e_eri_6_6_2_1(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, Kappa)
call auto2e_Normalize_6_6_2_1_BACD(G, H, NormC, NormD, NormB, NormA)
end subroutine auto2e_frontend_1_2_6_6

subroutine auto2e_frontend_Spher_1_2_6_6(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(2535) :: H
call auto2e_eri_Spher_6_6_2_1(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, Kappa)
call auto2e_Normalize_Spher_6_6_2_1_BACD(G, H, NormC, NormD, NormB, NormA)
end subroutine auto2e_frontend_Spher_1_2_6_6

subroutine auto2e_frontend_2_1_6_6(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(14112) :: H
call auto2e_eri_6_6_2_1(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_6_6_2_1_BADC(G, H, NormC, NormD, NormA, NormB)
end subroutine auto2e_frontend_2_1_6_6

subroutine auto2e_frontend_Spher_2_1_6_6(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(2535) :: H
call auto2e_eri_Spher_6_6_2_1(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_Spher_6_6_2_1_BADC(G, H, NormC, NormD, NormA, NormB)
end subroutine auto2e_frontend_Spher_2_1_6_6

subroutine auto2e_frontend_6_6_1_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(14112) :: H
call auto2e_eri_6_6_2_1(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, Kappa)
call auto2e_Normalize_6_6_2_1_CDBA(G, H, NormA, NormB, NormD, NormC)
end subroutine auto2e_frontend_6_6_1_2

subroutine auto2e_frontend_Spher_6_6_1_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(2535) :: H
call auto2e_eri_Spher_6_6_2_1(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, Kappa)
call auto2e_Normalize_Spher_6_6_2_1_CDBA(G, H, NormA, NormB, NormD, NormC)
end subroutine auto2e_frontend_Spher_6_6_1_2

subroutine auto2e_frontend_6_6_2_1(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(14112) :: H
call auto2e_eri_6_6_2_1(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_6_6_2_1_DCBA(G, H, NormA, NormB, NormC, NormD)
end subroutine auto2e_frontend_6_6_2_1

subroutine auto2e_frontend_Spher_6_6_2_1(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(2535) :: H
call auto2e_eri_Spher_6_6_2_1(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_Spher_6_6_2_1_DCBA(G, H, NormA, NormB, NormC, NormD)
end subroutine auto2e_frontend_Spher_6_6_2_1
end module auto2e_eri_iidp