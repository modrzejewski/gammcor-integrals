!
! Compute a batch of two-electron repulsion integrals, (AB|CD), in the Cartesian Gaussian
! and real spherical harmonic basis.
! Angular momenta: A=5, B=2, C=5, and D=1
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
module auto2e_eri_hdhp
use arithmetic
use auto2e_KetTransfer
use auto2e_SpherTransf
use auto2e_WMatrix_part1
use auto2e_WMatrix_part2
use auto2e_WMatrix_part3
implicit none

contains

subroutine auto2e_eri_5_2_5_1(G, Ra, CntrA, ExpA, NprimA, &
Rb, CntrB, ExpB, NprimB, Rc, CntrC, ExpC, NprimC, &
Rd, CntrD, ExpD, NprimD, Kappa)
!
! Two-electron integrals class: (hd|hp)
! L-dependent normalization constants are not applied.
! Output memory layout: CDAB
!
real(F64), dimension(126, *), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(85, 49) :: W
!
! Compute (LS|KS) integrals, L=5..7 and K=5..6
!
W = ZERO
call auto2e_W_5_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 1)
call auto2e_W_6_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 22, 1)
call auto2e_W_7_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 50, 1)
call auto2e_W_5_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 22)
call auto2e_W_6_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 22, 22)
call auto2e_W_7_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 50, 22)
call auto2e_BraKetTransfer_5_2_5_1(G, W, Ra, Rb, Rc, Rd)
end subroutine auto2e_eri_5_2_5_1

subroutine auto2e_eri_Spher_5_2_5_1(G, Ra, CntrA, ExpA, NprimA, &
Rb, CntrB, ExpB, NprimB, Rc, CntrC, ExpC, NprimC, &
Rd, CntrD, ExpD, NprimD, Kappa)
!
! Two-electron integrals class: (hd|hp)
! L-dependent normalization constants are not applied.
! Output memory layout: CDAB
!
real(F64), dimension(55, *), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(85, 49) :: W
!
! Compute (LS|KS) integrals, L=5..7 and K=5..6
!
W = ZERO
call auto2e_W_5_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 1)
call auto2e_W_6_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 22, 1)
call auto2e_W_7_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 50, 1)
call auto2e_W_5_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 22)
call auto2e_W_6_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 22, 22)
call auto2e_W_7_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 50, 22)
call auto2e_BraKetTransfer_Spher_5_2_5_1(G, W, Ra, Rb, Rc, Rd)
end subroutine auto2e_eri_Spher_5_2_5_1

subroutine auto2e_BraKetTransfer_5_2_5_1(G, W, Ra, Rb, Rc, Rd)
!
! Transform (LS|KS), L=5..7, K=5..6 using the horizontal recurrence relation
! of Eq. 18 in Head-Gordon, M. and Pople, J.A. J. Chem. Phys. 89, 5777 (1988);
! doi: 10.1063/1.455553
!
! The output integrals are expressed in the Cartesian GTOs basis.
!
real(F64), dimension(63, *), intent(out) :: G
real(F64), dimension(85, 49), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(3) :: AB, CD
real(F64), dimension(85, 63) :: T
real(F64), dimension(63, 85) :: U
AB = Ra - Rb
CD = Rc - Rd
call auto2e_KetTransfer_5_1(T, W, CD(1), CD(2), CD(3))
U = transpose(T)
call auto2e_KetTransfer_5_2(G(:, 1:126), U, AB(1), AB(2), AB(3))
end subroutine auto2e_BraKetTransfer_5_2_5_1

subroutine auto2e_BraKetTransfer_Spher_5_2_5_1(G, W, Ra, Rb, Rc, Rd)
!
! Transform (LS|KS), L=5..7, K=5..6 using the horizontal recurrence relation
! of Eq. 18 in Head-Gordon, M. and Pople, J.A. J. Chem. Phys. 89, 5777 (1988);
! doi: 10.1063/1.455553
!
! The output integrals are expressed in the real spherical harmonics basis.
!
real(F64), dimension(33, *), intent(out) :: G
real(F64), dimension(85, 49), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(3) :: AB, CD
real(F64), dimension(85, 63) :: Txyz
real(F64), dimension(85, 33) :: T
real(F64), dimension(33, 85) :: U
real(F64), dimension(33, 126) :: Gxyz
AB = Ra - Rb
CD = Rc - Rd
call auto2e_KetTransfer_5_1(Txyz, W, CD(1), CD(2), CD(3))
call auto2e_SpherTransf_U_Vector_5_1(T, Txyz, 85)
U = transpose(T)
call auto2e_KetTransfer_5_2(Gxyz, U, AB(1), AB(2), AB(3))
call auto2e_SpherTransf_U_Vector_5_2(G(:, 1:55), Gxyz, 33)
end subroutine auto2e_BraKetTransfer_Spher_5_2_5_1

subroutine auto2e_Normalize_5_2_5_1_ABCD(Gabcd, Gcdab, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: CDAB
! Output memory layout: ABCD
!
real(F64), dimension(*), intent(out) :: Gabcd
real(F64), dimension(*), intent(in) :: Gcdab
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do b = 1, 6
do a = 1, 21
do d = 1, 3
do c = 1, 21
w = a + (b-1)*21 + (c-1)*126 + (d-1)*2646
Gabcd(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_5_2_5_1_ABCD

subroutine auto2e_Normalize_5_2_5_1_ABDC(Gabdc, Gcdab, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: CDAB
! Output memory layout: ABDC
!
real(F64), dimension(*), intent(out) :: Gabdc
real(F64), dimension(*), intent(in) :: Gcdab
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do b = 1, 6
do a = 1, 21
do d = 1, 3
do c = 1, 21
w = a + (b-1)*21 + (d-1)*126 + (c-1)*378
Gabdc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_5_2_5_1_ABDC

subroutine auto2e_Normalize_5_2_5_1_BACD(Gbacd, Gcdab, NormA, NormB, NormC, NormD)
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
do a = 1, 21
do d = 1, 3
do c = 1, 21
w = b + (a-1)*6 + (c-1)*126 + (d-1)*2646
Gbacd(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_5_2_5_1_BACD

subroutine auto2e_Normalize_5_2_5_1_BADC(Gbadc, Gcdab, NormA, NormB, NormC, NormD)
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
do a = 1, 21
do d = 1, 3
do c = 1, 21
w = b + (a-1)*6 + (d-1)*126 + (c-1)*378
Gbadc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_5_2_5_1_BADC

subroutine auto2e_Normalize_5_2_5_1_CDAB(G, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants
!
real(F64), dimension(*), intent(inout) :: G
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v
v = 1
do b = 1, 6
do a = 1, 21
do d = 1, 3
do c = 1, 21
G(v) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*G(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_5_2_5_1_CDAB

subroutine auto2e_Normalize_5_2_5_1_CDBA(Gcdba, Gcdab, NormA, NormB, NormC, NormD)
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
do a = 1, 21
do d = 1, 3
do c = 1, 21
w = c + (d-1)*21 + (b-1)*63 + (a-1)*378
Gcdba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_5_2_5_1_CDBA

subroutine auto2e_Normalize_5_2_5_1_DCAB(Gdcab, Gcdab, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: CDAB
! Output memory layout: DCAB
!
real(F64), dimension(*), intent(out) :: Gdcab
real(F64), dimension(*), intent(in) :: Gcdab
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do b = 1, 6
do a = 1, 21
do d = 1, 3
do c = 1, 21
w = d + (c-1)*3 + (a-1)*63 + (b-1)*1323
Gdcab(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_5_2_5_1_DCAB

subroutine auto2e_Normalize_5_2_5_1_DCBA(Gdcba, Gcdab, NormA, NormB, NormC, NormD)
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
do a = 1, 21
do d = 1, 3
do c = 1, 21
w = d + (c-1)*3 + (b-1)*63 + (a-1)*378
Gdcba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_5_2_5_1_DCBA

subroutine auto2e_Normalize_Spher_5_2_5_1_ABCD(Gabcd, Gcdab, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: CDAB
! Output memory layout: ABCD
!
real(F64), dimension(*), intent(out) :: Gabcd
real(F64), dimension(*), intent(in) :: Gcdab
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do b = 1, 5
do a = 1, 11
do d = 1, 3
do c = 1, 11
w = a + (b-1)*11 + (c-1)*55 + (d-1)*605
Gabcd(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_5_2_5_1_ABCD

subroutine auto2e_Normalize_Spher_5_2_5_1_ABDC(Gabdc, Gcdab, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: CDAB
! Output memory layout: ABDC
!
real(F64), dimension(*), intent(out) :: Gabdc
real(F64), dimension(*), intent(in) :: Gcdab
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do b = 1, 5
do a = 1, 11
do d = 1, 3
do c = 1, 11
w = a + (b-1)*11 + (d-1)*55 + (c-1)*165
Gabdc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_5_2_5_1_ABDC

subroutine auto2e_Normalize_Spher_5_2_5_1_BACD(Gbacd, Gcdab, NormA, NormB, NormC, NormD)
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
do a = 1, 11
do d = 1, 3
do c = 1, 11
w = b + (a-1)*5 + (c-1)*55 + (d-1)*605
Gbacd(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_5_2_5_1_BACD

subroutine auto2e_Normalize_Spher_5_2_5_1_BADC(Gbadc, Gcdab, NormA, NormB, NormC, NormD)
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
do a = 1, 11
do d = 1, 3
do c = 1, 11
w = b + (a-1)*5 + (d-1)*55 + (c-1)*165
Gbadc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_5_2_5_1_BADC

subroutine auto2e_Normalize_Spher_5_2_5_1_CDAB(G, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants
!
real(F64), dimension(*), intent(inout) :: G
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v
v = 1
do b = 1, 5
do a = 1, 11
do d = 1, 3
do c = 1, 11
G(v) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*G(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_5_2_5_1_CDAB

subroutine auto2e_Normalize_Spher_5_2_5_1_CDBA(Gcdba, Gcdab, NormA, NormB, NormC, NormD)
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
do a = 1, 11
do d = 1, 3
do c = 1, 11
w = c + (d-1)*11 + (b-1)*33 + (a-1)*165
Gcdba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_5_2_5_1_CDBA

subroutine auto2e_Normalize_Spher_5_2_5_1_DCAB(Gdcab, Gcdab, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: CDAB
! Output memory layout: DCAB
!
real(F64), dimension(*), intent(out) :: Gdcab
real(F64), dimension(*), intent(in) :: Gcdab
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do b = 1, 5
do a = 1, 11
do d = 1, 3
do c = 1, 11
w = d + (c-1)*3 + (a-1)*33 + (b-1)*363
Gdcab(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_5_2_5_1_DCAB

subroutine auto2e_Normalize_Spher_5_2_5_1_DCBA(Gdcba, Gcdab, NormA, NormB, NormC, NormD)
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
do a = 1, 11
do d = 1, 3
do c = 1, 11
w = d + (c-1)*3 + (b-1)*33 + (a-1)*165
Gdcba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_5_2_5_1_DCBA

subroutine auto2e_frontend_5_2_5_1(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(7938) :: H
call auto2e_eri_5_2_5_1(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_5_2_5_1_DCBA(G, H, NormA, NormB, NormC, NormD)
end subroutine auto2e_frontend_5_2_5_1

subroutine auto2e_frontend_Spher_5_2_5_1(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(1815) :: H
call auto2e_eri_Spher_5_2_5_1(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_Spher_5_2_5_1_DCBA(G, H, NormA, NormB, NormC, NormD)
end subroutine auto2e_frontend_Spher_5_2_5_1

subroutine auto2e_frontend_1_5_5_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(7938) :: H
call auto2e_eri_5_2_5_1(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, Kappa)
call auto2e_Normalize_5_2_5_1_BACD(G, H, NormC, NormD, NormB, NormA)
end subroutine auto2e_frontend_1_5_5_2

subroutine auto2e_frontend_Spher_1_5_5_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(1815) :: H
call auto2e_eri_Spher_5_2_5_1(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, Kappa)
call auto2e_Normalize_Spher_5_2_5_1_BACD(G, H, NormC, NormD, NormB, NormA)
end subroutine auto2e_frontend_Spher_1_5_5_2

subroutine auto2e_frontend_5_1_5_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(7938) :: H
call auto2e_eri_5_2_5_1(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_5_2_5_1_BADC(G, H, NormC, NormD, NormA, NormB)
end subroutine auto2e_frontend_5_1_5_2

subroutine auto2e_frontend_Spher_5_1_5_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(1815) :: H
call auto2e_eri_Spher_5_2_5_1(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_Spher_5_2_5_1_BADC(G, H, NormC, NormD, NormA, NormB)
end subroutine auto2e_frontend_Spher_5_1_5_2

subroutine auto2e_frontend_2_5_1_5(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
call auto2e_eri_5_2_5_1(G, Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, &
   Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, Kappa)
call auto2e_Normalize_5_2_5_1_CDAB(G, NormB, NormA, NormD, NormC)
end subroutine auto2e_frontend_2_5_1_5

subroutine auto2e_frontend_Spher_2_5_1_5(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
call auto2e_eri_Spher_5_2_5_1(G, Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, &
   Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, Kappa)
call auto2e_Normalize_Spher_5_2_5_1_CDAB(G, NormB, NormA, NormD, NormC)
end subroutine auto2e_frontend_Spher_2_5_1_5

subroutine auto2e_frontend_5_1_2_5(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(7938) :: H
call auto2e_eri_5_2_5_1(H, Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_5_2_5_1_ABDC(G, H, NormD, NormC, NormA, NormB)
end subroutine auto2e_frontend_5_1_2_5

subroutine auto2e_frontend_Spher_5_1_2_5(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(1815) :: H
call auto2e_eri_Spher_5_2_5_1(H, Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_Spher_5_2_5_1_ABDC(G, H, NormD, NormC, NormA, NormB)
end subroutine auto2e_frontend_Spher_5_1_2_5

subroutine auto2e_frontend_1_5_2_5(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(7938) :: H
call auto2e_eri_5_2_5_1(H, Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, &
   Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, Kappa)
call auto2e_Normalize_5_2_5_1_ABCD(G, H, NormD, NormC, NormB, NormA)
end subroutine auto2e_frontend_1_5_2_5

subroutine auto2e_frontend_Spher_1_5_2_5(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(1815) :: H
call auto2e_eri_Spher_5_2_5_1(H, Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, &
   Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, Kappa)
call auto2e_Normalize_Spher_5_2_5_1_ABCD(G, H, NormD, NormC, NormB, NormA)
end subroutine auto2e_frontend_Spher_1_5_2_5

subroutine auto2e_frontend_5_2_1_5(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(7938) :: H
call auto2e_eri_5_2_5_1(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, Kappa)
call auto2e_Normalize_5_2_5_1_CDBA(G, H, NormA, NormB, NormD, NormC)
end subroutine auto2e_frontend_5_2_1_5

subroutine auto2e_frontend_Spher_5_2_1_5(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(1815) :: H
call auto2e_eri_Spher_5_2_5_1(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, Kappa)
call auto2e_Normalize_Spher_5_2_5_1_CDBA(G, H, NormA, NormB, NormD, NormC)
end subroutine auto2e_frontend_Spher_5_2_1_5

subroutine auto2e_frontend_2_5_5_1(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(7938) :: H
call auto2e_eri_5_2_5_1(H, Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_5_2_5_1_DCAB(G, H, NormB, NormA, NormC, NormD)
end subroutine auto2e_frontend_2_5_5_1

subroutine auto2e_frontend_Spher_2_5_5_1(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(1815) :: H
call auto2e_eri_Spher_5_2_5_1(H, Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_Spher_5_2_5_1_DCAB(G, H, NormB, NormA, NormC, NormD)
end subroutine auto2e_frontend_Spher_2_5_5_1
end module auto2e_eri_hdhp