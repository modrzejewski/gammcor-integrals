!
! Compute a batch of two-electron repulsion integrals, (AB|CD), in the Cartesian Gaussian
! and real spherical harmonic basis.
! Angular momenta: A=5, B=0, C=4, and D=2
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
module auto2e_eri_hsgd
use arithmetic
use auto2e_KetTransfer
use auto2e_SpherTransf
use auto2e_WMatrix_part1
use auto2e_WMatrix_part2
use auto2e_WMatrix_part3
implicit none

contains

subroutine auto2e_eri_5_0_4_2(G, Ra, CntrA, ExpA, NprimA, &
Rb, CntrB, ExpB, NprimB, Rc, CntrC, ExpC, NprimC, &
Rd, CntrD, ExpD, NprimD, Kappa)
!
! Two-electron integrals class: (hs|gd)
! L-dependent normalization constants are not applied.
! Output memory layout: ABCD
!
real(F64), dimension(21, *), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(21, 64) :: W
!
! Compute (LS|KS) integrals, L=5..5 and K=4..6
!
W = ZERO
call auto2e_W_5_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 1)
call auto2e_W_5_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 16)
call auto2e_W_5_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 37)
call auto2e_BraKetTransfer_5_0_4_2(G, W, Rc, Rd)
end subroutine auto2e_eri_5_0_4_2

subroutine auto2e_eri_Spher_5_0_4_2(G, Ra, CntrA, ExpA, NprimA, &
Rb, CntrB, ExpB, NprimB, Rc, CntrC, ExpC, NprimC, &
Rd, CntrD, ExpD, NprimD, Kappa)
!
! Two-electron integrals class: (hs|gd)
! L-dependent normalization constants are not applied.
! Output memory layout: CDAB
!
real(F64), dimension(11, *), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(21, 64) :: W
!
! Compute (LS|KS) integrals, L=5..5 and K=4..6
!
W = ZERO
call auto2e_W_5_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 1)
call auto2e_W_5_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 16)
call auto2e_W_5_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 37)
call auto2e_BraKetTransfer_Spher_5_0_4_2(G, W, Rc, Rd)
end subroutine auto2e_eri_Spher_5_0_4_2

subroutine auto2e_BraKetTransfer_5_0_4_2(G, W, Rc, Rd)
!
! Transform (LS|KS), L=5..5, K=4..6 using the horizontal recurrence relation
! of Eq. 18 in Head-Gordon, M. and Pople, J.A. J. Chem. Phys. 89, 5777 (1988);
! doi: 10.1063/1.455553
!
! The output integrals are expressed in the Cartesian GTOs basis.
!
real(F64), dimension(21, *), intent(out) :: G
real(F64), dimension(21, 64), intent(in) :: W
real(F64), dimension(3), intent(in) :: Rc, Rd
real(F64), dimension(3) :: CD
CD = Rc - Rd
call auto2e_KetTransfer_4_2(G(:, 1:90), W, CD(1), CD(2), CD(3))
end subroutine auto2e_BraKetTransfer_5_0_4_2

subroutine auto2e_BraKetTransfer_Spher_5_0_4_2(G, W, Rc, Rd)
!
! Transform (LS|KS), L=5..5, K=4..6 using the horizontal recurrence relation
! of Eq. 18 in Head-Gordon, M. and Pople, J.A. J. Chem. Phys. 89, 5777 (1988);
! doi: 10.1063/1.455553
!
! The output integrals are expressed in the real spherical harmonics basis.
!
real(F64), dimension(45, *), intent(out) :: G
real(F64), dimension(21, 64), intent(in) :: W
real(F64), dimension(3), intent(in) :: Rc, Rd
real(F64), dimension(3) :: CD
real(F64), dimension(21, 45) :: T
real(F64), dimension(21, 90) :: Txyz
real(F64), dimension(45, 21) :: Gxyz
CD = Rc - Rd
call auto2e_KetTransfer_4_2(Txyz, W, CD(1), CD(2), CD(3))
call auto2e_SpherTransf_U_Vector_4_2(T, Txyz, 21)
Gxyz = transpose(T)
call auto2e_SpherTransf_U_Vector_5_0(G(:, 1:11), Gxyz, 45)
end subroutine auto2e_BraKetTransfer_Spher_5_0_4_2

subroutine auto2e_Normalize_5_0_4_2_ABCD(G, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants
!
real(F64), dimension(*), intent(inout) :: G
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v
v = 1
do d = 1, 6
do c = 1, 15
do b = 1, 1
do a = 1, 21
G(v) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*G(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_5_0_4_2_ABCD

subroutine auto2e_Normalize_5_0_4_2_ABDC(Gabdc, Gabcd, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: ABCD
! Output memory layout: ABDC
!
real(F64), dimension(*), intent(out) :: Gabdc
real(F64), dimension(*), intent(in) :: Gabcd
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do d = 1, 6
do c = 1, 15
do b = 1, 1
do a = 1, 21
w = a + (b-1)*21 + (d-1)*21 + (c-1)*126
Gabdc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gabcd(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_5_0_4_2_ABDC

subroutine auto2e_Normalize_5_0_4_2_CDAB(Gcdab, Gabcd, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: ABCD
! Output memory layout: CDAB
!
real(F64), dimension(*), intent(out) :: Gcdab
real(F64), dimension(*), intent(in) :: Gabcd
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do d = 1, 6
do c = 1, 15
do b = 1, 1
do a = 1, 21
w = c + (d-1)*15 + (a-1)*90 + (b-1)*1890
Gcdab(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gabcd(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_5_0_4_2_CDAB

subroutine auto2e_Normalize_5_0_4_2_DCAB(Gdcab, Gabcd, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: ABCD
! Output memory layout: DCAB
!
real(F64), dimension(*), intent(out) :: Gdcab
real(F64), dimension(*), intent(in) :: Gabcd
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do d = 1, 6
do c = 1, 15
do b = 1, 1
do a = 1, 21
w = d + (c-1)*6 + (a-1)*90 + (b-1)*1890
Gdcab(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gabcd(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_5_0_4_2_DCAB

subroutine auto2e_Normalize_Spher_5_0_4_2_ABCD(Gabcd, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 1
do a = 1, 11
do d = 1, 5
do c = 1, 9
w = a + (b-1)*11 + (c-1)*11 + (d-1)*99
Gabcd(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_5_0_4_2_ABCD

subroutine auto2e_Normalize_Spher_5_0_4_2_ABDC(Gabdc, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 1
do a = 1, 11
do d = 1, 5
do c = 1, 9
w = a + (b-1)*11 + (d-1)*11 + (c-1)*55
Gabdc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_5_0_4_2_ABDC

subroutine auto2e_Normalize_Spher_5_0_4_2_CDAB(G, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants
!
real(F64), dimension(*), intent(inout) :: G
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v
v = 1
do b = 1, 1
do a = 1, 11
do d = 1, 5
do c = 1, 9
G(v) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*G(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_5_0_4_2_CDAB

subroutine auto2e_Normalize_Spher_5_0_4_2_DCAB(Gdcab, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 1
do a = 1, 11
do d = 1, 5
do c = 1, 9
w = d + (c-1)*5 + (a-1)*45 + (b-1)*495
Gdcab(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_5_0_4_2_DCAB

subroutine auto2e_frontend_2_4_5_0(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
call auto2e_eri_5_0_4_2(G, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, Kappa)
call auto2e_Normalize_5_0_4_2_ABCD(G, NormC, NormD, NormB, NormA)
end subroutine auto2e_frontend_2_4_5_0

subroutine auto2e_frontend_Spher_2_4_5_0(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(495) :: H
call auto2e_eri_Spher_5_0_4_2(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, Kappa)
call auto2e_Normalize_Spher_5_0_4_2_ABCD(G, H, NormC, NormD, NormB, NormA)
end subroutine auto2e_frontend_Spher_2_4_5_0

subroutine auto2e_frontend_2_4_0_5(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
call auto2e_eri_5_0_4_2(G, Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, &
   Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, Kappa)
call auto2e_Normalize_5_0_4_2_ABCD(G, NormD, NormC, NormB, NormA)
end subroutine auto2e_frontend_2_4_0_5

subroutine auto2e_frontend_Spher_2_4_0_5(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(495) :: H
call auto2e_eri_Spher_5_0_4_2(H, Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, &
   Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, Kappa)
call auto2e_Normalize_Spher_5_0_4_2_ABCD(G, H, NormD, NormC, NormB, NormA)
end subroutine auto2e_frontend_Spher_2_4_0_5

subroutine auto2e_frontend_0_5_4_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(1890) :: H
call auto2e_eri_5_0_4_2(H, Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_5_0_4_2_DCAB(G, H, NormB, NormA, NormC, NormD)
end subroutine auto2e_frontend_0_5_4_2

subroutine auto2e_frontend_Spher_0_5_4_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(495) :: H
call auto2e_eri_Spher_5_0_4_2(H, Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_Spher_5_0_4_2_DCAB(G, H, NormB, NormA, NormC, NormD)
end subroutine auto2e_frontend_Spher_0_5_4_2

subroutine auto2e_frontend_4_2_0_5(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(1890) :: H
call auto2e_eri_5_0_4_2(H, Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_5_0_4_2_ABDC(G, H, NormD, NormC, NormA, NormB)
end subroutine auto2e_frontend_4_2_0_5

subroutine auto2e_frontend_Spher_4_2_0_5(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(495) :: H
call auto2e_eri_Spher_5_0_4_2(H, Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_Spher_5_0_4_2_ABDC(G, H, NormD, NormC, NormA, NormB)
end subroutine auto2e_frontend_Spher_4_2_0_5

subroutine auto2e_frontend_5_0_4_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(1890) :: H
call auto2e_eri_5_0_4_2(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_5_0_4_2_DCAB(G, H, NormA, NormB, NormC, NormD)
end subroutine auto2e_frontend_5_0_4_2

subroutine auto2e_frontend_Spher_5_0_4_2(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(495) :: H
call auto2e_eri_Spher_5_0_4_2(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_Spher_5_0_4_2_DCAB(G, H, NormA, NormB, NormC, NormD)
end subroutine auto2e_frontend_Spher_5_0_4_2

subroutine auto2e_frontend_0_5_2_4(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(1890) :: H
call auto2e_eri_5_0_4_2(H, Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, &
   Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, Kappa)
call auto2e_Normalize_5_0_4_2_CDAB(G, H, NormB, NormA, NormD, NormC)
end subroutine auto2e_frontend_0_5_2_4

subroutine auto2e_frontend_Spher_0_5_2_4(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
call auto2e_eri_Spher_5_0_4_2(G, Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, &
   Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, Kappa)
call auto2e_Normalize_Spher_5_0_4_2_CDAB(G, NormB, NormA, NormD, NormC)
end subroutine auto2e_frontend_Spher_0_5_2_4

subroutine auto2e_frontend_5_0_2_4(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(1890) :: H
call auto2e_eri_5_0_4_2(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, Kappa)
call auto2e_Normalize_5_0_4_2_CDAB(G, H, NormA, NormB, NormD, NormC)
end subroutine auto2e_frontend_5_0_2_4

subroutine auto2e_frontend_Spher_5_0_2_4(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
call auto2e_eri_Spher_5_0_4_2(G, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, Kappa)
call auto2e_Normalize_Spher_5_0_4_2_CDAB(G, NormA, NormB, NormD, NormC)
end subroutine auto2e_frontend_Spher_5_0_2_4

subroutine auto2e_frontend_4_2_5_0(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(1890) :: H
call auto2e_eri_5_0_4_2(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_5_0_4_2_ABDC(G, H, NormC, NormD, NormA, NormB)
end subroutine auto2e_frontend_4_2_5_0

subroutine auto2e_frontend_Spher_4_2_5_0(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(495) :: H
call auto2e_eri_Spher_5_0_4_2(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_Spher_5_0_4_2_ABDC(G, H, NormC, NormD, NormA, NormB)
end subroutine auto2e_frontend_Spher_4_2_5_0
end module auto2e_eri_hsgd