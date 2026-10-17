!
! Compute a batch of two-electron repulsion integrals, (AB|CD), in the Cartesian Gaussian
! and real spherical harmonic basis.
! Angular momenta: A=6, B=3, C=3, and D=3
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
module auto2e_eri_ifff
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

subroutine auto2e_eri_6_3_3_3(G, Ra, CntrA, ExpA, NprimA, &
Rb, CntrB, ExpB, NprimB, Rc, CntrC, ExpC, NprimC, &
Rd, CntrD, ExpD, NprimD, Kappa)
!
! Two-electron integrals class: (if|ff)
! L-dependent normalization constants are not applied.
! Output memory layout: CDAB
!
real(F64), dimension(280, *), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(164, 74) :: W
!
! Compute (LS|KS) integrals, L=6..9 and K=3..6
!
W = ZERO
call auto2e_W_6_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 1)
call auto2e_W_7_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 1)
call auto2e_W_8_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 1)
call auto2e_W_9_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 1)
call auto2e_W_6_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 11)
call auto2e_W_7_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 11)
call auto2e_W_8_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 11)
call auto2e_W_9_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 11)
call auto2e_W_6_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 26)
call auto2e_W_7_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 26)
call auto2e_W_8_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 26)
call auto2e_W_9_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 26)
call auto2e_W_6_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 47)
call auto2e_W_7_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 47)
call auto2e_W_8_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 47)
call auto2e_W_9_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 47)
call auto2e_BraKetTransfer_6_3_3_3(G, W, Ra, Rb, Rc, Rd)
end subroutine auto2e_eri_6_3_3_3

subroutine auto2e_eri_Spher_6_3_3_3(G, Ra, CntrA, ExpA, NprimA, &
Rb, CntrB, ExpB, NprimB, Rc, CntrC, ExpC, NprimC, &
Rd, CntrD, ExpD, NprimD, Kappa)
!
! Two-electron integrals class: (if|ff)
! L-dependent normalization constants are not applied.
! Output memory layout: CDAB
!
real(F64), dimension(91, *), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(164, 74) :: W
!
! Compute (LS|KS) integrals, L=6..9 and K=3..6
!
W = ZERO
call auto2e_W_6_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 1)
call auto2e_W_7_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 1)
call auto2e_W_8_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 1)
call auto2e_W_9_3(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 1)
call auto2e_W_6_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 11)
call auto2e_W_7_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 11)
call auto2e_W_8_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 11)
call auto2e_W_9_4(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 11)
call auto2e_W_6_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 26)
call auto2e_W_7_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 26)
call auto2e_W_8_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 26)
call auto2e_W_9_5(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 26)
call auto2e_W_6_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 47)
call auto2e_W_7_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 47)
call auto2e_W_8_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 47)
call auto2e_W_9_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 47)
call auto2e_BraKetTransfer_Spher_6_3_3_3(G, W, Ra, Rb, Rc, Rd)
end subroutine auto2e_eri_Spher_6_3_3_3

subroutine auto2e_BraKetTransfer_6_3_3_3(G, W, Ra, Rb, Rc, Rd)
!
! Transform (LS|KS), L=6..9, K=3..6 using the horizontal recurrence relation
! of Eq. 18 in Head-Gordon, M. and Pople, J.A. J. Chem. Phys. 89, 5777 (1988);
! doi: 10.1063/1.455553
!
! The output integrals are expressed in the Cartesian GTOs basis.
!
real(F64), dimension(100, *), intent(out) :: G
real(F64), dimension(164, 74), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(3) :: AB, CD
real(F64), dimension(164, 100) :: T
real(F64), dimension(100, 164) :: U
AB = Ra - Rb
CD = Rc - Rd
call auto2e_KetTransfer_3_3(T, W, CD(1), CD(2), CD(3))
U = transpose(T)
call auto2e_KetTransfer_6_3(G(:, 1:280), U, AB(1), AB(2), AB(3))
end subroutine auto2e_BraKetTransfer_6_3_3_3

subroutine auto2e_BraKetTransfer_Spher_6_3_3_3(G, W, Ra, Rb, Rc, Rd)
!
! Transform (LS|KS), L=6..9, K=3..6 using the horizontal recurrence relation
! of Eq. 18 in Head-Gordon, M. and Pople, J.A. J. Chem. Phys. 89, 5777 (1988);
! doi: 10.1063/1.455553
!
! The output integrals are expressed in the real spherical harmonics basis.
!
real(F64), dimension(49, *), intent(out) :: G
real(F64), dimension(164, 74), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(3) :: AB, CD
real(F64), dimension(164, 100) :: Txyz
real(F64), dimension(164, 49) :: T
real(F64), dimension(49, 164) :: U
real(F64), dimension(49, 280) :: Gxyz
AB = Ra - Rb
CD = Rc - Rd
call auto2e_KetTransfer_3_3(Txyz, W, CD(1), CD(2), CD(3))
call auto2e_SpherTransf_U_Vector_3_3(T, Txyz, 164)
U = transpose(T)
call auto2e_KetTransfer_6_3(Gxyz, U, AB(1), AB(2), AB(3))
call auto2e_SpherTransf_U_Vector_6_3(G(:, 1:91), Gxyz, 49)
end subroutine auto2e_BraKetTransfer_Spher_6_3_3_3

subroutine auto2e_Normalize_6_3_3_3_ABDC(Gabdc, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 10
do a = 1, 28
do d = 1, 10
do c = 1, 10
w = a + (b-1)*28 + (d-1)*280 + (c-1)*2800
Gabdc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_6_3_3_3_ABDC

subroutine auto2e_Normalize_6_3_3_3_BADC(Gbadc, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 10
do a = 1, 28
do d = 1, 10
do c = 1, 10
w = b + (a-1)*10 + (d-1)*280 + (c-1)*2800
Gbadc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_6_3_3_3_BADC

subroutine auto2e_Normalize_6_3_3_3_DCAB(Gdcab, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 10
do a = 1, 28
do d = 1, 10
do c = 1, 10
w = d + (c-1)*10 + (a-1)*100 + (b-1)*2800
Gdcab(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_6_3_3_3_DCAB

subroutine auto2e_Normalize_6_3_3_3_DCBA(Gdcba, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 10
do a = 1, 28
do d = 1, 10
do c = 1, 10
w = d + (c-1)*10 + (b-1)*100 + (a-1)*1000
Gdcba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_6_3_3_3_DCBA

subroutine auto2e_Normalize_Spher_6_3_3_3_ABDC(Gabdc, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 7
do a = 1, 13
do d = 1, 7
do c = 1, 7
w = a + (b-1)*13 + (d-1)*91 + (c-1)*637
Gabdc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_6_3_3_3_ABDC

subroutine auto2e_Normalize_Spher_6_3_3_3_BADC(Gbadc, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 7
do a = 1, 13
do d = 1, 7
do c = 1, 7
w = b + (a-1)*7 + (d-1)*91 + (c-1)*637
Gbadc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_6_3_3_3_BADC

subroutine auto2e_Normalize_Spher_6_3_3_3_DCAB(Gdcab, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 7
do a = 1, 13
do d = 1, 7
do c = 1, 7
w = d + (c-1)*7 + (a-1)*49 + (b-1)*637
Gdcab(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_6_3_3_3_DCAB

subroutine auto2e_Normalize_Spher_6_3_3_3_DCBA(Gdcba, Gcdab, NormA, NormB, NormC, NormD)
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
do b = 1, 7
do a = 1, 13
do d = 1, 7
do c = 1, 7
w = d + (c-1)*7 + (b-1)*49 + (a-1)*343
Gdcba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_6_3_3_3_DCBA

subroutine auto2e_frontend_3_3_6_3(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(28000) :: H
call auto2e_eri_6_3_3_3(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_6_3_3_3_BADC(G, H, NormC, NormD, NormA, NormB)
end subroutine auto2e_frontend_3_3_6_3

subroutine auto2e_frontend_Spher_3_3_6_3(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(4459) :: H
call auto2e_eri_Spher_6_3_3_3(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_Spher_6_3_3_3_BADC(G, H, NormC, NormD, NormA, NormB)
end subroutine auto2e_frontend_Spher_3_3_6_3

subroutine auto2e_frontend_6_3_3_3(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(28000) :: H
call auto2e_eri_6_3_3_3(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_6_3_3_3_DCBA(G, H, NormA, NormB, NormC, NormD)
end subroutine auto2e_frontend_6_3_3_3

subroutine auto2e_frontend_Spher_6_3_3_3(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(4459) :: H
call auto2e_eri_Spher_6_3_3_3(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_Spher_6_3_3_3_DCBA(G, H, NormA, NormB, NormC, NormD)
end subroutine auto2e_frontend_Spher_6_3_3_3

subroutine auto2e_frontend_3_6_3_3(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(28000) :: H
call auto2e_eri_6_3_3_3(H, Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_6_3_3_3_DCAB(G, H, NormB, NormA, NormC, NormD)
end subroutine auto2e_frontend_3_6_3_3

subroutine auto2e_frontend_Spher_3_6_3_3(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(4459) :: H
call auto2e_eri_Spher_6_3_3_3(H, Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_Spher_6_3_3_3_DCAB(G, H, NormB, NormA, NormC, NormD)
end subroutine auto2e_frontend_Spher_3_6_3_3

subroutine auto2e_frontend_3_3_3_6(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(28000) :: H
call auto2e_eri_6_3_3_3(H, Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_6_3_3_3_ABDC(G, H, NormD, NormC, NormA, NormB)
end subroutine auto2e_frontend_3_3_3_6

subroutine auto2e_frontend_Spher_3_3_3_6(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(4459) :: H
call auto2e_eri_Spher_6_3_3_3(H, Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_Spher_6_3_3_3_ABDC(G, H, NormD, NormC, NormA, NormB)
end subroutine auto2e_frontend_Spher_3_3_3_6
end module auto2e_eri_ifff