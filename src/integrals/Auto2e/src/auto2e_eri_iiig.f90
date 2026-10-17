!
! Compute a batch of two-electron repulsion integrals, (AB|CD), in the Cartesian Gaussian
! and real spherical harmonic basis.
! Angular momenta: A=6, B=6, C=6, and D=4
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
module auto2e_eri_iiig
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

subroutine auto2e_eri_6_6_6_4(G, Ra, CntrA, ExpA, NprimA, &
Rb, CntrB, ExpB, NprimB, Rc, CntrC, ExpC, NprimC, &
Rd, CntrD, ExpD, NprimD, Kappa)
!
! Two-electron integrals class: (ii|ig)
! L-dependent normalization constants are not applied.
! Output memory layout: CDAB
!
real(F64), dimension(784, *), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(399, 230) :: W
!
! Compute (LS|KS) integrals, L=6..12 and K=6..10
!
W = ZERO
call auto2e_W_6_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 1)
call auto2e_W_7_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 1)
call auto2e_W_8_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 1)
call auto2e_W_9_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 1)
call auto2e_W_10_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 165, 1)
call auto2e_W_11_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 231, 1)
call auto2e_W_12_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 309, 1)
call auto2e_W_6_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 29)
call auto2e_W_7_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 29)
call auto2e_W_8_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 29)
call auto2e_W_9_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 29)
call auto2e_W_10_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 165, 29)
call auto2e_W_11_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 231, 29)
call auto2e_W_12_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 309, 29)
call auto2e_W_6_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 65)
call auto2e_W_7_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 65)
call auto2e_W_8_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 65)
call auto2e_W_9_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 65)
call auto2e_W_10_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 165, 65)
call auto2e_W_11_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 231, 65)
call auto2e_W_12_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 309, 65)
call auto2e_W_6_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 110)
call auto2e_W_7_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 110)
call auto2e_W_8_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 110)
call auto2e_W_9_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 110)
call auto2e_W_10_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 165, 110)
call auto2e_W_11_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 231, 110)
call auto2e_W_12_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 309, 110)
call auto2e_W_6_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 165)
call auto2e_W_7_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 165)
call auto2e_W_8_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 165)
call auto2e_W_9_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 165)
call auto2e_W_10_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 165, 165)
call auto2e_W_11_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 231, 165)
call auto2e_W_12_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 309, 165)
call auto2e_BraKetTransfer_6_6_6_4(G, W, Ra, Rb, Rc, Rd)
end subroutine auto2e_eri_6_6_6_4

subroutine auto2e_eri_Spher_6_6_6_4(G, Ra, CntrA, ExpA, NprimA, &
Rb, CntrB, ExpB, NprimB, Rc, CntrC, ExpC, NprimC, &
Rd, CntrD, ExpD, NprimD, Kappa)
!
! Two-electron integrals class: (ii|ig)
! L-dependent normalization constants are not applied.
! Output memory layout: CDAB
!
real(F64), dimension(169, *), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(399, 230) :: W
!
! Compute (LS|KS) integrals, L=6..12 and K=6..10
!
W = ZERO
call auto2e_W_6_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 1)
call auto2e_W_7_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 1)
call auto2e_W_8_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 1)
call auto2e_W_9_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 1)
call auto2e_W_10_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 165, 1)
call auto2e_W_11_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 231, 1)
call auto2e_W_12_6(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 309, 1)
call auto2e_W_6_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 29)
call auto2e_W_7_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 29)
call auto2e_W_8_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 29)
call auto2e_W_9_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 29)
call auto2e_W_10_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 165, 29)
call auto2e_W_11_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 231, 29)
call auto2e_W_12_7(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 309, 29)
call auto2e_W_6_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 65)
call auto2e_W_7_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 65)
call auto2e_W_8_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 65)
call auto2e_W_9_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 65)
call auto2e_W_10_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 165, 65)
call auto2e_W_11_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 231, 65)
call auto2e_W_12_8(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 309, 65)
call auto2e_W_6_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 110)
call auto2e_W_7_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 110)
call auto2e_W_8_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 110)
call auto2e_W_9_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 110)
call auto2e_W_10_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 165, 110)
call auto2e_W_11_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 231, 110)
call auto2e_W_12_9(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 309, 110)
call auto2e_W_6_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 1, 165)
call auto2e_W_7_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 29, 165)
call auto2e_W_8_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 65, 165)
call auto2e_W_9_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 110, 165)
call auto2e_W_10_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 165, 165)
call auto2e_W_11_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 231, 165)
call auto2e_W_12_10(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, 309, 165)
call auto2e_BraKetTransfer_Spher_6_6_6_4(G, W, Ra, Rb, Rc, Rd)
end subroutine auto2e_eri_Spher_6_6_6_4

subroutine auto2e_BraKetTransfer_6_6_6_4(G, W, Ra, Rb, Rc, Rd)
!
! Transform (LS|KS), L=6..12, K=6..10 using the horizontal recurrence relation
! of Eq. 18 in Head-Gordon, M. and Pople, J.A. J. Chem. Phys. 89, 5777 (1988);
! doi: 10.1063/1.455553
!
! The output integrals are expressed in the Cartesian GTOs basis.
!
real(F64), dimension(420, *), intent(out) :: G
real(F64), dimension(399, 230), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(3) :: AB, CD
real(F64), dimension(399, 420) :: T
real(F64), dimension(420, 399) :: U
AB = Ra - Rb
CD = Rc - Rd
call auto2e_KetTransfer_6_4(T, W, CD(1), CD(2), CD(3))
U = transpose(T)
call auto2e_KetTransfer_6_6(G(:, 1:784), U, AB(1), AB(2), AB(3))
end subroutine auto2e_BraKetTransfer_6_6_6_4

subroutine auto2e_BraKetTransfer_Spher_6_6_6_4(G, W, Ra, Rb, Rc, Rd)
!
! Transform (LS|KS), L=6..12, K=6..10 using the horizontal recurrence relation
! of Eq. 18 in Head-Gordon, M. and Pople, J.A. J. Chem. Phys. 89, 5777 (1988);
! doi: 10.1063/1.455553
!
! The output integrals are expressed in the real spherical harmonics basis.
!
real(F64), dimension(117, *), intent(out) :: G
real(F64), dimension(399, 230), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(3) :: AB, CD
real(F64), dimension(399, 420) :: Txyz
real(F64), dimension(399, 117) :: T
real(F64), dimension(117, 399) :: U
real(F64), dimension(117, 784) :: Gxyz
AB = Ra - Rb
CD = Rc - Rd
call auto2e_KetTransfer_6_4(Txyz, W, CD(1), CD(2), CD(3))
call auto2e_SpherTransf_U_Vector_6_4(T, Txyz, 399)
U = transpose(T)
call auto2e_KetTransfer_6_6(Gxyz, U, AB(1), AB(2), AB(3))
call auto2e_SpherTransf_U_Vector_6_6(G(:, 1:169), Gxyz, 117)
end subroutine auto2e_BraKetTransfer_Spher_6_6_6_4

subroutine auto2e_Normalize_6_6_6_4_BACD(Gbacd, Gcdab, NormA, NormB, NormC, NormD)
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
do d = 1, 15
do c = 1, 28
w = b + (a-1)*28 + (c-1)*784 + (d-1)*21952
Gbacd(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_6_6_6_4_BACD

subroutine auto2e_Normalize_6_6_6_4_BADC(Gbadc, Gcdab, NormA, NormB, NormC, NormD)
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
do d = 1, 15
do c = 1, 28
w = b + (a-1)*28 + (d-1)*784 + (c-1)*11760
Gbadc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_6_6_6_4_BADC

subroutine auto2e_Normalize_6_6_6_4_CDBA(Gcdba, Gcdab, NormA, NormB, NormC, NormD)
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
do d = 1, 15
do c = 1, 28
w = c + (d-1)*28 + (b-1)*420 + (a-1)*11760
Gcdba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_6_6_6_4_CDBA

subroutine auto2e_Normalize_6_6_6_4_DCBA(Gdcba, Gcdab, NormA, NormB, NormC, NormD)
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
do d = 1, 15
do c = 1, 28
w = d + (c-1)*15 + (b-1)*420 + (a-1)*11760
Gdcba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_6_6_6_4_DCBA

subroutine auto2e_Normalize_Spher_6_6_6_4_BACD(Gbacd, Gcdab, NormA, NormB, NormC, NormD)
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
do d = 1, 9
do c = 1, 13
w = b + (a-1)*13 + (c-1)*169 + (d-1)*2197
Gbacd(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_6_6_6_4_BACD

subroutine auto2e_Normalize_Spher_6_6_6_4_BADC(Gbadc, Gcdab, NormA, NormB, NormC, NormD)
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
do d = 1, 9
do c = 1, 13
w = b + (a-1)*13 + (d-1)*169 + (c-1)*1521
Gbadc(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_6_6_6_4_BADC

subroutine auto2e_Normalize_Spher_6_6_6_4_CDBA(Gcdba, Gcdab, NormA, NormB, NormC, NormD)
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
do d = 1, 9
do c = 1, 13
w = c + (d-1)*13 + (b-1)*117 + (a-1)*1521
Gcdba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_6_6_6_4_CDBA

subroutine auto2e_Normalize_Spher_6_6_6_4_DCBA(Gdcba, Gcdab, NormA, NormB, NormC, NormD)
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
do d = 1, 9
do c = 1, 13
w = d + (c-1)*9 + (b-1)*117 + (a-1)*1521
Gdcba(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*Gcdab(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize_Spher_6_6_6_4_DCBA

subroutine auto2e_frontend_4_6_6_6(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(329280) :: H
call auto2e_eri_6_6_6_4(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, Kappa)
call auto2e_Normalize_6_6_6_4_BACD(G, H, NormC, NormD, NormB, NormA)
end subroutine auto2e_frontend_4_6_6_6

subroutine auto2e_frontend_Spher_4_6_6_6(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(19773) :: H
call auto2e_eri_Spher_6_6_6_4(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Rb, CntrB, ExpB, NprimB, &
   Ra, CntrA, ExpA, NprimA, Kappa)
call auto2e_Normalize_Spher_6_6_6_4_BACD(G, H, NormC, NormD, NormB, NormA)
end subroutine auto2e_frontend_Spher_4_6_6_6

subroutine auto2e_frontend_6_4_6_6(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(329280) :: H
call auto2e_eri_6_6_6_4(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_6_6_6_4_BADC(G, H, NormC, NormD, NormA, NormB)
end subroutine auto2e_frontend_6_4_6_6

subroutine auto2e_frontend_Spher_6_4_6_6(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(19773) :: H
call auto2e_eri_Spher_6_6_6_4(H, Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, &
   Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, Kappa)
call auto2e_Normalize_Spher_6_6_6_4_BADC(G, H, NormC, NormD, NormA, NormB)
end subroutine auto2e_frontend_Spher_6_4_6_6

subroutine auto2e_frontend_6_6_4_6(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(329280) :: H
call auto2e_eri_6_6_6_4(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, Kappa)
call auto2e_Normalize_6_6_6_4_CDBA(G, H, NormA, NormB, NormD, NormC)
end subroutine auto2e_frontend_6_6_4_6

subroutine auto2e_frontend_Spher_6_6_4_6(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(19773) :: H
call auto2e_eri_Spher_6_6_6_4(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rd, CntrD, ExpD, NprimD, &
   Rc, CntrC, ExpC, NprimC, Kappa)
call auto2e_Normalize_Spher_6_6_6_4_CDBA(G, H, NormA, NormB, NormD, NormC)
end subroutine auto2e_frontend_Spher_6_6_4_6

subroutine auto2e_frontend_6_6_6_4(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(329280) :: H
call auto2e_eri_6_6_6_4(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_6_6_6_4_DCBA(G, H, NormA, NormB, NormC, NormD)
end subroutine auto2e_frontend_6_6_6_4

subroutine auto2e_frontend_Spher_6_6_6_4(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
real(F64), dimension(19773) :: H
call auto2e_eri_Spher_6_6_6_4(H, Ra, CntrA, ExpA, NprimA, &
   Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, &
   Rd, CntrD, ExpD, NprimD, Kappa)
call auto2e_Normalize_Spher_6_6_6_4_DCBA(G, H, NormA, NormB, NormC, NormD)
end subroutine auto2e_frontend_Spher_6_6_6_4
end module auto2e_eri_iiig