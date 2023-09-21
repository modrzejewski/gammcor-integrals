module SpherGTO
      use arithmetic
      use spherh
      use Auto2e

      implicit none

contains

      subroutine SpherGTO_TransformMatrix(ASpher, ACart, LmaxGTO, NormFactorsSpher, &
            NormFactorsCart, ShellLocSpher, ShellLocCart, ShellMomentum, ShellParamsIdx, &
            NAOSpher, NAOCart, NShells, TransfWork)
            !
            ! Transform a symmetric matrix A from the Cartesian to
            ! the spherical harmonic AO basis. This transformation should
            ! be applied to density matrices built from AO Cartesian GTO
            ! coefficients.
            !
            real(F64), dimension(NAOSpher, NAOSpher), intent(out) :: ASpher
            real(F64), dimension(NAOCart, NAOCart), intent(in)    :: ACart
            integer, intent(in)                                   :: LmaxGTO
            real(F64), dimension(:, :), intent(in)                :: NormFactorsSpher
            real(F64), dimension(:, :), intent(in)                :: NormFactorsCart
            integer, dimension(NShells), intent(in)               :: ShellLocSpher
            integer, dimension(NShells), intent(in)               :: ShellLocCart
            integer, dimension(:), intent(in)                     :: ShellMomentum
            integer, dimension(:), intent(in)                     :: ShellParamsIdx
            integer, intent(in)                                   :: NAOSpher
            integer, intent(in)                                   :: NAOCart
            integer, intent(in)                                   :: NShells
            real(F64), dimension(NAOCart, NAOSpher), intent(out)  :: TransfWork

            integer :: k

            do k = 1, NAOCart
                  call SpherGTO_TransformVectors(ASpher(:, 1), ACart(:, k), LmaxGTO, NormFactorsSpher, &
                        NormFactorsCart, ShellLocSpher, ShellLocCart, ShellMomentum, ShellParamsIdx, &
                        NAOSpher, NAOCart, NShells, 1)
                  TransfWork(k, :) = ASpher(:, 1)
            end do
            call SpherGTO_TransformVectors(ASpher, TransfWork, LmaxGTO, NormFactorsSpher, &
                  NormFactorsCart, ShellLocSpher, ShellLocCart, ShellMomentum, ShellParamsIdx, &
                  NAOSpher, NAOCart, NShells, NAOSpher)
      end subroutine SpherGTO_TransformMatrix


      subroutine SpherGTO_TransformMatrix_U(ASpher, ACart, LmaxGTO, NormFactorsSpher, &
            NormFactorsCart, ShellLocSpher, ShellLocCart, ShellMomentum, ShellParamsIdx, &
            NAOSpher, NAOCart, NShells, TransfWork)
            !
            ! Transform a symmetric matrix A from the Cartesian to
            ! the spherical harmonic AO basis. This transformation should be
            ! applied to components of the hamiltonian expressed via AO
            ! Cartesian GTO integrals.
            !
            real(F64), dimension(NAOSpher, NAOSpher), intent(out) :: ASpher
            real(F64), dimension(NAOCart, NAOCart), intent(in)    :: ACart
            integer, intent(in)                                   :: LmaxGTO
            real(F64), dimension(:, :), intent(in)                :: NormFactorsSpher
            real(F64), dimension(:, :), intent(in)                :: NormFactorsCart
            integer, dimension(NShells), intent(in)               :: ShellLocSpher
            integer, dimension(NShells), intent(in)               :: ShellLocCart
            integer, dimension(:), intent(in)                     :: ShellMomentum
            integer, dimension(:), intent(in)                     :: ShellParamsIdx
            integer, intent(in)                                   :: NAOSpher
            integer, intent(in)                                   :: NAOCart
            integer, intent(in)                                   :: NShells
            real(F64), dimension(NAOCart, NAOSpher), intent(out)  :: TransfWork

            integer :: k

            do k = 1, NAOCart
                  call SpherGTO_TransformVectors_U(ASpher(:, 1), ACart(:, k), LmaxGTO, NormFactorsSpher, &
                        NormFactorsCart, ShellLocSpher, ShellLocCart, ShellMomentum, ShellParamsIdx, &
                        NAOSpher, NAOCart, NShells, 1)
                  TransfWork(k, :) = ASpher(:, 1)
            end do
            call SpherGTO_TransformVectors_U(ASpher, TransfWork, LmaxGTO, NormFactorsSpher, &
                  NormFactorsCart, ShellLocSpher, ShellLocCart, ShellMomentum, ShellParamsIdx, &
                  NAOSpher, NAOCart, NShells, NAOSpher)
      end subroutine SpherGTO_TransformMatrix_U
      
      
      subroutine SpherGTO_TransformVectors(ASpher, ACart, LmaxGTO, NormFactorsSpher, &
            NormFactorsCart, ShellLocSpher, ShellLocCart, ShellMomentum, ShellParamsIdx, &
            NAOSpher, NAOCart, NShells, NVecs)
            real(F64), dimension(NAOSpher, NVecs), intent(out) :: ASpher
            real(F64), dimension(NAOCart, NVecs), intent(in)   :: ACart
            integer, intent(in)                                :: LmaxGTO
            real(F64), dimension(:, :), intent(in)             :: NormFactorsSpher
            real(F64), dimension(:, :), intent(in)             :: NormFactorsCart
            integer, dimension(NShells), intent(in)            :: ShellLocSpher
            integer, dimension(NShells), intent(in)            :: ShellLocCart
            integer, dimension(:), intent(in)                  :: ShellMomentum
            integer, dimension(:), intent(in)                  :: ShellParamsIdx
            integer, intent(in)                                :: NAOSpher
            integer, intent(in)                                :: NAOCart
            integer, intent(in)                                :: NShells
            integer, intent(in)                                :: NVecs
            
            integer :: NSpher, NCart
            integer :: c0, c1, s0, s1, k
            integer :: ShA, ShellParamsA, La
            real(F64), dimension(((LmaxGTO+1)*(LmaxGTO+2))/2) :: C
            real(F64), dimension(2*LmaxGTO+1) :: S

            do k = 1, NVecs
                  do ShA = 1, NShells
                        ShellParamsA = ShellParamsIdx(ShA)
                        La = ShellMomentum(ShellParamsA)
                        NCart = ((La + 1) * (La + 2)) / 2
                        NSpher = 2 * La + 1
                        c0 = ShellLocCart(ShA)
                        c1 = c0 + NCart - 1
                        s0 = ShellLocSpher(ShA)
                        s1 = s0 + NSpher - 1
                        C(1:NCart) = ACart(c0:c1, k) * NormFactorsCart(1:NCart, ShellParamsA)
                        call Auto1e_SpherTransf_V(La)%ptr(S(1:NSpher), C(1:NCart))
                        ASpher(s0:s1, k) = S(1:NSpher) / NormFactorsSpher(1:NSpher, ShellParamsA)
                  end do
            end do
      end subroutine SpherGTO_TransformVectors


      subroutine SpherGTO_TransformVectors_U(ASpher, ACart, LmaxGTO, NormFactorsSpher, &
            NormFactorsCart, ShellLocSpher, ShellLocCart, ShellMomentum, ShellParamsIdx, &
            NAOSpher, NAOCart, NShells, NVecs)
            real(F64), dimension(NAOSpher, NVecs), intent(out) :: ASpher
            real(F64), dimension(NAOCart, NVecs), intent(in)   :: ACart
            integer, intent(in)                                :: LmaxGTO
            real(F64), dimension(:, :), intent(in)             :: NormFactorsSpher
            real(F64), dimension(:, :), intent(in)             :: NormFactorsCart
            integer, dimension(NShells), intent(in)            :: ShellLocSpher
            integer, dimension(NShells), intent(in)            :: ShellLocCart
            integer, dimension(:), intent(in)                  :: ShellMomentum
            integer, dimension(:), intent(in)                  :: ShellParamsIdx
            integer, intent(in)                                :: NAOSpher
            integer, intent(in)                                :: NAOCart
            integer, intent(in)                                :: NShells
            integer, intent(in)                                :: NVecs
            
            integer :: NSpher, NCart
            integer :: c0, c1, s0, s1, k
            integer :: ShA, ShellParamsA, La
            real(F64), dimension(((LmaxGTO+1)*(LmaxGTO+2))/2) :: C
            real(F64), dimension(2*LmaxGTO+1) :: S

            do k = 1, NVecs
                  do ShA = 1, NShells
                        ShellParamsA = ShellParamsIdx(ShA)
                        La = ShellMomentum(ShellParamsA)
                        NCart = ((La + 1) * (La + 2)) / 2
                        NSpher = 2 * La + 1
                        c0 = ShellLocCart(ShA)
                        c1 = c0 + NCart - 1
                        s0 = ShellLocSpher(ShA)
                        s1 = s0 + NSpher - 1
                        C(1:NCart) = ACart(c0:c1, k) / NormFactorsCart(1:NCart, ShellParamsA)
                        call Auto1e_SpherTransf_U(La)%ptr(S(1:NSpher), C(1:NCart))
                        ASpher(s0:s1, k) = S(1:NSpher) * NormFactorsSpher(1:NSpher, ShellParamsA)
                  end do
            end do
      end subroutine SpherGTO_TransformVectors_U
      

      subroutine SpherGTO_NormFactors(NormFactors, ShellMomentum, &
            NPrimitives, CntrCoeffs, Exponents, NShellParams)
            real(F64), dimension(:, :), intent(out) :: NormFactors
            integer, dimension(:), intent(in)       :: ShellMomentum
            integer, dimension(:), intent(in)       :: NPrimitives
            real(F64), dimension(:, :), intent(in)  :: CntrCoeffs
            real(F64), dimension(:, :), intent(in)  :: Exponents
            integer, intent(in)                     :: NShellParams

            integer :: k

            do k = 1, NShellParams
                  call SpherGTO_NormFactors_Shell(NormFactors(:, k), ShellMomentum(k), &
                        CntrCoeffs(:, k), Exponents(:, k), NPrimitives(k))
            end do
      end subroutine SpherGTO_NormFactors


      subroutine SpherGTO_NormFactors_Shell(NormFactors, La, CntrA, ExpA, NprimA)
            real(F64), dimension(:), intent(out) :: NormFactors
            integer, intent(in)                  :: La
            real(F64), dimension(*), intent(in)  :: CntrA, ExpA
            integer, intent(in)                  :: NprimA

            real(F64), dimension(((La+1)*(La+2))/2, ((La+1)*(La+2))/2) :: Tcc
            real(F64), dimension(((La+1)*(La+2))/2, 2*La+1) :: Tcs
            real(F64), dimension(2*La+1, 2*La+1) :: Tss
            integer :: i, NCart

            NCart = ((La + 1) * (La + 2)) / 2
            NormFactors = ZERO
            !
            ! Compute overlap integrals between Cartesian orbitals
            !
            call ShellOverlap(Tcc, La, CntrA, ExpA, NprimA)
            !
            ! Cartesian -> spherical transformation
            !
            call Auto1e_SpherTransf_U_Vector(La)%ptr(Tcs, Tcc, NCart)
            do i = 1, 2 * La + 1
                  call Auto1e_SpherTransf_U(La)%ptr(Tss(:, i), Tcs(:, i))
            end do
            do i = 1, 2 * La + 1
                  NormFactors(i) = ONE / sqrt(Tss(i, i))
            end do
      end subroutine SpherGTO_NormFactors_Shell

      
      subroutine ShellOverlap(S, La, CntrA, ExpA, NprimA)
            !
            ! Compute overlap integrals for a pair of Cartesian Gaussian shells
            ! of contracted atomic orbitals centered on the same atom. 
            !
            real(F64), dimension(((La+1)*(La+2))/2, *), intent(out) :: S
            integer, intent(in)                       :: La
            real(F64), dimension(*), intent(in)       :: CntrA, ExpA
            integer, intent(in)                       :: NprimA

            real(F64) :: AlphaA, AlphaB, AlphaP
            real(F64), dimension(0:2*La) :: OneDInts
            integer :: t
            integer :: k, l
            integer :: a, b
            integer :: Na
            integer :: lx, ly, lz
            integer, dimension(((La+1)*(La+2))/2) :: ax, ay, az

            Na = ((La + 1) * (La + 2)) / 2
            a = 1
            do lx = La, 0, -1
                  do ly = La-lx, 0, -1
                        lz = La - lx - ly
                        ax(a) = lx
                        ay(a) = ly
                        az(a) = lz
                        a = a + 1
                  end do
            end do
            OneDInts = ZERO
            S(:, 1:Na) = ZERO
            do l = 1, NprimA
                  do k = 1, NprimA
                        AlphaA = ExpA(k)
                        AlphaB = ExpA(l)
                        AlphaP = AlphaA + AlphaB
                        do t = 0, La
                              !
                              ! Integrate(-Inf,+Inf) x**(2*t) * exp(-AlphaP * x**2) dx
                              !
                              OneDInts(2*t) = dblfact(2*t-1) / (TWO * AlphaP)**t * Sqrt(PI/AlphaP)
                        end do
                        do b = 1, Na
                              do a = 1, Na
                                    S(a, b) = S(a, b) + &
                                          CntrA(k) * CntrA(l) & 
                                          * OneDInts(ax(a)+ax(b)) &
                                          * OneDInts(ay(a)+ay(b)) &
                                          * OneDInts(az(a)+az(b))
                              end do
                        end do
                  end do
            end do
      end subroutine ShellOverlap
end module SpherGTO
