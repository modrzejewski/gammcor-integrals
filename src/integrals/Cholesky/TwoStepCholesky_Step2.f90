module TwoStepCholesky_Step2
      use arithmetic
      use math_constants
      use Auto2e
      use real_linalg
      use clock
      use string
      use sort
      use basis_sets
      use TwoStepCholesky_definitions
      
      implicit none
      
contains

      subroutine chol2_write_Vpqrs(M, T, V, Na, Nb, Nc, Nd, Npq, Nrs, Ipq, Irs)
            real(F64), dimension(:, :), intent(out)        :: M
            real(F64), dimension(Npq, Nrs), intent(out)    :: T
            real(F64), dimension(Na*Nb, Nc*Nd), intent(in) :: V
            integer, intent(in)                            :: Na, Nb, Nc, Nd
            integer, intent(in)                            :: Npq, Nrs
            integer, dimension(:), intent(in)              :: Ipq, Irs

            integer :: pq, rs

            do rs = 1, Nrs
                  do pq = 1, Npq
                        T(pq, rs) = V(Ipq(pq), Irs(rs))
                  end do
            end do
            M = T
      end subroutine chol2_write_Vpqrs

      
      subroutine chol2_write_Vabrs(M, T, V, Na, Nb, Nc, Nd, Nrs, Irs)
            real(F64), dimension(:, :), intent(out)        :: M
            real(F64), dimension(Na*Nb, Nrs), intent(out)  :: T
            real(F64), dimension(Na*Nb, Nc*Nd), intent(in) :: V
            integer, intent(in)                            :: Na, Nb, Nc, Nd
            integer, intent(in)                            :: Nrs
            integer, dimension(:), intent(in)              :: Irs

            integer :: rs

            do rs = 1, Nrs
                  T(:, rs) = V(:, Irs(rs))
            end do
            M = T
      end subroutine chol2_write_Vabrs


      subroutine chol2_write_Vaars(M, T, V, Na, Nc, Nd, Nab, Nrs, Irs)
            real(F64), dimension(:, :), intent(out)         :: M
            real(F64), dimension(Nab, Nrs), intent(out)     :: T
            real(F64), dimension(Na, Na, Nc*Nd), intent(in) :: V
            integer, intent(in)                             :: Na, Nc, Nd
            integer, intent(in)                             :: Nab, Nrs
            integer, dimension(:), intent(in)               :: Irs

            integer :: rs
            integer :: a, b, q

            do rs = 1, Nrs
                  q = 1
                  do b = 1, Na
                        do a = b, Na
                              T(q, rs) = V(a, b, Irs(rs))
                              q = q + 1
                        end do
                  end do
            end do
            M = T
      end subroutine chol2_write_Vaars

      
      subroutine chol2_Vpqrs_FullInterface(Vpqrs, PivotShellPairs, PivotShellPairLoc, PivotShellPairDim, &
            PivotOrbPairs, NPivotShellPairs, ShellPairs, ShellCenters, AtomCoords, ShellParamsIdx, &
            ShellMomentum, NAngFunc, NPrimitives, CntrCoeffs, Exponents, NormFactors, Kappa, &
            MaxNAngFunc, PtrOffset)

            real(F64), dimension(:, :), intent(out)  :: Vpqrs
            integer, dimension(:), intent(in)        :: PivotShellPairs
            integer, dimension(:), intent(in)        :: PivotShellPairLoc
            integer, dimension(:), intent(in)        :: PivotShellPairDim
            integer, dimension(:), intent(in)        :: PivotOrbPairs
            integer, intent(in)                      :: NPivotShellPairs
            integer, dimension(:, :), intent(in)     :: ShellPairs
            integer, dimension(:), intent(in)        :: ShellCenters
            real(F64), dimension(:, :), intent(in)   :: AtomCoords
            integer, dimension(:), intent(in)        :: ShellParamsIdx
            integer, dimension(:), intent(in)        :: ShellMomentum
            integer, dimension(:), intent(in)        :: NAngFunc
            integer, dimension(:), intent(in)        :: NPrimitives
            real(F64), dimension(:, :), intent(in)   :: CntrCoeffs
            real(F64), dimension(:, :), intent(in)   :: Exponents
            real(F64), dimension(:, :), intent(in)   :: NormFactors
            real(F64), intent(in)                    :: Kappa
            integer, intent(in)                      :: MaxNAngFunc
            integer, intent(in)                      :: PtrOffset

            integer :: KL, K, L, AtomA, AtomB, AtomC, AtomD
            integer :: ShellParamsA, ShellParamsB, ShellParamsC, ShellParamsD
            integer :: ShAB, ShA, ShB, pq0, pq1, Npq
            integer :: ShCD, ShC, ShD, rs0, rs1, Nrs
            integer :: La, Lb, Lc, Ld
            integer :: Na, Nb, Nc, Nd
            real(F64), dimension(MaxNAngFunc**4) :: V, T

            Vpqrs = ZERO
            !$omp parallel do schedule(guided) &
            !$omp private(L, K, KL) &
            !$omp private(ShAB, ShCD, ShA, ShB, ShC, ShD, AtomA, AtomB, AtomC, AtomD) &
            !$omp private(ShellParamsA, ShellParamsB, ShellParamsC, ShellParamsD) &
            !$omp private(La, Lb, Lc, Ld, Na, Nb, Nc, Nd, Npq, Nrs, pq0, pq1, rs0, rs1) &
            !$omp private(V, T) &
            !$omp default(shared)
            do KL = 1, (NPivotShellPairs * (NPivotShellPairs + 1)) / 2
                  call chol2_pq2p_ge_q(K, L, KL, NPivotShellPairs)
                  ShAB = PivotShellPairs(K)
                  ShCD = PivotShellPairs(L)

                  ShA = ShellPairs(1, ShAB)
                  ShB = ShellPairs(2, ShAB)
                  ShC = ShellPairs(1, ShCD)
                  ShD = ShellPairs(2, ShCD)

                  AtomA = ShellCenters(ShA)
                  AtomB = ShellCenters(ShB)
                  AtomC = ShellCenters(ShC)
                  AtomD = ShellCenters(ShD)

                  ShellParamsB = ShellParamsIdx(ShB)
                  ShellParamsA = ShellParamsIdx(ShA)
                  ShellParamsC = ShellParamsIdx(ShC)
                  ShellParamsD = ShellParamsIdx(ShD)

                  La = ShellMomentum(ShellParamsA)
                  Lb = ShellMomentum(ShellParamsB)
                  Lc = ShellMomentum(ShellParamsC)
                  Ld = ShellMomentum(ShellParamsD)

                  Na = NAngFunc(ShellParamsA)
                  Nb = NAngFunc(ShellParamsB)
                  Nc = NAngFunc(ShellParamsC)                        
                  Nd = NAngFunc(ShellParamsD)

                  Npq = PivotShellPairDim(K)
                  Nrs = PivotShellPairDim(L)

                  pq0 = PivotShellPairLoc(K)
                  pq1 = PivotShellPairLoc(K) + Npq - 1
                  rs0 = PivotShellPairLoc(L)
                  rs1 = PivotShellPairLoc(L) + Nrs - 1

                  call Auto2eERI(PtrOffset+auto2e_idx(Ld, Lc, Lb, La))%ptr( &
                        V, &
                        !
                        ! ShellD
                        !
                        AtomCoords(:, AtomD), CntrCoeffs(:, ShellParamsD), &
                        NormFactors(:, ShellParamsD), Exponents(:, ShellParamsD), &
                        NPrimitives(ShellParamsD), &
                        !
                        ! ShellC
                        !
                        AtomCoords(:, AtomC), CntrCoeffs(:, ShellParamsC), &
                        NormFactors(:, ShellParamsC), Exponents(:, ShellParamsC), &
                        NPrimitives(ShellParamsC), &
                        !
                        ! ShellB
                        !
                        AtomCoords(:, AtomB), CntrCoeffs(:, ShellParamsB), &
                        NormFactors(:, ShellParamsB), Exponents(:, ShellParamsB), &
                        NPrimitives(ShellParamsB), &
                        !
                        ! ShellA
                        !
                        AtomCoords(:, AtomA), CntrCoeffs(:, ShellParamsA), &
                        NormFactors(:, ShellParamsA), Exponents(:, ShellParamsA), &
                        NPrimitives(ShellParamsA), &
                        Kappa)

                  call chol2_write_Vpqrs(Vpqrs(pq0:pq1, rs0:rs1), T, V, Na, Nb, Nc, Nd, Npq, Nrs, &
                        PivotOrbPairs(pq0:pq1), PivotOrbPairs(rs0:rs1))
            end do
            !$omp end parallel do
      end subroutine chol2_Vpqrs_FullInterface


      subroutine chol2_Vpqrs(Vpqrs, PivotShellPairs, PivotShellPairLoc, PivotShellPairDim, &
            PivotOrbPairs, NPivotShellPairs, ShellPairs, AOBasis, Kappa)

            real(F64), dimension(:, :), intent(out)  :: Vpqrs
            integer, dimension(:), intent(in)        :: PivotShellPairs
            integer, dimension(:), intent(in)        :: PivotShellPairLoc
            integer, dimension(:), intent(in)        :: PivotShellPairDim
            integer, dimension(:), intent(in)        :: PivotOrbPairs
            integer, intent(in)                      :: NPivotShellPairs
            integer, dimension(:, :), intent(in)     :: ShellPairs
            type(TAOBasis), intent(in)               :: AOBasis
            real(F64), intent(in)                    :: Kappa

            integer :: PtrOffset, MaxNAngFunc
            
            associate ( &
                  LmaxGTO => AOBasis%LmaxGTO, &
                  SpherAO => AOBasis%SpherAO, &
                  ShellCenters => AOBasis%ShellCenters, &
                  AtomCoords => AOBasis%AtomCoords, &
                  ShellParamsIdx => AOBasis%ShellParamsIdx, &
                  ShellMomentum => AOBasis%ShellMomentum, &
                  NPrimitives => AOBasis%NPrimitives, &
                  CntrCoeffs => AOBasis%CntrCoeffs, &
                  Exponents => AOBasis%Exponents, &
                  NAngFuncSpher => AOBasis%NAngFuncSpher, &
                  NAngFuncCart => AOBasis%NAngFuncCart, &
                  NormFactorsSpher => AOBasis%NormFactorsSpher, &
                  NormFactorsCart => AOBasis%NormFactorsCart &
                  )
                  if (SpherAO) then
                        PtrOffset = AUTO2E_SPHER_OFFSET
                        MaxNAngFunc = 2 * LmaxGTO + 1
                  else
                        PtrOffset = 0
                        MaxNAngFunc = ((LmaxGTO + 1) * (LmaxGTO + 2)) / 2
                  end if
                  if (SpherAO) then
                        call chol2_Vpqrs_FullInterface(Vpqrs, PivotShellPairs, PivotShellPairLoc, PivotShellPairDim, &
                              PivotOrbPairs, NPivotShellPairs, ShellPairs, ShellCenters, AtomCoords, ShellParamsIdx, &
                              ShellMomentum, NAngFuncSpher, NPrimitives, CntrCoeffs, Exponents, NormFactorsSpher, Kappa, &
                              MaxNAngFunc, PtrOffset)
                  else
                        call chol2_Vpqrs_FullInterface(Vpqrs, PivotShellPairs, PivotShellPairLoc, PivotShellPairDim, &
                              PivotOrbPairs, NPivotShellPairs, ShellPairs, ShellCenters, AtomCoords, ShellParamsIdx, &
                              ShellMomentum, NAngFuncCart, NPrimitives, CntrCoeffs, Exponents, NormFactorsCart, Kappa, &
                              MaxNAngFunc, PtrOffset)
                  end if
            end associate
      end subroutine chol2_Vpqrs


      subroutine chol2_Vabrs_FullInterface(Vabrs, SubsetBounds, PivotShellPairs, PivotShellPairLoc, PivotShellPairDim, &
            PivotOrbPairs, NPivotShellPairs, ShellPairLoc, ShellPairs, ShellPairDim, &
            ShellCenters, AtomCoords, ShellParamsIdx, ShellMomentum, NAngFunc, &
            NPrimitives, CntrCoeffs, Exponents, NormFactors, Kappa, MaxNAngFunc, PtrOffset)

            real(F64), dimension(:, :), intent(out)  :: Vabrs
            integer, dimension(2), intent(in)        :: SubsetBounds
            integer, dimension(:), intent(in)        :: PivotShellPairs
            integer, dimension(:), intent(in)        :: PivotShellPairLoc
            integer, dimension(:), intent(in)        :: PivotShellPairDim
            integer, dimension(:), intent(in)        :: PivotOrbPairs
            integer, intent(in)                      :: NPivotShellPairs
            integer, dimension(:, :), intent(in)     :: ShellPairLoc
            integer, dimension(:, :), intent(in)     :: ShellPairs
            integer, dimension(:), intent(in)        :: ShellPairDim
            integer, dimension(:), intent(in)        :: ShellCenters
            real(F64), dimension(:, :), intent(in)   :: AtomCoords
            integer, dimension(:), intent(in)        :: ShellParamsIdx
            integer, dimension(:), intent(in)        :: ShellMomentum
            integer, dimension(:), intent(in)        :: NAngFunc
            integer, dimension(:), intent(in)        :: NPrimitives
            real(F64), dimension(:, :), intent(in)   :: CntrCoeffs
            real(F64), dimension(:, :), intent(in)   :: Exponents
            real(F64), dimension(:, :), intent(in)   :: NormFactors
            real(F64), intent(in)                    :: Kappa
            integer, intent(in)                      :: MaxNAngFunc
            integer, intent(in)                      :: PtrOffset

            integer :: L, AtomA, AtomB, AtomC, AtomD
            integer :: ShellParamsA, ShellParamsB, ShellParamsC, ShellParamsD
            integer :: ShAB, ShA, ShB, ab0, ab1, Nab
            integer :: ShCD, ShC, ShD, rs0, rs1, Nrs
            integer :: La, Lb, Lc, Ld
            integer :: Na, Nb, Nc, Nd
            real(F64), dimension(MaxNAngFunc**4) :: V
            real(F64), dimension(MaxNAngFunc**4) :: T

            Vabrs = ZERO
            !$omp parallel do collapse(2) schedule(guided) &
            !$omp private(L) &
            !$omp private(ShAB, ShCD, ShA, ShB, ShC, ShD, AtomA, AtomB, AtomC, AtomD) &
            !$omp private(ShellParamsA, ShellParamsB, ShellParamsC, ShellParamsD) &
            !$omp private(La, Lb, Lc, Ld, Na, Nb, Nc, Nd, ab0, ab1, rs0, rs1) &
            !$omp private(Nab, Nrs) &
            !$omp private(V, T) &
            !$omp default(shared)
            do L = 1, NPivotShellPairs
                  do ShAB = SubsetBounds(1), SubsetBounds(2)
                        ShCD = PivotShellPairs(L)

                        ShA = ShellPairs(1, ShAB)
                        ShB = ShellPairs(2, ShAB)
                        ShC = ShellPairs(1, ShCD)
                        ShD = ShellPairs(2, ShCD)

                        AtomA = ShellCenters(ShA)
                        AtomB = ShellCenters(ShB)
                        AtomC = ShellCenters(ShC)
                        AtomD = ShellCenters(ShD)

                        ShellParamsB = ShellParamsIdx(ShB)
                        ShellParamsA = ShellParamsIdx(ShA)
                        ShellParamsC = ShellParamsIdx(ShC)
                        ShellParamsD = ShellParamsIdx(ShD)

                        La = ShellMomentum(ShellParamsA)
                        Lb = ShellMomentum(ShellParamsB)
                        Lc = ShellMomentum(ShellParamsC)
                        Ld = ShellMomentum(ShellParamsD)

                        Na = NAngFunc(ShellParamsA)
                        Nb = NAngFunc(ShellParamsB)
                        Nc = NAngFunc(ShellParamsC)                        
                        Nd = NAngFunc(ShellParamsD)

                        Nab = ShellPairDim(ShAB)
                        Nrs = PivotShellPairDim(L)

                        ab0 = ShellPairLoc(CHOL2_SUBSET_STORAGE, ShAB)
                        ab1 = ShellPairLoc(CHOL2_SUBSET_STORAGE, ShAB) + Nab - 1
                        rs0 = PivotShellPairLoc(L)
                        rs1 = PivotShellPairLoc(L) + Nrs - 1

                        call Auto2eERI(PtrOffset+auto2e_idx(Ld, Lc, Lb, La))%ptr( &
                              V, &
                              !
                              ! ShellD
                              !
                              AtomCoords(:, AtomD), CntrCoeffs(:, ShellParamsD), &
                              NormFactors(:, ShellParamsD), Exponents(:, ShellParamsD), &
                              NPrimitives(ShellParamsD), &
                              !
                              ! ShellC
                              !
                              AtomCoords(:, AtomC), CntrCoeffs(:, ShellParamsC), &
                              NormFactors(:, ShellParamsC), Exponents(:, ShellParamsC), &
                              NPrimitives(ShellParamsC), &
                              !
                              ! ShellB
                              !
                              AtomCoords(:, AtomB), CntrCoeffs(:, ShellParamsB), &
                              NormFactors(:, ShellParamsB), Exponents(:, ShellParamsB), &
                              NPrimitives(ShellParamsB), &
                              !
                              ! ShellA
                              !
                              AtomCoords(:, AtomA), CntrCoeffs(:, ShellParamsA), &
                              NormFactors(:, ShellParamsA), Exponents(:, ShellParamsA), &
                              NPrimitives(ShellParamsA), &
                              Kappa)

                        if (ShA /= ShB) then
                              call chol2_write_Vabrs(Vabrs(ab0:ab1, rs0:rs1), T, V, &
                                    Na, Nb, Nc, Nd, Nrs, PivotOrbPairs(rs0:rs1))
                        else
                              call chol2_write_Vaars(Vabrs(ab0:ab1, rs0:rs1), T, V, &
                                    Na, Nc, Nd, Nab, Nrs, PivotOrbPairs(rs0:rs1))
                        end  if
                  end do
            end do
            !$omp end parallel do
      end subroutine chol2_Vabrs_FullInterface

      
      subroutine chol2_Vabrs(Vabrs, SubsetBounds, PivotShellPairs, PivotShellPairLoc, PivotShellPairDim, &
            PivotOrbPairs, NPivotShellPairs, ShellPairLoc, ShellPairs, ShellPairDim, AOBasis, Kappa)

            real(F64), dimension(:, :), intent(out)  :: Vabrs
            integer, dimension(2), intent(in)        :: SubsetBounds
            integer, dimension(:), intent(in)        :: PivotShellPairs
            integer, dimension(:), intent(in)        :: PivotShellPairLoc
            integer, dimension(:), intent(in)        :: PivotShellPairDim
            integer, dimension(:), intent(in)        :: PivotOrbPairs
            integer, intent(in)                      :: NPivotShellPairs
            integer, dimension(:, :), intent(in)     :: ShellPairLoc
            integer, dimension(:, :), intent(in)     :: ShellPairs
            integer, dimension(:), intent(in)        :: ShellPairDim
            type(TAOBasis), intent(in)               :: AOBasis
            real(F64), intent(in)                    :: Kappa
            
            integer :: PtrOffset, MaxNAngFunc
            
            associate ( &
                  LmaxGTO => AOBasis%LmaxGTO, &
                  SpherAO => AOBasis%SpherAO, &
                  ShellCenters => AOBasis%ShellCenters, &
                  AtomCoords => AOBasis%AtomCoords, &
                  ShellParamsIdx => AOBasis%ShellParamsIdx, &
                  ShellMomentum => AOBasis%ShellMomentum, &
                  NPrimitives => AOBasis%NPrimitives, &
                  CntrCoeffs => AOBasis%CntrCoeffs, &
                  Exponents => AOBasis%Exponents, &
                  NAngFuncSpher => AOBasis%NAngFuncSpher, &
                  NAngFuncCart => AOBasis%NAngFuncCart, &
                  NormFactorsSpher => AOBasis%NormFactorsSpher, &
                  NormFactorsCart => AOBasis%NormFactorsCart &
                  )
                  if (SpherAO) then
                        PtrOffset = AUTO2E_SPHER_OFFSET
                        MaxNAngFunc = 2 * LmaxGTO + 1
                  else
                        PtrOffset = 0
                        MaxNAngFunc = ((LmaxGTO + 1) * (LmaxGTO + 2)) / 2
                  end if

                  if (SpherAO) then
                        call chol2_Vabrs_FullInterface(Vabrs, SubsetBounds, PivotShellPairs, PivotShellPairLoc, PivotShellPairDim, &
                              PivotOrbPairs, NPivotShellPairs, ShellPairLoc, ShellPairs, ShellPairDim, &
                              ShellCenters, AtomCoords, ShellParamsIdx, ShellMomentum, NAngFuncSpher, &
                              NPrimitives, CntrCoeffs, Exponents, NormFactorsSpher, Kappa, MaxNAngFunc, PtrOffset)
                  else
                        call chol2_Vabrs_FullInterface(Vabrs, SubsetBounds, PivotShellPairs, PivotShellPairLoc, PivotShellPairDim, &
                              PivotOrbPairs, NPivotShellPairs, ShellPairLoc, ShellPairs, ShellPairDim, &
                              ShellCenters, AtomCoords, ShellParamsIdx, ShellMomentum, NAngFuncCart, &
                              NPrimitives, CntrCoeffs, Exponents, NormFactorsCart, Kappa, MaxNAngFunc, PtrOffset)
                  end if
            end associate
      end subroutine chol2_Vabrs
      
      
      subroutine chol2_pq2p_ge_q(p, q, pq, n)
            !
            ! Decode a lower-triangle compound index into individual
            ! indices:
            ! PQ -> (P, Q)
            ! Assumptions:
            ! 0) P = 1, 2, ..., N,
            !    Q = 1, 2, ..., N,
            ! 1) P >= Q (diagonal indices admissible)
            !
            ! An example of how this algorithm traverses an N=3 triangle:
            !
            ! 1
            ! 2 5
            ! 3 6 4
            !
            integer, intent(out) :: p
            integer, intent(out) :: q
            integer, intent(in)  :: pq
            integer, intent(in)  :: n

            integer :: q_base
            integer :: v
            integer :: interval1
            integer :: in1, in2
            !
            ! pq = (q_base - 1) * (n + 1) + v
            !
            q_base = (pq - 1) / (n + 1) + 1
            v = pq - (n + 1) * (q_base - 1)
            !
            ! Decide if v is in interval_1 or interval_2:
            ! in1 == 1 and in2 == 0 if v <= INTERVAL1
            ! in1 == 0 and in2 == 1 if v > INTERVAL1
            !
            interval1 = n - q_base + 1
            in2 = v / (interval1 + 1)
            !
            ! 1 -> 0, 0 -> 1
            !
            in1 = ieor((in2), 1)

            p = in1 * (q_base + v - 1) + in2 * (v - interval1 + n - q_base)          
            q = in1 * q_base + in2 * interval1
      end subroutine chol2_pq2p_ge_q
      

      subroutine chol2_PivotShellPairs(PivotShellPairs, PivotShellPairLoc, PivotShellPairDim, &
            PivotOrbPairs, NPivotShellPairs, Pivots, ShellPairs, ShellPairDim, ShellPairLoc, &
            ShellParamsIdx, NAngFunc, NShellPairs, NCholesky)
            
            integer, dimension(:), intent(out)       :: PivotShellPairs
            integer, dimension(:), intent(out)       :: PivotShellPairLoc
            integer, dimension(:), intent(out)       :: PivotShellPairDim
            integer, dimension(:), intent(out)       :: PivotOrbPairs
            integer, intent(out)                     :: NPivotShellPairs
            integer, dimension(:), intent(in)        :: Pivots
            integer, dimension(:, :), intent(in)     :: ShellPairs
            integer, dimension(:), intent(in)        :: ShellPairDim
            integer, dimension(:, :), intent(in)     :: ShellPairLoc
            integer, dimension(:), intent(in)        :: ShellParamsIdx
            integer, dimension(:), intent(in)        :: NAngFunc
            integer, intent(in)                      :: NShellPairs  
            integer, intent(in)                      :: NCholesky

            integer :: p, p0, p1, q, Nab, s
            integer :: ShAB, ShA, ShB
            integer :: ab, a, b
            integer :: Na, ShellParamsA
            integer :: N

            N = sum(Pivots)
            if (N /= NCholesky) then
                  call msg("Invalid Pivots array", MSG_ERROR)                  
                  error stop
            end if
            NPivotShellPairs = 0
            q = 0
            do ShAB = 1, NShellPairs
                  ShA = ShellPairs(1, ShAB)
                  ShB = ShellPairs(2, ShAB)
                  Nab = ShellPairDim(ShAB)
                  p0 = ShellPairLoc(CHOL2_FULL_STORAGE, ShAB)
                  p1 = ShellPairLoc(CHOL2_FULL_STORAGE, ShAB) + Nab - 1
                  s = sum(Pivots(p0:p1))
                  if (s > 0) then
                        NPivotShellPairs = NPivotShellPairs + 1
                        PivotShellPairs(NPivotShellPairs) = ShAB
                        PivotShellPairLoc(NPivotShellPairs) = q + 1
                        PivotShellPairDim(NPivotShellPairs) = s
                        if (ShA /= ShB) then
                              do p = p0, p1
                                    if (Pivots(p) == 1) then
                                          q = q + 1
                                          ab = p - p0 + 1
                                          PivotOrbPairs(q) = ab
                                    end if
                              end do
                        else
                              ShellParamsA = ShellParamsIdx(ShA)
                              Na = NAngFunc(ShellParamsA)
                              do p = p0, p1
                                    if (Pivots(p) == 1) then
                                          q = q + 1
                                          ab = p - p0 + 1
                                          call chol2_pq2p_ge_q(a, b, ab, Na)
                                          PivotOrbPairs(q) = a + (b - 1) * Na
                                    end if
                              end do
                        end if
                  end if
            end do
      end subroutine chol2_PivotShellPairs
end module TwoStepCholesky_Step2
