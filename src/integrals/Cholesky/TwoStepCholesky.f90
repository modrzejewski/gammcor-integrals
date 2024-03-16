module TwoStepCholesky
      ! -------------------------------------------------------------------------------------
      ! Two-step integral-direct Cholesky decomposition of the Coulomb matrix
      ! -------------------------------------------------------------------------------------
      ! Compute the Cholesky factorization V=R**T*R of the Coulomb matrix V(pq,rs).
      ! The Cholesky vectors,
      !
      !                 R(k,pq) for k=1, ..., NCholesky, p>=q
      !                 decomposition threshold: diagonal Dpq > Tau
      !
      ! include only the AO index pairs pq which passed the initial Schwarz inequality test.
      ! This subset of indices is referred to as the BASE subset throughout this program.
      !
      ! The implementation is based on the two step algorithm based on Ref. 2, which is
      ! an improved version of the single step algorithm described in detail in Ref. 1.
      ! In the first step, only the pivots are computed without full Rk's.
      ! The pivot AO pairs are then used in the second step by applying an efficient,
      ! non-pivoted Cholesky subroutine to a reduced-size Coulomb matrix of
      ! dimension (NCholesky, NCholesky).
      !
      ! 1. Aquilante et al. Subsection 13.5, Cholesky Decomposition Techniques in Electronic
      !    Structure Theory in Linear-Scaling Techniques in Computational Chemistry and
      !    Physics: Methods and Applications, 301-343, Springer 2011;
      !    doi: 10.1007/978-90-481-2853-2_13
      ! 2. Sarai D. Folkestad, Eirik F. Kjønstad and Henrik Koch,
      !    J. Chem. Phys. 150, 194112 (2019);
      !    doi: 10.1063/1.5083802
      !
      use arithmetic
      use math_constants
      use Auto2e
      use real_linalg
      use clock
      use string
      use sort
      use basis_sets
      use TwoStepCholesky_definitions
      use TwoStepCholesky_Step1
      use TwoStepCholesky_Step2

      implicit none

contains

      subroutine chol2_Algo_Koch_JCP2019(Chol2Vecs, AOBasis, Chol2Params)
            type(TChol2Vecs), intent(out)  :: Chol2Vecs
            type(TAOBasis), intent(in)     :: AOBasis
            type(TChol2Params), intent(in) :: Chol2Params
            !
            ! Generate the pivots of the Cholesky factorization V = R**T * R.
            !
            ! Small diagonal Coulomb integrals are removed from the computations
            ! at the prescreening phase using the Schwarz inequality. The indices (ShA, ShB)
            ! of the remaining shell pairs AB=1,2,...,NShellPairs are stored in ShellPairs(AB).
            ! The location of the first matrix element of each shell pair is stored in ShellPairLoc(:, AB).
            ! NOrbPairs is the total number of orbital pairs remaining after the prescreening step,
            ! which is different from NAO*(NAO+1)/2.
            !
            ! The algorithm for the pivot generation is based on Refs. 1 and 2.
            !
            ! 1. S.D. Folkestad, E.F. Kjonstad, and H. Koch.,
            !    An efficient algorithm for Cholesky decomposition of electron
            !    repulsion integrals, J. Chem. Phys. 150, 194112 (2019);
            !    doi: 10.1063/1.5083802
            ! 2. Aquilante et al. Subsection 13.5, Cholesky Decomposition Techniques in Electronic Structure Theory in
            !    Linear-Scaling Techniques in Computational Chemistry and Physics: Methods and Applications,
            !    301-343, Springer 2011; doi: 10.1007/978-90-481-2853-2_13
            !
            call chol2_Step1(Chol2Vecs, AOBasis, Chol2Params)
            call chol2_Step2(Chol2Vecs, AOBasis, Chol2Params)
      end subroutine chol2_Algo_Koch_JCP2019
      

      subroutine chol2_Step1_FullInterface(Pivots, NCholesky, ShellPairs, ShellPairLoc, ShellPairDim, NShellPairs, &
            NOrbPairs, SubsetDim, SubsetBounds, NSubsets, TauThresh, ShellCenters, AtomCoords, ShellParamsIdx, &
            ShellMomentum, NAngFunc, NPrimitives, CntrCoeffs, Exponents, NormFactors, Kappa, LmaxGTO, &
            NAO, NShells, SpherAO, OrbPairsBlock, CholVecsBlock)

            integer, dimension(:), allocatable, intent(out)             :: Pivots
            integer, intent(out)                                        :: NCholesky
            integer, dimension(2, (NShells*(NShells+1))/2), intent(out) :: ShellPairs
            integer, dimension(3, (NShells*(NShells+1))/2), intent(out) :: ShellPairLoc
            integer, dimension((NShells*(NShells+1))/2), intent(out)    :: ShellPairDim
            integer, intent(out)                                        :: NShellPairs
            integer, intent(out)                                        :: NOrbPairs
            integer, dimension(:), allocatable, intent(out)             :: SubsetDim
            integer, dimension(:, :), allocatable, intent(out)          :: SubsetBounds
            integer, dimension(2), intent(out)                          :: NSubsets
            real(F64), intent(in)                                       :: TauThresh
            integer, dimension(NShells), intent(in)                     :: ShellCenters
            real(F64), dimension(:, :), intent(in)                      :: AtomCoords
            integer, dimension(:), intent(in)                           :: ShellParamsIdx
            integer, dimension(:), intent(in)                           :: ShellMomentum
            integer, dimension(:), intent(in)                           :: NAngFunc
            integer, dimension(:), intent(in)                           :: NPrimitives
            real(F64), dimension(:, :), intent(in)                      :: CntrCoeffs
            real(F64), dimension(:, :), intent(in)                      :: Exponents
            real(F64), dimension(:, :), intent(in)                      :: NormFactors
            real(F64), intent(in)                                       :: Kappa
            integer, intent(in)                                         :: LmaxGTO
            integer, intent(in)                                         :: NAO
            integer, intent(in)                                         :: NShells
            logical, intent(in)                                         :: SpherAO
            integer, intent(in)                                         :: OrbPairsBlock
            integer, intent(in)                                         :: CholVecsBlock

            real(F64) :: PrescreenError
            integer :: PtrOffset
            integer :: MaxNOrbPairs, MaxNCholesky, MaxNAngFunc, MaxSubsetDim
            integer, dimension(:), allocatable :: NShellPairs0, NOrbPairs0
            integer, dimension(:, :), allocatable :: ShellPairLoc0
            real(F64), dimension(:), allocatable :: Vdiag
            real(F64), dimension(:), allocatable :: D
            integer :: ThisImage
            !
            ! Maximum number of Cholesky vectors expressed as a multiply
            ! of the number of AOs
            !
            real(F64), parameter :: MaxNAOMult = 15.0_F64

            ThisImage = this_image()
            !
            ! Select the set of two-electron integral subroutines:
            ! Cartesian or spherical AOs
            !
            if (SpherAO) then
                  PtrOffset = AUTO2E_SPHER_OFFSET
                  MaxNAngFunc = 2 * LmaxGTO + 1
            else
                  PtrOffset = 0
                  MaxNAngFunc = ((LmaxGTO + 1) * (LmaxGTO + 2)) / 2
            end if
            !
            ! Complete all diagonal two-electron itegrals, without prescreening
            !
            MaxNOrbPairs = (NAO*(NAO+1))/2
            allocate(Vdiag(MaxNOrbPairs))
            call chol2_Vdiag(Vdiag, ShellCenters, AtomCoords, ShellParamsIdx, &
                  ShellMomentum, NAngFunc, NPrimitives, CntrCoeffs, Exponents, NormFactors, Kappa, &
                  MaxNAngFunc, NAO, NShells, PtrOffset)
            !
            ! Use the Schwarz inequality to prescreen orbital pairs related to small integrals.
            ! Shell pairs comprised of small integrals will be removed from the subsequent
            ! storage and computational steps
            !
            allocate(NShellPairs0(4))
            allocate(NOrbPairs0(4))
            allocate(ShellPairLoc0(4, (NShells*(NShells+1))/2))
            call chol2_DefineBase(ShellPairs, ShellPairLoc0, ShellPairDim, NShellPairs0, NOrbPairs0, D, PrescreenError, &
                  TauThresh, Vdiag, ShellParamsIdx, ShellMomentum, NAngFunc, NShells, NAO)
            call chol2_DefineSubsets(SubsetDim, SubsetBounds, ShellPairLoc0, NSubsets, &
                  ShellPairDim, NShellPairs0, NOrbPairs0, OrbPairsBlock)
            NShellPairs = NShellPairs0(CHOL2_BASE)
            NOrbPairs = NOrbPairs0(CHOL2_BASE)
            MaxNCholesky = min(ceiling(NAO * MaxNAOMult), NOrbPairs)
            MaxSubsetDim = maxval(SubsetDim)
            allocate(Pivots(NOrbPairs))
            if (ThisImage == 1) then
                  call chol2_Pivots(Pivots, NCholesky, NShellPairs0, NOrbPairs0, &
                        D, MaxNCholesky, TauThresh, &
                        ShellPairs, ShellPairLoc0, ShellPairDim, ShellCenters, &
                        AtomCoords, ShellParamsIdx, ShellMomentum, NAngFunc, NPrimitives, CntrCoeffs, Exponents, &
                        NormFactors, Kappa, MaxNAngFunc, NAO, PtrOffset, CholVecsBlock)
            end if
            ShellPairLoc(:, :) = ShellPairLoc0(1:3, :)
            call co_broadcast(Pivots, source_image=1)
      end subroutine chol2_Step1_FullInterface


      subroutine chol2_Step1(Chol2Vecs, AOBasis, Chol2Params)
            type(TChol2Vecs), intent(out)  :: Chol2Vecs
            type(TAOBasis), intent(in)     :: AOBasis
            type(TChol2Params), intent(in) :: Chol2Params

            integer :: MaxNShellPairs
            
            MaxNShellPairs = (AOBasis%NShells * (AOBasis%NShells + 1)) / 2
            allocate(Chol2Vecs%ShellPairs(2, MaxNShellPairs))
            allocate(Chol2Vecs%ShellPairLoc(3, MaxNShellPairs))
            allocate(Chol2Vecs%ShellPairDim(MaxNShellPairs))
            if (AOBasis%SpherAO) then
                  call chol2_Step1_FullInterface( &
                        Chol2Vecs%Pivots, &
                        Chol2Vecs%NVecs, &
                        Chol2Vecs%ShellPairs, &
                        Chol2Vecs%ShellPairLoc, &
                        Chol2Vecs%ShellPairDim, &
                        Chol2Vecs%NShellPairs, &
                        Chol2Vecs%NOrbPairs, &
                        Chol2Vecs%SubsetDim, &
                        Chol2Vecs%SubsetBounds, &
                        Chol2Vecs%NSubsets, &
                        Chol2Params%CholeskyTauThresh, & 
                        AOBasis%ShellCenters, &
                        AOBasis%AtomCoords, &
                        AOBasis%ShellParamsIdx, &
                        AOBasis%ShellMomentum, &
                        AOBasis%NAngFuncSpher, &
                        AOBasis%NPrimitives, &
                        AOBasis%CntrCoeffs, &
                        AOBasis%Exponents, &
                        AOBasis%NormFactorsSpher, &
                        Chol2Params%Kappa, &
                        AOBasis%LmaxGTO, &
                        AOBasis%NAOSpher, &
                        AOBasis%NShells, &
                        AOBasis%SpherAO, &
                        Chol2Params%MaxBlockDim, &
                        Chol2Params%CholVecsBlock)
            else
                  call chol2_Step1_FullInterface( &
                        Chol2Vecs%Pivots, &
                        Chol2Vecs%NVecs, &
                        Chol2Vecs%ShellPairs, &
                        Chol2Vecs%ShellPairLoc, &
                        Chol2Vecs%ShellPairDim, &
                        Chol2Vecs%NShellPairs, &
                        Chol2Vecs%NOrbPairs, &
                        Chol2Vecs%SubsetDim, &
                        Chol2Vecs%SubsetBounds, &
                        Chol2Vecs%NSubsets, &
                        Chol2Params%CholeskyTauThresh, &
                        AOBasis%ShellCenters, &
                        AOBasis%AtomCoords, &
                        AOBasis%ShellParamsIdx, &
                        AOBasis%ShellMomentum, &
                        AOBasis%NAngFuncCart, &
                        AOBasis%NPrimitives, &
                        AOBasis%CntrCoeffs, &
                        AOBasis%Exponents, &
                        AOBasis%NormFactorsCart, &
                        Chol2Params%Kappa, &
                        AOBasis%LmaxGTO, &
                        AOBasis%NAOCart, &
                        AOBasis%NShells, &
                        AOBasis%SpherAO, &
                        Chol2Params%MaxBlockDim, &
                        Chol2Params%CholVecsBlock)
            end if
      end subroutine chol2_Step1


      subroutine chol2_Step2(Chol2Vecs, AOBasis, Chol2Params)
            !
            ! Perform the second step of the pivoted Cholesky decomposition
            ! according to the algorithm of Koch et al.
            !
            ! 1. Build a reduced-dimension Coulomb matrix V(pq,rs) where
            !    pq, rs are pivot AO pairs obtained after step 1.
            ! 2. Factorize V using non-pivoted Cholesky decomposistion
            !                    V = L * L**T
            !    The dimension of L is (NCholesky, NCholesky).
            ! 3. Invert the lower triangular matrix L.  
            !
            ! From this point, L**(-1) can be used to generate full-dimension
            ! Cholesky vectors using Eq. 3 of Ref. 1.
            !
            ! 1. S.D. Folkestad, E.F. Kjonstad, and H. Koch.,
            !    An efficient algorithm for Cholesky decomposition of electron
            !    repulsion integrals, J. Chem. Phys. 150, 194112 (2019);
            !    doi: 10.1063/1.5083802
            !
            type(TChol2Vecs), intent(inout) :: Chol2Vecs
            type(TAOBasis), intent(in)      :: AOBasis
            type(TChol2Params), intent(in)  :: Chol2Params

            real(F64), dimension(:, :), allocatable :: Vpqrs
            type(TClock) :: timer, DeltaT
            real(F64) :: t_Vpqrs, t_Cholesky, t_Inverse, t_Total
            
            associate ( &
                  NShellPairs => Chol2Vecs%NShellPairs, &
                  NOrbPairs => Chol2Vecs%NOrbPairs, &
                  NCholesky => Chol2Vecs%NVecs, &
                  ShellPairDim => Chol2Vecs%ShellPairDim, &
                  ShellPairLoc => Chol2Vecs%ShellPairLoc, &
                  ShellPairs => Chol2Vecs%ShellPairs, &
                  Pivots => Chol2Vecs%Pivots, &
                  ShellParamsIdx => AOBasis%ShellParamsIdx, &
                  NAngFuncSpher => AOBasis%NAngFuncSpher, &
                  NAngFuncCart => AOBasis%NAngFuncCart, &
                  SpherAO => AOBasis%SpherAO &
                  )
                  call clock_start(timer)
                  allocate(Chol2Vecs%PivotShellPairs(NShellPairs))
                  allocate(Chol2Vecs%PivotShellPairLoc(NShellPairs))
                  allocate(Chol2Vecs%PivotShellPairDim(NShellPairs))
                  allocate(Chol2Vecs%PivotOrbPairs(NOrbPairs))
                  !
                  ! Process the array of pivot AO pairs from Step 1
                  ! into a computationally efficient format
                  !
                  if (SpherAO) then
                        call chol2_PivotShellPairs( &
                              Chol2Vecs%PivotShellPairs, &
                              Chol2Vecs%PivotShellPairLoc, &
                              Chol2Vecs%PivotShellPairDim, &
                              Chol2Vecs%PivotOrbPairs, &
                              Chol2Vecs%NPivotShellPairs, &
                              Pivots, &
                              ShellPairs, &
                              ShellPairDim, &
                              ShellPairLoc, &
                              ShellParamsIdx, &
                              NAngFuncSpher, &
                              NShellPairs, &
                              NCholesky)
                  else
                        call chol2_PivotShellPairs( &
                              Chol2Vecs%PivotShellPairs, &
                              Chol2Vecs%PivotShellPairLoc, &
                              Chol2Vecs%PivotShellPairDim, &
                              Chol2Vecs%PivotOrbPairs, &
                              Chol2Vecs%NPivotShellPairs, &
                              Pivots, &
                              ShellPairs, &
                              ShellPairDim, &
                              ShellPairLoc, &
                              ShellParamsIdx, &
                              NAngFuncCart, &
                              NShellPairs, &
                              NCholesky)
                  end if
                  !
                  ! Reduced-dimension Coulomb matrix including only the pivot AO orbital pairs pq, rs
                  !
                  call clock_start(DeltaT)
                  allocate(Vpqrs(NCholesky, NCholesky))
                  call chol2_Vpqrs( &
                        Vpqrs, &
                        Chol2Vecs%PivotShellPairs, &
                        Chol2Vecs%PivotShellPairLoc, &
                        Chol2Vecs%PivotShellPairDim, &
                        Chol2Vecs%PivotOrbPairs, &
                        Chol2Vecs%NPivotShellPairs, &
                        ShellPairs, &
                        AOBasis, &
                        Chol2Params%Kappa)
                  t_Vpqrs = clock_readwall(DeltaT)
                  !
                  ! V=L*L**T
                  !
                  call clock_start(DeltaT)
                  call real_Cholesky(Vpqrs)
                  t_Cholesky = clock_readwall(DeltaT)
                  !
                  ! L**(-1), used in Eq. 3 of Ref. 1.
                  !
                  call clock_start(DeltaT)
                  call real_LowerTriangularInverse(Vpqrs)
                  call move_alloc(from=Vpqrs, to=Chol2Vecs%Inv_L)
                  t_Inverse = clock_readwall(DeltaT)
                  t_Total = clock_readwall(timer)
                  !
                  ! Timings
                  !
                  call msg("Step 2 of Cholesky decomposition completed in " // str(t_Total,d=1) // " seconds")
                  call msg("Detailed timings in seconds")
                  call msg(lfield("Reduced-dimension matrix V", 30) // str(t_Vpqrs,d=1))
                  call msg(lfield("V=L*L**T", 30) // str(t_Cholesky,d=1))
                  call msg(lfield("L**(-1)", 30) // str(t_Inverse,d=1))
                  call blankline()
            end associate
      end subroutine chol2_Step2


      subroutine chol2_AllocWorkspace(Wabrs, Chol2Vecs)
            real(F64), dimension(:, :), allocatable, intent(out) :: Wabrs
            type(TChol2Vecs), intent(in)                         :: Chol2Vecs

            integer :: MaxSubsetDim
            
            MaxSubsetDim = maxval(Chol2Vecs%SubsetDim)
            allocate(Wabrs(MaxSubsetDim, Chol2Vecs%NVecs))
      end subroutine chol2_AllocWorkspace
      

      subroutine chol2_FullDimVectors_Batch(Rkpq, Wabrs, SubsetIdx, Chol2Vecs, AOBasis, Chol2Params)
            !
            ! Compute the full dimension Cholesky vectors for a given subset of AO indices pq
            ! according to Eq. 3 in Ref. 1.
            ! 
            ! 1. Sarai D. Folkestad, Eirik F. Kjønstad and Henrik Koch,
            !    J. Chem. Phys. 150, 194112 (2019);
            !    doi: 10.1063/1.5083802
            !
            real(F64), dimension(:, :), intent(out) :: Rkpq
            real(F64), dimension(:, :), intent(out) :: Wabrs
            integer, intent(in)                     :: SubsetIdx
            type(TChol2Vecs), intent(in)            :: Chol2Vecs
            type(TAOBasis), intent(in)              :: AOBasis
            type(TChol2Params), intent(in)          :: Chol2Params
            
            call chol2_Vabrs(Wabrs, &
                  Chol2Vecs%SubsetBounds(:, SubsetIdx), &
                  Chol2Vecs%PivotShellPairs, &
                  Chol2Vecs%PivotShellPairLoc, &
                  Chol2Vecs%PivotShellPairDim, &
                  Chol2Vecs%PivotOrbPairs, &
                  Chol2Vecs%NPivotShellPairs, &
                  Chol2Vecs%ShellPairLoc, &
                  Chol2Vecs%ShellPairs, &
                  Chol2Vecs%ShellPairDim, &
                  AOBasis, &
                  Chol2Params%Kappa)
            !
            ! Eq. 3 in Ref. 1
            ! R(k,pq) <- InvL(k,rs)*V(pq,rs)
            !
            associate ( &
                  InvL => Chol2Vecs%Inv_L, &
                  Vabrs => Wabrs &
                  )
                  call real_abT(Rkpq, InvL, Wabrs)
            end associate
      end subroutine chol2_FullDimVectors_Batch


      subroutine chol2_FullDimVectors(Rkpq, Chol2Vecs, AOBasis, Chol2Params)                    
            real(F64), dimension(:, :, :), allocatable, intent(inout) :: Rkpq[:]
            type(TChol2Vecs), intent(in)                              :: Chol2Vecs
            type(TAOBasis), intent(in)                                :: AOBasis
            type(TChol2Params), intent(in)                            :: Chol2Params
            

            real(F64), dimension(:, :), allocatable :: Wabrs
            integer :: SubsetIdx, X, Y
            integer :: MaxSubsetDim
            type(TClock) :: timer

            call clock_start(timer)
            call msg("Full-dimension Cholesky vectors of the Coulomb matrix will be stored in memory")            
            MaxSubsetDim = maxval(Chol2Vecs%SubsetDim)
            if (allocated(Rkpq)) deallocate(Rkpq)
            allocate(Rkpq(Chol2Vecs%NVecs, MaxSubsetDim, Chol2Vecs%NSubsets(1))[*])
            call chol2_AllocWorkspace(Wabrs, Chol2Vecs)
            Y = this_image()
            do X = 1, Chol2Vecs%NSubsets(1)
                  SubsetIdx = X + (Y - 1) * Chol2Vecs%NSubsets(1)
                  call chol2_FullDimVectors_Batch(Rkpq(:, :, X), Wabrs, SubsetIdx, &
                        Chol2Vecs, AOBasis, Chol2Params)
            end do
            call msg("Cholesky vectors computed in " // str(clock_readwall(timer),d=1) // " seconds")
            sync all
      end subroutine chol2_FullDimVectors
end module TwoStepCholesky
