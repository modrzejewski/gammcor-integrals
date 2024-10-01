module TensorHypercontraction
      use arithmetic
      use BeckeGrid
      use grid_definitions
      use GridFunctions
      use real_linalg
      use basis_sets
      use TwoStepCholesky
      use sys_definitions
      use thc_definitions
      use TwoStepCholesky_definitions
      
      implicit none

      type TCompressedVecs
            real(F64), dimension(:, :), allocatable :: L
      end type TCompressedVecs

contains

      subroutine thc_Test(Zgk, Xgp, AOBasis)            
            real(F64), dimension(:, :), intent(in) :: Zgk
            real(F64), dimension(:, :), intent(in) :: Xgp
            type(TAOBasis), intent(in)             :: AOBasis

            real(F64), parameter :: Kappa = ZERO
            integer :: PtrOffset, MaxNAngFunc
            real(F64), dimension(:, :), allocatable :: Zgh
            integer :: NGrid

            NGrid = size(Zgk, dim=1)
            allocate(Zgh(NGrid, NGrid))
            call real_abT(Zgh, Zgk, Zgk)
            !
            ! Select the set of two-electron integral subroutines:
            ! Cartesian or spherical AOs
            !
            if (AOBasis%SpherAO) then
                  PtrOffset = AUTO2E_SPHER_OFFSET
                  MaxNAngFunc = 2 * AOBasis%LmaxGTO + 1
            else
                  PtrOffset = 0
                  MaxNAngFunc = ((AOBasis%LmaxGTO + 1) * (AOBasis%LmaxGTO + 2)) / 2
            end if
            if (AOBasis%SpherAO) then
                  associate ( &
                        NAngFunc => AOBasis%NAngFuncSpher, &
                        NormFactors => AOBasis%NormFactorsSpher, &
                        ShellLoc => AOBasis%ShellLocSpher &
                        )                        
                        call thc_Test_pqpq(Zgh, Xgp, &
                              AOBasis%ShellCenters, &
                              AOBasis%AtomCoords, &
                              AOBasis%ShellParamsIdx, &
                              AOBasis%ShellMomentum, &
                              NAngFunc, &
                              AOBasis%NPrimitives, &
                              AOBasis%CntrCoeffs, &
                              AOBasis%Exponents, &
                              NormFactors, Kappa, &
                              MaxNAngFunc, &
                              AOBasis%NShells, &
                              ShellLoc, &
                              PtrOffset)
                  end associate
            else
                  associate ( &
                        NAngFunc => AOBasis%NAngFuncCart, &
                        NormFactors => AOBasis%NormFactorsCart, &
                        ShellLoc => AOBasis%ShellLocCart &
                        )
                        call thc_Test_pqpq(Zgh, Xgp, &
                              AOBasis%ShellCenters, &
                              AOBasis%AtomCoords, &
                              AOBasis%ShellParamsIdx, &
                              AOBasis%ShellMomentum, &
                              NAngFunc, &
                              AOBasis%NPrimitives, &
                              AOBasis%CntrCoeffs, &
                              AOBasis%Exponents, &
                              NormFactors, Kappa, &
                              MaxNAngFunc, &
                              AOBasis%NShells, &
                              ShellLoc, &
                              PtrOffset)
                  end associate
            end if
      end subroutine thc_Test


      subroutine thc_Test_pqpq(Zgh, Xgp, ShellCenters, AtomCoords, ShellParamsIdx, &
            ShellMomentum, NAngFunc, NPrimitives, CntrCoeffs, Exponents, NormFactors, Kappa, &
            MaxNAngFunc, NShells, ShellLoc, PtrOffset)
            
            real(F64), dimension(:, :), intent(in)             :: Zgh
            real(F64), dimension(:, :), intent(in)             :: Xgp
            integer, dimension(NShells), intent(in)            :: ShellCenters
            real(F64), dimension(:, :), intent(in)             :: AtomCoords
            integer, dimension(:), intent(in)                  :: ShellParamsIdx
            integer, dimension(:), intent(in)                  :: ShellMomentum
            integer, dimension(:), intent(in)                  :: NAngFunc
            integer, dimension(:), intent(in)                  :: NPrimitives
            real(F64), dimension(:, :), intent(in)             :: CntrCoeffs
            real(F64), dimension(:, :), intent(in)             :: Exponents
            real(F64), dimension(:, :), intent(in)             :: NormFactors
            real(F64), intent(in)                              :: Kappa
            integer, intent(in)                                :: MaxNAngFunc
            integer, intent(in)                                :: NShells
            integer, dimension(:), intent(in)                  :: ShellLoc
            integer, intent(in)                                :: PtrOffset

            integer :: n
            integer :: AtomA, AtomB, ShA, ShB, Lb, Nb, La, Na
            real(F64) :: Diff, RelDiff
            integer :: p0, p1, q0, q1
            integer :: ShellParamsA, ShellParamsB
            real(F64), dimension(MaxNAngFunc**4) :: Vexact
            real(F64), dimension(MaxNAngFUnc**4) :: Vapprox
            real(F64), dimension(3) :: MaxAbsError, MaxRelError

            MaxAbsError = ZERO
            MaxRelError = ZERO
            ShellB: do ShB = 1, NShells
                  AtomB = ShellCenters(ShB)
                  ShellParamsB = ShellParamsIdx(ShB)
                  Lb = ShellMomentum(ShellParamsB)
                  Nb = NAngFunc(ShellParamsB)
                  ShellA: do ShA = ShB, NShells
                        AtomA = ShellCenters(ShA)
                        ShellParamsA = ShellParamsIdx(ShA)
                        La = ShellMomentum(ShellParamsA)
                        Na = NAngFunc(ShellParamsA)
                        
                        call Auto2eERI(PtrOffset+auto2e_idx(Lb, La, Lb, La))%ptr( &
                              Vexact, &
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

                        p0 = ShellLoc(ShA)
                        p1 = ShellLoc(ShA) + Na - 1
                        q0 = ShellLoc(ShB)
                        q1 = ShellLoc(ShB) + Nb - 1
                        call thc_Test_Vapprox(Vapprox, Zgh, &
                              Xgp(:, p0:p1), Xgp(:, q0:q1), Xgp(:, p0:p1), Xgp(:, q0:q1), &
                              Na, Nb, Na, Nb)

                        do n = 1, Na*Nb*Na*Nb
                              Diff = abs(Vapprox(n)-Vexact(n))
                              if (Diff > MaxAbsError(1)) then
                                    MaxAbsError(1) = Diff
                                    MaxAbsError(2) = Vapprox(n)
                                    MaxAbsError(3) = Vexact(n)
                              end if
                              if (abs(Vexact(n)) > 1.0E-6_F64) then
                                    RelDiff = abs(Vapprox(n)-Vexact(n))/abs(Vexact(n))
                                    if (RelDiff > MaxRelError(1)) then
                                          MaxRelError(1) = RelDiff
                                          MaxRelError(2) = Vapprox(n)
                                          MaxRelError(3) = Vexact(n)
                                    end if
                              end if
                        end do
                        
                  end do ShellA
            end do ShellB

            call msg("Performed tests of diagonal Coulomb integrals (pq|pq)")
            call msg("Largest absolute error |Vexact-Vapprox|           " // str(MaxAbsError(1),d=1))
            call msg("Vapprox                                           " // str(MaxAbsError(2),d=8))
            call msg("Vexact                                            " // str(MaxAbsError(3),d=8))
            call msg("Largest relative error |Vexact-Vapprox|/|Vexact|  " // str(MaxRelError(1),d=1))
            call msg("Vapprox                                           " // str(MaxRelError(2),d=8))
            call msg("Vexact                                            " // str(MaxRelError(3),d=8))
      end subroutine thc_Test_pqpq


      subroutine thc_Test_Vapprox(V, Zgh, Xgp, Xgq, Xhr, Xhs, Na, Nb, Nc, Nd)
            real(F64), dimension(Na, Nb, Nc, Nd), intent(out) :: V
            real(F64), dimension(:, :), intent(in)            :: Zgh
            real(F64), dimension(:, :), intent(in)            :: Xgp
            real(F64), dimension(:, :), intent(in)            :: Xgq
            real(F64), dimension(:, :), intent(in)            :: Xhr
            real(F64), dimension(:, :), intent(in)            :: Xhs
            integer, intent(in)                               :: Na, Nb
            integer, intent(in)                               :: Nc, Nd

            integer :: p, q, r, s, g, h
            integer :: NGrid

            NGrid = size(Zgh, dim=1)
            V = ZERO
            do s = 1, Nd
                  do r = 1, Nc
                        do q = 1, Nb
                              do p = 1, Na
                                    do h = 1, NGrid
                                          do g = 1, NGrid
                                                V(p, q, r, s) = V(p, q, r, s) + Xgp(g, p) * Xgq(g, q) * &
                                                      Xhr(h, r) * Xhs(h, s) * Zgh(g, h)
                                          end do
                                    end do
                              end do
                        end do
                  end do
            end do
      end subroutine thc_Test_Vapprox


      subroutine thc_CoulombMatrix_QuadraticMemory(THCGrid, AOBasis, System, THCParams, Chol2Params)
            type(TCoulTHCGrid), intent(out) :: THCGrid
            type(TAOBasis), intent(in)      :: AOBasis
            type(TSystem), intent(in)       :: System
            type(TTHCParams), intent(in)    :: THCParams
            type(TChol2Params), intent(in)  :: Chol2Params

            type(TChol2Vecs) :: Chol2Vecs
            real(F64), dimension(1, 1, 1) :: Rkpq

            if (.not. THCParams%THC_QuadraticMemory) then
                  call msg("Invalid value of QuadraticMemory", MSG_ERROR)
                  error stop
            end if
            !
            ! Locate pivots of the Coulomb matrix. This step is required before
            ! generating the Z matrix.
            !
            call chol2_Algo_Koch_JCP2019(Chol2Vecs, AOBasis, Chol2Params)            
            call thc_Grid( &
                  THCGrid%Xgp, &
                  THCGrid%NGrid, &
                  THCGrid%NGridReduced, &
                  THCParams%THC_BeckeGridKind, &  ! parent molecular grid
                  THCParams%PhiSquaredThresh, &   ! threshold for small values of atomic orbitals
                  THCParams%QRThresh, &           ! Threshold for rank-revealing QR/Cholesky
                  THCParams%QRThreshReduced, &    ! threshold for the reduced-size grid
                  THCParams%THC_BlockDim, &       ! block dimension for the on the fly THC/Cholesky
                  AOBasis, System)
            call thc_Z( &
                  THCGrid%Zgk, &
                  THCGrid%ZgkReduced, &
                  THCGrid%NGrid, &
                  THCGrid%NGridReduced, &
                  THCGrid%Xgp, &
                  Rkpq, Chol2Vecs, Chol2Params, AOBasis, THCParams)
            allocate(THCGrid%Zgh(THCGrid%NGrid, THCGrid%NGrid))
            call real_abT(THCGrid%Zgh, THCGrid%Zgk, THCGrid%Zgk)
      end subroutine thc_CoulombMatrix_QuadraticMemory


      subroutine thc_ReduceGrid(THCGrid)
            type(TCoulTHCGrid), intent(inout) :: THCGrid
            
            real(F64), dimension(:, :), allocatable :: XgpFull
            integer :: NCholesky, NAO
            integer :: NGridReduced

            if (THCGrid%NGridReduced < THCGrid%NGrid) then
                  NGridReduced = THCGrid%NGridReduced
                  NCholesky = size(THCGrid%Zgk, dim=2)
                  NAO = size(THCGrid%Xgp, dim=2)
                  THCGrid%NGrid = NGridReduced
                  !
                  ! Remove the points outside of the reduced-size grid
                  !
                  deallocate(THCGrid%Zgh)
                  call move_alloc(From=THCGrid%ZgkReduced, to=THCGrid%Zgk)
                  call move_alloc(From=THCGrid%Xgp, to=XgpFull)
                  allocate(THCGrid%Xgp(NGridReduced, NAO))
                  THCGrid%Xgp(:, :) = XgpFull(1:NGridReduced, :)
                  deallocate(XgpFull)
                  allocate(THCGrid%Zgh(NGridReduced, NGridReduced))
                  call real_abT(THCGrid%Zgh, THCGrid%Zgk, THCGrid%Zgk)
            end if
      end subroutine thc_ReduceGrid
      

      subroutine thc_Grid(Xgp, NGrid, NGridReduced, BeckeGridKind, PhiSquaredThresh, &
            QRThresh, QRThreshReduced, BlockDim, AOBasis, System)
            !
            ! Cholesky/rank-revealing QR pruned molecular grid for the THC decomposition
            !
            ! (1) Generate initial Becke grid
            ! (2) Reduce the number of grid points by computing the pivots of
            !
            !     S(g,h) = Sum(pq) X(g,pq)X(h,pq)
            !
            ! The pivots of S correspond to the linearly independent columns
            ! of X(g,pq) in rank-revealing QR decomposition of X (Ref. 1).
            ! The matrix S for the parent grid is too big to be kept in memory
            ! so it's recomputed on the fly (M. Lesiuk, private communication).
            ! The Cholesky decomposition of S is equivalent to the pivot-only
            ! first step of the algorithm of Koch et al. (Ref. 2). The grid
            ! points below the pivot threshold are removed by dynamic
            ! reallocation of arrays between macro iterations of Cholesky.
            !
            ! This algorithm only generates pivots and full Cholesky vectors
            ! are never computed.
            !
            !
            ! 1. Robert M. Parrish, Edward G. Hohenstein, Todd J. Martinez,
            !    and C. David Sherrill, Tensor hypercontraction. II. Least-squares
            !    renormalization, J. Chem. Phys. 137, 224106 (2012);
            !     doi: 10.1063/1.4768233
            !
            ! 2. Devin A. Matthews, Improved Grid Optimization and
            !    Fitting in Least Squares Tensor Hypercontraction
            !    J. Chem. Theory Comput. 16, 1382 (2020);
            !    doi: 10.1021/acs.jctc.9b01205
            !
            ! 3. S.D. Folkestad, E.F. Kjonstad, and H. Koch.,
            !    An efficient algorithm for Cholesky decomposition of electron
            !    repulsion integrals, J. Chem. Phys. 150, 194112 (2019);
            !    doi: 10.1063/1.5083802
            !
            real(F64), dimension(:, :), allocatable, intent(out) :: Xgp
            integer, intent(out)                                 :: NGrid
            integer, intent(out)                                 :: NGridReduced
            integer, intent(in)                                  :: BeckeGridKind
            real(F64), intent(in)                                :: PhiSquaredThresh
            real(F64), intent(in)                                :: QRThresh
            real(F64), intent(in)                                :: QRThreshReduced
            integer, intent(in)                                  :: BlockDim
            type(TAOBasis), intent(in)                           :: AOBasis
            type(TSystem), intent(in)                            :: System

            real(F64), dimension(:), allocatable :: X, Y, Z, W
            integer :: NAO
            integer :: ThisImage
            real(F64) :: t_Grid
            type(TClock) :: timer_Grid
            
            ThisImage = this_image()
            call clock_start(timer_Grid)
            call blankline()
            call toprule()
            call msg("Tensor hypercontraction")
            call midrule()
            !
            ! Atomic orbital values at the nodes of the base grid. The orbitals are transformed
            ! to the spherical basis if AOBasis%SpherAO==.true.
            !
            if (AOBasis%SpherAO) then
                  NAO = AOBasis%NAOSpher
            else
                  NAO = AOBasis%NAOCart
            end if
            call becke_MolecularGrid(X, Y, Z, W, NGrid, BeckeGridKind, &
                  System, AOBasis, PhiSquaredThresh)
            deallocate(W)
            !
            ! Molecular grid pruning by rank-revealing Cholesky decomposition.
            !
            ! The Cholesky decomposition is used here instead of QR becuase
            ! it produces the same permuation matrix P but does not require
            ! the construction of full matrix X(pq,g). Instead of X, this
            ! variant employs the Cholesky decomposition of a smaller oveverlap
            ! matrix S(g,h)=Sum(pq)X(pq,g)*X(pq,h). S(g,h) is built on the fly,
            ! in blocks, without the need for memory storage.
            !
            call thc_Chol_Pivots(X, Y, Z, NGrid, NGridReduced, AOBasis, QRThresh, &
                  QRThreshReduced, BlockDim)
            allocate(Xgp(NGrid, NAO))
            call gridfunc_Orbitals(Xgp, X, Y, Z, NGrid, NAO, AOBasis)
            call thc_normalize_Xgp(Xgp)
            t_Grid = clock_readwall(timer_Grid)
            call msg("THC grid completed in " // str(t_Grid,d=1) // " seconds")
      end subroutine thc_Grid

      
      subroutine thc_Z(Zgk, ZgkReduced, NGrid, NGridReduced, Xgp, Rkpq, Chol2Vecs, &
            Chol2Params, AOBasis, THCParams)
            !
            ! Compute the Z factor of the THC decomposition
            ! of Cholesky-decomposed Coulomb integrals (Eq. 13 in Ref. 1)
            !
            ! (pq|rs) = Sum(k) R(k,pq)*R(k,rs)
            ! = Sum(gh) X(g,p)*X(g,q)*Z(g,h)*X(h,r)*X(h,s)
            ! = Sum(gh)Sum(k) X(g,p)*X(g,q)*Z'(g,k)*Z'(h,k)*X(h,r)*X(h,s)
            !
            ! where X is the vector of weighted orbital values on the molecular grid, also
            ! referred to as the collocation vector. The matrix Z'(g,k) is smaller
            ! than Z(g,h)
            !
            ! Z(g,h) = Sum(k=1...NCholesky) Z'(g,k)*Z'(h,k)
            ! g,h = 1...NGrid,  k = 1...NCholesky,
            ! in numerical tests NCholesky < NGrid
            !
            ! Z'(g,k) is obtained using the LDL**T decomposition.
            !
            ! 1. Robert M. Parrish, Edward G. Hohenstein, Todd J. Martinez,
            !    and C. David Sherrill, Tensor hypercontraction. II. Least-squares
            !    renormalization, J. Chem. Phys. 137, 224106 (2012);
            !    doi: 10.1063/1.4768233
            !
            ! 2. Devin A. Matthews, Improved Grid Optimization and
            !    Fitting in Least Squares Tensor Hypercontraction
            !    J. Chem. Theory Comput. 16, 1382 (2020);
            !    doi: 10.1021/acs.jctc.9b01205
            !
            real(F64), dimension(:, :), allocatable, intent(out) :: Zgk
            real(F64), dimension(:, :), allocatable, intent(out) :: ZgkReduced
            integer, intent(in)                                  :: NGrid
            integer, intent(in)                                  :: NGridReduced
            real(F64), dimension(:, :), intent(in)               :: Xgp
            real(F64), dimension(:, :, :), intent(in)            :: Rkpq
            type(TChol2Vecs), intent(in)                         :: Chol2Vecs
            type(TChol2Params), intent(in)                       :: Chol2Params
            type(TAOBasis), intent(in)                           :: AOBasis
            type(TTHCParams), intent(in)                         :: THCParams
            
            real(F64), dimension(:, :), allocatable :: Sgh
            integer :: ThisImage
            integer :: NAO, NCholesky
            real(F64) :: t_Z
            type(TClock) :: timer_Z

            ThisImage = this_image()
            call msg("Least squares fitting of Z: LDL**T linear system solver")
            if (THCParams%THC_QuadraticMemory) then
                  call msg("thc_Z will generate full-dimension Cholesky vecs on the fly")
            else
                  call msg("thc_Z will use precomputed full set of Cholesky vectors")
            end if
            call clock_start(timer_Z)
            NAO = size(Xgp, dim=2)
            NCholesky = Chol2Vecs%NVecs
            allocate(Zgk(NGrid, NCholesky))
            call thc_XR(Zgk, Rkpq, Xgp, NGrid, AOBasis, THCParams%THC_QuadraticMemory, &
                  Chol2Vecs, Chol2Params)
            if (NGridReduced < NGrid) then
                  allocate(ZgkReduced(NGridReduced, NCholesky))
                  ZgkReduced(:, :) = Zgk(1:NGridReduced, :)
            end if
            if (ThisImage == 1) then
                  allocate(Sgh(NGrid, NGrid))
                  call thc_S(Sgh, Xgp, NGrid)
                  call real_Axb_symmetric_sysv(Zgk, Sgh)
                  if (NGridReduced < NGrid) then
                        deallocate(Sgh)
                        allocate(Sgh(NGridReduced, NGridReduced))
                        call thc_S(Sgh, Xgp, NGridReduced)
                        call real_Axb_symmetric_sysv(ZgkReduced, Sgh)
                  end if
            end if
            call co_broadcast(Zgk, source_image=1)
            if (NGridReduced < NGrid) then
                  call co_broadcast(ZgkReduced, source_image=1)
            end if
            t_Z = clock_readwall(timer_Z)
            call msg("THC least squares completed in " // str(t_Z,d=1) // " seconds")
            call blankline()
            sync all
      end subroutine thc_Z


      subroutine thc_XR(XRgk, Rkpq, Xgp, NGrid, AOBasis, QuadraticMemory, Chol2Vecs, Chol2Params)
            !
            ! Compute metric matrix E(1:NGrid, 1:NGrid) = XR(1:NGrid,1:NCholesky)*XR(1:NGrid,1:NCholesky)
            ! (Eq. 35 in Ref. 1).
            !
            ! 1. Robert M. Parrish, Edward G. Hohenstein, Todd J. Martinez,
            !    and C. David Sherrill, Tensor hypercontraction. II. Least-squares
            !    renormalization, J. Chem. Phys. 137, 224106 (2012); doi: 10.1063/1.4768233
            !
            real(F64), dimension(:, :), intent(out)              :: XRgk
            real(F64), dimension(:, :, :), intent(in)            :: Rkpq
            real(F64), dimension(:, :), intent(in)               :: Xgp
            integer, intent(in)                                  :: NGrid
            type(TAOBasis), intent(in)                           :: AOBasis
            logical, intent(in)                                  :: QuadraticMemory
            type(TChol2Vecs), intent(in)                         :: Chol2Vecs
            type(TChol2Params), intent(in)                       :: Chol2Params

            integer :: X, Y, L
            integer :: Npq
            integer :: MaxSubsetDim, ldR
            integer :: ThisImage
            real(F64), dimension(:, :), allocatable :: Xgpq
            real(F64), dimension(:, :), allocatable :: Wabrs
            real(F64), dimension(:, :), allocatable :: RkpqBatch

            associate ( &
                  NCholesky => Chol2Vecs%NVecs, &
                  NSubsets => Chol2Vecs%NSubsets, &
                  ShellPairs => Chol2Vecs%ShellPairs, &
                  ShellPairLoc => Chol2Vecs%ShellPairLoc, &
                  ShellPairDim => Chol2Vecs%ShellPairDim, &
                  SubsetDim => Chol2Vecs%SubsetDim, &
                  SubsetBounds => Chol2Vecs%SubsetBounds &
                  )
                  ThisImage = this_image()
                  Y = ThisImage
                  MaxSubsetDim = maxval(SubsetDim)
                  allocate(Xgpq(NGrid, MaxSubsetDim))
                  if (QuadraticMemory) then
                        ldR = 0
                        call chol2_AllocWorkspace(Wabrs, Chol2Vecs)
                        allocate(RkpqBatch(NCholesky, MaxSubsetDim))
                  else
                        ldR = size(Rkpq, dim=1)
                  end if
                  XRgk = ZERO
                  do X = 1, NSubsets(1)
                        L = X + (Y - 1) * NSubsets(1)
                        Npq = SubsetDim(L)
                        !
                        ! Collocation matrices X(g,pq) for the current subset
                        ! of orbital pair indices pq. These need to be precomputed
                        ! to enable fast matrix multiplication.
                        !
                        if (AOBasis%SpherAO) then
                              associate ( &
                                    NAO => AOBasis%NAOSpher, &
                                    ShellLoc => AOBasis%ShellLocSpher, &
                                    NAngFunc => AOBasis%NAngFuncSpher, &
                                    ShellParamsIdx => AOBasis%ShellParamsIdx &
                                    )
                                    call thc_Xgpq(Xgpq, Xgp, SubsetBounds(:, L), NGrid, Npq, &
                                          ShellPairs, ShellPairLoc, ShellPairDim, ShellLoc, ShellParamsIdx, &
                                          NAngFunc, NAO)                                    
                              end associate
                        else
                              associate ( &
                                    NAO => AOBasis%NAOCart, &
                                    ShellLoc => AOBasis%ShellLocCart, &
                                    NAngFunc => AOBasis%NAngFuncCart, &
                                    ShellParamsIdx => AOBasis%ShellParamsIdx &
                                    )
                                    call thc_Xgpq(Xgpq, Xgp, SubsetBounds(:, L), NGrid, Npq, &
                                          ShellPairs, ShellPairLoc, ShellPairDim, ShellLoc, ShellParamsIdx, &
                                          NAngFunc, NAO)
                              end associate
                        end if
                        !
                        ! XR(1:NGrid,1:NCholesky) = 2 * X(1:NGrid,1:Npq)*R(1:NCholesky,1:Npq)**T
                        !
                        ! Note the scaling factor of 2. This takes into account the permutational symmetry
                        ! Sum(pq) X(g,pq)*R(k,pq) = 2 * Sum(p>q) X(g,pq)*R(k,pq) + diagonal terms
                        ! The diagonal terms of X(g,pq) are scaled by 1/2.
                        !
                        if (QuadraticMemory) then
                              call chol2_FullDimVectors_Batch(RkpqBatch, Wabrs, L, &
                                    Chol2Vecs, AOBasis, Chol2Params)
                              call real_abT_x(XRgk, NGrid, Xgpq, NGrid, RkpqBatch, NCholesky, &
                                    NGrid, NCholesky, Npq, TWO, ONE)      
                        else
                              call real_abT_x(XRgk, NGrid, Xgpq, NGrid, Rkpq(:, :, X), ldR, &
                                    NGrid, NCholesky, Npq, TWO, ONE)
                        end if
                  end do
                  call co_sum(XRgk, result_image=1)
                  sync all
            end associate
      end subroutine thc_XR


      subroutine thc_S(Sgh, Xgp, NGrid)
            real(F64), dimension(:, :), intent(out) :: Sgh
            real(F64), dimension(:, :), intent(in)  :: Xgp
            integer, intent(in)                     :: NGrid

            integer :: ldS, ldX, NAO
            integer :: g, h

            ldS = size(Sgh, dim=1)
            ldX = size(Xgp, dim=1)
            NAO = size(Xgp, dim=2)
            call real_abT_x(Sgh, ldS, Xgp, ldX, Xgp, ldX, NGrid, NGrid, NAO, ONE, ZERO)
            do h = 1, NGrid
                  do g = 1, NGrid
                        Sgh(g, h) = Sgh(g, h)**2
                  end do
            end do
      end subroutine thc_S


      subroutine thc_normalize_Xgp(Xgp)
            !
            ! Normalize collocation matrices as in Ref. 1. The sum of squares
            ! can't be zero because all zero grid points were screened out during
            ! Becke grid generation:
            !
            ! Rg: Max|PhiP(Rg)|**2>=PhiSquaredThresh
            !
            !
            ! 1. Robert M. Parrish, Edward G. Hohenstein, Todd J. Martinez,
            !    and C. David Sherrill, Tensor hypercontraction. II. Least-squares
            !    renormalization, J. Chem. Phys. 137, 224106 (2012);
            !    doi: 10.1063/1.4768233
            !
            real(F64), dimension(:, :), intent(inout) :: Xgp
            
            integer :: NAO, NGrid
            integer :: p, g
            real(F64), dimension(:), allocatable :: N
            
            NGrid = size(Xgp, dim=1)
            NAO = size(Xgp, dim=2)
            allocate(N(NGrid))
            !
            ! Numerically stable evaluation of L2 norm
            !
            N = norm2(Xgp, dim=2)
            !$omp parallel do private(g)
            do g = 1, NGrid
                  N(g) = ONE / N(g)
            end do
            !$omp end parallel do
            !$omp parallel do private(p)
            do p = 1, NAO
                  Xgp(:, p) = N(:) * Xgp(:, p)
            end do
            !$omp end parallel do
      end subroutine thc_normalize_Xgp


      subroutine thc_Xgpq(Xgpq, Xgp, SubsetBounds, NPoints, Npq, &
            ShellPairs, ShellPairLoc, ShellPairDim, ShellLoc, ShellParamsIdx, &
            NAngFunc, NAO)
            !
            ! Generate products of collocation vectors with compound
            ! orbital index pq, see Eq. 22 in Ref. 1
            !
            ! X(g,pq) = X(g,p) * X(x,q)         p>=q
            !
            ! where g is the grid point, p and q are orbital indices, and
            ! the compound index pq refers to nonredundant orbital pairs
            ! p>=q. The index pq is compatible with the pq of the Cholesky vectors.
            ! The arrray X(g,pq) can be matrix multiplied with Cholesky vectors
            ! using efficient linear algebra subroutines.
            !
            ! The array SubsetBounds defines the current shell pair subset
            !
            ! SubsetBounds(1) <= Included shell pairs <= SubsetBounds(2)
            !
            ! Xgpq is computed in subsets in order to avoid keeping in
            ! memory the whole three-index array.
            !
            ! 1. Robert M. Parrish, Edward G. Hohenstein, Todd J. Martinez,
            !    and C. David Sherrill, Tensor hypercontraction. II. Least-squares
            !    renormalization, J. Chem. Phys. 137, 224106 (2012); doi: 10.1063/1.4768233
            !
            real(F64), dimension(NPoints, Npq), intent(out)    :: Xgpq
            real(F64), dimension(NPoints, NAO), intent(in)     :: Xgp
            integer, dimension(2), intent(in)                  :: SubsetBounds
            integer, intent(in)                                :: NPoints
            integer, intent(in)                                :: Npq
            integer, dimension(:, :), intent(in)               :: ShellPairs
            integer, dimension(:, :), intent(in)               :: ShellPairLoc
            integer, dimension(:), intent(in)                  :: ShellPairDim
            integer, dimension(:), intent(in)                  :: ShellLoc
            integer, dimension(:), intent(in)                  :: ShellParamsIdx
            integer, dimension(:), intent(in)                  :: NAngFunc
            integer, intent(in)                                :: NAO

            integer :: ShAB, LocAB, Nab, pq0, pq1, p0, p1, q0, q1
            integer :: ShA, ShellParamsA, Na, LocA
            integer :: ShB, ShellParamsB, Nb, LocB

            !$omp parallel do schedule(guided) &
            !$omp default(shared) &
            !$omp private(ShA, Na, ShellParamsA, LocA) &
            !$omp private(ShB, Nb, ShellParamsB, LocB) &
            !$omp private(Nab, LocAB) &
            !$omp private(pq0, pq1, p0, p1, q0, q1) &
            !$omp private(ShAB)
            do ShAB = SubsetBounds(1), SubsetBounds(2)
                  LocAB = ShellPairLoc(CHOL2_SUBSET_STORAGE, ShAB)
                  Nab = ShellPairDim(ShAB)
                  pq0 = LocAB
                  pq1 = LocAB + Nab - 1
                  
                  ShA = ShellPairs(1, ShAB)
                  ShellParamsA = ShellParamsIdx(ShA)
                  Na = NAngFunc(ShellParamsA)
                  LocA = ShellLoc(ShA)
                  p0 = LocA
                  p1 = LocA + Na - 1
                  
                  ShB = ShellPairs(2, ShAB)
                  ShellParamsB = ShellParamsIdx(ShB)
                  Nb = NAngFunc(ShellParamsB)
                  LocB = ShellLoc(ShB)
                  q0 = LocB
                  q1 = LocB + Nb - 1

                  if (ShA /= ShB) then
                        call thc_Xgpq_ab(Xgpq(:, pq0:pq1), Xgp(:, p0:p1), Xgp(:, q0:q1), Na, Nb, NPoints)
                  else
                        call thc_Xgpq_aa(Xgpq(:, pq0:pq1), Xgp(:, p0:p1), Na, NPoints)
                  end if
            end do
            !$omp end parallel do
      end subroutine thc_Xgpq


      subroutine thc_Xgpq_ab(Xgpq, Xgp, Xgq, Np, Nq, NPoints)
            !
            ! ShA /= ShB
            ! ShC /= ShD
            !
            real(F64), dimension(NPoints, Np, Nq), intent(out) :: Xgpq
            real(F64), dimension(NPoints, Np), intent(in)      :: Xgp
            real(F64), dimension(NPoints, Nq), intent(in)      :: Xgq
            integer, intent(in)                                :: Np
            integer, intent(in)                                :: Nq
            integer, intent(in)                                :: NPoints
            
            integer :: p, q

            do q = 1, Nq
                  do p = 1, Np
                        Xgpq(:, p, q) = Xgp(:, p) * Xgq(:, q)
                  end do
            end do
      end subroutine thc_Xgpq_ab


      subroutine thc_Xgpq_aa(Xgpq, Xgp, Np, NPoints)
            !
            ! ShA == ShB
            !
            real(F64), dimension(NPoints, (Np*(Np+1))/2), intent(out) :: Xgpq
            real(F64), dimension(NPoints, Np), intent(in)             :: Xgp
            integer, intent(in)                                       :: Np
            integer, intent(in)                                       :: NPoints
            
            integer :: p, q, pq

            pq = 1
            do q = 1, Np
                  !
                  ! Diagonal
                  !
                  p = q
                  Xgpq(:, pq) = (ONE/TWO)*Xgp(:, p)**2
                  pq = pq + 1
                  !
                  ! Off-diagonal
                  !
                  do p = q + 1, Np
                        Xgpq(:, pq) = Xgp(:, p) * Xgp(:, q)                        
                        pq = pq + 1
                  end do
            end do
      end subroutine thc_Xgpq_aa


      subroutine thc_Chol_Pivots(X, Y, Z, NGrid, NGridReduced, AOBasis, QRThresh, &
            QRThreshReduced, BlockDim)
            !
            ! Perform grid pruning by the Cholesky decomposition of
            !
            ! S(g,h) = Sum(pq) X(g,pq)X(h,pq)
            !
            ! The pivots of S correspond to the linearly independent columns
            ! of X(g,pq) in rank-revealing QR decomposition of X (Ref. 1).
            ! The matrix S for the parent grid is too big to be kept in memory
            ! so it's recomputed on the fly (M. Lesiuk, private communication).
            ! The Cholesky decomposition of S is equivalent to the pivot-only
            ! first step of the algorithm of Koch et al. (Ref. 2). The grid
            ! points below the pivot threshold are removed by dynamic
            ! reallocation of arrays between macro iterations of Cholesky.
            !
            ! This algorithm only generates pivots and full Cholesky vectors
            ! are never computed.
            !
            ! 1. Robert M. Parrish, Edward G. Hohenstein, Todd J. Martinez,
            !    and C. David Sherrill, Tensor hypercontraction. II. Least-squares
            !    renormalization, J. Chem. Phys. 137, 224106 (2012); doi: 10.1063/1.4768233
            !
            ! 2. S.D. Folkestad, E.F. Kjonstad, and H. Koch.,
            !    An efficient algorithm for Cholesky decomposition of electron
            !    repulsion integrals, J. Chem. Phys. 150, 194112 (2019);
            !    doi: 10.1063/1.5083802
            !
            ! 3. Devin A. Matthews, Improved Grid Optimization and
            !    Fitting in Least Squares Tensor Hypercontraction
            !    J. Chem. Theory Comput. 16, 1382 (2020);
            !    doi: 10.1021/acs.jctc.9b01205
            !
            real(F64), dimension(:), allocatable, intent(inout) :: X
            real(F64), dimension(:), allocatable, intent(inout) :: Y
            real(F64), dimension(:), allocatable, intent(inout) :: Z
            integer, intent(inout)                              :: NGrid
            integer, intent(out)                                :: NGridReduced            
            type(TAOBasis), intent(in)                          :: AOBasis
            real(F64), intent(in)                               :: QRThresh
            real(F64), intent(in)                               :: QRThreshReduced
            integer, intent(in)                                 :: BlockDim

            integer :: NPivots, NCandidates, NAO, NPivotsReduced
            real(F64) :: PivotThresh, PivotThreshReduced
            real(F64), dimension(:), allocatable :: D
            real(F64), dimension(:, :), allocatable :: Xgp
            real(F64), dimension(:), allocatable :: XhJ, YhJ, ZhJ
            integer :: b, MaxIters, J
            type(TCompressedVecs), dimension(:), allocatable :: V
            logical :: Converged
            type(TClock) :: timer_Chol, timer_Iter, timer_NextBlock, timer_Compress
            real(F64) :: t_Chol, t_Iter, t_Xgp, t_NextBlock, t_Compress
            real(F64) :: Dmax
            integer :: MaxJ
            real(F64) :: MemoryMB
            integer :: Percentage
            character(:), allocatable :: line
            integer :: ThisImage

            ThisImage = this_image()
            if (ThisImage == 1) then
                  call clock_start(timer_Chol)
                  if (AOBasis%SpherAO) then
                        NAO = AOBasis%NAOSpher
                  else
                        NAO = AOBasis%NAOCart
                  end if
                  call blankline()
                  call msg("Pivoted Cholesky decomposition with on-the-fly generation of matrix S", underline=.true.)
                  call msg("QR relative pivot threshold Eps=" // str(QRThresh,d=1))
                  if (QRThreshReduced > ZERO) call msg("QR relative pivot threshold (reduced grid) Eps=" // str(QRThreshReduced,d=1))
                  call msg("Block dimension: " // str(BlockDim) // " Cholesky vectors per macro iteration")
                  call msg("Starting from Becke grid composed of " // str(NGrid) // " points")
                  !
                  ! Due to normalization of Xgp, the initial diagonals are exactly equal to 1.0.
                  !
                  allocate(D(NGrid))
                  D = ONE
                  NCandidates = NGrid
                  !
                  ! QRThresh is the rank-revealing QR decomposition threshold
                  ! Eps defined in Section 1 of Ref. 2.
                  ! QRThresh is translated into PivotThresh of the Cholesky
                  ! decomposition. See the definition of
                  ! "function prune_grid" in the supplementary info of Ref. 3.
                  !
                  PivotThresh = QRThresh**2
                  if (QRThreshReduced > ZERO) then
                        PivotThreshReduced = QRThreshReduced**2
                  else
                        PivotThreshReduced = -ONE
                  end if
                  call msg("Cholesky absolute pivot threshold Eps**2=" // str(PivotThresh,d=1))
                  if (QRThreshReduced > ZERO) call msg("Cholesky absolute pivot threshold (reduced grid) Eps**2=" &
                        // str(PivotThreshReduced,d=1))
                  MaxIters = NCandidates / BlockDim
                  if (modulo(NCandidates, BlockDim) > 0) MaxIters = MaxIters + 1
                  allocate(XhJ(NCandidates))
                  allocate(YhJ(NCandidates))
                  allocate(Zhj(NCandidates))
                  allocate(V(MaxIters))
                  t_Xgp = ZERO
                  t_NextBlock = ZERO
                  t_Compress = ZERO
                  Converged = .false.
                  NPivots = 0
                  NPivotsReduced = 0
                  call blankline()
                  call toprule()
                  line = lfield("#", 5) // lfield("Dmax", 12) // lfield("Memory (MB)",15) //  &
                        lfield("NPivots", 12) // lfield("%", 5) // lfield("Time", 12)
                  call msg(line)
                  call midrule()
                  do b = 1, MaxIters
                        call clock_start(timer_Iter)
                        allocate(Xgp(NCandidates, NAO))
                        call gridfunc_Orbitals(Xgp, X, Y, Z, NCandidates, NAO, AOBasis)
                        call thc_normalize_Xgp(Xgp)
                        t_Xgp = t_Xgp + clock_readwall(timer_Iter)
                        !
                        ! Generate a batch of J <= BlockDim Cholesky vectors.
                        ! If on exit J < BlockDim then all remaining residuals are below
                        ! the pivot threshold and convergence is achieved.
                        !
                        MaxJ = min(BlockDim, size(XhJ)-NPivots)
                        call clock_start(timer_NextBlock)
                        call thc_Chol_NextBlock(V, XhJ, YhJ, ZhJ, J, NPivotsReduced, D, X, Y, Z, NPivots, Xgp, &
                              NCandidates, BlockDim, PivotThresh, PivotThreshReduced, NAO, MaxJ)
                        t_NextBlock = t_NextBlock + clock_readwall(timer_NextBlock)
                        NPivots = NPivots + J
                        MemoryMB = (storage_size(Xgp,kind=I64)*(NCandidates*NAO+NCandidates*NPivots)) &
                              / (8.0_F64*1024.0_F64*1024.0_F64)
                        !
                        ! Remove points corresponding to diagonal elements
                        ! below the pivot threshold. Reallocate arrays to new size.
                        ! This is a crucial part of the two-step Cholesky algorithm
                        ! of Koch et al. which lowers memory and computational cost.
                        !
                        call clock_start(timer_Compress)
                        call thc_Chol_Compress(Dmax, X, Y, Z, D, V, NCandidates, NPivots, &
                              BlockDim, PivotThresh)
                        t_Compress = t_Compress + clock_readwall(timer_Compress)
                        deallocate(Xgp)
                        t_Iter = clock_readwall(timer_Iter)
                        Percentage = nint(100*(real(NGrid-NCandidates,F64)/real(NGrid,F64)))
                        line = lfield(str(b), 5) // lfield(str(Dmax,d=1), 12) // lfield(str(MemoryMB,d=1),15) // &
                              lfield(str(NPivots), 12) // lfield(str(Percentage), 5) // lfield(str(t_Iter,d=1), 12)
                        call msg(line)
                        if (NCandidates == 0) then
                              Converged = .true.
                              exit
                        end if
                  end do
                  call midrule()
                  if (Converged) then
                        call msg("Rank-revealing decomposition converged with " // str(NPivots) // " pivots")
                        call msg("Average " // str(NPivots/AOBasis%NAtoms) // " points per atom")
                        call msg("NGridTHC/NAO=" // str(real(NPivots,F64)/NAO, d=1))
                        if (QRThreshReduced > ZERO) then
                              call msg("Reduced grid: " // str(NPivotsReduced) // " pivots")
                              call msg("Reduced grid: " // "NGridTHC/NAO=" // str(real(NPivotsReduced,F64)/NAO, d=1))
                        end if
                  else
                        call msg("Cholesky decomposistion failed to converge", MSG_ERROR)
                        error stop
                  end if
                  NGrid = NPivots
                  if (QRThreshReduced > ZERO) then
                        NGridReduced = NPivotsReduced
                  else
                        NGridReduced = NGrid
                  end if
                  t_Chol = clock_readwall(timer_Chol)
                  call msg("Total time " // str(t_Chol,d=1) // " seconds")
                  call blankline()
                  call msg("Detailed timings", underline=.true.)
                  call msg(lfield("Orbital values Xgp", 35) // str(t_Xgp,d=1))
                  call msg(lfield("Overlap matrix S + Cholesky vecs", 35) // str(t_NextBlock,d=1))
                  call msg(lfield("Reallocation", 35) // str(t_Compress,d=1))
                  call blankline()
                  deallocate(X, Y, Z)
                  allocate(X(NGrid), Y(NGrid), Z(NGrid))
                  X(:) = Xhj(1:NGrid)
                  Y(:) = Yhj(1:NGrid)
                  Z(:) = Zhj(1:NGrid)
            end if
            call co_broadcast(NGrid, source_image=1)
            if (ThisImage /= 1) then
                  deallocate(X, Y, Z)
                  allocate(X(NGrid), Y(NGrid), Z(NGrid))
            end if
            call co_broadcast(X, source_image=1)
            call co_broadcast(Y, source_image=1)
            call co_broadcast(Z, source_image=1)
            sync all
      end subroutine thc_Chol_Pivots


      subroutine thc_Chol_NextBlock(V, XhJ, YhJ, ZhJ, J, NPivotsReduced, D, X, Y, Z, NPivots, Xgp, &
            NCandidates, BlockDim, PivotThresh, PivotThreshReduced, NAO, MaxJ)

            type(TCompressedVecs), dimension(:), intent(inout)          :: V
            real(F64), dimension(:), intent(inout)                      :: XhJ
            real(F64), dimension(:), intent(inout)                      :: YhJ
            real(F64), dimension(:), intent(inout)                      :: ZhJ
            integer, intent(out)                                        :: J
            integer, intent(inout)                                      :: NPivotsReduced
            real(F64), dimension(NCandidates), intent(inout)            :: D
            real(F64), dimension(NCandidates), intent(in)               :: X
            real(F64), dimension(NCandidates), intent(in)               :: Y
            real(F64), dimension(NCandidates), intent(in)               :: Z
            integer, intent(in)                                         :: NPivots
            real(F64), dimension(NCandidates, NAO), intent(in)          :: Xgp
            integer, intent(in)                                         :: NCandidates
            integer, intent(in)                                         :: BlockDim
            real(F64), intent(in)                                       :: PivotThresh
            real(F64), intent(in)                                       :: PivotThreshReduced
            integer, intent(in)                                         :: NAO
            integer, intent(in)                                         :: MaxJ

            integer :: k, hJ
            real(F64) :: DhJ
            integer :: b, NBlocks, CurrentBlock
            real(F64), dimension(:), allocatable :: XphJ, LJ
            !
            ! Allocate storage space for the next block of Cholesky vectors
            !
            NBlocks = NPivots / BlockDim
            CurrentBlock = NBlocks + 1
            if (allocated(V(CurrentBlock)%L)) then
                  call msg("Error on entry to thc_Chol_NextBlock: invalid number of blocks", MSG_ERROR)
                  error stop
            end if
            allocate(V(CurrentBlock)%L(BlockDim, NCandidates))
            allocate(XphJ(NAO))
            allocate(LJ(NCandidates))
            !
            ! J is the counter of Cholesky vectors
            ! in the new block
            !
            J = 0
            do k = 1, MaxJ
                  hJ = maxloc(D, dim=1)
                  DhJ = D(hJ)
                  if (DhJ > PivotThresh) then
                        J = J + 1
                        if (DhJ > PivotThreshReduced) then
                              NPivotsReduced = NPivots + J
                        end if
                        !
                        ! XYZ coordinates of the current pivot
                        !
                        XhJ(NPivots+J) = X(hJ)
                        YhJ(NPivots+J) = Y(hJ)
                        ZhJ(NPivots+J) = Z(hJ)
                        !
                        ! LJ(g) <- S(g,hJ), g is arbitrary
                        !
                        XphJ = Xgp(hJ, :)
                        call real_av_x(LJ, Xgp, NCandidates, XphJ, NCandidates, NAO, ONE, ZERO)
                        LJ(:) = LJ(:)**2
                        !
                        ! Subtract the contribution of the Cholesky vectors from
                        ! the previous blocks
                        !
                        ! LJ(g) <- LJ(g) - Sum(k=1,NPivots) L(k,g)*L(k,hJ)
                        !
                        do b = 1, NBlocks
                              call real_aTv_x(LJ, V(b)%L, BlockDim, V(b)%L(:, hJ), &
                                    BlockDim, NCandidates, -ONE, ONE)
                        end do
                        !
                        ! Subract the contribution of the Cholesky Vectors
                        ! from the current block
                        !
                        ! LJ(g) <- LJ(g) - Sum(k=1,J-1) L(NPivots+k,g)*L(NPivots+k,hJ)
                        !
                        if (J > 1) then
                              call real_aTv_x(LJ, V(CurrentBlock)%L, BlockDim, V(CurrentBlock)%L(:, hJ), &
                                    J-1, NCandidates, -ONE, ONE)
                        end if
                        LJ(:) = LJ(:) / Sqrt(DhJ)
                        V(CurrentBlock)%L(J, :) = LJ
                        !
                        ! D(g) <- D(g) - L(g,J)*L(g,J)
                        !
                        D(:) = D(:) - LJ(:)**2
                  else
                        exit
                  end if
            end do
      end subroutine thc_Chol_NextBlock


      subroutine thc_Chol_Compress(Dmax, X, Y, Z, D, V, NCandidates, NPivots, BlockDim, PivotThresh)
            real(F64), intent(out)                               :: Dmax
            real(F64), dimension(:), intent(inout)               :: X
            real(F64), dimension(:), intent(inout)               :: Y
            real(F64), dimension(:), intent(inout)               :: Z
            real(F64), dimension(:), intent(inout)               :: D
            type(TCompressedVecs), dimension(:), intent(inout)   :: V
            integer, intent(inout)                               :: NCandidates
            integer, intent(in)                                  :: NPivots
            integer, intent(in)                                  :: BlockDim
            real(F64), intent(in)                                :: PivotThresh

            integer :: g, b
            integer :: NCandidatesNew, NBlocks
            real(F64), dimension(:, :), allocatable :: T
            integer, dimension(:), allocatable :: Map

            allocate(Map(NCandidates))
            NCandidatesNew = 0
            Dmax = ZERO
            do g = 1, NCandidates
                  Dmax = max(Dmax, D(g))
                  if (D(g) > PivotThresh) then
                        NCandidatesNew = NCandidatesNew + 1
                        Map(NCandidatesNew) = g
                  end if
            end do
            NCandidates = NCandidatesNew
            if (NCandidates > 0) then
                  !
                  ! Move pivot candidates to the front
                  !
                  do g = 1, NCandidates
                        D(g) = D(Map(g))
                        X(g) = X(Map(g))
                        Y(g) = Y(Map(g))
                        Z(g) = Z(Map(g))
                  end do
                  !
                  ! Move pivot candidate columns to the front
                  ! and remove non-candidate columns from
                  ! the storage.
                  !
                  NBlocks = NPivots / BlockDim
                  do b = 1, NBlocks
                        call move_alloc(from=V(b)%L, to=T)
                        !
                        ! Reallocation of L with smaller number of columns saves memory
                        !
                        allocate(V(b)%L(BlockDim, NCandidates))
                        associate (LNew => V(b)%L)
                              !$omp parallel do private(g)
                              do g = 1, NCandidates
                                    LNew(:, g) = T(:, Map(g))
                              end do
                              !$omp end parallel do
                        end associate
                  end do
                  if (NBlocks >= 1) deallocate(T)                  
            end if
      end subroutine thc_Chol_Compress
end module TensorHypercontraction
