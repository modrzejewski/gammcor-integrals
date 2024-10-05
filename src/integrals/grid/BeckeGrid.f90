module BeckeGrid
      use arithmetic
      use LebedevGrid
      use sys_definitions
      use grid_definitions
      use math_constants
      use periodic
      use string
      use display
      use GridFunctions
      
      implicit none
      
contains
      
      subroutine becke_MolecularGrid(Xk, Yk, Zk, Wk, NPoints, GridType, System, AOBasis, Phi2Thresh)
            real(F64), dimension(:), allocatable, intent(out)  :: Xk
            real(F64), dimension(:), allocatable, intent(out)  :: Yk
            real(F64), dimension(:), allocatable, intent(out)  :: Zk
            real(F64), dimension(:), allocatable, intent(out)  :: Wk
            integer, intent(out)                               :: NPoints            
            integer, intent(in)                                :: GridType
            type(TSystem), intent(in)                          :: System
            type(TAOBasis), intent(in)                         :: AOBasis
            real(F64), optional, intent(in)                    :: Phi2Thresh
            
            type(TLebedevGrid), dimension(LEBEDEV_NGRIDS) :: LebedevGrids
            integer, dimension(LEBEDEV_LMAX) :: LebedevLMap
            integer :: k, l, m, Z
            integer, dimension(:), allocatable :: ZList
            integer, dimension(:), allocatable :: ZCount
            integer, dimension(:), allocatable :: AtomElementMap
            integer :: NElements
            integer :: NAtoms
            integer :: NPointsInit, NPointsK, NPointsK_Init, NCentroidsK, NPointsBecke
            type(TBeckeGrid), dimension(:), allocatable :: IsolatedAtomGrids
            real(F64), dimension(:), allocatable :: AtomGridMap
            integer, dimension(:), allocatable :: Nk
            real(F64) :: WkBecke, WmBecke, WkBeckeNorm
            integer :: g0, g1, h0, h1
            real(F64), dimension(3) :: Rg
            integer :: ThisImage, NImages
            real(F64), dimension(:), allocatable :: Xl, Yl, Zl, Wl
            integer :: MaxNPoints
            integer :: NDiscarded1, NDiscarded2
            real(F64), parameter :: WkBeckeThresh = ZERO
            real(F64), parameter :: Phi2ThreshDefault = 1.0E-24_F64
            real(F64) :: PhiSquaredThresh

            ThisImage = this_image()
            NImages = num_images()
            NAtoms = System%NAtoms
            allocate(AtomElementMap(NAtoms))
            call sys_ElementsList(ZList, ZCount, AtomElementMap, NElements, System, SYS_ALL_ATOMS)
            call lebedev_init(LebedevGrids, LebedevLMap)
            if (present(Phi2Thresh)) then
                  PhiSquaredThresh = Phi2Thresh
            else
                  PhiSquaredThresh = Phi2ThreshDefault
            end if
            MaxNPoints = 0
            allocate(IsolatedAtomGrids(NElements))
            do k = 1, NElements
                  Z = ZList(k)
                  call becke_AtomicGrid( &
                        IsolatedAtomGrids(k)%X, &
                        IsolatedAtomGrids(k)%Y, &
                        IsolatedAtomGrids(k)%Z, &
                        IsolatedAtomGrids(k)%W, &
                        IsolatedAtomGrids(k)%NPoints, &
                        Z, GridType, LebedevGrids, LebedevLMap)
                  MaxNPoints = max(MaxNPoints, IsolatedAtomGrids(k)%NPoints)
            end do                        
            !
            ! Compute the number of grid points for each atom.
            ! This is the initial number of points which might be
            ! reduced after total weight computation and grid pruning.
            !            
            allocate(AtomGridMap(NAtoms))
            NPointsBecke = 0
            NPointsInit = 0
            do k = 1, NAtoms
                  NPointsK = IsolatedAtomGrids(AtomElementMap(k))%NPoints
                  AtomGridMap(k) = NPointsInit + 1
                  NPointsInit = NPointsInit + NPointsK
                  NPointsBecke = NPointsBecke + NPointsK
            end do
            allocate(Xk(NPointsInit))
            allocate(Yk(NPointsInit))
            allocate(Zk(NPointsInit))
            allocate(Wk(NPointsInit))
            allocate(Nk(NAtoms))
            Xk = ZERO
            Yk = ZERO
            Zk = ZERO
            Wk = ZERO
            Nk = 0
            associate ( &
                  AtomCoords => System%AtomCoords, &
                  SortedDistances => System%SortedDistances, &
                  SortedDistancesIdx => System%SortedDistancesIdx &
                  )
                  !$omp parallel &
                  !$omp private(k, l, m, NPointsK, NPointsK_Init, NCentroidsK, Rg, g0, g1, WkBecke, WmBecke, WkBeckeNorm) &
                  !$omp private(Xl, Yl, Zl, Wl) & 
                  !$omp default(shared)
                  allocate(Xl(MaxNPoints))
                  allocate(Yl(MaxNPoints))
                  allocate(Zl(MaxNPoints))
                  allocate(Wl(MaxNPoints))
                  !$omp do schedule(dynamic)
                  do k = ThisImage, NAtoms, NImages
                        NPointsK_Init = IsolatedAtomGrids(AtomElementMap(k))%NPoints
                        NPointsK = 0
                        do l = 1, NPointsK_Init
                              Rg(1) = AtomCoords(1, k) + IsolatedAtomGrids(AtomElementMap(k))%X(l)
                              Rg(2) = AtomCoords(2, k) + IsolatedAtomGrids(AtomElementMap(k))%Y(l)
                              Rg(3) = AtomCoords(3, k) + IsolatedAtomGrids(AtomElementMap(k))%Z(l)
                              WkBeckeNorm = ZERO
                              WkBecke = ZERO
                              do m = 1, NAtoms
                                    call becke_BeckeWeight(WmBecke, Rg, m, AtomCoords, &
                                          SortedDistances, SortedDistancesIdx, NAtoms)
                                    WkBeckeNorm = WkBeckeNorm + WmBecke
                                    if (m == k) WkBecke = WmBecke
                              end do
                              WkBecke = WkBecke / WkBeckeNorm
                              if (WkBecke > WkBeckeThresh) then
                                    NPointsK = NPointsK + 1
                                    Wl(NPointsK) = WkBecke * IsolatedAtomGrids(AtomElementMap(k))%W(l)
                                    Xl(NPointsK) = Rg(1)
                                    Yl(NPointsK) = Rg(2)
                                    Zl(NPointsK) = Rg(3)
                              end if
                        end do
                        g0 = AtomGridMap(k)
                        g1 = AtomGridMap(k) + NPointsK - 1
                        Nk(k) = NPointsK
                        Wk(g0:g1) = Wl(1:NPointsK)
                        Xk(g0:g1) = Xl(1:NPointsK)
                        Yk(g0:g1) = Yl(1:NPointsK)
                        Zk(g0:g1) = Zl(1:NPointsK)
                  end do
                  !$omp end do
                  deallocate(Xl, Yl, Zl, Wl)
                  !$omp end parallel
                  call co_sum(Wk, result_image=1)
                  call co_sum(Xk, result_image=1)
                  call co_sum(Yk, result_image=1)
                  call co_sum(Zk, result_image=1)
                  call co_sum(Nk, result_image=1)
                  if (ThisImage == 1) then
                        !
                        ! Remove gaps between atomic grids
                        !
                        allocate(Xl(MaxNPoints))
                        allocate(Yl(MaxNPoints))
                        allocate(Zl(MaxNPoints))
                        allocate(Wl(MaxNPoints))
                        NPoints = 0
                        do k = 1, NAtoms
                              NPointsK = Nk(k)
                              g0 = NPoints + 1
                              g1 = NPoints + NPointsK
                              if (g0 /= AtomGridMap(k)) then
                                    h0 = AtomGridMap(k) 
                                    h1 = AtomGridMap(k) + NPointsK - 1
                                    Xl(1:NPointsK) = Xk(h0:h1)
                                    Yl(1:NPointsK) = Yk(h0:h1) 
                                    Zl(1:NPointsK) = Zk(h0:h1)
                                    Wl(1:NPointsK) = Wk(h0:h1)
                                    Xk(g0:g1) = Xl(1:NPointsK)
                                    Yk(g0:g1) = Yl(1:NPointsK) 
                                    Zk(g0:g1) = Zl(1:NPointsK)
                                    Wk(g0:g1) = Wl(1:NPointsK)
                              end if
                              NPoints = NPoints + NPointsK
                        end do
                        if (NPoints < NPointsInit) then
                              Wk(NPoints+1:) = ZERO
                              Xk(NPoints+1:) = ZERO
                              Yk(NPoints+1:) = ZERO
                              Zk(NPoints+1:) = ZERO
                        end if
                        NDiscarded1 = NPointsInit - NPoints
                        !
                        ! Remove points for which all atomic orbitals are almost zero.
                        ! Removing those points is important for grids employed in tensor
                        ! hypercontraction because the THC subroutine assumes that at least
                        ! one single atomic orbital is nonnegligible at each grid point.
                        !
                        call becke_RemoveDistantPoints(Xk, Yk, Zk, Wk, NPoints, AOBasis, PhiSquaredThresh)
                        NDiscarded2 = NPointsInit - NDiscarded1 - NPoints
                  end if
                  call co_broadcast(Xk, source_image=1)
                  call co_broadcast(Yk, source_image=1)
                  call co_broadcast(Zk, source_image=1)
                  call co_broadcast(Wk, source_image=1)
                  call co_broadcast(NPoints, source_image=1)
            end associate
            if (ThisImage == 1) then
                  call msg("Initial Becke grid: " // trim(BECKE_PARAMS_LABEL(GridType)) // " " // str(NPointsBecke) // " points")
                  call msg("Discarded " // str(NDiscarded1) // " points with Becke weight <= " // str(WkBeckeThresh,d=1))
                  call msg("Discarded " // str(NDiscarded2) // " points with Max|PhiP(r)|**2 <= " // str(PhiSquaredThresh,d=1))
                  call msg("Final grid: " // str(NPoints))
                  call msg("Average " // str(NPoints/NAtoms) // " points per atom")
            end if
            sync all
      end subroutine becke_MolecularGrid


      subroutine becke_RemoveDistantPoints(X, Y, Z, W, NPoints, AOBasis, PhiSquaredThresh)
            real(F64), dimension(:), intent(inout) :: X
            real(F64), dimension(:), intent(inout) :: Y
            real(F64), dimension(:), intent(inout) :: Z
            real(F64), dimension(:), intent(inout) :: W
            integer, intent(inout)                 :: NPoints
            type(TAOBasis), intent(in)             :: AOBasis
            real(F64), intent(in)                  :: PhiSquaredThresh

            integer :: ShA, ShellParamsA, AtomA, La, NaCart, NaSpher
            integer :: MaxNAngFuncsCart, MaxNAngFuncsSpher
            real(F64), dimension(:, :), allocatable :: PhiCart, PhiSpher
            real(F64), dimension(:), allocatable :: MaxPhi2
            integer :: N, g

            associate (NShells => AOBasis%NShells, &
                  ShellCenters => AOBasis%ShellCenters, &
                  NormFactorsSpher => AOBasis%NormFactorsSpher, &
                  NormFactorsCart => AOBasis%NormFactorsCart, &
                  CartPolyX => AOBasis%CartPolyX, &
                  CartPolyY => AOBasis%CartPolyY, &
                  CartPolyZ => AOBasis%CartPolyZ, &
                  AtomCoords => AOBasis%AtomCoords, &
                  ShellParamsIdx => AOBasis%ShellParamsIdx, &
                  CntrCoeffs => AOBasis%CntrCoeffs, &
                  Exponents => AOBasis%Exponents, &
                  NPrimitives => AOBasis%NPrimitives, &
                  NAngFuncSpher => AOBasis%NAngFuncSpher, &
                  NAngFuncCart => AOBasis%NAngFuncCart, &
                  ShellMomentum => AOBasis%ShellMomentum, &
                  SpherAO => AOBasis%SpherAO, &
                  LmaxGTO => AOBasis%LmaxGTO &
                  )
                  MaxNAngFuncsCart = ((LmaxGTO + 1) * (LmaxGTO + 2)) / 2
                  MaxNAngFuncsSpher = 2 * LmaxGTO + 1
                  allocate(PhiCart(NPoints, MaxNAngFuncsCart))
                  if (SpherAO) then
                        allocate(PhiSpher(NPoints, MaxNAngFuncsSpher))
                  else
                        allocate(PhiSpher(1, 1))
                  end if
                  allocate(MaxPhi2(NPoints))
                  MaxPhi2 = ZERO
                  do ShA = 1, NShells
                        AtomA = ShellCenters(ShA)
                        ShellParamsA = ShellParamsIdx(ShA)
                        La = ShellMomentum(ShellParamsA)
                        NaCart = NAngFuncCart(ShellParamsA)
                        NaSpher = NAngFuncSpher(ShellParamsA)
                        if (SpherAO) then
                              call gridfunc_Orbitals_Shell(PhiSpher(:, 1:NaSpher), PhiCart(:, 1:NaCart), &
                                    X, Y, Z, NPoints, &
                                    NormFactorsSpher(:, ShellParamsA), &
                                    NormFactorsCart(:, ShellParamsA), &
                                    CntrCoeffs(:, ShellParamsA), &
                                    Exponents(:, ShellParamsA), &
                                    NPrimitives(ShellParamsA), &
                                    La, &
                                    AtomCoords(:, AtomA), &
                                    CartPolyX(:, La), &
                                    CartPolyY(:, La), &
                                    CartPolyZ(:, La), &
                                    NaSpher, &
                                    NaCart, &
                                    SpherAO)
                              call becke_MaxPhi2(MaxPhi2, PhiSpher, NPoints, NaSpher)
                        else
                              call gridfunc_Orbitals_Shell(PhiSpher, PhiCart(:, 1:NaCart), &
                                    X, Y, Z, NPoints, &
                                    NormFactorsSpher(:, ShellParamsA), &
                                    NormFactorsCart(:, ShellParamsA), &
                                    CntrCoeffs(:, ShellParamsA), &
                                    Exponents(:, ShellParamsA), &
                                    NPrimitives(ShellParamsA), &
                                    La, &
                                    AtomCoords(:, AtomA), &
                                    CartPolyX(:, La), &
                                    CartPolyY(:, La), &
                                    CartPolyZ(:, La), &
                                    NaSpher, &
                                    NaCart, &
                                    SpherAO)
                              call becke_MaxPhi2(MaxPhi2, PhiCart, NPoints, NaCart)
                        end if
                  end do
            end associate
            N = 0
            do g = 1, NPoints
                  if (MaxPhi2(g) > PhiSquaredThresh) then
                        N = N + 1
                        X(N) = X(g)
                        Y(N) = Y(g)
                        Z(N) = Z(g)
                        W(N) = W(g)
                  end if
            end do
            if (NPoints > N) then
                  X(N+1:NPoints) = ZERO
                  Y(N+1:NPoints) = ZERO
                  Z(N+1:NPoints) = ZERO
                  W(N+1:NPoints) = ZERO
            end if
            NPoints = N
      end subroutine becke_RemoveDistantPoints


      subroutine becke_MaxPhi2(MaxPhi2, Phi, NPoints, M)
            real(F64), dimension(NPoints), intent(inout)    :: MaxPhi2
            real(F64), dimension(NPoints, M), intent(inout) :: Phi
            integer, intent(in)                             :: NPoints
            integer, intent(in)                             :: M

            integer :: p
            
            Phi = Phi**2
            do p = 1, M
                  MaxPhi2 = max(MaxPhi2, Phi(:, p))
            end do
      end subroutine becke_MaxPhi2
      
      
      subroutine becke_AtomicGrid(Xk, Yk, Zk, Wk, NPoints, Z, GridType, LebedevGrids, LebedevLMap)
            real(F64), dimension(:), allocatable, intent(out) :: Xk
            real(F64), dimension(:), allocatable, intent(out) :: Yk
            real(F64), dimension(:), allocatable, intent(out) :: Zk
            real(F64), dimension(:), allocatable, intent(out) :: Wk
            integer, intent(out)                              :: NPoints
            integer, intent(in)                               :: Z
            integer, intent(in)                               :: GridType
            type(TLebedevGrid), dimension(:), intent(in)      :: LebedevGrids
            integer, dimension(:), intent(in)                 :: LebedevLMap

            integer :: k
            integer :: Lk, Nk
            integer :: NRadial
            integer :: N
            integer, dimension(:), allocatable :: LebedevIdx
            real(F64), dimension(:), allocatable :: Rk
            real(F64), dimension(:), allocatable :: WkRadial
            integer, dimension(:), allocatable :: s0, s1
            integer :: s0k, s1k
            
            NRadial = BECKE_PARAMS_NRadial(GridType)
            allocate(LebedevIdx(NRadial))
            allocate(Rk(NRadial))
            allocate(WkRadial(NRadial))
            allocate(s0(NRadial))
            allocate(s1(NRadial))
            call becke_EulerMaclaurin(Rk, WkRadial, Z, NRadial)
            select case (GridType)
            case (BECKE_PARAMS_THC)
                  do k = 1, NRadial
                        call becke_Pruning_Lk_Sherill2013(Lk, BECKE_PARAMS_THC_Lmax, Rk(k), Z)
                        LebedevIdx(k) = LebedevLMap(Lk)
                  end do
            case (BECKE_PARAMS_SG1, BECKE_PARAMS_MEDIUM, BECKE_PARAMS_FINE, BECKE_PARAMS_XFINE)
                  do k = 1, NRadial
                        call becke_Pruning_Lk_Pople1993(LebedevIdx(k), BECKE_PARAMS_LebedevIdx(:, GridType), Rk(k), Z)
                  end do
            case default
                  call msg("Invalid quadrature type", MSG_ERROR)
                  error stop
            end select
            N = 0
            do k = 1, NRadial
                  Nk = LebedevGrids(LebedevIdx(k))%NPoints
                  s0(k) = N + 1
                  s1(k) = N + Nk
                  N = N + Nk
            end do
            allocate(Xk(N))
            allocate(Yk(N))
            allocate(Zk(N))
            allocate(Wk(N))
            do k = 1, NRadial
                  s0k = s0(k)
                  s1k = s1(k)
                  Nk = s1k - s0k + 1
                  Xk(s0k:s1k) = Rk(k) * LebedevGrids(LebedevIdx(k))%x(1:Nk)
                  Yk(s0k:s1k) = Rk(k) * LebedevGrids(LebedevIdx(k))%y(1:Nk)
                  Zk(s0k:s1k) = Rk(k) * LebedevGrids(LebedevIdx(k))%z(1:Nk)
                  Wk(s0k:s1k) = WkRadial(k) * LebedevGrids(LebedevIdx(k))%w(1:Nk)                  
            end do
            NPoints = N
      end subroutine becke_AtomicGrid


      subroutine becke_EulerMaclaurin(Rk, Wk, Z, NRadial)
            !
            ! The Euler-Maclaurin radial grid for a single atom.
            !
            ! Gill, P., Johnson, B., Pople, J., A standard grid for density
            ! functional calculations, Chem. Phys. Lett. 209, 506(1993);
            ! doi: 10.1016/0009-2614(93)80125-9
            !
            real(F64), dimension(:), intent(out) :: Rk
            real(F64), dimension(:), intent(out) :: Wk
            integer, intent(in)                  :: Z
            integer, intent(in)                  :: NRadial

            real(F64) :: AtomicRadius
            integer :: k

            AtomicRadius = ATOMIC_RADII(Z)
            do k = 1, NRadial
                  !
                  ! Eq. 7 in Gill et al.
                  !
                  Rk(k) = AtomicRadius * k**2 / (NRadial + 1 - k)**2
                  !
                  ! The weights of Euler-Maclaurin quadrature are scaled by
                  ! 4Pi factor due to the unit normalization of the Lebedev grid
                  !
                  Wk(k) = (FOUR*Pi*Rk(k)**2)*((NRadial+1)*2*AtomicRadius*k)/(NRadial+1-k)**3
            end do
      end subroutine becke_EulerMaclaurin


      subroutine becke_BeckeWeight(w, Rg, i, AtomCoords, SortedDistances, SortedDistancesIdx, NAtoms)
            !
            ! Unnormalized Becke weight: equation 9 in Ref. 1.
            !
            ! 1. Stratmann, R.E., Scuseria, G.E., Frisch, M.J., Achieving linear scaling
            !    in exchange-correlation density functional quadratures,
            !    Chem. Phys. Lett. 257, 213 (1996); doi: 10.1016/0009-2614(96)00600-8
            !
            real(F64), intent(out)                 :: w
            real(F64), dimension(3), intent(in)    :: Rg
            integer, intent(in)                    :: i
            real(F64), dimension(:, :), intent(in) :: AtomCoords
            real(F64), dimension(:, :), intent(in) :: SortedDistances
            integer, dimension(:, :), intent(in)   :: SortedDistancesIdx
            integer, intent(in)                    :: NAtoms

            real(F64) :: mu, s
            integer :: k
            real(F64), dimension(3) :: Ri, Rj
            real(F64) :: Dig, Dij, Djg
            !
            ! Numerical constant used to compute the step function.
            ! Defined in Eq. 14 in Ref. 1
            !
            real(F64), parameter :: a = 0.64d+0

            Ri = AtomCoords(:, i)
            Dig = norm2(Ri-Rg)
            w = ONE
            !
            ! Loop over ith atom's neighbours, sorted
            ! according to increased distance from the ith atom
            !
            atomloop: do k = 2, NAtoms
                  Dij = SortedDistances(k, i)
                  Rj = AtomCoords(:, SortedDistancesIdx(k, i))
                  Djg = norm2(Rj-Rg)
                  !
                  ! Eq. 4 in Ref. 1
                  !
                  mu = (Dig - Djg) / (a * Dij)
                  if (mu <= -ONE) then
                        !
                        ! Condition for w=1.0
                        ! Eq. 15 in Ref. 1.
                        !
                        if (Dig <= (ONE-a)/TWO*Dij) then
                              exit atomloop
                        else
                              cycle atomloop
                        end if
                  else if (mu >= ONE) then
                        w = ZERO
                        exit atomloop
                  else
                        call becke_StepFunction(s, mu)
                        w = w * s
                  end if
            end do atomloop
      end subroutine becke_BeckeWeight


      pure subroutine becke_StepFunction(s, mu)
            ! -------------------------------------------------------------------------------
            ! Compute the step function (Eq. 8 in Ref. 1) for a given value of mu_ij.
            ! The polynomial employed to represent g(mu_ij) is defined in Eq. 14 of Ref. 1.
            ! It is assumed that the numerical constant a is already included in mu.
            ! This subroutine should only be use for for -1 <= mu <= 1.
            !
            ! The total step function of Stratmann et al. is defined as
            ! s(mu) = 1                     for mu <= -1
            ! s(mu) = becke_StepFunction(s, mu)  for -1 < mu < 1
            ! s(mu) = 0                     for 1 >= mu
            ! (See Fig. 1 in Ref. 1)
            !
            ! -------------------------------------------------------------------------------
            ! Mu 
            !      The input for the step function, mu_ij defined in Eq. 4 scaled by the
            !      numerical constant a.
            ! S
            !      Value of the step function (Eqs. 8 and 14)
            ! -------------------------------------------------------------------------------
            ! 1. Stratmann, R.E., Scuseria, G.E., Frisch, M.J., Achieving linear scaling
            !    in exchange-correlation density functional quadratures,
            !    Chem. Phys. Lett. 257, 213 (1996); doi: 10.1016/0009-2614(96)00600-8
            !
            real(F64), intent(out) :: s
            real(F64), intent(in) :: mu

            real(F64), parameter :: c1 = 35.0_F64 / 32.0_F64
            real(F64), parameter :: c3 = -35.0_F64 / 32.0_F64
            real(F64), parameter :: c5 = 21.0_F64 / 32.0_F64
            real(F64), parameter :: c7 = -5.0_F64 / 32.0_F64

            s = ONE/TWO - (c1 * mu + c3 * mu**3 + c5 * mu**5 + c7 * mu**7)
      end subroutine becke_StepFunction


      subroutine becke_Pruning_Lk_Sherill2013(Lk, Lmax, Rk, Z) 
            !
            ! Pruning function defined by Sherill et al. in the tensor
            ! hypercontraction paper.
            !
            ! 1. Robert M. Parrish, Edward G. Hohenstein, Todd J. Martínez,
            !    and C. David Sherrill1, Discrete variable representation in electronic
            !    structure theory: Quadrature grids for least-squares tensor
            !    hypercontraction, J. Chem. Phys. 138, 194107 (2013);
            !    doi: 10.1063/1.4802773
            !
            ! 2. Devin A. Matthews, Improved Grid Optimization and
            !    Fitting in Least Squares Tensor Hypercontraction
            !    J. Chem. Theory Comput. 16, 1382 (2020);
            !    doi: 10.1021/acs.jctc.9b01205
            !
            integer, intent(out)      :: Lk
            integer, intent(in)       :: Lmax
            real(F64), intent(in)     :: Rk
            integer, intent(in)       :: Z

            real(F64) :: AtomicRadius, R0

            AtomicRadius = ATOMIC_RADII(Z)
            R0 = Rk / AtomicRadius
            !
            ! Function defined at the bottom of the right column of page 7 in Ref. 1.
            ! Parameters taken from Matthews (supporting info of Ref. 2).
            !
            Lk = min(Lmax, ceiling(Lmax * (erf(1.2_F64 * R0) - 0.5_F64 * erf(1.2_F64 * R0 - 3.0_F64) - 0.5_F64)))
            Lk = min(LEBEDEV_LMAX, Lk)
            Lk = max(LEBEDEV_LMIN, Lk)
      end subroutine becke_Pruning_Lk_Sherill2013


      subroutine becke_Pruning_Lk_Pople1993(LebedevIdxK, LebedevIdx, Rk, Z)
            !
            ! 1. Gill, P., Johnson, B., Pople, J., A standard grid for density
            ! functional calculations, Chem. Phys. Lett. 209, 506(1993);
            ! doi: 10.1016/0009-2614(93)80125-9
            !
            integer, intent(out)              :: LebedevIdxK
            integer, dimension(5), intent(in) :: LebedevIdx
            real(F64), intent(in)             :: Rk
            integer, intent(in)               :: Z
            
            integer :: row
            real(F64) :: AtomicRadius, R0
            real(F64) :: r1, r2, r3, r4
            logical :: sg1avail
            !
            ! Pruning parameters defined in Table 4 of Ref. 1.
            !
            real(F64), dimension(0:2), parameter :: sg1_alpha1 = [ &
                  ! --------------------
                  ! VALUE      | ROW
                  ! --------------------
                  0.2500d+0, & ! 0 H-He
                  0.1667d+0, & ! 1 Li-Ne
                  0.1000d+0  & ! 2 Na-Ar
                  ]
            real(F64), dimension(0:2), parameter :: sg1_alpha2 = [ &
                  0.5000d+0, & ! 0 H-He
                  0.5000d+0, & ! 1 Li-Ne
                  0.4000d+0  & ! 2 Na-Ar
                  ]
            real(F64), dimension(0:2), parameter :: sg1_alpha3 = [ &
                  1.0000d+0, & ! 0 H-He
                  0.9000d+0, & ! 1 Li-Ne
                  0.8000d+0  & ! 2 Na-ar
                  ]
            real(F64), dimension(0:2), parameter :: sg1_alpha4 = [ &
                  4.5000d+0, & ! 0 H-He
                  3.5000d+0, & ! 1 Li-Ne
                  2.5000d+0  & ! 2 Na-Ar
                  ]

            AtomicRadius = ATOMIC_RADII(Z)
            R0 = Rk / AtomicRadius
            row = periodic_table_row(Z)
            if (row+1 > size(sg1_alpha1)) then
                  sg1avail = .false.
            else
                  sg1avail = .true.
            end if
            if (sg1avail) then
                  r1 = sg1_alpha1(row) * AtomicRadius
                  r2 = sg1_alpha2(row) * AtomicRadius
                  r3 = sg1_alpha3(row) * AtomicRadius
                  r4 = sg1_alpha4(row) * AtomicRadius
                  if (R0 <= r1) then
                        LebedevIdxK = LebedevIdx(1)
                  else if (R0 <= r2) then
                        LebedevIdxK = LebedevIdx(2)
                  else if (R0 <= r3) then
                        LebedevIdxK = LebedevIdx(3)
                  else if (R0 <= r4) then
                        LebedevIdxK = LebedevIdx(4)
                  else
                        LebedevIdxK = LebedevIdx(5)
                  end if
            else
                  LebedevIdxK = maxval(LebedevIdx)
            end if
      end subroutine becke_pruning_Lk_Pople1993
end module BeckeGrid
