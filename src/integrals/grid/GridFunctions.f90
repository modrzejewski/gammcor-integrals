module GridFunctions
      use arithmetic
      use basis_sets
      use Auto2e_SpherTransf
      use real_linalg

      implicit none

      integer, parameter :: GRID_ORBITAL_VALUE = 1
      integer, parameter :: GRID_ORBITAL_GRAD_X = 2
      integer, parameter :: GRID_ORBITAL_GRAD_Y = 3
      integer, parameter :: GRID_ORBITAL_GRAD_Z = 4

contains

      subroutine gridfunc_GGA_Variables(Rho, Phi, C, OccNumbers)
            real(F64), dimension(:, :), intent(out)   :: Rho
            real(F64), dimension(:, :, :), intent(in) :: Phi
            real(F64), dimension(:, :), intent(in)    :: C
            real(F64), dimension(:), intent(in)       :: OccNumbers

            integer :: NPoints, NMO
            integer :: i
            real(F64), dimension(:, :), allocatable :: W, Wx

            NPoints = size(Phi, dim=1)
            NMO = size(C, dim=2)
            allocate(W(NPoints, NMO))
            allocate(Wx(NPoints, NMO))
            call real_ab(W, Phi(:, :, 1), C)
            call multiply_W_W(Rho(:, 1), W, W, OccNumbers)
            do i = 2, 4
                  call real_ab(Wx, Phi(:, :, i), C)
                  call multiply_W_W(Rho(:, i), W, Wx, OccNumbers)
                  Rho(:, i) = TWO * Rho(:, i)
            end do
            
      contains
            
            subroutine multiply_W_W(Rho, W1, W2, OccNumbers)
                  real(F64), dimension(:), intent(out)   :: Rho
                  real(F64), dimension(:, :), intent(in) :: W1
                  real(F64), dimension(:, :), intent(in) :: W2
                  real(F64), dimension(:), intent(in)    :: OccNumbers

                  integer :: p
                  
                  Rho = ZERO
                  do p = 1, NMO
                        Rho(:) = Rho(:) + OccNumbers(p) * W1(:, p) * W2(:, p)
                  end do
            end subroutine multiply_W_W
      end subroutine gridfunc_GGA_Variables

      
      subroutine gridfunc_RhoIntegral(RhoIntegral, Xki, Wk, OccNumber)
            real(F64), intent(out)                 :: RhoIntegral
            real(F64), dimension(:, :), intent(in) :: Xki
            real(F64), dimension(:), intent(in)    :: Wk
            real(F64), intent(in)                  :: OccNumber

            integer :: NOcc, NPoints, i, k
            real(F64) :: t

            NPoints = size(Xki, dim=1)
            NOcc = size(Xki, dim=2)
            RhoIntegral = ZERO
            do i = 1, NOcc
                  t = ZERO
                  !$omp parallel do private(k) reduction(+:t)
                  do k = 1, NPoints
                        t = t + Wk(k) * Xki(k, i)**2
                  end do
                  !$omp end parallel do
                  RhoIntegral = RhoIntegral + t
            end do
            RhoIntegral = OccNumber * RhoIntegral
      end subroutine gridfunc_RhoIntegral


      subroutine gridfunc_Rho(Rho, Xg, Yg, Zg, EigenVecs, NVecs, AOBasis, MaxBatchDim)
            real(F64), dimension(:), intent(out)      :: Rho
            real(F64), dimension(:), intent(in)       :: Xg
            real(F64), dimension(:), intent(in)       :: Yg
            real(F64), dimension(:), intent(in)       :: Zg
            real(F64), dimension(:, :, :), intent(in) :: EigenVecs
            integer, dimension(2), intent(in)         :: NVecs
            type(TAOBasis), intent(in)                :: AOBasis
            integer, intent(in)                       :: MaxBatchDim

            integer :: NBatch, NGrid, g0, g1, NAO, N
            integer :: s, t, Nt, NSpins, MaxNVecs
            real(F64), dimension(:), allocatable :: Xgp, Xgi

            NSpins = size(EigenVecs, dim=3)
            NGrid = size(Rho, dim=1)
            NBatch = min(NGrid, MaxBatchDim)
            MaxNVecs = max(NVecs(1), NVecs(2))
            if (AOBasis%SpherAO) then
                  NAO = AOBasis%NAOSpher
            else
                  NAO = AOBasis%NAOCart
            end if
            allocate(Xgp(NBatch*NAO))
            allocate(Xgi(NBatch*MaxNVecs))
            Rho = ZERO
            N = NGrid / NBatch
            if (NGrid > NBatch .and. modulo(NGrid, NBatch) > 0) N = N + 1
            do t = 1, N
                  g0 = 1 + (t - 1) * NBatch
                  g1 = min(NGrid, t * NBatch)
                  Nt = g1 - g0 + 1
                  !
                  ! AO orbitals
                  !
                  call gridfunc_Orbitals(Xgp(1:Nt*NAO), Xg(g0:g1), Yg(g0:g1), Zg(g0:g1), Nt, NAO, AOBasis)
                  do s = 1, NSpins                        
                        call gridfunc_Orbitals_Transform(Xgi(1:Nt*NVecs(s)), Xgp, EigenVecs(:, :, s), &
                              Nt, NVecs(s), NAO)
                        call update_Rho(Rho(g0:g1), Xgi(1:Nt*NVecs(s)), Nt, NVecs(s))
                  end do
            end do
            deallocate(Xgp)
            deallocate(Xgi)
            
      contains
            subroutine update_Rho(Rho, Xgi, NGrid, NVecs)
                  real(F64), dimension(NGrid), intent(inout)     :: Rho
                  real(F64), dimension(NGrid, NVecs), intent(in) :: Xgi
                  integer, intent(in)                            :: NGrid
                  integer, intent(in)                            :: NVecs
                  
                  integer :: i
                  
                  do i = 1, NVecs
                        Rho(:) = Rho(:) + Xgi(:, i)**2
                  end do
            end subroutine update_Rho
      end subroutine gridfunc_Rho

      
      subroutine gridfunc_Orbitals_Transform(Xki, Xkp, Cpi, NGrid, NMO, NAO)
            real(F64), dimension(NGrid, NMO), intent(out) :: Xki
            real(F64), dimension(NGrid, NAO), intent(in)  :: Xkp
            real(F64), dimension(NAO, NMO), intent(in)    :: Cpi
            integer, intent(in)                           :: NGrid
            integer, intent(in)                           :: NMO
            integer, intent(in)                           :: NAO

            call real_ab(Xki, Xkp, Cpi)
      end subroutine gridfunc_Orbitals_Transform

      
      subroutine gridfunc_Orbitals(Phi, Xg, Yg, Zg, NPoints, NAO, AOBasis)
            real(F64), dimension(NPoints, NAO), intent(out) :: Phi
            real(F64), dimension(NPoints), intent(in)       :: Xg
            real(F64), dimension(NPoints), intent(in)       :: Yg
            real(F64), dimension(NPoints), intent(in)       :: Zg
            integer, intent(in)                             :: NPoints
            integer, intent(in)                             :: NAO
            type(TAOBasis), intent(in)                      :: AOBasis

            integer :: ShA, ShellParamsA, AtomA, La, NaCart, NaSpher
            integer :: p0, p1
            integer :: MaxNAngFuncsCart, MaxNAngFuncsSpher
            real(F64), dimension(:, :), allocatable :: PhiCart, PhiSpher

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
                  ShellLocCart => AOBasis%ShellLocCart, &
                  ShellLocSpher => AOBasis%ShellLocSpher, &
                  SpherAO => AOBasis%SpherAO, &
                  LmaxGTO => AOBasis%LmaxGTO &
                  )
                  MaxNAngFuncsCart = ((LmaxGTO + 1) * (LmaxGTO + 2)) / 2
                  MaxNAngFuncsSpher = 2 * LmaxGTO + 1
                  if (SpherAO) then
                        allocate(PhiCart(NPoints, MaxNAngFuncsCart))
                  else
                        allocate(PhiSpher(1, 1))
                  end if                  
                  do ShA = 1, NShells
                        AtomA = ShellCenters(ShA)
                        ShellParamsA = ShellParamsIdx(ShA)
                        La = ShellMomentum(ShellParamsA)
                        NaCart = NAngFuncCart(ShellParamsA)
                        NaSpher = NAngFuncSpher(ShellParamsA)
                        if (SpherAO) then
                              p0 = ShellLocSpher(ShA)
                              p1 = ShellLocSpher(ShA) + NaSpher - 1
                              call gridfunc_Orbitals_Shell(Phi(:, p0:p1), PhiCart(:, 1:NaCart), &
                                    Xg, Yg, Zg, NPoints, &
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
                        else
                              p0 = ShellLocCart(ShA)
                              p1 = ShellLocCart(ShA) + NaCart - 1
                              call gridfunc_Orbitals_Shell(PhiSpher, Phi(:, p0:p1), &
                                    Xg, Yg, Zg, NPoints, &
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
                        end if
                  end do
            end associate
      end subroutine gridfunc_Orbitals


      subroutine gridfunc_Orbitals_Grad(Phi, Xg, Yg, Zg, AOBasis)
            real(F64), dimension(:, :, :), intent(out) :: Phi
            real(F64), dimension(:), intent(in)        :: Xg
            real(F64), dimension(:), intent(in)        :: Yg
            real(F64), dimension(:), intent(in)        :: Zg
            type(TAOBasis), intent(in)                 :: AOBasis

            integer :: ShA, ShellParamsA, AtomA, La, NaCart, NaSpher
            integer :: p0
            integer :: MaxNAngFuncsCart, MaxNAngFuncsSpher
            real(F64), dimension(:, :, :), allocatable :: PhiCart, PhiSpher
            integer :: NPoints, NAO

            NPoints = size(Xg)
            NAO = size(Phi, dim=3)            
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
                  ShellLocCart => AOBasis%ShellLocCart, &
                  ShellLocSpher => AOBasis%ShellLocSpher, &
                  SpherAO => AOBasis%SpherAO, &
                  LmaxGTO => AOBasis%LmaxGTO &
                  )
                  MaxNAngFuncsCart = ((LmaxGTO + 1) * (LmaxGTO + 2)) / 2
                  MaxNAngFuncsSpher = 2 * LmaxGTO + 1
                  allocate(PhiCart(NPoints, MaxNAngFuncsCart, 4))
                  allocate(PhiSpher(NPoints, MaxNAngFuncsSpher, 4))
                  do ShA = 1, NShells
                        AtomA = ShellCenters(ShA)
                        ShellParamsA = ShellParamsIdx(ShA)
                        La = ShellMomentum(ShellParamsA)
                        NaCart = NAngFuncCart(ShellParamsA)
                        NaSpher = NAngFuncSpher(ShellParamsA)
                        if (SpherAO) then
                              p0 = ShellLocSpher(ShA)
                        else
                              p0 = ShellLocCart(ShA)
                        end if
                        call gridfunc_Orbitals_Grad_Shell(Phi, PhiSpher, PhiCart, &
                              Xg, Yg, Zg, NPoints, &
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
                              SpherAO, p0)
                  end do
            end associate
      end subroutine gridfunc_Orbitals_Grad
      
      
      subroutine gridfunc_Orbitals_Shell(PhiSpher, PhiCart, Xg, Yg, Zg, NPoints, NormFactorsSpher, &
            NormFactorsCart, CntrCoeffs, Exponents, NPrimitives, L, Ra, CartPolyX, CartPolyY, &
            CartPolyZ, NFuncSpher, NFuncCart, SpherAO)

            real(F64), dimension(:, :), intent(out)       :: PhiSpher
            real(F64), dimension(:, :), intent(out)       :: PhiCart
            real(F64), dimension(:), intent(in)           :: Xg
            real(F64), dimension(:), intent(in)           :: Yg
            real(F64), dimension(:), intent(in)           :: Zg
            integer, intent(in)                           :: NPoints
            real(F64), dimension(NFuncCart), intent(in)   :: NormFactorsCart
            real(F64), dimension(NFuncSpher), intent(in)  :: NormFactorsSpher
            real(F64), dimension(NPrimitives), intent(in) :: CntrCoeffs
            real(F64), dimension(NPrimitives), intent(in) :: Exponents
            integer, intent(in)                           :: NPrimitives
            integer, intent(in)                           :: L
            real(F64), dimension(3), intent(in)           :: Ra
            integer, dimension(NFuncCart), intent(in)     :: CartPolyX
            integer, dimension(NFuncCart), intent(in)     :: CartPolyY
            integer, dimension(NFuncCart), intent(in)     :: CartPolyZ
            integer, intent(in)                           :: NFuncSpher
            integer, intent(in)                           :: NFuncCart
            logical, intent(in)                           :: SpherAO
            
            real(F64) :: Xga, Yga, Zga, XgaLx, YgaLy, ZgaLz
            integer :: Lx, Ly, Lz
            integer :: k, i, v
            real(F64) :: R2, PhiRadial

            !$omp parallel do private(k, i, Xga, Yga, Zga, R2, PhiRadial) &
            !$omp shared(PhiCart) &
            !$omp default(shared)
            do k = 1, NPoints
                  Xga = Xg(k) - Ra(1)
                  Yga = Yg(k) - Ra(2)
                  Zga = Zg(k) - Ra(3)
                  R2 = Xga**2 + Yga**2 + Zga**2
                  PhiRadial = ZERO
                  do i = 1, NPrimitives
                        PhiRadial = PhiRadial + CntrCoeffs(i) * exp(-Exponents(i) * R2)
                  end do
                  PhiCart(k, 1) = PhiRadial
            end do
            !$omp end parallel do
            do v = 2, NFuncCart
                  PhiCart(:, v) = PhiCart(:, 1)
            end do
            !$omp parallel do private(v, k) collapse(2) &
            !$omp private(Lx, Ly, Lz, Xga, Yga, Zga, XgaLx, YgaLy, ZgaLz) &
            !$omp shared(PhiCart) &
            !$omp default(shared)
            do v = 1, NFuncCart
                  do k = 1, NPoints
                        Lx = CartPolyX(v)
                        Ly = CartPolyY(v)
                        Lz = CartPolyZ(v)
                        Xga = Xg(k) - Ra(1)
                        Yga = Yg(k) - Ra(2)
                        Zga = Zg(k) - Ra(3)
                        XgaLx = Xga**Lx
                        YgaLy = Yga**Ly
                        ZgaLz = Zga**Lz
                        PhiCart(k, v) = XgaLx * YgaLy * ZgaLz * PhiCart(k, v)
                  end do
            end do
            !$omp end parallel do
            if (SpherAO) then
                  call Auto1e_SpherTransf_U_Vector(L)%ptr(PhiSpher, PhiCart, NPoints)
                  do v = 1, NFuncSpher
                        PhiSpher(:, v) = NormFactorsSpher(v) * PhiSpher(:, v)
                  end do
            else
                  do v = 1, NFuncCart
                        PhiCart(:, v) = NormFactorsCart(v) * PhiCart(:, v)
                  end do
                  PhiSpher = ZERO
            end if
      end subroutine gridfunc_Orbitals_Shell


      subroutine gridfunc_Orbitals_Grad_Shell(Phi, PhiSpher, PhiCart, Xg, Yg, Zg, NPoints, NormFactorsSpher, &
            NormFactorsCart, CntrCoeffs, Exponents, NPrimitives, L, Ra, CartPolyX, CartPolyY, &
            CartPolyZ, NFuncSpher, NFuncCart, SpherAO, p0)

            real(F64), dimension(:, :, :), intent(inout)              :: Phi
            real(F64), dimension(NPoints, NFuncSpher, 4), intent(out) :: PhiSpher
            real(F64), dimension(NPoints, NFuncCart, 4), intent(out)  :: PhiCart
            real(F64), dimension(:), intent(in)                       :: Xg
            real(F64), dimension(:), intent(in)                       :: Yg
            real(F64), dimension(:), intent(in)                       :: Zg
            integer, intent(in)                                       :: NPoints
            real(F64), dimension(NFuncCart), intent(in)               :: NormFactorsCart
            real(F64), dimension(NFuncSpher), intent(in)              :: NormFactorsSpher
            real(F64), dimension(NPrimitives), intent(in)             :: CntrCoeffs
            real(F64), dimension(NPrimitives), intent(in)             :: Exponents
            integer, intent(in)                                       :: NPrimitives
            integer, intent(in)                                       :: L
            real(F64), dimension(3), intent(in)                       :: Ra
            integer, dimension(NFuncCart), intent(in)                 :: CartPolyX
            integer, dimension(NFuncCart), intent(in)                 :: CartPolyY
            integer, dimension(NFuncCart), intent(in)                 :: CartPolyZ
            integer, intent(in)                                       :: NFuncSpher
            integer, intent(in)                                       :: NFuncCart
            logical, intent(in)                                       :: SpherAO
            integer, intent(in)                                       :: p0
            
            real(F64) :: Xga, Yga, Zga, XgaLx, YgaLy, ZgaLz
            real(F64) :: DXgaLx, DYgaLy, DZgaLz, XYZ
            integer :: Lx, Ly, Lz
            integer :: k, i, v
            real(F64) :: CExpR2
            real(F64) :: R2, PhiRadial, PhiRadial_x, PhiRadial_y, PhiRadial_z

            !$omp parallel do private(k, i, Xga, Yga, Zga, R2, PhiRadial) &
            !$omp private(PhiRadial_x, PhiRadial_y, PhiRadial_z, CExpR2) &
            !$omp shared(PhiCart) &
            !$omp default(shared)
            do k = 1, NPoints
                  Xga = Xg(k) - Ra(1)
                  Yga = Yg(k) - Ra(2)
                  Zga = Zg(k) - Ra(3)
                  R2 = Xga**2 + Yga**2 + Zga**2
                  PhiRadial = ZERO
                  PhiRadial_x = ZERO
                  PhiRadial_y = ZERO
                  PhiRadial_z = ZERO
                  do i = 1, NPrimitives
                        CExpR2 = CntrCoeffs(i) * exp(-Exponents(i) * R2)
                        PhiRadial = PhiRadial + CExpR2
                        PhiRadial_x = PhiRadial_x - TWO * Exponents(i) * Xga * CExpR2
                        PhiRadial_y = PhiRadial_y - TWO * Exponents(i) * Yga * CExpR2
                        PhiRadial_z = PhiRadial_z - TWO * Exponents(i) * Zga * CExpR2
                  end do
                  PhiCart(k, 1, GRID_ORBITAL_VALUE) = PhiRadial
                  PhiCart(k, 1, GRID_ORBITAL_GRAD_X) = PhiRadial_x
                  PhiCart(k, 1, GRID_ORBITAL_GRAD_Y) = PhiRadial_y
                  PhiCart(k, 1, GRID_ORBITAL_GRAD_Z) = PhiRadial_z
            end do
            !$omp end parallel do
            do i = 1, 4
                  do v = 2, NFuncCart
                        PhiCart(:, v, i) = PhiCart(:, 1, i)
                  end do
            end do
            !$omp parallel do private(v, k) collapse(2) &
            !$omp private(Lx, Ly, Lz, Xga, Yga, Zga, XgaLx, YgaLy, ZgaLz) &
            !$omp private(DXgaLx, DYgaLy, DZgaLz, PhiRadial, XYZ) &
            !$omp private(PhiRadial_x, PhiRadial_y, PhiRadial_z) &
            !$omp shared(PhiCart) &
            !$omp default(shared)
            do v = 1, NFuncCart
                  do k = 1, NPoints
                        Lx = CartPolyX(v)
                        Ly = CartPolyY(v)
                        Lz = CartPolyZ(v)
                        Xga = Xg(k) - Ra(1)
                        Yga = Yg(k) - Ra(2)
                        Zga = Zg(k) - Ra(3)
                        XgaLx = Xga**Lx
                        YgaLy = Yga**Ly
                        ZgaLz = Zga**Lz
                        DXgaLx = Lx * Xga**(max(0, Lx-1))
                        DYgaLy = Ly * Yga**(max(0, Ly-1))
                        DZgaLz = Lz * Zga**(max(0, Lz-1))
                        PhiRadial = PhiCart(k, v, GRID_ORBITAL_VALUE)
                        PhiRadial_x = PhiCart(k, v, GRID_ORBITAL_GRAD_X)
                        PhiRadial_y = PhiCart(k, v, GRID_ORBITAL_GRAD_Y)
                        PhiRadial_z = PhiCart(k, v, GRID_ORBITAL_GRAD_Z)
                        XYZ = XgaLx * YgaLy * ZgaLz
                        PhiCart(k, v, GRID_ORBITAL_VALUE) = XYZ * PhiRadial
                        PhiCart(k, v, GRID_ORBITAL_GRAD_X) = XYZ * PhiRadial_x + DXgaLx*YgaLy*ZgaLz * PhiRadial
                        PhiCart(k, v, GRID_ORBITAL_GRAD_Y) = XYZ * PhiRadial_y + XgaLx*DYgaLy*ZgaLz * PhiRadial
                        PhiCart(k, v, GRID_ORBITAL_GRAD_Z) = XYZ * PhiRadial_z + XgaLx*YgaLy*DZgaLz * PhiRadial
                  end do
            end do
            !$omp end parallel do
            do i = 1, 4
                  if (SpherAO) then
                        call Auto1e_SpherTransf_U_Vector(L)%ptr(PhiSpher(:, :, i), PhiCart(:, :, i), NPoints)
                        do v = 1, NFuncSpher
                              Phi(:, p0+v-1, i) = NormFactorsSpher(v) * PhiSpher(:, v, i)
                        end do
                  else
                        do v = 1, NFuncCart
                              Phi(:, p0+v-1, i) = NormFactorsCart(v) * PhiCart(:, v, i)
                        end do
                  end if
            end do
      end subroutine gridfunc_Orbitals_Grad_Shell
end module GridFunctions
