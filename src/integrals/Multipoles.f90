module Multipoles
      use arithmetic
      use math_constants
      use real_linalg
      use basis_sets
      use sys_definitions
      use OneElectronInts
      use string
      use display

      implicit none
      !
      ! Conventions of the quadrupole operator
      ! ---
      ! 1. Buckingham, A. D., Permanent and induced molecular moments
      !    and long-range intermolecular forces. Adv. Chem. Phys., 107, 12 (1967)
      ! 2. Shortley, G., The Computation of Quadrupole and Magnetic-Dipole
      !    Transition probabilities, Phys. Rev., 225, 57 (1940)
      !
      ! Primitive Cartesian:
      ! <a| Qxx |b> = <a| X**2 |b>
      ! <a| Qxy |b> = <a| X * Y |b>
      ! 
      integer, parameter :: MULTI_QUAD_PRIMITIVE = 0
      !
      ! Traceless tensor. Buckingham's traceless quadrupole definition [1]
      ! <a| Qxx |b> = 1/2 <a| 3 * X**2 - (X**2+Y**2+Z**2)|b>
      ! <a| Qxy |b> = 1/2 <a| 3 * X * Y |b>
      ! Note that the minus sign is not included in the definition of
      ! this integral.
      !
      integer, parameter :: MULTI_QUAD_TRACELESS_BUCKINGHAM = 1
      !
      ! Traceless tensor. Shortley's traceless quadrupole definition [2]
      ! eq. (1)
      ! <a| Qxx |b> = <a| X**2 - 1/3*(X**2+Y**2+Z**2)|b>
      ! <a| Qxy |b> = <a| X * Y |b>
      ! Note that the minus sign is not included in the definition of
      ! this integral.
      !
      integer, parameter :: MULTI_QUAD_TRACELESS_SHORTLEY = 2

contains

      subroutine multi_Display(D, Q)
            real(F64), dimension(3), intent(in)    :: D
            real(F64), dimension(3, 3), intent(in) :: Q

            character(:), allocatable :: line
            real(F64) :: x, y, z
            character(:), allocatable :: sx, sy, sz
            integer :: i
            integer, parameter :: w = 15
            character(1), dimension(3), parameter :: Row = ["x", "y", "z"]

            x = todebye(D(1))
            y = todebye(D(2))
            z = todebye(D(3))
            sx = rfield(str(x, d=3), w)
            sy = rfield(str(y, d=3), w)
            sz = rfield(str(z, d=3), w)
            call blankline()
            call msg("dipole moment (Debye)")
            call msg(cfield("x", w) // cfield("y", w) // cfield("z", w))
            line = sx // sy // sz
            call msg(line)
            call blankline()
            call msg("traceless quadrupole moment (Debye*Angs)")
            call msg(cfield("", w) // cfield("x", w) // cfield("y", w) // cfield("z", w))
            do i = 1, 3
                  x = toang(todebye(Q(i, 1)))
                  y = toang(todebye(Q(i, 2)))
                  z = toang(todebye(Q(i, 3)))
                  sx = rfield(str(x, d=3), w)
                  sy = rfield(str(y, d=3), w)
                  sz = rfield(str(z, d=3), w)
                  line = lfield(Row(i), w) // sx // sy // sz
                  call msg(line)
            end do
            call blankline()
      end subroutine multi_Display
      

      subroutine multi_TotalMultipoles(D, Q, Rho, System, AOBasis)
            real(F64), dimension(3), intent(out)      :: D
            real(F64), dimension(3, 3), intent(out)   :: Q
            real(F64), dimension(:, :, :), intent(in) :: Rho
            type(TSystem), intent(in)                 :: System
            type(TAOBasis), intent(in)                :: AOBasis

            real(F64) :: Nx, Ny, Nz
            real(F64) :: Ex, Ey, Ez
            real(F64) :: Nyx, Nzx, Nzy, Nxx, Nyy, Nzz
            real(F64) :: Eyx, Ezx, Ezy, Exx, Eyy, Ezz
            real(F64), dimension(:, :), allocatable :: Dx, Dy, Dz
            real(F64), dimension(:, :), allocatable :: Qyx, Qzx, Qzy, Qxx, Qyy, Qzz
            real(F64), dimension(3) :: Rc
            integer :: NAO, NSpins, s

            NAO = AOBasis%NAOSpher
            NSpins = size(Rho, dim=3)
            call sys_ChargeCenter(Rc, System)
            call sys_NuclearMultipoles(Nx, Ny, Nz, Nyx, Nzx, &
                  Nzy, Nxx, Nyy, Nzz, Rc, System)
            !
            ! Electronic dipole
            !
            allocate(Dx(NAO, NAO))
            allocate(Dy(NAO, NAO))
            allocate(Dz(NAO, NAO))
            call multi_ElectronicDipole(Dx, Dy, Dz, Rc, AOBasis)
            D(1) = Nx
            D(2) = Ny
            D(3) = Nz            
            do s = 1, NSpins
                  call real_vw_x(Ex, Rho(:, :, s), Dx, NAO**2)
                  call real_vw_x(Ey, Rho(:, :, s), Dy, NAO**2)
                  call real_vw_x(Ez, Rho(:, :, s), Dz, NAO**2)
                  D(1) = D(1) - Ex
                  D(2) = D(2) - Ey
                  D(3) = D(3) - Ez
            end do
            !
            ! Electronic quadrupole
            !
            allocate(Qyx(NAO, NAO))
            allocate(Qzx(NAO, NAO))
            allocate(Qzy(NAO, NAO))
            allocate(Qxx(NAO, NAO))
            allocate(Qyy(NAO, NAO))
            allocate(Qzz(NAO, NAO))
            Q(1, 1) = Nxx
            Q(2, 1) = Nyx
            Q(3, 1) = Nzx
            Q(2, 2) = Nyy
            Q(3, 2) = Nzy
            Q(3, 3) = Nzz
            call multi_ElectronicQuadrupole(QYX, QZX, QZY, QXX, QYY, QZZ, Rc, AOBasis)
            do s = 1, NSpins
                  call real_vw_x(Eyx, Rho(:, :, s), QYX, NAO**2)
                  call real_vw_x(Ezx, Rho(:, :, s), QZX, NAO**2)
                  call real_vw_x(Ezy, Rho(:, :, s), QZY, NAO**2)
                  call real_vw_x(Exx, Rho(:, :, s), QXX, NAO**2)
                  call real_vw_x(Eyy, Rho(:, :, s), QYY, NAO**2)
                  call real_vw_x(Ezz, Rho(:, :, s), QZZ, NAO**2)
                  Q(1, 1) = Q(1, 1) - Exx
                  Q(2, 1) = Q(2, 1) - Eyx
                  Q(3, 1) = Q(3, 1) - Ezx
                  Q(2, 2) = Q(2, 2) - Eyy
                  Q(3, 2) = Q(3, 2) - Ezy
                  Q(3, 3) = Q(3, 3) - Ezz
            end do
            Q(1, 2) = Q(2, 1)
            Q(1, 3) = Q(3, 1)
            Q(2, 3) = Q(3, 2)
      end subroutine multi_TotalMultipoles


      subroutine multi_TracelessQuadrupole(QTraceless, QPrimitive, QuadDefinition)
            real(F64), dimension(3, 3), intent(out) :: QTraceless
            real(F64), dimension(3, 3), intent(in)  :: QPrimitive
            integer, intent(in)                     :: QuadDefinition

            real(F64) :: Trace
            integer :: r

            Trace = QPrimitive(1, 1) + QPrimitive(2, 2) + QPrimitive(3, 3)
            if (QuadDefinition == MULTI_QUAD_TRACELESS_BUCKINGHAM) then
                  QTraceless = (THREE/TWO) * QPrimitive
                  do r = 1, 3
                        QTraceless(r, r) = QTraceless(r, r) - (ONE/TWO) * Trace
                  end do
            else if (QuadDefinition == MULTI_QUAD_TRACELESS_SHORTLEY) then
                  QTraceless = QPrimitive
                  do r = 1, 3
                        QTraceless(r, r) = QTraceless(r, r) - (ONE/THREE) * Trace
                  end do                  
            end if
      end subroutine multi_TracelessQuadrupole

      
      subroutine multi_ElectronicDipole_AB(DXab, DYab, DZab, Rpc, ExAB, EyAB, EzAB, Prefactor, La, Lb, &
            Na, Nb, NormA, NormB, CartPolyX, CartPolyY, CartPolyZ)

            real(F64), dimension(Na, Nb), intent(inout) :: DXab, DYab, DZab
            real(F64), dimension(3), intent(in)         :: Rpc
            real(F64), dimension(:), intent(in)         :: ExAB
            real(F64), dimension(:), intent(in)         :: EyAB
            real(F64), dimension(:), intent(in)         :: EzAB
            real(F64), intent(in)                       :: Prefactor
            integer, intent(in)                         :: La
            integer, intent(in)                         :: Lb
            integer, intent(in)                         :: Na
            integer, intent(in)                         :: Nb
            real(F64), dimension(:), intent(in)         :: NormA
            real(F64), dimension(:), intent(in)         :: NormB
            integer, dimension(:, 0:), intent(in)       :: CartPolyX
            integer, dimension(:, 0:), intent(in)       :: CartPolyY
            integer, dimension(:, 0:), intent(in)       :: CartPolyZ

            integer :: a, b
            integer :: x0, y0, z0
            integer :: lxA, lyA, lzA, lxB, lyB, lzB
            real(F64) :: F

            do b = 1, Nb
                  do a = 1, Na
                        lxA = CartPolyX(a, La)
                        lyA = CartPolyY(a, La)
                        lzA = CartPolyZ(a, La)
                        lxB = CartPolyX(b, Lb)
                        lyB = CartPolyY(b, Lb)
                        lzB = CartPolyZ(b, Lb)

                        x0 = ints1e_ELoc(lxA, lxB, 0)
                        y0 = ints1e_ELoc(lyA, lyB, 0)
                        z0 = ints1e_ELoc(lzA, lzB, 0)

                        F = Prefactor * NormA(a) * NormB(b)
                        !
                        ! Helgaker, T., Jorgensen, P., Olsen, J., Section 9.5.3 in 
                        ! Molecular Electronic-Structure Theory, John Wiley and Sons 2000.
                        !
                        if (lxA+lxB > 0) then
                              DXab(a, b) = DXab(a, b) + F * (ExAB(x0+1) + Rpc(1) * ExAB(x0)) * EyAB(y0) * EzAB(z0)
                        else
                              DXab(a, b) = DXab(a, b) + F * Rpc(1) * ExAB(x0) * EyAB(y0) * EzAB(z0)
                        end if
                        if (lyA+lyB > 0) then
                              DYab(a, b) = DYab(a, b) + F * ExAB(x0) * (EyAB(y0+1) + Rpc(2) * EyAB(y0)) * EzAB(z0)
                        else
                              DYab(a, b) = DYab(a, b) + F * ExAB(x0) * Rpc(2) * EyAB(y0) * EzAB(z0)
                        end if
                        if (lzA+lzB > 0) then
                              DZab(a, b) = DZab(a, b) + F * ExAB(x0) * EyAB(y0) * (EzAB(z0+1) + Rpc(3) * EzAB(z0))
                        else
                              DZab(a, b) = DZab(a, b) + F * ExAB(x0) * EyAB(y0) * Rpc(3) * EzAB(z0)
                        end if
                  end do
            end do
      end subroutine multi_ElectronicDipole_AB


      subroutine multi_ElectronicDipole_Cartesian(Dx, Dy, Dz, Rc, AOBasis)
            !
            real(F64), dimension(:, :), intent(out) :: Dx, Dy, Dz
            real(F64), dimension(3), intent(in)     :: Rc
            type(TAOBasis), intent(in)              :: AOBasis

            integer :: ShellA, ShellB, ShellParamsA, ShellParamsB, ShellAB
            integer :: La, Na, a0, a1
            integer :: Lb, Nb, b0, b1
            integer :: i, j
            integer :: Lab
            real(F64), dimension(3) :: Ra, Rb, Rp, Rpa, Rpb, Rpc, Rab
            real(F64) :: Prefactor, AlphaA, AlphaB, AlphaAB, Mu
            integer, parameter :: MaxNFunc = ((AUTO2E_MAXL + 1) * (AUTO2E_MAXL + 2)) / 2
            real(F64), dimension(MaxNFunc**2) :: DXab, DYab, DZab
            integer, parameter :: MaxIndex = 2 * AUTO2E_MAXL
            !
            ! Dimension = Sum(k = 0, MaxIndex) (k + 1)**2
            !
            real(F64), dimension((2*MaxIndex**3+9*MaxIndex**2+13*MaxIndex+6)/6) :: ExAB, EyAB, EzAB

            Dx = ZERO
            Dy = ZERO
            Dz = ZERO
            associate ( &
                  ShellParamsIdx => AOBasis%ShellParamsIdx, &
                  ShellMomentum => AOBasis%ShellMomentum, &
                  NAngFunc => AOBasis%NAngFuncCart, &
                  ShellLoc => AOBasis%ShellLocCart, &
                  NShells => AOBasis%NShells, &
                  CartPolyX => AOBasis%CartPolyX, &
                  CartPolyY => AOBasis%CartPolyY, &
                  CartPolyZ => AOBasis%CartPolyZ, &
                  LmaxGTO => AOBasis%LmaxGTO, &
                  ShellCenters => AOBasis%ShellCenters, &
                  AtomCoords => AOBasis%AtomCoords, &
                  CntrCoeffs => AOBasis%CntrCoeffs, &
                  Exponents => AOBasis%Exponents, &
                  NPrimitives => AOBasis%NPrimitives, &
                  NormFactors => AOBasis%NormFactorsCart &
                  )
                  !$omp parallel &
                  !$omp private(ShellA, ShellParamsA, La, Na, a0, a1, Ra) &
                  !$omp private(ShellB, ShellParamsB, Lb, Nb, b0, b1, Rb) &
                  !$omp private(ShellAB) &
                  !$omp private(AlphaA, AlphaB, Mu, Prefactor) &
                  !$omp private(ExAB, EyAB, EzAB) &
                  !$omp private(i, j) &
                  !$omp private(Lab, AlphaAB) &
                  !$omp private(DXab, DYab, DZab) &
                  !$omp private(Rp, Rpa, Rpb, Rpc, Rab) &
                  !$omp shared(DX, DY, DZ) &
                  !$omp default(shared)
                  !$omp do schedule(dynamic)
                  ShellsAB: do ShellAB = 1, (NShells * (NShells + 1)) / 2
                        call ints1e_decode_pq(ShellAB, NShells, ShellA, ShellB)
                        ShellParamsA = ShellParamsIdx(ShellA)
                        La = ShellMomentum(ShellParamsA)
                        Na = NAngFunc(ShellParamsA)                        
                        a0 = ShellLoc(ShellA)
                        a1 = ShellLoc(ShellA) + Na - 1

                        ShellParamsB = ShellParamsIdx(ShellB)
                        Lb = ShellMomentum(ShellParamsB)
                        Nb = NAngFunc(ShellParamsB)
                        b0 = ShellLoc(ShellB)
                        b1 = ShellLoc(ShellB) + Nb - 1

                        Ra = AtomCoords(:, ShellCenters(ShellA))
                        Rb = AtomCoords(:, ShellCenters(ShellB))
                        Rab = Ra - Rb

                        Lab = La + Lb
                        DXab = ZERO
                        DYab = ZERO
                        DZab = ZERO
                        PrimitivesA: do j = 1, NPrimitives(ShellParamsB)
                              PrimitivesB: do i = 1, NPrimitives(ShellParamsA)
                                    AlphaA = Exponents(i, ShellParamsA)
                                    AlphaB = Exponents(j, ShellParamsB)
                                    AlphaAB = AlphaA + AlphaB
                                    Mu = AlphaA * AlphaB / (AlphaA + AlphaB)
                                    Prefactor = CntrCoeffs(i, ShellParamsA) * CntrCoeffs(j, ShellParamsB) &
                                          * exp(-Mu * dot_product(Rab, Rab)) * (PI / AlphaAB)**(THREE/TWO)
                                    Rp = (AlphaA*Ra + AlphaB*Rb) / AlphaAB
                                    Rpa = Rp - Ra
                                    Rpb = Rp - Rb
                                    Rpc = Rp - Rc                                    
                                    !
                                    ! Seed for calculating E^{ij}_t coeffs is 1.0 instead of
                                    ! exp(-Mu * Xab**2) as in Helgaker's textbook
                                    ! because the coeffs are linear functions of the seed so
                                    ! it may be incorporated into the prefactor.
                                    !
                                    call ints1e_eijmatrix(Lab, ONE, AlphaAB, Rpa(1), Rpb(1), ExAB)
                                    call ints1e_eijmatrix(Lab, ONE, AlphaAB, Rpa(2), Rpb(2), EyAB)
                                    call ints1e_eijmatrix(Lab, ONE, AlphaAB, Rpa(3), Rpb(3), EzAB)
                                    call multi_ElectronicDipole_AB(DXab, DYab, DZab, Rpc, &
                                          ExAB, EyAB, EzAB, Prefactor, La, Lb, &
                                          Na, Nb, NormFactors(:, ShellParamsA), NormFactors(:, ShellParamsB), &
                                          CartPolyX, CartPolyY, CartPolyZ)
                              end do PrimitivesB
                        end do PrimitivesA
                        DX(a0:a1, b0:b1) = reshape(DXab(1:Na*Nb), [Na, Nb])
                        DY(a0:a1, b0:b1) = reshape(DYab(1:Na*Nb), [Na, Nb])
                        DZ(a0:a1, b0:b1) = reshape(DZab(1:Na*Nb), [Na, Nb])
                  end do ShellsAB
                  !$omp end do
                  !$omp end parallel
            end associate
      end subroutine multi_ElectronicDipole_Cartesian


      subroutine multi_ElectronicDipole(Dx, Dy, Dz, Rc, AOBasis)
            !
            ! Electronic dipole moment matrix in the basis of
            ! spherical Gaussian atomic orbitals.
            !
            ! Both upper and lower triangles of the output matrix
            ! are filled with data.
            !
            real(F64), dimension(:, :), intent(out) :: Dx, Dy, Dz
            real(F64), dimension(3), intent(in)     :: Rc
            type(TAOBasis), intent(in)              :: AOBasis

            real(F64), dimension(:, :), allocatable :: Dx_cao, Dy_cao, Dz_cao
            integer :: NAOCart

            NAOCart = AOBasis%NAOCart
            allocate(Dx_cao(NAOCart, NAOCart))
            allocate(Dy_cao(NAOCart, NAOCart))
            allocate(Dz_cao(NAOCart, NAOCart))
            call multi_ElectronicDipole_Cartesian(Dx_cao, Dy_cao, Dz_cao, Rc, AOBasis)
            call real_smfill(Dx_cao)
            call real_smfill(Dy_cao)
            call real_smfill(Dz_cao)
            call ints1e_SpherAOTransf(Dx, Dx_cao, AOBasis)
            call ints1e_SpherAOTransf(Dy, Dy_cao, AOBasis)
            call ints1e_SpherAOTransf(Dz, Dz_cao, AOBasis)
      end subroutine multi_ElectronicDipole


      subroutine multi_ElectronicQuadrupole_AB(QYXab, QZXab, QZYab, QXXab, QYYab, QZZab, &
            Rpc, ExAB, EyAB, EzAB, Prefactor, La, Lb, &
            Na, Nb, NormA, NormB, CartPolyX, CartPolyY, CartPolyZ, &
            AlphaAB)

            real(F64), dimension(Na, Nb), intent(inout) :: QYXab, QZXab, QZYab, QXXab, QYYab, QZZab
            real(F64), dimension(3), intent(in)         :: Rpc
            real(F64), dimension(:), intent(in)         :: ExAB
            real(F64), dimension(:), intent(in)         :: EyAB
            real(F64), dimension(:), intent(in)         :: EzAB
            real(F64), intent(in)                       :: Prefactor
            integer, intent(in)                         :: La
            integer, intent(in)                         :: Lb
            integer, intent(in)                         :: Na
            integer, intent(in)                         :: Nb
            real(F64), dimension(:), intent(in)         :: NormA
            real(F64), dimension(:), intent(in)         :: NormB
            integer, dimension(:, 0:), intent(in)       :: CartPolyX
            integer, dimension(:, 0:), intent(in)       :: CartPolyY
            integer, dimension(:, 0:), intent(in)       :: CartPolyZ
            real(F64), intent(in)                       :: AlphaAB

            integer :: a, b
            integer :: x0, y0, z0
            integer :: lxA, lyA, lzA, lxB, lyB, lzB
            real(F64) :: F
            real(F64) :: ExAB_0, ExAB_1, ExAB_2
            real(F64) :: EyAB_0, EyAB_1, EyAB_2
            real(F64) :: EzAB_0, EzAB_1, EzAB_2

            do b = 1, Nb
                  do a = 1, Na
                        lxA = CartPolyX(a, La)
                        lyA = CartPolyY(a, La)
                        lzA = CartPolyZ(a, La)
                        lxB = CartPolyX(b, Lb)
                        lyB = CartPolyY(b, Lb)
                        lzB = CartPolyZ(b, Lb)

                        x0 = ints1e_ELoc(lxA, lxB, 0)
                        y0 = ints1e_ELoc(lyA, lyB, 0)
                        z0 = ints1e_ELoc(lzA, lzB, 0)

                        F = Prefactor * NormA(a) * NormB(b)
                        !
                        ! Helgaker, T., Jorgensen, P., Olsen, J., Section 9.5.3 in 
                        ! Molecular Electronic-Structure Theory, John Wiley and Sons 2000.
                        !
                        ExAB_0 = ExAB(x0)
                        EyAB_0 = EyAB(y0)
                        EzAB_0 = EzAB(z0)
                        if (lxA+lxB >= 2) then
                              ExAB_1 = ExAB(x0+1)
                              ExAB_2 = ExAB(x0+2)
                        else if (lxA+lxB == 1) then
                              ExAB_1 = ExAB(x0+1)
                              ExAB_2 = ZERO
                        else
                              ExAB_1 = ZERO
                              ExAB_2 = ZERO
                        end if
                        if (lyA+lyB >= 2) then
                              EyAB_1 = EyAB(y0+1)
                              EyAB_2 = EyAB(y0+1)
                        else if (lyA+lyB == 1) then
                              EyAB_1 = EyAB(y0+1)
                              EyAB_2 = ZERO
                        else
                              EyAB_1 = ZERO
                              EyAB_2 = ZERO
                        end if
                        if (lzA+lzB >= 2) then
                              EzAB_1 = EzAB(z0+1)
                              EzAB_2 = EzAB(z0+2)
                        else if (lzA+lzB == 1) then
                              EzAB_1 = EzAB(z0+1)
                              EzAB_2 = ZERO
                        else
                              EzAB_1 = ZERO
                              EzAB_2 = ZERO
                        end if

                        QYXab(a, b) = QYXab(a, b) + F * (ExAB_1 + Rpc(1)*ExAB_0) * (EyAB_1 + Rpc(2)*EyAB_0) * EzAB_0
                        QZXab(a, b) = QZXab(a, b) + F * (ExAB_1 + Rpc(1)*ExAB_0) * EyAB_0 * (EzAB_1 + Rpc(3)*EzAB_0)
                        QZYab(a, b) = QZYab(a, b) + F * ExAB_0 * (EyAB_1 + Rpc(2)*EyAB_0) * (EzAB_1 + Rpc(3)*EzAB_0)
                        QXXab(a, b) = QXXab(a, b) &
                              + F * ((Rpc(1)**2 + (ONE/TWO)/AlphaAB)*ExAB_0 + TWO*Rpc(1)*ExAB_1 + TWO*ExAB_2) &
                              * EyAB_0 &
                              * EzAB_0
                        QYYab(a, b) = QYYab(a, b) &
                              + F * ExAB_0 &
                              * ((Rpc(2)**2 + (ONE/TWO)/AlphaAB)*EyAB_0 + TWO*Rpc(2)*EyAB_1 + TWO*EyAB_2) &
                              * EzAB_0
                        QZZab(a, b) = QZZab(a, b) &
                              + F * ExAB_0 &
                              * EyAB_0 &
                              * ((Rpc(3)**2 + (ONE/TWO)/AlphaAB)*EzAB_0 + TWO*Rpc(3)*EzAB_1 + TWO*EzAB_2)
                  end do
            end do
      end subroutine multi_ElectronicQuadrupole_AB


      subroutine multi_ElectronicQuadrupole_Cartesian(QYX, QZX, QZY, QXX, QYY, QZZ, &
            Rc, AOBasis)
            
            real(F64), dimension(:, :), intent(out) :: QYX, QZX, QZY, QXX, QYY, QZZ
            real(F64), dimension(3), intent(in)     :: Rc
            type(TAOBasis), intent(in)              :: AOBasis

            integer :: ShellA, ShellB, ShellParamsA, ShellParamsB, ShellAB
            integer :: La, Na, a0, a1
            integer :: Lb, Nb, b0, b1
            integer :: i, j
            integer :: Lab
            real(F64), dimension(3) :: Ra, Rb, Rp, Rpa, Rpb, Rpc, Rab
            real(F64) :: Prefactor, AlphaA, AlphaB, AlphaAB, Mu
            integer, parameter :: MaxNFunc = ((AUTO2E_MAXL + 1) * (AUTO2E_MAXL + 2)) / 2
            real(F64), dimension(MaxNFunc**2) :: QYXab, QZXab, QZYab, QXXab, QYYab, QZZab
            integer, parameter :: MaxIndex = 2 * AUTO2E_MAXL
            !
            ! Dimension = Sum(k = 0, MaxIndex) (k + 1)**2
            !
            real(F64), dimension((2*MaxIndex**3+9*MaxIndex**2+13*MaxIndex+6)/6) :: ExAB, EyAB, EzAB

            QYX = ZERO
            QZX = ZERO
            QZY = ZERO
            QXX = ZERO
            QYY = ZERO
            QZZ = ZERO
            associate ( &
                  ShellParamsIdx => AOBasis%ShellParamsIdx, &
                  ShellMomentum => AOBasis%ShellMomentum, &
                  NAngFunc => AOBasis%NAngFuncCart, &
                  ShellLoc => AOBasis%ShellLocCart, &
                  NShells => AOBasis%NShells, &
                  CartPolyX => AOBasis%CartPolyX, &
                  CartPolyY => AOBasis%CartPolyY, &
                  CartPolyZ => AOBasis%CartPolyZ, &
                  LmaxGTO => AOBasis%LmaxGTO, &
                  ShellCenters => AOBasis%ShellCenters, &
                  AtomCoords => AOBasis%AtomCoords, &
                  CntrCoeffs => AOBasis%CntrCoeffs, &
                  Exponents => AOBasis%Exponents, &
                  NPrimitives => AOBasis%NPrimitives, &
                  NormFactors => AOBasis%NormFactorsCart &
                  )
                  !$omp parallel &
                  !$omp private(ShellA, ShellParamsA, La, Na, a0, a1, Ra) &
                  !$omp private(ShellB, ShellParamsB, Lb, Nb, b0, b1, Rb) &
                  !$omp private(ShellAB) &
                  !$omp private(AlphaA, AlphaB, Mu, Prefactor) &
                  !$omp private(ExAB, EyAB, EzAB) &
                  !$omp private(i, j) &
                  !$omp private(Lab, AlphaAB) &
                  !$omp private(QYXab, QZXab, QZYab, QXXab, QYYab, QZZab) &
                  !$omp private(Rp, Rpa, Rpb, Rpc, Rab) &
                  !$omp shared(QYX, QZX, QZY, QXX, QYY, QZZ) &
                  !$omp default(shared)
                  !$omp do schedule(dynamic)
                  ShellsAB: do ShellAB = 1, (NShells * (NShells + 1)) / 2
                        call ints1e_decode_pq(ShellAB, NShells, ShellA, ShellB)
                        ShellParamsA = ShellParamsIdx(ShellA)
                        La = ShellMomentum(ShellParamsA)
                        Na = NAngFunc(ShellParamsA)                        
                        a0 = ShellLoc(ShellA)
                        a1 = ShellLoc(ShellA) + Na - 1

                        ShellParamsB = ShellParamsIdx(ShellB)
                        Lb = ShellMomentum(ShellParamsB)
                        Nb = NAngFunc(ShellParamsB)
                        b0 = ShellLoc(ShellB)
                        b1 = ShellLoc(ShellB) + Nb - 1

                        Ra = AtomCoords(:, ShellCenters(ShellA))
                        Rb = AtomCoords(:, ShellCenters(ShellB))
                        Rab = Ra - Rb

                        Lab = La + Lb
                        QYXab = ZERO
                        QZXab = ZERO
                        QZYab = ZERO
                        QXXab = ZERO
                        QYYab = ZERO
                        QZZab = ZERO
                        PrimitivesA: do j = 1, NPrimitives(ShellParamsB)
                              PrimitivesB: do i = 1, NPrimitives(ShellParamsA)
                                    AlphaA = Exponents(i, ShellParamsA)
                                    AlphaB = Exponents(j, ShellParamsB)
                                    AlphaAB = AlphaA + AlphaB
                                    Mu = AlphaA * AlphaB / (AlphaA + AlphaB)
                                    Prefactor = CntrCoeffs(i, ShellParamsA) * CntrCoeffs(j, ShellParamsB) &
                                          * exp(-Mu * dot_product(Rab, Rab)) * (PI / AlphaAB)**(THREE/TWO)
                                    Rp = (AlphaA*Ra + AlphaB*Rb) / AlphaAB
                                    Rpa = Rp - Ra
                                    Rpb = Rp - Rb
                                    Rpc = Rp - Rc
                                    !
                                    ! Seed for calculating E^{ij}_t coeffs is 1.d+0 instead of
                                    ! exp(-alpha_reduced * xabsq) as in Helgaker's textbook
                                    ! because the coeffs are linear functions of the seed so
                                    ! it may be incorporated into the prefactor.
                                    !
                                    call ints1e_eijmatrix(Lab, ONE, AlphaAB, Rpa(1), Rpb(1), ExAB)
                                    call ints1e_eijmatrix(Lab, ONE, AlphaAB, Rpa(2), Rpb(2), EyAB)
                                    call ints1e_eijmatrix(Lab, ONE, AlphaAB, Rpa(3), Rpb(3), EzAB)
                                    call multi_ElectronicQuadrupole_AB(QYXab, QZXab, QZYab, QXXab, QYYab, QZZab, Rpc, &
                                          ExAB, EyAB, EzAB, Prefactor, La, Lb, &
                                          Na, Nb, NormFactors(:, ShellParamsA), NormFactors(:, ShellParamsB), &
                                          CartPolyX, CartPolyY, CartPolyZ, AlphaAB)
                              end do PrimitivesB
                        end do PrimitivesA
                        QYX(a0:a1, b0:b1) = reshape(QYXab(1:Na*Nb), [Na, Nb])
                        QZX(a0:a1, b0:b1) = reshape(QZXab(1:Na*Nb), [Na, Nb])
                        QZY(a0:a1, b0:b1) = reshape(QZYab(1:Na*Nb), [Na, Nb])
                        QXX(a0:a1, b0:b1) = reshape(QXXab(1:Na*Nb), [Na, Nb])
                        QYY(a0:a1, b0:b1) = reshape(QYYab(1:Na*Nb), [Na, Nb])
                        QZZ(a0:a1, b0:b1) = reshape(QZZab(1:Na*Nb), [Na, Nb])
                  end do ShellsAB
                  !$omp end do
                  !$omp end parallel
            end associate
      end subroutine multi_ElectronicQuadrupole_Cartesian


      subroutine multi_ElectronicQuadrupole(QYX, QZX, QZY, QXX, QYY, QZZ, Rc, AOBasis)
            !
            ! Electronic quadrupole moment matrix in the basis of
            ! spherical Gaussian atomic orbitals.
            !
            ! Both upper and lower triangles of the output matrix
            ! are filled with data.
            !
            real(F64), dimension(:, :), intent(out) :: QYX, QZX, QZY, QXX, QYY, QZZ
            real(F64), dimension(3), intent(in)     :: Rc
            type(TAOBasis), intent(in)              :: AOBasis

            real(F64), dimension(:, :), allocatable :: QYX_cao, QZX_cao, QZY_cao
            real(F64), dimension(:, :), allocatable :: QXX_cao, QYY_cao, QZZ_cao
            integer :: NAOCart, NAO

            NAOCart = AOBasis%NAOCart
            NAO = AOBasis%NAOSpher
            allocate(QYX_cao(NAOCart, NAOCart))
            allocate(QZX_cao(NAOCart, NAOCart))
            allocate(QZY_cao(NAOCart, NAOCart))
            allocate(QXX_cao(NAOCart, NAOCart))
            allocate(QYY_cao(NAOCart, NAOCart))
            allocate(QZZ_cao(NAOCart, NAOCart))
            call multi_ElectronicQuadrupole_Cartesian(QYX_cao, QZX_cao, QZY_cao, &
                  QXX_cao, QYY_cao, QZZ_cao, &
                  Rc, AOBasis)
            call real_smfill(QYX_cao)
            call real_smfill(QZX_cao)
            call real_smfill(QZY_cao)
            call real_smfill(QXX_cao)
            call real_smfill(QYY_cao)
            call real_smfill(QZZ_cao)
            call ints1e_SpherAOTransf(QYX, QYX_cao, AOBasis)
            call ints1e_SpherAOTransf(QZX, QZX_cao, AOBasis)
            call ints1e_SpherAOTransf(QZY, QZY_cao, AOBasis)
            call ints1e_SpherAOTransf(QXX, QXX_cao, AOBasis)
            call ints1e_SpherAOTransf(QYY, QYY_cao, AOBasis)
            call ints1e_SpherAOTransf(QZZ, QZZ_cao, AOBasis)
      end subroutine multi_ElectronicQuadrupole
end module Multipoles
