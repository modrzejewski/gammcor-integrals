module Auto2eInterface
      use arithmetic
      use basis_sets
      use Auto2e
      use display
      use string
!      use OneElectronInts
      use sys_definitions
      use linalg
      use spherh
      
      implicit none

      integer, parameter :: ORBITAL_ORDERING_AUTO2E = 0
      integer, parameter :: ORBITAL_ORDERING_MOLPRO = 1
      integer, parameter :: ORBITAL_ORDERING_ORCA = 2
      integer, parameter :: ORBITAL_ORDERING_DALTON = 3

      integer, parameter, private :: Auto2e_MaxAngFuncSpher = 2 * Auto2e_MaxL + 1
      integer, parameter, private :: Auto2e_MaxAngFuncCart = ((Auto2e_MaxL+1)*(Auto2e_MaxL+2))/2

      integer, private :: k
      ! --------------------------------------------------------------------------------------
      ! Ordering of solid harmonics in Orca
      ! --------------------------------------------------------------------------------------
      !                                                                  s
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: ORCA_S = [0, (-1, k=1,10)] + 1
      !
      ! Note the exception: p orbitals are not transformed to the spherical basis by the subroutines
      ! in the Auto2e module, see the SPHER_TRANSF_LMIN parameter. The Cartesian ordering, (px, py, pz),
      ! is the same as in Molpro, but different from Orca's. 
      !                                                                  pz px  py
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: ORCA_P = [3, 1, 2, (0,k=1,8)]
      !                                                                  d0  d1+  d1- d2+  d2-
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: ORCA_D = [0,  1,   -1,  2,   -2,  (-3,k=1,6)] + 3
      !                                                                 f0  f1+  f1-  f2+ f2-  f3+ f3-
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: ORCA_F = [0,  1,   -1,   2,  -2,  3,  -3, (-4,k=1,4)] + 4
      !                                                                  g0  g1+  g1-  g2+  g2-  g3+  g3- g4+   g4-
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: ORCA_G = [0,  1,   -1,  2,   -2,  3,   -3,  4,   -4, (-5,k=1,2)] + 5            
      !                                                                  h0   h1+  h1-  h2+  h2-  h3+  h3-  h4+  h4-  h5+  h5-
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: ORCA_H = [0,   1,  -1,   2,   -2,  3,   -3,  4,   -4,  5,   -5] + 6
      ! --------------------------------------------------------------------------------------
      ! Ordering of solid harmonics in Molpro
      ! --------------------------------------------------------------------------------------
      !                                                                    s
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: MOLPRO_S = [0, (-1, k=1,10)] + 1
      !
      ! Note the exception: p orbitals are not transformed to the spherical basis by the subroutines
      ! in the Auto2e module, see the SPHER_TRANSF_LMIN parameter. The Cartesian ordering, (px, py, pz),
      ! is the same as in Molpro. (Otherwise, the pz orbital corresponding to m=0 would be stored in
      ! the second element of the array instead of the third.)
      !
      !                                                                    px py pz
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: MOLPRO_P = [1, 2, 3, (0,k=1,8)]
      !                                                                    d0  d2- d1+ d2+ d1-
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: MOLPRO_D = [0,  -2,  1, 2, -1, (-3,k=1,6)] + 3
      !                                                                    f1+  f1-  f0  f3+  f2-  f3- f2+
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: MOLPRO_F = [1,   -1,   0,  3,  -2,  -3,  2, (-4,k=1,4)] + 4
      !                                                                    g0  g2-  g1+  g4+  g1-  g2+  g4-  g3+  g3-
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: MOLPRO_G = [0,  -2,  1,   4,   -1,  2,   -4,  3,   -3, (-5,k=1,2)] + 5
      !                                                                    h1+  h1-  h2+  h3+  h4-  h3-  h4+  h5-  h0  h5+  h2-
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: MOLPRO_H = [1,   -1,  2,   3,   -4,  -3,   4,  -5,   0,  5,   -2] + 6
      !
      ! ------------------------------------------------------------------------------------------------------------
      ! Ordering of solid harmonics in Dalton
      ! ------------------------------------------------------------------------------------------------------------
      !                                                                    s
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: DALTON_S = [0, (-1, k=1,10)] + 1
      !
      ! Note the exception: p orbitals are not transformed to the spherical basis by the subroutines
      ! in the Auto2e module, see the SPHER_TRANSF_LMIN parameter. The Cartesian ordering, (px, py, pz),
      ! is the same as in Molpro. (Otherwise, the pz orbital corresponding to m=0 would be stored in
      ! the second element of the array instead of the third.)
      !                                                                    px py pz
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: DALTON_P = [1, 2, 3, (0,k=1,8)]
      !
      ! Higher-order solid harmonics in Dalton are ordered from -L to L
      !
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: DALTON_D = [(k,k=-2,2), (-3,k=1,6)] + 3
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: DALTON_F = [(k,k=-3,3), (-4,k=1,4)] + 4
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: DALTON_G = [(k,k=-4,4), (-5,k=1,2)] + 5
      integer, dimension(Auto2e_MaxAngFuncSpher), parameter :: DALTON_H = [(k,k=-5,5)] + 6
      !
      ! ------------------------------------------------------------------------------------------------------------
      ! Ordering of Cartesian functions in Molpro
      ! ------------------------------------------------------------------------------------------------------------
      !                                                                        s
      integer, dimension(Auto2e_MaxAngFuncCart), parameter :: MOLPRO_CART_S = [1,(0,k=1,Auto2e_MaxAngFuncCart-1)]
      !                                                                        x  y  z
      integer, dimension(Auto2e_MaxAngFuncCart), parameter :: MOLPRO_CART_P = [1, 2, 3, (0,k=1,Auto2e_MaxAngFuncCart-3)]
      !
      ! Auto2e ordering
      ! 1  xx
      ! 2  xy
      ! 3  xz
      ! 4  yy
      ! 5  yz
      ! 6  zz
      !                                                                        xx yy zz xy xz yz
      integer, dimension(Auto2e_MaxAngFuncCart), parameter :: MOLPRO_CART_D = [1, 4, 6, 2, 3, 5, (0,k=1,Auto2e_MaxAngFuncCart-6)]
      !
      ! Auto2e ordering
      !
      ! 1  xxx
      ! 2  xxy
      ! 3  xxz
      ! 4  xyy
      ! 5  xyz
      ! 6  xzz
      ! 7  yyy
      ! 8  yyz
      ! 9  yzz
      ! 10 zzz
      !
      integer, dimension(Auto2e_MaxAngFuncCart), parameter :: MOLPRO_CART_F = &
            ! xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz
            [ 1,  7,  10, 2,  3,  4,  8,  6,  9,  5, (0,k=1,Auto2e_MaxAngFuncCart-10)]
      !
      ! Auto2e ordering
      !
      ! 1  xxxx
      ! 2  xxxy
      ! 3  xxxz
      ! 4  xxyy
      ! 5  xxyz
      ! 6  xxzz
      ! 7  xyyy
      ! 8  xyyz
      ! 9  xyzz
      ! 10 xzzz
      ! 11 yyyy
      ! 12 yyyz
      ! 13 yyzz
      ! 14 yzzz
      ! 15 zzzz
      !
      integer, dimension(Auto2e_MaxAngFuncCart), parameter :: MOLPRO_CART_G = &
            ! xxxx yyyy zzzz xxxy xxxz xyyy yyyz xzzz yzzz xxyy xxzz yyzz xxyz xyyz xyzz
            [ 1,   11,  15,  2,   3,   7,   12,  10,  14,  4,   6,   13,  5,   8,   9, (0,k=1,Auto2e_MaxAngFuncCart-15)]
      !
      ! And finally, the ordring of h functions in Molpro is systematic and the same as in Auto2e.
      !
      integer, dimension(Auto2e_MaxAngFuncCart), parameter :: MOLPRO_CART_H = [(k,k=1,Auto2e_MaxAngFuncCart)]
      
contains

      subroutine auto2e_interface_AngFuncTransf(X_out, X_in, FromExternalAO, TwoIndexTransf, &
            AOBasis, ExternalOrdering)
            !
            ! Convert AO matrices/vectors between external AO and Auto2e formats.
            ! The direction the conversion depends on the parameter FromExternalAO
            !
            ! X_out
            !                  Output, converted matrix.
            ! X_in
            !                  Input, matrix before conversion.
            ! FromExternalAO
            !                  If false: X_in(Auto2e) -> X_out(External program)
            !                  If true: X_in(External program) -> X_out(Auto2e)
            ! TwoIndexTransf
            !                  If true, both row and column indices are transformed.
            !                  This should be used for the one-electron hamiltonian,
            !                  overlap matrix, etc. If false, only the row index of X
            !                  is transformed. This option should be used for MO
            !                  coefficients in the AO basis stored in columns.
            ! AOBasis
            !                  Basis set data.
            ! ExternalOrdering
            !                  Speficies the external program, e.g., Molpro, Dalton, Orca.
            ! 
            !
            real(F64), dimension(:, :), intent(out) :: X_out
            real(F64), dimension(:, :), intent(in)  :: X_in
            logical, intent(in)                     :: FromExternalAO            
            logical, intent(in)                     :: TwoIndexTransf
            type(TAOBasis), intent(in)              :: AOBasis
            integer, intent(in)                     :: ExternalOrdering

            integer, dimension(:), allocatable :: Map
            integer, dimension(:), allocatable :: Map_ao2extao
            integer :: p, n
            integer :: NAO

            NAO = size(X_in, dim=1)
            allocate(Map(NAO))
            allocate(Map_ao2extao(NAO))
            call auto2e_interface_AOMap(Map, Map_ao2extao, AOBasis, ExternalOrdering)
            if (.not. FromExternalAO) then
                  Map = Map_ao2extao
            end if            
            if (TwoIndexTransf) then
                  !
                  ! Examples: overlap matrix, one-electron matrices of the Fock operator, ...
                  !
                  call auto2e_interface_TransfMatrix(X_out, X_in, Map)
            else
                  !
                  ! Examples: MO coefficients stored in columns, e.g., C(p,i), where p is AO and i is MO
                  !
                  n = size(X_in, dim=2)
                  do p = 1, n
                        call auto2e_interface_TransfVector(X_out(:, p), X_in(:, p), Map)
                  end do
            end if
      end subroutine auto2e_interface_AngFuncTransf
      

      subroutine auto2e_interface_TransfMatrix(X_out, X_in, Map)
            real(F64), dimension(:, :), intent(out) :: X_out
            real(F64), dimension(:, :), intent(in)  :: X_in
            integer, dimension(:), intent(in)       :: Map
            
            integer :: q
            integer :: NAO
            real(F64), dimension(:, :), allocatable :: W
            
            NAO = size(Map)
            allocate(W(NAO, NAO))
            do q = 1, NAO
                  call auto2e_interface_TransfVector(W(:, q), X_in(:, q), Map)
            end do
            do q = 1, NAO
                  X_out(:, Map(q)) = W(:, q)
            end do
      end subroutine auto2e_interface_TransfMatrix


      subroutine auto2e_interface_TransfVector(X_out, X_in, Map)
            real(F64), dimension(:), intent(out) :: X_out
            real(F64), dimension(:), intent(in)  :: X_in
            integer, dimension(:), intent(in)    :: Map
            
            integer :: p
            integer :: NAO
            
            NAO = size(Map)
            do p = 1, NAO
                  X_out(Map(p)) = X_in(p)
            end do
      end subroutine auto2e_interface_TransfVector


      subroutine auto2e_interface_C(C_ao, C_extao, AOBasis, ExternalOrdering)
            real(F64), dimension(:, :), intent(out) :: C_ao
            real(F64), dimension(:, :), intent(in)  :: C_extao
            type(TAOBasis), intent(in)              :: AOBasis
            integer, intent(in)                     :: ExternalOrdering
            
            logical :: FromExternalAO
            logical :: TwoIndexTransf

            FromExternalAO = .true. ! AOs from external program -> AOs in the Auto2e format
            TwoIndexTransf = .false. ! Transform only the index p of C(p,k)
            call auto2e_interface_AngFuncTransf(C_ao, C_extao, FromExternalAO, TwoIndexTransf, AOBasis, ExternalOrdering)
            if (ExternalOrdering == ORBITAL_ORDERING_ORCA) then
                  call auto2e_interface_ApplyOrcaPhases_Matrix(C_ao, AOBasis, TwoIndexTransf)
            end if
      end subroutine auto2e_interface_C


      subroutine auto2e_interface_TransfMOCoeffs(X_out, X_in, Map)
            real(F64), dimension(:, :), intent(out) :: X_out
            real(F64), dimension(:, :), intent(in)  :: X_in
            integer, dimension(:), intent(in)       :: Map

            integer :: m, k
            
            m = size(X_in, dim=1)
            do k = 1, m
                  call auto2e_interface_TransfVector(X_out(:, k), X_in(:, k), Map)
            end do
      end subroutine auto2e_interface_TransfMOCoeffs
      
      
      subroutine auto2e_interface_AOMap(Map_extao2ao, Map_ao2extao, AOBasis, ExternalOrdering)
            integer, dimension(:), intent(out) :: Map_extao2ao
            integer, dimension(:), intent(out) :: Map_ao2extao
            type(TAOBasis), intent(in)         :: AOBasis
            integer, intent(in)                :: ExternalOrdering

            integer :: k
            !
            ! Ordering of real solid harmonics in Molpro
            !
            integer, parameter :: Auto2e_MaxAngFuncSpher = 2 * Auto2e_MaxL + 1
            integer, parameter :: Auto2e_MaxAngFuncCart = ((Auto2e_MaxL+1)*(Auto2e_MaxL+2))/2

            integer, dimension(Auto2e_MaxAngFuncSpher, 0:Auto2e_MaxL), parameter :: Molpro_AngFuncMap &
                  = reshape([MOLPRO_S, MOLPRO_P, MOLPRO_D,  MOLPRO_F, MOLPRO_G, MOLPRO_H], [Auto2e_MaxAngFuncSpher,Auto2e_MaxL+1])
            
            integer, dimension(Auto2e_MaxAngFuncSpher, 0:Auto2e_MaxL), parameter :: Orca_AngFuncMap &
                  = reshape([ORCA_S, ORCA_P, ORCA_D,  ORCA_F, ORCA_G, ORCA_H], [Auto2e_MaxAngFuncSpher,Auto2e_MaxL+1])

            integer, dimension(Auto2e_MaxAngFuncCart, 0:Auto2e_MaxL), parameter :: Molpro_CartAngFuncMap &
                  = reshape([MOLPRO_CART_S, MOLPRO_CART_P, MOLPRO_CART_D,  MOLPRO_CART_F, MOLPRO_CART_G, MOLPRO_CART_H], &
                  [Auto2e_MaxAngFuncCart,Auto2e_MaxL+1])

            integer, dimension(Auto2e_MaxAngFuncSpher, 0:Auto2e_MaxL), parameter :: Dalton_AngFuncMap &
                  = reshape([DALTON_S, DALTON_P, DALTON_D,  DALTON_F, DALTON_G, DALTON_H], [Auto2e_MaxAngFuncSpher,Auto2e_MaxL+1])

            integer, dimension(Auto2e_MaxAngFuncSpher, 0:Auto2e_MaxL) :: AngFuncMap
            integer, dimension(Auto2e_MaxAngFuncCart, 0:Auto2e_MaxL) :: CartAngFuncMap
            integer :: kk, L, n
            integer :: p0, p1, p
            integer :: NAO

            if (AOBasis%SpherAO) then
                  NAO = AOBasis%NAOSpher
            else
                  NAO = AOBasis%NAOCart
            end if
            Map_extao2ao = -1
            if (AOBasis%SpherAO) then
                  select case (ExternalOrdering)
                  case (ORBITAL_ORDERING_MOLPRO)
                        AngFuncMap = Molpro_AngFuncMap
                  case (ORBITAL_ORDERING_ORCA)
                        AngFuncMap = Orca_AngFuncMap
                  case (ORBITAL_ORDERING_DALTON)
                        AngFuncMap = Dalton_AngFuncMap
                  case default
                        call msg("Invalid external orbital ordering")
                        error stop
                  end select
                  associate ( &
                        ShellLoc => AOBasis%ShellLocSpher, &
                        ShellParamsIdx => AOBasis%ShellParamsIdx, &
                        ShellMomentum => AOBasis%ShellMomentum, &
                        NAngFunc => AOBasis%NAngFuncSpher, &
                        NShells => AOBasis%NShells, &
                        NAO => AOBasis%NAOSpher &
                        )
                        do kk = 1, NShells
                              k = ShellParamsIdx(kk)
                              L = ShellMomentum(k)
                              n = NAngFunc(k)
                              p0 = ShellLoc(kk)
                              p1 = ShellLoc(kk) + n - 1
                              do p = p0, p1
                                    Map_extao2ao(p) = p0 + AngFuncMap(p-p0+1, L) - 1
                              end do
                        end do
                  end associate
            else
                  select case (ExternalOrdering)
                  case (ORBITAL_ORDERING_MOLPRO)
                        CartAngFuncMap = Molpro_CartAngFuncMap
                  case default
                        call msg("Invalid external orbital ordering")
                        error stop
                  end select
                  associate ( &
                        ShellLoc => AOBasis%ShellLocCart, &
                        ShellParamsIdx => AOBasis%ShellParamsIdx, &
                        ShellMomentum => AOBasis%ShellMomentum, &
                        NAngFunc => AOBasis%NAngFuncCart, &
                        NShells => AOBasis%NShells, &
                        NAO => AOBasis%NAOCart &
                        )
                        do kk = 1, NShells
                              k = ShellParamsIdx(kk)
                              L = ShellMomentum(k)
                              n = NAngFunc(k)
                              p0 = ShellLoc(kk)
                              p1 = ShellLoc(kk) + n - 1
                              do p = p0, p1
                                    Map_extao2ao(p) = p0 + CartAngFuncMap(p-p0+1, L) - 1
                              end do
                        end do
                  end associate
            end if
            do p = 1, NAO
                  Map_ao2extao(Map_extao2ao(p)) = p
            end do
      end subroutine auto2e_interface_AOMap


      subroutine auto2e_interface_ApplyOrcaPhases_Vector(W, AOBasis)
            !
            ! Change phases of the real combinations of the spherical
            ! f, g, and h functions to match the Orca convention.
            ! The nonstandard phase for f and g functions, m=+-3 and +-4,
            ! is described in Orca's gto1.cpp file. We guessed that the
            ! same convention is applied for the h functions m=+-3 and
            ! +-4, but not for m=+-5.
            !
            real(F64), dimension(:), intent(inout) :: W
            type(TAOBasis), intent(in)             :: AOBasis
            
            integer :: k, kk
            integer :: L
            integer :: p0, p1

            do kk = 1, AOBasis%NShells
                  k = AOBasis%ShellParamsIdx(kk)
                  L = AOBasis%ShellMomentum(k)
                  select case (L)
                  case (3) ! f functions
                        p0 = AOBasis%ShellLocSpher(kk) + L + 3
                        p1 = AOBasis%ShellLocSpher(kk) + L - 3
                        W(p0) = -W(p0)
                        W(p1) = -W(p1)
                  case (4) ! g functions
                        p0 = AOBasis%ShellLocSpher(kk) + L + 3
                        p1 = AOBasis%ShellLocSpher(kk) + L - 3
                        W(p0) = -W(p0)
                        W(p1) = -W(p1)
                        p0 = AOBasis%ShellLocSpher(kk) + L + 4
                        p1 = AOBasis%ShellLocSpher(kk) + L - 4
                        W(p0) = -W(p0)
                        W(p1) = -W(p1)
                  case (5) ! h functions
                        p0 = AOBasis%ShellLocSpher(kk) + L + 3
                        p1 = AOBasis%ShellLocSpher(kk) + L - 3
                        W(p0) = -W(p0)
                        W(p1) = -W(p1)
                        p0 = AOBasis%ShellLocSpher(kk) + L + 4
                        p1 = AOBasis%ShellLocSpher(kk) + L - 4
                        W(p0) = -W(p0)
                        W(p1) = -W(p1)
                  end select
            end do
      end subroutine auto2e_interface_ApplyOrcaPhases_Vector


      subroutine auto2e_interface_ApplyOrcaPhases_Matrix(X, AOBasis, TwoIndexTransf)
            real(F64), dimension(:, :), intent(inout) :: X
            type(TAOBasis), intent(in)                :: AOBasis
            logical, intent(in)                       :: TwoIndexTransf

            integer :: NAO, k
            real(F64), dimension(:), allocatable :: W

            NAO = size(X, dim=1)
            do k = 1, NAO
                  call auto2e_interface_ApplyOrcaPhases_Vector(X(:, k), AOBasis)
            end do
            if (TwoIndexTransf) then
                  allocate(W(NAO))
                  do k = 1, NAO
                        W = X(k, :)
                        call auto2e_interface_ApplyOrcaPhases_Vector(W, AOBasis)
                        X(k, :) = W
                  end do
            end if
      end subroutine auto2e_interface_ApplyOrcaPhases_Matrix
end module Auto2eInterface
