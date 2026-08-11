module basis_sets
   use arithmetic
   use SpherGTO
   use spherh
   use specf
   use sort
   use math_constants
   use sys_definitions
   use periodic
   use io
   use string
   use Auto2e
   use grid_definitions, only: GRID_AOTHRESH
   use real_linalg
   use basis_definitions

   implicit none


   type TBasisConfig
      !
      ! Basis set configuration for an atom or an element
      ! ---
      !
      ! Atomic number
      !
      integer :: Z
      !
      ! Path to the text file with basis set parameters
      !
      character(:), allocatable :: PathToParams
      !
      ! Path to guess density matrix
      !
      character(:), allocatable :: PathToGuess
      logical :: GuessAvailable
      !
      ! Number of shell parameter sets per atom
      !
      integer :: NShellParams
      !
      ! Offset in global shell parameter arrays
      !
      integer :: Offset
   end type TBasisConfig

contains

   subroutine basis_CreateConfigs(Configs, AtomConfigMap, NConfigs, System, BasisAssign)
      type(TBasisConfig), allocatable, intent(out)    :: Configs(:)
      integer, dimension(:), allocatable, intent(out) :: AtomConfigMap
      integer, intent(out)                            :: NConfigs
      type(TSystem), intent(in)                       :: System
      type(TBasisAssignment), intent(in)              :: BasisAssign

      integer :: a, c, Z
      character(:), allocatable :: PathToParams
      character(:), allocatable :: PathToGuessDir
      character(:), allocatable :: PathToGuess
      logical :: GuessAvailable
      logical :: found
      type(TBasisRule) :: Rule

      allocate(Configs(System%NAtoms))
      allocate(AtomConfigMap(System%NAtoms))
      NConfigs = 0

      do a = 1, System%NAtoms
         Z = System%ZNumbers(a)
         PathToGuess = ""
         GuessAvailable = .false.

         Rule = BasisAssign%get_atom_rule(a, Z)
         PathToParams = Rule%PathToParams

         if (Rule%GuessAvailable) then
            PathToGuessDir = Rule%PathToGuessDir
            PathToGuess = PathToGuessDir // trim(lowercase(elname_short(Z))) // ".txt"
            if (io_exists(PathToGuess)) then
               GuessAvailable = .true.
            else
               PathToGuess = ""
            end if
         end if

         found = .false.
         do c = 1, NConfigs
            if (Configs(c)%Z == Z .and. Configs(c)%PathToParams == PathToParams) then
               AtomConfigMap(a) = c
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            NConfigs = NConfigs + 1
            Configs(NConfigs)%Z = Z
            Configs(NConfigs)%PathToParams = PathToParams
            Configs(NConfigs)%GuessAvailable = GuessAvailable
            if (GuessAvailable) then
               Configs(NConfigs)%PathToGuess = PathToGuess
            end if
            AtomConfigMap(a) = NConfigs
         end if
      end do
   end subroutine basis_CreateConfigs


   subroutine basis_AtomicRhoGuess(Rho_ao, AOBasis, System, SpherAO)
      !
      ! Generate guess density matrix as a superposition of atomic block-diagonal densities.
      ! Read the correct guess density for each atom from disk using the AOBasis%Assignment rules.
      !
      ! Fused basis sets made by composition of two basis sets with basis_FuseBasisSets are not supported.
      !
      ! The guess density is generated only for non-ghost (real) atoms.
      !
      real(F64), dimension(:, :), intent(out) :: Rho_ao
      type(TAOBasis), intent(in)              :: AOBasis
      type(TSystem), intent(in)               :: System
      logical, intent(in)                     :: SpherAO

      type(TBasisConfig), allocatable :: Configs(:)
      integer, dimension(:), allocatable :: AtomConfigMap
      integer :: NConfigs
      integer :: a, c, s
      integer :: i0, i1, p0, p1
      integer :: sh1, sh2
      logical :: rholoaded

      if (AOBasis%Fused) then
         call msg("basis_AtomicRhoGuess: Fused basis sets made by composition of two basis sets with basis_FuseBasisSets are not supported.", MSG_ERROR)
         stop
      end if

      Rho_ao = ZERO

      call basis_CreateConfigs(Configs, AtomConfigMap, NConfigs, System, AOBasis%Assignment)

      associate(AtomShellMap => AOBasis%AtomShellMap, &
         ShellLocSpher => AOBasis%ShellLocSpher, &
         ShellLocCart => AOBasis%ShellLocCart, &
         ShellParamsIdx => AOBasis%ShellParamsIdx)
         do c = 1, NConfigs
            if (Configs(c)%GuessAvailable) then
               rholoaded = .false.
               do s = 1, 2
                  do a = System%RealAtoms(1, s), System%RealAtoms(2, s)
                     if (AtomConfigMap(a) == c) then
                        !
                        ! Determine the basis function boundaries for atom
                        ! Since fused basis sets are not supported, there is only 1 segment
                        !
                        sh1 = AtomShellMap(1, 1, a)
                        sh2 = AtomShellMap(2, 1, a)

                        if (SpherAO) then
                           i0 = ShellLocSpher(sh1)
                           i1 = ShellLocSpher(sh2) + AOBasis%NAngFuncSpher(ShellParamsIdx(sh2)) - 1
                        else
                           i0 = ShellLocCart(sh1)
                           i1 = ShellLocCart(sh2) + AOBasis%NAngFuncCart(ShellParamsIdx(sh2)) - 1
                        end if

                        if (.not. rholoaded) then
                           call io_text_read(Rho_ao(i0:i1, i0:i1), Configs(c)%PathToGuess)
                           rholoaded = .true.
                           p0 = i0
                           p1 = i1
                        else
                           Rho_ao(i0:i1, i0:i1) = Rho_ao(p0:p1, p0:p1)
                        end if
                     end if
                  end do
               end do
            else
               call msg("Note: Incomplete SCF guess. Atomic density missing for " // trim(elname_short(Configs(c)%Z)))
            end if
         end do
      end associate
   end subroutine basis_AtomicRhoGuess


   subroutine basis_Init(AOBasis, System, FilePath, SpherAO, ShellOrder, BasisAssign)
      !
      ! Create a new instance of a user-defined type which encapsulates all necessary
      ! basis-set data used to calculate integrals in the selected Gaussian-type
      ! atomic orbital basis set. The coefficients are loaded from a text file
      ! in the EMSL database format for the GAMESS program.
      !
      ! The data includes contraction coefficients, exponents, number of primitives,
      ! angular momenta, centers on which the functions are located, location of a particular
      ! orbital shell within a vector of AO indices, etc.
      !
      ! The same instance of TAOBasis can be used for subroutines which work in Cartesian
      ! and spherical Gaussian basis sets. That is because the relevant quantities
      ! are computed for both spherical harmonic and Cartesian basis sets, regardless
      ! of the value of SpherAO. The SpherAO parameter can be safely changed
      ! without invoking this subroutine again.
      !
      ! Arguments:
      !   AOBasis              Output TAOBasis object to be initialized.
      !
      !   System               TSystem object containing atomic coordinates and numbers.
      !
      !   FilePath             (Optional) Path to EMSL-format basis set parameters file.
      !                        Use this for a global basis set applied to all atoms.
      !                        For per-atom assignments, use BasisAssign instead.
      !                        Either FilePath or BasisAssign must be provided.
      !   SpherAO              (Optional) Logical flag indicating if spherical harmonics
      !                        are used (default: .true.).
      !
      !   ShellOrder           (Optional) Mode of shell ordering (see basis_definitions).
      !                        Defaults to SHELL_ORDER_BY_RADIUS.
      !
      !   BasisAssign          (Optional) TBasisAssignment object that provides a detailed
      !                        basis set-to-atom map. Either FilePath or BasisAssign
      !                        must be provided, but not both.
      !
      type(TAOBasis), intent(out)                  :: AOBasis
      type(TSystem), intent(in)                    :: System
      character(*), intent(in), optional           :: FilePath
      logical, intent(in), optional                :: SpherAO
      integer, optional, intent(in)                :: ShellOrder
      type(TBasisAssignment), intent(in), optional :: BasisAssign

      type(TBasisConfig), allocatable :: Configs(:)
      integer, dimension(:), allocatable :: AtomConfigMap
      integer :: a, NConfigs, NShells

      integer :: L, k, p, p0, p1, q, q0, q1
      integer :: NShellParamsTotal, MaxNPrimitives, LmaxGTO
      integer :: MaxNPrimitives_k, LmaxGTO_k
      integer :: MaxNAngFuncCart
      integer :: n
      integer, dimension(:), allocatable :: ShellMomentum
      integer, dimension(:), allocatable :: NPrimitives
      integer, dimension(:), allocatable :: ShellParamsIdx, ShellCenters
      integer, dimension(:, :), allocatable :: ConfigShellsMap
      integer, dimension(:), allocatable :: S
      real(F64), dimension(:), allocatable :: W, R2Max
      real(F64), dimension(:, :), allocatable :: CntrCoeffs, Exponents, NormFactorsCart
      integer :: Order
      logical :: Spherical
      type(TBasisRule)       :: ResolvedRule
      type(TBasisAssignment) :: Assignment

      if (present(FilePath) .and. present(BasisAssign)) then
         call msg("basis_Init: FilePath and BasisAssign cannot be provided simultaneously.", MSG_ERROR)
         error stop
      end if
      if (.not. present(FilePath) .and. .not. present(BasisAssign)) then
         call msg("basis_Init: At least one of FilePath or BasisAssign must be provided", MSG_ERROR)
         error stop
      end if

      if (present(BasisAssign)) then
         Assignment = BasisAssign
      else if (present(FilePath)) then
         call basis_ResolvePath(ResolvedRule, Assignment, "FILE " // FilePath)
         ResolvedRule%id = 0
         call Assignment%add_global_fallback(ResolvedRule)
      end if

      if (present(SpherAO)) then
         Spherical = SpherAO
      else
         Spherical = .true.
      end if

      if (present(ShellOrder)) then
         Order = ShellOrder
      else
         Order = SHELL_ORDER_BY_RADIUS
      end if
      !
      ! 1. Group Unique Blocks by (Z, PathToParams) & Determine Paths
      !
      call basis_CreateConfigs(Configs, AtomConfigMap, NConfigs, &
         System, Assignment)
      !
      ! 3. Query properties
      !
      MaxNPrimitives = -1
      LmaxGTO = -1
      do k = 1, NConfigs
         call basis_Query(Configs(k)%NShellParams, MaxNPrimitives_k, LmaxGTO_k, &
            Configs(k)%Z, Configs(k)%PathToParams)
         MaxNPrimitives = max(MaxNPrimitives, MaxNPrimitives_k)
         LmaxGTO = max(LmaxGTO, LmaxGTO_k)
      end do
      NShellParamsTotal = 0
      do k = 1, NConfigs
         NShellParamsTotal = NShellParamsTotal + Configs(k)%NShellParams
      end do
      allocate(NPrimitives(NShellParamsTotal))
      allocate(ShellMomentum(NShellParamsTotal))
      allocate(CntrCoeffs(MaxNPrimitives, NShellParamsTotal))
      allocate(Exponents(MaxNPrimitives, NShellParamsTotal))
      allocate(ConfigShellsMap(2, NConfigs))

      p0 = 1
      p1 = 1
      do k = 1, NConfigs
         p1 = p0 + Configs(k)%NShellParams - 1
         call basis_ReadElementData(CntrCoeffs(:, p0:p1), Exponents(:, p0:p1), &
            ShellMomentum(p0:p1), NPrimitives(p0:p1), Configs(k)%Z, Configs(k)%PathToParams)
         ConfigShellsMap(1, k) = p0
         ConfigShellsMap(2, k) = p1
         Configs(k)%Offset = p0 - 1
         p0 = p0 + Configs(k)%NShellParams
      end do
      MaxNAngFuncCart = ((LmaxGTO + 1) * (LmaxGTO + 2)) / 2
      allocate(NormFactorsCart(MaxNAngFuncCart, NShellParamsTotal))
      call basis_NormalizeCntrCoeffs(NormFactorsCart, CntrCoeffs, &
         Exponents, NPrimitives, ShellMomentum)
      select case(Order)
       case(SHELL_ORDER_BY_RADIUS)
         !
         ! Compute the radii (R2MAX) beyond which the absolute values of atomic orbitals
         ! fall below some small value, e.g., eps~10**(-12). The shells within each atom
         ! are then sorted according to decreasing effective radii. Sorting is done
         ! every time a basis set is loaded, but it is actually used only for orbital
         ! evalation on the numerical grid to exit the loop over quadrature points
         ! as soon as possible. It will have no effect on the code
         ! that does not employ numerical integration on the molecular grid.
         !
         ! Note that shell sorting makes the shell order different than that defined
         ! in the basis set coefficients file.
         !
         allocate(R2Max(NShellParamsTotal))
         call basis_R2Max(R2Max, GRID_AOTHRESH, ShellMomentum, NPrimitives, &
            CntrCoeffs, Exponents, NormFactorsCart, NShellParamsTotal)
         !
         ! Sort shell indices for every atom according to R2MAX,
         ! in decreasing order. Sorting is performed only within
         ! atomic orbitals belonging to one atom at a time
         !
         allocate(W(NShellParamsTotal))
         allocate(S(NShellParamsTotal))
         do p = 1, NShellParamsTotal
            S(p) = p
         end do
         do k = 1, NConfigs
            p0 = ConfigShellsMap(1, k)
            p1 = ConfigShellsMap(2, k)
            n = p1 - p0 + 1
            W(1:n) = -R2Max(p0:p1)
            call dsort(W(1:n), S(p0:p1), n)
         end do
       case(SHELL_ORDER_BY_MOMENTUM)
         !
         ! Sort shells within each atom according to increasing angular momentum.
         ! This will disable the computation of orbitals spatial extent (R2Max)
         ! and grid screening in numerical integrals on the molecular grid.
         ! The order of different shells of the same angular momentum will be
         ! the same as in the text file with basis set parameters.
         !
         allocate(R2Max(NShellParamsTotal))
         R2Max = huge(ONE)
         allocate(S(NShellParamsTotal))
         S = -1
         do k = 1, NConfigs
            p0 = ConfigShellsMap(1, k)
            p1 = ConfigShellsMap(2, k)
            q = 0
            do L = 0, LmaxGTO
               do p = p0, p1
                  if (ShellMomentum(p) == L) then
                     S(p0+q) = p
                     q = q + 1
                  end if
               end do
            end do
         end do
       case(SHELL_ORDER_FIXED)
         allocate(R2Max(NShellParamsTotal))
         R2Max = huge(ONE)
         allocate(S(NShellParamsTotal))
         do p = 1, NShellParamsTotal
            S(p) = p
         end do
       case default
         call msg("basis_Init: Invalid ShellOrder value.", MSG_ERROR)
         error stop
      end select
      !
      ! Assign shells to each atom in the system
      !
      NShells = 0
      do a = 1, System%NAtoms
         k = AtomConfigMap(a)
         p0 = ConfigShellsMap(1, k)
         p1 = ConfigShellsMap(2, k)
         NShells = NShells + p1 - p0 + 1
      end do
      allocate(ShellParamsIdx(NShells))
      allocate(ShellCenters(NShells))
      q = 0
      do a = 1, System%NAtoms
         k = AtomConfigMap(a)
         p0 = ConfigShellsMap(1, k)
         p1 = ConfigShellsMap(2, k)
         n = p1 - p0 + 1
         q0 = q + 1
         q1 = q + n
         ShellParamsIdx(q0:q1) = S(p0:p1)
         ShellCenters(q0:q1) = a
         q = q + n
      end do
      call basis_Init_2(AOBasis, System%AtomCoords, ShellCenters, ShellParamsIdx, ShellMomentum, &
         NPrimitives, CntrCoeffs, Exponents, NormFactorsCart, R2Max, Spherical)
      AOBasis%Assignment = Assignment
   end subroutine basis_Init


   function basis_Ang2Int(Momentum)
      integer                  :: basis_Ang2Int
      character(1), intent(in) :: Momentum

      select case (uppercase(Momentum))
       case ("S")
         basis_Ang2Int = 0
       case ("P")
         basis_Ang2Int = 1
       case ("D")
         basis_Ang2Int = 2
       case ("F")
         basis_Ang2Int = 3
       case ("G")
         basis_Ang2Int = 4
       case ("H")
         basis_Ang2Int = 5
       case ("I")
         basis_Ang2Int = 6
       case ("K")
         basis_Ang2Int = 7
       case ("L") ! S+P shell
         basis_Ang2Int = -1
       case default
         basis_Ang2Int = -2
      end select
   end function basis_Ang2Int


   subroutine basis_Query(NShellParams, MaxNPrimitives, LmaxGTO, Z, FilePath)
      integer, intent(out)     :: NShellParams
      integer, intent(out)     :: MaxNPrimitives
      integer, intent(out)     :: LmaxGTO
      integer, intent(in)      :: Z
      character(*), intent(in) :: FilePath

      integer :: u
      logical :: eof
      logical :: found_element, end_of_data
      logical :: blank, comment, endkey
      character(1) :: Momentum
      integer :: NPrimitives, L, k
      character(:), allocatable :: ElementLabel
      character(:), allocatable :: raw_line, line

      if (Z < 1 .or. Z > KNOWN_ELEMENTS) then
         call msg("Elements list includes invalid element", MSG_ERROR)
         error stop
      end if
      u = io_text_open(FilePath, "OLD")
      ElementLabel = uppercase(trim(ELNAME_LONG(Z)))
      eof = .false.
      found_element = .false.
      do while (.not. (eof .or. found_element))
         call io_text_readline(raw_line, u, eof)
         if (.not. eof) then
            line = uppercase(trim(raw_line))
            if (isblank(line) .or. iscomment(line)) then
               cycle
            else if (line=="$END") then
               eof = .true.
            else
               if (line==ElementLabel) then
                  found_element = .true.
               end if
            end if
         end if
      end do
      NShellParams = 0
      LmaxGTO = -1
      MaxNPrimitives = -1
      if (found_element .and. .not. eof) then
         end_of_data = .false.
         do while (.not. end_of_data)
            call io_text_readline(raw_line, u, eof)
            line = uppercase(trim(raw_line))
            blank = isblank(line)
            comment = iscomment(line)
            endkey = (line=="$END")
            if (eof .or. blank .or. comment .or. endkey) end_of_data = .true.
            if (.not. end_of_data) then
               read(line, *) Momentum, NPrimitives
               L = basis_Ang2Int(Momentum)
               if (L == -2) then
                  call msg("Unknown angular function symbol", MSG_ERROR)
                  error stop
               end if
               if (L == -1) then
                  !
                  ! S+P shells
                  !
                  NShellParams = NShellParams + 2
                  L = 1
               else
                  NShellParams = NShellParams + 1
               end if
               if (L > AUTO2E_MAXL) then
                  call msg("Angular momentum exceeds the upper limit of the Auto2e module", MSG_ERROR)
                  error stop
               end if
               LmaxGTO = max(LmaxGTO, L)
               MaxNPrimitives = max(MaxNPrimitives, NPrimitives)
               do k = 1, NPrimitives
                  call io_text_readline(raw_line, u, eof)
                  line = uppercase(trim(raw_line))
                  blank = isblank(line)
                  comment = iscomment(line)
                  endkey = (line=="$END")
                  if (eof .or. blank .or. comment .or. endkey) then
                     call msg("Unexpected end of basis set parameters", MSG_ERROR)
                     error stop
                  end if
               end do
            end if
         end do
      end if
      if (NShellParams == 0 .or. LmaxGTO < 0 .or. MaxNPrimitives <= 0) then
         call msg("Invalid basis set parameters for " // ElementLabel, MSG_ERROR)
         error stop
      end if
      close(u)
   end subroutine basis_Query


   subroutine basis_ReadElementData(CntrCoeffs, Exponents, ShellMomentum, NPrimitives, Z, FilePath)
      real(F64), dimension(:, :), intent(out) :: CntrCoeffs
      real(F64), dimension(:, :), intent(out) :: Exponents
      integer, dimension(:), intent(out)      :: ShellMomentum
      integer, dimension(:), intent(out)      :: NPrimitives
      integer, intent(in)                     :: Z
      character(*), intent(in)                :: FilePath

      integer :: u
      logical :: eof
      logical :: found_element, end_of_data
      logical :: blank, comment, endkey
      character(1) :: Momentum
      integer :: NP, L, k, v
      integer :: NShellParams
      integer :: LmaxGTO, MaxNPrimitives
      real(F64) :: Ev, Cv, C1v, C2v
      character(:), allocatable :: ElementLabel
      character(:), allocatable :: raw_line, line

      if (Z < 1 .or. Z > KNOWN_ELEMENTS) then
         call msg("Elements list includes invalid element", MSG_ERROR)
         error stop
      end if
      u = io_text_open(FilePath, "OLD")
      ElementLabel = uppercase(trim(ELNAME_LONG(Z)))
      eof = .false.
      found_element = .false.
      do while (.not. (eof .or. found_element))
         call io_text_readline(raw_line, u, eof)
         if (.not. eof) then
            line = uppercase(trim(raw_line))
            if (isblank(line) .or. iscomment(line)) then
               cycle
            else if (line=="$END") then
               eof = .true.
            else
               if (line==ElementLabel) then
                  found_element = .true.
               end if
            end if
         end if
      end do
      NShellParams = 0
      LmaxGTO = -1
      MaxNPrimitives = -1
      if (found_element .and. .not. eof) then
         end_of_data = .false.
         do while (.not. end_of_data)
            call io_text_readline(raw_line, u, eof)
            line = uppercase(trim(raw_line))
            blank = isblank(line)
            comment = iscomment(line)
            endkey = (line=="$END")
            if (eof .or. blank .or. comment .or. endkey) end_of_data = .true.
            if (.not. end_of_data) then
               read(line, *) Momentum, NP
               L = basis_Ang2Int(Momentum)
               if (L == -2) then
                  call msg("Unknown angular function symbol", MSG_ERROR)
                  error stop
               else if (L == -1) then
                  !
                  ! S+P shell
                  !
                  do k = 1, NP
                     call io_text_readline(raw_line, u, eof)
                     line = uppercase(trim(raw_line))
                     blank = isblank(line)
                     comment = iscomment(line)
                     endkey = (line=="$END")
                     if (eof .or. blank .or. comment .or. endkey) then
                        call msg("Unexpected end of basis set parameters", MSG_ERROR)
                        error stop
                     end if
                     read(line, *) v, Ev, C1v, C2v
                     CntrCoeffs(v, NShellParams+1) = C1v
                     Exponents(v, NShellParams+1) = Ev
                     CntrCoeffs(v, NShellParams+2) = C2v
                     Exponents(v, NShellParams+2) = Ev
                  end do
                  NPrimitives(NShellParams+1) = NP
                  ShellMomentum(NShellParams+1) = 0
                  NPrimitives(NShellParams+2) = NP
                  ShellMomentum(NShellParams+2) = 1
                  NShellParams = NShellParams + 2
                  L = 1
               else
                  do k = 1, NP
                     call io_text_readline(raw_line, u, eof)
                     line = uppercase(trim(raw_line))
                     blank = isblank(line)
                     comment = iscomment(line)
                     endkey = (line=="$END")
                     if (eof .or. blank .or. comment .or. endkey) then
                        call msg("Unexpected end of basis set parameters", MSG_ERROR)
                        error stop
                     end if
                     read(line, *) v, Ev, Cv
                     CntrCoeffs(v, NShellParams+1) = Cv
                     Exponents(v, NShellParams+1) = Ev
                  end do
                  NPrimitives(NShellParams+1) = NP
                  ShellMomentum(NShellParams+1) = L
                  NShellParams = NShellParams + 1
               end if
               LmaxGTO = max(L, LmaxGTO)
               MaxNPrimitives = max(NP, MaxNPrimitives)
            end if
         end do
      end if
      if (NShellParams == 0 .or. LmaxGTO < 0 .or. MaxNPrimitives <= 0) then
         call msg("Invalid basis set parameters for " // ElementLabel, MSG_ERROR)
         error stop
      end if
      close(u)
   end subroutine basis_ReadElementData


   subroutine basis_Init_2(AOBasis, AtomCoords, ShellCenters, ShellParamsIdx, ShellMomentum, &
      NPrimitives, CntrCoeffs, Exponents, NormFactorsCart, R2Max, SpherAO)

      type(TAOBasis), intent(out)                :: AOBasis
      real(F64), dimension(:, :), intent(in)     :: AtomCoords
      integer, dimension(:), intent(in)          :: ShellCenters
      integer, dimension(:), intent(in)          :: ShellParamsIdx
      integer, dimension(:), intent(in)          :: ShellMomentum
      integer, dimension(:), intent(in)          :: NPrimitives
      real(F64), dimension(:, :), intent(in)     :: CntrCoeffs
      real(F64), dimension(:, :), intent(in)     :: Exponents
      real(F64), dimension(:, :), intent(in)     :: NormFactorsCart
      real(F64), dimension(:), intent(in)        :: R2Max
      logical, intent(in)                        :: SpherAO

      integer :: MaxNAngFuncSpher, MaxNAngFuncCart
      integer :: L, lx, ly, i, j, a, b, j0
      integer :: k
      integer :: NShellsI

      associate (NShellParams=>AOBasis%NShellParams, NShells=>AOBasis%NShells, &
         LmaxGTO=>AOBasis%LmaxGTO, NAOSpher=>AOBasis%NAOSpher, NAOCart=>AOBasis%NAOCart, &
         MaxNPrimitives=>AOBasis%MaxNPrimitives, NAtoms=>AOBasis%NAtoms, &
         MaxNShells=>AOBasis%MaxNShells)

         NShellParams = size(ShellMomentum)
         NShells = size(ShellCenters)
         NAtoms = size(AtomCoords, dim=2)
         LmaxGTO = maxval(ShellMomentum(1:NShellParams))
         MaxNPrimitives = maxval(NPrimitives(1:NShellParams))
         MaxNAngFuncSpher = 2 * LmaxGTO + 1
         MaxNAngFuncCart = ((LmaxGTO + 1) * (LmaxGTO + 2)) / 2
         allocate(AOBasis%AtomCoords(3, NAtoms))
         allocate(AOBasis%ShellCenters(NShells))
         allocate(AOBasis%ShellParamsIdx(NShells))
         allocate(AOBasis%ShellMomentum(NShellParams))
         allocate(AOBasis%NPrimitives(NShellParams))
         allocate(AOBasis%CntrCoeffs(MaxNPrimitives, NShellParams))
         allocate(AOBasis%Exponents(MaxNPrimitives, NShellParams))
         allocate(AOBasis%NormFactorsCart(MaxNAngFuncCart, NShellParams))
         allocate(AOBasis%NormFactorsSpher(MaxNAngFuncSpher, NShellParams))
         allocate(AOBasis%NAngFuncCart(NShellParams))
         allocate(AOBasis%NAngFuncSpher(NShellParams))
         allocate(AOBasis%ShellLocSpher(NShells))
         allocate(AOBasis%ShellLocCart(NShells))
         allocate(AOBasis%R2Max(NShellParams))
         AOBasis%R2Max = R2Max
         AOBasis%AtomCoords = AtomCoords
         AOBasis%ShellCenters = ShellCenters
         allocate(AOBasis%AtomShellN(NAtoms))
         AOBasis%AtomShellN = 0
         do i = 1, NAtoms
            do j = 1, NShells
               if (j > 1) then
                  a = ShellCenters(j - 1)
                  b = ShellCenters(j)
                  if (b == i .and. a /= i) then
                     AOBasis%AtomShellN(i) = AOBasis%AtomShellN(i) + 1
                  end if
               else
                  b = ShellCenters(j)
                  if (b == i) then
                     AOBasis%AtomShellN(i) = 1
                  end if
               end if
            end do
         end do
         allocate(AOBasis%AtomShellMap(2, maxval(AOBasis%AtomShellN), NAtoms))
         AOBasis%AtomShellMap = -1
         do i = 1, NAtoms
            !
            ! Loop over shell segments, i.e., contiguous ranges of shells belonging
            ! to the same atom
            !
            kloop: do k = 1, AOBasis%AtomShellN(i)
               if (k == 1) then
                  j0 = 1
               else
                  j0 = AOBasis%AtomShellMap(2, k-1, i) + 1
               end if
               do j = j0, NShells
                  if (ShellCenters(j) == i) then
                     if (AOBasis%AtomShellMap(1, k, i) == -1) then
                        !
                        ! First shell in the current segment
                        !
                        AOBasis%AtomShellMap(1, k, i) = j
                        AOBasis%AtomShellMap(2, k, i) = j
                     else
                        !
                        ! Next shell in the current segment
                        !
                        AOBasis%AtomShellMap(2, k, i) = j
                     end if
                  else
                     if (AOBasis%AtomShellMap(2, k, i) /= -1) then
                        !
                        ! End of segment
                        !
                        cycle kloop
                     end if
                  end if
               end do
            end do kloop
         end do
         MaxNShells = 0
         do i = 1, NAtoms
            NShellsI = 0
            do k = 1, AOBasis%AtomShellN(i)
               a = AOBasis%AtomShellMap(1, k, i)
               b = AOBasis%AtomShellMap(2, k, i)
               NShellsI = NShellsI + b - a + 1
            end do
            MaxNShells = max(MaxNShells, NShellsI)
         end do
         AOBasis%ShellParamsIdx = ShellParamsIdx
         AOBasis%ShellMomentum = ShellMomentum
         AOBasis%NPrimitives = NPrimitives
         AOBasis%CntrCoeffs = CntrCoeffs(1:MaxNPrimitives, :)
         AOBasis%Exponents = Exponents(1:MaxNPrimitives, :)
         AOBasis%NormFactorsCart = NormFactorsCart(1:MaxNAngFuncCart, :)
         AOBasis%SpherAO = SpherAO
         call basis_CountOrbitals(NAOSpher, NAOCart, &
            AOBasis%NAngFuncSpher, &
            AOBasis%NAngFuncCart, &
            AOBasis%ShellLocSpher, &
            AOBasis%ShellLocCart, &
            ShellMomentum, ShellParamsIdx, NShellParams, NShells)
         call SpherGTO_NormFactors(AOBasis%NormFactorsSpher, &
            ShellMomentum, NPrimitives, CntrCoeffs, Exponents, NShellParams)
         !
         ! Exponents defining the sequence of Cartesian angular functions (x**l)*(y**m)*(z**n).
         ! The sequence of Cartesian polynomials has to be consistent with the angular functions
         ! defined in the automatic two-electron integrals code (Auto2e).
         !
         allocate(AOBasis%CartPolyX(MaxNAngFuncCart, 0:LmaxGTO))
         allocate(AOBasis%CartPolyY(MaxNAngFuncCart, 0:LmaxGTO))
         allocate(AOBasis%CartPolyZ(MaxNAngFuncCart, 0:LmaxGTO))
         do L = 0, LmaxGTO
            i = 1
            do lx = L, 0, -1
               do ly = L - lx, 0, -1
                  AOBasis%CartPolyX(i, L) = lx
                  AOBasis%CartPolyY(i, L) = ly
                  AOBasis%CartPolyZ(i, L) = L - lx - ly
                  i = i + 1
               end do
            end do
         end do
         !
         ! Maximum orbital angular momentum for each atom
         !
         allocate(AOBasis%MaxAtomL(NAtoms))
         AOBasis%MaxAtomL = 0
         do i = 1, NAtoms
            do k = 1, AOBasis%AtomShellN(i)
               a = AOBasis%AtomShellMap(1, k, i)
               b = AOBasis%AtomShellMap(2, k, i)
               do j = a, b
                  AOBasis%MaxAtomL(i) = max(AOBasis%MaxAtomL(i), ShellMomentum(ShellParamsIdx(j)))
               end do
            end do
         end do
      end associate
   end subroutine basis_Init_2


   subroutine basis_CountOrbitals(NAOSpher, NAOCart, NAngFuncSpher, NAngFuncCart, &
      ShellLocSpher, ShellLocCart, ShellMomentum, ShellParamsIdx, NShellParams, NShells)

      integer, intent(out)                          :: NAOSpher, NAOCart
      integer, dimension(NShellParams), intent(out) :: NAngFuncSpher, NAngFuncCart
      integer, dimension(NShells), intent(out)      :: ShellLocSpher, ShellLocCart
      integer, dimension(NShellParams), intent(in)  :: ShellMomentum
      integer, dimension(NShells)                   :: ShellParamsIdx
      integer, intent(in)                           :: NShellParams
      integer, intent(in)                           :: NShells

      integer :: s, a
      integer :: L
      integer :: ShellParamsA

      do s = 1, NShellParams
         L = ShellMomentum(s)
         NAngFuncSpher(s) = 2  * L + 1
         NAngFuncCart(s) = ((L + 1) * (L + 2)) / 2
      end do
      NAOSpher = 0
      NAOCart = 0
      do a = 1, NShells
         ShellParamsA = ShellParamsIdx(a)
         ShellLocSpher(a) = NAOSpher + 1
         ShellLocCart(a) = NAOCart + 1
         NAOSpher = NAOSpher + NAngFuncSpher(ShellParamsA)
         NAOCart = NAOCart + NAngFuncCart(ShellParamsA)
      end do
   end subroutine basis_CountOrbitals


   subroutine basis_CntrCoeffs_DataBase_Format(CntrCoeffsEMSL, AOBasis)
      !
      ! Translate contraction coefficients from the internal storage format
      ! to the EMSL database format.
      !
      ! In the EMSL format it is assumed that (i) the contraction coefficients
      ! multiply normalized gaussian primitives, (ii) the contracted Cartesian
      ! Gaussian corresponding to x**l is normalized
      !
      real(F64), dimension(:, :), intent(out) :: CntrCoeffsEMSL
      type(TAOBasis), intent(in) :: AOBasis

      integer :: u, v, L
      real(F64) :: N, Alpha

      associate (NShellParams => AOBasis%NShellParams, &
         NPrimitives => AOBasis%NPrimitives, &
         NormFactors => AOBasis%NormFactorsCart, &
         ShellMomentum => AOBasis%ShellMomentum, &
         Exponents => AOBasis%Exponents, &
         CntrCoeffs => AOBasis%CntrCoeffs &
         )
         CntrCoeffsEMSL = ZERO
         do u = 1, NShellParams
            L = ShellMomentum(u)
            do v = 1, NPrimitives(u)
               Alpha = Exponents(v, u)
               N = TWO**(L+THREE/FOUR) * Alpha**((TWO*L+THREE)/FOUR) &
                  / (PI**(THREE/FOUR) * dblfactorial(L))
               CntrCoeffsEMSL(v, u) = NormFactors(1, u) * CntrCoeffs(v, u) / N
            end do
         end do
      end associate
   end subroutine basis_CntrCoeffs_DataBase_Format


   subroutine basis_NormalizeCntrCoeffs(NormFactorsCart, CntrCoeffs, &
      Exponents, NPrimitives, ShellMomentum)
      !
      ! Rescale contraction coefficients and compute normalization factors
      ! for Cartesian Gaussians.
      !
      real(F64), dimension(:, :), intent(out)   :: NormFactorsCart
      real(F64), dimension(:, :), intent(inout) :: CntrCoeffs
      real(F64), dimension(:, :), intent(in)    :: Exponents
      integer, dimension(:), intent(in)         :: NPrimitives
      integer, dimension(:), intent(in)         :: ShellMomentum

      integer :: NShellParams
      integer :: a, NP, i, j, v
      real(F64) :: AlphaI, AlphaJ, CI, CJ
      real(F64) :: x, y, z
      integer :: lx, ly, lz
      integer :: L

      NShellParams = size(CntrCoeffs, dim=2)
      NormFactorsCart = ZERO
      do a = 1, NShellParams
         NP = NPrimitives(a)
         L = ShellMomentum(a)
         do i = 1, NP
            AlphaI = Exponents(i, a)
            CntrCoeffs(i, a) = CntrCoeffs(i, a) * (TWO/PI)**(THREE/FOUR)*TWO**L*(AlphaI)**((TWO*L+THREE)/FOUR)
         end do
         z = ZERO
         do i = 1, NP
            do j = 1, NP
               AlphaI = Exponents(i, a)
               AlphaJ = Exponents(j, a)
               CI = CntrCoeffs(i, a)
               CJ = CntrCoeffs(j, a)
               z = z + CI * CJ / (AlphaI + AlphaJ)**(L+THREE/TWO)
            end do
         end do
         y = Sqrt(PI**(THREE/TWO) / TWO**L)
         z = Sqrt(z)
         v = 1
         do lx = L, 0, -1
            do ly = L - lx, 0, -1
               lz = L - lx - ly
               x = Sqrt(dblfact(2*lx-1) * dblfact(2*ly-1) * dblfact(2*lz-1))
               NormFactorsCart(v, a) = ONE / (x * y * z)
               v = v + 1
            end do
         end do
      end do
   end subroutine basis_NormalizeCntrCoeffs


   subroutine basis_R2Max(R2Max, AOThresh, ShellMomentum, NPrimitives, &
      CntrCoeffs, Exponents, NormFactorsCart, NShellParams)

      real(F64), dimension(:), intent(out)   :: R2Max
      real(F64), intent(in)                  :: AOThresh
      integer, dimension(:), intent(in)      :: ShellMomentum
      integer, dimension(:), intent(in)      :: NPrimitives
      real(F64), dimension(:, :), intent(in) :: CntrCoeffs
      real(F64), dimension(:, :), intent(in) :: Exponents
      real(F64), dimension(:, :), intent(in) :: NormFactorsCart
      integer, intent(in)                    :: NShellParams

      real(F64) :: Ra
      integer :: a

      do a = 1, NShellParams
         call find_radius(Ra, AOThresh, ShellMomentum(a), NPrimitives(a), &
            CntrCoeffs(1:NPrimitives(a), a), Exponents(1:NPrimitives(a), a), &
            NormFactorsCart(1, a))
         R2Max(a) = Ra**2
      end do

   contains

      function phi_radial(r, L, NPrimitives, CntrCoeffs, Exponents, NormFactor)
         real(F64) :: phi_radial
         real(F64), intent(in)               :: r
         integer, intent(in)                 :: L
         integer, intent(in)                 :: NPrimitives
         real(F64), dimension(:), intent(in) :: CntrCoeffs
         real(F64), dimension(:), intent(in) :: Exponents
         real(F64), intent(in)               :: NormFactor

         real(F64) :: r2, rL
         integer :: i

         r2 = r**2
         rL = r**L
         phi_radial = ZERO
         do i = 1, NPrimitives
            phi_radial = phi_radial + CntrCoeffs(i) * rL * exp(-Exponents(i)*r2)
         end do
         phi_radial = NormFactor * phi_radial
      end function phi_radial


      subroutine find_radius(radius, AOThresh, L, NPrimitives, CntrCoeffs, Exponents, NormFactor)
         !
         ! Solve Phi(r) = AOThresh
         !
         real(F64)                           :: radius
         real(F64), intent(in)               :: AOThresh
         integer, intent(in)                 :: L
         integer, intent(in)                 :: NPrimitives
         real(F64), dimension(:), intent(in) :: CntrCoeffs
         real(F64), dimension(:), intent(in) :: Exponents
         real(F64), intent(in)               :: NormFactor

         integer :: i
         real(F64), dimension(:), allocatable :: AbsCntrCoeffs
         real(F64) :: AlphaMin
         real(F64) :: x, y, xleft, xright, deltax
         real(F64), parameter :: XThresh = 1.0E-12_F64
         logical :: converged
         integer, parameter :: MaxIters = 128

         allocate(AbsCntrCoeffs(NPrimitives))
         AbsCntrCoeffs = abs(CntrCoeffs(1:NPrimitives))
         !
         ! Look for x for which |Sum(i) a_i x**L exp(-alpha_i x**2)| is guaranteed
         ! to monotonically decrease
         ! ---
         ! Compute left bracketing value
         !
         AlphaMin = minval(Exponents(1:NPrimitives))
         if (L == 0) then
            xleft = ZERO
         else
            xleft = Sqrt(real(L, F64)/(TWO*AlphaMin))
         end if
         !
         ! Compute right bracketing value
         !
         deltax = ONE/AlphaMin
         converged = .false.
         iters1: do i = 1, MaxIters
            xright = xleft + deltax
            y = phi_radial(xright, L, NPrimitives, AbsCntrCoeffs, Exponents, NormFactor)
            if (y < AOThresh) then
               converged = .true.
               exit iters1
            else
               xleft = xright
            end if
         end do iters1
         if (.not. converged) then
            call msg("Subroutine find_radius failed", MSG_ERROR)
            error stop
         end if
         !
         ! Bisection algorithm
         !
         converged = .false.
         iters2: do i = 1, MaxIters
            x = (xleft + xright) / TWO
            y = phi_radial(x, L, NPrimitives, AbsCntrCoeffs, Exponents, NormFactor)
            if (y > AOThresh) then
               xleft = x
            else
               xright = x
            end if
            if (abs(xleft - xright) < XThresh) then
               converged = .true.
               exit iters2
            end if
         end do iters2
         if (.not. converged) then
            call msg("Subroutine find_radius failed", MSG_ERROR)
            error stop
         end if
         radius = (xleft + xright) / TWO
      end subroutine find_radius
   end subroutine basis_R2Max


   pure function basis_NAngFuncCart(L)
      integer :: basis_NAngFuncCart
      integer, intent(in) :: L
      basis_NAngFuncCart =  ((L + 1) * (L + 2)) / 2
   end function basis_NAngFuncCart


   subroutine basis_FuseBasisSets(AOBasisAB, AOBasisA, AOBasisB, System, SpherAO)
      type(TAOBasis), intent(out) :: AOBasisAB
      type(TAOBasis), intent(in)  :: AOBasisA
      type(TAOBasis), intent(in)  :: AOBasisB
      type(TSystem), intent(in)   :: System
      logical, intent(in)         :: SpherAO

      integer :: NShellParamsTotal, MaxNPrimitives, LmaxGTO
      integer :: MaxNAngFuncCart
      integer :: NShells
      integer, dimension(:), allocatable :: ShellMomentum
      integer, dimension(:), allocatable :: NPrimitives
      integer, dimension(:), allocatable :: ShellParamsIdx, ShellCenters
      real(F64), dimension(:), allocatable :: R2Max
      real(F64), dimension(:, :), allocatable :: CntrCoeffs, Exponents, NormFactorsCart
      integer :: Na, Nb, Ma, Mb

      MaxNPrimitives = max(AOBasisA%MaxNPrimitives, AOBasisB%MaxNPrimitives)
      LmaxGTO = max(AOBasisA%LmaxGTO, AOBasisB%LmaxGTO)
      NShellParamsTotal = AOBasisA%NShellParams + AOBasisB%NShellParams

      allocate(NPrimitives(NShellParamsTotal))
      Na = AOBasisA%NShellParams
      Nb = AOBasisB%NShellParams
      NPrimitives(1:Na) = AOBasisA%NPrimitives(1:Na)
      NPrimitives(Na+1:Na+Nb) = AOBasisB%NPrimitives(1:Nb)

      allocate(ShellMomentum(NShellParamsTotal))
      ShellMomentum(1:Na) = AOBasisA%ShellMomentum(1:Na)
      ShellMomentum(Na+1:Na+Nb) = AOBasisB%ShellMomentum(1:Nb)

      allocate(CntrCoeffs(MaxNPrimitives, NShellParamsTotal))
      Ma = AOBasisA%MaxNPrimitives
      Mb = AOBasisB%MaxNPrimitives
      CntrCoeffs = ZERO
      CntrCoeffs(1:Ma, 1:Na) = AOBasisA%CntrCoeffs(1:Ma, 1:Na)
      CntrCoeffs(1:Mb, Na+1:Na+Nb) = AOBasisB%CntrCoeffs(1:Mb, 1:Nb)

      allocate(Exponents(MaxNPrimitives, NShellParamsTotal))
      Exponents = ZERO
      Exponents(1:Ma, 1:Na) = AOBasisA%Exponents(1:Ma, 1:Na)
      Exponents(1:Mb, Na+1:Na+Nb) = AOBasisB%Exponents(1:Mb, 1:Nb)

      MaxNAngFuncCart = ((LmaxGTO + 1) * (LmaxGTO + 2)) / 2
      allocate(NormFactorsCart(MaxNAngFuncCart, NShellParamsTotal))
      NormFactorsCart = ZERO
      Ma = ((AOBasisA%LmaxGTO + 1) * (AOBasisA%LmaxGTO + 2)) / 2
      Mb = ((AOBasisB%LmaxGTO + 1) * (AOBasisB%LmaxGTO + 2)) / 2
      NormFactorsCart(1:Ma, 1:Na) = AOBasisA%NormFactorsCart(1:Ma, 1:Na)
      NormFactorsCart(1:Mb, Na+1:Na+Nb) = AOBasisB%NormFactorsCart(1:Mb, 1:Nb)

      allocate(R2Max(NShellParamsTotal))
      R2Max(1:Na) = AOBasisA%R2Max(1:Na)
      R2Max(Na+1:Na+Nb) = AOBasisB%R2Max(1:Nb)

      NShells = AOBasisA%NShells + AOBasisB%NShells
      allocate(ShellParamsIdx(NShells))
      allocate(ShellCenters(NShells))
      Ma = AOBasisA%NShells
      Mb = AOBasisB%NShells
      ShellParamsIdx(1:Ma) = AOBasisA%ShellParamsIdx(1:Ma)
      ShellParamsIdx(Ma+1:Ma+Mb) = AOBasisB%ShellParamsIdx(1:Mb) + Na
      ShellCenters(1:Ma) = AOBasisA%ShellCenters(1:Ma)
      ShellCenters(Ma+1:Ma+Mb) = AOBasisB%ShellCenters(1:Mb)

      call basis_Init_2(AOBasisAB, System%AtomCoords, ShellCenters, ShellParamsIdx, ShellMomentum, &
         NPrimitives, CntrCoeffs, Exponents, NormFactorsCart, R2Max, SpherAO)
      AOBasisAB%Fused = .true.
      !
      ! The metadata about the origin of basis set parameters
      ! is not carried to the fused basis object.
      !
      AOBasisAB%Assignment%Initialized = .false.
   end subroutine basis_FuseBasisSets


   subroutine basis_OAO(MOBasisVecsCart, MOBasisVecsSpher, S_cao, AOBasis, LinDepThresh, Silent)
      !
      ! Build the matrix of linearly-independent molecular-orbital (MO) vectors
      ! in Cartesian and spherical bases.
      !
      ! The matrix of MO vectors serves as the transformation matrix from
      ! the atomic orbital basis (AO) to ortogonalized orbital basis (OAO):
      ! 1. The resulting orbitals are linear combinations of real spherical harmonics
      !    or Cartesian AOs
      ! 2. The linear dependencies are removed, i.e., the eigenvectors of the overlap
      !    matrix corresponding to sk <= LinDepThresh are removed from the final basis set.
      !
      ! The number of columns of MOBasisVecs is the number of linearly independent
      ! vectors of the OAO basis.
      !
      ! This subroutine assumes that the vectors of the Cartesian Gaussian AO basis
      ! are normalized to unity.
      !
      ! Definitions
      ! ------------
      ! V Coao = Cao
      ! V = W U s**(-1/2)
      ! W = rectangular matrix of the transformation from Cartesian Gaussians
      !     to spherical Gaussians
      ! U = eigenvectors of the overlap matrix S in the spherical Gaussian basis.
      !     The columns are the eigenvectors uk for which sk > LinDepThresh.
      ! s**(-1/2) = square matrix with sk**(-1/2) on the diagonal.
      !
      ! How to do AO <-> OAO transformations
      ! ------------------------------------
      ! Transformation of the Kohn-Sham matrix:
      ! Foao = V**T Fao V
      ! Transformation of the density matrix:
      ! Dao = V Doao V**T
      ! Coao and Doao can be computed by diagonalization of Foao.
      !
      ! 1. Lopata, K. and Govind, N. J. Chem. Theory Comput. 7, 1344 (2011);
      !    doi: 10.1021/ct200137z
      !
      real(F64), dimension(:, :), allocatable, intent(out)  :: MOBasisVecsCart
      real(F64), dimension(:, :), allocatable, intent(out)  :: MOBasisVecsSpher
      real(F64), dimension(:, :), contiguous, intent(in)    :: S_cao
      type(TAOBasis), intent(in)                            :: AOBasis
      real(F64), intent(in)                                 :: LinDepThresh
      logical, optional, intent(in)                         :: Silent

      real(F64) :: NormInt, NormFact
      integer :: u, v
      integer :: ss, s, l, p, p0, p1, q, m
      integer :: Nexcl
      integer :: Noao
      integer :: lx, ly, lz
      real(F64), dimension(:, :), allocatable :: w
      real(F64), dimension(:), allocatable :: c_spher
      real(F64), dimension(:, :), allocatable :: work_transf
      real(F64), dimension(:, :), allocatable :: st, W_cao
      real(F64), dimension(:), allocatable :: eig
      logical :: DisplaySummary
      integer :: ThisImage, NImages

      ThisImage = this_image()
      NImages = num_images()
      if (present(Silent)) then
         DisplaySummary = .not. Silent
      else
         DisplaySummary = .true.
      end if
      associate (SpherAO => AOBasis%SpherAO, &
         NAOCart => AOBasis%NAOCart, &
         NAOSpher => AOBasis%NAOSpher, &
         NShells => AOBasis%NShells, &
         LmaxGTO => AOBasis%LmaxGTO, &
         ShellParamsIdx => AOBasis%ShellParamsIdx, &
         ShellMomentum => AOBasis%ShellMomentum, &
         ShellLocCart => AOBasis%ShellLocCart, &
         ShellLocSpher => AOBasis%ShellLocSpher, &
         NAngFuncCart => AOBasis%NAngFuncCart, &
         CartPolyX => AOBasis%CartPolyX, &
         CartPolyY => AOBasis%CartPolyY, &
         CartPolyZ => AOBasis%CartPolyZ, &
         NormFactorsCart => AOBasis%NormFactorsCart, &
         NormFactorsSpher => AOBasis%NormFactorsSpher &
         )
         if (SpherAO) then
            if (ThisImage == 1) then
               allocate(w(NAOCart, NAOSpher))
               allocate(c_spher(numxyz(LmaxGTO)))
               w = ZERO
               q = 1
               do ss = 1, NShells
                  s = ShellParamsIdx(ss)
                  l = ShellMomentum(s)
                  p0 = ShellLocCart(ss)
                  p1 = ShellLocCart(ss) + NAngFuncCart(s) - 1
                  do m = -l, l
                     call rshu(c_spher, l, m)
                     do p = p0, p1
                        lx = CartPolyX(p - p0 + 1, l)
                        ly = CartPolyY(p - p0 + 1, l)
                        lz = CartPolyZ(p - p0 + 1, l)
                        !
                        ! Dividing by NormFactorsCart cancels out the LMN-dependent
                        ! normalization of the AO basis functions.
                        ! (The LMN dependence means that the normalization
                        ! factor of a Cartesian Gaussian function is different for every
                        ! angular function X**l Y**m Z**n.)
                        ! Then the transformation proceeds identically as for
                        ! uncontracted primitives without any LMN-dependent normalization.
                        !
                        w(p, q) = c_spher(lxlylzpos(lx, ly, lz)) / NormFactorsCart(p-p0+1, s)
                     end do
                     !
                     ! Compute the normalization integral to properly normalize solid harmonic
                     ! orbitals to unity. This is important becase the threshold for linear
                     ! dependencies implicitly assumes normalized orbitals.
                     !
                     NormInt = ZERO
                     do u = p0, p1
                        do v = p0, p1
                           NormInt = NormInt + w(u, q) * S_cao(u, v) * w(v, q)
                        end do
                     end do
                     NormFact = ONE / Sqrt(NormInt)
                     w(p0:p1, q) = NormFact * w(p0:p1, q)
                     q = q + 1
                  end do
               end do
               !
               ! Transform the overlap matrix to the spherical basis
               ! and compute the eigenvalues and eigenvectors
               !
               allocate(work_transf(NAOSpher, NAOCart))
               allocate(st(NAOSpher, NAOSpher))
               call real_aTba(st, w, S_cao, work_transf)
               allocate(eig(NAOSpher))
               call symmetric_eigenproblem(eig, st, NAOSpher, .true.)
               !
               ! Find the number of eigenvectors corresponding to small eigenvalues
               ! of the overlap matrix. Note that the eigenvalues are arranged in
               ! increasing order.
               !
               Nexcl = 0
               do p = 1, NAOSpher
                  if (eig(p) <= LinDepThresh) then
                     Nexcl = Nexcl + 1
                  else
                     st(:, p) = st(:, p) / sqrt(eig(p))
                  end if
               end do
               Noao = NAOSpher - Nexcl
               if (DisplaySummary) then
                  call midrule()
                  call msg("Spherical harmonic basis contains " // str(NAOSpher) // " vectors")
                  call msg("Linear independence thresh: Eig(S) > " // str(LinDepThresh, d=1))
               end if
            end if
            if (NImages > 1) then
               call co_broadcast(Noao, source_image=1)
            end if
            !
            ! MOBasisVecsCart <- W (U s**(-1/2))
            !
            allocate(MOBasisVecsCart(NAOCart, Noao))
            if (ThisImage == 1) then
               call real_ab(MOBasisVecsCart, w, st(:, Nexcl+1:Nexcl+Noao))
            end if
            if (NImages > 1) then
               call co_broadcast(MOBasisVecsCart, source_image=1)
            end if
            allocate(MOBasisVecsSpher(NAOSpher, Noao))
            call SpherGTO_TransformVectors(MOBasisVecsSpher, MOBasisVecsCart, LmaxGTO, NormFactorsSpher, &
               NormFactorsCart, ShellLocSpher, ShellLocCart, ShellMomentum, ShellParamsIdx, &
               NAOSpher, NAOCart, NShells, Noao)
         else
            if (ThisImage == 1) then
               allocate(eig(NAOCart))
               allocate(W_cao(NAOCart, NAOCart))
               W_cao = S_cao
               call symmetric_eigenproblem(eig, W_cao, NAOCart, .true.)
               Nexcl = 0
               do p = 1, NAOCart
                  if (eig(p) <= LinDepThresh) then
                     Nexcl = Nexcl + 1
                  else
                     W_cao(:, p) = W_cao(:, p) / sqrt(eig(p))
                  end if
               end do
               Noao = NAOCart - Nexcl
            end if
            if (NImages > 1) then
               call co_broadcast(Noao, source_image=1)
            end if
            allocate(MOBasisVecsCart(NAOCart, Noao))
            allocate(MOBasisVecsSpher(1, 1))
            if (ThisImage == 1) then
               MOBasisVecsCart = W_cao(:, Nexcl+1:Nexcl+Noao)
            end if
            if (NImages > 1) then
               call co_broadcast(MOBasisVecsCart, source_image=1)
            end if
            if (ThisImage == 1) then
               if (DisplaySummary) then
                  call midrule()
                  call msg("Cartesian basis contains " // str(NAOCart) // " vectors")
               end if
            end if
         end if
         if (ThisImage == 1) then
            if (DisplaySummary) then
               if (Nexcl > 0) then
                  call msg("Removed " // str(Nexcl) // " vectors to make MOs non-redundant")
               end if
               call msg("Computations will be carried out in the space of " // str(Noao) // " MOs")
            end if
         end if
      end associate
   end subroutine basis_OAO


   subroutine basis_NonredundantOrthogonal(Qpk, NOAO, Spq, LinDepThresh)
      real(F64), dimension(:, :), allocatable, intent(out) :: Qpk
      integer, intent(out)                                 :: NOAO
      real(F64), dimension(:, :), intent(in)               :: Spq
      real(F64), intent(in)                                :: LinDepThresh

      integer :: NAO, k
      integer :: ErrorCode
      real(F64), dimension(:, :), allocatable :: Wpq
      real(F64), dimension(:), allocatable :: Lambda

      NAO = size(Spq, dim=1)
      allocate(Wpq(NAO, NAO))
      Wpq(:, :) = Spq(:, :)
      allocate(Lambda(NAO))
      call real_EVD(Lambda, Wpq, NAO, .true., &
         Algorithm=LINALG_EVD_ALGORITHM_1, &
         Info=ErrorCode)
      if (ErrorCode > 0) then
         call msg("Eigensolver 1 did not converge in basis_NonredundantOrthogonal", &
            MSG_WARNING)
         call msg("Attempting diagonalization with eigensolver 2", MSG_WARNING)
         !
         ! Overlap matrix has been overwritten by the previous eigensolver
         ! and needs to be recreated before another attempt
         !
         Wpq(:, :) = Spq(:, :)
         call real_EVD(Lambda, Wpq, NAO, .true., &
            Algorithm=LINALG_EVD_ALGORITHM_2, &
            Info=ErrorCode)
         if (ErrorCode > 0) then
            call msg("Eigensolver 2 did not converge in basis_NonredundantOrthogonal", &
               MSG_ERROR)
            error stop
         end if
      end if
      NOAO = 0
      do k = NAO, 1, -1
         if (Lambda(k) > LinDepThresh) then
            NOAO = NOAO + 1
         else
            exit
         end if
      end do
      if (NOAO == 0) then
         call msg("Failed to generate OAO vectors for LinDepThresh=" // str(LinDepThresh,d=1), MSG_ERROR)
         error stop
      end if
      allocate(Qpk(NAO, NOAO))
      do k = 1, NOAO
         Qpk(:, k) = Wpq(:, NAO-k+1) / Sqrt(Lambda(NAO-k+1))
      end do
   end subroutine basis_NonredundantOrthogonal
end module basis_sets
