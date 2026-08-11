module basis_definitions
   use arithmetic
   use string
   use periodic
   use display
   use io, only: io_exists, DIRSEP

   implicit none

   !
   ! Order orbital shells by their effective orbital radius
   !
   integer, parameter :: SHELL_ORDER_BY_RADIUS = 0
   !
   ! Order orbital shells by their angular momentum
   !
   integer, parameter :: SHELL_ORDER_BY_MOMENTUM = 1
   !
   ! Keep the original order of shells as defined in the basis set parameters file
   !
   integer, parameter :: SHELL_ORDER_FIXED = 2

   type TBasisRule
      !
      ! Atom index, element atomic number,
      ! or zero for a global fallback rule
      ! applied to the entire system
      !
      integer :: id
      !
      ! Path to the basis set parameters file,
      ! which can be part of the library or a
      ! user-provided file
      !
      character(:), allocatable :: PathToParams
      !
      ! Internal basis set name or filename base
      !
      character(:), allocatable :: BaseName
      !
      ! Human-readable name of the basis set
      !
      character(:), allocatable :: DisplayedName
      !
      ! Path to the directory containing guess density
      ! files for the given basis set. Allocated only
      ! if guesses are available for this basis set.
      !
      character(:), allocatable :: PathToGuessDir
      logical :: GuessAvailable = .false.
      !
      ! Flag indicating if the basis set is loaded
      ! from the library directory where all basis
      ! sets are stored
      !
      logical :: FromLibrary = .false.
   end type TBasisRule

   type TBasisAssignment
      !
      ! Mapping between atoms or elements and basis set files
      !
      ! Rules are evaluated in the following order of priority:
      ! 1. AtomRules: Explicit atom indices (e.g., atom 3)
      ! 2. ElementRules: Element atomic numbers (e.g., Z=8 for O)
      ! 3. GlobalFallback: Defined via '*' in the `basis_assignment` block
      !
      logical :: Initialized = .false.

      type(TBasisRule), allocatable :: AtomRules(:)
      integer :: NAtomRules = 0

      type(TBasisRule), allocatable :: ElementRules(:)
      integer :: NElementRules = 0

      type(TBasisRule) :: GlobalFallback
      logical :: FallbackAvailable = .false.

      character(:), allocatable :: LibraryDir
      character(:), allocatable :: GuessDir
   contains
      procedure :: set_library_dir => basis_set_library_dir
      procedure :: set_guess_dir => basis_set_guess_dir
      procedure :: add_atom_rule => basis_add_atom_rule
      procedure :: add_element_rule => basis_add_element_rule
      procedure :: add_global_fallback => basis_add_global_fallback
      procedure :: read_line => basis_ReadAssignLine
      procedure :: assign => basis_assignment_copy
      procedure :: is_uniform => basis_is_uniform
      procedure :: get_atom_rule => basis_get_atom_rule
      generic :: assignment(=) => assign
   end type TBasisAssignment

   type TAOBasis
      !
      ! Gaussian atomic orbital basis set parameters and mapping arrays.
      !
      ! A contracted Cartesian Gaussian orbital centered on atom c at Rc has the form:
      !   Phi(r) = N * (x-Xc)**lx * (y-Yc)**ly * (z-Zc)**lz * Sum_i [ c_i * exp(-alpha_i * |r-Rc|**2) ]
      !
      ! Spherical AOs are formed by contracting Cartesian Gaussians with solid
      ! harmonic coefficients. The normalization factor N depends on whether
      ! the final AO is Cartesian or spherical (see docstrings in integrals/Auto2e).
      !
      ! Parameters:
      !   N       : Normalization constant from NormFactorsCart or NormFactorsSpher.
      !   lx,ly,lz: Cartesian polynomial exponents from CartPolyX/Y/Z.
      !   c_i     : Contraction coefficients from CntrCoeffs.
      !   alpha_i : Gaussian exponents from Exponents.
      !   Rc      : (Xc, Yc, Zc) coordinates of atom c from AtomCoords.
      !
      ! Scalars:
      !   SpherAO          : True if spherical AOs are preferred, but parameters for both spherical
      !                      and Cartesian AOs are present in this structure.
      !   NShellParams     : Total number of unique shell parameter sets. Using atom-specific basis
      !                      assignments can result in different basis sets for atoms of the same
      !                      element, increasing this number.
      !   NShells          : Total number of shells instantiated on atoms.
      !   LmaxGTO          : Maximum angular momentum in the basis.
      !   MaxNPrimitives   : Maximum primitives in any single shell.
      !   NAOSpher         : Total number of spherical basis functions.
      !   NAOCart          : Total number of Cartesian basis functions.
      !   NAtoms           : Number of atoms.
      !   MaxNShells       : Maximum shells assigned to any single atom.
      !   Fused            : True if this basis set was created by joining two TAOBasis objects.
      !
      ! Arrays:
      !   AtomCoords       : (3, NAtoms) Cartesian coordinates of atoms.
      !   ShellCenters     : (NShells) Atom index each shell is centered on.
      !   ShellParamsIdx   : (NShells) Index mapping a shell to its unique parameter set.
      !   ShellMomentum    : (NShellParams) Angular momentum quantum number (L).
      !   AtomShellMap     : (2, MaxNShellSegments, NAtoms) Start and end shell indices for
      !                      contiguous segments.
      !   AtomShellN       : (NAtoms) Number of contiguous shell segments per atom.
      !                      This can be greater than one only if it's a fused basis set
      !                      created by joining two TAOBasis objects with basis_FuseBasisSets.
      !   NPrimitives      : (NShellParams) Number of Gaussian primitives per shell.
      !   CntrCoeffs       : (MaxNPrimitives, NShellParams) Contraction coefficients.
      !   Exponents        : (MaxNPrimitives, NShellParams) Gaussian exponents.
      !   NormFactorsCart  : (MaxNAngFuncCart, NShellParams) Normalization constants for
      !                      Cartesian functions.
      !   NormFactorsSpher : (MaxNAngFuncSpher, NShellParams) Normalization constants for
      !                      spherical functions.
      !   NAngFuncSpher    : (NShellParams) Number of spherical functions (2L+1).
      !   NAngFuncCart     : (NShellParams) Number of Cartesian functions.
      !   ShellLocSpher    : (NShells) Global starting index in the spherical basis.
      !   ShellLocCart     : (NShells) Global starting index in the Cartesian basis.
      !   CartPolyX/Y/Z    : (MaxNAngFuncCart, 0:LmaxGTO) Cartesian polynomial exponents (lx, ly, lz)
      !                      defining the angular part x**lx * y**ly * z**lz.
      !   R2Max            : (NShellParams) Effective squared radial extent of an orbital used
      !                      for grid screening.
      !   MaxAtomL         : (NAtoms) Maximum angular momentum on a given atom.
      !
      ! Objects:
      !   Assignment       : Stores file paths to basis set parameters for each atom.
      !                      Allows different basis sets for atoms of the same element.
      !                      Used as metadata to track the origin of basis set parameters.
      !
      real(F64), dimension(:, :), allocatable :: AtomCoords
      integer, dimension(:), allocatable :: ShellCenters
      integer, dimension(:), allocatable :: ShellParamsIdx
      integer, dimension(:), allocatable :: ShellMomentum
      integer, dimension(:, :, :), allocatable :: AtomShellMap
      integer, dimension(:), allocatable :: AtomShellN
      integer, dimension(:), allocatable :: NPrimitives
      real(F64), dimension(:, :), allocatable :: CntrCoeffs
      real(F64), dimension(:, :), allocatable :: Exponents
      real(F64), dimension(:, :), allocatable :: NormFactorsCart
      real(F64), dimension(:, :), allocatable :: NormFactorsSpher
      integer, dimension(:), allocatable :: NAngFuncSpher
      integer, dimension(:), allocatable :: NAngFuncCart
      integer, dimension(:), allocatable :: ShellLocSpher
      integer, dimension(:), allocatable :: ShellLocCart
      integer, dimension(:, :), allocatable :: CartPolyX
      integer, dimension(:, :), allocatable :: CartPolyY
      integer, dimension(:, :), allocatable :: CartPolyZ
      real(F64), dimension(:), allocatable :: R2Max
      integer, dimension(:), allocatable :: MaxAtomL
      type(TBasisAssignment) :: Assignment
      logical :: SpherAO
      logical :: Fused = .false.
      integer :: NShellParams
      integer :: NShells
      integer :: LmaxGTO
      integer :: MaxNPrimitives
      integer :: NAOSpher
      integer :: NAOCart
      integer :: NAtoms
      integer :: MaxNShells
   contains
      procedure :: display => basis_display
   end type TAOBasis

contains

   subroutine basis_display(this)
      class(TAOBasis), intent(in) :: this

      integer :: i, z_num, atom_idx
      integer :: n_items, current_item
      character(:), allocatable :: param_str
      character(2)  :: sym
      integer, parameter :: w = 45 ! Control the width of the label column

      call msg("Atomic orbital basis set")
      call dotted_separator(76)

      call msg(lfield("atomic centers", w) // str(this%NAtoms))
      call msg(lfield("shell parameter sets", w) // str(this%NShellParams))
      call msg(lfield("spherical AOs", w) // str(this%NAOSpher))
      call msg(lfield("max angular momentum", w) // str(this%LmaxGTO))

      if (this%Assignment%is_uniform()) then
         if (this%Assignment%FallbackAvailable) then
            call msg(lfield("parameter set", w) // this%Assignment%GlobalFallback%DisplayedName)
         else
            if (this%Assignment%NElementRules > 0) then
               call msg(lfield("parameter set", w) // this%Assignment%ElementRules(1)%DisplayedName)
            else if (this%Assignment%NAtomRules > 0) then
               call msg(lfield("parameter set", w) // this%Assignment%AtomRules(1)%DisplayedName)
            end if
         end if
      else
         call msg("parameter set")

         n_items = 0
         if (this%Assignment%FallbackAvailable) n_items = n_items + 1
         n_items = n_items + this%Assignment%NElementRules + this%Assignment%NAtomRules

         current_item = 0

         if (this%Assignment%FallbackAvailable) then
            current_item = current_item + 1
            if (current_item < n_items) then
               call msg(lfield("   ├─ default", w) // this%Assignment%GlobalFallback%DisplayedName)
            else
               call msg(lfield("   └─ default", w) // this%Assignment%GlobalFallback%DisplayedName)
            end if
         end if

         do i = 1, this%Assignment%NElementRules
            current_item = current_item + 1
            z_num = this%Assignment%ElementRules(i)%id
            sym = ELNAME_SHORT(z_num)
            if (current_item < n_items) then
               param_str = "   ├─ " // trim(sym)
            else
               param_str = "   └─ " // trim(sym)
            end if
            call msg(lfield(param_str, w) // this%Assignment%ElementRules(i)%DisplayedName)
         end do

         do i = 1, this%Assignment%NAtomRules
            current_item = current_item + 1
            atom_idx = this%Assignment%AtomRules(i)%id
            if (current_item < n_items) then
               param_str = "   ├─ atom #" // str(atom_idx)
            else
               param_str = "   └─ atom #" // str(atom_idx)
            end if
            call msg(lfield(param_str, w) // this%Assignment%AtomRules(i)%DisplayedName)
         end do
      end if
      call dotted_separator(76)
   end subroutine basis_display


   subroutine basis_add_atom_rule(BasisAssign, a, Rule)
      class(TBasisAssignment), intent(inout) :: BasisAssign
      integer, intent(in)                    :: a
      type(TBasisRule), intent(in)           :: Rule
      type(TBasisRule), allocatable :: temp(:)
      integer :: i

      if (allocated(BasisAssign%AtomRules)) then
         allocate(temp(BasisAssign%NAtomRules + 1))
         do i = 1, BasisAssign%NAtomRules
            temp(i) = BasisAssign%AtomRules(i)
         end do
         temp(BasisAssign%NAtomRules + 1) = Rule
         temp(BasisAssign%NAtomRules + 1)%id = a
         call move_alloc(temp, BasisAssign%AtomRules)
      else
         allocate(BasisAssign%AtomRules(1))
         BasisAssign%AtomRules(1) = Rule
         BasisAssign%AtomRules(1)%id = a
      end if
      BasisAssign%NAtomRules = BasisAssign%NAtomRules + 1
      BasisAssign%Initialized = .true.
   end subroutine basis_add_atom_rule


   subroutine basis_add_element_rule(BasisAssign, Z, Rule)
      class(TBasisAssignment), intent(inout) :: BasisAssign
      integer, intent(in)                    :: Z
      type(TBasisRule), intent(in)           :: Rule
      type(TBasisRule), allocatable :: temp(:)
      integer :: i

      if (allocated(BasisAssign%ElementRules)) then
         allocate(temp(BasisAssign%NElementRules + 1))
         do i = 1, BasisAssign%NElementRules
            temp(i) = BasisAssign%ElementRules(i)
         end do
         temp(BasisAssign%NElementRules + 1) = Rule
         temp(BasisAssign%NElementRules + 1)%id = Z
         call move_alloc(temp, BasisAssign%ElementRules)
      else
         allocate(BasisAssign%ElementRules(1))
         BasisAssign%ElementRules(1) = Rule
         BasisAssign%ElementRules(1)%id = Z
      end if
      BasisAssign%NElementRules = BasisAssign%NElementRules + 1
      BasisAssign%Initialized = .true.
   end subroutine basis_add_element_rule

   subroutine basis_set_library_dir(BasisAssign, path)
      class(TBasisAssignment), intent(inout) :: BasisAssign
      character(*), intent(in)               :: path
      BasisAssign%LibraryDir = trim(path)
      if (len(BasisAssign%LibraryDir) > 0) then
         if (.not. endswith(BasisAssign%LibraryDir, DIRSEP)) &
            BasisAssign%LibraryDir = BasisAssign%LibraryDir // DIRSEP
      end if
   end subroutine basis_set_library_dir

   subroutine basis_set_guess_dir(BasisAssign, path)
      class(TBasisAssignment), intent(inout) :: BasisAssign
      character(*), intent(in)               :: path
      BasisAssign%GuessDir = trim(path)
      if (len(BasisAssign%GuessDir) > 0) then
         if (.not. endswith(BasisAssign%GuessDir, DIRSEP)) &
            BasisAssign%GuessDir = BasisAssign%GuessDir // DIRSEP
      end if
   end subroutine basis_set_guess_dir

   subroutine basis_add_global_fallback(BasisAssign, Rule)
      class(TBasisAssignment), intent(inout) :: BasisAssign
      type(TBasisRule), intent(in)           :: Rule

      BasisAssign%GlobalFallback = Rule
      BasisAssign%FallbackAvailable = .true.
      BasisAssign%Initialized = .true.
   end subroutine basis_add_global_fallback

   subroutine basis_ReadAssignLine(BasisAssign, line)
      class(TBasisAssignment), intent(inout) :: BasisAssign
      character(*), intent(in)               :: line

      character(:), allocatable :: key, val
      type(TBasisRule)          :: ResolvedRule
      integer :: start_idx, end_idx, a, Z

      call split(line, key, val)

      call basis_ResolvePath(ResolvedRule, BasisAssign, val)

      if (key == "*" .or. uppercase(key) == "BASIS") then
         ResolvedRule%id = 0
         call BasisAssign%add_global_fallback(ResolvedRule)
      else if (index(key, "-") > 0) then
         ! Atom range (e.g., 1-3)
         read(key(1:index(key,"-")-1), *) start_idx
         read(key(index(key,"-")+1:), *) end_idx
         do a = start_idx, end_idx
            call BasisAssign%add_atom_rule(a, ResolvedRule)
         end do
      else if (isinteger(key)) then
         ! Single Atom ID
         read(key, *) a
         call BasisAssign%add_atom_rule(a, ResolvedRule)
      else
         ! Element symbol
         Z = znumber_short(key)
         if (Z > 0) then
            call BasisAssign%add_element_rule(Z, ResolvedRule)
         else
            call msg("Unknown element symbol '" // key // &
               "' in basis_assignment block.", priority=MSG_ERROR)
            error stop
         end if
      end if
   end subroutine basis_ReadAssignLine

   subroutine basis_assignment_copy(lhs, rhs)
      class(TBasisAssignment), intent(out) :: lhs
      class(TBasisAssignment), intent(in)  :: rhs

      lhs%Initialized = rhs%Initialized
      lhs%NAtomRules = rhs%NAtomRules
      lhs%NElementRules = rhs%NElementRules

      if (rhs%NAtomRules > 0) then
         allocate(lhs%AtomRules, source=rhs%AtomRules)
      end if

      if (rhs%NElementRules > 0) then
         allocate(lhs%ElementRules, source=rhs%ElementRules)
      end if

      if (rhs%FallbackAvailable) then
         lhs%GlobalFallback = rhs%GlobalFallback
         lhs%FallbackAvailable = .true.
      end if

      if (allocated(rhs%LibraryDir)) then
         lhs%LibraryDir = rhs%LibraryDir
      end if

      if (allocated(rhs%GuessDir)) then
         lhs%GuessDir = rhs%GuessDir
      end if
   end subroutine basis_assignment_copy


   function basis_is_uniform(this)
      !
      ! Check if the basis assignment is globally uniform, that is,
      ! without special cases for selected elements or atoms
      !
      logical                             :: basis_is_uniform
      class(TBasisAssignment), intent(in) :: this
      integer                             :: k
      character(:), allocatable           :: RefPath

      if (.not. this%Initialized) then
         call msg("basis_is_uniform: TBasisAssignment is not initialized", MSG_ERROR)
         error stop
      end if

      basis_is_uniform = .true.
      if (this%NAtomRules == 0 .and. this%NElementRules == 0) return

      if (this%FallbackAvailable) then
         RefPath = this%GlobalFallback%PathToParams
      else if (this%NElementRules > 0) then
         RefPath = this%ElementRules(1)%PathToParams
      else
         RefPath = this%AtomRules(1)%PathToParams
      end if

      do k = 1, this%NAtomRules
         if (this%AtomRules(k)%PathToParams /= RefPath) basis_is_uniform = .false.
      end do
      do k = 1, this%NElementRules
         if (this%ElementRules(k)%PathToParams /= RefPath) basis_is_uniform = .false.
      end do
   end function basis_is_uniform


   function basis_get_atom_rule(this, atom_idx, Z)
      !
      ! Extract basis rule for atom.
      !
      ! Resolution priorities:
      ! 1. Atom-specific rules.
      ! 2. Element-specific rules.
      ! 3. Global fallback.
      !
      ! Note that in most cases there is only a single uniform basis set
      ! assigned to the entire system, which corresponds to the global
      ! fallback without any other atom-specific or element-specific exceptions.
      !
      type(TBasisRule)                    :: basis_get_atom_rule
      class(TBasisAssignment), intent(in) :: this
      integer, intent(in)                 :: atom_idx
      integer, intent(in)                 :: Z
      integer                             :: c

      if (.not. this%Initialized) then
         call msg("basis_get_atom_rule: TBasisAssignment is not initialized.", MSG_ERROR)
         error stop
      end if

      do c = 1, this%NAtomRules
         if (this%AtomRules(c)%id == atom_idx) then
            basis_get_atom_rule = this%AtomRules(c)
            return
         end if
      end do

      do c = 1, this%NElementRules
         if (this%ElementRules(c)%id == Z) then
            basis_get_atom_rule = this%ElementRules(c)
            return
         end if
      end do

      if (this%FallbackAvailable) then
         basis_get_atom_rule = this%GlobalFallback
      else
         call msg("No basis set assigned for atom " // str(atom_idx) // &
            " (" // trim(ELNAME_SHORT(Z)) // ")", MSG_ERROR)
         error stop
      end if
   end function basis_get_atom_rule


   subroutine basis_ResolvePath(Rule, BasisAssign, ValString)
      !
      ! Map an EMSL basis set name or custom file path to an absolute path.
      !
      ! Parameters:
      !   Rule                : Output basis rule object.
      !   BasisAssign         : Input TBasisAssignment object containing LibraryDir and GuessDir state.
      !   ValString           : Input basis string, either an EMSL alias or 'FILE path'.
      !
      type(TBasisRule), intent(out)       :: Rule
      class(TBasisAssignment), intent(in) :: BasisAssign
      character(*), intent(in)            :: ValString

      character(:), allocatable :: a1, a2, p, n

      call split(ValString, a1, a2)
      if (uppercase(a1) == "FILE") then
         if (io_exists(a2)) then
            p = a2
            n = a2
            Rule%BaseName = a2
            Rule%FromLibrary = .false.
         else
            call msg("Basis set coefficients file is inaccessible: " // a2, MSG_ERROR)
            error stop
         end if
      else
         ! Standard EMSL basis sets mapped here
         select case(uppercase(ValString))
          case ("AHLRICHS-VDZ")
            p = "ahlrichs-vdz"
            n = "Ahlrichs VDZ"
          case ("3-21G")
            p = "3-21g"
            n = "3-21G"
          case ("3-21++G")
            p = "3-21g_diff_all"
            n = "3-21++G"
          case ("6-31G")
            p = "6-31g"
            n = "6-31G"
          case ("6-31G**")
            p = "6-31g_pol_all"
            n = "6-31G**"
          case ("6-31++G")
            p = "6-31g_diff_all"
            n = "6-31++G"
          case ("6-31+G**")
            p = "6-31g_diff_pol_all"
            n = "6-31+G**"
          case ("6-31++G**")
            p = "6-31g_diff_all_pol_all"
            n = "6-31++G**"
          case ("6-31G(3DF,3PD)", "6-31G(3DF, 3PD)")
            p = "6-31g_3df_3pd"
            n = "6-31G(3df,3dp)"
          case ("6-311G")
            p = "6-311g"
            n = "6-311G"
          case ("6-311G**")
            p = "6-311g_pol_all"
            n = "6-311G**"
          case ("6-311+G*")
            p = "6-311g_diff_pol"
            n = "6-311+G*"
          case ("6-311++G**")
            p = "6-311g_diff_all_pol_all"
            n = "6-311++G**"
          case ("6-311++G(2D,2P)", "6-311++G(2D, 2P)")
            p = "6-311g_diff_all_2d_2p"
            n = "6-311++G(2d,2p)"
          case ("6-311++G(3DF,3PD)", "6-311++G(3DF, 3PD)")
            p = "6-311g_diff_all_3df_3pd"
            n = "6-311++G(3df,3pd)"
          case ("6-311(3+,3+)G**", "6-311(3+, 3+)G**", "6-311(3+,3+)G(d,p)", &
             "6-311(3+, 3+)G(d, p)")
            p = "6-311g_triple_diff_all_pol_all"
            n = "6-311(3+,3+)G**"
          case ("CC-PVDZ", "CC-PV(D+D)Z")
            p = "cc-pvddz"
            n = "cc-pVDZ"
          case ("CC-PVDZ_OLD", "CC-PVDZ-OLD")
            p = "cc-pvdz"
            n = "cc-pVDZ (old)"
          case ("CC-PVDZ-PP")
            p = "cc-pvdz-pp"
            n = "cc-pVDZ-PP"
          case ("AUG-CC-PVDZ", "AUG-CC-PV(D+D)Z")
            p = "aug-cc-pvddz"
            n = "aug-cc-pVDZ"
          case ("AUG-CC-PVDZ_OLD", "AUG-CC-PVDZ-OLD")
            p = "aug-cc-pvdz"
            n = "aug-cc-pVDZ (old)"
          case ("AUG-CC-PVDZ-PP")
            p = "aug-cc-pvdz-pp"
            n = "aug-cc-pVDZ-PP"
          case ("D-AUG-CC-PVDZ")
            p = "d-aug-cc-pvdz"
            n = "d-aug-cc-pVDZ"
          case ("CC-PVTZ", "CC-PV(T+D)Z")
            p = "cc-pvtdz"
            n = "cc-pVTZ"
          case ("CC-PVTZ_OLD", "CC-PVTZ-OLD")
            p = "cc-pvtz"
            n = "cc-pVTZ (old)"
          case ("CC-PVTZ-PP")
            p = "cc-pvtz-pp"
            n = "cc-pVTZ-PP"
          case ("AUG-CC-PVTZ", "AUG-CC-PV(T+D)Z")
            p = "aug-cc-pvtdz"
            n = "aug-cc-pVTZ"
          case ("AUG-CC-PVTZ_OLD", "AUG-CC-PVTZ-OLD")
            p = "aug-cc-pvtz"
            n = "aug-cc-pVTZ (old)"
          case ("AUG-CC-PVTZ-PP")
            p = "aug-cc-pvtz-pp"
            n = "aug-cc-pVTZ-PP"
          case ("D-AUG-CC-PVTZ")
            p = "d-aug-cc-pvtz"
            n = "d-aug-cc-pVTZ"
          case ("CC-PVQZ", "CC-PV(Q+D)Z")
            p = "cc-pvqdz"
            n = "cc-pVQZ"
          case ("CC-PVQZ_OLD", "CC-PVQZ-OLD")
            p = "cc-pvqz"
            n = "cc-pVQZ (old)"
          case ("AUG-CC-PVQZ", "AUG-CC-PV(Q+D)Z")
            p = "aug-cc-pvqdz"
            n = "aug-cc-pVQZ"
          case ("AUG-CC-PVQZ_OLD", "AUG-CC-PVQZ-OLD")
            p = "aug-cc-pvqz"
            n = "aug-cc-pVQZ (old)"
          case ("D-AUG-CC-PVQZ")
            p = "d-aug-cc-pvqz"
            n = "d-aug-cc-pVQZ"
          case ("CC-PV5Z", "CC-PV(5+D)Z")
            p = "cc-pv5dz"
            n = "cc-pV5Z"
          case ("CC-PV5Z_OLD", "CC-PV5Z-OLD")
            p = "cc-pv5z"
            n = "cc-pV5Z (old)"
          case ("AUG-CC-PV5Z", "AUG-CC-PV(5+D)Z")
            p = "aug-cc-pv5dz"
            n = "aug-cc-pV5Z"
          case ("AUG-CC-PV5Z_OLD", "AUG-CC-PV5Z-OLD")
            p = "aug-cc-pv5z"
            n = "aug-cc-pV5Z (old)"
          case ("D-AUG-CC-PV5Z")
            p = "d-aug-cc-pv5z"
            n = "d-aug-cc-pV5Z"
          case ("CC-PCVTZ")
            p = "cc-pcvtz"
            n = "cc-pCVTZ"
          case ("AUG-CC-PCVTZ")
            p = "aug-cc-pcvtz"
            n = "aug-cc-pCVTZ"
          case ("CC-PCVQZ")
            p = "cc-pcvqz"
            n = "cc-pCVQZ"
          case ("AUG-CC-PCVQZ")
            p = "aug-cc-pcvqz"
            n = "aug-cc-pCVQZ"
            ! ----------------------------------------------
            !               (aug-)cc-pwCXZ
            ! ----------------------------------------------
          case ("CC-PWCVQZ")
            p = "cc-pwcvqz"
            n = "cc-pwCVQZ"
          case ("AUG-CC-PWCVQZ")
            p = "aug-cc-pwcvqz"
            n = "aug-cc-pwCVQZ"
          case ("CC-PWCV5Z")
            p = "cc-pwcv5z"
            n = "cc-pwCV5Z"
          case ("AUG-CC-PWCV5Z")
            p = "aug-cc-pwcv5z"
            n = "aug-cc-pwCV5Z"
          case ("DEF2-QZVP")
            p = "def2-qzvp"
            n = "Def2-QZVP"
          case ("DEF2-QZVPD")
            p = "def2-qzvpd"
            n = "Def2-QZVPD"
          case ("DEF2-QZVPP")
            p = "def2-qzvpp"
            n = "Def2-QZVPP"
          case ("DEF2-QZVPPD")
            p = "def2-qzvppd"
            n = "Def2-QZVPPD"
          case ("DEF2-SV(P)")
            p = "def2-sv_p"
            n = "Def2-SV(P)"
          case ("DEF2-SVP")
            p = "def2-svp"
            n = "Def2-SVP"
          case ("DEF2-SVPD")
            p = "def2-svpd"
            n = "Def2-SVPD"
          case ("DEF2-TZVP")
            p = "def2-tzvp"
            n = "Def2-TZVP"
          case ("DEF2-TZVPD")
            p = "def2-tzvpd"
            n = "Def2-TZVPD"
          case ("DEF2-TZVPP")
            p = "def2-tzvpp"
            n = "Def2-TZVPP"
          case ("DEF2-TZVPPD")
            p = "def2-tzvppd"
            n = "Def2-TZVPPD"
          case ("SADLEJ-PVTZ")
            p = "sadlej-pvtz"
            n = "Sadlej-pVTZ"
          case ("SAPPORO-QZP-2012")
            p = "sapporo-qzp-2012"
            n = "Sapporo-QZP-2012"
          case ("CRENBL")
            p = "crenbl"
            n = "CRENBL"
          case default
            call msg("Unknown basis set", MSG_ERROR)
            error stop
         end select

         Rule%BaseName = p
         Rule%FromLibrary = .true.

         if (allocated(BasisAssign%LibraryDir)) then
            p = BasisAssign%LibraryDir // p // ".txt"
         else
            call msg("Basis library directory is not specified", MSG_ERROR)
            error stop
         end if
      end if

      if (allocated(BasisAssign%GuessDir) .and. Rule%FromLibrary) then
         Rule%GuessAvailable = .true.
         !
         ! GuessDir points to a root directory containing guess density files.
         ! Each basis set has its own subdirectory under GuessDir named after
         ! the BaseName of the basis set. The guess files for individual elements
         ! will reside inside this subdirectory.
         !
         Rule%PathToGuessDir = BasisAssign%GuessDir // Rule%BaseName // "/"
      else
         Rule%GuessAvailable = .false.
      end if

      Rule%PathToParams = p
      Rule%DisplayedName = n
   end subroutine basis_ResolvePath
end module basis_definitions
