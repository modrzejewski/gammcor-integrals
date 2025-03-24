module sys_definitions
      use arithmetic
      use periodic
      use string
      use sort
      use display
      use io
      
      implicit none
      !
      ! Units used in the xyz file
      !
      integer, parameter :: SYS_UNITS_ANGSTROM = 1
      integer, parameter :: SYS_UNITS_BOHR = 2      
      !
      ! Subsystems required for single points, two-body interaction energies,
      ! and three-body nonadditive interaction energies
      !
      integer, parameter :: SYS_TOTAL = 1
      integer, parameter :: SYS_MONO_A = 2
      integer, parameter :: SYS_MONO_B = 3
      integer, parameter :: SYS_MONO_C = 4
      integer, parameter :: SYS_DIMER_AB = 5
      integer, parameter :: SYS_DIMER_BC = 6
      integer, parameter :: SYS_DIMER_AC = 7
      !
      ! Subsystems required for four-body nonadditive interaction energies
      !
      integer, parameter :: SYS_MONO_D = 8
      integer, parameter :: SYS_DIMER_AD = 9
      integer, parameter :: SYS_DIMER_BD = 10
      integer, parameter :: SYS_DIMER_CD = 11
      integer, parameter :: SYS_TRIMER_ABC = 12
      integer, parameter :: SYS_TRIMER_ABD = 13
      integer, parameter :: SYS_TRIMER_ACD = 14
      integer, parameter :: SYS_TRIMER_BCD = 15

      integer, parameter :: SYS_NONE = 0
      integer, parameter :: SYS_MOLECULE = 1
      integer, parameter :: SYS_DIMER = 2
      integer, parameter :: SYS_TRIMER = 3
      integer, parameter :: SYS_TETRAMER = 4

      integer, parameter :: SYS_ALL_ATOMS = 1
      integer, parameter :: SYS_REAL_ATOMS = 2
      integer, parameter :: SYS_GHOST_ATOMS = 3

      type TSystem
            real(F64), dimension(:, :), allocatable :: AtomCoords
            !
            ! Nuclear charges (correspond to physical nuclei,
            ! not affected by pseudopotentials)
            !
            integer, dimension(:), allocatable :: ZNumbers
            !
            ! Effective nuclear charges applied when
            ! a pseudopotential is present
            !
            integer, dimension(:), allocatable :: ZNumbersECP
            logical :: ECPCharges = .false.
            !
            ! Additional point charges used to generate an electrostatic
            ! field for the molecule. Can be used, e.g., to remove energy
            ! level degeneracy.
            !
            integer :: NPointCharges = 0
            real(F64), dimension(:), allocatable :: PointCharges
            real(F64), dimension(:, :), allocatable :: PointChargeCoords
            !
            integer :: Mult = 1
            integer :: Charge = 0
            integer :: NElectrons
            integer :: NAtoms = 0
            integer, dimension(2, 2) :: RealAtoms
            integer, dimension(4) :: SubsystemAtoms = [0, 0, 0, 0]
            integer, dimension(4) :: SubsystemCharges = [0, 0, 0, 0]
            integer, dimension(15) :: SubsystemMult = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            integer :: SystemKind = SYS_NONE
            integer :: SubsystemKind = SYS_NONE
            real(F64), dimension(:, :), allocatable :: SortedDistances
            integer, dimension(:, :), allocatable :: SortedDistancesIdx
      end type TSystem

contains

      subroutine sys_Init(System, i)
            type(TSystem), intent(inout) :: System
            integer, intent(in)          :: i

            integer :: p0, p1, q0, q1, s
            integer :: a0, a1, b0, b1, c0, c1, d0, d1
            
            associate (RealAtoms=>System%RealAtoms, NAtoms=>System%NAtoms, Charge=>System%Charge, &
                  SubsystemAtoms=>System%SubsystemAtoms, Mult=>System%Mult, SubsystemMult=>System%SubsystemMult, &
                  SubsystemCharges=>System%SubsystemCharges, SystemKind=>System%SystemKind, ZNumbers=>System%ZNumbers, &
                  NElectrons=>System%NElectrons, SubsystemKind=>System%SubsystemKind, ECPCharges=>System%ECPCharges, &
                  ZNumbersECP=>System%ZNumbersECP)
                  SubsystemKind = i
                  Charge = sum(SubsystemCharges)
                  NAtoms = sum(SubsystemAtoms)
                  p0 = 1
                  p1 = NAtoms
                  !
                  ! When q1<q0, the loops over atoms q0...q1 perform zero cycles.
                  ! That's why q0 and q1 are modified only for non-contiguous
                  ! ranges of atoms.
                  !
                  q0 = 1
                  q1 = 0
                  select case (SystemKind)
                  case (SYS_DIMER)
                        Charge = sum(SubsystemCharges(1:2))
                        select case (i)
                        case (SYS_MONO_A)
                              p0 = 1
                              p1 = SubsystemAtoms(1)
                              Charge = SubsystemCharges(1)
                        case (SYS_MONO_B)
                              p0 = SubsystemAtoms(1) + 1
                              p1 = NAtoms
                              Charge = SubsystemCharges(2)
                        end select
                  case (SYS_TRIMER)
                        Charge = sum(SubsystemCharges(1:3))
                        select case (i)
                        case (SYS_MONO_A)
                              p0 = 1
                              p1 = SubsystemAtoms(1)
                              Charge = SubsystemCharges(1)
                        case (SYS_MONO_B)
                              p0 = SubsystemAtoms(1) + 1
                              p1 = SubsystemAtoms(1) + SubsystemAtoms(2)
                              Charge = SubsystemCharges(2)
                        case (SYS_MONO_C)
                              p0 = SubsystemAtoms(1) + SubsystemAtoms(2) + 1
                              p1 = NAtoms
                              Charge = SubsystemCharges(3)
                        case (SYS_DIMER_AB)
                              p0 = 1
                              p1 = SubsystemAtoms(1) + SubsystemAtoms(2)
                              Charge = SubsystemCharges(1) + SubsystemCharges(2)
                        case (SYS_DIMER_BC)
                              p0 = SubsystemAtoms(1) + 1
                              p1 = NAtoms
                              Charge = SubsystemCharges(2) + SubsystemCharges(3)
                        case (SYS_DIMER_AC)
                              p0 = 1
                              p1 = SubsystemAtoms(1)
                              q0 = SubsystemAtoms(1) + SubsystemAtoms(2) + 1
                              q1 = NAtoms
                              Charge = SubsystemCharges(1) + SubsystemCharges(3)
                        end select
                  case (SYS_TETRAMER)
                        Charge = sum(SubsystemCharges(1:4))
                        a0 = 1
                        a1 = SubsystemAtoms(1)
                        b0 = SubsystemAtoms(1) + 1
                        b1 = SubsystemAtoms(1) + SubsystemAtoms(2)
                        c0 = SubsystemAtoms(1) + SubsystemAtoms(2) + 1
                        c1 = SubsystemAtoms(1) + SubsystemAtoms(2) + SubsystemAtoms(3)
                        d0 = SubsystemAtoms(1) + SubsystemAtoms(2) + SubsystemAtoms(3) + 1
                        d1 = SubsystemAtoms(1) + SubsystemAtoms(2) + SubsystemAtoms(3) + SubsystemAtoms(4)
                        select case (i)
                        case (SYS_MONO_A)
                              p0 = a0
                              p1 = a1
                              Charge = SubsystemCharges(1)
                        case (SYS_MONO_B)
                              p0 = b0
                              p1 = b1
                              Charge = SubsystemCharges(2)
                        case (SYS_MONO_C)
                              p0 = c0
                              p1 = c1
                              Charge = SubsystemCharges(3)
                        case (SYS_DIMER_AB)
                              p0 = a0
                              p1 = b1
                              Charge = SubsystemCharges(1) + SubsystemCharges(2)
                        case (SYS_DIMER_BC)
                              p0 = b0
                              p1 = c1
                              Charge = SubsystemCharges(2) + SubsystemCharges(3)
                        case (SYS_DIMER_AC)
                              p0 = a0
                              p1 = a1
                              q0 = c0
                              q1 = c1
                              Charge = SubsystemCharges(1) + SubsystemCharges(3)
                        case (SYS_MONO_D)
                              p0 = d0
                              p1 = d1
                              Charge = SubsystemCharges(4)
                        case (SYS_DIMER_AD)
                              p0 = a0
                              p1 = a1
                              q0 = d0
                              q1 = d1
                              Charge = SubsystemCharges(1) + SubsystemCharges(4)
                        case (SYS_DIMER_BD)
                              p0 = b0
                              p1 = b1
                              q0 = d0
                              q1 = d1
                              Charge = SubsystemCharges(2) + SubsystemCharges(4)
                        case (SYS_DIMER_CD)
                              p0 = c0
                              p1 = d1
                              Charge = SubsystemCharges(3) + SubsystemCharges(4)
                        case (SYS_TRIMER_ABC)
                              p0 = a0
                              p1 = c1
                              Charge = SubsystemCharges(1) + SubsystemCharges(2) + SubsystemCharges(3)
                        case (SYS_TRIMER_ABD)
                              p0 = a0
                              p1 = b1
                              q0 = d0
                              q1 = d1
                              Charge = SubsystemCharges(1) + SubsystemCharges(2) + SubsystemCharges(4)
                        case (SYS_TRIMER_ACD)
                              p0 = a0
                              p1 = a1
                              q0 = c0
                              q1 = d1
                              Charge = SubsystemCharges(1) + SubsystemCharges(3) + SubsystemCharges(4)
                        case (SYS_TRIMER_BCD)
                              p0 = b0
                              p1 = d1
                              Charge = SubsystemCharges(2) + SubsystemCharges(3) + SubsystemCharges(4)
                        end select
                  end select
                  RealAtoms(1, 1) = p0
                  RealAtoms(2, 1) = p1
                  RealAtoms(1, 2) = q0
                  RealAtoms(2, 2) = q1
                  NElectrons = 0
                  do s = 1, 2
                        if (RealAtoms(2, s) >= RealAtoms(1, s)) then
                              p0 = RealAtoms(1, s)
                              p1 = RealAtoms(2, s)
                              if (ECPCharges) then
                                    NElectrons = NElectrons + sum(ZNumbersECP(p0:p1))
                              else
                                    NElectrons = NElectrons + sum(ZNumbers(p0:p1))
                              end if
                        end if
                  end do
                  NElectrons = NElectrons - Charge
                  Mult = SubsystemMult(i)
            end associate
      end subroutine sys_Init

      
      subroutine sys_ElementsList(ZList, ZCount, AtomElementMap, NElements, System, AtomsType)
            integer, dimension(:), allocatable, intent(out) :: ZList
            integer, dimension(:), allocatable, intent(out) :: ZCount
            integer, dimension(:), intent(out)              :: AtomElementMap
            integer, intent(out)                            :: NElements
            type(TSystem), intent(in)                       :: System
            integer, intent(in)                             :: AtomsType

            integer, dimension(KNOWN_ELEMENTS) :: AllElementsCount
            integer, dimension(KNOWN_ELEMENTS) :: RealElementsCount
            integer, dimension(KNOWN_ELEMENTS) :: GhostElementsCount
            integer, dimension(KNOWN_ELEMENTS) :: ZElementMap
            integer :: z, k, l, s

            if (AtomsType == SYS_ALL_ATOMS .or. AtomsType == SYS_GHOST_ATOMS) then
                  AllElementsCount = 0
                  do k = 1, System%NAtoms
                        z = System%ZNumbers(k)
                        AllElementsCount(z) = AllElementsCount(z) + 1
                  end do
            end if
            
            if (AtomsType == SYS_REAL_ATOMS .or. AtomsType == SYS_GHOST_ATOMS) then
                  RealElementsCount = 0
                  do s = 1, 2
                        do k = System%RealAtoms(1, s), System%RealAtoms(2, s)
                              z = System%ZNumbers(k)
                              RealElementsCount(z) = RealElementsCount(z) + 1
                        end do
                  end do
            end if
            
            if (AtomsType == SYS_GHOST_ATOMS) then
                  GhostElementsCount = AllElementsCount - RealElementsCount
            end if
            
            if (AtomsType == SYS_REAL_ATOMS) then
                  AllElementsCount = RealElementsCount
            else if (AtomsType == SYS_GHOST_ATOMS) then
                  AllElementsCount = GhostElementsCount
            end if
            
            NElements = 0
            do z = 1, KNOWN_ELEMENTS
                  if (AllElementsCount(z) > 0) NElements = NElements + 1
            end do
            allocate(ZList(NElements))
            allocate(ZCount(NElements))
            if (NElements > 0) then
                  l = 1
                  do z = 1, KNOWN_ELEMENTS
                        if (AllElementsCount(z) > 0) then
                              ZList(l) = z
                              if (AtomsType == SYS_ALL_ATOMS) then
                                    ZElementMap(z) = l
                              end if
                              ZCount(l) = AllElementsCount(z)
                              l = l + 1
                        end if
                  end do
                  if (AtomsType == SYS_ALL_ATOMS) then
                        do k = 1, System%NAtoms
                              AtomElementMap(k) = ZElementMap(System%ZNumbers(k))
                        end do
                  else
                        AtomElementMap = 0
                  end if
            end if
      end subroutine sys_ElementsList


      function sys_ChemicalFormula(System)
            !
            ! Compute chemical formula of a given molecule.
            !
            character(:), allocatable         :: sys_ChemicalFormula
            type(TSystem), intent(in)         :: System

            integer, dimension(:), allocatable :: ZListReal, ZListGhost
            integer, dimension(:), allocatable :: ZCountReal, ZCountGhost
            integer :: NElementsReal, NElementsGhost
            integer, dimension(KNOWN_ELEMENTS) :: AtomElementMap
            
            character(:), allocatable :: s1, s2
            integer :: k

            s1 = ""
            s2 = ""
            if (System%SystemKind == SYS_MOLECULE .or. System%SubsystemKind == SYS_TOTAL) then
                  call sys_ElementsList(ZListReal, ZCountReal, AtomElementMap, NElementsReal, System, SYS_ALL_ATOMS)
                  do k = 1, NElementsReal
                        s1 = s1 // trim(ELNAME_SHORT(ZListReal(k)))
                        if (ZCountReal(k) > 1) then
                              s1 = s1 // str(ZCountReal(k))
                        end if
                  end do
            else
                  select case (System%SubsystemKind)
                  case (SYS_MONO_A)     ! --- 1 ---
                        s1 = "Monomer A: "
                  case (SYS_MONO_B)     ! --- 2 ---
                        s1 = "Monomer B: "
                  case (SYS_MONO_C)     ! --- 3 ---
                        s1 = "Monomer C: "
                  case (SYS_DIMER_AB)   ! --- 4 ---
                        s1 = "Dimer AB: "
                  case (SYS_DIMER_BC)   ! --- 5 ---
                        s1 = "Dimer BC: "
                  case (SYS_DIMER_AC)   ! --- 6 ---
                        s1 = "Dimer AC: "
                  case (SYS_MONO_D)     ! --- 7 ---
                        s1 = "Monomer D: "
                  case (SYS_DIMER_AD)   ! --- 8 ---
                        s1 = "Dimer AD: "
                  case (SYS_DIMER_BD)   ! --- 9 ---
                        s1 = "Dimer BD: "
                  case (SYS_DIMER_CD)   ! --- 10 ---
                        s1 = "Dimer CD: "
                  case (SYS_TRIMER_ABC) ! --- 11 ---
                        s1 = "Trimer ABC: "
                  case (SYS_TRIMER_ABD) ! --- 12 ---
                        s1 = "Trimer ABD: "
                  case (SYS_TRIMER_ACD) ! --- 13 ---
                        s1 = "Trimer ACD: "
                  case (SYS_TRIMER_BCD) ! --- 14 ---
                        s1 = "Trimer BCD: "
                  end select
                  call sys_ElementsList(ZListReal, ZCountReal, AtomElementMap, NElementsReal, System, SYS_REAL_ATOMS)
                  call sys_ElementsList(ZListGhost, ZCountGhost, AtomElementMap, NElementsGhost, System, SYS_GHOST_ATOMS)
                  do k = 1, NElementsReal
                        s1 = s1 // trim(ELNAME_SHORT(ZListReal(k)))
                        if (ZCountReal(k) > 1) then
                              s1 = s1 // str(ZCountReal(k))
                        end if
                  end do
                  do k = 1, NElementsGhost
                        s2 = s2 // trim(ELNAME_SHORT(ZListGhost(k)))
                        if (ZCountGhost(k) > 1) then
                              s2 = s2 // str(ZCountGhost(k))
                        end if
                  end do
            end if
            if (System%Charge > 0 .and. System%Mult > 1) then
                  s1 = s1 // " (charge=" // str(System%Charge) // ", 2S+1=" // str(System%Mult) // ")"
            else if (System%Charge > 0) then
                  s1 = s1 // " (charge=" // str(System%Charge) // ")"
            else if (System%Mult > 1) then
                  s1 = s1 // " (2S+1=" // str(System%Mult) // ")"
            end if
            if (len(s2) > 0) then
                  sys_ChemicalFormula = s1 // " + ghost centers " // s2
            else 
                  sys_ChemicalFormula = s1
            end if
      end function sys_ChemicalFormula


      function sys_IsDummyAtom(System, k)
            logical :: sys_IsDummyAtom
            type(TSystem), intent(in) :: System
            integer, intent(in)       :: k

            sys_IsDummyAtom = .true.
            if (k >= System%RealAtoms(1, 1) .and. k <= System%RealAtoms(2, 1)) sys_IsDummyAtom = .false.
            if (k >= System%RealAtoms(1, 2) .and. k <= System%RealAtoms(2, 2)) sys_IsDummyAtom = .false.
      end function sys_IsDummyAtom


      subroutine sys_SortDistances(System)
            type(TSystem), intent(inout) :: System

            real(F64), dimension(3) :: ri, rj, d
            real(F64) :: dr
            integer :: i, j
            
            if (allocated(System%SortedDistances)) deallocate(System%SortedDistances)
            if (allocated(System%SortedDistancesIdx)) deallocate(System%SortedDistancesIdx)
            allocate(System%SortedDistances(System%NAtoms, System%NAtoms))
            allocate(System%SortedDistancesIdx(System%NAtoms, System%NAtoms))

            associate ( &
                  NAtoms => System%NAtoms, &
                  AtomCoords => System%AtomCoords, &
                  SortedDistances => System%SortedDistances, &
                  SortedDistancesIdx => System%SortedDistancesIdx &
                  )
                  if (NAtoms > 1) then
                        do i = 1, NAtoms
                              ri = AtomCoords(:, i)
                              SortedDistances(i, i) = ZERO
                              do j = 1, i - 1
                                    rj = AtomCoords(:, j)
                                    d = ri - rj
                                    dr = norm2(d)
                                    SortedDistances(i, j) = dr
                                    SortedDistances(j, i) = dr
                              end do
                        end do
                        !
                        ! For every nucleus J sort all
                        ! remaining nuclei according
                        ! to SORTEDDISTANCES(I, J)
                        !
                        do j = 1, NAtoms
                              do i = 1, NAtoms
                                    SortedDistancesIdx(i, j) = i
                              end do
                              call dsort(SortedDistances(:, j), SortedDistancesIdx(:, j), NAtoms)
                        end do
                  else
                        SortedDistances(1, 1) = ZERO
                        SortedDistancesIdx(1, 1) = 1
                  end if
            end associate
      end subroutine sys_SortDistances


      subroutine sys_NuclearRepulsion(Enucl, System)
            real(F64), intent(out)    :: Enucl
            type(TSystem), intent(in) :: System

            integer :: i, j, s, t
            real(F64), dimension(3) :: Ra, Rb, Rab
            real(F64) :: Dab, Qa, Qb

            associate ( &
                  AtomCoords => System%AtomCoords, &
                  RealAtoms => System%RealAtoms, &
                  ECPCharges => System%ECPCharges, &
                  ZNumbers => System%ZNumbers, &
                  ZNumbersECP => System%ZNumbersECP &                  
                  )
                  Enucl = ZERO
                  do s = 1, 2
                        do i = RealAtoms(1, s), RealAtoms(2, s)
                              Ra = AtomCoords(:, i)
                              if (ECPCharges) then
                                    Qa = real(ZNumbersECP(i), F64)
                              else
                                    Qa = real(ZNumbers(i), F64)
                              end if
                              do t = 1, 2
                                    do j = max(i+1, RealAtoms(1, t)), RealAtoms(2, t)
                                          Rb = AtomCoords(:, j)
                                          if (ECPCharges) then
                                                Qb = real(ZNumbersECP(j), F64)
                                          else
                                                Qb = real(ZNumbers(j), F64)
                                          end if
                                          Rab = Ra - Rb
                                          Dab = norm2(Rab)
                                          Enucl = Enucl + Qa * Qb / Dab
                                    end do
                              end do
                        end do
                  end do
            end associate
      end subroutine sys_NuclearRepulsion


      subroutine sys_NuclearMultipoles(Dx, Dy, Dz, Qyx, Qzx, &
            Qzy, Qxx, Qyy, Qzz, Rc, System)
            
            real(F64), intent(out)    :: Dx, Dy, Dz
            real(F64), intent(out)    :: Qyx, Qzx, Qzy, Qxx, Qyy, Qzz
            real(F64), dimension(3)   :: Rc
            type(TSystem), intent(in) :: System

            integer :: s, p0, p1, p
            integer :: Z
            real(F64), dimension(3) :: Rp, Rpc

            Dx = ZERO
            Dy = ZERO
            Dz = ZERO
            Qyx = ZERO
            Qzx = ZERO
            Qzy = ZERO
            Qxx = ZERO
            Qyy = ZERO
            Qzz = ZERO
            do s = 1, 2
                  if (System%RealAtoms(2, s) >= System%RealAtoms(1, s)) then
                        p0 = System%RealAtoms(1, s)
                        p1 = System%RealAtoms(2, s)
                        do p = p0, p1
                              if (System%ECPCharges) then
                                    Z = System%ZNumbersECP(p) 
                              else
                                    Z = System%ZNumbers(p)
                              end if
                              Rp(:) = System%AtomCoords(:, p)
                              Rpc(:) = Rp(:) - Rc(:)
                              Dx = Dx + Rpc(1) * Z
                              Dy = Dy + Rpc(2) * Z
                              Dz = Dz + Rpc(3) * Z
                              Qyx = Qyx + Rpc(1)*Rpc(2) * Z
                              Qzx = Qzx + Rpc(1)*Rpc(3) * Z
                              Qzy = Qzy + Rpc(2)*Rpc(3) * Z
                              Qxx = Qxx + Rpc(1)**2 * Z
                              Qyy = Qyy + Rpc(2)**2 * Z
                              Qzz = Qzz + Rpc(3)**2 * Z
                        end do
                  end if
            end do
      end subroutine sys_NuclearMultipoles


      subroutine sys_ChargeCenter(Rc, System)
            real(F64), dimension(3), intent(out) :: Rc
            type(TSystem), intent(in)            :: System

            integer :: SumZ, Z
            integer :: s, p, p0, p1

            Rc = ZERO
            SumZ = System%NElectrons + System%Charge
            do s = 1, 2
                  if (System%RealAtoms(2, s) >= System%RealAtoms(1, s)) then
                        p0 = System%RealAtoms(1, s)
                        p1 = System%RealAtoms(2, s)
                        do p = p0, p1
                              if (System%ECPCharges) then
                                    Z = System%ZNumbersECP(p) 
                              else
                                    Z = System%ZNumbers(p)
                              end if
                              Rc(:) = Rc(:) + real(Z, F64)/SumZ * System%AtomCoords(:, p)
                        end do
                  end if
            end do
      end subroutine sys_ChargeCenter

      
      subroutine sys_Read_XYZ(System, FilePath, Units)
            type(TSystem), intent(out)    :: System
            character(*), intent(in)      :: FilePath
            integer, optional, intent(in) :: Units

            logical :: XYZDefined, ReadingXYZBlock
            integer :: AtomIdx
            integer :: u
            character(:), allocatable :: key, val
            character(:), allocatable :: line
            logical :: eof
            integer :: Units0

            if (present(Units)) then
                  Units0 = Units
            else
                  Units0 = SYS_UNITS_ANGSTROM
            end if
            
            u = io_text_open(FilePath, "OLD")
            XYZDefined = .false.
            ReadingXYZBlock = .false.
            AtomIdx = -1
            lines: do
                  call io_text_readline(line, u, eof)
                  if (eof) then
                        exit lines
                  end if                  
                  if (isblank(line) .or. iscomment(line)) then
                        cycle lines
                  end if
                  call split(line, key, val)
                  key = uppercase(key)
                  if (key == "XYZ") then
                        XYZDefined = .true.
                        ReadingXYZBlock = .true.
                        cycle lines
                  else if (key == "END") then
                        if (ReadingXYZBlock) then
                              call sys_Init(System, SYS_TOTAL)
                              ReadingXYZBlock = .false.
                              exit lines
                        else
                              cycle lines
                        end if
                  else
                        if (ReadingXYZBlock) then
                              call sys_Read_XYZ_NextLine(System, AtomIdx, line, Units0)
                        end if
                  end  if
            end do lines
            if (.not. XYZDefined) then
                  call msg("XYZ coordinates not defined in file " // FilePath, MSG_ERROR)
                  error stop
            end if
            call sys_SortDistances(System)
            close(u)
      end subroutine sys_Read_XYZ


      subroutine sys_Read_XYZ_NextLine(System, AtomIdx, line, Units)
            type(TSystem), intent(inout) :: System
            integer, intent(inout)       :: AtomIdx
            character(*), intent(in)     :: line
            integer, intent(in)          :: Units

            character(:), allocatable :: key, val
            character(:), allocatable :: element, coords
            integer :: k, z
            integer :: NSubsystems

            call split(line, key, val)
            key = uppercase(key)
            if (System%SystemKind == SYS_NONE) then
                  NSubsystems = IntListLength(line)
                  select case (NSubsystems)
                  case (1)
                        System%SystemKind = SYS_MOLECULE
                  case (2)
                        System%SystemKind = SYS_DIMER
                  case (3)
                        System%SystemKind = SYS_TRIMER
                  case (4)
                        System%SystemKind = SYS_TETRAMER
                  case default
                        call msg("First line of the XYZ block has an invalid format", MSG_ERROR)
                        error stop
                  end select
                  AtomIdx = 0
                  read(line, *) (System%SubsystemAtoms(k), k=1,System%SystemKind)
                  System%NAtoms = sum(System%SubsystemAtoms)
                  System%RealAtoms(:, 1) = [1, System%NAtoms]
                  System%RealAtoms(:, 2) = [1, 0]
                  allocate(System%AtomCoords(3, System%NAtoms))
                  allocate(System%ZNumbers(System%NAtoms))
            else if (key == "CHARGE" .or. key == "CHARGES") then
                  read(val, *) (System%SubsystemCharges(k), k=1,System%SystemKind)
                  System%Charge = sum(System%SubsystemCharges)
            else if (key == "MULT" .or. key == "MULTIPLICITY") then
                  if (System%SystemKind==SYS_MOLECULE) then
                        read(val, *) System%SubsystemMult(1)
                  else if (System%SystemKind==SYS_DIMER) then
                        read(val, *) (System%SubsystemMult(k), k=1,3)
                  else if (System%SystemKind==SYS_TRIMER) then
                        read(val, *) (System%SubsystemMult(k), k=1,7)
                  else ! Multiplicities of subsystems in a tetramer
                        read(val, *) (System%SubsystemMult(k), k=1,15)
                  end if
                  System%Mult = System%SubsystemMult(1)
            else
                  if (AtomIdx > -1) then
                        AtomIdx = AtomIdx + 1
                        if (AtomIdx <= System%NAtoms) then
                              element = key
                              coords = val
                              z = znumber_short(element)
                              if (z > 0) then
                                    System%ZNumbers(AtomIdx) = z
                                    read(coords, *) (System%AtomCoords(k, AtomIdx), k=1,3)
                                    if (Units == SYS_UNITS_ANGSTROM) then
                                          System%AtomCoords(:, AtomIdx) = tobohr(System%AtomCoords(:, AtomIdx))
                                    end if
                              else
                                    call msg("Unknown element: " // element, MSG_ERROR)
                                    error stop
                              end if
                        else
                              call msg("Inconsistent number of atoms specified", MSG_ERROR)
                              error stop
                        end if
                  else
                        call msg("Unspecified number of atoms", MSG_ERROR)
                        error stop
                  end if
            end if
      end subroutine sys_Read_XYZ_NextLine
end module sys_definitions
