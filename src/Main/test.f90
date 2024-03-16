program test
      use sorter_Cholesky
      use Cholesky
      use basis_sets
      use sys_definitions
      use chol_definitions
      use CholeskyOTF, only: chol_CoulombMatrix_OTF, TCholeskyVecsOTF
      use arithmetic
      use display
      use Cholesky_driver
      use string
      use Auto2e
      use Auto2eInterface
      use OneElectronInts
      use sphergto
      use BeckeGrid
      use GridFunctions
      use grid_definitions
      
      implicit none

      character(:), allocatable :: RawIntegralsPath, BasisSetPath, XYZPath, NaturalOrbitalsPath
      character(:), allocatable :: Example, SBinaryFilePath, HFOrbitalsPath
      integer :: AOSource
      logical :: SpherAO
      integer :: Accuracy
      integer :: ExternalOrdering
      integer :: Units
      logical :: SortAngularMomenta

      Example = "Dalton/water-cc-pVTZ"
      RawIntegralsPath = "./examples/" // Example // "/results/AOTWOINT"
      BasisSetPath = "./examples/" // Example // "/basis.txt"
      XYZPath = "./examples/" // Example // "/molecule.xyz"
      SBinaryFilePath = "./examples/" // Example // "/S.bin"
      NaturalOrbitalsPath = "./examples/" // Example // "/results/CAONO.mol"
      HFOrbitalsPath = "./examples/" // Example // "/results/CAOMO.bin"
      !AOSource = 2 ! Molpro binary file
      AOSource = 1 ! Dalton binary file
      SpherAO = .true.
      Accuracy = CHOL_ACCURACY_DEFAULT
      ExternalOrdering = ORBITAL_ORDERING_DALTON
      Units = SYS_UNITS_BOHR
      if (ExternalOrdering == ORBITAL_ORDERING_ORCA) then
            SortAngularMomenta = .false.
      else if (ExternalOrdering == ORBITAL_ORDERING_MOLPRO) then
            SortAngularMomenta = .true.
      else if (ExternalOrdering == ORBITAL_ORDERING_DALTON) then
            SortAngularMomenta = .false.
      else
            call msg("Check angular momenta ordering")
      end if

!      call test_molcas(XYZPath, Units, BasisSetPath)
      


      
      ! call test_Grid(BasisSetPath, XYZPath, HFOrbitalsPath, &
      !       AOSource, SpherAO, ExternalOrdering, SortAngularMomenta, Units)
      
      
      !call test_Fock_Molpro(BasisSetPath, XYZPath, NaturalOrbitalsPath, SpherAO, Accuracy)

      if (ExternalOrdering == ORBITAL_ORDERING_DALTON) then
            call test_Dalton(XYZPath, Units, BasisSetPath)
            call test_Transform_OTF(RawIntegralsPath, BasisSetPath, XYZPath, HFOrbitalsPath, &
                  AOSource, SpherAO, Accuracy, ExternalOrdering, SortAngularMomenta, Units)
      else if (ExternalOrdering == ORBITAL_ORDERING_MOLPRO) then
            call test_Transform_OTF(RawIntegralsPath, BasisSetPath, XYZPath, NaturalOrbitalsPath, &
                  AOSource, SpherAO, Accuracy, ExternalOrdering, SortAngularMomenta, Units)
      end if

  !    if (ExternalOrdering == ORBITAL_ORDERING_ORCA) then
  !          call msg("-- Testing Orca overlap integral ---")
  !          call test_orca(XYZPath, Units, BasisSetPath, SBinaryFilePath)
  !    end if
            
contains

      subroutine test_orca(XYZPath, Units, BasisSetPath, SBinaryFilePath)
            character(*), intent(in) :: XYZPath
            integer, intent(in)      :: Units
            character(*), intent(in) :: BasisSetPath
            character(*), intent(in) :: SBinaryFilePath

            type(TSystem) :: System
            type(TAOBasis) :: AOBasis
            real(F64), dimension(:, :), allocatable :: S_cao, S_sao, S_extao, S_extao_2
            real(F64), dimension(:), allocatable :: TransfWork
            integer :: NAOSpher
            integer :: NAOCart
            integer :: NAO
            logical, parameter :: SpherAO = .true.
            logical, parameter :: SortAngularMomenta = .false.
            
            
            call auto2e_init()
            call sys_Read_XYZ(System, XYZPath, Units)
            call basis_NewAOBasis(AOBasis, System, BasisSetPath, SpherAO, SortAngularMomenta)

            NAOSpher = AOBasis%NAOSpher
            NAOCart = AOBasis%NAOCart
            NAO = NAOSpher
            
            allocate(S_cao(NAOCart, NAOCart))
            allocate(S_sao(NAO, NAO))
            allocate(S_extao(NAO, NAO))
            allocate(S_extao_2(NAO, NAO))

            call ints1e_OverlapMatrix(S_cao, AOBasis)            
            call linalg_smfill(S_cao)

            allocate(TransfWork(NAOSpher*NAOCart))
            call SpherGTO_TransformMatrix_U(S_sao, S_cao, &
                  AOBasis%LmaxGTO, &
                  AOBasis%NormFactorsSpher, &
                  AOBasis%NormFactorsCart, &
                  AOBasis%ShellLocSpher, &
                  AOBasis%ShellLocCart, &
                  AOBasis%ShellMomentum, &
                  AOBasis%ShellParamsIdx, &
                  AOBasis%NAOSpher, &
                  AOBasis%NAOCart, &
                  AOBasis%NShells, TransfWork)

            call auto2e_interface_ApplyOrcaPhases_Matrix(S_sao, AOBasis, .true.)
            call auto2e_interface_AngFuncTransf(S_extao_2, S_sao, .false., .true., AOBasis, ORBITAL_ORDERING_ORCA)
            call orca_read_symmetric_matrix(S_extao, SBinaryFilePath)

            call msg("--- Auto2e overlap matrix (Auto2e-Orca) ---")
            call geprn(S_extao_2)
            call msg("--- end of Auto2e overlap matrix ---")

            call msg("--- Orca overlap matrix ---")
            call geprn(S_extao)
            call msg("--- end of Orca overlap matrix ---")
      end subroutine test_orca


      subroutine test_molcas(XYZPath, Units, BasisSetPath)
            character(*), intent(in) :: XYZPath
            integer, intent(in)      :: Units
            character(*), intent(in) :: BasisSetPath

            type(TSystem) :: System
            type(TAOBasis) :: AOBasis
            real(F64), dimension(:, :), allocatable :: S_cao, S_sao, S_extao, S_extao_2
            real(F64), dimension(:), allocatable :: TransfWork
            integer :: NAOSpher
            integer :: NAOCart
            integer :: NAO
            logical, parameter :: SpherAO = .true.
            logical, parameter :: SortAngularMomenta = .false.
            
            
            call auto2e_init()
            call sys_Read_XYZ(System, XYZPath, Units)
            call basis_NewAOBasis(AOBasis, System, BasisSetPath, SpherAO, SortAngularMomenta)

            NAOSpher = AOBasis%NAOSpher
            NAOCart = AOBasis%NAOCart
            NAO = NAOSpher
            
            allocate(S_cao(NAOCart, NAOCart))
            allocate(S_sao(NAO, NAO))
            allocate(S_extao(NAO, NAO))
            allocate(S_extao_2(NAO, NAO))

            call ints1e_OverlapMatrix(S_cao, AOBasis)            
            call linalg_smfill(S_cao)

            allocate(TransfWork(NAOSpher*NAOCart))
            call SpherGTO_TransformMatrix_U(S_sao, S_cao, &
                  AOBasis%LmaxGTO, &
                  AOBasis%NormFactorsSpher, &
                  AOBasis%NormFactorsCart, &
                  AOBasis%ShellLocSpher, &
                  AOBasis%ShellLocCart, &
                  AOBasis%ShellMomentum, &
                  AOBasis%ShellParamsIdx, &
                  AOBasis%NAOSpher, &
                  AOBasis%NAOCart, &
                  AOBasis%NShells, TransfWork)

            call auto2e_interface_AngFuncTransf(S_extao_2, S_sao, .false., .true., AOBasis, ORBITAL_ORDERING_OPEN_MOLCAS)
            !            call orca_read_symmetric_matrix(S_extao, SBinaryFilePath)
!            S_extao_2 = S_sao

            call msg("--- Auto2e overlap matrix (Auto2e-Orca) ---")
            call geprn(S_extao_2)
            call msg("--- end of Auto2e overlap matrix ---")

            ! call msg("--- Orca overlap matrix ---")
            ! call geprn(S_extao)
            ! call msg("--- end of Orca overlap matrix ---")
      end subroutine test_molcas


      subroutine orca_read_symmetric_matrix(S, BinaryFilePath)
            real(F64), dimension(:, :), intent(out) :: S
            character(*), intent(in)                :: BinaryFilePath

            integer :: NAO
            integer :: irow, icol, ityp
            integer :: p, q, r
            real(F64), dimension(:), allocatable :: Work
            integer :: u
            
            NAO = size(S, dim=1)
            allocate(Work((NAO*(NAO+1))/2))
            
            open(newunit=u,File=BinaryFilePath,form='unformatted',access='stream', status='Old')
            read(u) IROW,ICOL,ITYP
            read(u) Work
            close(u)

            r = 1
            do p = 1, NAO
                  do q = 1, p
                        S(p, q) = Work(r)
                        r = r + 1
                  end do
            end do
            call linalg_smfill(S)
      end subroutine orca_read_symmetric_matrix


      
      

      subroutine test_Transform_OTF(RawIntegralsPath, BasisSetPath, XYZPath, NaturalOrbitalsPath, &
            AOSource, SpherAO, Accuracy, ExternalOrdering, SortAngularMomenta, Units)
            
            character(*), intent(in) :: RawIntegralsPath
            character(*), intent(in) :: BasisSetPath
            character(*), intent(in) :: XYZPath
            character(*), intent(in) :: NaturalOrbitalsPath
            integer, intent(in)      :: AOSource
            logical, intent(in)      :: SpherAO
            integer, intent(in)      :: Accuracy
            integer, intent(in)      :: ExternalOrdering
            logical, intent(in)      :: SortAngularMomenta
            integer, intent(in)      :: Units

            integer :: NAO, Na, Nb, a0, a1, b0, b1
            real(F64), dimension(:, :), allocatable :: CA, CB
            real(F64), dimension(:, :), allocatable :: Rkab, RkabOTF, RkabTHC
            integer :: iunit
            integer :: nrep, nbas(8)
            logical :: OnTheFly
            real(F64), dimension(:, :), allocatable :: CAONO, CSAONO
            integer, parameter :: MaxBufferDimMB = 10 ! buffer for mo transf, in megabytes
                        
            ! open(newunit=iunit,file=RawIntegralsPath,form='unformatted')
            ! read(iunit) nrep
            ! read(iunit) nbas(1:nrep)
            ! close(iunit)
            ! NAO = sum(nbas(1:nrep))

            open(newunit=iunit,file=NaturalOrbitalsPath,form='UNFORMATTED')
            read(iunit) NAO
            allocate(CAONO(1:NAO,1:NAO))
            allocate(CSAONO(1:NAO,1:NAO))
            read(iunit) CAONO(1:NAO,1:NAO)
            read(iunit) CSAONO(1:NAO,1:NAO)
            close(iunit)

            a0 = 1
            a1 = NAO
            b0 = 1
            b1 = NAO
            na = a1-a0+1
            nb = b1-b0+1

            OnTheFly = .true.
            call chol_MOVecs_v2(RkabOTF, CAONO, a0, a1, CAONO, b0, b1, MaxBufferDimMB, OnTheFly, &
                  RawIntegralsPath, AOSource, BasisSetPath, XYZPath, Accuracy, SpherAO, ExternalOrdering, &
                  SortAngularMomenta, Units)

            OnTheFly = .false.
            call chol_MOVecs_v2(Rkab, CAONO, a0, a1, CAONO, b0, b1, MaxBufferDimMB, OnTheFly, &
                  RawIntegralsPath, AOSource, BasisSetPath, XYZPath, Accuracy, SpherAO, ExternalOrdering, &
                  SortAngularMomenta, Units)

            call chol_MOVecs_THC(RkabTHC, CAONO, a0, a1, CAONO, b0, b1, BasisSetPath, XYZPath, Accuracy, SpherAO, &
                  ExternalOrdering, SortAngularMomenta, Units)

            call msg("------------- checking Cholesky binary vs Cholesky on the fly ----------------------")            
            call test_CheckRkab(Rkab, RkabOTF, NA, NB, size(Rkab,dim=1), size(RkabOTF,dim=1))
            call msg("------------- checking Cholesky binary vs Cholesky THC  ----------------------")
            call test_CheckRkab(Rkab, RkabTHC, NA, NB, size(Rkab,dim=1), size(RkabTHC,dim=1))

            call msg("------------- checking Cholesky OTF vs Cholesky THC  ----------------------")
            call test_CheckRkab(RkabOTF, RkabTHC, NA, NB, size(RkabOTF,dim=1), size(RkabTHC,dim=1))
      end subroutine test_Transform_OTF


      subroutine test_CheckRkab(Rkab, RkabOTF, NA, NB, NCholesky, NCholeskyOTF)
            real(F64), dimension(NCholesky, NA, NB), intent(in)    :: Rkab
            real(F64), dimension(NCholeskyOTF, NA, NB), intent(in) :: RkabOTF
            integer, intent(in)                                    :: NA, NB
            integer, intent(in)                                    :: NCholesky
            integer, intent(in)                                    :: NCholeskyOTF

            integer :: a, b, c, d
            real(F64) :: Vabcd, VabcdOTF
            real(F64) :: AbsErr, MaxAbsErr, RelErr, MaxRelErr
            real(F64) :: Val1, Val2, Val1OTF, Val2OTF

            call msg("Testing on the fly integrals")
            MaxAbsErr = ZERO
            MaxRelErr = ZERO
            do b = 1, NB
                  do a = 1, NA
                        c = a
                        d = b
                        Vabcd = dot_product(Rkab(:, a, b), Rkab(:, c, d))
                        VabcdOTF = dot_product(RkabOTF(:, a, b), RkabOTF(:, c, d))                                    
                        AbsErr = abs(Vabcd-VabcdOTF)
                        if (abs(Vabcd) > 1.0E-5_F64) then
                              RelErr = abs(Vabcd-VabcdOTF)/abs(Vabcd)
                        else
                              RelErr = ZERO
                        end if
                        if (AbsErr > MaxAbsErr) then
                              MaxAbsErr = AbsErr
                              Val1 = Vabcd
                              Val1OTF = VabcdOTF
                        end if
                        if (RelErr > MaxRelErr) then
                              MaxRelErr = RelErr
                              Val2 = Vabcd
                              Val2OTF = VabcdOTF
                        end if
                  end do
            end do
            call msg("Maximum absolute error: " // str(MaxAbsErr,d=1) // " at (ab|cd) = " // str(Val1,d=10))
            call msg("                               (ab|cd)Approx = " // str(Val1OTF,d=10))
            call msg("Maximum relative error: " // str(MaxRelErr,d=1) // " at (ab|cd) = " // str(Val2,d=10))
            call msg("                               (ab|cd)Approx = " // str(Val2OTF,d=10))
      end  subroutine test_CheckRkab



      subroutine test_Fock_Molpro(BasisSetPath, XYZPath, HFOrbitalsPath, SpherAO, Accuracy)            
            character(*), intent(in) :: BasisSetPath
            character(*), intent(in) :: XYZPath
            character(*), intent(in) :: HFOrbitalsPath
            logical, intent(in)      :: SpherAO
            integer, intent(in)      :: Accuracy

            integer :: NAO, NMO
            integer :: iunit
            integer :: NOcc
            integer :: p
            logical :: OnTheFly
            type(TSystem) :: System
            type(TAOBasis) :: AOBasis
            type(TCholeskyVecsOTF) :: CholeskyVecs
            real(F64), dimension(:, :), allocatable :: Rho_mo, F_mo, H0_mo
            real(F64), dimension(:), allocatable :: OrbEnergies
            integer, parameter :: Units = SYS_UNITS_ANGSTROM
            real(F64), dimension(:, :), allocatable :: CAONO, CSAONO
            integer, parameter :: MaxBufferDimMB = 10 ! buffer for mo transf, in megabytes
            logical, parameter :: SortAngularMomenta = .false.
            real(F64), parameter :: KScal = 1.0_F64
            integer, parameter :: ExternalOrdering = ORBITAL_ORDERING_MOLPRO
                        
            open(newunit=iunit,file=HFOrbitalsPath,form='UNFORMATTED')
            read(iunit) NAO
            allocate(CAONO(1:NAO,1:NAO))
            allocate(CSAONO(1:NAO,1:NAO))
            read(iunit) CAONO(1:NAO,1:NAO)
            read(iunit) CSAONO(1:NAO,1:NAO)
            close(iunit)

            call auto2e_init()
            call sys_Read_XYZ(System, XYZPath, Units)
            call basis_NewAOBasis(AOBasis, System, BasisSetPath, SpherAO, SortAngularMomenta)

            if ((SpherAO .and. NAO /= AOBasis%NAOSpher) .or. &
                  (.not. SpherAO .and. NAO /= AOBAsis%NAOCart)) then
                  call msg("Invalid number of orbitals")
                  error stop
            end if
            
            NMO = NAO
            NOcc = System%NElectrons/2
            allocate(Rho_mo(NMO, NMO))
            Rho_mo = ZERO
            do p = 1, NOcc
                  Rho_mo(p, p) = TWO
            end do
            allocate(F_mo(NMO, NMO))
            allocate(H0_mo(NMO, NMO))

            call chol_Rkpq_OTF(CholeskyVecs, AOBasis, Accuracy)

            call chol_F(F_mo, H0_mo, Rho_mo, CAONO, KScal, CholeskyVecs, AOBasis, &
                  System, ExternalOrdering, MaxBufferDimMB)

            print *, "------------------------- F_MO ----------------------------"            
            call geprn(F_mo)
            print *, "------------------------- H0_MO ---------------------------"
            call geprn(H0_mo)                       
      end subroutine test_Fock_Molpro

      
      subroutine test_Dalton(XYZPath, Units, BasisSetPath)
            character(*), intent(in) :: XYZPath
            integer, intent(in)      :: Units
            character(*), intent(in) :: BasisSetPath

            type(TSystem) :: System
            type(TAOBasis) :: AOBasis
            real(F64), dimension(:, :), allocatable :: S_cao, S_sao, S_extao, T_cao, T_sao, T_extao, V_cao, V_sao, V_extao
            real(F64), dimension(:), allocatable :: TransfWork
            integer :: NAOSpher
            integer :: NAOCart
            integer :: NAO
            logical, parameter :: SpherAO = .true.
            logical, parameter :: SortAngularMomenta = .false.
            
            
            call auto2e_init()
            !
            ! Initialize the Boys function interpolation table
            ! (used for Coulomb integrals evaluation).
            !
            call boys_init(4 * AUTO2E_MAXL)
            call sys_Read_XYZ(System, XYZPath, Units)
            call basis_NewAOBasis(AOBasis, System, BasisSetPath, SpherAO, SortAngularMomenta)

            NAOSpher = AOBasis%NAOSpher
            NAOCart = AOBasis%NAOCart
            NAO = NAOSpher
            allocate(S_cao(NAOCart, NAOCart))
            allocate(S_sao(NAO, NAO))
            allocate(S_extao(NAO, NAO))

            allocate(T_cao(NAOCart, NAOCart))
            allocate(T_sao(NAO, NAO))
            allocate(T_extao(NAO, NAO))

            allocate(V_cao(NAOCart, NAOCart))
            allocate(V_sao(NAO, NAO))
            allocate(V_extao(NAO, NAO))

            call ints1e_OverlapMatrix(S_cao, AOBasis)
            call ints1e_Coulomb(V_cao, AOBasis, System)
            call ints1e_Kinetic(T_cao, AOBasis)
            
            call linalg_smfill(S_cao)
            call linalg_smfill(T_cao)
            call linalg_smfill(V_cao)
            
            allocate(TransfWork(NAOSpher*NAOCart))
            call SpherGTO_TransformMatrix_U(S_sao, S_cao, &
                  AOBasis%LmaxGTO, &
                  AOBasis%NormFactorsSpher, &
                  AOBasis%NormFactorsCart, &
                  AOBasis%ShellLocSpher, &
                  AOBasis%ShellLocCart, &
                  AOBasis%ShellMomentum, &
                  AOBasis%ShellParamsIdx, &
                  AOBasis%NAOSpher, &
                  AOBasis%NAOCart, &
                  AOBasis%NShells, TransfWork)
            call SpherGTO_TransformMatrix_U(T_sao, T_cao, &
                  AOBasis%LmaxGTO, &
                  AOBasis%NormFactorsSpher, &
                  AOBasis%NormFactorsCart, &
                  AOBasis%ShellLocSpher, &
                  AOBasis%ShellLocCart, &
                  AOBasis%ShellMomentum, &
                  AOBasis%ShellParamsIdx, &
                  AOBasis%NAOSpher, &
                  AOBasis%NAOCart, &
                  AOBasis%NShells, TransfWork)
            call SpherGTO_TransformMatrix_U(V_sao, V_cao, &
                  AOBasis%LmaxGTO, &
                  AOBasis%NormFactorsSpher, &
                  AOBasis%NormFactorsCart, &
                  AOBasis%ShellLocSpher, &
                  AOBasis%ShellLocCart, &
                  AOBasis%ShellMomentum, &
                  AOBasis%ShellParamsIdx, &
                  AOBasis%NAOSpher, &
                  AOBasis%NAOCart, &
                  AOBasis%NShells, TransfWork)

            
            call auto2e_interface_AngFuncTransf(S_extao, S_sao, .false., .true., AOBasis, ORBITAL_ORDERING_DALTON)
            call auto2e_interface_AngFuncTransf(T_extao, T_sao, .false., .true., AOBasis, ORBITAL_ORDERING_DALTON)
            call auto2e_interface_AngFuncTransf(V_extao, V_sao, .false., .true., AOBasis, ORBITAL_ORDERING_DALTON)

            call msg("--- Auto2e overlap matrix (Auto2e->Dalton) ---")
            call geprn(S_extao)
            call msg("--- end of Auto2e overlap matrix ---")
            
            call msg("--- Auto2e Vne matrix (Auto2e->Dalton) ---")
            call geprn(V_extao)
            call msg("--- end of Auto2e Vne matrix ---")


            call msg("--- Auto2e T matrix (Auto2e->Dalton) ---")
            call geprn(T_extao)
            call msg("--- end of Auto2e T matrix ---")

            call boys_free()
      end subroutine test_Dalton


      subroutine test_Grid(BasisSetPath, XYZPath, NaturalOrbitalsPath, &
            AOSource, SpherAO, ExternalOrdering, SortAngularMomenta, Units)
            
            character(*), intent(in) :: BasisSetPath
            character(*), intent(in) :: XYZPath
            character(*), intent(in) :: NaturalOrbitalsPath
            integer, intent(in)      :: AOSource
            logical, intent(in)      :: SpherAO
            integer, intent(in)      :: ExternalOrdering
            logical, intent(in)      :: SortAngularMomenta
            integer, intent(in)      :: Units

            integer :: NAO
            integer :: iunit
            real(F64), dimension(:, :), allocatable :: CAONO, CSAONO
            integer, parameter :: MaxBufferDimMB = 10 ! buffer for mo transf, in megabytes

            integer :: NPoints
            real(F64), dimension(:), allocatable :: Xg, Yg, Zg, Wg
            real(F64), dimension(:, :, :), allocatable :: Phi
            real(F64), dimension(:, :), allocatable :: C_ao
            integer, parameter :: GridType = BECKE_PARAMS_FINE
            integer :: NOcc
            real(F64) :: OccNumber
            real(F64) :: RhoIntegral, DivRhoIntegral
            type(TAOBasis) :: AOBasis
            type(TSystem) :: System
            real(F64), dimension(:), allocatable :: OccNumbers
            real(F64), dimension(:, :), allocatable :: Rho
            integer :: k
                        
            call auto2e_init()
            call sys_Read_XYZ(System, XYZPath, Units)
            call basis_NewAOBasis(AOBasis, System, BasisSetPath, SpherAO, SortAngularMomenta)
            if (AOBasis%SpherAO) then
                  NAO = AOBasis%NAOSpher
            else
                  NAO = AOBasis%NAOCart
            end if

            open(newunit=iunit,file=NaturalOrbitalsPath,form='UNFORMATTED')
            read(iunit) NAO
            allocate(CAONO(1:NAO,1:NAO))
            allocate(CSAONO(1:NAO,1:NAO))
            read(iunit) CAONO(1:NAO,1:NAO)
            read(iunit) CSAONO(1:NAO,1:NAO)
            close(iunit)
            !
            ! Molecular grid
            !
            call becke_MolecularGrid(Xg, Yg, Zg, Wg, NPoints, GridType, System, AOBasis)
            !
            ! Atomic orbitals and gradient components on the grid
            ! Phi(k, p, i)
            ! k grid point
            ! p index of AO
            ! i = 1 (value), 2 (d/dx), 3 (d/dy), 4 (d/dz)
            !
            allocate(Phi(NPoints, NAO, 4))
            call gridfunc_Orbitals_Grad(Phi, Xg, Yg, Zg, AOBasis)
            !
            ! MO coefficients (external format) -> MO coefficients (native format)
            !
            allocate(C_ao(NAO, NAO))
            call auto2e_interface_C(C_ao, CAONO, AOBasis, ExternalOrdering)
            !
            ! Occupation numbers (fake)
            !
            NOcc = System%NElectrons / 2
            allocate(OccNumbers(NAO))
            OccNumbers = ZERO
            OccNumbers(1:NOcc) = TWO
            !
            ! Rho and GradRho values on the grid
            ! Rho(k, i)
            ! k grid point
            ! i = 1 (rho value), 2 (dRho/dx), 3 (dRho/dy), 4 (dRho/dz)
            !
            allocate(Rho(NPoints, 4))
            call gridfunc_GGA_Variables(Rho, Phi, C_ao, OccNumbers)
            !
            ! Tests:
            ! * integral of Rho over the whole space (exact value: number of electrons)
            ! * integral of div(Rho) over the whole space (exact value: zero)
            !
            DivRhoIntegral = ZERO
            RhoIntegral = ZERO
            do k = 1, NPoints
                  RhoIntegral = RhoIntegral + Wg(k) * Rho(k, 1)
                  DivRhoIntegral = DivRhoIntegral + Wg(k) * (Rho(k, 2) + Rho(k, 3) + Rho(k, 4))
            end do
            call msg("Rho integral    " // str(RhoIntegral,d=6))
            call msg("DivRho integral " // str(DivRhoIntegral,d=6))
      end subroutine test_Grid
end program test
