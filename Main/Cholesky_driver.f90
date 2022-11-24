module Cholesky_driver
      use arithmetic
      use Cholesky, only: chol_CoulombMatrix, TCholeskyVecs, chol_Rkab_ExternalBinary, chol_MOTransf_TwoStep
      use CholeskyOTF, only: chol_CoulombMatrix_OTF, TCholeskyVecsOTF, chol_MOTransf_TwoStep_OTF
      use CholeskyFock, only: chf_H0, chf_JK
      use OneElectronInts
      use boys
      use Auto2e
      use basis_sets
      use sys_definitions
      use chol_definitions
      use display
      use Auto2eInterface
      
      implicit none

contains

      subroutine chol_MOVecs(Rkab, CA, a0, a1, CB, b0, b1, MaxBufferDimMB, OnTheFly, &
            RawIntegralsPath, AOSource, BasisSetPath, XYZPath, Accuracy, SpherAO, ExternalOrdering, &
            SortAngularMomenta)

            real(F64), dimension(:, :), allocatable, intent(out) :: Rkab
            real(F64), dimension(:, :), intent(in)               :: CA
            integer, intent(in)                                  :: a0, a1
            real(F64), dimension(:, :), intent(in)               :: CB
            integer, intent(in)                                  :: b0, b1
            integer, intent(in)                                  :: MaxBufferDimMB
            logical, intent(in)                                  :: OnTheFly
            character(*), intent(in)                             :: RawIntegralsPath
            integer, intent(in)                                  :: AOSource
            character(*), intent(in)                             :: BasisSetPath
            character(*), intent(in)                             :: XYZPath
            integer, intent(in)                                  :: Accuracy
            logical, intent(in)                                  :: SpherAO
            integer, intent(in)                                  :: ExternalOrdering
            logical, intent(in)                                  :: SortAngularMomenta

            type(TCholeskyVecs) :: CholeskyVecs
            type(TCholeskyVecsOTF) :: CholeskyVecsOTF
            type(TAOBasis) :: AOBasis
            type(TSystem) :: System
            integer :: NCholesky, NA, NB
            real(F64), dimension(:, :), allocatable :: CA_ao, CB_ao
            logical :: TwoIndexTransf, FromExternalAO

            NA = a1 - a0 + 1
            NB = b1 - b0 + 1
            if (OnTheFly) then
                  ! ----------------------------------------------------------------
                  ! Cholesky with on the fly integrals calculation
                  ! Compatible with MO coeffs from the external programs defined
                  ! in the Auto2eInterface module.
                  ! ----------------------------------------------------------------
                  !
                  ! Initialize the molecular integrals library.
                  ! The init subroutine binds subroutines to subroutine pointers.
                  !
                  call auto2e_init()
                  !
                  ! Initialize the Boys function interpolation table
                  ! (used for Coulomb integrals evaluation).
                  !
                  call boys_init(4 * AUTO2E_MAXL)
                  !
                  ! Read the XYZ coordinates and atom types
                  !
                  call sys_Read_XYZ(System, XYZPath)
                  !
                  ! Read the basis set parameters from an EMSL text file
                  ! (GAMESS-US format, no need for any edits, just download it straight from the website)
                  !
                  call basis_NewAOBasis(AOBasis, System, BasisSetPath, SpherAO, SortAngularMomenta)
                  if (AOBasis%LmaxGTO > AUTO2E_MAXL) then
                        call msg("Basis set includes angular momenta unsupported by the Auto2e subroutine")
                        error stop
                  end if
                  if (SpherAO) then
                        if (size(CA, dim=1) /= AOBasis%NAOSpher .or. size(CB, dim=1) /= AOBasis%NAOSpher) then
                              call msg("Invalid number of atomic orbitals in CA/CB")
                              error stop
                        end if
                  else
                        if (size(CA, dim=1) /= AOBasis%NAOCart .or. size(CB, dim=1) /= AOBasis%NAOCart) then
                              call msg("Invalid number of atomic orbitals in CA/CB")
                              error stop
                        end if
                  end if
                  !
                  ! Change the ordering of angular functions from the external program
                  ! to the internal ordering of the Auto2e module.
                  ! The ordering of whole shells is assumed to be already the same
                  ! as in the external program.
                  !
                  allocate(CA_ao, mold=CA)
                  allocate(CB_ao, mold=CB)
                  FromExternalAO = .true. ! AOs in external format, e.g., Molpro, Orca -> AOs in the format of Auto2e
                  TwoIndexTransf = .false. ! Transform only the index p of C(p,k)
                  call auto2e_interface_AngFuncTransf(CA_ao, CA, FromExternalAO, TwoIndexTransf, AOBasis, ExternalOrdering)
                  call auto2e_interface_AngFuncTransf(CB_ao, CB, FromExternalAO, TwoIndexTransf, AOBasis, ExternalOrdering)
                  if (ExternalOrdering == ORBITAL_ORDERING_ORCA) then
                        call auto2e_interface_ApplyOrcaPhases_Matrix(CA_ao, AOBasis, TwoIndexTransf)
                        call auto2e_interface_ApplyOrcaPhases_Matrix(CB_ao, AOBasis, TwoIndexTransf)
                  end if
                  !
                  ! Compute the Cholesky vectors is AO basis
                  !
                  call chol_CoulombMatrix_OTF(CholeskyVecsOTF, Accuracy, AOBasis)
                  NCholesky = CholeskyVecsOTF%NVecs
                  !
                  ! Transform the Cholesky vectors to the MO basis
                  !
                  allocate(Rkab(NCholesky, NA*NB))
                  call chol_MOTransf_TwoStep_OTF(Rkab, CholeskyVecsOTF, CA_ao, CB_ao, a0, a1, b0, b1, &
                        MaxBufferDimMB, AOBasis)
                  !
                  ! Deallocate the Boys function interpolation tables
                  ! to avoid allocation of already allocated arrays
                  ! if this subroutine is called again
                  !
                  call boys_free()
            else
                  !
                  ! Cholesky with two-electron integrals read from the disk
                  !
                  block
                        integer :: iunit
                        integer :: NAO
                        integer :: nres, nrep, nbas(8)
                        
                        open(newunit=iunit,file=RawIntegralsPath,form='unformatted')
                        read(iunit) nrep
                        read(iunit) nbas(1:nrep)
                        close(iunit)
                        NAO = sum(nbas(1:nrep))
                        call chol_CoulombMatrix(CholeskyVecs, NAO, RawIntegralsPath, AOSource, Accuracy)
                        NCholesky = CholeskyVecs%NCholesky
                        allocate(Rkab(NCholesky, NA*NB))
                        call chol_MOTransf_TwoStep(Rkab, CholeskyVecs, CA, a0, a1, CB, b0, b1, 10)
                  end block
            end if
      end subroutine chol_MOVecs


      subroutine chol_MOVecs_v2(Rkab, CA, a0, a1, CB, b0, b1, MaxBufferDimMB, OnTheFly, &
            RawIntegralsPath, AOSource, BasisSetPath, XYZPath, Accuracy, SpherAO, ExternalOrdering, &
            SortAngularMomenta)

            real(F64), dimension(:, :), allocatable, intent(out) :: Rkab
            real(F64), dimension(:, :), intent(in)               :: CA
            integer, intent(in)                                  :: a0, a1
            real(F64), dimension(:, :), intent(in)               :: CB
            integer, intent(in)                                  :: b0, b1
            integer, intent(in)                                  :: MaxBufferDimMB
            logical, intent(in)                                  :: OnTheFly
            character(*), intent(in)                             :: RawIntegralsPath
            integer, intent(in)                                  :: AOSource
            character(*), intent(in)                             :: BasisSetPath
            character(*), intent(in)                             :: XYZPath
            integer, intent(in)                                  :: Accuracy
            logical, intent(in)                                  :: SpherAO
            integer, intent(in)                                  :: ExternalOrdering
            logical, intent(in)                                  :: SortAngularMomenta

            type(TCholeskyVecs) :: CholeskyVecs
            type(TCholeskyVecsOTF) :: CholeskyVecsOTF
            type(TAOBasis) :: AOBasis
            type(TSystem) :: System
            integer :: NCholesky, NA, NB

            NA = a1 - a0 + 1
            NB = b1 - b0 + 1
            if (OnTheFly) then
                  !
                  ! Read the XYZ coordinates and atom types
                  !
                  call sys_Read_XYZ(System, XYZPath)
                  !
                  ! Read the basis set parameters from an EMSL text file
                  ! (GAMESS-US format, no need for any edits, just download it straight from the website)
                  !
                  call basis_NewAOBasis(AOBasis, System, BasisSetPath, SpherAO, SortAngularMomenta)
                  !
                  ! Compute Cholesky vectors in AO basis
                  !
                  call chol_Rkpq_OTF(CholeskyVecsOTF, AOBasis, Accuracy)
                  !
                  ! Transform Cholesky vectors to MO basis
                  !
                  NCholesky = CholeskyVecsOTF%NVecs
                  allocate(Rkab(NCholesky,NA*NB))
                  call chol_Rkab_OTF(Rkab, CA, a0, a1, CB, b0, b1, MaxBufferDimMB, &
                        CholeskyVecsOTF, AOBasis, ExternalOrdering)
            else
                  !
                  ! Cholesky with two-electron integrals read from the disk
                  !
                  call chol_Rkpq_ExternalBinary(CholeskyVecs, RawIntegralsPath, AOSource, Accuracy)
                  !
                  ! Transform Cholesky vectors to MO basis                  
                  !
                  NCholesky = CholeskyVecs%NCholesky
                  allocate(Rkab(NCholesky,NA*NB))
                  call chol_Rkab_ExternalBinary(Rkab, CholeskyVecs, CA, a0, a1, CB, b0, b1, MaxBufferDimMB)
            end if
      end subroutine chol_MOVecs_v2
      

      subroutine chol_Rkpq_ExternalBinary(CholeskyVecs, RawIntegralsPath, AOSource, Accuracy)
            !
            ! Cholesky vectors in AO basis computed with two-electron integrals read from disk.
            !
            type(TCholeskyVecs), intent(out) :: CholeskyVecs
            character(*), intent(in)         :: RawIntegralsPath
            integer, intent(in)              :: AOSource
            integer, intent(in)              :: Accuracy

            integer :: iunit
            integer :: NAO
            integer :: nrep, nbas(8)
            
            open(newunit=iunit,file=RawIntegralsPath,form='unformatted')
            read(iunit) nrep
            read(iunit) nbas(1:nrep)
            close(iunit)
            NAO = sum(nbas(1:nrep))
            call chol_CoulombMatrix(CholeskyVecs, NAO, RawIntegralsPath, AOSource, Accuracy)
      end subroutine chol_Rkpq_ExternalBinary


      subroutine chol_Rkpq_OTF(CholeskyVecsOTF, AOBasis, Accuracy)
            !
            ! Cholesky vectors in AO basis computed on the fly.
            !
            type(TCholeskyVecsOTF), intent(out)                  :: CholeskyVecsOTF
            type(TAOBasis), intent(in)                           :: AOBasis
            integer, intent(in)                                  :: Accuracy
            !
            ! Initialize the molecular integrals library.
            ! The init subroutine binds subroutines to subroutine pointers.
            !
            call auto2e_init()
            !
            ! Initialize the Boys function interpolation table
            ! (used for Coulomb integrals evaluation).
            !
            call boys_init(4 * AUTO2E_MAXL)
            if (AOBasis%LmaxGTO > AUTO2E_MAXL) then
                  call msg("Basis set includes angular momenta unsupported by the Auto2e subroutine")
                  error stop
            end if
            !
            ! Compute the Cholesky vectors is AO basis
            !
            call chol_CoulombMatrix_OTF(CholeskyVecsOTF, Accuracy, AOBasis)
            !
            ! Deallocate the Boys function interpolation tables
            ! to avoid allocation of already allocated arrays
            ! if this subroutine is called again
            !
            call boys_free()
      end subroutine chol_Rkpq_OTF


      subroutine chol_Rkab_OTF(Rkab, CA, a0, a1, CB, b0, b1, MaxBufferDimMB, CholeskyVecsOTF, AOBasis, ExternalOrdering)
            !
            ! Transformation of Cholesky vectors computed in AO basis with the on the fly
            ! algorithm.
            !
            ! Compatible with MO coeffs generated with the programs listed in the Auto2eInterface
            ! module.
            !
            real(F64), dimension(:, :), intent(out) :: Rkab
            real(F64), dimension(:, :), intent(in)  :: CA
            integer, intent(in)                     :: a0, a1
            real(F64), dimension(:, :), intent(in)  :: CB
            integer, intent(in)                     :: b0, b1
            integer, intent(in)                     :: MaxBufferDimMB
            type(TCholeskyVecsOTF), intent(in)      :: CholeskyVecsOTF
            type(TAOBasis), intent(in)              :: AOBasis
            integer, intent(in)                     :: ExternalOrdering

            integer :: NCholesky, NA, NB
            real(F64), dimension(:, :), allocatable :: CA_ao, CB_ao
            logical :: TwoIndexTransf, FromExternalAO, SpherAO

            NA = a1 - a0 + 1
            NB = b1 - b0 + 1
            SpherAO = AOBasis%SpherAO
            NCholesky = CholeskyVecsOTF%NVecs
            if (SpherAO) then
                  if (size(CA, dim=1) /= AOBasis%NAOSpher .or. size(CB, dim=1) /= AOBasis%NAOSpher) then
                        call msg("Invalid number of atomic orbitals in CA/CB")
                        error stop
                  end if
            else
                  if (size(CA, dim=1) /= AOBasis%NAOCart .or. size(CB, dim=1) /= AOBasis%NAOCart) then
                        call msg("Invalid number of atomic orbitals in CA/CB")
                        error stop
                  end if
            end if
            if (size(Rkab,dim=1) /= NCholesky .or. size(Rkab,dim=2) /= NA*NB) then
                  call msg("Invalid dimension of Rkab")
                  error stop
            end if
            !
            ! Change the ordering of angular functions from Molpro
            ! to the internal ordering of the Auto2e module.
            ! The ordering of whole shells is assumed to be already
            ! the same as in the external program.
            !
            allocate(CA_ao, mold=CA)
            allocate(CB_ao, mold=CB)                  
            FromExternalAO = .true. ! AOs from external program -> AOs in the Auto2e format
            TwoIndexTransf = .false. ! Transform only the index p of C(p,k)
            call auto2e_interface_AngFuncTransf(CA_ao, CA, FromExternalAO, TwoIndexTransf, AOBasis, ExternalOrdering)
            call auto2e_interface_AngFuncTransf(CB_ao, CB, FromExternalAO, TwoIndexTransf, AOBasis, ExternalOrdering)
            if (ExternalOrdering == ORBITAL_ORDERING_ORCA) then
                  call auto2e_interface_ApplyOrcaPhases_Matrix(CA_ao, AOBasis, TwoIndexTransf)
                  call auto2e_interface_ApplyOrcaPhases_Matrix(CB_ao, AOBasis, TwoIndexTransf)
            end if
            !
            ! Transform the Cholesky vectors to the MO basis
            !
            call chol_MOTransf_TwoStep_OTF(Rkab, CholeskyVecsOTF, CA_ao, CB_ao, a0, a1, b0, b1, &
                  MaxBufferDimMB, AOBasis)
      end subroutine chol_Rkab_OTF


      subroutine chol_F(F_mo, H0_mo, Rho_mo, C, KScal, CholeskyVecs, AOBasis, &
            System, ExternalOrdering, MaxBufferDimMB, J_mo, K_mo)
            
            real(F64), dimension(:, :), intent(out) :: F_mo
            real(F64), dimension(:, :), intent(out) :: H0_mo
            real(F64), dimension(:, :), intent(in)  :: Rho_mo
            real(F64), dimension(:, :), intent(in)  :: C
            real(F64), intent(in)                   :: KScal
            type(TCholeskyVecsOTF), intent(in)      :: CholeskyVecs
            type(TAOBasis), intent(in)              :: AOBasis
            type(TSystem), intent(in)               :: System
            integer, intent(in)                     :: ExternalOrdering
            integer, intent(in)                     :: MaxBufferDimMB
            real(F64), dimension(:, :), optional, intent(out) :: J_mo
            real(F64), dimension(:, :), optional, intent(out) :: K_mo

            real(F64), dimension(:, :), allocatable :: C_ao
            real(F64), dimension(:, :), allocatable :: Ceig_mo
            real(F64), dimension(:, :, :), allocatable :: Ceig_ao
            real(F64), dimension(:, :), allocatable :: OccNum
            real(F64), dimension(:, :), allocatable :: H0_ao
            real(F64), dimension(:, :), allocatable :: TransfWork
            real(F64), dimension(:, :, :), allocatable :: Rho_ao, JK_ao, J_ao, K_ao
            integer :: NMO, NAO
            logical :: FromExternalAO, TwoIndexTransf
            integer :: p
            integer, dimension(2) :: NOcc
            logical :: CoulContrib, ExchContrib

            call auto2e_init()
            call boys_init(4*AUTO2E_MAXL)
            associate (NAOCart => AOBasis%NAOCart, &
                  NAOSpher => AOBasis%NAOSpher, &
                  SpherAO => AOBasis%SpherAO)

                  NMO = size(C, dim=2)
                  NAO = size(C, dim=1)
                  if ((SpherAO .and. NAO /= NAOSpher) .or. (.not.SpherAO .and. NAO /= NAOCart)) then
                        call msg("Invalid dimension of matrix C")
                        error stop
                  end if
                  if (NMO /= size(Rho_mo,dim=1)) then
                        call msg("Invalid dimension of matrix Rho_mo")
                        error stop
                  end if
                  !
                  ! Bare nuclei hamiltonian (atomic orbital basis, auto2e ordering)
                  !
                  allocate(H0_ao(NAO, NAO))
                  call chol_H0(H0_ao, AOBasis, System)
                  !
                  ! MO coefficients in AO basis (auto2e ordering)
                  !
                  allocate(C_ao(NAO, NMO))
                  FromExternalAO = .true. ! AOs from external program -> AOs in the Auto2e format
                  TwoIndexTransf = .false. ! Transform only the index p of C(p,k)
                  call auto2e_interface_AngFuncTransf(C_ao, C, FromExternalAO, TwoIndexTransf, AOBasis, ExternalOrdering)
                  if (ExternalOrdering == ORBITAL_ORDERING_ORCA) then
                        call auto2e_interface_ApplyOrcaPhases_Matrix(C_ao, AOBasis, TwoIndexTransf)
                  end if
                  !
                  ! Bare nuclei hamiltonian (MO basis)
                  !
                  allocate(TransfWork(NAO, NMO))
                  call linalg_ab(TransfWork, H0_ao, C_ao)
                  call linalg_aTb(H0_mo, C_ao, TransfWork)
                  !
                  ! Eigenvectors of the density matrix
                  ! Occupation numbers
                  !
                  allocate(Ceig_mo(NMO, NMO))
                  allocate(OccNum(NMO, 1))
                  Ceig_mo = Rho_mo
                  call symmetric_eigenproblem(OccNum(:, 1), Ceig_mo, NMO, .true.)
                  do p = 1, NMO
                        OccNum(p, 1) = max(ZERO, OccNum(p, 1))
                  end do
                  allocate(Ceig_ao(NAO, NMO, 1))
                  call linalg_ab(Ceig_ao(:, :, 1), C_ao, Ceig_mo)
                  !
                  ! 1-RDM in AO basis (auto2e ordering)
                  !
                  allocate(Rho_ao(NAO, NAO, 1))
                  call linalg_ab(TransfWork, C_ao, Rho_mo)
                  call linalg_abT(Rho_ao(:, :, 1), TransfWork, C_ao)
                  !
                  ! Coulomb (J) and exchange (K) matrices
                  ! (AO basis, auto2e ordering)
                  !
                  NOcc(1) = NMO
                  NOcc(2) = 0
                  allocate(JK_ao(NAO, NAO, 1))
                  if ((.not. present(J_mo)) .or. (.not. present(K_mo))) then
                        CoulContrib = .true.
                        ExchContrib = .true.
                        call chf_JK(JK_ao, Rho_ao, Ceig_ao, CholeskyVecs, KScal, OccNum, NOcc, AOBasis, &
                              MaxBufferDimMB, CoulContrib, ExchContrib)
                  else
                        allocate(J_ao(NAO, NAO, 1))
                        allocate(K_ao(NAO, NAO, 1))
                        CoulContrib = .true.
                        ExchContrib = .false.
                        call chf_JK(J_ao, Rho_ao, Ceig_ao, CholeskyVecs, KScal, OccNum, NOcc, AOBasis, &
                              MaxBufferDimMB, CoulContrib, ExchContrib)
                        CoulContrib = .false.
                        ExchContrib = .true.
                        call chf_JK(K_ao, Rho_ao, Ceig_ao, CholeskyVecs, KScal, OccNum, NOcc, AOBasis, &
                              MaxBufferDimMB, CoulContrib, ExchContrib)
                        call linalg_smfill(J_ao(:, :, 1))
                        call linalg_smfill(K_ao(:, :, 1))
                        JK_ao = J_ao + K_ao
                        call linalg_ab(TransfWork, J_ao(:, :, 1), C_ao)
                        call linalg_aTb(J_mo, C_ao, TransfWork)
                        call linalg_ab(TransfWork, K_ao(:, :, 1), C_ao)
                        call linalg_aTb(K_mo, C_ao, TransfWork)
                  end if
                  call linalg_smfill(JK_ao(:, :, 1))
                  call linalg_ab(TransfWork, JK_ao(:, :, 1), C_ao)
                  call linalg_aTb(F_mo, C_ao, TransfWork)
                  F_mo = F_mo + H0_mo
            end associate
            call boys_free()
      end subroutine chol_F

      
      subroutine chol_H0_mo(H0_mo, C, AOBasis, System, ExternalOrdering)            
            real(F64), dimension(:, :), intent(out) :: H0_mo
            real(F64), dimension(:, :), intent(in)  :: C
            type(TAOBasis), intent(in)              :: AOBasis
            type(TSystem), intent(in)               :: System
            integer, intent(in)                     :: ExternalOrdering

            real(F64), dimension(:, :), allocatable :: C_ao
            real(F64), dimension(:, :), allocatable :: H0_ao
            real(F64), dimension(:, :), allocatable :: TransfWork
            integer :: NMO, NAO
            logical :: FromExternalAO, TwoIndexTransf

            call auto2e_init()
            call boys_init(4*AUTO2E_MAXL)
            associate (NAOCart => AOBasis%NAOCart, &
                  NAOSpher => AOBasis%NAOSpher, &
                  SpherAO => AOBasis%SpherAO)

                  NMO = size(C, dim=2)
                  NAO = size(C, dim=1)
                  if ((SpherAO .and. NAO /= NAOSpher) .or. (.not.SpherAO .and. NAO /= NAOCart)) then
                        call msg("Invalid dimension of matrix C")
                        error stop
                  end if
                  !
                  ! Bare nuclei hamiltonian (atomic orbital basis, auto2e ordering)
                  !
                  allocate(H0_ao(NAO, NAO))
                  call chol_H0(H0_ao, AOBasis, System)
                  !
                  ! MO coefficients in AO basis (auto2e ordering)
                  !
                  allocate(C_ao(NAO, NMO))
                  FromExternalAO = .true. ! AOs from external program -> AOs in the Auto2e format
                  TwoIndexTransf = .false. ! Transform only the index p of C(p,k)
                  call auto2e_interface_AngFuncTransf(C_ao, C, FromExternalAO, TwoIndexTransf, AOBasis, ExternalOrdering)
                  if (ExternalOrdering == ORBITAL_ORDERING_ORCA) then
                        call auto2e_interface_ApplyOrcaPhases_Matrix(C_ao, AOBasis, TwoIndexTransf)
                  end if
                  !
                  ! Bare nuclei hamiltonian (MO basis)
                  !
                  allocate(TransfWork(NAO, NMO))
                  call linalg_ab(TransfWork, H0_ao, C_ao)
                  call linalg_aTb(H0_mo, C_ao, TransfWork)
            end associate
            call boys_free()
      end subroutine chol_H0_mo
      

      subroutine chol_H0(H0_ao, AOBasis, System)
            real(F64), dimension(:, :), intent(out) :: H0_ao
            type(TAOBasis), intent(in)              :: AOBasis
            type(TSystem), intent(in)               :: System

            integer :: NAO
            real(F64), dimension(:, :), allocatable :: H0_cao
            real(F64), dimension(:), allocatable :: TransfWork

            associate (NAOCart => AOBasis%NAOCart, &
                  NAOSpher => AOBasis%NAOSpher, &
                  SpherAO => AOBasis%SpherAO)                  
                  NAO = size(H0_ao, dim=1)
                  if ((SpherAO .and. NAO /= NAOSpher) .or. (.not.SpherAO .and. NAO /= NAOCart)) then
                        call msg("Invalid dimension of matrix C_ao")
                        error stop
                  end if
                  if (SpherAO) then
                        allocate(H0_cao(NAOCart, NAOCart))
                        call chf_H0(H0_cao, AOBasis, System)
                        call linalg_smfill(H0_cao)
                        !
                        ! Transform the bare-nuclei hamiltonian to the spherical
                        ! AO basis.
                        !
                        allocate(TransfWork(NAOSpher*NAOCart))
                        call SpherGTO_TransformMatrix_U(H0_ao, H0_cao, &
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
                  else
                        call chf_H0(H0_ao, AOBasis, System)
                        call linalg_smfill(H0_ao)
                  end if
            end associate
      end subroutine chol_H0
end module Cholesky_driver
