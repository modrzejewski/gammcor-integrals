module THC_Gammcor
      use arithmetic
      use thc_definitions
      use sys_definitions
      use basis_sets
      use display
      use boys
      use Auto2e
      use Auto2eInterface
      use real_linalg
      use TensorHypercontraction
      use THCFock
      use clock
      
      implicit none

      integer, parameter :: THC_ACCURACY_DEFAULT = 1
      integer, parameter :: THC_ACCURACY_TIGHT = 2
      integer, parameter :: THC_ACCURACY_LUDICROUS = 3
      integer, parameter :: THC_ACCURACY_DEBUG = 4

contains

      subroutine thc_gammcor_Rkab(Rkab, CA, a0, a1, CB, b0, b1, BasisSetPath, XYZPath, Accuracy, &
            SpherAO, ExternalOrdering, SortAngularMomenta, Units)

            real(F64), dimension(:, :), allocatable, intent(out) :: Rkab
            real(F64), dimension(:, :), intent(in)               :: CA
            integer, intent(in)                                  :: a0, a1
            real(F64), dimension(:, :), intent(in)               :: CB
            integer, intent(in)                                  :: b0, b1
            character(*), intent(in)                             :: BasisSetPath
            character(*), intent(in)                             :: XYZPath
            integer, intent(in)                                  :: Accuracy
            logical, intent(in)                                  :: SpherAO
            integer, intent(in)                                  :: ExternalOrdering
            logical, intent(in)                                  :: SortAngularMomenta
            integer, intent(in)                                  :: Units

            type(TAOBasis) :: AOBasis
            type(TSystem) :: System
            real(F64), dimension(:, :), allocatable :: Xgp, Zgk
            real(F64), dimension(:, :), allocatable :: Xga, Xgb            
            integer :: NGridTHC, NCholesky, NA, NB
            integer :: NAO
            !
            ! Read the XYZ coordinates and atom types
            !
            call sys_Read_XYZ(System, XYZPath, Units)
            !
            ! Read the basis set parameters from an EMSL text file
            ! (GAMESS-US format, no need for any edits, just download it straight from the website)
            !
            call basis_NewAOBasis(AOBasis, System, BasisSetPath, SpherAO, SortAngularMomenta)
            if (AOBasis%SpherAO) then
                  NAO = AOBasis%NAOSpher
            else
                  NAO = AOBasis%NAOCart
            end if
            NA = a1 - a0 + 1
            NB = b1 - b0 + 1
            call thc_gammcor_XZ(Xgp, Zgk, AOBasis, System, Accuracy)
            NGridTHC = size(Zgk, dim=1)
            NCholesky = size(Zgk, dim=2)
            allocate(Xga(NGridTHC, NA))
            allocate(Xgb(NGridTHC, NB))
            allocate(Rkab(NCholesky, NA*NB))
            call thc_gammcor_Xga(Xga, Xgp, CA(:, a0:a1), AOBasis, ExternalOrdering)
            call thc_gammcor_Xga(Xgb, Xgp, CB(:, b0:b1), AOBasis, ExternalOrdering)
            call thc_gammcor_Rkab_2(Rkab, Xga, Xgb, Zgk, Na, Nb, NCholesky, NGridTHC)        
      end subroutine thc_gammcor_Rkab
      

      subroutine thc_gammcor_Rkab_2(Rkab, Xga, Xgb, Zgk, Na, Nb, NCholesky, NGridTHC)
            integer, intent(in)                                   :: Na, Nb
            integer, intent(in)                                   :: NCholesky, NGridTHC
            real(F64), dimension(NCholesky, Na, Nb), intent(out)  :: Rkab
            real(F64), dimension(NGridTHC, Na), intent(in)        :: Xga
            real(F64), dimension(NGridTHC, Nb), intent(in)        :: Xgb
            real(F64), dimension(NGridTHC, NCholesky), intent(in) :: Zgk

            integer :: b, k
            real(F64), dimension(:, :), allocatable :: XZgk

            allocate(XZgk(NGridTHC, NCholesky))
            do b = 1, Nb
                  !$omp parallel do private(k)
                  do k = 1, NCholesky
                        XZgk(:, k) = Xgb(:, b) * Zgk(:, k)
                  end do
                  !$omp end parallel do
                  !
                  ! R(k,a;b) = Sum(g) [XZ](g,k;b)*X(g,a)
                  !
                  call real_aTb(Rkab(:, :, b), XZgk, Xga)
            end do
      end subroutine thc_gammcor_Rkab_2


      subroutine thc_gammcor_Rkab_Batch_a_Fixed_b(Rkab, XXga, Xga, Xgb, Zgk, Na, NCholesky, NGridTHC)
            integer, intent(in)                                   :: Na
            integer, intent(in)                                   :: NCholesky, NGridTHC
            real(F64), dimension(NCholesky, Na), intent(out)      :: Rkab
            real(F64), dimension(NGridTHC, Na), intent(out)       :: XXga
            real(F64), dimension(NGridTHC, Na), intent(in)        :: Xga
            real(F64), dimension(NGridTHC), intent(in)            :: Xgb
            real(F64), dimension(NGridTHC, NCholesky), intent(in) :: Zgk

            integer :: a

            !$omp parallel do private(a)
            do a = 1, Na
                  XXga(:, a) = Xga(:, a) * Xgb(:)
            end do
            !$omp end parallel do
            call real_aTb(Rkab, Zgk, XXga)
      end subroutine thc_gammcor_Rkab_Batch_a_Fixed_b


      subroutine thc_gammcor_Vabcd_Batch_ab_Fixed_cd(Vabcd, XXgab, ZXXkab, ZXXkcd, &
            Xga, Xgb, Xgc, Xgd, Zgk, Na, Nb, NCholesky, NGridTHC)
            
            integer, intent(in)                                   :: Na, Nb
            integer, intent(in)                                   :: NCholesky, NGridTHC
            real(F64), dimension(Na, Nb), intent(out)             :: Vabcd
            real(F64), dimension(NGridTHC, Na), intent(out)       :: XXgab
            real(F64), dimension(NCholesky, Na), intent(out)      :: ZXXkab
            real(F64), dimension(NCholesky), intent(out)          :: ZXXkcd
            real(F64), dimension(NGridTHC, Na), intent(in)        :: Xga
            real(F64), dimension(NGridTHC, Nb), intent(in)        :: Xgb
            real(F64), dimension(NGridTHC), intent(in)            :: Xgc
            real(F64), dimension(NGridTHC), intent(in)            :: Xgd
            real(F64), dimension(NGridTHC, NCholesky), intent(in) :: Zgk

            integer :: a, b

            XXgab(:, 1) = Xgc(:) * Xgd(:)
            call real_ATv(ZXXkcd, Zgk, XXgab(:, 1))
            do b = 1, Nb
                  !$omp parallel do private(a)
                  do a = 1, Na
                        XXgab(:, a) = Xga(:, a) * Xgb(:, b)
                  end do
                  !$omp end parallel do
                  call real_aTb(ZXXkab, Zgk, XXgab)
                  call real_ATv(Vabcd(:, b), ZXXkab, ZXXkcd)
            end do
      end subroutine thc_gammcor_Vabcd_Batch_ab_Fixed_cd


      subroutine thc_gammcor_Xga(Xga, Xgp, C_extao, AOBasis, ExternalOrdering)
            real(F64), dimension(:, :), intent(out) :: Xga
            real(F64), dimension(:, :), intent(in)  :: Xgp
            real(F64), dimension(:, :), intent(in)  :: C_extao
            type(TAOBasis), intent(in)              :: AOBasis
            integer, intent(in)                     :: ExternalOrdering

            real(F64), dimension(:, :), allocatable :: C_ao
            !
            ! Change the ordering of angular functions from the external convention
            ! to the internal ordering of the Auto2e module.
            ! The ordering of whole shells is assumed to be already
            ! the same as in the external program.
            !
            allocate(C_ao, mold=C_extao)
            call auto2e_interface_C(C_ao, C_extao, AOBasis, ExternalOrdering)
            call real_ab(Xga, Xgp, C_ao)
      end subroutine thc_gammcor_Xga

      
      subroutine thc_gammcor_XZ(Xgp, Zgk, AOBasis, System, Accuracy, &
            CholeskyThresh, THCThresh, Omega)
            
            real(F64), dimension(:, :), allocatable, intent(out) :: Xgp
            real(F64), dimension(:, :), allocatable, intent(out) :: Zgk
            type(TAOBasis), intent(in)                           :: AOBasis
            type(TSystem), intent(in)                            :: System
            integer, intent(in)                                  :: Accuracy
            real(F64), optional, intent(in)                      :: CholeskyThresh
            real(F64), optional, intent(in)                      :: THCThresh
            real(F64), optional, intent(in)                      :: Omega

            type(TCoulTHCGrid) :: THCGrid
            type(TTHCParams) :: THCParams
            type(TChol2Params) :: Chol2Params
            real(F64) :: Kappa

            THCParams%THC_QuadraticMemory = .true.
            select case (Accuracy)
            case (THC_ACCURACY_DEFAULT)
                  Chol2Params%CholeskyTauThresh = 1.0E-5_F64
                  THCParams%QRThresh = 1.0E-3_F64
            case (THC_ACCURACY_TIGHT)
                  Chol2Params%CholeskyTauThresh = 1.0E-6_F64
                  THCParams%QRThresh = 1.0E-4_F64
            case (THC_ACCURACY_LUDICROUS)
                  Chol2Params%CholeskyTauThresh = 1.0E-7_F64
                  THCParams%QRThresh = 1.0E-5_F64
            case (THC_ACCURACY_DEBUG)
                  Chol2Params%CholeskyTauThresh = 1.0E-10_F64
                  THCParams%QRThresh = 1.0E-6_F64
            case default
                  call msg("Invalid THC accuracy value", MSG_ERROR)
                  error stop
            end select
            if (present(CholeskyThresh)) Chol2Params%CholeskyTauThresh = CholeskyThresh
            if (present(THCThresh)) THCParams%QRThresh = THCThresh
            THCParams%QRThreshReduced = THCParams%QRThresh
            if (present(Omega)) then
                  !
                  ! If Omega is present, long-range Coulomb integrals
                  ! will be computed with the Erf(Omega*r)/r operator
                  !
                  if (Omega > ZERO) then
                        Kappa = ONE / Omega**2
                  else
                        Kappa = ZERO
                  end if
            else
                  Kappa = ZERO
            end if
            Chol2Params%Kappa = Kappa
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
            call thc_CoulombMatrix_QuadraticMemory(THCGrid, AOBasis, System, THCParams, Chol2Params)
            call move_alloc(to=Zgk, from=THCGrid%Zgk)
            call move_alloc(to=Xgp, from=THCGrid%Xgp)
            !
            ! Deallocate the Boys function interpolation tables
            ! to avoid allocation of already allocated arrays
            ! if this subroutine is called again
            !
            call boys_free()
      end subroutine thc_gammcor_XZ


      subroutine thc_gammcor_F(Fij, Fvw, Cpi_extao, Cpa_extao, Cpv_extao, &
            OccNum, Zgk, Xgp, AOBasis, System, ExternalOrdering)
            !
            ! Compute the Fock matrix, F, using THC-decomposed Coulomb integrals
            !
            ! F(p,q) = T(p,q) + Vne(p,q)
            !        + Sum(k=1,...,NInactive+NActive) OccNum(i) (2 * (pq|kk) - (pk|kq))
            !
            ! The THC Fock matrix shoud be used to compute the canonical
            ! inactive and virtual orbitals.
            !
            ! Fij
            ! Fvw
            !              Inactive-inactive (ij) and virtual-virtual (vw) blocks of F
            !
            ! Cpi_extao
            ! Cpa_extao
            ! Cpv_extao
            !              Inactive (i), active (a), and virtual (v) orbital
            !              coefficients in the AO basis (AO ordering of the external
            !              program)
            ! OccNum
            !              Occupation numbers
            !              0 <= OccNum(k) <= 1  for k = 1, ..., NInactive+NActive
            !
            ! Zgk, Xgp
            !              Z and X matrices of the THC decomposition
            !
            ! AOBasis
            !              Basis set definition
            ! System
            !              System definition
            !
            ! ExternalOrdering
            !              Ordering of the AO funtions in the external program
            !
            real(F64), dimension(:, :), intent(out) :: Fij
            real(F64), dimension(:, :), intent(out) :: Fvw
            real(F64), dimension(:, :), intent(in)  :: Cpi_extao
            real(F64), dimension(:, :), intent(in)  :: Cpa_extao
            real(F64), dimension(:, :), intent(in)  :: Cpv_extao
            real(F64), dimension(:), intent(in )    :: OccNum
            real(F64), dimension(:, :), intent(in)  :: Zgk
            real(F64), dimension(:, :), intent(in)  :: Xgp
            type(TAOBasis), intent(in)              :: AOBasis
            type(TSystem), intent(in)               :: System
            integer, intent(in)                     :: ExternalOrdering

            integer :: NOccupied, NActive, NInactive, NVirtual, NMO, NAO, NGridTHC
            integer :: i0, i1, a0, a1, k
            real(F64), dimension(:, :), allocatable :: Cpi, Cpv
            real(F64), dimension(:, :, :), allocatable :: Fpq, Cpo
            real(F64), dimension(:, :), allocatable :: Zgh
            real(F64), dimension(:, :), allocatable :: Fpi, Fpv
            real(F64) :: Nk
            integer, dimension(2) :: NOcc
            type(TClock) :: timer
            
            call clock_start(timer)
            call msg("Started THC Fock matrix build")
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
            NAO = AOBasis%NAOSpher
            NActive = size(Cpa_extao, dim=2)
            NInactive = size(Cpi_extao, dim=2)
            NVirtual = size(Cpv_extao, dim=2)
            NOccupied = NInactive + NActive
            NMO = NInactive + NActive + NVirtual
            NGridTHC = size(Xgp, dim=1)
            if (size(OccNum) /= NOccupied) then
                  call msg("Invalid size of the OccNum array", MSG_ERROR)
                  error stop
            end if
            allocate(Cpi(NAO, NInactive))
            allocate(Cpv(NAO, NVirtual))
            allocate(Cpo(NAO, NOccupied, 1))
            allocate(Fpq(NAO, NAO, 1))
            allocate(Zgh(NGridTHC, NGridTHC))
            allocate(Fpi(NAO, NInactive))
            allocate(Fpv(NAO, NVirtual))
            i0 = 1
            i1 = NInactive
            a0 = NInactive + 1
            a1 = NInactive + NActive
            call auto2e_interface_C(Cpi, Cpi_extao, AOBasis, ExternalOrdering)
            call auto2e_interface_C(Cpv, Cpv_extao, AOBasis, ExternalOrdering)
            Cpo(:, i0:i1, 1) = Cpi(:, :)
            call auto2e_interface_C(Cpo(:, a0:a1, 1), Cpa_extao, AOBasis, ExternalOrdering)
            do k = 1, NOccupied
                  Nk = max(ZERO, OccNum(k))
                  Cpo(:, k, 1) = Sqrt(Nk) * Cpo(:, k, 1)
            end do
            call real_abT(Zgh, Zgk, Zgk)
            NOcc(1) = NOccupied
            NOcc(2) = 0
            call thc_Fock_F(Fpq, Cpo, NOcc, Zgh, Xgp, AOBasis, System)
            call real_ab(Fpi, Fpq(:, :, 1), Cpi)
            call real_aTb(Fij, Cpi, Fpi)
            call real_ab(Fpv, Fpq(:, :, 1), Cpv)
            call real_aTb(Fvw, Cpv, Fpv)
            call msg("Fock matrix completed in " // str(clock_readwall(timer),d=1) // " seconds")
            !
            ! Deallocate the Boys function interpolation tables
            ! to avoid allocation of already allocated arrays
            ! if this subroutine is called again
            !
            call boys_free()
      end subroutine thc_gammcor_F
end module THC_Gammcor
