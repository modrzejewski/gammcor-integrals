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

      
      subroutine thc_gammcor_XZ(Xgp, Zgk, AOBasis, System, Accuracy)
            real(F64), dimension(:, :), allocatable, intent(out) :: Xgp
            real(F64), dimension(:, :), allocatable, intent(out) :: Zgk
            type(TAOBasis), intent(in)                           :: AOBasis
            type(TSystem), intent(in)                            :: System
            integer, intent(in)                                  :: Accuracy

            type(TCoulTHCGrid) :: THCGrid
            type(TTHCParams) :: THCParams
            type(TChol2Params) :: Chol2Params

            THCParams%THC_QuadraticMemory = .true.
            select case (Accuracy)
            case (THC_ACCURACY_DEFAULT)
                  Chol2Params%CholeskyTauThresh = 1.0E-7_F64
                  THCParams%QRThresh = 1.0E-4_F64
            case (THC_ACCURACY_TIGHT)
                  Chol2Params%CholeskyTauThresh = 1.0E-7_F64
                  THCParams%QRThresh = 1.0E-5_F64
            case (THC_ACCURACY_LUDICROUS)
                  Chol2Params%CholeskyTauThresh = 1.0E-7_F64
                  THCParams%QRThresh = 1.0E-6_F64
            case (THC_ACCURACY_DEBUG)
                  Chol2Params%CholeskyTauThresh = 1.0E-10_F64
                  THCParams%QRThresh = 1.0E-7_F64
            case default
                  call msg("Invalid THC accuracy value", MSG_ERROR)
                  error stop
            end select
            THCParams%QRThreshReduced = THCParams%QRThresh
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
end module THC_Gammcor
