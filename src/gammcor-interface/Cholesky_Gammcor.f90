module Cholesky_Gammcor
      use arithmetic
      use TwoStepCholesky
      use TwoStepCholesky_definitions
      use real_linalg
      use clock
      use Auto2eInterface
      use sys_definitions
      use basis_sets
      use OneElectronInts_Gammcor

      implicit none
      
      integer, parameter :: CHOL_ACCURACY_DEFAULT = 1
      integer, parameter :: CHOL_ACCURACY_TIGHT = 2
      integer, parameter :: CHOL_ACCURACY_LUDICROUS = 3
      integer, parameter :: CHOL_ACCURACY_DEBUG = 4

      type TCholeskyVecsOTF
            real(F64), dimension(:, :, :), allocatable :: Rkpq
            type(TChol2Vecs) :: Chol2Data
      end type TCholeskyVecsOTF
      
contains

      subroutine chol_gammcor_Rkab(Rkab, CA, a0, a1, CB, b0, b1, MaxBufferDimMB, CholeskyVecsOTF, &
            AOBasis, ExternalOrdering)
            !
            ! Transformation of Cholesky vectors computed in AO basis with the on the fly
            ! algorithm.
            !
            ! Compatible with MO coeffs generated with the programs listed in the Auto2eInterface
            ! module.
            !
            real(F64), dimension(:, :), intent(out)    :: Rkab
            real(F64), dimension(:, :), intent(in)     :: CA
            integer, intent(in)                        :: a0, a1
            real(F64), dimension(:, :), intent(in)     :: CB
            integer, intent(in)                        :: b0, b1
            integer, intent(in)                        :: MaxBufferDimMB
            type(TCholeskyVecsOTF), intent(in)         :: CholeskyVecsOTF
            type(TAOBasis), intent(in)                 :: AOBasis
            integer, intent(in)                        :: ExternalOrdering

            integer :: NCholesky, NA, NB
            real(F64), dimension(:, :), allocatable :: CA_ao, CB_ao
            logical :: TwoIndexTransf, FromExternalAO, SpherAO

            NA = a1 - a0 + 1
            NB = b1 - b0 + 1
            SpherAO = AOBasis%SpherAO
            NCholesky = CholeskyVecsOTF%Chol2Data%NVecs
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
            ! Change the ordering of angular functions from the external convention
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
            call chol_gammcor_MOTransf(Rkab, &
                  CholeskyVecsOTF%Rkpq, &
                  CholeskyVecsOTF%Chol2Data, &
                  CA_ao, CB_ao, a0, a1, b0, b1, &
                  MaxBufferDimMB, AOBasis)
      end subroutine chol_gammcor_Rkab
      

      subroutine chol_gammcor_MOTransf(Rkab, Rkpq, CholeskyVecs, CA, CB, a0, a1, b0, b1, &
            MaxBufferDimMB, AOBasis)
            
            real(F64), dimension(:, :), intent(out)                    :: Rkab
            real(F64), dimension(:, :, :), intent(in)                  :: Rkpq
            type(TChol2Vecs), intent(in)                               :: CholeskyVecs
            real(F64), dimension(:, :), intent(in)                     :: CA
            real(F64), dimension(:, :), intent(in)                     :: CB
            integer, intent(in)                                        :: a0, a1
            integer, intent(in)                                        :: b0, b1
            integer, intent(in)                                        :: MaxBufferDimMB
            type(TAOBasis), intent(in)                                 :: AOBasis

            integer(I64) :: MaxBufferDim

            MaxBufferDim = (int(MaxBufferDimMB,I64)*1024*1024)/(storage_size(Rkab)/8)
            if (AOBasis%SpherAO) then
                  call chol_gammcor_MOTransf_2(Rkab, &
                        Rkpq, &
                        CA, CB, a0, a1, b0, b1, &
                        AOBasis%NAOSpher, &
                        CholeskyVecs%NVecs, &
                        MaxBufferDim, &
                        CholeskyVecs%SubsetBounds, &
                        CholeskyVecs%ShellPairs, &
                        CholeskyVecs%ShellPairLoc, &
                        AOBasis%ShellLocSpher, &
                        AOBasis%ShellParamsIdx, &
                        AOBasis%NAngFuncSpher, &
                        CholeskyVecs%NSubsets)
            else
                  call chol_gammcor_MOTransf_2(Rkab, &
                        Rkpq, &
                        CA, CB, a0, a1, b0, b1, &
                        AOBasis%NAOCart, &
                        CholeskyVecs%NVecs, &
                        MaxBufferDim, &
                        CholeskyVecs%SubsetBounds, &
                        CholeskyVecs%ShellPairs, &
                        CholeskyVecs%ShellPairLoc, &
                        AOBasis%ShellLocCart, &
                        AOBasis%ShellParamsIdx, &
                        AOBasis%NAngFuncCart, &
                        CholeskyVecs%NSubsets)
            end if
      end subroutine chol_gammcor_MOTransf
      

      subroutine chol_gammcor_MOTransf_2(Rab, Rpq, CA, CB, a0, a1, b0, b1, NAO, NCholesky, MaxBufferDim, &
            SubsetBounds, ShellPairs, ShellPairLoc, ShellLoc, ShellParamsIdx, NAngFunc, NSubsets)
            
            real(F64), dimension(:, :), intent(out)                    :: Rab
            real(F64), dimension(:, :, :), intent(in)                  :: Rpq
            real(F64), dimension(:, :), intent(in)                     :: CA
            real(F64), dimension(:, :), intent(in)                     :: CB
            integer, intent(in)                                        :: a0, a1
            integer, intent(in)                                        :: b0, b1
            integer, intent(in)                                        :: NAO
            integer, intent(in)                                        :: NCholesky
            integer(I64), intent(in)                                   :: MaxBufferDim
            integer, dimension(:, :), intent(in)                       :: SubsetBounds
            integer, dimension(:, :), intent(in)                       :: ShellPairs
            integer, dimension(:, :), intent(in)                       :: ShellPairLoc
            integer, dimension(:), intent(in)                          :: ShellLoc
            integer, dimension(:), intent(in)                          :: ShellParamsIdx
            integer, dimension(:), intent(in)                          :: NAngFunc
            integer, dimension(2), intent(in)                          :: NSubsets

            integer :: NA, NB
            integer :: MaxNk, NPasses, Nk
            integer :: k, a0k, a1k
            real(F64), dimension(:, :), allocatable :: CAT
            real(F64), dimension(:), allocatable :: W
            integer :: X, Y, SubsetIdx
            type(TClock) :: timer, timer_total

            Y = 1
            NA = a1 - a0 + 1
            NB = b1 - b0 + 1
            allocate(CAT(NA, NAO))
            CAT = transpose(CA(:, a0:a1))
            MaxNk = min(NA, int(MaxBufferDim/(NCholesky*NAO)))
            if (.not. MaxNk > 0) then
                  call msg("Buffer size too small for Rpq->Rab")
                  error stop
            end if
            NPasses = NA / MaxNk
            if (modulo(NA, MaxNk) > 0) then
                  NPasses = NPasses + 1
            end if
            allocate(W(NCholesky*MaxNk*NAO))
            call msg("Two step transformation of Cholesky vectors: R(k,pq)->R(k,ab)")
            call msg("Will perform " // str(NPasses) // " passes over R(k,pq)")
            call clock_start(timer_total)
            Rab = ZERO
            do k = 1, NPasses
                  call clock_start(timer)
                  ! ----------------------------------------------------------------
                  ! First index transformation R(:, pq) -> R(:, aq) for a=a0k...a1k
                  ! ----------------------------------------------------------------
                  W = ZERO
                  a0k = (k - 1) * MaxNk + a0
                  a1k = min(k * MaxNk, a1)
                  Nk = a1k - a0k + 1
                  do X = 1, NSubsets(1)
                        SubsetIdx = X + (Y - 1) * NSubsets(1)
                        call chol_gammcor_MOTransf_Step1(W, a0k, a1k, Rpq(:, :, X), NCholesky, CAT, a0, a1, &
                              SubsetBounds(:, SubsetIdx), ShellPairs, ShellPairLoc, ShellLoc, ShellParamsIdx, &
                              NAngFunc, NAO)
                  end do
                  ! -----------------------------------------------------------------------------
                  ! Second index transformation R(:, aq) -> R(:, ab) for a=a0k...a1k, b=b0...b1
                  ! -----------------------------------------------------------------------------
                  call chol_gammcor_MOTransf_Step2(Rab, W, CB, a0k, a1k, Nk, a0, a1, b0, b1, NCholesky, NAO)
                  call msg("Completed pass " // str(k) // " in " // str(clock_readwall(timer),d=1) // " seconds")
            end do
            call msg("Completed transformation in " // str(clock_readwall(timer_total),d=1) // " seconds")
      end subroutine chol_gammcor_MOTransf_2
      
      
      subroutine chol_gammcor_MOTransf_Step1(Wkap, a0k, a1k, Rkpq, NCholesky, CAT, a0, a1, &
            SubsetBounds, ShellPairs, ShellPairLoc, ShellLoc, ShellParamsIdx, &
            NAngFunc, NAO)

            integer, intent(in)                                            :: a0, a1
            integer, intent(in)                                            :: a0k, a1k
            real(F64), dimension(1:NCholesky, a0k:a1k, NAO), intent(inout) :: Wkap
            real(F64), dimension(:, :), intent(in)                         :: Rkpq
            integer, intent(in)                                            :: NCholesky
            real(F64), dimension(a0:a1, NAO), intent(in)                   :: CAT
            integer, dimension(2), intent(in)                              :: SubsetBounds
            integer, dimension(:, :), intent(in)                           :: ShellPairs
            integer, dimension(:, :), intent(in)                           :: ShellPairLoc
            integer, dimension(:), intent(in)                              :: ShellLoc
            integer, dimension(:), intent(in)                              :: ShellParamsIdx
            integer, dimension(:), intent(in)                              :: NAngFunc
            integer, intent(in)                                            :: NAO

            integer :: Nk
            integer :: ShAB, LocAB
            integer :: ShA, ShellParamsA, Na, LocA
            integer :: ShB, ShellParamsB, Nb, LocB
            integer :: a, b, p, q, pq

            Nk = a1k - a0k + 1
            do ShAB = SubsetBounds(1), SubsetBounds(2)
                  LocAB = ShellPairLoc(CHOL2_SUBSET_STORAGE, ShAB)

                  ShA = ShellPairs(1, ShAB)
                  ShellParamsA = ShellParamsIdx(ShA)
                  Na = NAngFunc(ShellParamsA)
                  LocA = ShellLoc(ShA)

                  ShB = ShellPairs(2, ShAB)
                  ShellParamsB = ShellParamsIdx(ShB)
                  Nb = NAngFunc(ShellParamsB)
                  LocB = ShellLoc(ShB)

                  if (ShA /= ShB) then
                        pq = LocAB - 1
                        do b = 1, Nb
                              do a = 1, Na
                                    p = LocA + a - 1
                                    q = LocB + b - 1
                                    pq = pq + 1
                                    !
                                    ! W(:, :, p) <- W(:, :, p) + R(1:NVecs, pq) * CA(q, a0k:a1k)**T
                                    ! W(:, :, q) <- W(:, :, q) + R(1:NVecs, pq) * CA(p, a0k:a1k)**T
                                    !
                                    call real_vwT_x(Wkap(:, :, p), NCholesky, Rkpq(:, pq), CAT(a0k:a1k, q), NCholesky, Nk, ONE)
                                    call real_vwT_x(Wkap(:, :, q), NCholesky, Rkpq(:, pq), CAT(a0k:a1k, p), NCholesky, Nk, ONE)
                              end do
                        end do
                  else
                        pq = LocAB - 1
                        do b = 1, Nb
                              pq = pq + 1
                              q = LocB + b - 1
                              p = q
                              !
                              ! W(:, :, p) <- W(:, :, p) + R(1:NVecs, pp) * CA(p, a0k:a1k)**T
                              !
                              call real_vwT_x(Wkap(:, :, p), NCholesky, Rkpq(:, pq), CAT(a0k:a1k, q), NCholesky, Nk, ONE)
                              do a = b+1, Na
                                    p = LocA + a - 1
                                    q = LocB + b - 1
                                    pq = pq + 1
                                    !
                                    ! W(:, :, p) <- W(:, :, p) + R(1:NVecs, pq) * CA(q, a0k:a1k)**T
                                    ! W(:, :, q) <- W(:, :, q) + R(1:NVecs, pq) * CA(p, a0k:a1k)**T
                                    !
                                    call real_vwT_x(Wkap(:, :, p), NCholesky, Rkpq(:, pq), CAT(a0k:a1k, q), NCholesky, Nk, ONE)
                                    call real_vwT_x(Wkap(:, :, q), NCholesky, Rkpq(:, pq), CAT(a0k:a1k, p), NCholesky, Nk, ONE)
                              end do
                        end do
                  end if
            end do
      end subroutine chol_gammcor_MOTransf_Step1


      subroutine chol_gammcor_MOTransf_Step2(Rab, W, CB, a0k, a1k, Nk, a0, a1, b0, b1, NCholesky, NAO)
            integer, intent(in)                                           :: a0, a1
            integer, intent(in)                                           :: b0, b1
            real(F64), dimension(NCholesky, a0:a1, b0:b1), intent(inout)  :: Rab
            real(F64), dimension(NCholesky, Nk, NAO), intent(in)          :: W
            real(F64), dimension(:, :), intent(in)                        :: CB
            integer, intent(in)                                           :: a0k, a1k
            integer, intent(in)                                           :: Nk
            integer, intent(in)                                           :: NCholesky
            integer, intent(in)                                           :: NAO

            integer :: m, n, ldW, b
            !
            ! Transform the second AO index by a series of vector-matrix multiplications
            ! --------------------------------------------------------------------------
            !
            ! Explicit loop:
            ! do b = b0, b1
            !       do p = 1, NAO
            !             Rab(:, a0k:a1k, b) = Rab(:, a0k:a1k, b) + W(:, :, p) * CB(p, b)
            !       end do
            ! end do
            !
            ! Call a linear algebra subroutine for each of individual
            ! vector-matrix multiplications
            !
            ! Rab(:, a0k:a1k, b) <- Rab(:, a0k:a1k, b) + Sum(p=1...NAO) W(:, :, p) CB(p, b)
            !
            ! Note that the above cannot be done as a single matrix multiplication due to
            ! the non-contiguous memory layout of Rab.
            !
            m = NCholesky * Nk
            n = NAO
            ldW = NCholesky * Nk
            do b = b0, b1
                  call real_av_x(Rab(:, a0k:a1k, b), W, ldW, CB(:, b), m, n, ONE, ONE)
            end do
      end subroutine chol_gammcor_MOTransf_Step2


      subroutine chol_gammcor_Rkpq(CholeskyVecsOTF, AOBasis, Accuracy, Omega)
            !
            ! Cholesky vectors in AO basis computed on the fly.
            !
            type(TCholeskyVecsOTF), intent(out)                     :: CholeskyVecsOTF
            type(TAOBasis), intent(in)                              :: AOBasis
            integer, intent(in)                                     :: Accuracy
            real(F64), optional, intent(in)                         :: Omega
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
            if (present(Omega)) then
                  call chol_gammcor_CoulombMatrix( &
                        CholeskyVecsOTF%Rkpq, &
                        CholeskyVecsOTF%Chol2Data, &
                        Accuracy, AOBasis, Omega)
            else
                  call chol_gammcor_CoulombMatrix( &
                        CholeskyVecsOTF%Rkpq, &
                        CholeskyVecsOTF%Chol2Data, &
                        Accuracy, AOBasis)
            end if
            !
            ! Deallocate the Boys function interpolation tables
            ! to avoid allocation of already allocated arrays
            ! if this subroutine is called again
            !
            call boys_free()
      end subroutine chol_gammcor_Rkpq
      
      
      subroutine chol_gammcor_CoulombMatrix(Rkpq, Chol2Vecs, Accuracy, AOBasis, Omega)
            !
            ! Generate the Cholesky vectors matrix R: V = R**T * R, where V is the matrix
            ! of two electron Coulomb integrals (pq|rs).
            !
            real(F64), dimension(:, :, :), allocatable, intent(out) :: Rkpq
            type(TChol2Vecs), intent(out)                           :: Chol2Vecs
            integer, intent(in)                                     :: Accuracy
            type(TAOBasis), intent(in)                              :: AOBasis
            real(F64), optional, intent(in)                         :: Omega

            real(F64) :: Kappa
            type(TChol2Params) :: Chol2Params
            real(F64), dimension(:, :), allocatable :: Wabrs
            integer :: SubsetIdx, X, Y
            integer :: MaxSubsetDim
            type(TClock) :: timer
            
            select case (Accuracy)
            case (CHOL_ACCURACY_DEFAULT)
                  Chol2Params%CholeskyTauThresh = 1.0E-5_F64
            case (CHOL_ACCURACY_TIGHT)
                  Chol2Params%CholeskyTauThresh = 1.0E-6_F64
            case (CHOL_ACCURACY_LUDICROUS)
                  Chol2Params%CholeskyTauThresh = 1.0E-7_F64
            case (CHOL_ACCURACY_DEBUG)
                  Chol2Params%CholeskyTauThresh = 1.0E-8_F64
            case default
                  call msg("Invalid accuracy setting", MSG_ERROR)
                  error stop
            end select
            
            if (present(Omega)) then
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
            ! Step 1: locate pivots
            !
            call chol2_Algo_Koch_JCP2019(Chol2Vecs, AOBasis, Chol2Params)
            !
            ! Step 2: compute full-dimensional Cholesky vectors
            !
            call clock_start(timer)
            call msg("Step 2: Computing full-dimensional vectors")
            MaxSubsetDim = maxval(Chol2Vecs%SubsetDim)
            allocate(Rkpq(Chol2Vecs%NVecs, MaxSubsetDim, Chol2Vecs%NSubsets(1)))
            call chol2_AllocWorkspace(Wabrs, Chol2Vecs)
            Y = this_image()
            do X = 1, Chol2Vecs%NSubsets(1)
                  SubsetIdx = X + (Y - 1) * Chol2Vecs%NSubsets(1)
                  call chol2_FullDimVectors_Batch(Rkpq(:, :, X), Wabrs, SubsetIdx, &
                        Chol2Vecs, AOBasis, Chol2Params)
            end do
            call msg("Step 2 completed in " // str(clock_readwall(timer),d=1) // " seconds")            
      end subroutine chol_gammcor_CoulombMatrix


      ! subroutine chol_gammcor_F(F_mo, H0_mo, Rho_mo, C, KScal, CholeskyVecs, AOBasis, &
      !       System, ExternalOrdering, MaxBufferDimMB, J_mo, K_mo)
            
      !       real(F64), dimension(:, :), intent(out)   :: F_mo
      !       real(F64), dimension(:, :), intent(out)   :: H0_mo
      !       real(F64), dimension(:, :), intent(in)    :: Rho_mo
      !       real(F64), dimension(:, :), intent(in)    :: C
      !       real(F64), intent(in)                     :: KScal
      !       type(TCholeskyVecsOTF), intent(in)        :: CholeskyVecs
      !       type(TAOBasis), intent(in)                :: AOBasis
      !       type(TSystem), intent(in)                 :: System
      !       integer, intent(in)                       :: ExternalOrdering
      !       integer, intent(in)                       :: MaxBufferDimMB
      !       real(F64), dimension(:, :), optional, intent(out) :: J_mo
      !       real(F64), dimension(:, :), optional, intent(out) :: K_mo

      !       real(F64), dimension(:, :), allocatable :: C_ao
      !       real(F64), dimension(:, :), allocatable :: Ceig_mo
      !       real(F64), dimension(:, :, :), allocatable :: Ceig_ao
      !       real(F64), dimension(:, :), allocatable :: OccNum
      !       real(F64), dimension(:, :), allocatable :: H0_ao
      !       real(F64), dimension(:, :), allocatable :: TransfWork
      !       real(F64), dimension(:, :, :), allocatable :: Rho_ao, JK_ao, J_ao, K_ao
      !       integer :: NMO, NAO
      !       logical :: FromExternalAO, TwoIndexTransf
      !       integer :: p
      !       integer, dimension(2) :: NOcc
      !       logical :: CoulContrib, ExchContrib

      !       call auto2e_init()
      !       call boys_init(4*AUTO2E_MAXL)
      !       associate (NAOCart => AOBasis%NAOCart, &
      !             NAOSpher => AOBasis%NAOSpher, &
      !             SpherAO => AOBasis%SpherAO)

      !             NMO = size(C, dim=2)
      !             NAO = size(C, dim=1)
      !             if ((SpherAO .and. NAO /= NAOSpher) .or. (.not.SpherAO .and. NAO /= NAOCart)) then
      !                   call msg("Invalid dimension of matrix C")
      !                   error stop
      !             end if
      !             if (NMO /= size(Rho_mo,dim=1)) then
      !                   call msg("Invalid dimension of matrix Rho_mo")
      !                   error stop
      !             end if
      !             !
      !             ! Bare nuclei hamiltonian (atomic orbital basis, auto2e ordering)
      !             !
      !             allocate(H0_ao(NAO, NAO))
      !             call ints1e_gammcor_H0(H0_ao, AOBasis, System)
      !             !
      !             ! MO coefficients in AO basis (auto2e ordering)
      !             !
      !             allocate(C_ao(NAO, NMO))
      !             FromExternalAO = .true. ! AOs from external program -> AOs in the Auto2e format
      !             TwoIndexTransf = .false. ! Transform only the index p of C(p,k)
      !             call auto2e_interface_AngFuncTransf(C_ao, C, FromExternalAO, TwoIndexTransf, AOBasis, ExternalOrdering)
      !             if (ExternalOrdering == ORBITAL_ORDERING_ORCA) then
      !                   call auto2e_interface_ApplyOrcaPhases_Matrix(C_ao, AOBasis, TwoIndexTransf)
      !             end if
      !             !
      !             ! Bare nuclei hamiltonian (MO basis)
      !             !
      !             allocate(TransfWork(NAO, NMO))
      !             call real_ab(TransfWork, H0_ao, C_ao)
      !             call real_aTb(H0_mo, C_ao, TransfWork)
      !             !
      !             ! Eigenvectors of the density matrix
      !             ! Occupation numbers
      !             !
      !             allocate(Ceig_mo(NMO, NMO))
      !             allocate(OccNum(NMO, 1))
      !             Ceig_mo = Rho_mo
      !             call symmetric_eigenproblem(OccNum(:, 1), Ceig_mo, NMO, .true.)
      !             do p = 1, NMO
      !                   OccNum(p, 1) = max(ZERO, OccNum(p, 1))
      !             end do
      !             allocate(Ceig_ao(NAO, NMO, 1))
      !             call real_ab(Ceig_ao(:, :, 1), C_ao, Ceig_mo)
      !             !
      !             ! 1-RDM in AO basis (auto2e ordering)
      !             !
      !             allocate(Rho_ao(NAO, NAO, 1))
      !             call real_ab(TransfWork, C_ao, Rho_mo)
      !             call real_abT(Rho_ao(:, :, 1), TransfWork, C_ao)
      !             !
      !             ! Coulomb (J) and exchange (K) matrices
      !             ! (AO basis, auto2e ordering)
      !             !
      !             NOcc(1) = NMO
      !             NOcc(2) = 0
      !             allocate(JK_ao(NAO, NAO, 1))
      !             if ((.not. present(J_mo)) .or. (.not. present(K_mo))) then
      !                   CoulContrib = .true.
      !                   ExchContrib = .true.
      !                   call chf_JK(JK_ao, Rho_ao, Ceig_ao, Rkpq, CholeskyVecs, KScal, OccNum, NOcc, AOBasis, &
      !                         MaxBufferDimMB, CoulContrib, ExchContrib)
      !             else
      !                   allocate(J_ao(NAO, NAO, 1))
      !                   allocate(K_ao(NAO, NAO, 1))
      !                   CoulContrib = .true.
      !                   ExchContrib = .false.
      !                   call chf_JK(J_ao, Rho_ao, Ceig_ao, &
      !                         CholeskyVecs%Rkpq, &
      !                         CholeskyVecs%Chol2Data, &
      !                         KScal, OccNum, NOcc, AOBasis, &
      !                         MaxBufferDimMB, CoulContrib, ExchContrib)
      !                   CoulContrib = .false.
      !                   ExchContrib = .true.
      !                   call chf_JK(K_ao, Rho_ao, Ceig_ao, &
      !                         CholeskyVecs%Rkpq, &
      !                         CholeskyVecs%Chol2Data, &
      !                         KScal, OccNum, NOcc, AOBasis, &
      !                         MaxBufferDimMB, CoulContrib, ExchContrib)
      !                   call real_smfill(J_ao(:, :, 1))
      !                   call real_smfill(K_ao(:, :, 1))
      !                   JK_ao = J_ao + K_ao
      !                   call real_ab(TransfWork, J_ao(:, :, 1), C_ao)
      !                   call real_aTb(J_mo, C_ao, TransfWork)
      !                   call real_ab(TransfWork, K_ao(:, :, 1), C_ao)
      !                   call real_aTb(K_mo, C_ao, TransfWork)
      !             end if
      !             call real_smfill(JK_ao(:, :, 1))
      !             call real_ab(TransfWork, JK_ao(:, :, 1), C_ao)
      !             call real_aTb(F_mo, C_ao, TransfWork)
      !             F_mo = F_mo + H0_mo
      !       end associate
      !       call boys_free()
      ! end subroutine chol_gammcor_F
end module Cholesky_Gammcor
