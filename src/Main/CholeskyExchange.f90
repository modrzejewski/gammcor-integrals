module CholeskyExchange
      use arithmetic
      use linalg
      use CholeskyOTF
      use basis_sets

      implicit none

contains

      subroutine chf_EfficientBufferDim(EfficientBufferDim, NCholesky, NAO)
            !
            ! Compute the preffered buffer dimension for the single-index
            ! Cholesky vector transformation. The preferred dimension depends
            ! depends on the loop parallelization:
            ! larger buffer -> better work partition between concurrent threads
            !
            integer(I64), intent(out) :: EfficientBufferDim
            integer, intent(in)       :: NCholesky
            integer, intent(in)       :: NAO

            integer :: NThreads
            
            !$omp parallel default(shared)
            !$omp master
            NThreads = 1
            !$ NThreads = omp_get_num_threads()
            !$omp end master
            !$omp end parallel

            EfficientBufferDim = 2 * int(NThreads,I64) * int(NCholesky,I64) * int(NAO,I64)
      end subroutine chf_EfficientBufferDim
      

      pure function chf_p_ge_q2pq(p, q, m)
            !
            ! Compute compound 2-index, (pq), assuming that p is always
            ! greater or equal q. Min. index: 1, max index: M.
            ! Example: P and Q enumerate, respectively, rows and columns
            ! of a symmetric matrix. The compound 2-index enumerates 
            ! consecutive elements of the lower triangle of the 5x5 matrix.
            !
            ! M = 5
            !             Q 
            !    | 1                |
            !    | 2  6             |
            ! P  | 3  7  10         |
            !    | 4  8  11  13     |
            !    | 5  9  12  14  15 |
            !
            integer             :: chf_p_ge_q2pq
            integer, intent(in) :: p
            integer, intent(in) :: q
            integer, intent(in) :: m

            integer :: i1, i2

            i1 = ((2 * m - q + 2) * (q - 1)) / 2
            i2 = p - q + 1
            chf_p_ge_q2pq = i1 + i2
      end function chf_p_ge_q2pq


      pure subroutine chf_pq2p_ge_q(p, q, pq, n)
            !
            ! Decode a lower-triangle compound index into individual
            ! indices:
            ! PQ -> (P, Q)
            ! Assumptions:
            ! 0) P = 1, 2, ..., N,
            !    Q = 1, 2, ..., N,
            ! 1) P >= Q (diagonal indices admissible)
            !
            ! An example of how this algorithm traverses an N=3 triangle:
            !
            ! 1
            ! 2 5
            ! 3 6 4
            !
            integer, intent(out) :: p
            integer, intent(out) :: q
            integer, intent(in)  :: pq
            integer, intent(in)  :: n

            integer :: q_base
            integer :: v
            integer :: interval1
            integer :: in1, in2
            !
            ! pq = (q_base - 1) * (n + 1) + v
            !
            q_base = (pq - 1) / (n + 1) + 1
            v = pq - (n + 1) * (q_base - 1)
            !
            ! Decide if v is in interval_1 or interval_2:
            ! in1 == 1 and in2 == 0 if v <= INTERVAL1
            ! in1 == 0 and in2 == 1 if v > INTERVAL1
            !
            interval1 = n - q_base + 1
            in2 = v / (interval1 + 1)
            !
            ! 1 -> 0, 0 -> 1
            !
            in1 = ieor((in2), 1)

            p = in1 * (q_base + v - 1) + in2 * (v - interval1 + n - q_base)          
            q = in1 * q_base + in2 * interval1
      end subroutine chf_pq2p_ge_q


      subroutine chf_Kpq_Cip(Cip, Cpi, i0k, i1k, NAO)
            real(F64), dimension(i0k:i1k, NAO), intent(out) :: Cip
            real(F64), dimension(:, :), intent(in)          :: Cpi
            integer, intent(in)                             :: i0k
            integer, intent(in)                             :: i1k
            integer, intent(in)                             :: NAO

            Cip = transpose(Cpi(:, i0k:i1k))
      end subroutine chf_Kpq_Cip


      subroutine chf_Kpq_Transf(Wkip, Rkpq, NCholesky, Cip, Nk, &
            SubsetBounds, ShellPairs, ShellPairLoc, ShellLoc, ShellParamsIdx, &
            NAngFunc, NAO)

            real(F64), dimension(NCholesky, Nk, NAO), intent(inout)        :: Wkip
            real(F64), dimension(:, :), intent(in)                         :: Rkpq
            integer, intent(in)                                            :: NCholesky
            real(F64), dimension(Nk, NAO), intent(in)                      :: Cip
            integer, intent(in)                                            :: Nk
            integer, dimension(2), intent(in)                              :: SubsetBounds
            integer, dimension(:, :), intent(in)                           :: ShellPairs
            integer, dimension(:, :), intent(in)                           :: ShellPairLoc
            integer, dimension(:), intent(in)                              :: ShellLoc
            integer, dimension(:), intent(in)                              :: ShellParamsIdx
            integer, dimension(:), intent(in)                              :: NAngFunc
            integer, intent(in)                                            :: NAO

            integer :: ShAB, LocAB
            integer :: ShA, ShellParamsA, Na, LocA
            integer :: ShB, ShellParamsB, Nb, LocB
            integer :: a, b, p, q, pq
            integer :: i, kappa

            do ShAB = SubsetBounds(1), SubsetBounds(2)
                  LocAB = ShellPairLoc(SUBSET_STORAGE, ShAB)

                  ShA = ShellPairs(1, ShAB)
                  ShellParamsA = ShellParamsIdx(ShA)
                  Na = NAngFunc(ShellParamsA)
                  LocA = ShellLoc(ShA)

                  ShB = ShellPairs(2, ShAB)
                  ShellParamsB = ShellParamsIdx(ShB)
                  Nb = NAngFunc(ShellParamsB)
                  LocB = ShellLoc(ShB)
                  !
                  ! W(:, :, p) <- W(:, :, p) + R(1:NVecs, pq) * C(q, i0k:i1k)**T
                  ! W(:, :, q) <- W(:, :, q) + R(1:NVecs, pq) * C(p, i0k:i1k)**T
                  !
                  if (ShA /= ShB) then
                        do b = 1, Nb
                              q = LocB + b - 1
                              !$omp parallel do collapse(3) &
                              !$omp private(a, i, p, pq, kappa) &
                              !$omp default(shared)
                              do a = 1, Na
                                    do i = 1, Nk
                                          do kappa = 1, NCholesky
                                                p = LocA + a - 1
                                                pq = LocAB + a-1 + Na * (b - 1)
                                                Wkip(kappa, i, p) = Wkip(kappa, i, p) + Rkpq(kappa, pq) * Cip(i, q)
                                          end do
                                    end do
                              end do
                              !$omp end parallel do
                        end do
                        
                        do a = 1, Na
                              p = LocA + a - 1
                              !$omp parallel do collapse(3) &
                              !$omp private(b, i, q, pq, kappa) &
                              !$omp default(shared)
                              do b = 1, Nb
                                    do i = 1, Nk
                                          do kappa = 1, NCholesky
                                                q = LocB + b - 1
                                                pq = LocAB + a-1 + Na * (b - 1)
                                                Wkip(kappa, i, q) = Wkip(kappa, i, q) + Rkpq(kappa, pq) * Cip(i, p)
                                          end do
                                    end do
                              end do
                              !$omp end parallel do
                        end do
                  else
                        do b = 1, Na
                              q = LocB + b - 1
                              !$omp parallel do collapse(3) &
                              !$omp private(a, i, p, pq, kappa) &
                              !$omp default(shared)
                              do a = b, Na
                                    do i = 1, Nk
                                          do kappa = 1, NCholesky
                                                p = LocA + a - 1
                                                pq = LocAB + chf_p_ge_q2pq(a, b, Na) - 1
                                                Wkip(kappa, i, p) = Wkip(kappa, i, p) + Rkpq(kappa, pq) * Cip(i, q)
                                          end do
                                    end do
                              end do
                              !$omp end parallel do
                        end do

                        do a = 1, Na
                              p = LocA + a - 1
                              !$omp parallel do collapse(3) &
                              !$omp private(b, i, q, pq, kappa) &
                              !$omp default(shared)
                              do b = 1, a - 1
                                    do i = 1, Nk
                                          do kappa = 1, NCholesky
                                                q = LocB + b - 1
                                                pq = LocAB + chf_p_ge_q2pq(a, b, Na) - 1                                          
                                                Wkip(kappa, i, q) = Wkip(kappa, i, q) + Rkpq(kappa, pq) * Cip(i, p)
                                          end do
                                    end do
                              end do
                              !$omp end parallel do
                        end do
                  end if
            end do
      end subroutine chf_Kpq_Transf


      subroutine chf_CopyKpqBlock(Kpq, KpqBlock, p0, p1, q0, q1, Np, Nq)
            real(F64), dimension(:, :), intent(inout) :: Kpq
            real(F64), dimension(Np, Nq), intent(in)  :: KpqBlock
            integer, intent(in)                       :: p0, p1
            integer, intent(in)                       :: q0, q1
            integer, intent(in)                       :: Np, Nq

            Kpq(p0:p1, q0:q1) = Kpq(p0:p1, q0:q1) + KpqBlock
      end subroutine chf_CopyKpqBlock
      

      subroutine chf_RkipRkiq(Kpq, Wkip, Ni, NCholesky, TargetBlockDim, NAO, KScal)
            real(F64), dimension(:, :), intent(inout) :: Kpq
            real(F64), dimension(NCholesky*Ni, NAO)   :: Wkip
            integer, intent(in)                       :: Ni
            integer, intent(in)                       :: NCholesky
            integer, intent(in)                       :: TargetBlockDim
            integer, intent(in)                       :: NAO
            real(F64), intent(in)                     :: KScal
            
            integer :: mn, m, n
            integer :: p0, p1, q0, q1, Np, Nq
            real(F64), dimension(:), allocatable :: KpqBlock
            integer :: ThisImage, NImages
            integer :: NBlocks

            ! NImages = num_images()
            ! ThisImage = this_image()
            NImages = 1
            ThisImage = 1
            NBlocks = max(1, NAO / TargetBlockDim)
            if (NBlocks * TargetBlockDim < NAO) then
                  NBlocks = NBlocks + 1
            end if
            allocate(KpqBlock(TargetBlockDim**2))
            KpqBlock = ZERO
            do mn = ThisImage, (NBlocks * (NBlocks + 1)) / 2, NImages
                  call chf_pq2p_ge_q(m, n, mn, NBlocks)
                  p0 = TargetBlockDim * (m - 1) + 1
                  p1 = min(NAO, TargetBlockDim * m)
                  q0 = TargetBlockDim * (n - 1) + 1
                  q1 = min(NAO, TargetBlockDim * n)
                  Np = p1 - p0 + 1
                  Nq = q1 - q0 + 1
                  !
                  ! Kpq = Kpq + Sum(ki) Rkip*Rkiq (multiplied by -KScal)
                  !
                  call linalg_aTb_x(KpqBlock, Np, Wkip(:, p0:p1), NCholesky*Ni, &
                        Wkip(:, q0:q1), NCholesky*Ni, Np, Nq, NCholesky*Ni, -KScal, ZERO)
                  call chf_CopyKpqBlock(Kpq, KpqBlock, p0, p1, q0, q1, Np, Nq)
            end do
      end subroutine chf_RkipRkiq


      subroutine chf_NiWkip(Wkip, NCholesky, Ni, NAO, OccNum)
            real(F64), dimension(NCholesky, Ni, NAO), intent(inout) :: Wkip
            integer, intent(in)                                     :: NCholesky
            integer, intent(in)                                     :: Ni
            integer, intent(in)                                     :: NAO
            real(F64), dimension(Ni), intent(in)                    :: OccNum

            integer :: p, i

            !$omp parallel do collapse(2) &
            !$omp private(p,i) &
            !$omp default(shared)
            do p = 1, NAO
                  do i = 1, Ni
                        Wkip(:, i, p) = Sqrt(OccNum(i)) * Wkip(:, i, p)
                  end do
            end do
            !$omp end parallel do
      end subroutine chf_NiWkip
      
      
      subroutine chf_Kpq(Kpq, Rkpq, Cpi, NAO, NCholesky, MaxBufferDim, TargetBlockDim, &
            SubsetBounds, ShellPairs, ShellPairLoc, ShellLoc, ShellParamsIdx, &
            NAngFunc, NSubsets, KScal, OccNum)

            real(F64), dimension(:, :), intent(inout)                  :: Kpq
            real(F64), dimension(:, :, :), intent(in)                  :: Rkpq
            real(F64), dimension(:, :), intent(in)                     :: Cpi
            integer, intent(in)                                        :: NAO
            integer, intent(in)                                        :: NCholesky
            integer(I64), intent(in)                                   :: MaxBufferDim
            integer, intent(in)                                        :: TargetBlockDim
            integer, dimension(:, :), intent(in)                       :: SubsetBounds
            integer, dimension(:, :), intent(in)                       :: ShellPairs
            integer, dimension(:, :), intent(in)                       :: ShellPairLoc
            integer, dimension(:), intent(in)                          :: ShellLoc
            integer, dimension(:), intent(in)                          :: ShellParamsIdx
            integer, dimension(:), intent(in)                          :: NAngFunc
            integer, dimension(2), intent(in)                          :: NSubsets
            real(F64), intent(in)                                      :: KScal
            real(F64), dimension(:), intent(in)                        :: OccNum

            integer :: Ni
            integer :: MaxNik, NPasses, Nik
            integer :: k, i0k, i1k
            integer :: i0, i1
            real(F64), dimension(:, :), allocatable :: Cip
            real(F64), dimension(:), allocatable :: Wkip
            integer :: X, Y, SubsetIdx
            integer :: ThisImage

            ! ThisImage = this_image()
            ThisImage = 1
            Y = ThisImage
            i0 = 1
            i1 = size(Cpi, dim=2)
            Ni = i1 - i0 + 1
            MaxNik = min(Ni, int(MaxBufferDim/(NCholesky*NAO)))
            if (.not. MaxNik > 0) then
                  call msg("Buffer size too small for Rkpq->Rkiq")
                  error stop
            end if
            NPasses = NI / MaxNik
            if (modulo(NI, MaxNik) > 0) then
                  NPasses = NPasses + 1
            end if
            allocate(Cip(MaxNik, NAO))
            allocate(Wkip(NCholesky*MaxNik*NAO))
            do k = 1, NPasses
                  !
                  ! First index transformation R(:, pq) -> R(:, iq) for i=i0k...i1k
                  !
                  Wkip = ZERO
                  i0k = (k - 1) * MaxNik + i0
                  i1k = min(k * MaxNik, i1)
                  call chf_Kpq_Cip(Cip, Cpi, i0k, i1k, NAO)
                  Nik = i1k - i0k + 1
                  do X = 1, NSubsets(1)
                        SubsetIdx = X + (Y - 1) * NSubsets(1)
                        call chf_Kpq_Transf(Wkip, Rkpq(:, :, X), NCholesky, Cip, Nik, &
                              SubsetBounds(:, SubsetIdx), ShellPairs, ShellPairLoc, ShellLoc, &
                              ShellParamsIdx, NAngFunc, NAO)
                  end do
                  !
                  ! Sum over concurrent images
                  !
                  ! call co_sum(Wkip)
                  !
                  ! Kpq = Kpq + Sum(ki) OccNum(i)*Rkip*Rkiq
                  !
                  call chf_NiWkip(Wkip, NCholesky, Nik, NAO, OccNum(i0k:i1k))
                  call chf_RkipRkiq(Kpq, Wkip, Nik, NCholesky, TargetBlockDim, NAO, KScal)
            end do
      end subroutine chf_Kpq


      subroutine chf_K(F, C_ao, OccNum, NOcc, AOBasis, CholeskyVecs, &
            MaxBufferDimMB, TargetBlockDim, KScal)
            !
            ! Compute the exchange component of the Fock matrix F(p,q)
            ! using a single-index transformation of the Cholesky
            ! vectors in AO basis (Rkpq->Rkiq).
            !
            ! The matrix K is defined as
            ! K(p,q) = Sum(i) (pi|iq) = -KScal*Sum(ki) OccNum(i)*Rkip*Rkiq
            !
            ! Closed-shell systems:
            ! F(p,q) <- F(p,q) + K(p,q)
            !
            ! Open-shell systems:
            ! F(p,q,s) <- F(p,q) + K(p,q,s)
            !
            ! Only the leftmost NOcc columns of C_ao are referenced.
            ! 
            !
            real(F64), dimension(:, :, :), intent(inout) :: F
            real(F64), dimension(:, :, :), intent(in)    :: C_ao
            real(F64), dimension(:, :), intent(in)       :: OccNum
            integer, dimension(:), intent(in)            :: NOcc
            type(TAOBasis), intent(in)                   :: AOBasis
            type(TCholeskyVecsOTF), intent(in)           :: CholeskyVecs
            integer, intent(in)                          :: MaxBufferDimMB
            integer, intent(in)                          :: TargetBlockDim
            real(F64), intent(in)                        :: KScal

            integer :: NSpins, s
            integer(I64) :: MaxBufferDim

            associate ( &
                  Rkpq => CholeskyVecs%R, &
                  NCholesky => CholeskyVecs%NVecs, &
                  ShellPairs => CholeskyVecs%ShellPairs, &
                  ShellPairLoc => CholeskyVecs%ShellPairLoc, &
                  ShellPairDim => CholeskyVecs%ShellPairDim, &
                  SubsetBounds => CholeskyVecs%SubsetBounds, &
                  SubsetDim => CholeskyVecs%SubsetDim, &
                  NSubsets => CholeskyVecs%NSubsets, &
                  ShellParamsIdx => AOBasis%ShellParamsIdx, &
                  NAOCart => AOBasis%NAOCart, &
                  NAOSpher => AOBasis%NAOSpher, &
                  NAngFuncCart => AOBasis%NAngFuncCart, &
                  NAngFuncSpher => AOBasis%NAngFuncSpher, &
                  ShellLocCart => AOBasis%ShellLocCart, &
                  ShellLocSpher => AOBasis%ShellLocSpher, &
                  SpherAO => AOBasis%SpherAO &
                  )
                  !
                  ! Maximum scratch array length for the single-index transformation
                  ! of the Cholesky vectors
                  !
                  MaxBufferDim = (int(MaxBufferDimMB,I64)*1024*1024)/(storage_size(Rkpq)/8)
                  NSpins = size(C_ao, dim=3)
                  do s = 1, NSpins
                        if (NOcc(s) > 0) then
                              if (SpherAO) then
                                    call chf_Kpq(F(:, :, s), Rkpq, C_ao(:, 1:NOcc(s), s), NAOSpher, NCholesky, MaxBufferDim, &
                                          TargetBlockDim, SubsetBounds, ShellPairs, ShellPairLoc, ShellLocSpher, ShellParamsIdx, &
                                          NAngFuncSpher, NSubsets, KScal, OccNum(1:NOcc(s),s))
                              else
                                    call chf_Kpq(F(:, :, s), Rkpq, C_ao(:, 1:NOcc(s), s), NAOCart, NCholesky, MaxBufferDim, &
                                          TargetBlockDim, SubsetBounds, ShellPairs, ShellPairLoc, ShellLocCart, ShellParamsIdx, &
                                          NAngFuncCart, NSubsets, KScal, OccNum(1:NOcc(s),s))
                              end if
                        end if
                  end do
            end associate
      end subroutine chf_K
end module CholeskyExchange
