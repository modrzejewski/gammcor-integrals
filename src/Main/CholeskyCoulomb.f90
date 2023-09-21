module CholeskyCoulomb
      use arithmetic
      use linalg
      use CholeskyOTF
      use basis_sets
      
      implicit none

contains

      subroutine coul_Cab(C, Rho, Na, Nb, LocA, LocB, NAO)
            !
            ! ShA /= ShB
            ! ShC /= ShD
            !
            real(F64), dimension(Na, Nb), intent(out)  :: C
            real(F64), dimension(NAO, NAO), intent(in) :: Rho
            integer, intent(in)                        :: Na
            integer, intent(in)                        :: Nb
            integer, intent(in)                        :: LocA
            integer, intent(in)                        :: LocB
            integer, intent(in)                        :: NAO
            
            integer :: pp, qq
            integer :: p, q
            integer :: OffsetP, OffsetQ

            OffsetP = LocA - 1
            OffsetQ = LocB - 1
            do q = 1, Nb
                  do p = 1, Na
                        pp = OffsetP + p
                        qq = OffsetQ + q
                        C(p, q) = TWO * Rho(pp, qq)
                  end do
            end do
      end subroutine coul_Cab


      subroutine coul_Jab(J, RQ, Na, Nb, LocA, LocB, NAO)
            !
            ! ShA /= ShB
            ! ShC /= ShD
            !
            real(F64), dimension(NAO, NAO), intent(inout) :: J
            real(F64), dimension(Na, Nb), intent(in)      :: RQ
            integer, intent(in)                           :: Na
            integer, intent(in)                           :: Nb
            integer, intent(in)                           :: LocA
            integer, intent(in)                           :: LocB
            integer, intent(in)                           :: NAO
            
            integer :: p0, p1, q0, q1

            p0 = LocA
            p1 = LocA + Na - 1
            q0 = LocB
            q1 = LocB + Nb - 1
            J(p0:p1, q0:q1) = RQ
      end subroutine coul_Jab


      subroutine coul_Caa(C, Rho, Na, LocA, NAO)
            !
            ! ShA == ShB
            !
            real(F64), dimension(*), intent(out)       :: C
            real(F64), dimension(NAO, NAO), intent(in) :: Rho
            integer, intent(in)                        :: Na
            integer, intent(in)                        :: LocA
            integer, intent(in)                        :: NAO
            
            integer :: pp, qq
            integer :: p, q, pq
            integer :: OffsetP, OffsetQ

            OffsetP = LocA - 1
            OffsetQ = LocA - 1
            pq = 1
            do q = 1, Na
                  !
                  ! Diagonal
                  !
                  p = q
                  pp = OffsetP + p
                  qq = OffsetQ + q
                  C(pq) = Rho(pp, qq)
                  pq = pq + 1
                  !
                  ! Off-diagonal
                  !
                  do p = q + 1, Na
                        pp = OffsetP + p
                        qq = OffsetQ + q
                        C(pq) = TWO * Rho(pp, qq)
                        pq = pq + 1
                  end do
            end do
      end subroutine coul_Caa


      subroutine coul_Jaa(J, RQ, Na, LocA, NAO)
            !
            ! ShA == ShB
            !
            real(F64), dimension(NAO, NAO), intent(inout) :: J
            real(F64), dimension(*), intent(in)           :: RQ
            integer, intent(in)                           :: Na
            integer, intent(in)                           :: LocA
            integer, intent(in)                           :: NAO
            
            integer :: pp, qq
            integer :: p, q, pq
            integer :: OffsetP, OffsetQ

            OffsetP = LocA - 1
            OffsetQ = LocA - 1
            pq = 1
            do q = 1, Na
                  !
                  ! Diagonal
                  !
                  p = q
                  pp = OffsetP + p
                  qq = OffsetQ + q
                  J(pp, qq) = RQ(pq)
                  pq = pq + 1
                  !
                  ! Off-diagonal
                  !
                  do p = q + 1, Na
                        pp = OffsetP + p
                        qq = OffsetQ + q
                        J(pp, qq) = RQ(pq)
                        pq = pq + 1
                  end do
            end do
      end subroutine coul_Jaa
      

      subroutine coul_C(C, SubsetBounds, Rho, Npq, ShellPairs, ShellPairLoc, &
            ShellPairDim, ShellLoc, ShellParamsIdx, NAngFunc, NAO)
            
            real(F64), dimension(Npq), intent(out)             :: C
            integer, dimension(2), intent(in)                  :: SubsetBounds
            real(F64), dimension(NAO, NAO), intent(in)         :: Rho
            integer, intent(in)                                :: Npq
            integer, dimension(:, :), intent(in)               :: ShellPairs
            integer, dimension(:, :), intent(in)               :: ShellPairLoc
            integer, dimension(:), intent(in)                  :: ShellPairDim
            integer, dimension(:), intent(in)                  :: ShellLoc
            integer, dimension(:), intent(in)                  :: ShellParamsIdx
            integer, dimension(:), intent(in)                  :: NAngFunc
            integer, intent(in)                                :: NAO

            integer :: ShAB, LocAB, Nab, pq0, pq1
            integer :: ShA, ShellParamsA, Na, LocA
            integer :: ShB, ShellParamsB, Nb, LocB

            !$omp parallel do schedule(guided) &
            !$omp default(shared) &
            !$omp private(ShA, Na, ShellParamsA, LocA) &
            !$omp private(ShB, Nb, ShellParamsB, LocB) &
            !$omp private(Nab, LocAB) &
            !$omp private(pq0, pq1) &
            !$omp private(ShAB)
            do ShAB = SubsetBounds(1), SubsetBounds(2)
                  LocAB = ShellPairLoc(SUBSET_STORAGE, ShAB)
                  Nab = ShellPairDim(ShAB)
                  pq0 = LocAB
                  pq1 = LocAB + Nab - 1
                  
                  ShA = ShellPairs(1, ShAB)
                  ShellParamsA = ShellParamsIdx(ShA)
                  Na = NAngFunc(ShellParamsA)
                  LocA = ShellLoc(ShA)
                  
                  ShB = ShellPairs(2, ShAB)
                  ShellParamsB = ShellParamsIdx(ShB)
                  Nb = NAngFunc(ShellParamsB)
                  LocB = ShellLoc(ShB)

                  if (ShA /= ShB) then
                        call coul_Cab(C(pq0:pq1), Rho, Na, Nb, LocA, LocB, NAO)
                  else
                        call coul_Caa(C(pq0:pq1), Rho, Na, LocA, NAO)
                  end if
            end do
            !$omp end parallel do
      end subroutine coul_C


      subroutine coul_J_3(J, RQ, SubsetBounds, Npq, ShellPairs, ShellPairLoc, &
            ShellPairDim, ShellLoc, ShellParamsIdx, NAngFunc)

            real(F64), dimension(:, :), intent(inout)          :: J
            real(F64), dimension(Npq), intent(in)              :: RQ
            integer, dimension(2), intent(in)                  :: SubsetBounds
            integer, intent(in)                                :: Npq
            integer, dimension(:, :), intent(in)               :: ShellPairs
            integer, dimension(:, :), intent(in)               :: ShellPairLoc
            integer, dimension(:), intent(in)                  :: ShellPairDim
            integer, dimension(:), intent(in)                  :: ShellLoc
            integer, dimension(:), intent(in)                  :: ShellParamsIdx
            integer, dimension(:), intent(in)                  :: NAngFunc

            integer :: ShAB, LocAB, Nab, pq0, pq1
            integer :: ShA, ShellParamsA, Na, LocA
            integer :: ShB, ShellParamsB, Nb, LocB
            integer :: NAO

            NAO = size(J, dim=1)
            !$omp parallel do schedule(guided) &
            !$omp default(shared) &
            !$omp private(ShA, Na, ShellParamsA, LocA) &
            !$omp private(ShB, Nb, ShellParamsB, LocB) &
            !$omp private(Nab, LocAB) &
            !$omp private(pq0, pq1) &
            !$omp private(ShAB)
            do ShAB = SubsetBounds(1), SubsetBounds(2)
                  LocAB = ShellPairLoc(SUBSET_STORAGE, ShAB)
                  Nab = ShellPairDim(ShAB)
                  pq0 = LocAB
                  pq1 = LocAB + Nab - 1
                  
                  ShA = ShellPairs(1, ShAB)
                  ShellParamsA = ShellParamsIdx(ShA)
                  Na = NAngFunc(ShellParamsA)
                  LocA = ShellLoc(ShA)
                  
                  ShB = ShellPairs(2, ShAB)
                  ShellParamsB = ShellParamsIdx(ShB)
                  Nb = NAngFunc(ShellParamsB)
                  LocB = ShellLoc(ShB)

                  if (ShA /= ShB) then
                        call coul_Jab(J, RQ(pq0:pq1), Na, Nb, LocA, LocB, NAO)
                  else
                        call coul_Jaa(J, RQ(pq0:pq1), Na, LocA, NAO)
                  end if
            end do
            !$omp end parallel do
      end subroutine coul_J_3


      subroutine coul_J_2(J, R, Rho, ShellPairs, ShellPairLoc, &
            ShellPairDim, ShellLoc, ShellParamsIdx, SubsetBounds, SubsetDim, NSubsets, &
            NAngFunc, NAO, NVecs)

            real(F64), dimension(:, :), intent(out)    :: J
            real(F64), dimension(:, :, :), intent(in)  :: R
            real(F64), dimension(:, :, :), intent(in)  :: Rho
            integer, dimension(:, :), intent(in)       :: ShellPairs
            integer, dimension(:, :), intent(in)       :: ShellPairLoc
            integer, dimension(:), intent(in)          :: ShellPairDim
            integer, dimension(:), intent(in)          :: ShellLoc
            integer, dimension(:), intent(in)          :: ShellParamsIdx
            integer, dimension(:, :), intent(in)       :: SubsetBounds
            integer, dimension(:), intent(in)          :: SubsetDim
            integer, dimension(2), intent(in)          :: NSubsets
            integer, dimension(:), intent(in)          :: NAngFunc
            integer, intent(in)                        :: NAO
            integer, intent(in)                        :: NVecs

            integer :: NSpins
            integer :: x, y, s
            integer :: SubsetIdx, Npq
            integer :: MaxSubsetDim
            integer :: ThisImage
            integer :: ldR
            real(F64), dimension(:, :), allocatable :: C
            real(F64), dimension(:), allocatable :: Q, RQ

            NSpins = size(Rho, dim=3)
            ! ThisImage = this_image()
            ThisImage = 1
            allocate(Q(NVecs))
            Q = ZERO
            MaxSubsetDim = maxval(SubsetDim)
            ldR = size(R, dim=1)
            allocate(C(MaxSubsetDim, NSpins))
            Y = ThisImage
            do X = 1, NSubsets(1)
                  SubsetIdx = X + (Y - 1) * NSubsets(1)                        
                  Npq = SubsetDim(SubsetIdx)
                  if (Npq > 0) then
                        !
                        ! C stores the same information as Rho, but the matrix elements
                        ! are permuted so that a simple matrix-vector multiplication
                        ! R*C can be done. C is summed over spins.
                        !
                        do s = 1, NSpins
                              call coul_C(C(1:Npq, s), SubsetBounds(:, SubsetIdx), Rho(:, :, s), &
                                    Npq, ShellPairs, ShellPairLoc, ShellPairDim, ShellLoc, &
                                    ShellParamsIdx, NAngFunc, NAO)
                        end do
                        if (NSpins > 1) then
                              C(1:Npq, 1) = C(1:Npq, 1) + C(1:Npq, 2)
                        end if
                        !
                        ! Q <- Q + R*C
                        !
                        ! Indices of R: R(k, pq)
                        ! Indices of Q: Q(k)
                        ! k=Cholesky vector
                        !
                        call linalg_av_x(Q, R(:, :, X), ldR, C(:, 1), NVecs, Npq, ONE, ONE)
                  end if
            end do
            !
            ! Accumulate the vector Q and send back
            ! the sum to every image
            !
            ! call co_sum(Q)
            !
            ! RQ <- R**T*Q
            !
            allocate(RQ(MaxSubsetDim))
            J = ZERO
            Y = ThisImage
            do X = 1, NSubsets(1)
                  SubsetIdx = X + (Y - 1) * NSubsets(1)                        
                  Npq = SubsetDim(SubsetIdx)
                  if (Npq > 0) then
                        call linalg_aTv_x(RQ, R(:, :, X), ldR, Q, NVecs, Npq, ONE, ZERO)
                        call coul_J_3(J, RQ(1:Npq), SubsetBounds(:, SubsetIdx), Npq, ShellPairs, ShellPairLoc, &
                              ShellPairDim, ShellLoc, ShellParamsIdx, NAngFunc)
                  end if
            end do
      end subroutine coul_J_2

      
      subroutine coul_J(J, Rho, AOBasis, CholeskyVecs)
            real(F64), dimension(:, :), intent(out)    :: J
            real(F64), dimension(:, :, :), intent(in)  :: Rho            
            type(TAOBasis), intent(in)                 :: AOBasis
            type(TCholeskyVecsOTF), intent(in)         :: CholeskyVecs

            associate ( &
                  R => CholeskyVecs%R, &
                  NVecs => CholeskyVecs%NVecs, &
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
                  if (SpherAO) then
                        call coul_J_2(J, R, Rho, ShellPairs, ShellPairLoc, ShellPairDim, ShellLocSpher, &
                              ShellParamsIdx, SubsetBounds, SubsetDim, NSubsets, NAngFuncSpher, NAOSpher, NVecs)
                  else
                        call coul_J_2(J, R, Rho, ShellPairs, ShellPairLoc, ShellPairDim, ShellLocCart, &
                              ShellParamsIdx, SubsetBounds, SubsetDim, NSubsets, NAngFuncCart, NAOCart, NVecs)
                  end if
            end associate
      end subroutine coul_J
end module CholeskyCoulomb
