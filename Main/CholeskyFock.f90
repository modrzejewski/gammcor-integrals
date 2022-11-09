module CholeskyFock
      use arithmetic
      use CholeskyCoulomb
      use CholeskyExchange
      use CholeskyOTF
      use OneElectronInts
      use basis_sets
      use sys_definitions
      
      implicit none

contains

      subroutine chf_JK(JK_ao, Rho_ao, C_ao, CholeskyVecs, KScal, OccNum, NOcc, AOBasis, &
            MaxBufferDimMB, CoulContrib, ExchContrib)
            !
            ! Closed-shell case
            ! -----------------
            ! JK(p,q,1) = Sum(rs) (pq|rs)*Rho(rs,1) - KScal * Sum(i) OccNum(i,1) (pi|iq)
            ! (index i denotes occupied orbitals, i=1...NOcc(1))
            !
            ! Open-shell case
            ! ---------------
            ! JK(p,q,sigma) = Sum(rs) (pq|rs)*Rho(rs,sigma) - KScal*Sum(i)OccNum(i,s)(pi|iq) 
            ! (index i denotes occupied sigma spin-orbitals, i=1...NOcc(s))
            !
            real(F64), dimension(:, :, :), intent(out) :: JK_ao
            real(F64), dimension(:, :, :), intent(in)  :: Rho_ao
            real(F64), dimension(:, :, :), intent(in)  :: C_ao
            type(TCholeskyVecsOTF), intent(in)         :: CholeskyVecs
            real(F64), intent(in)                      :: KScal
            real(F64), dimension(:, :), intent(in)     :: OccNum
            integer, dimension(:), intent(in)          :: NOcc
            type(TAOBasis), intent(in)                 :: AOBasis
            integer, intent(in)                        :: MaxBufferDimMB
            logical, intent(in)                        :: CoulContrib
            logical, intent(in)                        :: ExchContrib

            integer :: NSpins
            integer, parameter :: TargetBlockDim = 100
            
            ! ThisImage = this_image()
            
            NSpins = size(Rho_ao, dim=3)
            JK_ao = ZERO
            if (CoulContrib) then
                  call coul_J(JK_ao(:, :, 1), Rho_ao, AOBasis, CholeskyVecs)
                  if (NSpins == 2) JK_ao(:, :, 2) = JK_ao(:, :, 1)
            end if
            if (ExchContrib) then
                  call chf_K(JK_ao, C_ao, OccNum, NOcc, AOBasis, CholeskyVecs, &
                        MaxBufferDimMB, TargetBlockDim, KScal)
            end if
      end subroutine chf_JK


      subroutine chf_H0(H0_cao, AOBasis, System)
            real(F64), dimension(:, :), intent(out) :: H0_cao
            type(TAOBasis), intent(in)              :: AOBasis
            type(TSystem), intent(in)               :: System

            integer :: NAOCart
            real(F64), dimension(:, :), allocatable :: Ts_cao

            NAOCart = AOBasis%NAOCart
            allocate(Ts_cao(NAOCart, NAOCart))
            call ints1e_Kinetic(Ts_cao, AOBasis)
            call ints1e_Coulomb(H0_cao, AOBasis, System)
            H0_cao = H0_cao + Ts_cao
      end subroutine chf_H0
end module CholeskyFock
