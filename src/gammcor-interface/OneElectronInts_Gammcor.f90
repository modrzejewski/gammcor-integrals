module OneElectronInts_Gammcor
      use arithmetic
      use OneElectronInts
      use sys_definitions
      use basis_sets
      use real_linalg
      use Auto2eInterface

      implicit none

contains

      subroutine ints1e_gammcor_H0(H0_ao, AOBasis, System)
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
                        call ints1e_gammcor_H0_cao(H0_cao, AOBasis, System)
                        call real_smfill(H0_cao)
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
                        call ints1e_gammcor_H0_cao(H0_ao, AOBasis, System)
                        call real_smfill(H0_ao)
                  end if
            end associate
      end subroutine ints1e_gammcor_H0
      

      subroutine ints1e_gammcor_H0_cao(H0_cao, AOBasis, System)
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
      end subroutine ints1e_gammcor_H0_cao

      
      subroutine ints1e_gammcor_H0_mo(H0_mo, C, AOBasis, System, ExternalOrdering)            
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
                  call ints1e_gammcor_H0(H0_ao, AOBasis, System)
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
                  call real_ab(TransfWork, H0_ao, C_ao)
                  call real_aTb(H0_mo, C_ao, TransfWork)
            end associate
            call boys_free()
      end subroutine ints1e_gammcor_H0_mo


      subroutine ints1e_gammcor_H0_extao(H0_extao, AOBasis, System, ExternalOrdering)            
            real(F64), dimension(:, :), intent(out) :: H0_extao
            type(TAOBasis), intent(in)              :: AOBasis
            type(TSystem), intent(in)               :: System
            integer, intent(in)                     :: ExternalOrdering

            real(F64), dimension(:, :), allocatable :: H0_ao
            integer :: NAO
            logical :: FromExternalAO, TwoIndexTransf

            call auto2e_init()
            call boys_init(4*AUTO2E_MAXL)
            associate (NAOCart => AOBasis%NAOCart, &
                  NAOSpher => AOBasis%NAOSpher, &
                  SpherAO => AOBasis%SpherAO)

                  NAO = size(H0_extao, dim=1)
                  if ((SpherAO .and. NAO /= NAOSpher) .or. (.not.SpherAO .and. NAO /= NAOCart)) then
                        call msg("Invalid dimension of matrix H0_extao")
                        error stop
                  end if
                  !
                  ! Bare nuclei hamiltonian (atomic orbital basis, auto2e ordering)
                  !
                  allocate(H0_ao(NAO, NAO))
                  call ints1e_gammcor_H0(H0_ao, AOBasis, System)
                  !
                  ! MO coefficients in AO basis (auto2e ordering)
                  !
                  FromExternalAO = .false.
                  TwoIndexTransf = .true.
                  call auto2e_interface_AngFuncTransf(H0_extao, H0_ao, FromExternalAO, TwoIndexTransf, AOBasis, ExternalOrdering)
                  if (ExternalOrdering == ORBITAL_ORDERING_ORCA) then
                        call auto2e_interface_ApplyOrcaPhases_Matrix(H0_extao, AOBasis, TwoIndexTransf)
                  end if
            end associate
            call boys_free()
      end subroutine ints1e_gammcor_H0_extao
end module OneElectronInts_Gammcor
