! module gammcor_interface
!       implicit none

! contains

      use basis_sets
      use sys_definitions
      use CholeskyOTF
      use Cholesky_driver: 

      type(TCholeskyVecsOTF) :: CholeskyVecsOTF     
      type(TAOBasis) :: AOBasis
      
      XYZPath = "/home/michalhapka/pr-dmft/cholesky_test/water.xyz"
      BasisSetPath = "/home/michalhapka/pr-dmft/cholesky_test/cc-pVDZ"
      call auto2e_init()
      call cholesky_ao_vectors(CholeskyVecsOTF, AOBasis, XYZPath, BasisSetPath, Accuracy, MaxBufferDimMB)
      !
      ! Transform Cholesky vectors to MO basis
      !
      NA = NBasis
      NB = NBasis
      a0 = 1
      a1 = NBasis
      b0 = 1
      b1 = NBasis
      NCholesky = CholeskyVecsOTF%NVecs
      allocate(Rkab(NCholesky,NA*NB))
      call chol_Rkab_OTF(Rkab, UAux, a0, a1, UAux, b0, b1, MaxBufferDimMB, CholeskyVecsOTF, AOBasis)

      
      subroutine cholesky_ao_vectors(CholeskyVecsOTF, AOBasis, XYZPath, BasisSetPath, Accuracy, MaxBufferDimMB)
            use arithmetic
            use auto2e
            use Cholesky, only: chol_CoulombMatrix, TCholeskyVecs, chol_Rkab_ExternalBinary, chol_MOTransf_TwoStep
            use CholeskyOTF, only: chol_CoulombMatrix_OTF, TCholeskyVecsOTF, chol_MOTransf_TwoStep_OTF
            use Cholesky_driver
            use basis_sets
            use sys_definitions
            use chol_definitions
            
            implicit none


            type(TCholeskyVecsOTF), intent(out) :: CholeskyVecsOTF
            type(TAOBasis), intent(out)         :: AOBasis
            character(*), intent(in)            :: XYZPath
            character(*), intent(in)            :: BasisSetPath
            integer, intent(in)                 :: Accuracy
            integer, intent(in)                 :: MaxBufferDimMB

            type(TSystem) :: System
            logical, parameter :: SpherAO = .true.
            !
            ! Initialize the two-electron intergrals library
            !
            call auto2e_init()
            !
            ! Read the XYZ coordinates and atom types
            !
            call sys_Read_XYZ(System, XYZPath)
            !
            ! Read the basis set parameters from an EMSL text file
            ! (GAMESS-US format, no need for any edits, just download it straight from the website)
            !
            call basis_NewAOBasis(AOBasis, System, BasisSetPath, SpherAO)
            !
            ! Compute Cholesky vectors in AO basis
            !
            call chol_Rkpq_OTF(CholeskyVecsOTF, AOBasis, Accuracy)
      end subroutine cholesky_ao_vectors
end module gammcor_interface
