module TwoStepCholesky_definitions
      use arithmetic
      
      implicit none
      !
      ! Orbital pair subsets      
      !
      integer, parameter :: CHOL2_BASE = 1        ! All orbital pairs accepted after the initial Schwarz prescreening
      integer, parameter :: CHOL2_CANDIDATES = 2  ! Pivot candidates in a macroiteration of the pivot finding algorithm
      integer, parameter :: CHOL2_BATCH = 3       ! Subset of pivot candidates in a microiteration of the pivot finding algorithm
      !
      ! Orbital pair storage modes
      ! used in the ShellPairLoc array
      !
      integer, parameter :: CHOL2_FULL_STORAGE = 1        ! All orbital pairs of the BASE subset are stored
      integer, parameter :: CHOL2_SUBSET_STORAGE = 2      ! Orbital pairs of the BASE subset are divided into subsets
      integer, parameter :: CHOL2_SUBSET_INDEX = 3        ! of size at most equal to MaxBlockDim. Used for parallelization.
      integer, parameter :: CHOL2_COMPRESSED_STORAGE = 4  ! Only the orbital pairs of the CANDIDATES subset are stored

      type TChol2Params
            !
            ! Threshold for screening diagonal elements,
            ! Eq. 17 in J. Chem. Phys. 150, 194112 (2019);
            ! doi: 10.1063/1.5083802
            !
            real(F64) :: CholeskyTauThresh = 1.0E-7_F64
            !
            ! Range separation parameter defining the screened
            ! interaction potential Erf(Omega*r12)/r12:
            !
            ! Kappa=1/Omega**2
            !
            ! Kappa=0 means that the full range Coulomb
            ! potential is used.
            !
            real(F64) :: Kappa = 0.0_F64
            !
            ! Maximum size of a block of the matrix W(pq,rs). W is partitioned into blocks,
            ! which are then passed to the matrix multiplication subroutine to compute
            ! matrix products such as RW(:, rs) <- R(:, pq) * W(pq, rs).
            !
            ! The size of the block should be big enough
            ! to fit into the range where the matrix multiplication subroutine achieves
            ! high floating point operations per second. In addition, larger MaxBlockDim
            ! yields better shared-memory parallel computation of W(pq,rs) within each block.
            ! 
            integer :: MaxBlockDim = 4000
            !
            ! Block size of Cholesky vectors used during the two-step Cholesky
            ! factorization
            !
            integer :: CholVecsBlock = 500
      end type TChol2Params

      
      type TChol2Vecs
            real(F64), dimension(:, :), allocatable :: Inv_L
            integer :: NVecs
            integer :: MaxSubsetDim
            integer, dimension(2) :: NSubsets
            integer :: NOrbPairs
            integer :: NShellPairs
            integer, dimension(:, :), allocatable :: ShellPairs
            integer, dimension(:, :), allocatable :: ShellPairLoc
            integer, dimension(:), allocatable :: ShellPairDim
            integer, dimension(:), allocatable :: SubsetDim
            integer, dimension(:, :), allocatable :: SubsetBounds
            integer, dimension(:), allocatable :: Pivots
            integer, dimension(:), allocatable :: PivotShellPairs
            integer, dimension(:), allocatable :: PivotShellPairLoc
            integer, dimension(:), allocatable :: PivotShellPairDim
            integer, dimension(:), allocatable :: PivotOrbPairs
            integer :: NPivotShellPairs
      end type TChol2Vecs

      type TCompressedCholVecs
            real(F64), dimension(:, :), allocatable :: L
      end type TCompressedCholVecs
end module TwoStepCholesky_definitions
