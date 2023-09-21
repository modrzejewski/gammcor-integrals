module grid_definitions
      use arithmetic

      implicit none

      integer, parameter :: TESSELLATION_WEIGHT_DENSITY = 1
      integer, parameter :: TESSELLATION_WEIGHT_ORBITALS = 2      

      integer, parameter :: BECKE_PRUNING_NONE = 0
      integer, parameter :: BECKE_PRUNING_POPLE1993 = 1
      integer, parameter :: BECKE_PRUNING_SHERILL2013 = 2
      !
      !                    PRE-DEFINED GRID SETTINGS
      !
      integer, parameter :: BECKE_PARAMS_SG1 = 1
      integer, parameter :: BECKE_PARAMS_MEDIUM = 2
      integer, parameter :: BECKE_PARAMS_FINE = 3
      integer, parameter :: BECKE_PARAMS_XFINE = 4
      integer, parameter :: BECKE_PARAMS_THC = 5
      !
      ! Grid labels
      !
      character(6), dimension(5), parameter :: BECKE_PARAMS_LABEL = [ &
            "SG-1  ", &
            "medium", &
            "fine  ", &
            "xfine ", &
            "THC   " &
            ]
      !
      ! Tensor hypercontraction
      ! ---
      ! Grid used to construct the collocation matrices
      ! in the decomposition of electron repulsion integrals
      !
      ! 1. J. Chem. Theory Comput. 16, 1382 (2020); doi: 10.1021/acs.jctc.9b01205
      !
      integer, parameter :: BECKE_PARAMS_THC_NRadial = 29 ! supporting information for Ref. 1
      integer, parameter :: BECKE_PARAMS_THC_Lmax = 17 ! supporting information for Ref. 1
      ! --------------------------- SG-1 ------------------------------
      integer, parameter :: BECKE_PARAMS_SG1_NRadial = 50
      integer, dimension(5), parameter :: BECKE_PARAMS_SG1_LebedevIdx = [1, 4, 7,  11,  7]
      integer, dimension(5), parameter :: BECKE_PARAMS_SG1_LebedevL = [3, 9, 15, 23, 15]
      !
      ! -------------------------- MEDIUM  ----------------------------
      ! EULER-MACLAURIN / LEBEDEV (96; 14, 74, 302, 146, 302)
      !
      integer, parameter :: BECKE_PARAMS_MEDIUM_NRadial    = 96
      integer, dimension(5), parameter :: BECKE_PARAMS_MEDIUM_LebedevIdx = [2, 7, 9, 14, 9]
      integer, dimension(5), parameter :: BECKE_PARAMS_MEDIUM_LebedevL   = [5, 15, 19, 29, 19]
      !
      ! --------------------------- FINE ------------------------------
      ! EULER-MACLAURIN / LEBEDEV (150; 26, 146, 302, 590, 302)
      ! 
      integer, parameter :: BECKE_PARAMS_FINE_NRadial      = 150
      integer, dimension(5), parameter :: BECKE_PARAMS_FINE_LebedevIdx = [3, 9, 14, 17, 14]
      integer, dimension(5), parameter :: BECKE_PARAMS_FINE_LebedevL   = [7, 19, 29, 41, 29]
      !
      ! ------------------------ EXTRA FINE ---------------------------
      ! EULER-MACLAURIN / LEBEDEV (250; 38, 230, 590, 1202, 590)
      !
      integer, parameter :: BECKE_PARAMS_XFINE_NRadial     = 250
      integer, dimension(5), parameter :: BECKE_PARAMS_XFINE_LebedevIdx = [4, 14, 17, 20, 17]
      integer, dimension(5), parameter :: BECKE_PARAMS_XFINE_LebedevL   = [9, 29, 41, 59, 41]

      integer, parameter :: BECKE_NParamSets = 5
      integer, dimension(BECKE_NParamSets), parameter :: BECKE_PARAMS_NRadial = &
            [BECKE_PARAMS_SG1_NRadial, BECKE_PARAMS_MEDIUM_NRadial, BECKE_PARAMS_FINE_NRadial, &
            BECKE_PARAMS_XFINE_NRadial, BECKE_PARAMS_THC_NRadial]
      integer, dimension(5, BECKE_NParamSets), parameter :: BECKE_PARAMS_LebedevIdx = reshape( &
            [BECKE_PARAMS_SG1_LebedevIdx, &
            BECKE_PARAMS_MEDIUM_LebedevIdx, &
            BECKE_PARAMS_FINE_LebedevIdx, &
            BECKE_PARAMS_XFINE_LebedevIdx, &
            [0, 0, 0, 0, 0]], &
            [5,BECKE_NParamSets] &
            )
      integer, dimension(5, BECKE_NParamSets), parameter :: BECKE_PARAMS_LebedevL = reshape( &
            [BECKE_PARAMS_SG1_LebedevL, &
            BECKE_PARAMS_MEDIUM_LebedevL, &
            BECKE_PARAMS_FINE_LebedevL, &
            BECKE_PARAMS_XFINE_LebedevL, &
            [0, 0, 0, 0, 0]], &
            [5,BECKE_NParamSets] &
            )

      type TBeckeGrid
            real(F64), dimension(:), allocatable :: X
            real(F64), dimension(:), allocatable :: Y
            real(F64), dimension(:), allocatable :: Z
            real(F64), dimension(:), allocatable :: W
            integer :: NPoints = 0
            integer :: NCentroids = 0
      end type TBeckeGrid
end module grid_definitions
