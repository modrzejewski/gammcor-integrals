module periodic
      use arithmetic
      use string

      implicit none

      integer, parameter :: KNOWN_ELEMENTS = 86
      
      character(len=2), dimension(KNOWN_ELEMENTS) :: ELNAME_SHORT = [ &
            "H ", "HE", "LI", "BE", "B ", &
            "C ", "N ", "O ", "F ", "NE", &
            "NA", "MG", "AL", "SI", "P ", &
            "S ", "CL", "AR", "K ", "CA", &
            "SC", "TI", "V ", "CR", "MN", &
            "FE", "CO", "NI", "CU", "ZN", &
            "GA", "GE", "AS", "SE", "BR", &
            "KR", "RB", "SR", "Y ", "ZR", &
            "NB", "MO", "TC", "RU", "RH", &
            "PD", "AG", "CD", "IN", "SN", &
            "SB", "TE", "I ", "XE", "CS", &
            "BA", "LA", "CE", "PR", "ND", &
            "PM", "SM", "EU", "GD", "TB", &
            "DY", "HO", "ER", "TM", "YB", &
            "LU", "HF", "TA", "W ", "RE", &
            "OS", "IR", "PT", "AU", "HG", &
            "TL", "PB", "BI", "PO", "AT", &
            "RN"]

      character(len=12), dimension(KNOWN_ELEMENTS) :: ELNAME_LONG = [ &
            "HYDROGEN    ", "HELIUM      ", "LITHIUM     ", "BERYLLIUM   ", "BORON       ", &
            "CARBON      ", "NITROGEN    ", "OXYGEN      ", "FLUORINE    ", "NEON        ", &
            "SODIUM      ", "MAGNESIUM   ", "ALUMINUM    ", "SILICON     ", "PHOSPHORUS  ", &
            "SULFUR      ", "CHLORINE    ", "ARGON       ", "POTASSIUM   ", "CALCIUM     ", &
            "SCANDIUM    ", "TITANIUM    ", "VANADIUM    ", "CHROMIUM    ", "MANGANESE   ", &
            "IRON        ", "COBALT      ", "NICKEL      ", "COPPER      ", "ZINC        ", &
            "GALLIUM     ", "GERMANIUM   ", "ARSENIC     ", "SELENIUM    ", "BROMINE     ", &
            "KRYPTON     ", "RUBIDIUM    ", "STRONTIUM   ", "YTTRIUM     ", "ZIRCONIUM   ", &
            "NIOBIUM     ", "MOLYBDENUM  ", "TECHNETIUM  ", "RUTHENIUM   ", "RHODIUM     ", &
            "PALLADIUM   ", "SILVER      ", "CADMIUM     ", "INDIUM      ", "TIN         ", &
            "ANTIMONY    ", "TELLURIUM   ", "IODINE      ", "XENON       ", "CESIUM      ", &
            "BARIUM      ", "LANTHANUM   ", "CERIUM      ", "PRASEODYMIUM", "NEODYMIUM   ", &
            "PROMETHIUM  ", "SAMARIUM    ", "EUROPIUM    ", "GADOLINIUM  ", "TERBIUM     ", &
            "DYSPROSIUM  ", "HOLMIUM     ", "ERBIUM      ", "THULIUM     ", "YTTERBIUM   ", &
            "LUTETIUM    ", "HAFNIUM     ", "TANTALUM    ", "TUNGSTEN    ", "RHENIUM     ", &
            "OSMIUM      ", "IRIDIUM     ", "PLATINUM    ", "GOLD        ", "MERCURY     ", &
            "THALLIUM    ", "LEAD        ", "BISMUTH     ", "POLONIUM    ", "ASTATINE    ", &
            "RADON       "]

      ! ---------------------------------------------------------------
      ! ATOMIC RADII (USED FOR INTEGRATION ON THE MOLECULAR GRID)
      ! ---------------------------------------------------------------
      ! 1. Gill, P., Johnson, B.G., Pople, J.A., 
      !    A standard grid for density functional calculations,
      !    Chem. Phys. Lett. 209, 506 (1993)
      ! 2. Ghosh, D. C. and Biswas, R., Theoretical Calculation of
      !    Absolute Radii of Atoms and Ions. Part 1. The Atomic Radii,
      !    Int. J. Mol. Sci. 3, 87 (2002)
      ! ---------------------------------------------------------------
      ! All values are written in atomic units. Radii for H...Ar are
      ! taken from Table 1 in [1]. Radii for K...Lr are taken from
      ! Table 1 in [2].
      !
      real(F64), dimension(103), parameter   :: ATOMIC_RADII = [ &
            1.0000d+0, & ! 1  H    
            0.5882d+0, & ! 2  He
            
            3.0769d+0, & ! 3  Li
            2.0513d+0, & ! 4  Be
            1.5385d+0, & ! 5  B
            1.2308d+0, & ! 6  C
            1.0256d+0, & ! 7  N
            0.8791d+0, & ! 8  O
            0.7692d+0, & ! 9  F
            0.6838d+0, & ! 10 Ne
            
            4.0909d+0, & ! 11 Na
            3.1579d+0, & ! 12 Mg
            2.5714d+0, & ! 13 Al
            2.1687d+0, & ! 14 Si
            1.8750d+0, & ! 15 P
            1.6514d+0, & ! 16 S
            1.4754d+0, & ! 17 Cl
            1.3333d+0, & ! 18 Ar
            
            6.7270d+0, & ! 19 K
            5.1928d+0, & ! 20
            4.9333d+0, & ! 21
            4.6980d+0, & ! 22
            4.4847d+0, & ! 23
            4.2899d+0, & ! 24
            4.1109d+0, & ! 25
            3.9467d+0, & ! 26
            3.7946d+0, & ! 27
            3.6542d+0, & ! 28
            3.5240d+0, & ! 29
            3.4023d+0, & ! 30
            2.9599d+0, & ! 31
            2.6195d+0, & ! 32
            2.3491d+0, & ! 33
            2.1295d+0, & ! 34
            1.9474d+0, &
            1.7939d+0, &
            9.0907d+0, &
            7.0175d+0, &
            6.6666d+0, &
            6.3491d+0, &
            6.0605d+0, &
            5.7971d+0, &
            5.5554d+0, &
            5.3332d+0, &
            5.1281d+0, &
            4.9382d+0, &
            4.7619d+0, &
            4.5977d+0, &
            4.0000d+0, &
            3.5398d+0, &
            3.1746d+0, &
            2.8777d+0, &
            2.6316d+0, &
            2.4241d+0, &
            11.4546d+0, &
            8.8417d+0, &
            7.2002d+0, &
            6.0723d+0, &
            5.2497d+0, &
            4.6238d+0, &
            4.1311d+0, &
            3.7333d+0, &
            3.4053d+0, &
            3.1303d+0, &
            2.8966d+0, &
            2.6951d+0, &
            2.5199d+0, &
            2.3661d+0, &
            2.2301d+0, &
            2.1087d+0, &
            1.9999d+0, &
            1.9047d+0, &
            1.8130d+0, &
            1.7319d+0, &
            1.6579d+0, &
            1.5898d+0, &
            1.5462d+0, &
            1.4695d+0, &
            1.4158d+0, &
            5.7894d+0, &
            5.0399d+0, &
            4.4603d+0, &
            4.0000d+0, &
            3.6258d+0, &
            3.3157d+0, &
            3.0546d+0, &
            13.6824d+0, &
            10.5611d+0, &
            10.0327d+0, &
            9.5562d+0, &
            6.9999d+0, &
            6.0806d+0, &
            5.3749d+0, &
            4.4590d+0, &
            4.0676d+0, &
            3.9868d+0, &
            3.4597d+0, &
            3.2191d+0, &
            3.0100d+0, &
            2.8263d+0, &
            2.6638d+0, &
            2.5188d+0, &
            2.4876d+0 ] ! 103 Lr

contains

      function znumber_long(longname)
            ! -----------------------------------------------------
            ! Determine the atomic number, Z, of the element whose
            ! non-abbreviated name is LONGNAME:
            ! "HYDROGEN" -> 1
            ! -----------------------------------------------------
            ! ZNUMBER_LONG - Returned value, the atomic number if
            !                element corresponding to LONGNAME is
            !                found, 0 otherwise.
            ! SHORTNAME    - Non-abbreviated name of chemical element:
            !                "Helium", "Carbon", "Bromine", ...
            !                The result is not case-sensitive
            !
            integer                  :: znumber_long
            character(*), intent(in) :: longname

            character(:), allocatable :: s
            integer :: k

            s = uppercase(longname)
            !
            ! Special case: "PHOSPHOROUS" name is used
            ! in definitions of basis sets
            !
            if (s .eq. "PHOSPHOROUS") s = "PHOSPHORUS"
            znumber_long = 0

            kloop: do k = 1, KNOWN_ELEMENTS
                  if (elname_long(k) .eq. s) then
                        znumber_long = k
                  end if
            end do kloop
      end function znumber_long


      function znumber_short(shortname)
            ! -----------------------------------------------------------
            ! Determine the atomic number, Z, of the element whose
            ! abbreviated name is SHORTNAME:
            ! "H" -> 1
            ! -----------------------------------------------------------
            ! ZNUMBER_SHOR,   Output, the atomic number of the element
            !                 element corresponding to SHORTNAME. If
            !                 the element is not found, 0 is returned.
            ! SHORTNAME     - Abbreviated name of the chemical element:
            !                 "He", "C", "Br", ... 
            !
            integer                  :: znumber_short
            character(*), intent(in) :: shortname

            character(:), allocatable :: s
            integer :: k

            s = uppercase(shortname)
            znumber_short = 0

            kloop: do k = 1, KNOWN_ELEMENTS
                  if (elname_short(k) .eq. s) then
                        znumber_short = k
                  end if
            end do kloop
      end function znumber_short

      
      function periodic_table_row(z)
            integer :: periodic_table_row
            integer, intent(in) :: z
            
            if (z <= 2) then
                  periodic_table_row = 0
            else if (z <= 10) then
                  periodic_table_row = 1
            else if (z <= 18) then
                  periodic_table_row = 2
            else if (z <= 36) then
                  periodic_table_row = 3
            else if (z <= 54) then
                  periodic_table_row = 4
            else if (z <= 86) then
                  periodic_table_row = 5
            else
                  periodic_table_row = 6
            end if
      end function periodic_table_row
end module periodic
