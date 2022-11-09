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
end module periodic
