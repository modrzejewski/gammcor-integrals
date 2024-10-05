module math_constants
      use arithmetic

      implicit none

      real(F64), parameter :: pi = 3.141592653589793_F64
      real(F64), parameter :: pi32 = sqrt(pi**3)
      real(F64), parameter :: pi34 = sqrt(pi32)
      real(F64), parameter :: pi52 = sqrt(pi**5)
      real(F64), parameter :: pi12 = sqrt(pi)
      real(F64), parameter :: pi13 = pi**(1.0_F64 / 3.0_F64)
      real(F64), parameter :: pi23 = pi**(2.0_F64 / 3.0_F64)
      real(F64), parameter :: EULER_CONSTANT = 0.57721566490153286_F64

      real(F64), parameter :: zero = 0.0_F64
      real(F64), parameter :: one = 1.0_F64
      real(F64), parameter :: two = 2.0_F64
      real(F64), parameter :: three = 3.0_F64
      real(F64), parameter :: four = 4.0_F64
      real(F64), parameter :: five = 5.0_F64
      real(F64), parameter :: six = 6.0_F64
      real(F64), parameter :: seven = 7.0_F64
      real(F64), parameter :: eight = 8.0_F64
      real(F64), parameter :: nine = 9.0_F64
      real(F64), parameter :: ten = 10.0_F64
      real(F64), parameter :: ELEVEN = 11.0_F64
      real(F64), parameter :: TWELVE = 12.0_F64
      real(F64), parameter :: FIFTEEN = 15.0_F64
      real(F64), parameter :: SIXTEEN = 16.0_F64
      real(F64), parameter :: TWENTY = 20.0_F64
      real(F64), parameter :: FORTY = 40.0_F64

      real(F64), parameter :: five12 = sqrt(FIVE)

      real(F64), parameter :: FRAC12 = one / two
      real(F64), parameter :: FRAC13 = one / three
      real(F64), parameter :: FRAC14 = one / four
      real(F64), parameter :: FRAC15 = ONE / FIVE
      real(F64), parameter :: FRAC16 = one / six
      real(F64), parameter :: FRAC23 = two / three
      real(F64), parameter :: FRAC32 = three / two
      real(F64), parameter :: FRAC34 = three / four
      real(F64), parameter :: FRAC43 = four / three
      real(F64), parameter :: FRAC53 = five / three
      !
      ! Planck constant in J s. Source: CODATA 2010
      !
      real(F64), parameter :: planck = 6.62606957E-34_F64
      !
      ! Planck constant over 2 pi in J s. Source: CODATA 2010
      !
      real(F64), parameter :: planckbar = 1.054571726E-34_F64
      !
      ! Electric constant \epsilon_0 in F m^{-1}. Source: CODATA 2010
      !
      real(F64), parameter :: eps_0 = 8.854187817E-12_F64
      !
      ! Speed of light in vacuum in m s^{-1}. Source: CODATA 2010
      !
      real(F64), parameter :: light_speed = 299792458_F64
      !
      ! Elementary charge in C. Source: CODATA 2010
      !
      real(F64), parameter :: e_charge = 1.602176565E-19_F64
      !
      ! Electron mass in kg. Source: CODATA 2010
      !
      real(F64), parameter :: me = 9.10938291E-31_F64
      !
      ! Atomic unit of velocity a_0 E_h / hbar in m s^{-1}.
      ! Source: CODATA 2010.
      !
      real(F64), parameter :: velocity_atomic = 2.18769126379E06_F64
      !
      ! Bohr unit length in meters. Source: CODATA 2010
      !
      real(F64), parameter :: bohr = 0.52917721092E-10_F64
      !
      ! Binary multiples of the byte
      !
      integer(I64), parameter :: KIBIBYTE = 1024_I64
      integer(I64), parameter :: MEBIBYTE = KIBIBYTE * KIBIBYTE
      integer(I64), parameter :: GIBIBYTE = MEBIBYTE * KIBIBYTE
      integer(I64), parameter :: TEBIBYTE = GIBIBYTE * KIBIBYTE
      integer(I64), parameter :: PEBIBYTE = TEBIBYTE * KIBIBYTE
      integer(I64), parameter :: EXBIBYTE = PEBIBYTE * KIBIBYTE
      !
      ! Decimal multiples of the byte
      !
      integer(I64), parameter :: KILOBYTE = 1000_I64
      integer(I64), parameter :: MEGABYTE = KILOBYTE * KILOBYTE
      integer(I64), parameter :: GIGABYTE = MEGABYTE * KILOBYTE
      integer(I64), parameter :: TERABYTE = GIGABYTE * KILOBYTE
      integer(I64), parameter :: PETABYTE = TERABYTE * KILOBYTE
      integer(I64), parameter :: EXABYTE = PETABYTE * KILOBYTE

contains

      pure function compare(a, b, tol)
            !
            ! Compare two finite-precision numbers
            ! represented in finite-precision arithmetic.
            ! The tolarance for the fuzzy comparison is
            ! expressed in terms of spacing between real
            ! numbers in the finite precision aithmetic.
            ! TOL = 3.0 means that the tolerance is 3.0 * spacing.
            !
            logical               :: compare
            real(F64), intent(in) :: a
            real(F64), intent(in) :: b
            real(F64), intent(in) :: tol

            real(F64) :: sa, sb

            sa = spacing(a)
            sb = spacing(b)
            compare = (abs(a - b) < tol * max(sa, sb))
      end function compare


      pure function tobyte(s, from)
            !
            ! Convert units of information storage
            ! User-specified unit -> byte
            ! 
            integer(I64)             :: tobyte
            real(F64), intent(in)    :: s
            integer(I64), intent(in) :: from

            tobyte = nint(s * real(from, F64), I64)
      end function tobyte


      pure function frombyte(s, to)
            !
            ! Convert units of information storage
            ! Byte -> user-specified unit
            !
            real(F64)                :: frombyte
            integer(I64), intent(in) :: s
            integer(I64), intent(in) :: to
            
            frombyte = real(s, F64) / real(to, F64)
      end function frombyte


      elemental function tobohr(l_angstrom)
            !
            ! Angstrom -> bohr conversion using
            ! CODATA 2010 value of bohr unit:
            ! 1 bohr = 0.52917721092d-10 m
            ! ------------------------------------------
            ! Note: GAMESS value of bohr is  0.52917725.
            ! Use this value when checking the code
            ! against GAMESS.
            !
            real(F64)             :: tobohr
            real(F64), intent(in) :: l_angstrom

            tobohr = l_angstrom / 0.52917721092_F64
      end function tobohr

      
      pure function omega_to_lambda_nm(omega)
            !
            ! Convert angular frequency omega=2pi/T in atomic units
            ! to wavelength in nanometers: lambda=2pi*c/omega
            !
            real(F64)             :: omega_to_lambda_nm
            real(F64), intent(in) :: omega

            real(F64), parameter :: c_au = light_speed / velocity_atomic
            real(F64) :: lambda_au
            
            lambda_au = TWO * pi * c_au / omega
            omega_to_lambda_nm = lambda_au * bohr * 1.0E9_F64
      end function omega_to_lambda_nm
      

      pure function lambda_nm_to_omega(lambda)
            !
            ! Convert wavelength in nanometers to
            ! angular frequency omega in atomic units:
            ! omega = 2pi*c/lambda
            !
            real(F64)             :: lambda_nm_to_omega
            real(F64), intent(in) :: lambda

            real(F64), parameter :: c_au = light_speed / velocity_atomic
            real(F64) :: lambda_au

            lambda_au = lambda * 1.0E-9_F64 / bohr
            lambda_nm_to_omega = TWO * PI * c_au / lambda_au
      end function lambda_nm_to_omega


      pure function time_tofs(t_au)
            !
            ! Convert time in atomic units to femtoseconds.
            ! Source: CODATA 2010
            !
            real(F64) :: time_tofs
            real(F64), intent(in) :: t_au

            time_tofs = t_au * 2.418884326502E-2_F64
      end function time_tofs


      elemental function toang(l_bohr)
            !
            ! Convert length in atomic units (bohr) to length in ansgroms.
            !
            real(F64) :: toang
            real(F64), intent(in) :: l_bohr

            toang = l_bohr / tobohr(one)
      end function toang


      function tokcal(e_hartree)
            !
            ! Hartree -> kcal / mol conversion using
            ! CODATA 2006 values of the Hartree unit, 
            ! Avogadro number, and the definition
            ! of calorie, cal = 4.184 J
            !
            real(F64)             :: tokcal
            real(F64), intent(in) :: e_hartree

            tokcal = 627.5094688043E+0_F64 * e_hartree
      end function tokcal


      function toev(e_hartree)
            !
            ! Hartree -> electron volt conversion.
            ! Conversion factor is taken directly
            ! from CODATA 2006 tables.
            !
            real(F64)             :: toev
            real(F64), intent(in) :: e_hartree

            toev = 27.21138386E+0_F64 * e_hartree
      end function toev


      function tocm_1(e_hartree)
            !
            ! Hartree -> reciprocal centimeters.
            ! Conversion factor is taken directly
            ! from CODATA 2010 tables.                                 
            !                                
            real(F64)             :: tocm_1
            real(F64), intent(in) :: e_hartree

            tocm_1 = 2.194746313708E+5_F64 * e_hartree
      end function tocm_1


      function tom_1(e_hartree)
            !
            ! Hartree -> reciprocal meters. 
            ! Conversion factor is taken directly
            ! from CODATA 2010 tables.
            !
            real(F64)             :: tom_1
            real(F64), intent(in) :: e_hartree

            tom_1 = 2.194746313708E+7_F64 * e_hartree
      end function tom_1


      function tojoule(e_hartree)
            !
            ! Convert energy in Hartree
            ! to Joules. The conversion factor
            ! is obtained from CODATA 2010.
            !
            real(F64) :: tojoule
            real(F64), intent(in) :: e_hartree
            
            tojoule = 4.35974434E-18_F64 * e_hartree
      end function tojoule


      function tolambdam(e_hartree)
            !
            ! Convert energy in Hartree
            ! to wavelength in meters
            ! 
            real(F64) :: tolambdam
            real(F64), intent(in) :: e_hartree

            tolambdam = planck * light_speed / tojoule(e_hartree)            
      end function tolambdam


      function tolambdaangs(e_hartree)
            !
            ! Convert energy in Hartree
            ! to wavelength in meters                                
            !                                    
            real(F64) :: tolambdaangs
            real(F64), intent(in) :: e_hartree
            
            tolambdaangs = tolambdam(e_hartree) * 1.0E+10_F64
      end function tolambdaangs


      function todebye(d_au)
            real(F64)             :: todebye
            real(F64), intent(in) :: d_au
            !-----------------------------------------
            ! Atomic units of dipole moment -> debye
            ! Conversion factor calculated using
            ! the definition of debye unit,
            ! 1 D = 1/c * 10^{-21} C m 
            ! and CODATA 2006 values of the speed of 
            ! light and atomic unit of dipole moment,
            ! 1 unit charge * bohr radius = 
            ! 8.478 352 81 * 10^{-30} C m
            !
            todebye = 2.541746229_F64 * d_au
      end function todebye


      function totransprob_dipole(S, e_hartree, g)
            real(F64) :: totransprob_dipole
            real(F64), intent(in) :: S
            real(F64), intent(in) :: e_hartree
            integer,   intent(in) :: g

            totransprob_dipole = (sixteen * pi**three * stom2C2(S)) &
                  / (three * planck * eps_0 * real(g, F64) * tolambdam(e_hartree)**3)
      end function totransprob_dipole
      
      
      function totransprob_dipole_odw(e_hartree, S, g)
            real(F64) :: totransprob_dipole_odw
            real(F64), intent(in) :: S
            real(F64), intent(in) :: e_hartree
            integer,   intent(in) :: g

            totransprob_dipole_odw = (sixteen * pi**three * stom2C2(S)) &
                  / (three * planck * eps_0 * real(g, F64) * tolambdam(e_hartree)**3)
      end function totransprob_dipole_odw


      function totransprob_dipole2(S, e_hartree, g)
            real(F64) :: totransprob_dipole2
            real(F64), intent(in) :: S
            real(F64), intent(in) :: e_hartree
            integer,   intent(in) :: g

            totransprob_dipole2 = 303.8d+16 * S / (tolambdaangs(e_hartree)**3&
                  *g*1.499_F64)
      end function totransprob_dipole2

      
      function totransprob_dipole4(S, e_hartree, g)
            real(F64) :: totransprob_dipole4
            real(F64), intent(in) :: S
            real(F64), intent(in) :: e_hartree
            integer,   intent(in) :: g

            totransprob_dipole4 = 2.0261E+18_F64 * S / (tolambdaangs(e_hartree)**3*g)        
      end function totransprob_dipole4

      
      function totransprob_quad(S, e_hartree, g)
            real(F64) :: totransprob_quad
            real(F64), intent(in) :: S
            real(F64), intent(in) :: e_hartree
            integer,   intent(in) :: g

            ! totransprob_quad = (sixteen * pi**five * stom4C2(S)*tom_1(e_hartree)**5) &
            !       / (fifteen * planck * eps_0 * dble(g))

            totransprob_quad = (sixteen * pi**five * stom4C2(S)) &
                  / (fifteen * planck * eps_0 * dble(g) * tolambdam(e_hartree)**5)
      end function totransprob_quad


      function totransprob_quad_shortley(S, e_hartree, g)
            !
            ! Shortley, G., The Computation of Quadrupole and Magnetic-Dipole
            ! Transition probabilities, Phys. Rev., 225, 57 (1940) eq(1)
            ! Converterd from CGS metric system to SI metric system
            ! that is: multiplied by 1/(4pie_0)
            ! The final expression is not divided by (2J+1) due to the fact
            ! that the line strength computed in this program is equal to 
            ! one term from the eq(1) summation. No summation is executed, 
            ! i.e. no division is exequted.
            !
            real(F64) :: totransprob_quad_shortley
            real(F64), intent(in) :: S
            real(F64), intent(in) :: e_hartree
            integer,   intent(in) :: g

            totransprob_quad_shortley = (sixteen *two * pi**six * stom4C2(S)*tom_1(e_hartree)**5) &
                  / (four * pi * five * planck * eps_0)

      end function totransprob_quad_shortley

      function totransprob_quad_shortley2(e_hartree, S, g)
            real(F64) :: totransprob_quad_shortley2
            real(F64), intent(in) :: S
            real(F64), intent(in) :: e_hartree
            integer,   intent(in) :: g

            totransprob_quad_shortley2 = (sixteen *two * pi**six * stom4C2(S)*tom_1(e_hartree)**5) &
                  / (four * pi * five * planck * eps_0)
      end function totransprob_quad_shortley2



      function stom2C2(S)
            real(F64) :: stom2C2
            real(F64), intent(in) :: S
            real(F64) :: bohr

            bohr = four * pi * eps_0 * planckbar**two &
                  / (me * e_charge**two)

            stom2C2 = S * bohr**two * e_charge**2
      end function stom2C2


      function stom4C2(S)
            real(F64) :: stom4C2
            real(F64), intent(in) :: S
            real(F64) :: bohr

            bohr = four * pi * eps_0 * planckbar**two &
                  / (me * e_charge**two)

            stom4C2= S * bohr**four* e_charge**2
      end function stom4C2


      function stom2C(S)
            real(F64) :: stom2C
            real(F64), intent(in) :: S
            real(F64) :: bohr

            bohr = four * pi * eps_0 * planckbar**two &
                  / (me * e_charge**two)

            print*, 'bohr', bohr
            print*, 'e', e_charge

            stom2C = S * bohr**two * e_charge * 1.0E+40_F64
      end function stom2C
end module math_constants
