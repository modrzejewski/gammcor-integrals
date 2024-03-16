module boys
      use math_constants
      use arithmetic
      use chebinterp
      use display
      use clock
      ! ---------------------------------------------------------
      !                      THE BOYS FUNCTION
      ! ---------------------------------------------------------
      ! This module contains a set of subroutines for efficient
      ! computation of the Boys function
      !
      ! Fm(x) = Int(0,1) t**(2m) exp(-x t**2) dt
      !
      ! The accuracy of Fm(x) is controlled by
      ! the BOYS_ACCURACY parameter.
      !
      implicit none

      real(F64), dimension(4), parameter :: X_ASYMPTOTIC_CHEB5 = [21.0_F64, 23.0_F64, 25.0_F64, 27.0_F64]
      integer, dimension(4), parameter :: NINTERVALS_CHEB5 = [170, 300, 500, 850]
      real(F64), dimension(4), parameter :: X_ASYMPTOTIC_CHEB7 = [21.0_F64, 23.0_F64, 25.0_F64, 27.0_F64]
      integer, dimension(4), parameter :: NINTERVALS_CHEB7 = [40, 60, 85, 130]
      ! ------------------------------------------------------
      ! Levels of accuracy for 5-point Chebyshev interpolation
      ! ------------------------------------------------------
      ! BOYS_ACCURACY   MaxUE        TIME
      ! 1               4.044E-10    5.749E+00
      ! 2               3.668E-11    5.722E+00
      ! 3               4.306E-12    5.706E+00
      ! 4               4.414E-13    5.706E+00
      ! ------------------------------------------------------
      ! Levels of accuracy for 7-point Chebyshev interpolation
      ! ------------------------------------------------------
      ! 1                3.656E-10   6.840E+00
      ! 2                2.923E-11   6.848E+00
      ! 3                4.005E-12   6.849E+00
      ! 4                3.395E-13   7.026E+00
      ! ------------------------------------------------------
      ! MaxUE is the maximum unsigned error of F_m(X), with m <= 24;
      ! verified on the interval (0, 60).
      !
      integer, parameter :: BOYS_ACCURACY = 3
      real(F64), parameter, private :: X_ASYMPTOTIC = X_ASYMPTOTIC_CHEB5(BOYS_ACCURACY)
      integer, parameter, private :: NINTERVALS = NINTERVALS_CHEB5(BOYS_ACCURACY)
      !
      ! Number of points used in the Chebyshev interpolation.
      ! Remember to change both INTERP_N and the function FM_INTERP.
      !
      integer, parameter, private :: INTERP_N = 5
      real(F64), parameter, private :: DELTAX = X_ASYMPTOTIC / real(NINTERVALS, F64)
      real(F64), dimension(:, :, :), allocatable, private, save :: INTERPTABLE
      integer, private, save :: MAXDEG

contains

      subroutine boys_init(m)
            !
            ! Initialize this module for the maximum order of the Boys function
            ! equal to M
            !
            integer, intent(in) :: m

            integer :: i, j, k
            real(F64) :: a, b
            real(F64), dimension(INTERP_N) :: interp_x, interp_f

            maxdeg = m

            allocate(INTERPTABLE(INTERP_N, NINTERVALS, 0:MAXDEG))
            !
            ! Generate coefficients of the expansion of F_m(X) in terms
            ! of the Chebyshev polynomials
            !
            do k = 0, MAXDEG
                  do i = 1, NINTERVALS
                        a = real(i - 1, F64) * DELTAX
                        b = real(i, F64) * DELTAX
                        call cheb_nodes(interp_x, INTERP_N, a, b)
                        do j = 1, INTERP_N
                              interp_f(j) = fm_reference(k, interp_x(j))
                        end do
                        call cheb_coeffs(INTERPTABLE(:, i, k), interp_f, INTERP_N)
                  end do
            end do
      end subroutine boys_init
      
      
      subroutine boys_free()
            deallocate(INTERPTABLE)
      end subroutine boys_free


      pure function fm_reference(m, x)
            ! 
            ! Compute the reference F_m(X). The unsigned relative error
            ! is smaller 4.0E-15. Use this function to compute interpolation nodes.
            !
            real(F64)             :: fm_reference
            integer, intent(in)   :: m
            real(F64), intent(in) :: x

            real(F64) :: expmx, twox
            real(F64) :: f, beta
            integer :: l, a

            integer :: deltam
            !
            ! All the magical values below were manually calibrated
            ! to achieve the target unsigned relative error < 4.0E-15
            ! for all X >= 0.
            !
            if (x < 35.0_F64) then
                  !
                  ! Compute approximate F_{m+d}(X) for a sufficiently high d,
                  ! and refine this approximation by the downward recurrence.
                  !
                  if (x < 1.0_F64) then
                        deltam = 15
                  else if (x < 5.0_F64) then
                        deltam = 30
                  else if (x < 10.0_F64) then
                        deltam = 50
                  else if (x < 15.0_F64) then
                        deltam = 60
                  else if (x < 20.0_F64) then
                        deltam = 70
                  else if (x < 25.0_F64) then
                        deltam = 80
                  else if (x < 30.0_F64) then
                        deltam = 90
                  else
                        deltam = 100
                  end if

                  expmx = exp(-x)
                  twox = TWO * x
                  f = fm_shavitt_expansion(m+deltam, twox, expmx)
                  a = 2 * (m + deltam - 1) + 1

                  do l = 1, deltam
                        f = ONE / real(a, F64) * (twox * f + expmx)
                        a = a - 2
                  end do
            else
                  !
                  ! Large x. Use the asymptotic expression for F_0(X)
                  !
                  f = (ONE/TWO)  * sqrt(pi / x)
                  expmx = exp(-x)
                  a = 1
                  beta = ONE / (TWO * x)
                  !
                  ! Upward recursion
                  !
                  do l = 2, m + 1
                        f = beta * (real(a, F64) * f - expmx)
                        a = a + 2
                  end do
            end if

            fm_reference = f
      end function fm_reference


      pure function fm_shavitt_expansion(m, twox, expmx)
            !
            ! 1. Eq. 6 in J. Math. Chem. 36, 301 (2004);
            !    doi: 10.1023/B:JOMC.0000044226.49921.f5
            ! 2. Shavitt, I. Methods in Computational Physics, eds.
            !    B. Alder, S. Fernbach, and M. Rotenberg (Academic
            !    Press, New York, 1963)
            !
            real(F64) :: fm_shavitt_expansion

            integer, intent(in)   :: m
            real(F64), intent(in) :: twox
            real(F64), intent(in) :: expmx

            real(F64) :: xi, sum
            integer :: a, d
            integer :: i

            integer, parameter :: nterms = 6

            sum = ZERO
            a = 2 * m + 1
            d = a
            xi = ONE

            do i = 0, nterms - 1
                  sum = sum + xi / real(d, F64)
                  xi = xi * twox
                  a = a + 2
                  d = d * a
            end do

            fm_shavitt_expansion = expmx * sum
      end function fm_shavitt_expansion


      pure function fm_interp(x, m)
            !
            ! Return an interpolated value of the Boys function F_m(X)
            ! for a given X in (0, X_ASYMPTOTIC).
            !
            real(F64)             :: fm_interp
            real(F64), intent(in) :: x
            integer, intent(in)   :: m

            real(F64) :: ab_sum
            integer :: k

            k = floor(x / DELTAX)
            !
            ! a = real(k, F64) * DELTAX
            ! b = real(k+1, F64) * DELTAX
            !
            ab_sum = real(2 * k + 1, F64) * DELTAX
            call cheb_approx_5(fm_interp, x, INTERPTABLE(:, k+1, m), ab_sum, DELTAX)
      end function fm_interp
      

      pure subroutine fm(m, x, fmarray)
            ! -----------------------------------------------------------------
            ! Evaluate the Boys function F_n(X) for n=0, ..., m:
            ! F_n(x) = \int_0^1 t^{2n} \exp(-x t^2) dt
            ! -----------------------------------------------------------------
            ! 1. T. Helgaker, Molecular Electronic-Structure Theory,
            !    Eqs. 9.8.24, 9.8.13
            ! 2. Guseinov, I.I., Mamedov, B.A., J. Math. Chem., 40, 179 (2006)
            ! 3. Mamedov, B.A., J. Math. Chem. 36, 301 (2004)
            ! 4. Guseinov, I.I., Mamedov, B.A., J. Math. Chem., 36, 341 (2004)
            ! -----------------------------------------------------------------
            ! M       - Maximum order of F_m(x). The Boys function F_m(X) for
            !           m = 0, ..., m is stored in the FMARRAY array
            ! X       - Point at which function is evaluated
            ! FMARRAY - Output array: F_0(x) is written to FMARRAY(1),
            !           F_1(x) to FMARRAY(2), etc.
            !
            integer, intent(in) :: m
            real(F64), intent(in) :: x
            real(F64), dimension(:), intent(out) :: fmarray

            real(F64) :: beta, expmx, a, twox
            real(F64) :: f
            integer :: l

            if (x < X_ASYMPTOTIC) then
                  f = fm_interp(x, m)
                  fmarray(m + 1) = f
                  expmx = exp(-x)
                  twox = TWO * x
                  a = real(2 * (m - 1) + 1, F64)
                  !
                  ! Downward recursion
                  !
                  do l = m, 1, -1
                        f = (twox * f + expmx) / a
                        fmarray(l) = f
                        a = a - TWO
                  end do
            else
                  !
                  ! Large x. Use the asymptotic expression for F_0(X)
                  !
                  f = (ONE/TWO)  * sqrt(pi / x)
                  fmarray(1) = f
                  expmx = exp(-x)
                  a = ONE
                  beta = ONE / (TWO * x)
                  !
                  ! Upward recursion
                  !
                  do l = 2, m + 1
                        f = beta * (a * f - expmx)
                        a = a + TWO
                        fmarray(l) = f
                  end do
            end if
      end subroutine fm


      pure subroutine f0(x, fmarray)
            real(F64), intent(in) :: x
            real(F64), dimension(:), intent(out) :: fmarray

            if (x < X_ASYMPTOTIC) then
                  fmarray(1) = fm_interp(x, 0)
            else
                  fmarray(1) = (ONE/TWO)  * sqrt(pi / x)
            end if
      end subroutine f0

      
      pure subroutine f1(x, fmarray)
            real(F64), intent(in) :: x
            real(F64), dimension(:), intent(out) :: fmarray

            real(F64) :: f, beta, expmx, twox

            if (x < X_ASYMPTOTIC) then
                  f = fm_interp(x, 1)
                  fmarray(2) = f

                  expmx = exp(-x)
                  twox = TWO * x
                  f = twox * f + expmx
                  fmarray(1) = f
            else
                  f = (ONE/TWO)  * sqrt(pi / x)
                  fmarray(1) = f
                  expmx = exp(-x)

                  beta = ONE / (TWO * x)
                  f = beta * (f - expmx)
                  fmarray(2) = f
            end if
      end subroutine f1


      pure subroutine f2(x, fmarray)
            real(F64), intent(in) :: x
            real(F64), dimension(:), intent(out) :: fmarray

            real(F64) :: f, beta, expmx, twox

            if (x < X_ASYMPTOTIC) then
                  f = fm_interp(x, 2)
                  fmarray(3) = f


                  expmx = exp(-x)
                  twox = TWO * x

                  f = (twox * f + expmx) / 3.0_F64
                  fmarray(2) = f
                  f = twox * f + expmx
                  fmarray(1) = f
            else
                  f = (ONE/TWO)  * sqrt(pi / x)
                  fmarray(1) = f
                  expmx = exp(-x)
                  
                  beta = ONE / (TWO * x)
                  f = beta * (f - expmx)
                  fmarray(2) = f
                  f = beta * (3.0_F64 * f - expmx)
                  fmarray(3) = f
            end if
      end subroutine f2


      pure subroutine f3(x, fmarray)
            real(F64), intent(in) :: x
            real(F64), dimension(:), intent(out) :: fmarray

            real(F64) :: beta, expmx, twox
            real(F64) :: f

            if (x < X_ASYMPTOTIC) then
                  f = fm_interp(x, 3)
                  fmarray(4) = f
                  expmx = exp(-x)
                  twox = TWO * x

                  f = (twox * f + expmx) / 5.0_F64
                  fmarray(3) = f
                  f = (twox * f + expmx) / 3.0_F64
                  fmarray(2) = f
                  f = twox * f + expmx
                  fmarray(1) = f
            else
                  f = (ONE/TWO)  * sqrt(pi / x)
                  fmarray(1) = f
                  expmx = exp(-x)

                  beta = ONE / (TWO * x)
                  f = beta * (f - expmx)
                  fmarray(2) = f
                  f = beta * (3.0_F64 * f - expmx)
                  fmarray(3) = f
                  f = beta * (5.0_F64 * f - expmx)
                  fmarray(4) = f
            end if
      end subroutine f3


      pure subroutine f4(x, fmarray)
            real(F64), intent(in) :: x
            real(F64), dimension(:), intent(out) :: fmarray

            real(F64) :: beta, expmx, twox
            real(F64) :: f

            if (x < X_ASYMPTOTIC) then
                  f = fm_interp(x, 4)
                  fmarray(5) = f
                  expmx = exp(-x)
                  twox = TWO * x

                  f = (twox * f + expmx) / 7.0_F64
                  fmarray(4) = f
                  f = (twox * f + expmx) / 5.0_F64
                  fmarray(3) = f
                  f = (twox * f + expmx) / 3.0_F64
                  fmarray(2) = f
                  f = twox * f + expmx
                  fmarray(1) = f
            else
                  f = (ONE/TWO)  * sqrt(pi / x)
                  fmarray(1) = f
                  expmx = exp(-x)
                  beta = ONE / (TWO * x)

                  f = beta * (f - expmx)
                  fmarray(2) = f
                  f = beta * (3.0_F64 * f - expmx)
                  fmarray(3) = f
                  f = beta * (5.0_F64 * f - expmx)
                  fmarray(4) = f
                  f = beta * (7.0_F64 * f - expmx)
                  fmarray(5) = f
            end if
      end subroutine f4


      pure subroutine f5(x, fmarray)
            real(F64), intent(in) :: x
            real(F64), dimension(:), intent(out) :: fmarray

            real(F64) :: beta, expmx, twox
            real(F64) :: f

            if (x < X_ASYMPTOTIC) then
                  f = fm_interp(x, 5)
                  fmarray(6) = f
                  expmx = exp(-x)
                  twox = TWO * x

                  f = (twox * f + expmx) / 9.0_F64
                  fmarray(5) = f
                  f = (twox * f + expmx) / 7.0_F64
                  fmarray(4) = f
                  f = (twox * f + expmx) / 5.0_F64
                  fmarray(3) = f
                  f = (twox * f + expmx) / 3.0_F64
                  fmarray(2) = f
                  f = twox * f + expmx
                  fmarray(1) = f
            else
                  f = (ONE/TWO)  * sqrt(pi / x)
                  fmarray(1) = f
                  expmx = exp(-x)

                  beta = ONE / (TWO * x)
                  f = beta * (f - expmx)
                  fmarray(2) = f
                  f = beta * (3.0_F64 * f - expmx)
                  fmarray(3) = f
                  f = beta * (5.0_F64 * f - expmx)
                  fmarray(4) = f
                  f = beta * (7.0_F64 * f - expmx)
                  fmarray(5) = f
                  f = beta * (9.0_F64 * f - expmx)
                  fmarray(6) = f
            end if
      end subroutine f5


      pure subroutine f6(x, fmarray)
            real(F64), intent(in) :: x
            real(F64), dimension(:), intent(out) :: fmarray

            real(F64) :: beta, expmx, twox
            real(F64) :: f

            if (x < X_ASYMPTOTIC) then
                  f = fm_interp(x, 6)
                  fmarray(7) = f
                  expmx = exp(-x)
                  twox = TWO * x

                  f = (twox * f + expmx) / 11.0_F64
                  fmarray(6) = f
                  f = (twox * f + expmx) / 9.0_F64
                  fmarray(5) = f
                  f = (twox * f + expmx) / 7.0_F64
                  fmarray(4) = f
                  f = (twox * f + expmx) / 5.0_F64
                  fmarray(3) = f
                  f = (twox * f + expmx) / 3.0_F64
                  fmarray(2) = f
                  f = twox * f + expmx
                  fmarray(1) = f
            else
                  f = (ONE/TWO)  * sqrt(pi / x)
                  fmarray(1) = f
                  expmx = exp(-x)

                  beta = ONE / (TWO * x)
                  f = beta * (f - expmx)
                  fmarray(2) = f
                  f = beta * (3.0_F64 * f - expmx)
                  fmarray(3) = f
                  f = beta * (5.0_F64 * f - expmx)
                  fmarray(4) = f
                  f = beta * (7.0_F64 * f - expmx)
                  fmarray(5) = f
                  f = beta * (9.0_F64 * f - expmx)
                  fmarray(6) = f
                  f = beta * (11.0_F64 * f - expmx)
                  fmarray(7) = f
            end if
      end subroutine f6


      subroutine boys_unittest()
            integer, parameter :: m_max = 24
            real(F64), parameter :: dx = 0.0005_F64
            real(F64) :: x
            integer :: k, l
            real(F64) :: max_error, max_error_x
            real(F64) :: e, fm_ref
            real(F64), dimension(m_max+1) :: fmarray_m
            type(tclock) :: timer
            real(F64) :: t

            call boys_init(m_max)
            max_error = ZERO
            max_error_x = ZERO
            
            do k = 1, nint(50.0_F64 / dx)
                  x = real((k - 1), F64) * dx

                  call fm(m_max, x, fmarray_m)

                  do l = 1, m_max + 1
                        fm_ref = fm_reference(l-1, x)
                        e = abs(fmarray_m(l) - fm_ref) / fm_ref
                        if (e > max_error) then 
                              max_error = e
                              max_error_x = x
                        end if
                  end do
            end do

            do k = 1, NINTERVALS
                  x = real((k - 1), F64) * DELTAX

                  call fm(m_max, x, fmarray_m)
                  do l = 1, m_max + 1
                        fm_ref = fm_reference(l-1, x)
                        e = abs(fmarray_m(l) - fm_ref) / fm_ref
                        if (e > max_error) then 
                              max_error = e
                              max_error_x = x
                        end if
                  end do
            end do

            call clock_start(timer)

            !$omp parallel private(x, fmarray_m)
            !$omp do schedule(guided)
            do k = 1, 10**9
                  x = real(modulo(k, NINTERVALS), F64) * DELTAX
                  call fm(m_max, x, fmarray_m)
            end do 
            !$omp end do nowait
            !$omp end parallel

            t = clock_readwall(timer)

            call msg("TESTING COMPUTATION OF THE BOYS FUNCTION")
            call imsg("MAX. ORDER OF F_M(X)", m_max)
            call dmsg("MAX. UNSIGNED RELATIVE ERROR", max_error, fmt="ES10.3")
            call dmsg("MAX. ERROR OCCURED AT POINT X =", max_error_x, fmt="F10.6")
            call dmsg("PERFORMANCE MEASURE [s]", t, fmt="ES10.3")

            call boys_free()
      end subroutine boys_unittest
end module boys
