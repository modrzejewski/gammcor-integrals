module chebinterp
      use arithmetic
      use math_constants
      
      implicit none
      
contains
      
      pure subroutine cheb_nodes(x, n, a, b)
            !
            ! Compute the roots of the n-th order Chebyshev polynomial of the first kind
            ! and map them onto the interval (a, b).
            !
            real(F64), dimension(:), intent(out) :: x
            integer, intent(in)                  :: n
            real(F64), intent(in)                :: a
            real(F64), intent(in)                :: b

            integer :: k
            real(F64) :: t, w0, w1, w2
            
            w0 = PI / real(n, F64)
            w1 = (b - a) / TWO
            w2 = (b + a) / TWO

            do k = 1, n
                  t = w0 * (real(k, F64) - FRAC12)
                  x(k) = w1 * cos(t) + w2
            end do
      end subroutine cheb_nodes

      
      pure subroutine cheb_coeffs(c, f, n)
            !
            ! Compute the coefficients of the Chebyshev approximation to F(x)
            ! (sum of Chebyshev polynomials of order <= (n-1))
            !
            real(F64), dimension(:), intent(out) :: c
            real(F64), dimension(:), intent(in)  :: f
            integer, intent(in)                  :: n

            integer :: j, k
            real(F64) :: s
            real(F64) :: w0, w1, t

            w0 = PI / real(n, F64)
            do j = 1, n
                  s = ZERO
                  w1 = real(j - 1, F64)
                  do k = 1, n
                        t = w0 * w1 * (real(k, F64) - FRAC12)
                        s = s + f(k) * cos(t)
                  end do
                  c(j) = TWO / real(n, F64) * s
            end do
      end subroutine cheb_coeffs


      pure subroutine cheb_approx_n(f, x, c, n, ab_sum, ab_diff)
            !
            ! Compute the Chebyshev approximation to F(x)
            ! (sum of Chebyshev polynomials of order <= (n-1))
            ! Use Clenshaw's efficient summation. [1]
            !
            ! 1. Clenshaw, C.W., A note on the summation of
            !    Chebyshev series, Math. Comp. 9, 118 (1955);
            !    doi: 10.1090/S0025-5718-1955-0071856-0 
            !
            real(F64), intent(out)              :: f
            real(F64), intent(in)               :: x
            real(F64), dimension(:), intent(in) :: c
            integer, intent(in)                 :: n
            real(F64), intent(in)               :: ab_sum
            real(F64), intent(in)               :: ab_diff

            integer :: j
            real(F64) :: y, twoy
            real(F64) :: d0, d1, d2
            
            y = (TWO * x - ab_sum) / ab_diff
            twoy = TWO * y
            d1 = ZERO
            d2 = ZERO
            do j = n, 2, -1
                  d0 = twoy * d1 - d2 + c(j)
                  d2 = d1
                  d1 = d0
            end do
            f = y * d1 - d2 + FRAC12 * c(1)
      end subroutine cheb_approx_n


      pure subroutine cheb_approx_4(f, x, c, ab_sum, ab_diff)
            !
            ! Compute the Chebyshev approximation to F(x)
            ! (sum of Chebyshev polynomials of order <= 3)
            !
            real(F64), intent(out)              :: f
            real(F64), intent(in)               :: x
            real(F64), dimension(:), intent(in) :: c
            real(F64), intent(in)               :: ab_sum
            real(F64), intent(in)               :: ab_diff

            real(F64) :: y, twoy
            real(F64) :: d2, d3, d4
            
            y = (TWO * x - ab_sum) / ab_diff
            twoy = TWO * y
            !
            ! d(j) = 2 * y * d(j+1) - d(j+2) + c(j)
            !
            d4 = c(4)
            d3 = twoy * d4 + c(3)
            d2 = twoy * d3 - d4 + c(2)
            !
            ! F(x) = y * d(2) - d(3) + 1/2 * c(1)
            !
            f = y * d2 - d3 + FRAC12 * c(1)
      end subroutine cheb_approx_4


      pure subroutine cheb_approx_5(f, x, c, ab_sum, ab_diff)
            !
            ! Compute the Chebyshev approximation to F(x)
            ! (sum of Chebyshev polynomials of order <= 4)
            !
            real(F64), intent(out)              :: f
            real(F64), intent(in)               :: x
            real(F64), dimension(:), intent(in) :: c
            real(F64), intent(in)               :: ab_sum
            real(F64), intent(in)               :: ab_diff

            real(F64) :: y, twoy
            real(F64) :: d2, d3, d4, d5
            
            y = (TWO * x - ab_sum) / ab_diff
            twoy = TWO * y
            !
            ! d(j) = 2 * y * d(j+1) - d(j+2) + c(j)
            !
            d5 = c(5)
            d4 = twoy * d5 + c(4)
            d3 = twoy * d4 - d5 + c(3)
            d2 = twoy * d3 - d4 + c(2)
            !
            ! F(x) = y * d(2) - d(3) + 1/2 * c(1)
            !
            f = y * d2 - d3 + FRAC12 * c(1)
      end subroutine cheb_approx_5


      pure subroutine cheb_approx_6(f, x, c, ab_sum, ab_diff)
            !
            ! Compute the Chebyshev approximation to F(x)
            ! (sum of Chebyshev polynomials of order <= 5)
            !
            real(F64), intent(out)              :: f
            real(F64), intent(in)               :: x
            real(F64), dimension(:), intent(in) :: c
            real(F64), intent(in)               :: ab_sum
            real(F64), intent(in)               :: ab_diff

            real(F64) :: y, twoy
            real(F64) :: d2, d3, d4, d5, d6
            
            y = (TWO * x - ab_sum) / ab_diff
            twoy = TWO * y
            !
            ! d(j) = 2 * y * d(j+1) - d(j+2) + c(j)
            !
            d6 = c(6)
            d5 = twoy * d6 + c(5)
            d4 = twoy * d5 - d6 + c(4)
            d3 = twoy * d4 - d5 + c(3)
            d2 = twoy * d3 - d4 + c(2)
            !
            ! F(x) = y * d(2) - d(3) + 1/2 * c(1)
            !
            f = y * d2 - d3 + FRAC12 * c(1)
      end subroutine cheb_approx_6


      pure subroutine cheb_approx_7(f, x, c, ab_sum, ab_diff)
            !
            ! Compute the Chebyshev approximation to F(x)
            ! (sum of Chebyshev polynomials of order <= 6)
            !
            real(F64), intent(out)              :: f
            real(F64), intent(in)               :: x
            real(F64), dimension(:), intent(in) :: c
            real(F64), intent(in)               :: ab_sum
            real(F64), intent(in)               :: ab_diff

            real(F64) :: y, twoy
            real(F64) :: d2, d3, d4, d5, d6, d7
            
            y = (TWO * x - ab_sum) / ab_diff
            twoy = TWO * y
            !
            ! d(j) = 2 * y * d(j+1) - d(j+2) + c(j)
            !
            d7 = c(7)
            d6 = twoy * d7 + c(6)
            d5 = twoy * d6 - d7 + c(5)
            d4 = twoy * d5 - d6 + c(4)
            d3 = twoy * d4 - d5 + c(3)
            d2 = twoy * d3 - d4 + c(2)
            !
            ! F(x) = y * d(2) - d(3) + 1/2 * c(1)
            !
            f = y * d2 - d3 + FRAC12 * c(1)
      end subroutine cheb_approx_7


      pure subroutine cheb_approx_9(f, x, c, ab_sum, ab_diff)
            !
            ! Compute the Chebyshev approximation to F(x)
            ! (sum of Chebyshev polynomials of order <= 8)
            !
            real(F64), intent(out)              :: f
            real(F64), intent(in)               :: x
            real(F64), dimension(:), intent(in) :: c
            real(F64), intent(in)               :: ab_sum
            real(F64), intent(in)               :: ab_diff

            real(F64) :: y, twoy
            real(F64) :: d2, d3, d4, d5, d6, d7, d8, d9
            
            y = (TWO * x - ab_sum) / ab_diff
            twoy = TWO * y
            !
            ! d(j) = 2 * y * d(j+1) - d(j+2) + c(j)
            !
            d9 = c(9)
            d8 = twoy * d9 + c(8)
            d7 = twoy * d8 - d9 + c(7)
            d6 = twoy * d7 - d8 + c(6)
            d5 = twoy * d6 - d7 + c(5)
            d4 = twoy * d5 - d6 + c(4)
            d3 = twoy * d4 - d5 + c(3)
            d2 = twoy * d3 - d4 + c(2)
            !
            ! F(x) = y * d(2) - d(3) + 1/2 * c(1)
            !
            f = y * d2 - d3 + FRAC12 * c(1)
      end subroutine cheb_approx_9
end module chebinterp
