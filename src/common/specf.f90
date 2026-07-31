module specf
      use math_constants

      implicit none

      interface ipack
            module procedure ipack2
            module procedure ipack4
      end interface

      interface iunpack
            module procedure iunpack2
            module procedure iunpack4
      end interface

contains

      pure function iphase(k)
            !
            ! (-1)^K
            !
            integer :: iphase
            integer, intent(in) :: k
            
            integer :: k0

            k0 = abs(k)
            iphase = 1 - modulo(k0, 2) * 2
      end function iphase


      pure function ifactfrac(a, b)
            !
            ! Compute A!/B! (A >= B)
            !
            integer :: ifactfrac
            integer, intent(in) :: a
            integer, intent(in) :: b
            
            integer :: v, k

            v = 1
            do k = b+1, a
                  v = v * k
            end do
            ifactfrac = v
      end function ifactfrac


      pure function ibinom(n, k)
            !
            ! Compute binomial coefficient: N choose K.
            ! A recurrence formula is used,
            ! C(N, K+1) = [(N-K) * C(N, K)] / (K+1)
            !
            integer :: ibinom
            integer, intent(in) :: n
            integer, intent(in) :: k

            integer :: i
            integer :: v
            integer :: k0

            k0 = min(k, n-k)
            v = 1
            do i = 1, k0
                  v = ((n - i + 1) * v) / i
            end do
            ibinom = v
      end function ibinom


      pure function ifact(n)
            !
            ! Factorial function: N!
            !
            integer :: ifact
            integer, intent(in) :: n

            integer :: k

            ifact = 1
            do k = 1, n
                  ifact = ifact * k
            end do
      end function ifact


      pure function idblfact(n)
            !
            ! Calculate double factorial n!!
            !
            integer :: idblfact
            integer, intent(in) :: n

            integer :: v
            integer :: w, w0
            
            v = 1
            w0 = 2 - modulo(n, 2)
            do w = w0, n, 2
                  v = v * w
            end do
            idblfact = v
      end function idblfact


      pure function dblfactorial(l)
            !
            ! Calculate double factorial (2l - 1)!!
            !
            double precision :: dblfactorial

            integer, intent(in) :: l

            double precision :: s, di
            integer :: i

            s = one
            di = one
            do i = 1, l
                  s = s * di
                  di = di + two
            end do

            dblfactorial = s
      end function dblfactorial


      pure subroutine int_sincos(n, m, a, b, val)
            ! -----------------------------------
            ! Calculate definite integral
            ! \int_a^b \sin^n x \cos^m x dx
            ! -----------------------------------
            integer, intent(in)           :: n, m
            double precision, intent(in)  :: a, b
            double precision, intent(out) :: val

            double precision :: cosa, sina
            double precision :: cosb, sinb
            double precision :: c, dk, dm, t, cosam, cosbm
            integer :: k

            if (n .eq. 0) then
                  call int_cos(m, a, b, val)
            else if (m .eq. 0) then
                  call int_sin(n, a, b, val)
            else
                  cosa = cos(a)
                  sina = sin(a)
                  cosb = cos(b)
                  sinb = sin(b)
                  cosam = cosa**(m + 1)
                  cosbm = cosb**(m + 1)
                  dm = dble(m)

                  val = zero
                  c = one
                  dk = dble(n)

                  do k = n, 2, -2
                        t = sinb**(k - 1) * cosbm - sina**(k - 1) * cosam
                        val = val - c * t / (dk + dm)
                        c = c * (dk - one) / (dk + dm)
                        dk = dk - two
                  end do

                  if (modulo(n, 2) .eq. 0) then
                        call int_cos(m, a, b, t)
                        val = val + c * t
                  else
                        val = val + c * (-cosbm + cosam) / (dm + one)
                  end if
                  
            end if
      end subroutine int_sincos

      
      pure subroutine int_cos(n, a, b, val)
            ! ----------------------------------
            ! Calculate definite integral
            ! \int_a^b \cos^n x dx
            ! ----------------------------------
            integer, intent(in)           :: n
            double precision, intent(in)  :: a, b
            double precision, intent(out) :: val

            double precision :: cosa, sina
            double precision :: cosb, sinb
            double precision :: c, dk, t
            integer :: k

            if (n .eq. 0) then
                  val = b - a
            else if (n .eq. 1) then
                  sina = sin(a)
                  sinb = sin(b)
                  val = sinb - sina
            else
                  cosa = cos(a)
                  sina = sin(a)
                  cosb = cos(b)
                  sinb = sin(b)

                  val = zero
                  c = one
                  dk = dble(n)

                  do k = n, 2, -2
                        t = cosb**(k - 1) * sinb - cosa**(k - 1) * sina
                        val = val + c * t / dk 
                        c = c * (dk - one) / dk
                        dk = dk - two
                  end do

                  if (modulo(n, 2) .eq. 0) then
                        val = val + c * (b - a)
                  else
                        val = val + c * (sinb - sina)
                  end if
            end if
      end subroutine int_cos


      pure subroutine int_sin(n, a, b, val)
            ! ----------------------------------
            ! Calculate definite integral
            ! \int_a^b \sin^n x dx
            ! ----------------------------------
            integer, intent(in)           :: n
            double precision, intent(in)  :: a, b
            double precision, intent(out) :: val

            double precision :: cosa, sina
            double precision :: cosb, sinb
            double precision :: c, dk, t
            integer :: k

            if (n .eq. 0) then
                  val = b - a
            else if (n .eq. 1) then
                  cosa = cos(a)
                  cosb = cos(b)
                  val = -cosb + cosa
            else
                  cosa = cos(a)
                  sina = sin(a)
                  cosb = cos(b)
                  sinb = sin(b)

                  val = zero
                  c = one
                  dk = dble(n)

                  do k = n, 2, -2
                        t = sinb**(k - 1) * cosb - sina**(k - 1) * cosa
                        val = val - c * t / dk 
                        c = c * (dk - one) / dk
                        dk = dk - two
                  end do
                  
                  if (modulo(n, 2) .eq. 0) then
                        val = val + c * (b - a)
                  else
                        val = val + c * (-cosb + cosa)
                  end if
            end if
      end subroutine int_sin


      pure subroutine ipack2(a, b, c)
            integer, intent(in)  :: a, b
            integer, intent(out) :: c

            integer, parameter :: i0 = 0
            integer, parameter :: fieldlen = bit_size(i0) / 2
            integer, parameter :: pos1 = 0
            integer, parameter :: pos2 = fieldlen

            c = b
            call mvbits(a, pos1, fieldlen, c, pos2)
      end subroutine ipack2


      pure subroutine iunpack2(c, a, b)
            integer, intent(in)  :: c
            integer, intent(out) :: a, b

            integer, parameter :: i0 = 0
            integer, parameter :: fieldlen = bit_size(i0) / 2
            integer, parameter :: pos1 = 0
            integer, parameter :: pos2 = fieldlen

            a = ibits(c, pos2, fieldlen)
            b = ibits(c, pos1, fieldlen)
      end subroutine iunpack2


      pure subroutine ipack4(a, b, c, d, e)
            integer, intent(in)  :: a, b, c, d
            integer, intent(out) :: e

            integer, parameter :: i0 = 0
            integer, parameter :: fieldlen = bit_size(i0) / 4
            integer, parameter :: pos1 = 0
            integer, parameter :: pos2 = fieldlen
            integer, parameter :: pos3 = 2 * fieldlen
            integer, parameter :: pos4 = 3 * fieldlen

            e = d
            call mvbits(a, pos1, fieldlen, e, pos4)
            call mvbits(b, pos1, fieldlen, e, pos3)
            call mvbits(c, pos1, fieldlen, e, pos2)
      end subroutine ipack4


      pure subroutine iunpack4(e, a, b, c, d)
            integer, intent(in)  :: e
            integer, intent(out) :: a, b, c, d

            integer, parameter :: i0 = 0
            integer, parameter :: fieldlen = bit_size(i0) / 4
            integer, parameter :: pos1 = 0
            integer, parameter :: pos2 = fieldlen
            integer, parameter :: pos3 = 2 * fieldlen
            integer, parameter :: pos4 = 3 * fieldlen
            
            a = ibits(e, pos4, fieldlen)
            b = ibits(e, pos3, fieldlen)
            c = ibits(e, pos2, fieldlen)
            d = ibits(e, pos1, fieldlen)
      end subroutine iunpack4
end module specf
