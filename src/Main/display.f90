module display
      use, intrinsic :: iso_fortran_env, only : OUTPUT_UNIT
      use arithmetic

      implicit none      
      !
      ! Standard output
      !
      integer, parameter :: STDOUNIT = OUTPUT_UNIT
      ! -------------------------------------
      ! MESSAGE PRIORITIES
      ! -------------------------------------
      !
      integer, parameter :: MSG_DEBUG  = 0
      integer, parameter :: MSG_NORMAL = 10
      integer, parameter :: MSG_WARNING = 80
      integer, parameter :: MSG_ERROR  = 100
      !
      ! Messages whose priority is strictly
      ! below MSG_PRIORITY_THRESH will not
      ! be displayed
      !
      integer, save :: MSG_PRIORITY_THRESH = MSG_NORMAL

contains

      subroutine msg(s, priority)
            !
            ! Output character string to the standard output
            !
            character(len=*), intent(in)  :: s
            integer, optional, intent(in) :: priority

            write(STDOUNIT, "(1X,A)") s
            flush(STDOUNIT)
      end subroutine msg

      
      subroutine midrule(width)
            integer, intent(in) :: width

            write(STDOUNIT, "(1X,A)") repeat("-", width)
            flush(STDOUNIT)
      end subroutine midrule


      subroutine blankline()
            write(STDOUNIT, "(1X)")
            flush(STDOUNIT)
      end subroutine blankline


      subroutine geprn(a)
            real(F64), dimension(:, :), intent(in) :: a

            integer :: k, l, m, n

            m = size(a, dim=1)
            n = size(a, dim=2)
            
            do k = 1, m
                  write(STDOUNIT, "(999F20.10)") (a(k, l), l = 1, n)
            end do
            flush(STDOUNIT)
      end subroutine geprn
end module display
