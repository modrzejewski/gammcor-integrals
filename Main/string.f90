module string
      use arithmetic

      implicit none

      character(len=*), parameter :: STR_LETTER_UPPER = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      character(len=*), parameter :: STR_LETTER_LOWER = "abcdefghijklmnopqrstuvwxyz"

      interface str
            module procedure :: str_i32
            module procedure :: str_i64
            module procedure :: str_f64
      end interface str

contains

      pure function str_i32(i)
            !
            ! Convert integer to string. The result does not
            ! contain any blanks.
            !
            character(len=:), allocatable :: str_i32
            integer(I32), intent(in) :: i

            character(len=I32_MAXWIDTH) :: t
            integer :: l

            write(t, "(I0)") i
            l = len_trim(t)
            allocate(character(l) :: str_i32)
            str_i32(:) = t(1:l)
      end function str_i32

      
      pure function str_i64(i)
            !
            ! Convert integer to string. The result does not
            ! contain any blanks.
            !
            character(:), allocatable :: str_i64
            integer(I64), intent(in) :: i

            character(len=I64_MAXWIDTH) :: t
            integer :: l

            write(t, "(I0)") i
            l = len_trim(t)
            allocate(character(l) :: str_i64)
            str_i64(:) = t(1:l)
      end function str_i64


      pure function str_f64(f, d)
            !
            ! Convert floating point number to string.
            ! The result does not contain any blanks.
            !
            character(:), allocatable     :: str_f64
            real(F64), intent(in)         :: f
            integer, optional, intent(in) :: d

            character(F64_ES_W) :: t
            character(:), allocatable :: fmt

            if (present(d)) then
                  !
                  ! (ES{F64_ES_W}.dE{F64_ES_E})
                  !
                  fmt = "(ES" // str(F64_ES_W) // "." // &
                        str(min(F64_ES_D, d)) // "E" // str(F64_ES_E) // ")"
            else
                  !
                  ! (ES{F64_ES_W}.{F64_ES_D}E{F64_ES_E})
                  !
                  fmt = "(ES" // str(F64_ES_W) // "." // &
                        str(F64_ES_D) // "E" // str(F64_ES_E) // ")"
            end if
            write(t, fmt) f
            str_f64 = trim(adjustl(t))
      end function str_f64


      pure function iscomment(s)
            logical                  :: iscomment
            character(*), intent(in) :: s
            
            character(:), allocatable :: sl

            iscomment = .false.
            if (.not. isblank(s)) then
                  sl = adjustl(s)
                  if (sl(1:1) == "!") then
                        iscomment = .true.
                  end if
            end if            
      end function iscomment


      pure function isblank(l)
            logical                      :: isblank
            character(len=*), intent(in) :: l

            if (len_trim(l) .eq. 0) then
                  isblank = .true.
            else
                  isblank = .false.
            end if
      end function isblank


      pure function isinteger(s)
            !
            ! Test if S is a character string
            ! representing a single integer number.
            !
            logical :: isinteger
            character(len=*), intent(in) :: s
            
            character(len=:), allocatable :: adjs
            integer :: k, l

            adjs = adjustl(s)
            l = len_trim(adjs)
            
            isinteger = .true.
            kloop: do k = 1, l
                  select case (adjs(k:k))
                        case ("0":"9")
                              continue
                        case ("+", "-")
                              if (k > 1) then
                                    isinteger = .false.
                                    exit kloop
                              end if
                        case default
                              isinteger = .false.
                              exit kloop
                  end select
            end do kloop
      end function isinteger

      
      pure subroutine split(s, s1, s2, delimiter)
            !
            ! Split a list of words into two pieces:
            ! "keyword    value1 value2" -> "keyword" + "value1 value2".
            !
            character(*), intent(in)               :: s
            character(:), allocatable, intent(out) :: s1
            character(:), allocatable, intent(out) :: s2
            character(1), intent(in), optional :: delimiter

            integer :: k
            character(:), allocatable :: w
            character(1) :: delim

            if (present(delimiter)) then
                  delim = delimiter
            else
                  delim = " "
            end if
            
            w = trim(adjustl(s))
            if (len(w) == 0) then
                  s1 = ""
                  s2 = ""
            else
                  k = index(w, delim)
                  if (k == 0) then
                        s1 = w
                        s2 = ""
                  else
                        s1 = w(1:k-1)
                        s2 = trim(adjustl(w(k+1:)))
                  end if
            end if
      end subroutine split


      pure function lfield(s, w)
            !
            ! Create a left-aligned character string of width W.
            ! Left-aligned cells are useful for printing table headers.
            !
            character(:), allocatable :: lfield
            character(*), intent(in) :: s
            integer, intent(in) :: w
            integer :: k

            allocate(character(w) :: lfield)
            k = min(w, len_trim(s))
            lfield(1:w) = " "
            lfield(1:k) = s(1:k)
      end function lfield


      pure function rfield(s, w)
            !
            ! Create a rightt-aligned character string of width W.
            !
            character(:), allocatable :: rfield
            character(*), intent(in) :: s
            integer, intent(in) :: w
            integer :: k, ls

            ls = len_trim(s)
            if (ls >= w) then
                  rfield = s(1:ls)
            else
                  allocate(character(w) :: rfield)
                  k = min(w, len_trim(s))
                  rfield(1:w) = " "
                  rfield(w-k+1:w) = s(1:k)
            end if
      end function rfield


      pure function cfield(s, w)
            !
            ! Create a centered character string of width W.
            ! Centered cells are useful for printing table headers.
            !
            character(:), allocatable :: cfield
            character(*), intent(in) :: s
            integer, intent(in) :: w
            integer :: k0, k1, l, ls

            ls = len_trim(s)
            if (ls >= w) then
                  cfield = s(1:ls)
            else
                  allocate(character(w) :: cfield)
                  l = min(w, len_trim(s))
                  k0 = (w - l) / 2 + 1
                  k1 = k0 + l - 1
                  cfield(1:w) = " "
                  cfield(k0:k1) = s(1:l)
            end if
      end function cfield
      

      pure function uppercase(s)
            !
            ! Convert characters to uppercase.
            ! Numbers and special characters are ignored.
            !
            character(:), allocatable :: uppercase
            character(*), intent(in)  :: s
            integer :: idx, k

            uppercase = s
            do k = 1, len_trim(s)
                  idx = index(STR_LETTER_LOWER, s(k:k))
                  if (idx > 0) then
                        uppercase(k:k) = STR_LETTER_UPPER(idx:idx)
                  end if
            end do
      end function uppercase


      pure function lowercase(s)
            !
            ! Convert characters to lowercase.
            ! Numbers and special characters are ignored.
            !
            character(:), allocatable :: lowercase
            character(*), intent(in)  :: s
            integer :: idx, k

            lowercase = s
            do k = 1, len_trim(s)
                  idx = index(STR_LETTER_UPPER, s(k:k))
                  if (idx > 0) then
                        lowercase(k:k) = STR_LETTER_LOWER(idx:idx)
                  end if
            end do
      end function lowercase


      function IntListLength(s)
            !
            ! Get the number of integers represented 
            ! in a character string S.
            !
            integer :: IntListLength
            character(*), intent(in) :: s

            integer :: l, k, i
            integer :: maxints
            integer :: t
            integer :: ios

            l = len_trim(s)
            IntListLength = 0
            if (l > 0) then
                  maxints = l / 2 + 1
                  do k = 1, maxints + 1
                        read(s, iostat=ios, fmt=*) (t, i=1, k)
                        if (ios .ne. 0) then
                              IntListLength = k - 1
                              exit
                        end if
                  end do
            end if
      end function IntListLength


      subroutine readreal(a, s)
            !
            ! Read a string of an arbitrary number of real numbers,
            ! allocate an array of the corresponding number of elements,
            ! and store the data.
            !
            real(F64), dimension(:), allocatable :: a
            character(*), intent(in)             :: s

            integer :: n, l, k, i
            integer :: maxn
            real(F64) :: t
            integer :: ios

            l = len_trim(s)
            n = 0
            if (l > 0) then
                  maxn = l / 2 + 1
                  do k = 1, maxn
                        read(s, iostat=ios, fmt=*) (t, i=1, k)
                        if (ios .ne. 0) then
                              n = k - 1
                              exit
                        end if
                  end do
            end if

            if (n > 0) then
                  allocate(a(n))
                  read(s, iostat=ios, fmt=*) (a(i), i=1, n)
            else
                  allocate(a(0))
            end if
      end subroutine readreal
end module string
