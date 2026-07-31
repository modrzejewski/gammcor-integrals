module io
      use arithmetic
      use display
      use string

      implicit none

      interface io_size_byte
            module procedure :: io_size_byte_rank1_F64
            module procedure :: io_size_byte_rank2_F64
            module procedure :: io_size_byte_rank3_F64
            module procedure :: io_size_byte_rank4_F64
      end interface io_size_byte
      
      !
      ! The maximum length of an I/O error message
      !
      integer, parameter :: IO_MAX_MSGLEN = 256

      !
      ! The character used by the operating system
      ! to separate pathname components.
      !
      character(len=1) :: DIRSEP = "/"

contains

      function io_exists(s)
            !
            ! Check if the file referenced by S exists.
            !
            logical :: io_exists
            character(*), intent(in) :: s
            
            inquire(file=s, exist=io_exists)
      end function io_exists


      function io_size_byte_rank1_F64(a)
            !
            ! Size of an array A in bytes. The size in bytes
            ! is stored in an integer of the I64 kind to prevent
            ! an overflow for large arrays.
            !
            integer(I64)                        :: io_size_byte_rank1_F64
            real(F64), dimension(:), intent(in) :: a

            io_size_byte_rank1_F64 = storage_size(a, kind=I64) &
                  * size(a, kind=I64) / 8_I64
      end function io_size_byte_rank1_F64


      function io_size_byte_rank2_F64(a)
            integer(I64)                           :: io_size_byte_rank2_F64
            real(F64), dimension(:, :), intent(in) :: a

            io_size_byte_rank2_F64 = storage_size(a, kind=I64) &
                  * size(a, kind=I64) / 8_I64
      end function io_size_byte_rank2_F64


      function io_size_byte_rank3_F64(a)
            integer(I64)                              :: io_size_byte_rank3_F64
            real(F64), dimension(:, :, :), intent(in) :: a

            io_size_byte_rank3_F64 = storage_size(a, kind=I64) &
                  * size(a, kind=I64) / 8_I64
      end function io_size_byte_rank3_F64


      function io_size_byte_rank4_F64(a)
            integer(I64)                                 :: io_size_byte_rank4_F64
            real(F64), dimension(:, :, :, :), intent(in) :: a

            io_size_byte_rank4_F64 = storage_size(a, kind=I64) &
                  * size(a, kind=I64) / 8_I64
      end function io_size_byte_rank4_F64


      function io_text_open(filename, s)
            !
            ! Open a text file.
            !
            integer                            :: io_text_open
            character(*), intent(in)           :: filename
            character(*), intent(in)           :: s

            integer :: open_stat
            character(len=IO_MAX_MSGLEN) :: errmsg

            open(newunit=io_text_open, file=filename, status=s, &
                  access="SEQUENTIAL", iostat=open_stat, iomsg=errmsg)
            
            if (open_stat .ne. 0) then
                  call msg("Could not open file")
                  call msg(trim(adjustl(filename)))
                  call msg(trim(errmsg))
                  error stop
            end if
      end function io_text_open


      subroutine io_text_readline(line, u, eof)
            !
            ! Read a line from a text file. The limit for the line
            ! size is MAXCHUNKS * DEFLEN characters (see the code).
            !
            character(:), allocatable, intent(out) :: line
            integer, intent(in)                    :: u
            logical, optional, intent(out)         :: eof

            character(len=80) :: chunk
            character(len=IO_MAX_MSGLEN) :: errmsg
            integer :: s, ios
            integer :: n
            integer, parameter :: maxchunks = 2**10

            line = ""
            if (present(eof)) eof = .false.
            
            lineloop: do n = 1, maxchunks
                  read(u, "(A)", advance="NO", size=s, &
                        iostat=ios, iomsg=errmsg) chunk

                  if (s > 0) then
                        line = line // chunk(1:s)
                  end if

                  if (ios == iostat_end) then
                        if (present(eof)) eof = .true.
                        exit lineloop
                  else if (ios == iostat_eor) then
                        exit lineloop
                  else if (ios .ne. 0) then
                        call msg("COULD NOT READ LINE")
                        call msg(trim(errmsg))
                        stop
                  end if
            end do lineloop
      end subroutine io_text_readline


      subroutine io_text_write(a, filename)
            !
            ! Save a matrix of floating-point numbers into a text file.
            ! An appropriate format is chosen to retain the full precision
            ! of the input. Use this subroutine to write portable, human-
            ! readable files. Use IO_BINARY_WRITE if the data need not be
            ! read on other machines.
            !
            real(F64), dimension(:, :), intent(in) :: a
            character(*), intent(in)           :: filename

            integer :: u
            integer :: m, n
            integer :: i, j
            character(:), allocatable :: fmt

            u = io_text_open(filename, "REPLACE")
            m = size(a, dim=1)
            n = size(a, dim=2)
            !
            ! ({N}(ES{F64_ES_W}.{F64_ES_D}E{F64_ES_E},:,1X))
            !
            fmt = "(" // str(n) // "(ES" // str(F64_ES_W) // "." // &
                  str(F64_ES_D) // "E" // str(F64_ES_E) // ",:,1X))"

            do i = 1, m
                  write(u, fmt) (a(i, j), j = 1, n)
            end do

            close(u)
      end subroutine io_text_write


      subroutine io_text_read(a, filename)
            !
            ! Read a matrix of floating-point numbers from a text file.
            ! See the comments for IO_TEXT_WRITE.
            !
            real(F64), dimension(:, :), intent(out) :: a
            character(*), intent(in)            :: filename

            integer :: u
            integer :: m, n
            integer :: i, j

            u = io_text_open(filename, "OLD")
            m = size(a, dim=1)
            n = size(a, dim=2)

            do i = 1, m
                  read(u, *) (a(i, j), j = 1, n)
            end do

            close(u)
      end subroutine io_text_read
end module io
