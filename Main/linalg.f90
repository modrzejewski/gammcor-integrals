module linalg
      use arithmetic
      use display
      use string

      implicit none

contains

      subroutine symmetric_eigenproblem(w, a, n, compute_eigenvecs)
            real(F64), dimension(:), contiguous, intent(out)      :: w
            real(F64), dimension(:, :), contiguous, intent(inout) :: a
            integer, intent(in)                                   :: n
            logical, intent(in)                                   :: compute_eigenvecs
            
            integer :: lwork, liwork, info, lda
            real(F64), dimension(:), allocatable :: work
            integer, dimension(:), allocatable :: iwork
            real(F64), dimension(1) :: work0
            integer, dimension(1) :: iwork0
            character(1) :: jobz
            external :: dsyevd

            lda = size(a, dim=1)
            if (compute_eigenvecs) then
                  jobz = "V"
            else
                  jobz = "N"
            end if
            !
            ! Compute the optimal size of temporary storage
            !
            call dsyevd(jobz, "L", n, a, lda, w, work0, -1, iwork0, -1, info)
            lwork = ceiling(work0(1))
            liwork = iwork0(1)
            allocate(work(lwork))
            allocate(iwork(liwork))
            call dsyevd(jobz, "L", n, a, lda, w, work, lwork, iwork, liwork, info)
            if (info .ne. 0) then
                  call msg("Eigensolver for symmetric matrices returned error code (" // str(info) // ")")
                  error stop
            end if
            deallocate(work)
            deallocate(iwork)
      end subroutine symmetric_eigenproblem
      

      subroutine linalg_smfill(a)
            !
            ! Reconstruct the upper triangle of a symmetric matrix
            !
            real(F64), dimension(:, :), intent(inout) :: a

            integer :: i, j, n

            n = size(a, dim=1)
            do j = 1, n
                  do i = j + 1, n
                        a(j, i) = a(i, j)
                  end do
            end do
      end subroutine linalg_smfill
      

      subroutine linalg_av_x(w, a, lda, v, m, n, alpha, beta)
            !
            ! Perform matrix-vector multiplication w = alpha * Av + beta * w
            !
            real(F64), dimension(*), intent(inout)   :: w
            real(F64), dimension(lda, *), intent(in) :: a
            integer, intent(in)                      :: lda
            real(F64), dimension(*), intent(in)      :: v
            integer, intent(in)                      :: m
            integer, intent(in)                      :: n
            real(F64), intent(in)                    :: alpha
            real(F64), intent(in)                    :: beta

            external :: dgemv
            
            call dgemv("N", m, n, alpha, a, lda, v, 1, beta, w, 1)
      end subroutine linalg_av_x
      

      subroutine linalg_vwT_x(a, lda, v, w, m, n, alpha)
            !
            ! A <- alpha * v*w**T + A
            !
            real(F64), dimension(lda, *), intent(inout) :: a
            integer, intent(in)                         :: lda
            real(F64), dimension(*), intent(in)         :: v
            real(F64), dimension(*), intent(in)         :: w
            integer, intent(in)                         :: m
            integer, intent(in)                         :: n
            real(F64), intent(in)                       :: alpha

            external :: dger

            call dger(m, n, alpha, v, 1, w, 1, a, lda)
      end subroutine linalg_vwT_x


      subroutine linalg_ab(c, a, b)
            !
            ! Compute C(1:u, 1:v) <- A(1:k, 1:l) B(1:m, 1:n)
            !
            real(F64), dimension(:, :), intent(out) :: c
            real(F64), dimension(:, :), intent(in)  :: a
            real(F64), dimension(:, :), intent(in)  :: b
            
            integer :: m, n, k, l, u, v
            real(F64), parameter :: alpha = ONE
            real(F64), parameter :: beta = ZERO
            external :: dgemm

            c = ZERO
            k = size(a, dim=1)
            l = size(a, dim=2)
            m = size(b, dim=1)
            n = size(b, dim=2)
            u = size(c, dim=1)
            v = size(c, dim=2)
            if (l .ne. m) then
                  call msg("real_ab: inconsistent dimensions of matrices A and B")
                  stop
            end if
            if ((u .ne. k) .or. (v .ne. n)) then
                  call msg("real_ab: inconsistent dimensions of matrix C")
                  stop
            end if
            call dgemm("N", "N", u, v, l, alpha, a, k, b, m, beta, c, u)
      end subroutine linalg_ab


      subroutine linalg_abT(c, a, b, alpha, beta)
            !
            ! Compute C(1:u, 1:v) <- A(1:k, 1:l) B(1:m, 1:n)**T
            !
            real(F64), dimension(:, :), intent(out) :: c
            real(F64), dimension(:, :), intent(in)  :: a
            real(F64), dimension(:, :), intent(in)  :: b
            real(F64), optional, intent(in)         :: alpha
            real(F64), optional, intent(in)         :: beta
            
            integer :: m, n, k, l, u, v
            real(F64) :: alpha0, beta0
            external :: dgemm

            c = ZERO
            k = size(a, dim=1)
            l = size(a, dim=2)
            m = size(b, dim=1)
            n = size(b, dim=2)
            u = size(c, dim=1)
            v = size(c, dim=2)
            if (n .ne. l) then
                  call msg("real_abT: inconsistent dimensions of matrices A and B")
                  error stop
            end if
            if ((u .ne. k) .or. (v .ne. m)) then
                  call msg("real_abT: inconsistent dimensions of matrix C")
                  error stop
            end if
            if (present(alpha)) then
                  alpha0 = alpha
            else
                  alpha0 = ONE
            end if
            if (present(beta)) then
                  beta0 = beta
            else
                  beta0 = ZERO
            end if
            call dgemm("N", "C", u, v, l, alpha0, a, k, b, m, beta0, c, u)
      end subroutine linalg_abT
      

      subroutine linalg_aTb(c, a, b)
            !
            ! Compute C(1:u, 1:v) <- A(1:k, 1:l)**T B(1:m, 1:n)
            !
            real(F64), dimension(:, :), intent(out) :: c
            real(F64), dimension(:, :), intent(in)  :: a
            real(F64), dimension(:, :), intent(in)  :: b
            
            integer :: m, n, k, l, u, v
            real(F64), parameter :: alpha = ONE
            real(F64), parameter :: beta = ZERO
            external :: dgemm

            c = ZERO
            k = size(a, dim=1)
            l = size(a, dim=2)
            m = size(b, dim=1)
            n = size(b, dim=2)
            u = size(c, dim=1)
            v = size(c, dim=2)
            if (k .ne. m) then
                  call msg("real_aTb: inconsistent dimensions of matrices A and B")
                  stop
            end if
            if ((u .ne. l) .or. (v .ne. n)) then
                  call msg("real_aTb: inconsistent dimensions of matrix C")
                  stop
            end if
            call dgemm("C", "N", u, v, k, alpha, a, k, b, m, beta, c, u)
      end subroutine linalg_aTb

      
      subroutine linalg_abT_x(c, ldc, a, lda, b, ldb, m, n, k, alpha, beta)
            !
            ! Compute C(1:m, 1:n) <- alpha*A(1:m, 1:k) B(1:n, 1:k)**T + beta*C(1:m, 1:n)
            !
            real(F64), dimension(ldc, *), intent(out) :: c
            integer, intent(in)                       :: ldc
            real(F64), dimension(lda, *), intent(in)  :: a
            integer, intent(in)                       :: lda
            real(F64), dimension(ldb, *), intent(in)  :: b
            integer, intent(in)                       :: ldb
            integer, intent(in)                       :: m
            integer, intent(in)                       :: n
            integer, intent(in)                       :: k
            real(F64), optional, intent(in)           :: alpha
            real(F64), optional, intent(in)           :: beta
            
            real(F64) :: alpha0, beta0
            external :: dgemm

            if (present(alpha)) then
                  alpha0 = alpha
            else
                  alpha0 = 1.0_F64
            end if

            if (present(beta)) then
                  beta0 = beta
            else
                  beta0 = 0.0_F64
            end if

            call dgemm("N", "C", m, n, k, alpha0, a, lda, b, ldb, beta0, c, ldc)
      end subroutine linalg_abT_x


      subroutine linalg_aTb_x(c, ldc, a, lda, b, ldb, m, n, k, alpha, beta)
            !
            ! Compute C(1:m, 1:n) <- alpha*A(1:k, 1:m)**T B(1:k, 1:n) + beta*C(1:m, 1:n)
            !
            real(F64), dimension(ldc, *), intent(out) :: c
            integer, intent(in)                       :: ldc
            real(F64), dimension(lda, *), intent(in)  :: a
            integer, intent(in)                       :: lda
            real(F64), dimension(ldb, *), intent(in)  :: b
            integer, intent(in)                       :: ldb
            integer, intent(in)                       :: m
            integer, intent(in)                       :: n
            integer, intent(in)                       :: k
            real(F64), optional, intent(in)           :: alpha, beta
            
            real(F64), parameter :: alpha_default = 1.0_F64
            real(F64), parameter :: beta_default = 0.0_F64
            real(F64) :: alpha_param, beta_param
            external :: dgemm

            if (present(alpha)) then
                  alpha_param = alpha
            else
                  alpha_param = alpha_default
            end if
            if (present(beta)) then
                  beta_param = beta
            else
                  beta_param = beta_default
            end if
            call dgemm("C", "N", m, n, k, alpha_param, a, lda, b, ldb, beta_param, c, ldc)
      end subroutine linalg_aTb_x


      subroutine linalg_aTv_x(w, a, lda, v, m, n, alpha, beta)
            !
            ! Perform matrix-vector multiplication w = alpha * A**T v + beta * w
            !
            real(F64), dimension(*), intent(inout)   :: w
            real(F64), dimension(lda, *), intent(in) :: a
            integer, intent(in)                      :: lda
            real(F64), dimension(*), intent(in)      :: v
            integer, intent(in)                      :: m
            integer, intent(in)                      :: n
            real(F64), intent(in)                    :: alpha
            real(F64), intent(in)                    :: beta

            external :: dgemv
            
            call dgemv("T", m, n, alpha, a, lda, v, 1, beta, w, 1)
      end subroutine linalg_aTv_x
end module linalg
