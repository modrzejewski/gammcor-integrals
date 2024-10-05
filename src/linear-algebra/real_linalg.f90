! -------------------------------------------------------------
!            LINEAR ALGEBRA MODULE (REAL NUMBERS)
! -------------------------------------------------------------
!
! * SIMPLIFIED INTERFACES FOR LINEAR ALGEBRA OF REAL-VALUED
!   MATRICES
! * RECTANGULAR MATRICES
! * NONCONTIGUOUS MEMORY STORAGE
!
! The code for rectangular matrices is designed for
! two special cases of vector/matrix  operations:
! 1) One of the dimensions is orders-of-magnitude smaller
!    than the other one.
! 2) Arrays have noncontiguous elements.
!
! In subroutines where summation of possibly large number of
! floating point numbers is performed, partial sums are used
! to prevent the accumulation of numerical error.
!
module real_linalg
      use arithmetic
      use math_constants
      use display
      use string
      use clock                              

      implicit none

      interface real_vta_rect
            module procedure :: real_vta_rect_scal
            module procedure :: real_vta_rect_noscal
      end interface real_vta_rect

contains

      subroutine real_Cholesky(A)
            real(F64), dimension(:, :), intent(inout) :: A

            integer :: n, info
            integer :: i, j
            external :: dpotrf

            n = size(A, dim=1)
            call dpotrf("L", n, A, n, info)
            do j = 1, n
                  do i = 1, j - 1
                        A(i, j) = ZERO
                  end do
            end do
            if (info /= 0) then
                  call msg("Cholesky factorization failed with exit code info=" // str(info), MSG_ERROR)
                  error stop
            end if
      end subroutine real_Cholesky


      subroutine real_LowerTriangularInverse(L)
            real(F64), dimension(:, :), intent(inout) :: L

            integer :: n, info
            external :: dtrtri

            n = size(L, dim=1)
            call dtrtri("L", "N", n, L, n, info)
            if (info /= 0) then
                  call msg("Lower triangular inverse failed with exit code info=" // str(info), MSG_ERROR)
                  error stop
            end if
      end subroutine real_LowerTriangularInverse
      
            
      subroutine real_PivotedCholesky(A, P, Rank, Eps)
            !
            ! Perform pivoted Cholesky decompostion of A
            ! 
            ! Pi**T * A * Pi = U**T * U
            !
            ! where Pi is the permuation matrix
            !
            ! Pi(k,l) = KroneckerDelta(P(l),k)
            !
            ! The output array P is a compressed storage form
            ! of the permutation matrix Pi.
            !
            ! The algorithm terminates if the current diagonal
            ! element at U <= Eps.
            !
            real(F64), dimension(:, :), intent(inout) :: A
            integer, dimension(:), intent(out)        :: P
            integer, intent(out)                      :: Rank
            real(F64), intent(in)                     :: Eps

            integer :: N, info
            real(F64), dimension(:), allocatable :: Work
            external :: dpstrf
            
            N = size(A, dim=1)
            allocate(Work(2*N))
            call dpstrf("U", N, A, N, P, Rank, Eps, Work, info)
            if (info < 0) then
                  call msg("Pivoted Cholesky decomposition returned with info="//str(info), MSG_ERROR)
                  error stop
            end if
      end subroutine real_PivotedCholesky
      

      subroutine real_Lyapunov(X, B, A)
            !
            ! Solve the Lyapunov equation
            !
            ! B**T X + X B = A
            !
            ! The input arrays, A and B, will be used as working space
            ! and destroyed on exit.
            !
            real(F64), dimension(:, :), intent(out)   :: X
            real(F64), dimension(:, :), intent(inout) :: A
            real(F64), dimension(:, :), intent(inout) :: B

            integer :: n
            real(F64) :: scale
            integer :: info
            character(1), parameter :: trana = "T"
            character(1), parameter :: tranb = "N"
            real(F64), dimension(:, :), allocatable :: Z
            real(F64), dimension(:), allocatable :: wr, wi
            integer, parameter :: isgn = 1
            
            external :: dtrsyl

            n = size(A, dim=1)
            allocate(Z(n, n))
            allocate(wr(n))
            allocate(wi(n))
            call real_Schur(wr, wi, Z, B)
            !
            ! Transform A to the basis of Schur vectors
            ! A <- Z**T A Z
            ! B already has the block triangular form
            !
            call real_aTb(X, Z, A)
            call real_ab(A, X, Z)
            call dtrsyl(trana, tranb, isgn, n, n, B, n, B, n, A, n, scale, info)
            !
            ! Back transform the solution from the Schur basis
            !
            call real_abT(B, A, Z)
            call real_ab(X, Z, B)
            if (scale < ONE) then
                  call msg("Overflow in the Lyapunov equation", MSG_ERROR)
                  error stop
            end if
      end subroutine real_Lyapunov


      subroutine real_Schur(wr, wi, Z, A)
            !
            ! Compute the Schur factorization A = Z*T*(Z**T)
            ! where T is in the real Schur form.
            !
            real(F64), dimension(:), intent(out)      :: wr
            real(F64), dimension(:), intent(out)      :: wi
            real(F64), dimension(:, :), intent(out)   :: Z
            real(F64), dimension(:, :), intent(inout) :: A

            external :: dgees

            integer :: n, ldA, sdim, ldZ
            integer :: info
            integer :: lwork
            logical, dimension(1) :: bwork
            real(F64), dimension(:), allocatable :: work
            real(F64), dimension(1) :: work0
            character(1), parameter :: jobvs = "V"
            character(1), parameter :: sort = "N"
            
            n = size(A, dim=1)
            ldA = n
            ldZ = n
            !
            ! Workspace query
            !
            lwork = -1
            call dgees(jobvs, sort, select_func, n, A, ldA, sdim, &
                  wr, wi, Z, ldZ, work0, lwork, bwork, info)
            lwork = ceiling(work0(1))
            allocate(work(lwork))
            call dgees(jobvs, sort, select_func, n, A, ldA, sdim, &
                  wr, wi, Z, ldZ, work, lwork, bwork, info)
            if (info /= 0) then
                  call msg("Schur decomposition failed with info=" // str(info), MSG_ERROR)
                  error stop
            end if
            
      contains
            !
            ! Dummy function: not referenced, but required by the interface
            !
            function select_func(wr, wi)
                  logical :: select_func
                  real(F64), intent(in) :: wr
                  real(F64), intent(in) :: wi

                  select_func = .true.
            end function select_func
      end subroutine real_Schur
      

      subroutine real_SVD(U, V, Sigma, A)
            !
            ! Singular value decomposition of a real matrix A
            !
            ! A = U Diag(Sigma) V**T
            !
            ! Columns of the output matrix U/V are the left/right singular
            ! vectors of A.
            !
            ! Matrix     Dimensions
            ! A          (m, n)
            ! U          (m, m)
            ! V          (n, n)
            ! Sigma      array of size >= min(m, n)
            !
            real(F64), dimension(:, :), intent(out) :: U
            real(F64), dimension(:, :), intent(out) :: V
            real(F64), dimension(:), intent(out)    :: Sigma
            real(F64), dimension(:, :), intent(in)  :: A

            integer :: m, n
            integer :: ldA, ldVT, ldU, info
            real(F64), dimension(:, :), allocatable :: VT
            integer, dimension(:), allocatable :: iwork
            real(F64), dimension(:), allocatable :: work
            integer :: lwork

            external :: dgesdd

            m = size(A, dim=1)
            n = size(A, dim=2)
            ldA = m
            ldU = m
            if (size(U, dim=1) /= m .or. size(U, dim=2) /= m) then
                  call msg("Invalid dimensions of the left singular vectors matrix", MSG_ERROR)
                  error stop
            end if
            if (size(V, dim=1) /= n .or. size(V, dim=2) /= n) then
                  call msg("Invalid dimensions of the right singular vectors matrix", MSG_ERROR)
                  error stop
            end if
            if (size(Sigma) < min(m, n)) then
                  call msg("Invalid size of the singular values array", MSG_ERROR)
                  error stop
            end if
            allocate(VT(n, n))
            ldVT = n
            !
            ! Query the optimal scratch space
            !
            allocate(iwork(8*min(m, n)))
            allocate(work(1))
            call dgesdd("A", m, n, A, ldA, Sigma, U, ldU, VT, ldVT, work, -1, iwork, info)            
            !
            ! Proper SVD call
            !
            lwork = ceiling(work(1))
            deallocate(work)
            allocate(work(lwork))
            call dgesdd("A", m, n, A, ldA, Sigma, U, ldU, VT, ldVT, work, lwork, iwork, info)
            V = transpose(VT)
      end subroutine real_SVD
      

      subroutine real_LeastSquares(b, rank, A, rcond)
            !
            ! Solve the linear least squares problem min||Ax-b||.
            !
            real(F64), dimension(:, :), intent(inout) :: b
            integer, intent(out)                      :: rank
            real(F64), dimension(:, :), intent(inout) :: A
            real(F64), intent(in)                     :: rcond

            real(F64), dimension(1) :: work0
            real(F64), dimension(:), allocatable :: work
            real(F64), dimension(:), allocatable :: s
            integer :: m, n, nrhs
            integer :: lwork, info
            external :: dgelss

            m = size(A, dim=1)
            n = size(A, dim=2)
            nrhs = size(b, dim=2)
            allocate(s(max(1, min(m, n))))
            call dgelss(m, n, nrhs, A, m, b, m, s, rcond, rank, work0, -1, info)
            lwork = ceiling(work0(1))
            allocate(work(lwork))
            call dgelss(m, n, nrhs, A, m, b, m, s, rcond, rank, work, lwork, info)
            if (info /= 0) then
                  call msg("Linear least squares subroutine failed with Info="//str(info), MSG_ERROR)
                  error stop
            end if
      end subroutine real_LeastSquares


      subroutine real_LeastSquares_PivotedQR(b, A, P, Rank, QRThresh)
            !
            ! Solve linear least squares problem
            ! x = min(x') || Ax' - b ||
            ! by QR factorization with column pivoting
            ! of a rank-deficient matrix A.
            !
            real(F64), dimension(:, :), intent(inout) :: b
            real(F64), dimension(:, :), intent(inout) :: A
            integer, dimension(:), intent(out)        :: P
            integer, intent(out)                      :: Rank
            real(F64), intent(in)                     :: QRThresh

            real(F64), dimension(:), allocatable :: work
            real(F64), dimension(1) :: work0
            integer :: lwork, info
            integer :: M, N, NRHS, ldA, ldB
            
            external :: dgelsy

            P = 0
            M = size(A, dim=1)
            N = size(A, dim=2)
            NRHS = size(b, dim=2)
            ldA = M
            ldB = N
            if (size(b, dim=1) /= N) then
                  call msg("Invalid size of matrix b", MSG_ERROR)
                  error stop
            end if
            if (size(P) /= N) then
                  call msg("Invalid size of the pivots matrix P", MSG_ERROR)
                  error stop
            end if
            lwork = -1
            call dgelsy(M, N, NRHS, A, ldA, B, ldB, P, QRThresh, Rank, work0, lwork, info)
            lwork = ceiling(work0(1))
            allocate(work(lwork))
            call dgelsy(M, N, NRHS, A, ldA, B, ldB, P, QRThresh, Rank, work, lwork, info)            
      end  subroutine real_LeastSquares_PivotedQR


      subroutine real_invert_robust(a)
            !
            ! Invert a general matrix A using the Lapack subroutines for poorly
            ! conditioned linear systems.
            !
            real(F64), dimension(:, :), intent(inout) :: a

            real(F64), dimension(:, :), allocatable :: b
            integer :: n, k

            n = size(a, dim=1)
            allocate(b(n, n))
            b = ZERO
            do k = 1, n
                  b(k, k) = ONE
            end do
            call real_Axb_robust(b, a)
            a = b
      end subroutine real_invert_robust
      

      subroutine real_Axb_robust(b, a)
            !
            ! Solve a linear system Ax = b using the Lapack subroutines for poorly
            ! conditioned problems.
            !
            real(F64), dimension(:, :), intent(inout) :: b
            real(F64), dimension(:, :), intent(inout) :: a

            real(F64), dimension(:, :), allocatable :: af
            real(F64), dimension(:, :), allocatable :: x
            real(F64), dimension(:), allocatable :: r, c
            integer, dimension(:), allocatable :: ipiv
            real(F64), dimension(:), allocatable :: work
            integer, dimension(:), allocatable :: iwork
            character :: equed
            real(F64) :: rcond, rpvgrw
            integer :: n, lda
            integer :: info
            real(F64), dimension(:), allocatable :: berr
            integer :: nrhs = 1
            integer, parameter :: nparams = 0
            real(F64), dimension(:), allocatable :: params
            real(F64), dimension(:, :), allocatable :: err_bnds_comp, err_bnds_norm
            integer, parameter :: n_err_bnds=3
            external :: dgesvxx

            n = size(b, dim=1)
            nrhs = size(b, dim=2)
            lda = size(a, dim=1)
            allocate(af(n, n))
            allocate(ipiv(n))
            allocate(r(n))
            allocate(c(n))
            allocate(x(n, nrhs))
            allocate(berr(nrhs))
            allocate(err_bnds_comp(nrhs, n_err_bnds))
            allocate(err_bnds_norm(nrhs, n_err_bnds))
            allocate(params(nparams))
            allocate(work(4*n))
            allocate(iwork(n))
            
            call dgesvxx( &
                  "E", &           ! FACT
                  "N", &           ! TRANS
                  n, &             ! N
                  nrhs, &          ! NRHS
                  a, &             ! A
                  lda, &           ! LDA
                  af, &            ! AF
                  n, &             ! LDAF
                  ipiv, &          ! IPIV
                  equed, &         ! EQUED
                  r, &             ! R
                  c, &             ! C
                  b, &             ! B
                  n, &             ! ldb
                  x, &             ! X
                  n, &             ! LDX
                  rcond, &         ! RCOND
                  rpvgrw, &        ! RPVGRW
                  berr, &          ! BERR
                  n_err_bnds, &    ! N_ERR_BNDS
                  err_bnds_norm, & ! ERR_BNDS_NORM
                  err_bnds_comp, & ! ERR_BNDS_COMP
                  nparams, &       ! NPARAMS
                  params, &        ! PARAMS
                  work, &          ! WORK
                  iwork, &         ! IWORK
                  info)            ! INFO
            b = x
      end subroutine real_Axb_robust


      subroutine real_Axb_symmetric_sysv(b, A)
            !
            ! Solve a linear system Ax = b, where A is a symmetric
            ! nondefinite matrix.
            !
            real(F64), dimension(:, :), intent(inout) :: b
            real(F64), dimension(:, :), intent(inout) :: a

            integer :: N, NRHS, info, lwork
            integer, dimension(:), allocatable :: ipiv
            real(F64), dimension(:), allocatable :: work
            real(F64), dimension(1) :: work0
            
            external :: dsysv

            N = size(A, dim=1)
            NRHS = size(b, dim=2)
            allocate(ipiv(N))
            lwork = -1
            call dsysv("L", &
                  N, &
                  NRHS, &
                  A, &
                  N, &
                  ipiv, &
                  b, &
                  N, &
                  work0, &
                  lwork, &
                  info)
            lwork = max(1, ceiling(work0(1)))
            allocate(work(lwork))
            call dsysv("L", &
                  N, &
                  NRHS, &
                  A, &
                  N, &
                  ipiv, &
                  b, &
                  N, &
                  work, &
                  lwork, &
                  info)
            if (info /= 0) then
                  call msg("Linear system solver returned with info="//str(info), MSG_ERROR)
                  error stop
            end if
      end subroutine real_Axb_symmetric_sysv


      subroutine real_Axb_nonsymmetric_gesv(b, A)
            !
            ! Solve a linear system Ax = b, where A is a nonsymmetric
            ! matrix.
            !
            real(F64), dimension(:, :), intent(inout) :: b
            real(F64), dimension(:, :), intent(inout) :: a

            integer :: N, NRHS, info
            integer, dimension(:), allocatable :: ipiv
            
            external :: dgesv

            N = size(A, dim=1)
            NRHS = size(b, dim=2)
            allocate(ipiv(N))
            call dgesv(N, &
                  NRHS, &
                  A, &
                  N, &
                  ipiv, &
                  b, &
                  N, &
                  info)
            if (info /= 0) then
                  call msg("Linear system solver returned with info="//str(info), MSG_ERROR)
                  error stop
            end if
      end subroutine real_Axb_nonsymmetric_gesv
      
      
      subroutine real_Pack(MPacked, MFull, NVecs)
            !
            ! Copy a symmetric matrix into a packed one-dimensional array
            !
            real(F64), dimension(:), intent(out)   :: MPacked
            real(F64), dimension(:, :), intent(in) :: MFull
            integer, intent(in)                    :: NVecs

            integer :: ldM, info
            external :: dtrttp

            ldM = size(MFull, dim=1)
            call dtrttp("L", NVecs, MFull, ldM, MPacked, info)
      end subroutine real_Pack


      subroutine real_Unpack(MFull, MPacked, NVecs)
            !
            ! Unpack a packed one-dimensional array into a symmetric matrix.
            ! The upper triangle of the unpacked matrix will be filled
            ! with zeros.
            !
            real(F64), dimension(:, :), intent(out) :: MFull
            real(F64), dimension(:), intent(in)     :: MPacked
            integer, intent(in)                     :: NVecs

            integer :: ldM, info
            external :: dtpttr

            ldM = size(MFull, dim=1)
            !
            ! Fill the upper triangle of MFull with zeros because it's not referenced
            ! by the unpacking subroutine
            !
            MFull = ZERO
            call dtpttr("L", NVecs, MPacked, MFull, ldM, info)
      end subroutine real_Unpack
      
      
      subroutine real_LogDet_query(lwork, lda, n)
            !
            ! Query the size of the work array for the Log(Det(X)) subroutine.
            !
            integer, intent(out) :: lwork
            integer, intent(in) :: lda
            integer, intent(in) :: n

            real(F64), dimension(1) :: tau, a, work
            integer :: info
            external :: dgeqrf

            call dgeqrf(n, n, a, lda, tau, work, -1, info)
            lwork = ceiling(work(1)) + n
      end subroutine real_LogDet_query


      subroutine real_LogDet(LogDet, a, lda, n, work, lwork)
            !
            ! Compute Log(Det(A)) where A is a real symmetric matrix.
            ! Note that both the upper and lower triangle of A are referenced.
            !
            real(F128), intent(out) :: LogDet
            real(F64), dimension(lda, *), intent(inout) :: a
            integer, intent(in) :: lda
            integer, intent(in) :: n
            real(F64), dimension(*), intent(in) :: work
            integer, intent(in) :: lwork

            integer :: k
            integer :: info
            real(F128) :: t
            external :: dgeqrf
            !
            ! The TAU array of the scalar factors of elementary
            ! reflectors is stored in the first k elements of WORK.
            ! The size of WORK is still optimal because the query
            ! subroutine adds k to the base value of LWORK.
            !
            call dgeqrf(n, n, a, lda, work(1:n), work(n+1:lwork), lwork-n, info)
            if (info .ne. 0) then
                  call msg("DGEQRF returned error code " // str(info), MSG_ERROR)
                  stop
            end if
            LogDet = 0.0_F128
            do k = 1, n
                  t = real(abs(a(k, k)), F128)
                  LogDet = LogDet + log(t)
            end do
      end subroutine real_LogDet
      

      subroutine real_QR_query(lwork, lda, m, n)
            !
            ! Query the size of the work array for the QR factorization subroutine
            !
            integer, intent(out) :: lwork
            integer, intent(in) :: lda
            integer, intent(in) :: m, n

            real(F64), dimension(1) :: tau, a, work
            integer :: info, lwork1, lwork2
            integer :: k
            external :: dgeqrf, dorgqr

            call dgeqrf(m, n, a, lda, tau, work, -1, info)
            lwork1 = ceiling(work(1))
            k = min(m, n)
            call dorgqr(m, n, k, a, lda, tau, work, -1, info)
            lwork2 = ceiling(work(1))
            lwork = k + max(lwork1, lwork2)
      end subroutine real_QR_query

      
      subroutine real_QR(a, lda, m, n, qrwork, lqrwork)
            !
            ! Compute the Q matrix of the QR matrix factorization.
            ! The size of the work array *must* be computed using
            ! the query subroutine.
            !
            real(F64), dimension(lda, *), intent(inout) :: a
            integer, intent(in) :: lda
            integer, intent(in) :: m, n
            real(F64), dimension(*), intent(out) :: qrwork
            integer, intent(in) :: lqrwork
            
            integer :: k
            integer :: info, lwork
            external :: dgeqrf, dorgqr

            k = min(m, n)
            lwork = lqrwork - k
            !
            ! The TAU array of the scalar factors of elementary
            ! reflectors is stored in the first k elements of WORK.
            ! The size of WORK is still optimal because the query
            ! subroutine adds k to the base value of LWORK.
            !
            associate(tau => qrwork(1:k), work => qrwork(k+1:lqrwork))
                  call dgeqrf(m, n, a, lda, tau, work, lwork, info)
                  if (info .ne. 0) then
                        call msg("DGEQRF returned error code " // str(info), MSG_ERROR)
                        stop
                  end if
                  call dorgqr(m, n, k, a, lda, tau, work, lwork, info)
                  if (info .ne. 0) then
                        call msg("DORGQR returned error code " // str(info), MSG_ERROR)
                        stop
                  end if
            end associate
      end subroutine real_QR
      

      subroutine real_vwT(a, v, w, alpha)
            !
            ! A <- alpha * v*w**T + A
            !
            real(F64), dimension(:, :), intent(inout) :: a
            real(F64), dimension(:), intent(in)       :: v
            real(F64), dimension(:), intent(in)       :: w
            real(F64), intent(in)                     :: alpha

            integer :: m, n
            external :: dger

            m = size(a, dim=1)
            n = size(a, dim=2)
            call dger(m, n, alpha, v, 1, w, 1, a, m)
      end subroutine real_vwT


      subroutine real_vwT_x(a, lda, v, w, m, n, alpha)
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
      end subroutine real_vwT_x
      
      
      subroutine real_aTba(c, a, b, scratch)
            real(F64), dimension(:, :), contiguous, intent(out) :: c
            real(F64), dimension(:, :), contiguous, intent(in)  :: a
            real(F64), dimension(:, :), contiguous, intent(in)  :: b
            real(F64), dimension(:, :), contiguous, intent(out) :: scratch

            integer :: n, k
            integer :: lda, ldb, ldc
            external :: dgemm

            n = size(a, dim=2)
            k = size(b, dim=2)
            lda = size(a, dim=1)
            ldb = size(b, dim=1)
            ldc = size(c, dim=1)

            if ((size(c, dim=1) .ne. size(a, dim=2)) .or. &
                  (size(a, dim=1) .ne. size(b, dim=2)) .or. &
                  (size(b, dim=1) .ne. size(b, dim=2))) then
                  call msg("Inconsistent dimensions on entry to real_aTba", MSG_ERROR)
                  stop
            end if

            if (size(scratch) < n * k) then
                  call msg("Scratch too small for real_aTba", MSG_ERROR)
                  stop
            end if
            !
            ! SCRATCH <- A^T B
            !
            call dgemm("T", "N", n, k, k, ONE, a, lda, b, ldb, ZERO, scratch, n)
            !
            ! C <- SCRATCH A
            !
            call dgemm("N", "N", n, n, k, ONE, scratch, n, a, lda, ZERO, c, ldc)
      end subroutine real_aTba


      subroutine real_abaT(c, a, b, scratch)
            real(F64), dimension(:, :), contiguous, intent(out) :: c
            real(F64), dimension(:, :), contiguous, intent(in)  :: a
            real(F64), dimension(:, :), contiguous, intent(in)  :: b
            real(F64), dimension(:, :), contiguous, intent(out) :: scratch

            integer :: n, k
            integer :: lda, ldb, ldc
            external :: dgemm

            scratch = ZERO
            c = ZERO
            n = size(a, dim=1)
            k = size(b, dim=1)
            lda = size(a, dim=1)
            ldb = size(b, dim=1)
            ldc = size(c, dim=1)

            if ((size(c, dim=1) .ne. size(a, dim=1)) .or. &
                  (size(a, dim=2) .ne. size(b, dim=2)) .or. &
                  (size(b, dim=1) .ne. size(b, dim=2))) then
                  call msg("Inconsistent dimensions on entry to real_abaT", MSG_ERROR)
                  stop
            end if

            if (size(scratch) < n * k) then
                  call msg("Scratch too small for real_abaT", MSG_ERROR)
                  stop
            end if
            !
            ! SCRATCH <- B A^T
            !
            call dgemm("N", "T", k, n, k, ONE, b, ldb, a, lda, ZERO, scratch, k)
            !
            ! C <- A SCRATCH
            !
            call dgemm("N", "N", n, n, k, ONE, a, lda, scratch, k, ZERO, c, ldc)
      end subroutine real_abaT
      

      subroutine real_aTbc(d, a, b, c)
            !
            ! D <- A**T B C
            !
            real(F64), dimension(:, :), intent(out) :: d
            real(F64), dimension(:, :), intent(in)  :: a
            real(F64), dimension(:, :), intent(in)  :: b
            real(F64), dimension(:, :), intent(in)  :: c
            !
            ! Compute D <- A**H B C
            !
            integer :: u, v
            real(F64), dimension(:, :), allocatable :: w

            u = size(a, dim=2)
            v = size(b, dim=2)
            allocate(w(u, v))
            !
            ! W <- A**T B
            !
            call real_aTb(w, a, b)
            !
            ! D <- W C
            !
            call real_ab(d, w, c)
      end subroutine real_aTbc


      subroutine real_aTb(c, a, b)
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
                  call msg("real_aTb: inconsistent dimensions of matrices A and B", MSG_ERROR)
                  stop
            end if
            if ((u .ne. l) .or. (v .ne. n)) then
                  call msg("real_aTb: inconsistent dimensions of matrix C", MSG_ERROR)
                  stop
            end if
            call dgemm("C", "N", u, v, k, alpha, a, k, b, m, beta, c, u)
      end subroutine real_aTb


      subroutine real_aTb_x(c, ldc, a, lda, b, ldb, m, n, k, alpha, beta)
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
      end subroutine real_aTb_x


      subroutine real_Av(w, A, v)
            real(F64), dimension(:), intent(out)   :: w
            real(F64), dimension(:, :), intent(in) :: A
            real(F64), dimension(:), intent(in)    :: v

            real(F64), parameter :: alpha = ONE
            real(F64), parameter :: beta = ZERO
            integer :: m, n
            external :: dgemv

            m = size(A, dim=1)
            n = size(A, dim=2)
            
            call dgemv("N", m, n, alpha, A, m, v, 1, beta, w, 1)
      end subroutine real_Av


      subroutine real_ATv(w, A, v)
            real(F64), dimension(:), intent(out)   :: w
            real(F64), dimension(:, :), intent(in) :: A
            real(F64), dimension(:), intent(in)    :: v

            real(F64), parameter :: alpha = ONE
            real(F64), parameter :: beta = ZERO
            integer :: m, n
            external :: dgemv

            m = size(A, dim=1)
            n = size(A, dim=2)
            
            call dgemv("T", m, n, alpha, A, m, v, 1, beta, w, 1)
      end subroutine real_ATv
      
      
      subroutine real_av_x(w, a, lda, v, m, n, alpha, beta)
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
      end subroutine real_av_x

      
      subroutine real_aTv_x(w, a, lda, v, m, n, alpha, beta)
            !
            ! Perform matrix-vector multiplication w(1:n) = alpha * A(1:m,1:n)**T v(1:m) + beta * w(1:n)
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
      end subroutine real_aTv_x


      subroutine real_vw_x(d, v, w, n)
            !
            ! Perform vector dot product d = vw
            !
            real(F64), intent(out)              :: d
            real(F64), dimension(*), intent(in) :: v
            real(F64), dimension(*), intent(in) :: w
            integer, intent(in)                 :: n

            real(F64) :: ddot
            external :: ddot

            d = ddot(n, v, 1, w, 1)
      end subroutine real_vw_x


      subroutine real_ab(c, a, b)
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
                  call msg("real_ab: inconsistent dimensions of matrices A and B", MSG_ERROR)
                  stop
            end if
            if ((u .ne. k) .or. (v .ne. n)) then
                  call msg("real_ab: inconsistent dimensions of matrix C", MSG_ERROR)
                  stop
            end if
            call dgemm("N", "N", u, v, l, alpha, a, k, b, m, beta, c, u)
      end subroutine real_ab


      subroutine real_ab_x(c, ldc, a, lda, b, ldb, m, n, k, alpha, beta)
            !
            ! Compute C(1:m, 1:n) <- alpha*A(1:m, 1:k) B(1:k, 1:n) + beta*C(1:m, 1:n)
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
                  alpha0 = ONE
            end if

            if (present(beta)) then
                  beta0 = beta
            else
                  beta0 = ZERO
            end if

            call dgemm("N", "N", m, n, k, alpha0, a, lda, b, ldb, beta0, c, ldc)
      end subroutine real_ab_x


      subroutine real_aTbT(c, a, b, alpha, beta)
            real(F64), dimension(:, :), intent(inout) :: c
            real(F64), dimension(:, :), intent(in)    :: a
            real(F64), dimension(:, :), intent(in)    :: b
            real(F64), optional, intent(in)           :: alpha
            real(F64), optional, intent(in)           :: beta

            integer :: m, n, k
            real(F64) :: alpha0, beta0
            
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

            m = size(c, dim=1)
            n = size(c, dim=2)
            k = size(a, dim=1)
            call real_aTbT_x(c, m, a, k, b, n, m, n, k, alpha0, beta0)
      end subroutine real_aTbT
      

      subroutine real_aTbT_x(c, ldc, a, lda, b, ldb, m, n, k, alpha, beta)
            !
            ! Compute C(1:m, 1:n) <- alpha*A(1:k, 1:m)**T B(1:n, 1:k)**T + beta*C(1:m, 1:n)
            !
            real(F64), dimension(ldc, *), intent(inout) :: c
            integer, intent(in)                         :: ldc
            real(F64), dimension(lda, *), intent(in)    :: a
            integer, intent(in)                         :: lda
            real(F64), dimension(ldb, *), intent(in)    :: b
            integer, intent(in)                         :: ldb
            integer, intent(in)                         :: m
            integer, intent(in)                         :: n
            integer, intent(in)                         :: k
            real(F64), optional, intent(in)             :: alpha
            real(F64), optional, intent(in)             :: beta
            
            real(F64) :: alpha0, beta0
            external :: dgemm

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

            call dgemm("C", "C", m, n, k, alpha0, a, lda, b, ldb, beta0, c, ldc)
      end subroutine real_aTbT_x
      
      
      subroutine real_abT(c, a, b)
            !
            ! Compute C(1:u, 1:v) <- A(1:k, 1:l) B(1:m, 1:n)**T
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
            if (n .ne. l) then
                  call msg("real_abT: inconsistent dimensions of matrices A and B", MSG_ERROR)
                  error stop
            end if
            if ((u .ne. k) .or. (v .ne. m)) then
                  call msg("real_abT: inconsistent dimensions of matrix C", MSG_ERROR)
                  error stop
            end if
            call dgemm("N", "C", u, v, l, alpha, a, k, b, m, beta, c, u)
      end subroutine real_abT


      subroutine real_abT_x(c, ldc, a, lda, b, ldb, m, n, k, alpha, beta)
            !
            ! Compute C(1:m, 1:n) <- alpha * A(1:m, 1:k) B(1:k, 1:n)**T + beta * C(1:m,1:n)
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
      end subroutine real_abT_x

      
      subroutine linear_system(y, A)
            !
            ! Solve a linear system of equations Ax = y.
            ! The solution overwrites the right-hand side vector y.
            !
            real(F64), dimension(:), contiguous, intent(inout)    :: y
            real(F64), dimension(:, :), contiguous, intent(inout) :: A

            integer :: n, info
            integer, dimension(:), allocatable :: ipiv
            external :: dgetrf, dgetrs
            !
            ! 1. Compute LU factorization of A
            ! 2. Solve linear system A * x = y
            !
            n = size(A, dim=1)
            allocate(ipiv(n))
            call dgetrf(n, n, A, n, ipiv, info)
            call dgetrs("N", n, 1, A, n, ipiv, y, n, info)
            if (info .ne. 0) then
                  call msg("Error in linear system solver", MSG_ERROR)
                  stop
            end if
      end subroutine linear_system


      subroutine nonsymmetric_eigenproblem(wr, wi, vl, vr, A)
            real(F64), dimension(:), intent(out)      :: wr
            real(F64), dimension(:), intent(out)      :: wi
            real(F64), dimension(:, :), intent(out)   :: vl
            real(F64), dimension(:, :), intent(out)   :: vr
            real(F64), dimension(:, :), intent(inout) :: A
            
            real(F64), dimension(1) :: work0
            real(F64), dimension(:), allocatable :: work
            integer :: lwork, info
            integer :: n

            external :: dgeev
            
            n = size(A, dim=1)
            lwork = -1
            call dgeev("V", "V", n, A, n, wr, wi, vl, n, vr, n, work0, lwork, info)
            lwork = ceiling(work0(1))
            allocate(work(lwork))
            call dgeev("V", "V", n, A, n, wr, wi, vl, n, vr, n, work, lwork, info)
            if (info /= 0) then
                  call msg("Nonsymmetric matrix eigendecompositino failed with info=" // str(info), MSG_ERROR)
                  error stop
            end if
      end subroutine nonsymmetric_eigenproblem

      
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
            type(tclock) :: t_eigen
            external :: dsyevd

            call clock_start(t_eigen)
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
                  call msg("Eigensolver for symmetric matrices returned error code (" // str(info) // ")", MSG_ERROR)
                  stop
            end if
            deallocate(work)
            deallocate(iwork)
      end subroutine symmetric_eigenproblem
      

      subroutine real_bta_rect(a, b)
            ! ---------------------------------------------------------
            ! A(1:M, 1:N) <- B(1:K, 1:M)^T A(1:K, 1:N),
            ! M, K << N.
            ! Use this subroutine only if M and K are small numbers.
            ! Because the matrix A stores both input and output, in
            ! some scenarios this subroutine is more memory-efficient
            ! than GEMM. A small array is allocated to temporarily
            ! store a single column of A.
            ! ---------------------------------------------------------
            ! A
            !    On entry, a rectangular matrix of dimension K x N.
            !    On exit, a rectangular matrix of dimension M x N.
            ! B
            !    On entry, a matrix of dimension K x N.
            !
            real(F64), dimension(:, :), intent(inout) :: a
            real(F64), dimension(:, :), intent(in)    :: b

            integer :: m, n, k
            integer :: p, q, w
            real(F64), dimension(:), allocatable :: aw
            real(F64) :: apw
            
            n = size(a, dim=2)
            k = size(a, dim=1)
            m = size(b, dim=2)

            if (size(b, dim=1) .ne. k) then
                  call msg("ERROR: INVALID ARGUMENT PASSED TO REAL_BTA_RECT", &
                        priority=MSG_ERROR)
                  stop
            end if

            !$omp parallel private(w, aw, p, q, apw) &
            !$omp shared(a, b, m, n, k) &
            !$omp default(none)
            allocate(aw(k))
            !$omp do schedule(static)
            do w = 1, n
                  aw(:) = a(1:k, w)
                  do p = 1, m
                        apw = ZERO
                        do q = 1, k
                              apw = apw + b(q, p) * aw(q)
                        end do
                        a(p, w) = apw
                  end do
            end do
            !$omp end do nowait
            deallocate(aw)
            !$omp end parallel
      end subroutine real_bta_rect


      subroutine real_vta_rect_scal(w, a, v, alpha)
            ! -------------------------------------------------------
            ! W(1:N) <- ALPHA * V(1:K)^T A(1:K, 1:N),
            ! K << N.
            ! A is a rectangular matrix.
            ! Use this subroutine only if A is a K x N matrix
            ! where K IS A SMALL NUMBER!
            ! -------------------------------------------------------
            ! A
            !    A rectangular matrix of dimension K x N.
            ! V
            !    A vector of dimension K.
            ! W 
            !    On exit, W = V^T A. A vector of dimension N.
            ! ALPHA
            !    Input, a scalar scaling factor.
            !
            real(F64), dimension(:), intent(out)   :: w
            real(F64), dimension(:, :), intent(in) :: a
            real(F64), dimension(:), intent(in)    :: v
            real(F64), intent(in)                  :: alpha
            
            integer :: n, k
            integer :: q, r
            real(F64) :: wr
            
            n = size(a, dim=2)
            k = size(a, dim=1)

            if (size(w, dim=1) .ne. n .or. size(v, dim=1) .ne. k) then
                  call msg("ERROR: INVALID ARGUMENT PASSED TO REAL_VTA_RECT", &
                        priority=MSG_ERROR)
                  stop
            end if

            !$omp parallel private(r, wr, q) &
            !$omp shared(a, v, w, n, k, alpha) &
            !$omp default(none)
            !$omp do schedule(static)
            do r = 1, n
                  wr = ZERO
                  do q = 1, k
                        wr = wr + v(q) * a(q, r)
                  end do
                  w(r) = alpha * wr
            end do
            !$omp end do nowait
            !$omp end parallel
      end subroutine real_vta_rect_scal


      subroutine real_vta_rect_noscal(w, a, v)
            ! -------------------------------------------------------
            ! W(1:N) <- V(1:K)^T A(1:K, 1:N),
            ! K << N.
            ! A is a rectangular matrix.
            ! Use this subroutine only if A is a K x N matrix
            ! where K IS A SMALL NUMBER!
            ! -------------------------------------------------------
            ! A
            !    A rectangular matrix of dimension K x N.
            ! V
            !    A vector of dimension K.
            ! W 
            !    On exit, W = V^T A. A vector of dimension N.
            ! ALPHA
            !    Input, a scalar scaling factor.
            !
            real(F64), dimension(:), intent(out)   :: w
            real(F64), dimension(:, :), intent(in) :: a
            real(F64), dimension(:), intent(in)    :: v
            
            integer :: n, k
            integer :: q, r
            real(F64) :: wr
            
            n = size(a, dim=2)
            k = size(a, dim=1)

            if (size(w, dim=1) .ne. n .or. size(v, dim=1) .ne. k) then
                  call msg("ERROR: INVALID ARGUMENT PASSED TO REAL_VTA_RECT", &
                        priority=MSG_ERROR)
                  stop
            end if

            !$omp parallel private(r, wr, q) &
            !$omp shared(a, v, w, n, k) &
            !$omp default(none)
            !$omp do schedule(static)
            do r = 1, n
                  wr = ZERO
                  do q = 1, k
                        wr = wr + v(q) * a(q, r)
                  end do
                  w(r) = wr
            end do
            !$omp end do nowait
            !$omp end parallel
      end subroutine real_vta_rect_noscal


      subroutine real_vtapw_rect(w, a, v, alpha)
            ! -------------------------------------------------------
            ! W(1:N) <- W(1:N) + ALPHA * V(1:K)^T A(1:K, 1:N)
            ! K << N.
            ! A is a rectangular matrix.
            ! Use this subroutine only if A is a K x N matrix
            ! where K IS A SMALL NUMBER!
            ! -------------------------------------------------------
            ! A
            !    A rectangular matrix of dimension K x N.
            ! V
            !    A vector of dimension K.
            ! W 
            !    On exit, W = V^T A. A vector of dimension N.
            ! ALPHA
            !    Input, a scalar scaling factor.
            !
            real(F64), dimension(:), intent(inout) :: w
            real(F64), dimension(:, :), intent(in) :: a
            real(F64), dimension(:), intent(in)    :: v
            real(F64), intent(in)                  :: alpha
            
            integer :: n, k
            integer :: q, r
            real(F64) :: wr
            
            n = size(a, dim=2)
            k = size(a, dim=1)

            if (size(w, dim=1) .ne. n .or. size(v, dim=1) .ne. k) then
                  call msg("ERROR: INVALID ARGUMENT PASSED TO REAL_VTA_RECT", &
                        priority=MSG_ERROR)
                  stop
            end if

            !$omp parallel private(r, wr, q) &
            !$omp shared(a, v, w, n, k, alpha) &
            !$omp default(none)
            !$omp do schedule(static)
            do r = 1, n
                  wr = ZERO
                  do q = 1, k
                        wr = wr + v(q) * a(q, r)
                  end do
                  w(r) = w(r) + alpha * wr
            end do
            !$omp end do nowait
            !$omp end parallel
      end subroutine real_vtapw_rect


      subroutine real_l2norm(l2norm, a)
            !
            ! L2NORM <- SQRT(A^T A)
            !
            real(F64), intent(out)              :: l2norm
            real(F64), dimension(:), intent(in) :: a

            real(F64) :: l2norm_partial
            integer :: k, l, n, r
            integer, parameter :: blocksize = 2**10

            n = size(a, dim=1)
            r = modulo(n, blocksize)
            l2norm = ZERO
            if (r .ne. 0) then
                  l2norm_partial = ZERO
                  do l = 1, r
                        l2norm_partial = l2norm_partial + a(l)**2
                  end do
                  l2norm = l2norm + l2norm_partial
            end if
            if (n >= blocksize) then
                  !$omp parallel reduction(+:l2norm) &
                  !$omp private(k, l, l2norm_partial) &
                  !$omp shared(a, n, r) default(none)
                  !$omp do
                  do k = r+1, n, blocksize
                        l2norm_partial = ZERO
                        do l = 1, blocksize
                              l2norm_partial = l2norm_partial + a(k-1+l)**2
                        end do
                        l2norm = l2norm + l2norm_partial
                  end do
                  !$omp end do nowait
                  !$omp end parallel
            end if

            l2norm = sqrt(l2norm)
      end subroutine real_l2norm


      subroutine real_dot(dot, a, b)
            !
            ! DOT <- A^T B
            !
            real(F64), intent(out)              :: dot
            real(F64), dimension(:), intent(in) :: a
            real(F64), dimension(:), intent(in) :: b

            real(F64) :: dot_partial
            integer :: k, l, n, r
            integer, parameter :: blocksize = 2**10

            if (size(a) .ne. size(b)) then
                  call msg("INCONSISTENT DIMENSIONS ON ENTRY TO REAL_DOT", &
                        priority=MSG_ERROR)
                  call imsg("SIZE(A)", size(a), MSG_ERROR)
                  call imsg("SIZE(B)", size(b), MSG_ERROR)
                  stop
            end if

            n = size(a)
            r = modulo(n, blocksize)
            dot = ZERO

            if (r .ne. 0) then
                  dot_partial = ZERO
                  do l = 1, r
                        dot_partial = dot_partial + a(l) * b(l)
                  end do
                  dot = dot + dot_partial
            end if
            
            if (n >= blocksize) then
                  !$omp parallel reduction(+:dot) &
                  !$omp private(k, l, dot_partial) &
                  !$omp shared(a, b, n, r) default(none)
                  !$omp do
                  do k = r+1, n, blocksize
                        dot_partial = ZERO
                        do l = 1, blocksize
                              dot_partial = dot_partial + a(k-1+l) * b(k-1+l)
                        end do
                        dot = dot + dot_partial
                  end do
                  !$omp end do nowait
                  !$omp end parallel
            end if
      end subroutine real_dot


      subroutine real_scal(a, alpha)
            !
            ! A <- ALPHA * A
            !
            real(F64), dimension(:), intent(inout) :: a
            real(F64), intent(in)                  :: alpha

            integer :: k, n

            n = size(a, dim=1)

            !$omp parallel shared(n, a, alpha) private(k) default(none)
            !$omp do schedule(static)
            do k = 1, n
                  a(k) = alpha * a(k)
            end do
            !$omp end do nowait
            !$omp end parallel
      end subroutine real_scal


      subroutine real_axpy(y, alpha, x)
            !
            ! Y <- ALPHA * X + Y
            !
            real(F64), dimension(:), intent(inout) :: y
            real(F64), intent(in)                  :: alpha
            real(F64), dimension(:), intent(in)    :: x

            integer :: k, n

            if (size(x) .ne. size(y)) then
                  call msg("INCONSISTENT DIMENSIONS ON ENTRY TO REAL_AXPY", &
                        priority=MSG_ERROR)
                  call imsg("SIZE(X)", size(x), MSG_ERROR)
                  call imsg("SIZE(Y)", size(y), MSG_ERROR)
                  stop
            end if

            n = size(y)

            !$omp parallel shared(n, x, y, alpha) private(k) default(none)
            !$omp do schedule(static)
            do k = 1, n
                  y(k) = y(k) + alpha * x(k)
            end do
            !$omp end do nowait
            !$omp end parallel
      end subroutine real_axpy


      subroutine real_conditional_gramschmidt(colwise, a, thresh, max_overlap)
            logical, intent(in)                       :: colwise
            real(F64), dimension(:, :), intent(inout) :: a
            real(F64), intent(in)                     :: thresh
            real(F64), intent(out)                    :: max_overlap
            
            real(F64) :: s_ij
            integer :: m, i, j
            
            max_overlap = ZERO
            if (colwise) then
                  m = size(a, dim=2)
                  outer1: do j = 1, m
                        do i = j + 1, m
                              call real_dot(s_ij, a(:, i), a(:, j))
                              max_overlap = max(max_overlap, abs(s_ij))
                        end do
                  end do outer1
            else
                  m = size(a, dim=1)
                  outer2: do j = 1, m
                        do i = j + 1, m
                              call real_dot(s_ij, a(i, :), a(j, :))
                              max_overlap = max(max_overlap, abs(s_ij))
                        end do
                  end do outer2
            end if

            if (max_overlap > thresh) then
                  call real_gramschmidt(colwise, a, 0)
            end if
      end subroutine real_conditional_gramschmidt


      subroutine real_gramschmidt(colwise, a, k)
            ! -----------------------------------------------------------------------
            ! Gram-Schmidt orthonormalization.
            ! -----------------------------------------------------------------------
            ! COLWISE
            !     .TRUE. if vectors are stored in columns of A.
            !     .FALSE. if vectors are stored in rows of A.
            !
            ! A
            !     On entry, the columns of A are the vectors to be orthonormalized.
            !     On exit, the matrix A contains K + L orthonormal columns.
            ! K
            !     the first K vectors be assumed to be orthonormal.
            !     K==0 if A contains no orthonormal vectors.
            !
            logical, intent(in)                       :: colwise
            real(F64), dimension(:, :), intent(inout) :: a
            integer, intent(in)                       :: k

            integer :: n, l
            integer :: i, j
            real(F64) :: t, l2norm

            if (colwise) then
                  n = size(a, dim=1)
                  l = size(a, dim=2) - k
            else
                  n = size(a, dim=2)
                  l = size(a, dim=1) - k
            end if

            if (l < 0 .or. k < 0) then
                  call msg("INVALID ARGUMENT PASSED TO REAL_GRAMSCHMIDT", &
                        priority=MSG_ERROR)
                  stop
            end if

            if (l == 0) return
            !
            ! Normalize each vector which is not yet normalized
            !
            do i = k + 1, k + l
                  if (colwise) then
                        call real_l2norm(l2norm, a(:, i))
                        call real_scal(a(:, i), ONE/l2norm)
                  else
                        call real_l2norm(l2norm, a(i, :))
                        call real_scal(a(i, :), ONE/l2norm)
                  end if
            end do
            !
            ! Project out non-orthogonal components from each
            ! of the recently appended coulumns
            !
            do i = k + 1, k + l
                  do j = 1, i - 1
                        if (colwise) then
                              call real_dot(t, a(:, j), a(:, i))
                              call real_axpy(a(:, i), -t, a(:, j))
                        else
                              call real_dot(t, a(j, :), a(i, :))
                              call real_axpy(a(i, :), -t, a(j, :))
                        end if
                  end do

                  if (colwise) then
                        call real_l2norm(l2norm, a(:, i))
                        call real_scal(a(:, i), ONE/l2norm)
                  else
                        call real_l2norm(l2norm, a(i, :))
                        call real_scal(a(i, :), ONE/l2norm)
                  end if
            end do
      end subroutine real_gramschmidt


      subroutine real_biorth(colwise, abar, a, u)
            ! ----------------------------------------------------------
            ! First U columns/rows of A (ABAR) contain the old vectors
            ! (from previous iterations), which are assumed to be
            ! biorthonormalized on entry to this subroutine:
            !
            ! <ABAR(:,p)|A(:,q)> = \delta_{p, q} for every p, q <= U
            !
            ! X = [............... ...................]
            !     |<----- U ----->|<------- V ------->|
            !      Old, already      New, to be
            !      biorthonormal     biorthonormalized
            !
            ! X = A, ABAR
            !
            ! The next V columns/rows of A (ABAR) correspond
            ! to the recently added vectors which are not yet
            ! biorthonormal either amongst themselves or with
            ! respect to the old U vectors. After the
            ! biorthonormalization the following condition
            !
            ! <ABAR(:,p)|A(:,q)> = \delta_{p, q}
            ! for every p, q <= U+V
            !
            ! is satisfied.
            ! ----------------------------------------------------------
            ! 1. Hirao, K. and Nakatsuji, H., A generalization of the
            !    Davidson's Method to Large Nonsymmetric Eigenvalue
            !    Problems, J. Comp. Phys. 45, 246 (1982).
            ! ----------------------------------------------------------
            ! COLWISE
            !           .TRUE. if vectors are stored in columns 
            !           .FALSE. if vectors are stored in rows
            !
            ! ABAR    
            !           Left vectors.
            ! A       
            !           Right vectors.
            ! U      
            !           Number of vectors assumed biorthonormal.
            !
            logical, intent(in)                       :: colwise
            real(F64), dimension(:, :), intent(inout) :: a
            real(F64), dimension(:, :), intent(inout) :: abar
            integer, intent(in)                       :: u

            integer :: n, v
            integer :: p, q
            real(F64) :: t, scala, scalb, sab


            if ((size(a, dim=1) .ne. size(abar, dim=1)) .or. &
                  (size(a, dim=2) .ne. size(abar, dim=2))) then
                  call msg("INCONSISTENT DIMENSIONS ON ENTRY TO REAL_BIORTH", &
                        priority=MSG_ERROR)
                  stop
            end if

            if (colwise) then
                  n = size(a, dim=1)
                  v = size(a, dim=2) - u
            else
                  n = size(a, dim=2)
                  v = size(a, dim=1) - u
            end if

            if (u < 0 .or. v < 0) then
                  call msg("INVALID ARGUMENT PASSED TO REAL_BIORTH", &
                        priority=MSG_ERROR)
                  stop
            end if

            if (v == 0) return
            !
            ! Loop over the V new vectors 
            !
            do p = u + 1, u + v
                  !
                  ! Schmidt biorthogonalization: Eq. 11 in [1].
                  ! Project out q-th vector:
                  ! ABAR(:, P) <- (1 - ABAR(:, Q) A(:, Q)T^) ABAR(:, P)
                  !
                  do q = 1, p - 1
                        if (colwise) then
                              call real_dot(t, a(:, q), abar(:, p))
                              call real_axpy(abar(:, p), -t, abar(:, q))
                        else
                              call real_dot(t, a(q, :), abar(p, :))
                              call real_axpy(abar(p, :), -t, abar(q, :))
                        end if
                  end do
                  !
                  ! Schmidt biorthogonalization: Eq. 11 in [1].
                  ! Project out q-th vector:
                  ! A(:, P) <- (1 - A(:, Q) ABAR(:, Q)^T) A(:, P)
                  !
                  do q = 1, p - 1
                        if (colwise) then
                              call real_dot(t, abar(:, q), a(:, p))
                              call real_axpy(a(:, p), -t, a(:, q))
                        else
                              call real_dot(t, abar(q, :), a(p, :))
                              call real_axpy(a(p, :), -t, a(q, :))
                        end if
                  end do
                  !
                  ! Bi-normalize pair of left and right vector
                  ! according to Eq. 13 in [1]
                  !
                  if (colwise) then
                        call real_dot(sab, abar(:, p), a(:, p))
                        scalb = ONE / sqrt(abs(sab))
                        !
                        ! Sign of one of the pair of vectors is changed
                        ! if the dot product if negative (see the comment
                        ! below Eq. 15 in [1]).
                        !
                        scala = sign(scalb, sab)
                        call real_scal(abar(:, p), scala)
                        call real_scal(a(:, p), scalb)
                  else
                        call real_dot(sab, abar(p, :), a(p, :))
                        scalb = ONE / sqrt(abs(sab))
                        scala = sign(scalb, sab)
                        call real_scal(abar(p, :), scala)
                        call real_scal(a(p, :), scalb)
                  end if
            end do
      end subroutine real_biorth


      subroutine real_smfill(a)
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
      end subroutine real_smfill
end module real_linalg
