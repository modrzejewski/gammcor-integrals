module spherh
      use math_constants
      use arithmetic

      implicit none

contains

      pure subroutine xyzint(w, l, lx0, ly0, lz0)
            ! --------------------------------------------------------------
            ! Compute the integral of x^{lx+lx0} y^{ly+ly0} z^{lz+lz0} over
            ! the unit sphere:
            ! Int_Omega x^{lx+lx0} y^{ly+ly0} z^{lz+lz0} dOmega.
            ! (See Eq. 15 in [1]) The tntegrals are computed for the fixed
            ! lx0, ly0, lz0, and for all (LX,LY,LZ) for which LX+LY+LZ=L.
            ! COMMON PITFALL: Do not assume that this subroutine computes
            ! integrals for angular momenta lower than L.
            ! --------------------------------------------------------------
            ! 1. Moreno-Flores, R., Alvarez-Mendez, R., Vela, A., and
            !    Koster, A.M., Half-Numerical Evaluation of Pseudopotential
            !    Integrals, J. Comput. Chem. 27, 1009 (2006);
            !    doi: 10.1002/jcc.20410
            !
            real(F64), dimension(:), intent(out) :: w
            integer, intent(in)                  :: l
            integer, intent(in)                  :: lx0
            integer, intent(in)                  :: ly0
            integer, intent(in)                  :: lz0

            integer :: lx, ly, lz
            integer :: x, y, z
            integer :: l0
            integer :: idx
            real(F64), parameter :: fourpi = FOUR * PI
            real(F64) :: a, b, c, d
            
            l0 = lx0 + ly0 + lz0
            idx = 0
            d = dblfact(l+l0+1)
            do lx = 0, l
                  x = lx + lx0
                  a = fourpi * dblfact(x-1)
                  do ly = 0, l - lx
                        lz = l - lx - ly
                        y = ly + ly0
                        z = lz + lz0
                        idx = idx + 1
                        if ((modulo(x, 2) .ne. 0) .or. &
                              (modulo(y, 2) .ne. 0) .or. &
                              (modulo(z, 2) .ne. 0)) then
                              w(idx) = ZERO
                        else
                              b = dblfact(y-1)
                              c = dblfact(z-1)
                              w(idx) = a * b * c / d
                        end if
                  end do
            end do
      end subroutine xyzint


      pure function nslm(l, m)
            !
            ! Compute the normalization factor Nlm of the solid harmonic Slm
            ! in the basis of x**lx y**ly z**lz functions.
            !
            ! Nlm in Racah's normalization is defined in Eq. 9.1.10 of Ref. 1,
            ! but we multiply it by  Sqrt(2l+1/4Pi) to change the normalization
            ! integral to unity.
            !
            ! 1. Helgaker, T., Jorgensen, P., Olsen, J., Molecular
            !    Electronic-Structure Theory, Wiley & Sons Chichester
            !    2000. Eqs. 6.4.47-6.4.50, p. 215
            !
            real(F64)    :: nslm
            integer, intent(in) :: l, m
            
            real(F64) :: sqra, sqrb, sqrc
            real(F64), parameter :: sq2 = sqrt(TWO)
            real(F64), parameter :: fourpi = FOUR * PI
            real(F64) :: a
            
            sqra = sqrtfact(l+abs(m))
            sqrb = sqrtfact(l-abs(m))
            if (m .eq. 0) then
                  sqrc = ONE
            else
                  sqrc = sq2
            end if
            !
            ! Real solid harmonics are normalized to 1
            ! instead of Racah's normalization in Helgaker's
            ! textbook
            !
            a = sqrt(real(2*l+1, F64)/fourpi)
            nslm = a * sqra * sqrb * sqrc / (real(2**abs(m), F64) * fact(l))
      end function nslm
      
      
      pure function vm(m)
            integer             :: vm
            integer, intent(in) :: m

            if (m .ge. 0) then
                  vm = 0
            else
                  vm = 1
            end if
      end function vm


      pure function clmtuv(l, m, t, u, v)
            !
            ! Compute the expansion coefficeint Clm(t, u, v) of the solid harmonic Slm
            ! in the basis of x**lx y**ly z**lz functions. Clm(t, u, v) is defined in
            ! Eq. 9.1.11 of Ref. 1.
            !
            ! Note that the Slm functions in Helgaker's textbook are in Racah's
            ! normalization, whereas the harmonics computed in this module are
            ! normalized to unity.
            !
            ! 1. Helgaker, T., Jorgensen, P., Olsen, J., Molecular
            !    Electronic-Structure Theory, Wiley & Sons Chichester
            !    2000. Eqs. 6.4.47-6.4.50, p. 215
            !
            real(F64) :: clmtuv
            integer, intent(in) :: l, m, t, u, v
            
            real(F64) :: a, b
            integer :: c
            
            a = binom(l, t) * binom(l - t, abs(m) + t)
            b = binom(t, u) * binom(abs(m), v)
            c = (-1)**(t + (v - vm(m)) / 2) * 4**t
            clmtuv = a * b / real(c, F64)
      end function clmtuv


      pure subroutine rshu(wu, l, m)
            ! ---------------------------------------------------------------
            ! Calculate coefficient of the X^{lx}Y^{ly}Z^{lz} expansion
            ! of real spherical orthonormal harmonics: U^{lm}_{lxlylz}
            ! coefficient in Eq. 32 and 33 in [2].
            !
            ! The implementation uses the equations given in Helgaker's
            ! textbook (Eqs. 9.1.9-12 in Ref. 1), but rescaled by the factor
            ! of Sqrt(2l+1/4Pi) to change the normalization integral to
            ! unity.
            !
            ! The complete array of coefficients, for all lx + ly + lz = l
            ! is computed in a single call.
            !
            ! Definition of real spherical harmonics:
            ! For m > 0: Sl^m = (-1)^m/Sqrt(2) Yl^m + 1/Sqrt(2)       Yl^(-m)
            ! For m < 0: Sl^m =      i/Sqrt(2) Yl^m - i(-1)^m/Sqrt(2) Yl^(-m)
            ! For m = 0: Sl^m = Yl^m
            ! Yl^m are orthonormal, complex-valued spherical harmonics.
            ! ---------------------------------------------------------------
            ! 1. Helgaker, T., Jorgensen, P., Olsen, J., Molecular
            !    Electronic-Structure Theory, Wiley & Sons Chichester
            !    2000. Eqs. 6.4.47-6.4.50, p. 215
            ! 2. Moreno-Flores, R., Alvarez-Mendez, R., Vela, A., and
            !    Koster, A.M., Half-Numerical Evaluation of Pseudopotential
            !    Integrals, J. Comput. Chem. 27, 1009 (2006);
            !    doi: 10.1002/jcc.20410
            !
            real(F64), dimension(:), intent(out) :: wu
            integer, intent(in)                  :: l
            integer, intent(in)                  :: m
            
            real(F64) :: npref
            integer :: k, t, u, v
            integer :: lx, ly, lz
            integer :: i, n
            integer :: absm

            absm = abs(m)
            npref = nslm(l, m)
            n = lxlylzpos(l, 0, 0)
            wu(1:n) = ZERO
            do t = 0, (l - absm) / 2
                  do u = 0, t
                        do k = 0, (absm - vm(m)) / 2
                              v = 2 * k + vm(m)
                              lx = 2 * t + absm - 2 * u - v
                              ly = 2 * u + v
                              lz = l - 2 * t - absm
                              i = lxlylzpos(lx, ly, lz)
                              wu(i) = wu(i) + npref * clmtuv(l, m, t, u, v)
                        end do
                  end do
            end do
      end subroutine rshu


      subroutine ClmkCoeffs(Clmk, l, kappa, Clm)
            !
            ! Compute the coefficients of Slmk=r**(l+2kappa) Slm expanded
            ! in the Cartesian basis
            !
            ! Slmk = r**(l+2*kappa) Slm =
            ! Sum(u+v+w=l+2*kappa) Clmk(lxlylzpos(u,v,w)) x**u y**v z**w
            !
            ! Definition of real spherical harmonics:
            ! For m > 0: Slm = (-1)**m/Sqrt(2) Ylm + 1/Sqrt(2)       Yl(-m)
            ! For m < 0: Slm =       i/Sqrt(2) Ylm - i(-1)^m/Sqrt(2) Yl(-m)
            ! For m = 0: Slm = Ylm
            ! Yl^m are orthonormal, complex-valued spherical harmonics.
            !
            ! Definition of the expansion coefficients Clm:
            !
            ! r**l Slm = Sum(u+v+w=l) Clm(lxlylzpos(u,v,w)) x**u y**v z**w
            !
            real(F64), dimension(:), intent(out) :: Clmk
            integer, intent(in)                  :: l
            integer, intent(in)                  :: kappa
            real(F64), dimension(:), intent(in)  :: Clm

            integer :: u, v, w, mx, my, mz, lx, ly, lz, uvwidx

            Clmk = ZERO
            do lx = 0, l
                  do ly = 0, l - lx
                        do mx = 0, kappa
                              do my = 0, kappa - mx
                                    lz = l - lx - ly
                                    mz = kappa - mx - my
                                    u = lx + 2 * mx
                                    v = ly + 2 * my
                                    w = lz + 2 * mz
                                    uvwidx = lxlylzpos(u, v, w)
                                    Clmk(uvwidx) = Clmk(uvwidx) + fact(kappa)/(fact(mx)*fact(my)*fact(mz)) &
                                          * Clm(lxlylzpos(lx,ly,lz))
                              end do
                        end do
                  end do
            end do
      end subroutine ClmkCoeffs


      function rshv(l, m, lx, ly, lz, ulm, work)
            ! --------------------------------------------------------------
            ! Calculate V(lxlylz,lm) integral, Eq. 36 in [1]. Note that
            ! the RSHV matrix is the inverse of RSHU:
            ! Sum(0<=j<=l,-j<=m<=j) V(lxlylz,lm) U(mxmymz,lm) =
            ! delta(lx,mx) * delta(ly,my) * delta(lz,mz)
            ! --------------------------------------------------------------
            ! 1. Moreno-Flores, R., Alvarez-Mendez, R., Vela, A., and
            !    Koster, A.M., Half-Numerical Evaluation of Pseudopotential
            !    Integrals, J. Comput. Chem. 27, 1009 (2006)
            ! --------------------------------------------------------------
            real(F64)    :: rshv
            integer, intent(in) :: l
            integer, intent(in) :: m
            integer, intent(in) :: lx
            integer, intent(in) :: ly
            integer, intent(in) :: lz
            real(F64), dimension(:), intent(in) :: ulm
            real(F64), dimension(:), intent(out) :: work
            
            integer :: n
            integer :: k0, k
            !
            ! Generate Int_Omega X Y Z d Omega integrals
            !
            call xyzint(work, l, lx, ly, lz)
            n = numxyz(l)
            k0 = ulmpos(l, m)
            rshv = ZERO
            do k = 1, n
                  rshv = rshv + ulm(k0+k-1) * work(k)
            end do
      end function rshv

      
      pure subroutine rsh(wslm, lmax, xn, yn, zn)
            ! ---------------------------------------------------------------
            ! Compute real spherical orthonormal harmonics in Cartesian
            ! representation, S_{LM}(X_N, Y_N, Z_N), where
            ! X_N^2 + Y_N^2 + Z_N^2 = 1.
            !
            ! S_{LM} functions are computed for every L <= LMAX.
            !
            ! The functions computed with this subroutine are normalized to 1:
            ! Int d Omega S_{l'm'} S_{lm} = delta_{ll'}delta{mm'}
            !
            ! Definition of real spherical harmonics:
            ! For m > 0: Sl^m = (-1)^m/Sqrt(2) Yl^m + 1/Sqrt(2)       Yl^(-m)
            ! For m < 0: Sl^m =      i/Sqrt(2) Yl^m - i(-1)^m/Sqrt(2) Yl^(-m)
            ! For m = 0: Sl^m = Yl^m
            ! Yl^m are orthonormal, complex-valued spherical harmonics.
            ! ---------------------------------------------------------------
            ! 1. Helgaker, T., Jorgensen, P., Olsen, J., Molecular
            !    Electronic-Structure Theory, Wiley & Sons Chichester
            !    2000. Eqs. 6.4.47-6.4.50, p. 215
            !
            !
            real(F64), dimension(:), intent(out) :: wslm
            integer, intent(in)                  :: lmax
            real(F64), intent(in)                :: xn
            real(F64), intent(in)                :: yn
            real(F64), intent(in)                :: zn

            real(F64) :: a, b, c, d, e
            real(F64), parameter :: fourpi = FOUR * PI
            real(F64) :: norm
            integer :: l, m
            !
            ! Generate orthogonal solid harmonics from recurrence relations
            ! 6.4.70-73, page 218 in [1]
            !
            wslm(slmpos(0, 0)) = ONE
            if (lmax >= 1) then
                  wslm(slmpos(1, 1)) = xn * wslm(slmpos(0, 0))
                  wslm(slmpos(1, -1)) = yn * wslm(slmpos(0, 0))
                  wslm(slmpos(1, 0)) = zn * wslm(slmpos(0, 0))
            end if

            do l = 1, lmax - 1
                  a = real(2 * l + 1, F64)
                  b = real(2 * l + 2, F64)
                  e = sqrt(a / b)
                  wslm(slmpos(l+1, l+1)) = e * (xn * wslm(slmpos(l, l)) - yn * wslm(slmpos(l, -l)))
                  wslm(slmpos(l+1, -l-1)) = e * (yn * wslm(slmpos(l, l)) + xn * wslm(slmpos(l, -l)))
                  do m = -l, l
                        d = sqrt(real((l+m+1)*(l-m+1), F64))
                        if ((m > -l) .and. (m < l)) then
                              c = sqrt(real((l+m)*(l-m), F64))
                              wslm(slmpos(l+1, m)) = (a * zn * wslm(slmpos(l, m)) - c * wslm(slmpos(l-1, m))) / d
                        else
                              wslm(slmpos(l+1, m)) = (a * zn * wslm(slmpos(l, m))) / d
                        end if
                  end do
            end do
            !
            ! Make all functions biorhonormal:
            ! Int d Omega S_{l'm'} S_{lm} = delta_{ll'}delta{mm'}
            !
            do l = 0, lmax
                  norm = sqrt(real(2*l+1, F64)/fourpi)
                  do m = -l, l
                        wslm(slmpos(l, m)) = norm * wslm(slmpos(l, m))
                  end do
            end do
      end subroutine rsh


      pure subroutine AngularMomentum(lvec, l, mm, m)
            !
            ! Compute a matrix element of the angular momentum operator (vector)
            ! in the basis of real spherical harmonics:
            !
            ! 1/i <Sl^(mm)|L|Sl^m> = 1/i (<Sl^(mm)|Lx|Sl^m>, <Sl^(mm)|Ly|Sl^m>, <Sl^(mm)|Lz|Sl^m>)
            !
            ! The matrix of L in the basis of real spherical harmonics is purely
            ! imaginary. The output of this subroutine corresponds to 1/i L to avoid
            ! complex numbers.
            !
            ! Definition of real spherical harmonics:
            ! For m > 0: Sl^m = (-1)^m/Sqrt(2) Yl^m + 1/Sqrt(2)       Yl^(-m)
            ! For m < 0: Sl^m =      i/Sqrt(2) Yl^m - i(-1)^m/Sqrt(2) Yl^(-m)
            ! For m = 0: Sl^m = Yl^m
            ! Yl^m are orthonormal, complex-valued spherical harmonics.
            !
            real(F64), dimension(3), intent(out) :: lvec
            integer, intent(in)                  :: l
            integer, intent(in)                  :: mm
            integer, intent(in)                  :: m
            
            real(F64) :: c, x, y, z
            
            x = ZERO
            y = ZERO
            z = ZERO
            if (mm > 0 .and. m > 0) then
                  if (mm == m + 1) then
                        c = (1 + l - mm) * (l + mm)
                        y = ONE/TWO*sqrt(c)
                  else if (mm == m - 1) then
                        c = (l - mm) * (1 + l + mm)
                        y = -ONE/TWO*sqrt(c)
                  end if
            else if (mm < 0 .and. m < 0) then
                  if (mm == m + 1) then
                        c = (1 + l - mm) * (l + mm)
                        y = -ONE/TWO*sqrt(c)
                  else if (mm == m - 1) then
                        c = (l - mm) * (1 + l + mm)
                        y = ONE/TWO*sqrt(c)
                  end if
            else if (mm == 1 .and. m == 0) then
                  c = l * (l + 1)
                  y = sqrt(c)/sqrt(TWO)
            else if (mm == 0 .and. m == 1) then
                  c = l * (l + 1)
                  y = -sqrt(c)/sqrt(TWO)
            else if (mm == -1 .and. m == 0) then
                  c = l * (l + 1)
                  x = -sqrt(c)/sqrt(TWO)
            else if (mm == 0 .and. m == -1) then
                  c = l * (l + 1)
                  x = sqrt(c)/sqrt(TWO)
            else if (mm > 0 .and. m < 0) then
                  if (1 + m + mm == 0) then
                        c = (l - mm) * (1 + l + mm)
                        x = ONE/TWO*sqrt(c)
                  else if (m + mm == 1) then
                        c = (1 + l - mm) * (l + mm)
                        x = ONE/TWO*sqrt(c)
                  else if (m + mm == 0) then
                        z = -mm
                  end if
            else if (mm < 0 .and. m > 0) then
                  if (1 + m + mm == 0) then
                        c = (l - mm) * (1 + l + mm)
                        x = -ONE/TWO*sqrt(c)
                  else if (m + mm == 1) then
                        c = (1 + l - mm) * (l + mm)
                        x = -ONE/TWO*sqrt(c)
                  else if (m + mm == 0) then
                        z = -mm
                  end if
            end if
            lvec(1) = x
            lvec(2) = y
            lvec(3) = z
      end subroutine AngularMomentum
      

      pure function lxlylzpos(lx, ly, lz)
            !
            ! Index of the Cartesian polynomial x**lx y**ly z**lz, defined
            ! by the loop
            !
            ! idx = 1
            ! do lx = 0, l
            ! do ly = 0, l - lx
            ! lz = l - lx - ly
            ! lxlylzpos(lx, ly, lz) <- idx
            ! idx = idx + 1
            ! end do
            ! end do
            ! end do
            !
            integer :: lxlylzpos
            integer, intent(in) :: lx, ly, lz
            
            integer :: l

            l = lx + ly + lz
            lxlylzpos = ((2 * l - lx + 3) * lx) / 2 + ly + 1
      end function lxlylzpos


      pure function ulmpos(l, m)
            !
            ! Position of a first element of Ulm(LX,LY,LZ)
            ! (Ulm(0,0,L)) in an array generated by RSHU
            ! subroutine
            !
            integer             :: ulmpos
            integer, intent(in) :: l
            integer, intent(in) :: m

            ulmpos = numulm(l-1) + numxyz(l) * (m + l) + 1
      end function ulmpos


      pure function plmpos(l, m)
            integer             :: plmpos
            integer, intent(in) :: l
            integer, intent(in) :: m

            integer :: n
            
            n = numplm(l-1)
            plmpos = n + m + 1
      end function plmpos


      pure function slmpos(l, m)
            integer :: slmpos
            integer, intent(in) :: l
            integer, intent(in) :: m

            integer :: n
            
            n = numslm(l-1)
            slmpos = n + m + l + 1
      end function slmpos


      pure function numxyz(l)
            !
            ! Number of three-index polynomials X^{LX} Y^{LY} Z^{LZ}
            ! for which LX+LY+LZ=L. The result is equivalent to
            ! LXLYLZPOS(L, 0, 0).
            !
            integer :: numxyz
            integer, intent(in) :: l
            
            numxyz = ((l + 1) * (l + 2)) / 2
      end function numxyz

      
      pure function numplm(lmax)
            !
            ! Number of associated Legendre polynomials, P_l^{|m|}, for L <= LMAX
            ! LMAX >= -1
            ! NUMPLM(-1) = 0
            !
            integer :: numplm
            integer, intent(in) :: lmax

            numplm = ((lmax + 1) * (lmax + 2)) / 2
      end function numplm


      pure function numslm(lmax)
            !
            ! Number of real spherical harmonics, S_{lm}, for L <= LMAX
            ! LMAX >= -1
            ! NUMSLM(-1) = 0
            !
            integer :: numslm
            integer, intent(in) :: lmax

            numslm = 2 * numplm(lmax) - (lmax + 1)
      end function numslm


      pure function numulm(lmax)
            !
            ! Number of U^{LM}_{LX,LY,LZ} matrix elements for
            ! L = 0, 1, ..., LMAX
            !
            integer             :: numulm
            integer, intent(in) :: lmax

            numulm = ((1+lmax)*(2+lmax)*(3+lmax)*(2+3*lmax)) / 12
      end function numulm


      pure function lmfrac12(l, m)
            !
            ! SQRT((L-ABS(M))!/(L+ABS(M))!)
            !
            real(F64)    :: lmfrac12
            integer, intent(in) :: l
            integer, intent(in) :: m

            lmfrac12 = sqrtfact(l-abs(m)) / sqrtfact(l+abs(m))
      end function lmfrac12


      pure function sqrtfact(n)
            !
            ! SQRT(N!)
            !
            real(F64)    :: sqrtfact
            integer, intent(in) :: n

            real(F64), dimension(0:30), parameter :: sf = (/ &
                  1.d+0, &                 ! 0
                  1.d+0, &                 ! 1
                  1.414213562373095d+0, &  ! 2
                  2.449489742783178d+0, &  ! 3
                  4.898979485566356d+0, &  ! 4
                  10.95445115010332d+0, &  ! 5
                  26.83281572999748d+0, &  ! 6
                  70.9929573971954d+0,  &  ! 7
                  200.7984063681781d+0, &  ! 8
                  602.3952191045344d+0, &  ! 9
                  1904.940943966505d+0, &  ! 10
                  6317.974358922328d+0, &  ! 11
                  21886.10518114176d+0, &  ! 12
                  78911.47445080469d+0, &  ! 13
                  295259.7012800765d+0, &  ! 14
                  1.143535905863913d+6, &  ! 15
                  4.574143623455652d+6, &  ! 16
                  1.885967730625315d+7, &  ! 17
                  8.001483428544985d+7, &  ! 18
                  3.487765766344294d+8, &  ! 19
                  1.559776268628498d+9, &  ! 20
                  7.147792818185865d+9, &  ! 21
                  3.352612008237171d+10, & ! 22
                  1.607856235454059d+11, & ! 23
                  7.876854713229384d+11, & ! 24
                  3.938427356614691d+12, & ! 25
                  2.008211794424596d+13, & ! 26
                  1.04349745809074d+14, &  ! 27
                  5.521669535672285d+14, & ! 28
                  2.973510046012911d+15, & ! 29
                  1.628658527169496d+16/)  ! 30
                  
                  sqrtfact = sf(n)
      end function sqrtfact


      pure function fact(n)
            !
            ! Factorial function: N!
            !
            real(F64) :: fact
            integer, intent(in) :: n

            real(F64), dimension(0:50), parameter :: f = [ &
                  1.0_F64, & ! 0
                  1.0_F64, & ! 1
                  2.0_F64, & ! 2
                  6.0_F64, & ! 3
                  24.0_F64, & ! 4
                  120.0_F64, & ! 5
                  720.0_F64, & ! 6
                  5040.0_F64, & ! 7
                  40320.0_F64, & ! 8
                  362880.0_F64, & ! 9
                  3.6288E+6_F64, & ! 10
                  3.99168E+7_F64,  & ! 11
                  4.790016E+8_F64, & ! 12
                  6.2270208E+9_F64, & ! 13
                  8.71782912E+10_F64, & ! 14
                  1.307674368E+12_F64, & ! 15
                  2.0922789888E+13_F64, & ! 16
                  3.55687428096E+14_F64, & ! 17
                  6.402373705728E+15_F64, & ! 18
                  1.21645100408832E+17_F64, & ! 19
                  2.43290200817664E+18_F64, & ! 20
                  5.109094217170944e19_F64, & ! 21
                  1.12400072777760768e21_F64, & ! 22
                  2.585201673888497664e22_F64, & ! 23
                  6.2044840173323943936e23_F64, & ! 24
                  1.5511210043330985984e25_F64, & ! 25
                  4.03291461126605635584e26_F64, & ! 26
                  1.0888869450418352160768e28_F64, & ! 27
                  3.04888344611713860501504e29_F64, & ! 28
                  8.841761993739701954543616e30_F64, & ! 29
                  2.6525285981219105863630848e32_F64, & ! 30
                  8.22283865417792281772556288e33_F64, & ! 31
                  2.6313083693369353016721801216e35_F64, & ! 32
                  8.68331761881188649551819440128e36_F64, & ! 33
                  2.9523279903960414084761860964352e38_F64, & ! 34
                  1.03331479663861449296666513375232e40_F64, & ! 35
                  3.719933267899012174679994481508352e41_F64, & ! 36
                  1.37637530912263450463159795815809024e43_F64, & ! 37
                  5.230226174666011117600072241000742912e44_F64, & ! 38
                  2.03978820811974433586402817399028973568e46_F64, & ! 39
                  8.15915283247897734345611269596115894272e47_F64, & ! 40
                  3.3452526613163807108170062053440751665152e49_F64, & ! 41
                  1.405006117752879898543142606244511569936384e51_F64, & ! 42
                  6.0415263063373835637355132068513997507264512e52_F64, & ! 43
                  2.658271574788448768043625811014615890319638528e54_F64, & ! 44
                  1.1962222086548019456196316149565771506438373376e56_F64, & ! 45
                  5.50262215981208894985030542880025489296165175296e57_F64, & ! 46
                  2.5862324151116818064296435515361197996919763238912e59_F64, & ! 47
                  1.241391559253607267086228904737337503852148635467776e61_F64, & ! 48
                  6.0828186403426756087225216332129537688755283137921024e62_F64, & ! 49
                  3.0414093201713378043612608166064768844377641568960512e64_F64 & ! 50
                  ]
            
            fact = f(n)
      end function fact


      pure function binom(n, k)
            !
            ! Compute the binomial coefficient C(n,k) = n!/(k!(n-k)!)
            ! The computation is done using the recurrence formula
            ! C(n, k+1) = [(n-k) * C(n, k)] / (k+1)
            !
            real(F64) :: binom
            integer, intent(in) :: n
            integer, intent(in) :: k

            integer :: i
            real(F64) :: v
            real(F64) :: di, dn
            integer :: k0

            dn = real(n, F64)
            k0 = min(k, n-k)
            v = ONE
            do i = 1, k0
                  di = real(i, F64)
                  v = ((dn - di + ONE) * v) / di
            end do
            binom = v
      end function binom


      pure function dblfact(n)
            !
            ! N!! = n * (n - 2) ... 5 * 3 * 1      n > 0 odd
            !       n * (n - 2) ... 6 * 4 * 2      n > 0 even
            !       1                              n = -1, 0.
            ! 
            !
            real(F64) :: dblfact
            integer, intent(in) :: n

            real(F64), dimension(-1:30), parameter :: f = (/ &
                  1.0_F64, &                   ! -1
                  1.0_F64, &                   ! 0
                  1.0_F64, &                   ! 1
                  2.0_F64, &                   ! 2
                  3.0_F64, &                   ! 3
                  8.0_F64, &                   ! 4
                  15.0_F64, &                  ! 5
                  48.0_F64, &                  ! 6
                  105.0_F64, &                 ! 7
                  384.0_F64, &                 ! 8
                  945.0_F64, &                 ! 9
                  3840.0_F64, &                ! 10
                  10395.0_F64, &               ! 11
                  46080.0_F64, &               ! 12
                  135135.0_F64, &              ! 13
                  645120.0_F64, &              ! 14
                  2.027025E+6_F64, &           ! 15
                  1.032192E+7_F64, &           ! 16
                  3.4459425E+7_F64, &          ! 17
                  1.8579456E+8_F64, &          ! 18
                  6.54729075E+8_F64, &         ! 19
                  3.7158912E+9_F64, &          ! 20
                  1.3749310575E+10_F64, &      ! 21
                  8.17496064E+10_F64, &        ! 22
                  3.16234143225E+11_F64, &     ! 23
                  1.9619905536E+12_F64, &      ! 24
                  7.905853580625E+12_F64, &    ! 25
                  5.10117543936E+13_F64, &     ! 26
                  2.13458046676875E+14_F64, &  ! 27
                  1.4283291230208E+15_F64, &   ! 28
                  6.190283353629375E+15_F64, & ! 29
                  4.2849873690624E+16_F64 &    ! 30
            /)

            dblfact = f(n)
      end function dblfact


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

      
      pure function dphase(k)
            !
            ! (-1)^K
            !
            real(F64) :: dphase
            integer, intent(in) :: k

            dphase = real(iphase(k), F64)
      end function dphase
end module spherh
