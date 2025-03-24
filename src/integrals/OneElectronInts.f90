module OneElectronInts
      use arithmetic
      use auto2e
      use auto2e_Hermite
      use auto2e_BraTransform
      use boys
      use basis_sets
      use sys_definitions
      use display
      use real_linalg
      use sphergto

      implicit none

      type TAuto2eRtuv
            procedure(auto2e_Rtuv_1), pointer, nopass :: ptr
      end type TAuto2eRtuv

      type TAuto2eBoys
            procedure(f0), pointer, nopass :: ptr
      end type TAuto2eBoys

contains

      pure function ints1e_hposition(l, j, t)
            !
            ! Calculate the position of the E^{ij}_t element
            ! (i + j = l) in a linearized array
            !
            integer :: ints1e_hposition

            integer, intent(in) :: l, j, t

            ints1e_hposition = j * (l + 1) + t + 1
      end function ints1e_hposition


      pure subroutine ints1e_eijmatrix_new_level(oldlev, newlev, lold, p, xpa, xpb)
            ! ---------------------------------------------------------------
            ! Compute the coefficients of the transformation from a product
            ! of Gassian-type functions to Hermite-Gaussian basis.
            ! The E^{ij}_t coefficients are defined in Helgaker's textbook,
            ! see Eqs 9.5.6 and 9.5.7 in [1]. This subroutine performs a 
            ! single iteration of the vertical recurrence relation.
            ! ---------------------------------------------------------------
            ! 1. Helgaker, T., Jorgensen, P., Olsen, J., Molecular
            !    Electronic-Structure Theory, 2000.
            !
            real(F64), dimension(:), intent(in)  :: oldlev
            real(F64), dimension(:), intent(out) :: newlev
            integer, intent(in)                  :: lold
            real(F64), intent(in)                :: p
            real(F64), intent(in)                :: xpa
            real(F64), intent(in)                :: xpb

            integer :: lnew
            integer :: k, t, u, pos1, pos2
            real(F64) :: e1, e2
            real(F64) :: twop
            real(F64) :: dt, di, dk

            lnew = lold + 1
            u = 1
            twop = two * p

            !
            ! E^{lnew, 0}_t, t = 0, ... , lnew
            ! ---
            ! E^{lnew, 0}_0
            !
            pos1 = 1
            pos2 = 2

            e1 = oldlev(pos1)
            e2 = oldlev(pos2)
            newlev(u) = xpa * e1 + e2
            u = u + 1
            !
            ! E^{lnew, 0}_t, t = 1, ... , lnew
            !
            di = real(lnew, F64)
            pos1 = 1
            do t = 1, lnew
                  dt = real(t, F64)
                  e1 = oldlev(pos1)
                  newlev(u) = one / (twop * dt) * di * e1
                  u = u + 1
                  pos1 = pos1 + 1
            end do
            !
            ! E^{lnew - k, k}_t, k = 1, ... , lnew - 1
            !
            do k = 1, lnew - 1
                  !
                  ! E^{lnew - k, k}_0
                  !
                  di = real(lnew-k, F64)
                  dk = real(k, F64)
                  pos1 = ints1e_hposition(lold, k, 0)
                  pos2 = pos1 + 1
                  e1 = oldlev(pos1)
                  e2 = oldlev(pos2)
                  newlev(u) = xpa * e1 + e2
                  u = u + 1
                  !
                  ! E^{lnew - k, k}_t, t = 1, ... , lnew
                  !
                  pos2 = ints1e_hposition(lold, k - 1, 0)
                  do t = 1, lnew
                        dt = real(t, F64)
                        e1 = oldlev(pos1)
                        e2 = oldlev(pos2)
                        newlev(u) = one / (twop * dt) * (di * e1 + dk * e2)
                        u = u + 1
                        pos1 = pos1 + 1
                        pos2 = pos2 + 1
                  end do
            end do
            !
            ! E^{0, lnew}_0
            !
            pos1 = ints1e_hposition(lold, lold, 0)
            pos2 = pos1 + 1
            e1 = oldlev(pos1)
            e2 = oldlev(pos2)
            newlev(u) = xpb * e1 + e2
            u = u + 1
            !
            ! E^{0, lnew}_t, t = 1, ... , lnew
            !
            dk = real(lnew, F64)
            pos1 = ints1e_hposition(lold, lnew - 1, 0)
            do t = 1, lnew
                  dt = real(t, F64)
                  e1 = oldlev(pos1)
                  newlev(u) = one / (twop * dt) * dk * e1
                  u = u + 1
                  pos1 = pos1 + 1
            end do
      end subroutine ints1e_eijmatrix_new_level


      pure subroutine ints1e_eijmatrix(l, seed, p, xpa, xpb, eijm)
            ! ---------------------------------------------------------------
            ! Compute the coefficients of the transformation from a product
            ! of Gassian-type functions to Hermite-Gaussian basis.
            ! The E^{ij}_t coefficients are defined in Helgaker's textbook,
            ! see Eqs 9.5.6 and 9.5.7 in [1].
            ! ---------------------------------------------------------------
            ! 1. Helgaker, T., Jorgensen, P., Olsen, J., Molecular
            !    Electronic-Structure Theory, 2000.
            !
            integer, intent(in)                  :: l
            real(F64), intent(in)                :: seed
            real(F64), intent(in)                :: p
            real(F64), intent(in)                :: xpa
            real(F64), intent(in)                :: xpb
            real(F64), dimension(:), intent(out) :: eijm

            integer :: new_level, old_level, k
            real(F64), dimension(2) :: init_array

            eijm(1) = seed

            if (l > 0) then
                  init_array(1) = seed
                  init_array(2) = ZERO
                  new_level = 2
                  call ints1e_eijmatrix_new_level(init_array, eijm(new_level:), 0, p, xpa, xpb)
                  old_level = new_level
                  new_level = new_level + 4
                  do k = 2, l
                        call ints1e_eijmatrix_new_level(eijm(old_level:), eijm(new_level:), k - 1, p, xpa, xpb)
                        old_level = new_level
                        new_level = new_level + (k + 1)**2
                  end do
            end if
      end subroutine ints1e_eijmatrix
      

      pure function ints1e_ELoc(lA, lB, t)
            integer             :: ints1e_ELoc
            integer, intent(in) :: lA
            integer, intent(in) :: lB
            integer, intent(in) :: t

            integer :: Lab

            Lab = lA + lB
            ints1e_ELoc = (Lab*(Lab+1)*(2*Lab+1))/6 + lB*(Lab+1) + t + 1
      end function ints1e_ELoc
      

      pure subroutine ints1e_decode_pq(pq, n, p, q)
            !
            ! Decode a lower-triangle compound index into individual
            ! indices:
            ! PQ -> (P, Q)
            ! Assumptions:
            ! 0) P = 1, 2, ..., N,
            !    Q = 1, 2, ..., N,
            ! 1) P >= Q (diagonal indices admissible)
            !
            ! An example of how this algorithm traverses an N=3 triangle:
            !
            !      Q
            !    1
            ! P  2 5
            !    3 6 4
            !
            integer, intent(in)  :: pq
            integer, intent(in)  :: n
            integer, intent(out) :: p
            integer, intent(out) :: q

            integer :: q_base
            integer :: v
            integer :: interval1
            integer :: in1, in2
            !
            ! pq = (q_base - 1) * (n + 1) + v
            !
            q_base = (pq - 1) / (n + 1) + 1
            v = pq - (n + 1) * (q_base - 1)
            !
            ! Decide if v is in interval_1 or interval_2:
            ! in1 == 1 and in2 == 0 if v <= INTERVAL1
            ! in1 == 0 and in2 == 1 if v > INTERVAL1
            !
            interval1 = n - q_base + 1
            in2 = v / (interval1 + 1)
            !
            ! 1 -> 0, 0 -> 1
            !
            in1 = ieor((in2), 1)

            p = in1 * (q_base + v - 1) + in2 * (v - interval1 + n - q_base)          
            q = in1 * q_base + in2 * interval1
      end subroutine ints1e_decode_pq


      subroutine ints1e_overlap_AB(S, Na, Nb, Ra, La, CntrA, ExpA, NormA, NprimA, Rb, Lb, &
            CntrB, ExpB, NormB, NprimB, ax, ay, az, bx, by, bz, BinomTable)
            !
            ! Compute overlap integrals for a pair of Cartesian Gaussian shells
            ! of contracted atomic orbitals. This subroutine will work for any angular momenta
            ! La and Lb provided that the binomial coefficients C(max(La,Lb),k) are stored
            ! in BinomTable.
            !
            real(F64), dimension(Na, Nb), intent(out) :: S
            integer, intent(in)                       :: Na
            integer, intent(in)                       :: Nb
            integer, intent(in)                       :: La
            real(F64), dimension(3), intent(in)       :: Ra
            real(F64), dimension(*), intent(in)       :: CntrA, ExpA, NormA
            integer, intent(in)                       :: NprimA
            real(F64), dimension(3), intent(in)       :: Rb
            integer, intent(in)                       :: Lb
            real(F64), dimension(*), intent(in)       :: CntrB, ExpB, NormB
            integer, intent(in)                       :: NPrimB
            integer, dimension(:), intent(in)         :: ax, ay, az
            integer, dimension(:), intent(in)         :: bx, by, bz
            real(F64), dimension(0:, 0:), intent(in)  :: BinomTable

            real(F64), dimension(3) :: Rp, Rpa, Rpb, Rab
            real(F64), dimension(3, 0:La) :: PA
            real(F64), dimension(3, 0:Lb) :: PB
            real(F64) :: GaussFactor, Rab2
            real(F64) :: AlphaA, AlphaB, AlphaP
            real(F64), dimension(3) :: W
            real(F64), dimension(0:(La+Lb)/2) :: OneDInts
            real(F64), dimension(0:La, 0:Lb) :: XInts, YInts, ZInts
            integer :: t, u
            integer :: k, l
            integer :: a, b

            OneDInts = ZERO
            S = ZERO
            Rab = Ra - Rb
            Rab2 = dot_product(Rab, Rab)
            do l = 1, NprimB
                  do k = 1, NprimA
                        AlphaA = ExpA(k)
                        AlphaB = ExpB(l)
                        AlphaP = ExpA(k) + ExpB(l)
                        GaussFactor = exp(-AlphaA*AlphaB*Rab2/AlphaP)
                        Rp = (AlphaA*Ra + AlphaB*Rb) / AlphaP
                        Rpa = Rp - Ra
                        Rpb = Rp - Rb
                        do a = 0, La
                              PA(:, a) = Rpa**a
                        end do
                        do b = 0, Lb
                              PB(:, b) = Rpb**b
                        end do
                        do t = 0, (La+Lb)/2
                              !
                              ! Integrate(-Inf,+Inf) x**(2*t) * exp(-AlphaP * x**2) dx
                              !
                              OneDInts(t) = dblfact(2*t-1) / (TWO * AlphaP)**t * Sqrt(PI/AlphaP)
                        end do
                        do b = 0, Lb
                              do a = 0, La
                                    W = ZERO
                                    do u = 0, b
                                          do t = modulo(a+b-u,2), a, 2
                                                W = W + BinomTable(a, t) * PA(:, t) &
                                                      * BinomTable(b, u) * PB(:, u) &
                                                      * OneDInts((a+b-t-u)/2)
                                          end do
                                    end do
                                    XInts(a, b) = W(1)
                                    YInts(a, b) = W(2)
                                    ZInts(a, b) = W(3)
                              end do
                        end do
                        XInts = GaussFactor * CntrA(k) * CntrB(l) * XInts
                        do b = 1, Nb
                              do a = 1, Na
                                    S(a, b) = S(a, b) + &
                                          NormA(a) * NormB(b) &
                                          * XInts(ax(a), bx(b)) &
                                          * YInts(ay(a), by(b)) &
                                          * ZInts(az(a), bz(b))
                              end do
                        end do
                  end do
            end do
      end subroutine ints1e_overlap_AB
      

      subroutine ints1e_OverlapMatrix(S, AOBasis)
            !
            ! Overlap matrix in a Cartesian GTO basis.
            ! Only lower triangle is referenced.
            !
            real(F64), dimension(:, :), intent(out) :: S
            type(TAOBasis), intent(in)              :: AOBasis

            integer :: ShellA, ShellB, ShellParamsA, ShellParamsB, ShellAB
            integer :: La, Na, a0, a1
            integer :: Lb, Nb, b0, b1
            integer :: L, t
            integer, parameter :: MaxNFunc = ((AUTO2E_MAXL + 1) * (AUTO2E_MAXL + 2)) / 2
            real(F64), dimension(MaxNFunc**2) :: Sab
            real(F64), dimension(:, :), allocatable :: BinomTable

            S = ZERO
            associate ( &
                  ShellParamsIdx => AOBasis%ShellParamsIdx, &
                  ShellMomentum => AOBasis%ShellMomentum, &
                  NAngFunc => AOBasis%NAngFuncCart, &
                  ShellLoc => AOBasis%ShellLocCart, &
                  NShells => AOBasis%NShells, &
                  CartPolyX => AOBasis%CartPolyX, &
                  CartPolyY => AOBasis%CartPolyY, &
                  CartPolyZ => AOBasis%CartPolyZ, &
                  LmaxGTO => AOBasis%LmaxGTO, &
                  ShellCenters => AOBasis%ShellCenters, &
                  AtomCoords => AOBasis%AtomCoords, &
                  CntrCoeffs => AOBasis%CntrCoeffs, &
                  Exponents => AOBasis%Exponents, &
                  NPrimitives => AOBasis%NPrimitives, &
                  NormFactors => AOBasis%NormFactorsCart &
                  )
                  allocate(BinomTable(0:LmaxGTO, 0:LmaxGTO))
                  do L = 0, LmaxGTO
                        do t = 0, L
                              BinomTable(L, t) = binom(L, t)
                        end do
                  end do
                  !$omp parallel &
                  !$omp private(ShellA, ShellParamsA, La, Na, a0, a1) &
                  !$omp private(ShellB, ShellParamsB, Lb, Nb, b0, b1) &
                  !$omp private(Sab, ShellAB) &
                  !$omp shared(S, BinomTable) &
                  !$omp default(shared)
                  !$omp do schedule(dynamic)
                  do ShellAB = 1, (NShells * (NShells + 1)) / 2
                        call ints1e_decode_pq(ShellAB, NShells, ShellA, ShellB)
                        ShellParamsA = ShellParamsIdx(ShellA)
                        La = ShellMomentum(ShellParamsA)
                        Na = NAngFunc(ShellParamsA)
                        a0 = ShellLoc(ShellA)
                        a1 = ShellLoc(ShellA) + Na - 1


                        ShellParamsB = ShellParamsIdx(ShellB)
                        Lb = ShellMomentum(ShellParamsB)
                        Nb = NAngFunc(ShellParamsB)
                        b0 = ShellLoc(ShellB)
                        b1 = ShellLoc(ShellB) + Nb - 1

                        call ints1e_overlap_AB(Sab, Na, Nb, &
                              !
                              AtomCoords(:, ShellCenters(ShellA)), &
                              ShellMomentum(ShellParamsA), &
                              CntrCoeffs(:, ShellParamsA), &
                              Exponents(:, ShellParamsA), &
                              NormFactors(:, ShellParamsA), &
                              NPrimitives(ShellParamsA), &
                              !
                              AtomCoords(:, ShellCenters(ShellB)), &
                              ShellMomentum(ShellParamsB), &
                              CntrCoeffs(:, ShellParamsB), &
                              Exponents(:, ShellParamsB), &
                              NormFactors(:, ShellParamsB), &
                              NPrimitives(ShellParamsB), &
                              !
                              CartPolyX(:, La), &
                              CartPolyY(:, La), &
                              CartPolyZ(:, La), &
                              !
                              CartPolyX(:, Lb), &
                              CartPolyY(:, Lb), &
                              CartPolyZ(:, Lb), &
                              !
                              BinomTable)

                        S(a0:a1, b0:b1) = reshape(Sab(1:Na*Nb), [Na, Nb])
                  end do
                  !$omp end do
                  !$omp end parallel
            end associate
      end subroutine ints1e_OverlapMatrix


      subroutine ints1e_Coulomb_core(Rtuv, RtuvC, fmarray, x, Qc, Lab, AlphaAB, Rpc, Auto2eRtuv, Auto2eBoys)
            real(F64), dimension(:), intent(inout)      :: Rtuv
            real(F64), dimension(:), intent(out)        :: RtuvC
            real(F64), dimension(:), intent(out)        :: fmarray
            real(F64), intent(in)                       :: x
            real(F64), intent(in)                       :: Qc
            integer, intent(in)                         :: Lab
            real(F64), intent(in)                       :: AlphaAB
            real(F64), dimension(3), intent(in)         :: Rpc
            type(TAuto2eRtuv), dimension(:), intent(in) :: Auto2eRtuv
            type(TAuto2eBoys), dimension(:), intent(in) :: Auto2eBoys

            integer :: N

            if (Lab == 0) then
                  call f0(x, fmarray)
                  RtuvC(1) = fmarray(1)
            else if (Lab <= 6) then
                  call Auto2eBoys(Lab)%ptr(x, fmarray)
                  call Auto2eRtuv(Lab)%ptr(RtuvC, fmarray, AlphaAB, Rpc(1), Rpc(2), Rpc(3))
            else
                  call fm(Lab, x, fmarray)
                  call Auto2eRtuv(Lab)%ptr(RtuvC, fmarray, AlphaAB, Rpc(1), Rpc(2), Rpc(3))
            end if
            N = ((Lab+1)*(Lab**2+5*Lab+6))/6
            Rtuv(1:N) = Rtuv(1:N) - Qc * RtuvC(1:N)
      end subroutine ints1e_Coulomb_core

      
      subroutine ints1e_HermiteTransf(Vab, R, La, Lb, Na, Nb, ExAB, EyAB, EzAB, &
            Prefactor, NormA, NormB, CartPolyX, CartPolyY, CartPolyZ)
            
            real(F64), dimension(Na, Nb), intent(inout) :: Vab
            real(F64), dimension(:), intent(in)         :: R
            integer, intent(in)                         :: La
            integer, intent(in)                         :: Lb
            integer, intent(in)                         :: Na
            integer, intent(in)                         :: Nb
            real(F64), dimension(:), intent(in)         :: ExAB
            real(F64), dimension(:), intent(in)         :: EyAB
            real(F64), dimension(:), intent(in)         :: EzAB
            real(F64), intent(in)                       :: Prefactor
            real(F64), dimension(:), intent(in)         :: NormA
            real(F64), dimension(:), intent(in)         :: NormB
            integer, dimension(:, 0:), intent(in)       :: CartPolyX
            integer, dimension(:, 0:), intent(in)       :: CartPolyY
            integer, dimension(:, 0:), intent(in)       :: CartPolyZ

            integer :: tuv
            integer :: a, b, t, u, v
            integer :: lxA, lyA, lzA
            integer :: lxB, lyB, lzB
            integer :: x0, y0, z0
            integer :: Lab, ly
            real(F64) :: F

            Lab = La + Lb
            do b = 1, Nb
                  do a = 1, Na
                        lxA = CartPolyX(a, La)
                        lyA = CartPolyY(a, La)
                        lzA = CartPolyZ(a, La)
                        lxB = CartPolyX(b, Lb)
                        lyB = CartPolyY(b, Lb)
                        lzB = CartPolyZ(b, Lb)
                        x0 = ints1e_ELoc(lxA, lxB, 0)
                        y0 = ints1e_ELoc(lyA, lyB, 0)
                        z0 = ints1e_ELoc(lzA, lzB, 0)
                        F = Prefactor * NormA(a) * NormB(b)
                        ly = lyA + lyB
                        tuv = 1
                        do v = 0, lzA+lzB
                              do u = 0, lyA+lyB
                                    do t = 0, lxA+lxB
                                          Vab(a, b) = Vab(a, b) + F * ExAB(x0+t) * EyAB(y0+u) * EzAB(z0+v) * R(tuv)
                                          tuv = tuv + 1
                                    end do
                                    tuv = tuv + (Lab-u-v) - (lxA+lxB)
                              end do
                              !
                              ! Mathematica code:
                              ! delta = Sum[Sum[1, {t,0,Lab-u-v}], {u, lyA+lyB+1, Lab-v}]
                              !
                              tuv = tuv + ((Lab-ly-v)*(1+Lab-ly-v))/2
                        end do
                  end do
            end do
      end subroutine ints1e_HermiteTransf


      subroutine ints1e_Coulomb(V, AOBasis, System)
            !
            ! Coulomb matrix in a Cartesian GTO basis
            !
            ! 1. Nuclear attraction integral - T. Helgaker, Molecular
            !    Electronic-Structure Theory, eq. 9.9.32
            !
            real(F64), dimension(:, :), intent(out) :: V
            type(TAOBasis), intent(in)              :: AOBasis
            type(TSystem), intent(in)               :: System

            integer :: ShellA, ShellB, ShellParamsA, ShellParamsB, ShellAB
            integer :: La, Na, a0, a1
            integer :: Lb, Nb, b0, b1
            integer :: i, j, k, l
            integer :: Lab, Nab
            real(F64), dimension(3) :: Ra, Rb, Rp, Rc, Rpa, Rpb, Rpc, Rab
            real(F64) :: Prefactor, AlphaA, AlphaB, AlphaAB, Mu
            integer, parameter :: MaxNFunc = ((AUTO2E_MAXL + 1) * (AUTO2E_MAXL + 2)) / 2
            real(F64), dimension(MaxNFunc**2) :: Vab
            real(F64), dimension(2*AUTO2e_MAXL+1) :: fmarray
            integer, parameter :: MaxRtuv = 2 * AUTO2E_MAXL
            !
            ! Dimension = Sum(Sum(Sum(1, t=0,MaxRtuv-u-v), u=0, MaxRtuv-v), v=0, MaxRtuv)
            !
            real(F64), dimension(((MaxRtuv+1)*(MaxRtuv**2+5*MaxRtuv+6))/6) :: Rtuv, RtuvC
            !
            ! Dimension = Sum(k = 0, MaxRtuv) (k + 1)**2
            !
            real(F64), dimension((2*MaxRtuv**3+9*MaxRtuv**2+13*MaxRtuv+6)/6) :: ExAB, EyAB, EzAB
            integer :: NAtomicCharges
            real(F64) :: Qc, x
            type(TAuto2eRtuv), dimension(20) :: Auto2eRtuv
            type(TAuto2eBoys), dimension(6) :: Auto2eBoys

            Auto2eRtuv(1)%ptr => auto2e_Rtuv_1
            Auto2eRtuv(2)%ptr => auto2e_Rtuv_2
            Auto2eRtuv(3)%ptr => auto2e_Rtuv_3
            Auto2eRtuv(4)%ptr => auto2e_Rtuv_4
            Auto2eRtuv(5)%ptr => auto2e_Rtuv_5

            Auto2eRtuv(6)%ptr => auto2e_Rtuv_6
            Auto2eRtuv(7)%ptr => auto2e_Rtuv_7
            Auto2eRtuv(8)%ptr => auto2e_Rtuv_8
            Auto2eRtuv(9)%ptr => auto2e_Rtuv_9
            Auto2eRtuv(10)%ptr => auto2e_Rtuv_10

            Auto2eRtuv(11)%ptr => auto2e_Rtuv_11
            Auto2eRtuv(12)%ptr => auto2e_Rtuv_12
            Auto2eRtuv(13)%ptr => auto2e_Rtuv_13
            Auto2eRtuv(14)%ptr => auto2e_Rtuv_14
            Auto2eRtuv(15)%ptr => auto2e_Rtuv_15

            Auto2eRtuv(16)%ptr => auto2e_Rtuv_16
            Auto2eRtuv(17)%ptr => auto2e_Rtuv_17
            Auto2eRtuv(18)%ptr => auto2e_Rtuv_18
            Auto2eRtuv(19)%ptr => auto2e_Rtuv_19
            Auto2eRtuv(20)%ptr => auto2e_Rtuv_20

            Auto2eBoys(1)%ptr => f1
            Auto2eBoys(2)%ptr => f2
            Auto2eBoys(3)%ptr => f3
            Auto2eBoys(4)%ptr => f4
            Auto2eBoys(5)%ptr => f5
            Auto2eBoys(6)%ptr => f6

            if (2*AOBasis%LmaxGTO > 20) then
                  call msg("Angular momentum not supported in the one-electron integrals subroutine", MSG_ERROR)
                  error stop
            end if

            V = ZERO            
            associate ( &
                  ShellParamsIdx => AOBasis%ShellParamsIdx, &
                  ShellMomentum => AOBasis%ShellMomentum, &
                  NAngFunc => AOBasis%NAngFuncCart, &
                  ShellLoc => AOBasis%ShellLocCart, &
                  NShells => AOBasis%NShells, &
                  CartPolyX => AOBasis%CartPolyX, &
                  CartPolyY => AOBasis%CartPolyY, &
                  CartPolyZ => AOBasis%CartPolyZ, &
                  LmaxGTO => AOBasis%LmaxGTO, &
                  ShellCenters => AOBasis%ShellCenters, &
                  AtomCoords => AOBasis%AtomCoords, &
                  CntrCoeffs => AOBasis%CntrCoeffs, &
                  Exponents => AOBasis%Exponents, &
                  NPrimitives => AOBasis%NPrimitives, &
                  NormFactors => AOBasis%NormFactorsCart, &
                  RealAtoms => System%RealAtoms, &
                  ECPCharges => System%ECPCharges, &
                  ZNumbers => System%ZNumbers, &
                  ZNumbersECP => System%ZNumbersECP, &
                  PointCharges => System%PointCharges, &
                  PointChargeCoords => System%PointChargeCoords, &
                  NPointCharges => System%NPointCharges &
                  )
                  !
                  ! Number of non-dummy atoms
                  !
                  NAtomicCharges = 0
                  do i = 1, 2
                        NAtomicCharges = NAtomicCharges + RealAtoms(2, i) - RealAtoms(1, i) + 1
                  end do
                  !$omp parallel &
                  !$omp private(ShellA, ShellParamsA, La, Na, a0, a1, Ra) &
                  !$omp private(ShellB, ShellParamsB, Lb, Nb, b0, b1, Rb) &
                  !$omp private(ShellAB) &
                  !$omp private(AlphaA, AlphaB, Mu, Prefactor) &
                  !$omp private(ExAB, EyAB, EzAB, Qc) &
                  !$omp private(i, j, k, l) &
                  !$omp private(Lab, Nab, AlphaAB, x) &
                  !$omp private(Vab) &
                  !$omp private(Rp, Rc, Rpa, Rpb, Rpc, Rab) &
                  !$omp private(fmarray, Rtuv, RtuvC) &
                  !$omp shared(Auto2eRtuv, Auto2eBoys) &
                  !$omp shared(V) &
                  !$omp default(shared)
                  !$omp do schedule(dynamic)
                  ShellsAB: do ShellAB = 1, (NShells * (NShells + 1)) / 2
                        call ints1e_decode_pq(ShellAB, NShells, ShellA, ShellB)
                        ShellParamsA = ShellParamsIdx(ShellA)
                        La = ShellMomentum(ShellParamsA)
                        Na = NAngFunc(ShellParamsA)                        
                        a0 = ShellLoc(ShellA)
                        a1 = ShellLoc(ShellA) + Na - 1

                        ShellParamsB = ShellParamsIdx(ShellB)
                        Lb = ShellMomentum(ShellParamsB)
                        Nb = NAngFunc(ShellParamsB)
                        b0 = ShellLoc(ShellB)
                        b1 = ShellLoc(ShellB) + Nb - 1

                        Ra = AtomCoords(:, ShellCenters(ShellA))
                        Rb = AtomCoords(:, ShellCenters(ShellB))
                        Rab = Ra - Rb

                        Lab = La + Lb
                        Nab = Na * Nb
                        Vab = ZERO
                        PrimitivesA: do j = 1, NPrimitives(ShellParamsB)
                              PrimitivesB: do i = 1, NPrimitives(ShellParamsA)
                                    AlphaA = Exponents(i, ShellParamsA)
                                    AlphaB = Exponents(j, ShellParamsB)
                                    AlphaAB = AlphaA + AlphaB
                                    Mu = AlphaA * AlphaB / (AlphaA + AlphaB)
                                    Prefactor = CntrCoeffs(i, ShellParamsA) * CntrCoeffs(j, ShellParamsB) &
                                          * exp(-Mu * dot_product(Rab, Rab)) * TWO * PI / AlphaAB
                                    Rp = (AlphaA*Ra + AlphaB*Rb) / AlphaAB
                                    Rpa = Rp - Ra
                                    Rpb = Rp - Rb
                                    !
                                    ! Seed for calculating E^{ij}_t coeffs is 1.d+0 instead of
                                    ! exp(-alpha_reduces * xabsq) as in Helgaker's textbook
                                    ! because the coeffs are linear functions of the seed so
                                    ! it may be incorporated into the prefactor.
                                    !
                                    call ints1e_eijmatrix(Lab, ONE, AlphaAB, Rpa(1), Rpb(1), ExAB)
                                    call ints1e_eijmatrix(Lab, ONE, AlphaAB, Rpa(2), Rpb(2), EyAB)
                                    call ints1e_eijmatrix(Lab, ONE, AlphaAB, Rpa(3), Rpb(3), EzAB)
                                    Rtuv = ZERO
                                    do k = 1, 2
                                          do l = RealAtoms(1, k), RealAtoms(2, k)
                                                Rc = AtomCoords(:, l)
                                                if (ECPCharges) then
                                                      Qc = real(ZNumbersECP(l), F64)
                                                else
                                                      Qc = real(ZNumbers(l), F64)
                                                end if
                                                Rpc = Rp - Rc
                                                x = AlphaAB * dot_product(Rpc, Rpc)
                                                call ints1e_Coulomb_core(Rtuv, RtuvC, fmarray, x, Qc, &
                                                      Lab, AlphaAB, Rpc, Auto2eRtuv, Auto2eBoys)
                                          end do
                                    end do
                                    
                                    do l = 1, NPointCharges
                                          Rc = PointChargeCoords(:, l)
                                          Qc = PointCharges(l)
                                          Rpc = Rp - Rc
                                          x = AlphaAB * dot_product(Rpc, Rpc)
                                          call ints1e_Coulomb_core(Rtuv, RtuvC, fmarray, x, Qc, &
                                                Lab, AlphaAB, Rpc, Auto2eRtuv, Auto2eBoys)
                                    end do

                                    call ints1e_HermiteTransf(Vab, Rtuv, La, Lb, Na, Nb, ExAB, EyAB, EzAB, &
                                          Prefactor, NormFactors(:, ShellParamsA), NormFactors(:, ShellParamsB), &
                                          CartPolyX, CartPolyY, CartPolyZ)
                              end do PrimitivesB
                        end do PrimitivesA
                        V(a0:a1, b0:b1) = reshape(Vab(1:Nab), [Na, Nb])
                  end do ShellsAB
                  !$omp end do
                  !$omp end parallel
            end associate
      end subroutine ints1e_Coulomb


      subroutine ints1e_Kinetic_core(Tab, ExAB, EyAB, EzAB, Prefactor, La, Lb, &
            Na, Nb, NormA, NormB, AlphaB, CartPolyX, CartPolyY, CartPolyZ)
            !
            ! Calculate kinetic energy integral between unnormalized gauss orbitals.
            ! J.T. Ferman, E.F. Valeev, Fundamentals of Molecular Integral Evaluation, eq. 4.8
            !
            real(F64), dimension(Na, Nb), intent(inout) :: Tab
            real(F64), dimension(:), intent(in)         :: ExAB
            real(F64), dimension(:), intent(in)         :: EyAB
            real(F64), dimension(:), intent(in)         :: EzAB
            real(F64), intent(in)                       :: Prefactor
            integer, intent(in)                         :: La
            integer, intent(in)                         :: Lb
            integer, intent(in)                         :: Na
            integer, intent(in)                         :: Nb
            real(F64), dimension(:), intent(in)         :: NormA
            real(F64), dimension(:), intent(in)         :: NormB
            integer, dimension(:, 0:), intent(in)       :: CartPolyX
            integer, dimension(:, 0:), intent(in)       :: CartPolyY
            integer, dimension(:, 0:), intent(in)       :: CartPolyZ
            real(F64), intent(in)                       :: AlphaB

            real(F64) :: Ix, Iy, Iz
            real(F64) :: int0, int1
            real(F64) :: dlxB, dlyB, dlzB
            real(F64) :: xy, yz, zx, Sab
            integer :: a, b
            integer :: x0, y0, z0
            integer :: lxA, lyA, lzA, lxB, lyB, lzB
            real(F64) :: F

            do b = 1, Nb
                  do a = 1, Na
                        lxA = CartPolyX(a, La)
                        lyA = CartPolyY(a, La)
                        lzA = CartPolyZ(a, La)
                        lxB = CartPolyX(b, Lb)
                        lyB = CartPolyY(b, Lb)
                        lzB = CartPolyZ(b, Lb)

                        x0 = ints1e_ELoc(lxA, lxB, 0)
                        y0 = ints1e_ELoc(lyA, lyB, 0)
                        z0 = ints1e_ELoc(lzA, lzB, 0)

                        F = Prefactor * NormA(a) * NormB(b)
                        xy = ExAB(x0) * EyAB(y0)
                        yz = EyAB(y0) * EzAB(z0)
                        zx = EzAB(z0) * ExAB(x0)
                        Sab = ExAB(x0) * EyAB(y0) * EzAB(z0)
                        !
                        ! Eq. 4.8
                        !
                        ! T_{12} = I_x + I_y + I_z
                        ! <+2|_x = x^{l + 2} y^m z^n \exp(-\alpha r^2)
                        ! I_x = \alpha_2 * (2 * l2 + 1) <0|0> - 2 \alpha_2^2 <0|+2>_x - l2(l2-1)/2 <0|-2>_x
                        !
                        dlxB = real(lxB, F64)
                        dlyB = real(lyB, F64)
                        dlzB = real(lzB, F64)

                       if (lxB >= 2) then
                             int0 = ExAB(ints1e_ELoc(lxA, lxB-2, 0)) * yz
                             int1 = ExAB(ints1e_ELoc(lxA, lxB+2, 0)) * yz
                             Ix = (TWO * dlxB + ONE) * AlphaB*Sab - TWO*AlphaB**2 * int1 - dlxB * (dlxB - ONE) / TWO * int0
                       else
                             int1 = ExAB(ints1e_ELoc(lxA, lxB+2, 0)) * yz
                             Ix = (TWO * dlxB + ONE) * AlphaB*Sab - TWO*AlphaB**2 * int1
                       end if

                        if (lyB >= 2) then
                              int0 = EyAB(ints1e_ELoc(lyA, lyB-2, 0)) * zx
                              int1 = EyAB(ints1e_ELoc(lyA, lyB+2, 0)) * zx
                              Iy = (TWO * dlyB + ONE) * AlphaB*Sab - TWO*AlphaB**2 * int1 - dlyB * (dlyB - ONE) / TWO * int0
                        else
                              int1 = EyAB(ints1e_ELoc(lyA, lyB+2, 0)) * zx
                              Iy = (TWO * dlyB + ONE) * AlphaB*Sab - TWO*AlphaB**2 * int1
                        end if

                        if (lzB >= 2) then
                              int0 = EzAB(ints1e_ELoc(lzA, lzB-2, 0)) * xy
                              int1 = EzAB(ints1e_ELoc(lzA, lzB+2, 0)) * xy
                              Iz = (TWO * dlzB + ONE) * AlphaB*Sab - TWO*AlphaB**2 * int1 - dlzB * (dlzB - ONE) / TWO * int0
                        else
                              int1 = EzAB(ints1e_ELoc(lzA, lzB+2, 0)) * xy
                              Iz = (TWO * dlzB + ONE) * AlphaB*Sab - TWO*AlphaB**2 * int1
                        end if

                        Tab(a, b) = Tab(a, b) + F * (Ix + Iy + Iz)
                  end do
            end do
      end subroutine ints1e_Kinetic_core
      

      subroutine ints1e_Kinetic(T, AOBasis)
            !
            ! Kinetic matrix in a Cartesian GTO basis
            !
            real(F64), dimension(:, :), intent(out) :: T
            type(TAOBasis), intent(in)              :: AOBasis

            integer :: ShellA, ShellB, ShellParamsA, ShellParamsB, ShellAB
            integer :: La, Na, a0, a1
            integer :: Lb, Nb, b0, b1
            integer :: i, j
            integer :: Lab, Nab
            real(F64), dimension(3) :: Ra, Rb, Rp, Rc, Rpa, Rpb, Rab
            real(F64) :: Prefactor, AlphaA, AlphaB, AlphaAB, Mu
            integer, parameter :: MaxNFunc = ((AUTO2E_MAXL + 1) * (AUTO2E_MAXL + 2)) / 2
            real(F64), dimension(MaxNFunc**2) :: Tab
            integer, parameter :: MaxIndex = 2 * AUTO2E_MAXL + 2
            !
            ! Dimension = Sum(k = 0, MaxIndex) (k + 1)**2
            !
            real(F64), dimension((2*MaxIndex**3+9*MaxIndex**2+13*MaxIndex+6)/6) :: ExAB, EyAB, EzAB

            T = ZERO            
            associate ( &
                  ShellParamsIdx => AOBasis%ShellParamsIdx, &
                  ShellMomentum => AOBasis%ShellMomentum, &
                  NAngFunc => AOBasis%NAngFuncCart, &
                  ShellLoc => AOBasis%ShellLocCart, &
                  NShells => AOBasis%NShells, &
                  CartPolyX => AOBasis%CartPolyX, &
                  CartPolyY => AOBasis%CartPolyY, &
                  CartPolyZ => AOBasis%CartPolyZ, &
                  LmaxGTO => AOBasis%LmaxGTO, &
                  ShellCenters => AOBasis%ShellCenters, &
                  AtomCoords => AOBasis%AtomCoords, &
                  CntrCoeffs => AOBasis%CntrCoeffs, &
                  Exponents => AOBasis%Exponents, &
                  NPrimitives => AOBasis%NPrimitives, &
                  NormFactors => AOBasis%NormFactorsCart &
                  )
                  !$omp parallel &
                  !$omp private(ShellA, ShellParamsA, La, Na, a0, a1, Ra) &
                  !$omp private(ShellB, ShellParamsB, Lb, Nb, b0, b1, Rb) &
                  !$omp private(ShellAB) &
                  !$omp private(AlphaA, AlphaB, Mu, Prefactor) &
                  !$omp private(ExAB, EyAB, EzAB) &
                  !$omp private(i, j) &
                  !$omp private(Lab, Nab, AlphaAB) &
                  !$omp private(Tab) &
                  !$omp private(Rp, Rc, Rpa, Rpb, Rab) &
                  !$omp shared(T) &
                  !$omp default(shared)
                  !$omp do schedule(dynamic)
                  ShellsAB: do ShellAB = 1, (NShells * (NShells + 1)) / 2
                        call ints1e_decode_pq(ShellAB, NShells, ShellA, ShellB)
                        ShellParamsA = ShellParamsIdx(ShellA)
                        La = ShellMomentum(ShellParamsA)
                        Na = NAngFunc(ShellParamsA)                        
                        a0 = ShellLoc(ShellA)
                        a1 = ShellLoc(ShellA) + Na - 1

                        ShellParamsB = ShellParamsIdx(ShellB)
                        Lb = ShellMomentum(ShellParamsB)
                        Nb = NAngFunc(ShellParamsB)
                        b0 = ShellLoc(ShellB)
                        b1 = ShellLoc(ShellB) + Nb - 1

                        Ra = AtomCoords(:, ShellCenters(ShellA))
                        Rb = AtomCoords(:, ShellCenters(ShellB))
                        Rab = Ra - Rb

                        Lab = La + Lb
                        Nab = Na * Nb
                        Tab = ZERO
                        PrimitivesA: do j = 1, NPrimitives(ShellParamsB)
                              PrimitivesB: do i = 1, NPrimitives(ShellParamsA)
                                    AlphaA = Exponents(i, ShellParamsA)
                                    AlphaB = Exponents(j, ShellParamsB)
                                    AlphaAB = AlphaA + AlphaB
                                    Mu = AlphaA * AlphaB / (AlphaA + AlphaB)
                                    Prefactor = CntrCoeffs(i, ShellParamsA) * CntrCoeffs(j, ShellParamsB) &
                                          * exp(-Mu * dot_product(Rab, Rab)) * (PI / AlphaAB)**(THREE/TWO)
                                    Rp = (AlphaA*Ra + AlphaB*Rb) / AlphaAB
                                    Rpa = Rp - Ra
                                    Rpb = Rp - Rb                                    
                                    !
                                    ! Seed for calculating E^{ij}_t coeffs is 1.d+0 instead of
                                    ! exp(-alpha_reduced * xabsq) as in Helgaker's textbook
                                    ! because the coeffs are linear functions of the seed so
                                    ! it may be incorporated into the prefactor.
                                    !
                                    call ints1e_eijmatrix(Lab+2, ONE, AlphaAB, Rpa(1), Rpb(1), ExAB)
                                    call ints1e_eijmatrix(Lab+2, ONE, AlphaAB, Rpa(2), Rpb(2), EyAB)
                                    call ints1e_eijmatrix(Lab+2, ONE, AlphaAB, Rpa(3), Rpb(3), EzAB)
                                    call ints1e_Kinetic_core(Tab, ExAB, EyAB, EzAB, Prefactor, La, Lb, &
                                          Na, Nb, NormFactors(:, ShellParamsA), NormFactors(:, ShellParamsB), &
                                          AlphaB, CartPolyX, CartPolyY, CartPolyZ)
                              end do PrimitivesB
                        end do PrimitivesA
                        T(a0:a1, b0:b1) = reshape(Tab(1:Nab), [Na, Nb])
                  end do ShellsAB
                  !$omp end do
                  !$omp end parallel
            end associate
      end subroutine ints1e_Kinetic


      subroutine ints1e_S(S, AOBasis)
            !
            ! Overlap matrix of spherical Gaussian atomic orbitals.
            ! Both upper and lower triangles of the output matrix
            ! are filled with data.
            !
            real(F64), dimension(:, :), intent(out) :: S
            type(TAOBasis), intent(in)              :: AOBasis

            real(F64), dimension(:, :), allocatable :: S_cao
            integer :: NAOCart

            NAOCart = AOBasis%NAOCart
            allocate(S_cao(NAOCart, NAOCart))
            call ints1e_OverlapMatrix(S_cao, AOBasis)
            call real_smfill(S_cao)
            call ints1e_SpherAOTransf(S, S_cao, AOBasis)
      end subroutine ints1e_S


      subroutine ints1e_T(T, AOBasis)
            !
            ! Kinetic energy operator matrix in the basis of spherical
            ! Gaussian atomic orbitals.
            ! Both upper and lower triangles of the output matrix
            ! are filled with data.
            !
            real(F64), dimension(:, :), intent(out) :: T
            type(TAOBasis), intent(in)              :: AOBasis

            real(F64), dimension(:, :), allocatable :: T_cao
            integer :: NAOCart

            NAOCart = AOBasis%NAOCart
            allocate(T_cao(NAOCart, NAOCart))
            call ints1e_Kinetic(T_cao, AOBasis)
            call real_smfill(T_cao)
            call ints1e_SpherAOTransf(T, T_cao, AOBasis)
      end subroutine ints1e_T


      subroutine ints1e_Vne(Vne, AOBasis, System)
            !
            ! One-electron Coulomb matrix in the spherical Gaussian
            ! in the basis of spherical Gaussian atomic orbitals.
            ! Both upper and lower triangles of the output matrix
            ! are filled with data.
            !
            real(F64), dimension(:, :), intent(out) :: Vne
            type(TAOBasis), intent(in)              :: AOBasis
            type(TSystem), intent(in)               :: System

            real(F64), dimension(:, :), allocatable :: Vne_cao
            integer :: NAOCart

            NAOCart = AOBasis%NAOCart
            allocate(Vne_cao(NAOCart, NAOCart))
            call ints1e_Coulomb(Vne_cao, AOBasis, System)
            call real_smfill(Vne_cao)
            call ints1e_SpherAOTransf(Vne, Vne_cao, AOBasis)
      end subroutine ints1e_Vne
      
      
      subroutine ints1e_SpherAOTransf(X_sao, X_cao, AOBasis)
            real(F64), dimension(:, :), intent(out) :: X_sao
            real(F64), dimension(:, :), intent(in)  :: X_cao
            type(TAOBasis), intent(in)              :: AOBasis

            integer :: NAOSpher, NAOCart
            real(F64), dimension(:), allocatable :: TransfWork
            
            NAOSpher = AOBasis%NAOSpher
            NAOCart = AOBasis%NAOCart
            allocate(TransfWork(NAOSpher*NAOCart))
            call SpherGTO_TransformMatrix_U(X_sao, X_cao, &
                  AOBasis%LmaxGTO, &
                  AOBasis%NormFactorsSpher, &
                  AOBasis%NormFactorsCart, &
                  AOBasis%ShellLocSpher, &
                  AOBasis%ShellLocCart, &
                  AOBasis%ShellMomentum, &
                  AOBasis%ShellParamsIdx, &
                  AOBasis%NAOSpher, &
                  AOBasis%NAOCart, &
                  AOBasis%NShells, TransfWork)
      end subroutine ints1e_SpherAOTransf
end module OneElectronInts
