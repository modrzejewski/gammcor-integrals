#!/usr/bin/env python3
# -------------------------------------------------------------------
#    Subroutines for efficient transformation of atomic orbitals
#    from the Cartesian to the real spherical harmonics basis
# -------------------------------------------------------------------
from scipy.special import comb
from math import sqrt, factorial, pi
from scipy.special import factorial2 as dblfact
import numpy as np
from Auto2eGlobal import *

SMALL_COEFF_THRESH = 10 * np.finfo(float).eps

def XYZInt_Eq15(l, lx0, ly0, lz0):
    #
    # Compute the integral of x^{lx+lx0} y^{ly+ly0} z^{lz+lz0} over
    # the unit sphere:
    # Int_Omega x^{lx+lx0} y^{ly+ly0} z^{lz+lz0} dOmega.
    # (See Eq. 15 in [1]) The tntegrals are computed for the fixed
    # lx0, ly0, lz0, and for all (LX,LY,LZ) for which LX+LY+LZ=L.
    # COMMON PITFALL: Do not assume that this subroutine computes
    # integrals for angular momenta lower than L.
    #
    # 1. Moreno-Flores, R., Alvarez-Mendez, R., Vela, A., and
    #    Koster, A.M., Half-Numerical Evaluation of Pseudopotential
    #    Integrals, J. Comput. Chem. 27, 1009 (2006);
    #    doi: 10.1002/jcc.20410
    #
    l0 = lx0 + ly0 + lz0
    d = dblfact(l+l0+1, exact=True)
    xyz = np.zeros((l+1,l+1,l+1))
    for lx in range(l+1):
        x = lx + lx0
        a = 4 * pi * dblfact(x-1, exact=True)
        for ly in range(l - lx + 1):
            lz = l - lx - ly
            y = ly + ly0
            z = lz + lz0
            if (x % 2 == 1) or (y % 2 == 1) or (z % 2 == 1):
                xyz[lx][ly][lz] = 0
            else:
                b = dblfact(y-1, exact=True)
                c = dblfact(z-1, exact=True)
                xyz[lx][ly][lz] = a * b * c / d
    return xyz
                
                
def RHSV(l, m, j):
    Uxyz = RSHU(l, m, l)
    Vxyz = np.zeros((j+1,j+1,j+1))
    for lx in range(j+1):
        for ly in range(j-lx+1):
            lz = j - lx - ly
            xyz = XYZInt_Eq15(l, lx, ly, lz)
            for mx in range(l+1):
                for my in range(l-mx+1):
                    mz = l - mx - my
                    Vxyz[lx][ly][lz] += Uxyz[mx][my][mz] * xyz[mx][my][mz]
    return Vxyz
    
                
def Nslm(l, m):
    #
    # Compute the normalization factor Nlm of the solid harmonic Slm
    # in the basis of x**lx y**ly z**lz functions.
    #
    # Nlm in Racah's normalization is defined in Eq. 9.1.10 of Ref. 1,
    # but we multiply it by  Sqrt(2l+1/4Pi) to change the normalization
    # integral to unity.
    #
    # 1. Helgaker, T., Jorgensen, P., Olsen, J., Molecular
    #    Electronic-Structure Theory, Wiley & Sons Chichester
    #    2000. Eqs. 6.4.47-6.4.50, p. 215
    #
    sqra = sqrt(factorial(l+abs(m)))
    sqrb = sqrt(factorial(l-abs(m)))
    if m == 0:
        sqrc = 1
    else:
        sqrc = sqrt(2)
    #
    # Real solid harmonics are normalized to 1
    # instead of Racah's normalization in Helgaker's
    # textbook
    #
    a = sqrt((2*l+1)/(4*pi))
    return a * sqra * sqrb * sqrc / (2**abs(m) * factorial(l))

    
def Clmtuv(l, m, t, u, v):
    #
    # Compute the expansion coefficeint Clm(t, u, v) of the solid harmonic Slm
    # in the basis of x**lx y**ly z**lz functions. Clm(t, u, v) is defined in
    # Eq. 9.1.11 of Ref. 1.
    #
    # Note that the Slm functions in Helgaker's textbook are in Racah's
    # normalization, whereas the harmonics computed in this module are
    # normalized to unity.
    #
    # 1. Helgaker, T., Jorgensen, P., Olsen, J., Molecular
    #    Electronic-Structure Theory, Wiley & Sons Chichester
    #    2000. Eqs. 6.4.47-6.4.50, p. 215
    #
    if m >= 0:
        vm = 0
    else:
        vm = 1
    AbsM = abs(m)
    a = comb(l, t, exact=True) * comb(l-t, AbsM+t, exact=True)
    b = comb(t, u, exact=True) * comb(AbsM, v, exact=True)
    c = (-1)**(t+(v-vm)//2) * 4**t
    return a * b / c


def RSHU_l_eq_j(l, m):
    #
    # Calculate coefficient of the X^{lx}Y^{ly}Z^{lz} expansion
    # of real spherical orthonormal harmonics: U^{lm}_{lxlylz}
    # coefficient in Eq. 32 and 33 in [2].
    #
    # The implementation uses the equations given in Helgaker's
    # textbook (Eqs. 9.1.9-12 in Ref. 1), but rescaled by the factor
    # of Sqrt(2l+1/4Pi) to change the normalization integral to
    # unity.
    #
    # The complete array of coefficients, for all lx + ly + lz = l
    # is computed in a single call.
    #
    # Definition of real spherical harmonics:
    # For m > 0: Sl^m = (-1)^m/Sqrt(2) Yl^m + 1/Sqrt(2)       Yl^(-m)
    # For m < 0: Sl^m =      i/Sqrt(2) Yl^m - i(-1)^m/Sqrt(2) Yl^(-m)
    # For m = 0: Sl^m = Yl^m
    # Yl^m are orthonormal, complex-valued spherical harmonics.
    #
    # 1. Helgaker, T., Jorgensen, P., Olsen, J., Molecular
    #    Electronic-Structure Theory, Wiley & Sons Chichester
    #    2000. Eqs. 6.4.47-6.4.50, p. 215
    # 2. Moreno-Flores, R., Alvarez-Mendez, R., Vela, A., and
    #    Koster, A.M., Half-Numerical Evaluation of Pseudopotential
    #    Integrals, J. Comput. Chem. 27, 1009 (2006);
    #    doi: 10.1002/jcc.20410
    #
    AbsM = abs(m)
    npref = Nslm(l, m)
    Uxyz = np.zeros((l+1,l+1,l+1))
    if m >= 0:
        vm = 0
    else:
        vm = 1
    for t in range((l-AbsM)//2+1):
        for u in range(t+1):
            for k in range((AbsM-vm)//2+1):
                v = 2 * k + vm
                lx = 2 * t + AbsM - 2 * u - v
                ly = 2 * u + v
                lz = l - 2 * t - AbsM
                Uxyz[lx][ly][lz] += npref * Clmtuv(l, m, t, u, v)
    return Uxyz


def RSHU(l, m, j):
    Ulm = RSHU_l_eq_j(l, m)
    if j == l:
        return Ulm
    Uxyz = np.zeros((j+1,j+1,j+1))
    if (j - l) % 2 == 1:
        return Uxyz
    kappa = (j - l) // 2
    for lx in range(l+1):
        for ly in range(l-lx+1):
            for mx in range(kappa+1):
                for my in range(kappa-mx+1):
                    lz = l - lx - ly
                    mz = kappa - mx - my
                    u = lx + 2 * mx
                    v = ly + 2 * my
                    w = lz + 2 * mz
                    Uxyz[u][v][w] += (factorial(kappa)/(factorial(mx)*factorial(my)*factorial(mz)) * Ulm[lx][ly][lz])
    return Uxyz


def RSHUV(LMax):
    #
    # Compute the Cartesian->Spherical transformation matrices
    # up to angular momentum LMax.
    #
    # The matrix U[(l,m,j)] contains coefficients of the transformation
    # from Cartesian polynomials x**lx * y**ly * z**lz, lx+ly+lz=j to
    # real spherical harmonics Ylm normalized to one on
    # the unit sphere.
    #
    # Ylm = Sum(lx+ly+lz=j) U[(l,m,j)][lx][ly][lz] x**lx * y**ly * z**lz
    #
    # The definition of the coefficients U[lx][ly][lz] is numerically
    # equivalent to the definition in Ref. 1, but it's coded according
    # to the Helgaker's Molecular electronic structure theory. Additionally,
    # the subroutines provided in this module can compute the coefficients
    # of contaminant spherical harmonics, i.e., the functions for which
    # l < j. (Note that the coefficients U[(l,m,j)] are nonzero only if
    # (j-l)%2==0)
    #
    # The matrix V is the inverse of U and should be used to transform
    # the matrices which act on the atomic orbitals:
    #
    # Sum(l,m) U[(l,m,j)][mx][my][mz] * U[(l,m,j)][nx][nx][nz] =
    # delta(mx,nx) * delta(my,ny) * delta(mz,nz)
    #
    #
    # 1. Moreno-Flores, R., Alvarez-Mendez, R., Vela, A., and
    #    Koster, A.M., Half-Numerical Evaluation of Pseudopotential
    #    Integrals, J. Comput. Chem. 27, 1009 (2006);
    #    doi: 10.1002/jcc.20410
    #
    # 2. Helgaker, T., Jorgensen, P., Olsen, J., Molecular
    #    Electronic-Structure Theory, Wiley & Sons Chichester
    #    2000. Eqs. 6.4.47-6.4.50, p. 215
    #
    U = {}
    V = {}
    for j in range(LMax+1):
        for l in range(j+1):
            for m in range(-l,l+1):
                Uxyz = RSHU(l, m, j)
                Vxyz = RHSV(l, m, j)
                U[(l, m, j)] = Uxyz
                V[(l, m, j)] = Vxyz
                
    return U, V


def lxlylzpos(lx, ly, lz):
    l = lx + ly + lz
    return ((2 * l - lx + 3) * lx) // 2 + ly + 1


def test_UV(J):
    #
    # Test if the matrix product of the transformation
    # matrices U and V is the identity matrix.
    #
    U, V = RSHUV(J)
    Nxyz = ((J+1)*(J+2))//2
    AB = np.zeros((Nxyz, Nxyz))
    for mx in range(J+1):
        for my in range(J-mx+1):
            mz = J - mx - my
            for nx in range(J+1):
                for ny in range(J-nx+1):
                    nz = J - nx - ny
                    a = lxlylzpos(mx, my, mz) - 1
                    b = lxlylzpos(nx, ny, nz) - 1
                    uv = 0
                    for l in range(J+1):
                        for m in range(-l, l+1):
                            uv += V[(l, m, J)][mx][my][mz] * U[(l, m, J)][nx][ny][nz]
                    AB[a][b] = uv
                    
    np.set_printoptions(suppress=True, precision=8, linewidth=160)
    print(AB)


def SpherTransf_1e(Yout, Xin, l, UCoeffs=True, Vectorized=True):
    U, V = RSHUV(l)
    if UCoeffs:
        coeffs = U
    else:
        coeffs = V
    s = ""
    for m in range(-l, l+1):
        if Vectorized:
            yvec = "{Y}(:, {iy})".format(Y=Yout, iy=l+m+1)
        else:
            yvec = "{Y}({iy})".format(Y=Yout, iy=l+m+1)
        cvec = []
        xvec = []
        for x, y, z in AngularFunctions[l]:
            uxyz = coeffs[(l, m, l)][x][y][z]
            if abs(uxyz) > SMALL_COEFF_THRESH:
                cvec.append(uxyz)
                if Vectorized:
                    xvec.append("{X}(:, {ix})".format(ix=XYZIdx[(x,y,z)], X=Xin))
                else:
                    xvec.append("{X}({ix})".format(ix=XYZIdx[(x,y,z)], X=Xin))
        s += SplitLine(FormatSum(yvec, cvec, xvec), ["+", "-"]) + "\n"
    return s


def SpherTransf_1e_Subroutine(l, UCoeffs=True, Vectorized=True):
    if Vectorized:
        arguments = """real(F64), dimension(n, *), intent(out) :: S
real(F64), dimension(n, *), intent(in) :: C
integer, intent(in) :: n"""
    else:
        arguments = """real(F64), dimension(*), intent(out) :: S
real(F64), dimension(*), intent(in) :: C"""

    if l >= SPHER_TRANSF_LMIN:
        transformation_code = SpherTransf_1e("S", "C", l, UCoeffs, Vectorized)
    else:
        if Vectorized:
            transformation_code = "S(:, 1:{NSpher}) = C(:, 1:{NSpher})".format(NSpher=2*l+1)
        else:
            transformation_code = "S(1:{NSpher}) = C(1:{NSpher})".format(NSpher=2*l+1)

    s = """
subroutine auto1e_SpherTransf{U}{Vector}_{l}(S, C{VecLen})
{arguments}
{transf}
end subroutine auto1e_SpherTransf{U}{Vector}_{l}
""".format(l=l, transf=transformation_code.strip(),
           U="_U" if UCoeffs else "_V", Vector="_Vector" if Vectorized else "",
           arguments=arguments, VecLen = ", n" if Vectorized else "")
    return s


def Cartesian_Dummy_Transf(l):
    a = {}
    m = -l
    for x, y, z in AngularFunctions[l]:
        t = np.zeros((l+1, l+1, l+1))
        t[x][y][z] = 1
        a[(l, m, l)] = t
        m += 1
    return a
    

def SpherTransf_2e(Yout, Xin, lA, lB, UCoeffs=True, Vectorized=True):
    Ua, Va = RSHUV(lA)
    Ub, Vb = RSHUV(lB)
    if UCoeffs:
        coeffsA = Ua
        coeffsB = Ub
    else:
        coeffsA = Va
        coeffsB = Vb
    if lA < SPHER_TRANSF_LMIN:
        coeffsA = Cartesian_Dummy_Transf(lA)
    if lB < SPHER_TRANSF_LMIN:
        coeffsB = Cartesian_Dummy_Transf(lB)
    s = ""
    for mB in range(-lB, lB+1):
        for mA in range(-lA, lA+1):
            if Vectorized:
                yvec = "{Y}(:, {iy})".format(Y=Yout, iy=GabIndexSpher(mA, lA, mB, lB))
            else:
                yvec = "{Y}({iy})".format(Y=Yout, iy=GabIndexSpher(mA, lA, mB, lB))
            cvec = []
            xvec = []
            for xB, yB, zB in AngularFunctions[lB]:
                for xA, yA, zA in AngularFunctions[lA]:
                    uxyzA = coeffsA[(lA, mA, lA)][xA][yA][zA]
                    uxyzB = coeffsB[(lB, mB, lB)][xB][yB][zB]
                    if abs(uxyzA*uxyzB) > SMALL_COEFF_THRESH:
                        cvec.append(uxyzA*uxyzB)
                        if Vectorized:
                            xvec.append("{X}(:, {ix})".format(ix=GabIndex((xA,yA,zA), (xB, yB, zB)), X=Xin))
                        else:
                            xvec.append("{X}({ix})".format(ix=GabIndex((xA,yA,zA), (xB, yB, zB)), X=Xin))
            s += SplitLine(FormatSum(yvec, cvec, xvec), ["+", "-"], TargetSize=100) + "\n"
    return s


def SpherTransf_2e_Subroutine(lA, lB, UCoeffs=True, Vectorized=True):
    NSpherA = 2 * lA + 1
    NSpherB = 2 * lB + 1
    NCartA = ((lA+1)*(lA+2))//2
    NCartB = ((lB+1)*(lB+2))//2
    if Vectorized:
        arguments = """real(F64), dimension(n, *), intent(out) :: S
real(F64), dimension(n, *), intent(in) :: C
integer, intent(in) :: n""".format(NSpher=NSpherA*NSpherB, NCart=NCartA*NCartB)
    else:
        arguments = """real(F64), dimension(*), intent(out) :: S
real(F64), dimension(*), intent(in) :: C""".format(NSpher=NSpherA*NSpherB, NCart=NCartA*NCartB)
        
    s = """
subroutine auto2e_SpherTransf{U}{Vector}_{lA}_{lB}(S, C{VecLen})
{arguments}
{transf}
end subroutine auto2e_SpherTransf{U}{Vector}_{lA}_{lB}
""".format(lA=lA, lB=lB, transf=SpherTransf_2e("S", "C", lA, lB, UCoeffs, Vectorized).strip(),
           U="_U" if UCoeffs else "_V", Vector="_Vector" if Vectorized else "",
           arguments=arguments, VecLen=", n" if Vectorized else "")
    return s


def Init_Subroutine(MaxL):
    PtrAssignments = []
    for lA in range(MaxL+1):
        PtrAssignments.append("Auto1e_SpherTransf_V({l})%ptr => auto1e_SpherTransf_V_{l}".format(l=lA))
    for lA in range(MaxL+1):
        PtrAssignments.append("Auto1e_SpherTransf_U({l})%ptr => auto1e_SpherTransf_U_{l}".format(l=lA))
    for lA in range(MaxL+1):
        PtrAssignments.append("Auto1e_SpherTransf_U_Vector({l})%ptr => auto1e_SpherTransf_U_Vector_{l}".format(l=lA))
    
    Subroutine = """
subroutine auto2e_SpherTransf_init()
{}
end subroutine auto2e_SpherTransf_init
""".format("\n".join(PtrAssignments))
    return Subroutine


def SpherTransf_Module(Lmax):
    subroutines = ""
    for lB in range(Lmax+1):
        for lA in range(Lmax+1):
            if not (lA < SPHER_TRANSF_LMIN and lB < SPHER_TRANSF_LMIN):
                subroutines += SpherTransf_2e_Subroutine(lA, lB, Vectorized=True, UCoeffs=True)

    for l in range(Lmax+1):
        subroutines += SpherTransf_1e_Subroutine(l, Vectorized=True, UCoeffs=True)
        subroutines += SpherTransf_1e_Subroutine(l, Vectorized=False, UCoeffs=True)
        subroutines += SpherTransf_1e_Subroutine(l, Vectorized=False, UCoeffs=False)

    module = """!
! Automatically generated Cartesian -> spherical transformation.
! The real spherical harmonics are normalized to 1 on the unit sphere.
! No scaling and reordering is applied for angular momenta 0 <= L < {Lmin}.
!
module auto2e_SpherTransf
use arithmetic
implicit none

type TAutoTransfVectPtr
procedure(auto1e_SpherTransf_U_Vector_0), pointer, nopass :: ptr
end type TAutoTransfVectPtr

type TAutoTransfPtr
procedure(auto1e_SpherTransf_U_0), pointer, nopass :: ptr
end type TAutoTransfPtr

type(TAutoTransfPtr), dimension(0:{Lmax}), save :: Auto1e_SpherTransf_V
type(TAutoTransfPtr), dimension(0:{Lmax}), save :: Auto1e_SpherTransf_U
type(TAutoTransfVectPtr), dimension(0:{Lmax}), save :: Auto1e_SpherTransf_U_Vector

contains

{init}

{subroutines}
end module auto2e_SpherTransf""".format(subroutines=subroutines.strip(), init=Init_Subroutine(Lmax).strip(),
                                        Lmax=Lmax, Lmin=SPHER_TRANSF_LMIN)
    return module

    
                                    
                    
