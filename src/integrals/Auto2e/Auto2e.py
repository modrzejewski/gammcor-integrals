#!/usr/bin/env python3

from numpy import array
import sys
from functools import reduce
from sympy import Symbol
import os
from os import path
import ParallelBuild
from Auto2eGlobal import *
from AutoSpherTransf import *
#
# Threshold for enabling the "RCopy" variant of the Hermite->Cartesian
# transformation of the ket pair. The first step is copying
# the Hermite Coulomb integrals (Rtuv) to contiguous memory
# locations. After that, a vectorized transformation is carried out
# for all elements at once. The algorithm is fast for large bra
# angular momenta (and saves a lot of lines of code) but slow
# for small LBra. The following optimal parameters were determined
# on a laptop with an Intel i7-6560U CPU.
#
RCopyLBraThresh = 7
RCopyLKetThresh = 1

def tuvindices(MaxIndex):
    """ Indices of the Rtuv and Stuv matrices """
    TUVIdx = {}
    Idx = 1
    for v in range(MaxIndex+1):
        for u in range(MaxIndex-v+1):
            for t in range(MaxIndex-v-u+1):
                TUVIdx[(t, u, v)] = Idx
                Idx += 1
    return TUVIdx            

def NeedsSpherTransf(lA, lB, lC, lD):
    return (max(lA, lB, lC, lD) >= SPHER_TRANSF_LMIN)

def WMatrixDims(A, B, C, D):
    """ Compute the dimension of W. Example for the (dd|dd) case:
                       |(ds|ds) (ds|fs) (ds|gs)|
            W    =     |(fs|ds) (fs|fs) (fs|gs)|
                       |(gs|ds) (gs|fs) (gs|gs)|
    """
    NBra = 0
    for l in range(A, A+B+1):
        NBra += NFuncCart(l)
    NKet = 0
    for l in range(C, C+D+1):
        NKet += NFuncCart(l)

    return (NBra, NKet)

def RMatrixDim(n):
    """ Compute the dimension of Rtuv, 0<=t+u+v<=n. """
    Idx = 0
    for v in range(n+1):
        for u in range(n-v+1):
            for t in range(n-v-u+1):
                Idx += 1
    return Idx

def WIndex(C, lC0):
    """ Translate (lx, ly, lz) to the index of the corresponding W matrix row/column. """
    lC = sum(C)
    #
    # Count how many columns are kept in memory for lC0 <= l < lC
    #
    nPrev = 0
    for l in range(lC0, lC):
        nPrev += NFuncCart(l)
    return nPrev + XYZIdx[tuple(C)]

def EnableRCopy(LBra, LKet):
    if LBra >= RCopyLBraThresh and LKet >= RCopyLKetThresh:
        return True
    else:
        return False
    
def RtuvKSubroutine(MaxIndex):
    """ Recursively compute the Hermite Coulomb integrals, Eqs. 9.9.18-20 p. 375 in Helgaker's textbook."""
    X = Symbol("X")
    Y = Symbol("Y")
    Z = Symbol("Z")
    RtuvA = "RA"
    RtuvB = "RB"
    TUVIdxA = tuvindices(MaxIndex)
    if MaxIndex > 0:
        TUVIdxB = tuvindices(MaxIndex-1)
    s = ""
    for v in range(MaxIndex+1):
        for u in range(MaxIndex-v+1):
            for t in range(MaxIndex-v-u+1):
                if t >= 1:
                    rhs = (t-1) * RtuvElement(t-2, u, v, RtuvB, TUVIdxB) + X * RtuvElement(t-1, u, v, RtuvB, TUVIdxB)
                elif u >= 1:
                    rhs = (u-1) * RtuvElement(t, u-2, v, RtuvB, TUVIdxB) + Y * RtuvElement(t, u-1, v, RtuvB, TUVIdxB)
                elif v >= 1:
                    rhs = (v-1) * RtuvElement(t, u, v-2, RtuvB, TUVIdxB) + Z * RtuvElement(t, u, v-1, RtuvB, TUVIdxB)
                else:
                    rhs = "(-2*p)**n * Fm(n+1)"
                idx = TUVIdxA[(t, u, v)]
                NextLine = "{R}({i}) = {rhs}\n".format(R=RtuvA, i=idx, rhs=str(rhs))
                s = s + NextLine
    Code = """subroutine auto2e_RtuvK_{MaxIndex}(RA, RB, n, Fm, X, Y, Z, p)
!
! Hermite Coulomb integrals Rtuv, 0<=t+u+v<={MaxIndex}
!
real(F64), dimension(:), intent(out) :: RA
real(F64), dimension(:), intent(in) :: RB
integer, intent(in) :: n
real(F64), dimension(:), intent(in) :: Fm
real(F64), intent(in) :: X
real(F64), intent(in) :: Y
real(F64), intent(in) :: Z
real(F64), intent(in) :: p
{Recurrence}end subroutine auto2e_RtuvK_{MaxIndex}
""".format(MaxIndex=MaxIndex, Recurrence=s)
    return Code

def RtuvElement(t, u, v, sRtuv, TUVIdx):
    if t < 0 or u < 0 or v < 0:
        return 0
    else:
        Idx = TUVIdx[(t, u, v)]
        Rtuv = Symbol("{R}({I})".format(R=sRtuv, I=Idx))
        return Rtuv
                
def PwrString(Var, Pwr):
    if Pwr == 0:
        return ""
    elif Pwr == 1:
        return "{Var}*".format(Var=Var)
    else:
        return "{V}**{P}*".format(V=Var, P=Pwr)

def MomentumTransfer_Recursion(Prefactor, C, D, lC0):
    """
Recursive implementation of Eq. 18 in Head-Gordon, M. and Pople, J.A. J. Chem. Phys. 89, 5777 (1988);
doi: 10.1063/1.455553. Applied for both bra and ket momentum transfer.
    """
    lD = sum(D)
    if lD > 0:
        if D[0] > 0:
            Delta = array([1,0,0])
        elif D[1] > 0:
            Delta = array([0,1,0])
        elif D[2] > 0:
            Delta = array([0,0,1])
        Sum1 = MomentumTransfer_Recursion(Prefactor, C+Delta, D-Delta, lC0)
        Sum2 = MomentumTransfer_Recursion(Prefactor+Delta, C, D-Delta, lC0)
        #
        # Sum1 <- Sum1 + Sum2
        #
        for key2, val2 in Sum2.items():
            if key2 in Sum1:
                P1, I1 = Sum1[key2]
                P2, I2 = val2
                if I1 != I2:
                    #
                    # That shouldn't happen if the algorithm works as intended
                    #
                    print("Error: the same prefactor assigned to two different T/W matrix elements")
                    sys.exit()
                else:
                    Sum1[key2] = (P1+P2, I1)
            else:
                Sum1[key2] = val2
        return Sum1
    else:
        TIdx = WIndex(C, lC0)
        lx = int(Prefactor[0])
        ly = int(Prefactor[1])
        lz = int(Prefactor[2])
        Term = {}
        Term[(lx, ly, lz)] = (1, TIdx)
        return Term

def MomentumTransfer_IsContiguous(Sum0, Sum1, VectorLength):
    Contiguous = True
    if len(Sum0) == len(Sum1):
        for key0 in Sum0.keys():
            if key0 in Sum1:
                p0, i0 = Sum0[key0]
                p1, i1 = Sum1[key0]
                if p0 != p1 or i0 + VectorLength - 1 != i1:
                    Contiguous = False
            else:
                Contiguous = False
    else:
        Contiguous = False
    return Contiguous

def MomentumTransfer_Print(Sum, VectorLength, v0, BraTransfer):
    """
Print the code for the horizontal recurrence relation. The printed code corresponds
to the complete set of steps going from |A+B,S) to |AB). If VectorLength>1, the arithmetic
operations are done for a range of indices all at once using the array slice notation.
    """
    if BraTransfer:
        Coeffs = ["ABx", "ABy", "ABz"]
        Output = "G({Index})"
        Input = "T({Index})"
    else:
        Coeffs = ["CDx", "CDy", "CDz"]
        Output = "T(:, {Index})"
        Input = "W(:, {Index})"
    
    if VectorLength > 1:
        i = "{v0}:{v1}".format(v0=v0, v1=v0+VectorLength-1)
    else:
        i = str(v0)
    s = Output.format(Index=i) + " = "
    FirstTerm = True
    for key, val in sorted(Sum.items()):
        lx, ly, lz = key
        Pref, TIdx = val
        a = PwrString(Coeffs[0], lx)
        b = PwrString(Coeffs[1], ly)
        c = PwrString(Coeffs[2], lz)
        if Pref == 1:
            p = ""
        else:
            p = "{p}*".format(p=Pref)
        t0 = TIdx
        t1 = TIdx + VectorLength - 1
        if VectorLength == 1:
            i = str(t0)
        else:
            i = "{i0}:{i1}".format(i0=t0, i1=t1)
        InputArray = Input.format(Index=i)
        if FirstTerm:
            sign = ""
        else:
            sign = " + "
        s += sign + p + a + b + c + InputArray
        FirstTerm = False
    return SplitLine(s, ["+"])
            
def KetTransfer_CDLoop(lC, lD):
    Body = ""
    VectorLength = 0
    for D in AngularFunctions[lD]:
        for C in AngularFunctions[lC]:
            if VectorLength == 0:
                Sum0 = MomentumTransfer_Recursion(array([0,0,0]), C, D, lC)
                VectorLength = 1
                C0 = C
                D0 = D
            else:
                Sum1 = MomentumTransfer_Recursion(array([0,0,0]), C, D, lC)
                Contiguous = MomentumTransfer_IsContiguous(Sum0, Sum1, VectorLength+1)
                if Contiguous:
                    VectorLength += 1
                else:
                    Body += MomentumTransfer_Print(Sum0, VectorLength, GabIndex(C0, D0), False) + "\n"
                    Sum0 = Sum1
                    C0 = C
                    D0 = D
                    VectorLength = 1
    Body += MomentumTransfer_Print(Sum0, VectorLength, GabIndex(C0, D0), False) + "\n"
    return Body

def KetTransfer_Subroutine(lC, lD):
    Code = """
subroutine auto2e_KetTransfer_{lc}_{ld}(T, W, CDx, CDy, CDz)
real(F64), dimension(:, :), intent(out) :: T
real(F64), dimension(:, :), intent(in) :: W
real(F64), intent(in) :: CDx, CDy, CDz
{KetTransfer}end subroutine auto2e_KetTransfer_{lc}_{ld}
""".format(KetTransfer=KetTransfer_CDLoop(lC, lD), lc=lC, ld=lD)
    return Code    

def KetTransfer_Module(MaxL):
    """ Module for subroutines performing the horizontal recursion, i.e., the momentum transfer,
    for ket pair of oribtals: (..|C+D,S) -> (..|CD). Memory layout: PhiC changes faster than PhiD.

    The subroutine for lC>=0, lD=0 is trivial and, therefore, it's not included in the module.
    """
    Subroutines = ""
    for lD in range(1, MaxL+1):
        for lC in range(lD, MaxL+1):
            Subroutines += KetTransfer_Subroutine(lC, lD)
    Module = """!
! Transform (..|C+D,S) to (..|CD) using the horizontal recurrence relation
! of Eq. 18 in Head-Gordon, M. and Pople, J.A. J. Chem. Phys. 89, 5777 (1988);
! doi: 10.1063/1.455553
!
! Code generated automatically.
!
module auto2e_KetTransfer
use arithmetic
implicit none
contains{Subroutines}end module auto2e_KetTransfer
""".format(Subroutines=Subroutines)
    return Module    

def BraKetTransfer_Subroutine(lA, lB, lC, lD, SpherTransf):
    """ Carry out the horizontal momentum transfer for  bra and ket orbital pairs.
    If requested, combine the momentum transfer with the Cartesian->spherical orbital
    transformation so that the output is in the real spherical GTOs basis.

    The list Order returned by this method contains information on the memory alignment
    of the array of two-electron integrals right after momentum transfers and, possibly,
    the spherical transformation. (The leftmost index changes the fastest, the rightmost
    index changes the slowest.) Due to the way the spherical transformation works,
    the order of the output may not be the same depending on whether SpherTransf
    is requested.
    """

    ApplySpherTransf = (max(lA, lB, lC, lD) >= SPHER_TRANSF_LMIN and SpherTransf)
    Dict = {"lA":lA, "lB":lB, "lC":lC, "lD":lD}
    Dict["Comment"] =  """!
! Transform (LS|KS), L={l0}..{l1}, K={k0}..{k1} using the horizontal recurrence relation
! of Eq. 18 in Head-Gordon, M. and Pople, J.A. J. Chem. Phys. 89, 5777 (1988);
! doi: 10.1063/1.455553
!
! {OutputBasis}
!""".format(l0=lA, l1=lA+lB, k0=lC, k1=lC+lD,
            OutputBasis="The output integrals are expressed in the real spherical harmonics basis."
            if ApplySpherTransf else "The output integrals are expressed in the Cartesian GTOs basis.")
    Dict["Wab"] = WMatrixDims(lA, lB, lC, lD)[0]
    Dict["Wcd"] = WMatrixDims(lA, lB, lC, lD)[1]
    Dict["NabCart"] = NFuncCart(lA) * NFuncCart(lB)
    Dict["NcdCart"] = NFuncCart(lC) * NFuncCart(lD)
    Dict["NabSpher"] = NFuncSpher(lA) * NFuncSpher(lB)
    Dict["NcdSpher"] = NFuncSpher(lC) * NFuncSpher(lD)
    ABSpherTransf = (max(lA, lB) >= SPHER_TRANSF_LMIN and SpherTransf)
    CDSpherTransf = (max(lC, lD) >= SPHER_TRANSF_LMIN and SpherTransf)
    
    if lB > 0 and lD > 0:
        return BraKetTransfer_Subroutine_LLLL(ABSpherTransf, CDSpherTransf, Dict)
    elif lB > 0:
        return BraKetTransfer_Subroutine_LLLS(ABSpherTransf, CDSpherTransf, Dict)
    elif lD > 0:
        return BraKetTransfer_Subroutine_LSLL(ABSpherTransf, CDSpherTransf, Dict)
    else: # lB == 0 and lD == 0
        return BraKetTransfer_Subroutine_LSLS(ABSpherTransf, CDSpherTransf, Dict)

    
def BraKetTransfer_Subroutine_LLLL(ABSpherTransf, CDSpherTransf, Dict):
    Order = ["c", "d", "a", "b"]
    if CDSpherTransf and ABSpherTransf:
        Code = """
subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}(G, W, Ra, Rb, Rc, Rd)
{Comment}
real(F64), dimension({NcdSpher}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(3) :: AB, CD
real(F64), dimension({Wab}, {NcdCart}) :: Txyz
real(F64), dimension({Wab}, {NcdSpher}) :: T
real(F64), dimension({NcdSpher}, {Wab}) :: U
real(F64), dimension({NcdSpher}, {NabCart}) :: Gxyz
AB = Ra - Rb
CD = Rc - Rd
call auto2e_KetTransfer_{lC}_{lD}(Txyz, W, CD(1), CD(2), CD(3))
call auto2e_SpherTransf_U_Vector_{lC}_{lD}(T, Txyz, {Wab})
U = transpose(T)
call auto2e_KetTransfer_{lA}_{lB}(Gxyz, U, AB(1), AB(2), AB(3))
call auto2e_SpherTransf_U_Vector_{lA}_{lB}(G(:, 1:{NabSpher}), Gxyz, {NcdSpher})
end subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
    elif CDSpherTransf and not ABSpherTransf:
        Code = """
subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}(G, W, Ra, Rb, Rc, Rd)
{Comment}
real(F64), dimension({NcdSpher}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(3) :: AB, CD
real(F64), dimension({Wab}, {NcdSpher}) :: T
real(F64), dimension({Wab}, {NcdCart}) :: Txyz
real(F64), dimension({NcdSpher}, {Wab}) :: U
AB = Ra - Rb
CD = Rc - Rd
call auto2e_KetTransfer_{lC}_{lD}(Txyz, W, CD(1), CD(2), CD(3))
call auto2e_SpherTransf_U_Vector_{lC}_{lD}(T, Txyz, {Wab})
U = transpose(T)
call auto2e_KetTransfer_{lA}_{lB}(G(:, 1:{NabCart}), U, AB(1), AB(2), AB(3))
end subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
    elif not CDSpherTransf and ABSpherTransf:
        Code = """
subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}(G, W, Ra, Rb, Rc, Rd)
{Comment}
real(F64), dimension({NcdCart}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(3) :: AB, CD
real(F64), dimension({Wab}, {NcdCart}) :: T
real(F64), dimension({NcdCart}, {Wab}) :: U
real(F64), dimension({NcdCart}, {NabCart}) :: Gxyz
AB = Ra - Rb
CD = Rc - Rd
call auto2e_KetTransfer_{lC}_{lD}(T, W, CD(1), CD(2), CD(3))
U = transpose(T)
call auto2e_KetTransfer_{lA}_{lB}(Gxyz, U, AB(1), AB(2), AB(3))
call auto2e_SpherTransf_U_Vector_{lA}_{lB}(G(:, 1:{NabSpher}), Gxyz, {NcdCart})
end subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
    else:
        Code = """
subroutine auto2e_BraKetTransfer_{lA}_{lB}_{lC}_{lD}(G, W, Ra, Rb, Rc, Rd)
{Comment}
real(F64), dimension({NcdCart}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(3) :: AB, CD
real(F64), dimension({Wab}, {NcdCart}) :: T
real(F64), dimension({NcdCart}, {Wab}) :: U
AB = Ra - Rb
CD = Rc - Rd
call auto2e_KetTransfer_{lC}_{lD}(T, W, CD(1), CD(2), CD(3))
U = transpose(T)
call auto2e_KetTransfer_{lA}_{lB}(G(:, 1:{NabCart}), U, AB(1), AB(2), AB(3))
end subroutine auto2e_BraKetTransfer_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
    return Code, Order

def BraKetTransfer_Subroutine_LSLL(ABSpherTransf, CDSpherTransf, Dict):
    if CDSpherTransf and ABSpherTransf:
        Code = """
subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}(G, W, Rc, Rd)
{Comment}
real(F64), dimension({NcdSpher}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
real(F64), dimension(3), intent(in) :: Rc, Rd
real(F64), dimension(3) :: CD
real(F64), dimension({NabCart}, {NcdSpher}) :: T
real(F64), dimension({NabCart}, {NcdCart}) :: Txyz
real(F64), dimension({NcdSpher}, {NabCart}) :: Gxyz
CD = Rc - Rd
call auto2e_KetTransfer_{lC}_{lD}(Txyz, W, CD(1), CD(2), CD(3))
call auto2e_SpherTransf_U_Vector_{lC}_{lD}(T, Txyz, {NabCart})
Gxyz = transpose(T)
call auto2e_SpherTransf_U_Vector_{lA}_{lB}(G(:, 1:{NabSpher}), Gxyz, {NcdSpher})
end subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
        Order = ["c", "d", "a", "b"]
    elif CDSpherTransf and not ABSpherTransf:
        Code = """
subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}(G, W, Rc, Rd)
{Comment}
real(F64), dimension({NabCart}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
real(F64), dimension(3), intent(in) :: Rc, Rd
real(F64), dimension(3) :: CD
real(F64), dimension({NabCart}, {NcdCart}) :: Gxyz
CD = Rc - Rd
call auto2e_KetTransfer_{lC}_{lD}(Gxyz, W, CD(1), CD(2), CD(3))
call auto2e_SpherTransf_U_Vector_{lC}_{lD}(G(:, 1:{NcdSpher}), Gxyz, {NabCart})
end subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
        Order = ["a", "b", "c", "d"]
    elif not CDSpherTransf and ABSpherTransf:
        Code = """
subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}(G, W, Rc, Rd)
{Comment}
real(F64), dimension({NcdCart}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
real(F64), dimension(3), intent(in) :: Rc, Rd
real(F64), dimension(3) :: CD
real(F64), dimension({NabCart}, {NcdCart}) :: T
real(F64), dimension({NcdCart}, {NabCart}) :: Gxyz
CD = Rc - Rd
call auto2e_KetTransfer_{lC}_{lD}(T, W, CD(1), CD(2), CD(3))
Gxyz = transpose(T)
call auto2e_SpherTransf_U_Vector_{lA}_{lB}(G(:, 1:{NabSpher}), Gxyz, {NcdCart})
end subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
        Order = ["c", "d", "a", "b"]
    else:
        Code = """
subroutine auto2e_BraKetTransfer_{lA}_{lB}_{lC}_{lD}(G, W, Rc, Rd)
{Comment}
real(F64), dimension({NabCart}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
real(F64), dimension(3), intent(in) :: Rc, Rd
real(F64), dimension(3) :: CD
CD = Rc - Rd
call auto2e_KetTransfer_{lC}_{lD}(G(:, 1:{NcdCart}), W, CD(1), CD(2), CD(3))
end subroutine auto2e_BraKetTransfer_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
        Order = ["a", "b", "c", "d"]
    return Code, Order

def BraKetTransfer_Subroutine_LLLS(ABSpherTransf, CDSpherTransf, Dict):
    Order = ["c", "d", "a", "b"]
    if CDSpherTransf and ABSpherTransf:
        Code = """
subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}(G, W, Ra, Rb)
{Comment}
real(F64), dimension({NcdSpher}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb
real(F64), dimension(3) :: AB
real(F64), dimension({Wab}, {NcdSpher}) :: T
real(F64), dimension({NcdSpher}, {Wab}) :: U
real(F64), dimension({NcdSpher}, {NabCart}) :: Gxyz
AB = Ra - Rb
call auto2e_SpherTransf_U_Vector_{lC}_{lD}(T, W, {Wab})
U = transpose(T)
call auto2e_KetTransfer_{lA}_{lB}(Gxyz, U, AB(1), AB(2), AB(3))
call auto2e_SpherTransf_U_Vector_{lA}_{lB}(G(:, 1:{NabSpher}), Gxyz, {NcdSpher})
end subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
    elif CDSpherTransf and not ABSpherTransf:
        Code = """
subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}(G, W, Ra, Rb)
{Comment}
real(F64), dimension({NcdSpher}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb
real(F64), dimension(3) :: AB
real(F64), dimension({Wab}, {NcdSpher}) :: T
real(F64), dimension({NcdSpher}, {Wab}) :: Uxyz
AB = Ra - Rb
call auto2e_SpherTransf_U_Vector_{lC}_{lD}(T, W, {Wab})
Uxyz = transpose(T)
call auto2e_KetTransfer_{lA}_{lB}(G(:, 1:{NabCart}), Uxyz, AB(1), AB(2), AB(3))
end subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
    elif not CDSpherTransf and ABSpherTransf:
        Code = """
subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}(G, W, Ra, Rb)
{Comment}
real(F64), dimension({NcdCart}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb
real(F64), dimension(3) :: AB
real(F64), dimension({NcdCart}, {Wab}) :: U
real(F64), dimension({NcdCart}, {NabCart}) :: Gxyz
AB = Ra - Rb
U = transpose(W)
call auto2e_KetTransfer_{lA}_{lB}(Gxyz, U, AB(1), AB(2), AB(3))
call auto2e_SpherTransf_U_Vector_{lA}_{lB}(G(:, 1:{NabSpher}), Gxyz, {NcdSpher})
end subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
    else:
         Code = """
subroutine auto2e_BraKetTransfer_{lA}_{lB}_{lC}_{lD}(G, W, Ra, Rb)
{Comment}
real(F64), dimension({NcdCart}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb
real(F64), dimension(3) :: AB
real(F64), dimension({NcdCart}, {Wab}) :: U
AB = Ra - Rb
U = transpose(W)
call auto2e_KetTransfer_{lA}_{lB}(G(:, 1:{NabCart}), U, AB(1), AB(2), AB(3))
end subroutine auto2e_BraKetTransfer_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
    return Code, Order

def BraKetTransfer_Subroutine_LSLS(ABSpherTransf, CDSpherTransf, Dict):
    if CDSpherTransf and ABSpherTransf:
        Code = """
subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}(G, W)
{Comment}
real(F64), dimension({NcdSpher}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
real(F64), dimension({NabCart}, {NcdSpher}) :: T
real(F64), dimension({NcdSpher}, {NabCart}) :: Gxyz
call auto2e_SpherTransf_U_Vector_{lC}_{lD}(T, W, {NabCart})
Gxyz = transpose(T)
call auto2e_SpherTransf_U_Vector_{lA}_{lB}(G(:, 1:{NabSpher}), Gxyz, {NcdSpher})
end subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
        Order = ["c", "d", "a", "b"]
    elif CDSpherTransf and not ABSpherTransf:
        Code = """
subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}(G, W)
{Comment}
real(F64), dimension({NabCart}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
call auto2e_SpherTransf_U_Vector_{lC}_{lD}(G(:, 1:{NcdSpher}), W, {NabCart})
end subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
        Order = ["a", "b", "c", "d"]
    elif not CDSpherTransf and ABSpherTransf:
        Code = """
subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}(G, W)
{Comment}
real(F64), dimension({NcdCart}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
real(F64), dimension({NcdCart}, {NabCart}) :: Gxyz
Gxyz = transpose(W)
call auto2e_SpherTransf_U_Vector_{lA}_{lB}(G(:, 1:{NabSpher}), Gxyz, {NcdCart})
end subroutine auto2e_BraKetTransfer_Spher_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
        Order = ["c", "d", "a", "b"]
    else:
        Code = """
subroutine auto2e_BraKetTransfer_{lA}_{lB}_{lC}_{lD}(G, W)
{Comment}
real(F64), dimension({NabCart}, *), intent(out) :: G
real(F64), dimension({Wab}, {Wcd}), intent(in) :: W
G(:, 1:{NcdCart}) = W
end subroutine auto2e_BraKetTransfer_{lA}_{lB}_{lC}_{lD}
""".format(**Dict)
        Order = ["a", "b", "c", "d"]
    return Code, Order
    
def LargestCommonFactor(EFactors, k0, k1):
    """ Compute the largest common factor of terms E(K0)....E(K1). """
    if not k1 > k0:
        CommonFactor = set([])
    else:
        CommonFactor = EFactors[k0]
        for k in range(k0+1, k1+1):
            CommonFactor &= EFactors[k]
    return CommonFactor
            
def PartialSum(EFactors, RFactors, Signs, CommonFactor, k0, k1):
    if k0 == 0:
        FirstTerm = True
    else:
        FirstTerm = False
        
    if len(CommonFactor) == 2:
        Ea, Eb = list(CommonFactor)
        if FirstTerm:
            s = "{A}*{B}*(".format(A=Ea, B=Eb)
        else:
            s = " + {A}*{B}*(".format(A=Ea, B=Eb)
        #
        # First term in the parentheses
        #
        Ec = list(EFactors[k0]-CommonFactor)[0]
        if Signs[k0] == "+":
            s += "{C}*{R}".format(C=Ec, R=RFactors[k0])
        else:
            s += "-{C}*{R}".format(C=Ec, R=RFactors[k0])
        #
        # Consecutive terms
        #
        for k in range(k0+1, k1+1):
            Ec = list(EFactors[k]-CommonFactor)[0]
            s += " {Sign} {C}*{R}".format(C=Ec, R=RFactors[k], Sign=Signs[k])
        s += ")"
    else:
        #
        # First term in the partial sum
        #
        Ea, Eb, Ec = list(EFactors[k0])
        if FirstTerm:
            if Signs[k0] == "+":
                s = "{A}*{B}*{C}*{R}".format(A=Ea, B=Eb, C=Ec, R=RFactors[k0])
            else:
                s = "-{A}*{B}*{C}*{R}".format(A=Ea, B=Eb, C=Ec, R=RFactors[k0])
        else:
            s = "{Sign} {A}*{B}*{C}*{R}".format(A=Ea, B=Eb, C=Ec, R=RFactors[k0], Sign=Signs[k0])
        #
        # Consecutive terms
        #
        for k in range(k0+1, k1+1):
            Ea, Eb, Ec = list(EFactors[k])
            s += " {Sign} {A}*{B}*{C}*{R}".format(A=Ea, B=Eb, C=Ec, R=RFactors[k], Sign=Signs[k])
    return s

def HermiteToCartesian_Transform_RHS(k0, EFactors, RFactors, Signs):
    if k0 >= len(EFactors):
        return ""
    #
    # Determine how many terms share the same CommonFactor of length 2
    #
    CommonFactor = EFactors[k0]
    k1 = k0
    for k in range(k0+1, len(EFactors)):
        NewCommonFactor = CommonFactor & EFactors[k]
        if len(NewCommonFactor) > 2:
            print("Invalid Ex*Ex*Ez coefficients")
            sys.exit(1)
        if len(NewCommonFactor) == 2:
            CommonFactor = NewCommonFactor
            k1 += 1
        else:
            break
    if k0 == k1:
        CommonFactor = set([])
    s = PartialSum(EFactors, RFactors, Signs, CommonFactor, k0, k1)
    return s + HermiteToCartesian_Transform_RHS(k1+1, EFactors, RFactors, Signs)

def HermiteToCartesian_KetTransform_SingleTerm(t, u, v, lx, ly, lz, LBra):
    """ Transform ket orbital pair from the Hermite to Cartesian basis 
    (Eq. 9.9.47 in Helgaker's textbook. The pahse factor (-1)**(tau+nu+phi)
    is included. """
    
    LKet = lx + ly + lz
    IdxBra = tuvindices(LBra)
    IdxTUV = tuvindices(LBra+LKet)
    if t[0] == t[1]:
        s = "S({I}) = ".format(I=IdxBra[(t[0], u, v)])
    else:
        s = "S({I0}:{I1}) = ".format(I0=IdxBra[(t[0], u, v)], I1=IdxBra[(t[1], u, v)])

    Signs = []
    EFactors = []
    RFactors = [] 
    for phi in range(lz+1):
        for nu in range(ly+1):
            for tau in range(lx+1):
                Ex = "Ex({i})".format(i=EtIdx[(lx, tau)])
                Ey = "Ey({i})".format(i=EtIdx[(ly, nu)])
                Ez = "Ez({i})".format(i=EtIdx[(lz, phi)])
                if t[0] == t[1]:
                    R = "R({i})".format(i=IdxTUV[(t[0]+tau, u+nu, v+phi)])
                else:
                    i0 = IdxTUV[(t[0]+tau, u+nu, v+phi)]
                    i1 = IdxTUV[(t[1]+tau, u+nu, v+phi)]
                    R = "R({i0}:{i1})".format(i0=i0, i1=i1)
                if (tau+nu+phi) % 2 == 0:
                    Sign = "+"
                else:
                    Sign = "-"
                EFactors.append(set([Ex, Ey, Ez]))
                RFactors.append(R)
                Signs.append(Sign)
                
    s += HermiteToCartesian_Transform_RHS(0, EFactors, RFactors, Signs) + "\n"
    return s

def RCopy_KetTransform_Subroutine(LBra, LKet):
    """ Hermite -> Cartesian ket transform using a contiguous copy of the Rtuv matrix. """
    Code = """
subroutine auto2e_RCopy_KetTransform_{LBra}_{LKet}(S, T, Ex, Ey, Ez, R, lx, ly, lz)
!
! Transform the ket shell pair from Hermite to Cartesian Gaussian basis.
! This variant of the transformation algorithm starts by copying the Rtuv
! matrix elements into contiguous memory locations.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(out) :: T
real(F64), dimension(:), intent(in) :: Ex, Ey, Ez
real(F64), dimension(:), intent(in) :: R
integer, intent(in) :: lx, ly, lz
real(F64) :: c
integer :: tau, nu, phi, i
integer :: x0, y0, z0, x, y, z
S = ZERO
x0 = ((lx + 1) * lx) / 2 + 1
y0 = ((ly + 1) * ly) / 2 + 1
z0 = ((lz + 1) * lz) / 2 + 1
do phi = 0, lz
do nu = 0, ly
do tau = 0, lx
i = ((2*{LKet}+3-phi)*phi)/2+nu+1
call auto2e_Rcopy_{LBra}_{LKet}(T, R(tau+1:), RCopyIdx_{LBra}_{LKet}(:, i))
x = x0 + tau
y = y0 + nu
z = z0 + phi
c = (-1)**modulo(tau+nu+phi, 2)*Ex(x)*Ey(y)*Ez(z)
S = S + c * T
end do
end do
end do
end subroutine auto2e_RCopy_KetTransform_{LBra}_{LKet}
""".format(LBra=LBra, LKet=LKet)
    return Code

def RCopyIdx(LBra, LKet):
    RIdx = tuvindices(LBra+LKet)
    n = 0
    indices = []
    for phi in range(LKet+1):
        for nu in range(LKet-phi+1):
            n += 1
            m = 0
            for v in range(LBra+1):
                for u in range(LBra+1-v):
                    m += 1
                    indices.append(str(RIdx[(0, u+nu, v+phi)]))
                    
    s = SplitLine(("integer, dimension({m}, {n}), parameter :: RCopyIdx_{LBra}_{LKet} = reshape([{i}], [{m}, {n}])"
         .format(m=m, n=n,  i=", ".join(indices), LBra=LBra, LKet=LKet)), [" "])
    return s

def RCopy_Subroutine(LBra, LKet):
    SIdx = tuvindices(LBra)
    loop = ""
    k = 1
    for v in range(LBra+1):
        for u in range(LBra+1-v):
            MaxT = LBra - u - v
            i0 = SIdx[(0, u, v)]
            i1 = SIdx[(MaxT, u, v)]
            if MaxT >= 1:
                loop += "T({i0}:{i1}) = R(idx({k}):idx({k})+{MaxT})\n".format(i0=i0, i1=i1, k=k, MaxT=MaxT)
            else:
                loop += "T({i0}) = R(idx({k}))\n".format(i0=i0, k=k)
            k += 1
    Subroutine = """
subroutine auto2e_RCopy_{LBra}_{LKet}(T, R, idx)
real(F64), dimension(:), intent(out) :: T
real(F64), dimension(:), intent(in) :: R
integer, dimension(:), intent(in) :: idx
{loop}end subroutine auto2e_RCopy_{LBra}_{LKet}
""".format(loop=loop, LBra=LBra, LKet=LKet)
    return Subroutine

def RCopy_KetTransform_Module(LBra, LKet):
    ModuleName = "auto2e_KetTransform_{LBra}_{LKet}".format(LBra=LBra, LKet=LKet)
    Code = []
    Code.append((ModuleName, ("""module {ModuleName}
use arithmetic
implicit none

{IdxArray}

contains
{RCopy}{KetTransform}end module {ModuleName}"""
                              .format(RCopy=RCopy_Subroutine(LBra, LKet),
                                      KetTransform=RCopy_KetTransform_Subroutine(LBra, LKet),
                                      ModuleName=ModuleName,
                                      IdxArray=RCopyIdx(LBra, LKet)))))
    return Code

def HermiteToCartesian_KetTransformSubroutine(LBra, lx, ly, lz):
    KetTransform = ""
    for v in range(LBra+1):
        for u in range(LBra-v+1):
            #
            # Transformation is done in a vectorized way, for all values of t = 0, ..., LBra-v-u+1
            #
            t0 = 0
            t1 = LBra-v-u
            NewLine = HermiteToCartesian_KetTransform_SingleTerm([t0,t1], u, v, lx, ly, lz, LBra)
            KetTransform += SplitLine(NewLine, ["+", "-", "*"])
    Code = """
subroutine auto2e_KetTransform_{L}_{x}_{y}_{z}(S, Ex, Ey, Ez, R)
!
! Perform the Hermite -> Cartesian transformation for ket pairs of orbitals:
! S(t,u,v) = Sum(tau,nu,phi) Etau(kx)*Enu(ky)*Ephi(kz) * R(t+tau,nu+u,phi+v)
! 0<=t+u+v<={L}
! Class=(LS|KS), L={L}, K=({x},{y},{z})
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
real(F64), dimension(:), intent(in) :: R
{Transform}end subroutine auto2e_KetTransform_{L}_{x}_{y}_{z}
""".format(L=LBra, x=lx, y=ly, z=lz, Transform=KetTransform)
    return Code

def KetTransform_Module(LBra, LKet):
    MaxNLines = 2000
    Subroutines = [""]
    NModules = 1
    for lx, ly, lz in AngularFunctions[LKet]:
        Subroutines[NModules-1] += HermiteToCartesian_KetTransformSubroutine(LBra, lx, ly, lz)
        NLines = Subroutines[NModules-1].count("\n")
        if NLines > MaxNLines:
            Subroutines.append("")
            NModules += 1
    if Subroutines[NModules-1] == "":
        NModules += -1
    Code = []
    for k in range(NModules):
        if NModules > 1:
            ModuleName = "auto2e_KetTransform_{LBra}_{LKet}_part{k}".format(LBra=LBra, LKet=LKet, k=k+1)
        else:
            ModuleName = "auto2e_KetTransform_{LBra}_{LKet}".format(LBra=LBra, LKet=LKet)
        Code.append((ModuleName, """module {ModuleName}
use arithmetic
implicit none
contains
{Subs}end module {ModuleName}""".format(Subs=Subroutines[k], ModuleName=ModuleName)))
    return Code

def HermiteToCartesian_RCopy_BraKetTransform(LBra, LKet):
    loop = """k = 0
do lx = {LKet}, 0, -1
do ly = {LKet}-lx, 0, -1
lz = {LKet} - lx - ly
call auto2e_RCopy_KetTransform_{LBra}_{LKet}(S, T, ExCD, EyCD, EzCD, R, lx, ly, lz)
call auto2e_BraTransform_{LBra}(W(Row0:, Col0+k), S, ExAB, EyAB, EzAB)
k = k + 1
end do
end do
""".format(LBra=LBra, LKet=LKet)
    return loop

def HermiteToCartesian_BraKetTransform(LBra, LKet):
    """ Code for the Hermite -> Cartesian transform of two-electron integrals
    belonging to the (LBra S|LKet S) class.
    
    Special-case code is applied for S angular momenta.
    """
    s = ""
    if LKet == 0:
        s += "S = ExCD(1)*EyCD(1)*EzCD(1)*R\n"
        if LBra == 0:
            s += "W(Row0, Col0) = W(Row0, Col0) + ExAB(1)*EyAB(1)*EzAB(1)*S(1)\n"
        elif LBra == 1:
            s += """W(Row0+0, Col0) = W(Row0+0, Col0) + EyAB(1)*EzAB(1)*(ExAB(2)*S(1) + ExAB(3)*S(2))
W(Row0+1, Col0) = W(Row0+1, Col0) + ExAB(1)*EzAB(1)*(EyAB(2)*S(1) + EyAB(3)*S(3))
W(Row0+2, Col0) = W(Row0+2, Col0) + EyAB(1)*ExAB(1)*(EzAB(2)*S(1) + EzAB(3)*S(4))
"""
        else:
            s += "call auto2e_BraTransform_{L}(W(Row0:, Col0), S, ExAB, EyAB, EzAB)\n".format(L=LBra)
    else:
        for lx, ly, lz in AngularFunctions[LKet]:
            if LBra == 1 and LKet == 1:
                if (lx, ly, lz) == (1, 0, 0):
                    s += """S(1:2) = EyCD(1)*EzCD(1)*(ExCD(2)*R(1:2) - ExCD(3)*R(2:3))
S(3) = EyCD(1)*EzCD(1)*(ExCD(2)*R(4) - ExCD(3)*R(5))
S(4) = EyCD(1)*EzCD(1)*(ExCD(2)*R(7) - ExCD(3)*R(8))
"""
                elif (lx, ly, lz) == (0, 1, 0):
                    s += """S(1:2) = ExCD(1)*EzCD(1)*(EyCD(2)*R(1:2) - EyCD(3)*R(4:5))
S(3) = ExCD(1)*EzCD(1)*(EyCD(2)*R(4) - EyCD(3)*R(6))
S(4) = ExCD(1)*EzCD(1)*(EyCD(2)*R(7) - EyCD(3)*R(9))
"""
                elif (lx, ly, lz) == (0, 0, 1):
                    s += """S(1:2) = EyCD(1)*ExCD(1)*(EzCD(2)*R(1:2) - EzCD(3)*R(7:8))
S(3) = EyCD(1)*ExCD(1)*(EzCD(2)*R(4) - EzCD(3)*R(9))
S(4) = EyCD(1)*ExCD(1)*(EzCD(2)*R(7) - EzCD(3)*R(10))
"""
            else:
                s += "call auto2e_KetTransform_{L}_{x}_{y}_{z}(S, ExCD, EyCD, EzCD, R)\n".format(L=LBra, x=lx, y=ly, z=lz)
            if LBra == 0:
                s += "W(Row0, Col0+{I}) = W(Row0, Col0+{I}) + ExAB(1)*EyAB(1)*EzAB(1)*S(1)\n".format(I=XYZIdx[(lx,ly,lz)]-1)
            elif LBra == 1:
                s += """W(Row0+0, Col0+{I}) = W(Row0+0, Col0+{I}) + EyAB(1)*EzAB(1)*(ExAB(2)*S(1) + ExAB(3)*S(2))
W(Row0+1, Col0+{I}) = W(Row0+1, Col0+{I}) + ExAB(1)*EzAB(1)*(EyAB(2)*S(1) + EyAB(3)*S(3))
W(Row0+2, Col0+{I}) = W(Row0+2, Col0+{I}) + EyAB(1)*ExAB(1)*(EzAB(2)*S(1) + EzAB(3)*S(4))
""".format(I=XYZIdx[(lx,ly,lz)]-1)
            else:
                s += "call auto2e_BraTransform_{L}(W(Row0:, Col0+{I}), S, ExAB, EyAB, EzAB)\n".format(L=LBra, I=XYZIdx[(lx,ly,lz)]-1)
    return s

def HermiteToCartesian_BraTransform_SingleTerm(lx, ly, lz):
    """ Transform bra orbital pair from the Hermite to Cartesian basis 
    (Eq. 9.9.48 in Helgaker's textbook. """
    
    LBra = lx + ly + lz
    IdxBra = tuvindices(LBra)
    s = "W({i}) = W({i}) + ".format(i=XYZIdx[(lx,ly,lz)])
    EFactors = []
    SFactors = []
    Signs = []
    for v in range(lz+1):
        for u in range(ly+1):
            for t in range(lx+1):
                Ex = "Ex({i})".format(i=EtIdx[(lx, t)])
                Ey = "Ey({i})".format(i=EtIdx[(ly, u)])
                Ez = "Ez({i})".format(i=EtIdx[(lz, v)])
                S = "S({i})".format(i=IdxBra[(t, u, v)])
                EFactors.append(set([Ex, Ey, Ez]))
                SFactors.append(S)
                Signs.append("+")
    s += HermiteToCartesian_Transform_RHS(0, EFactors, SFactors, Signs) + "\n"
    return SplitLine(s, ["+"])

def HermiteToCartesian_BraTransformSubroutine(LBra):
    """ Carry out the Hermite -> Cartesian transform for all angular functions
    of the bra shell pair LS, L=(lx,ly,lz), S=(0,0,0). """

    BraTransform = ""
    for lx, ly, lz in AngularFunctions[LBra]:
        BraTransform += HermiteToCartesian_BraTransform_SingleTerm(lx, ly, lz)
    Code = """
subroutine auto2e_BraTransform_{L}(W, S, Ex, Ey, Ez)
!
! Perform the Hermite -> Cartesian transformation for bra pairs of orbitals.
! Class=(LS|..), L={L}
! Code generated automatically.
!
real(F64), dimension(:), intent(inout) :: W
real(F64), dimension(:), intent(in) :: S
real(F64), dimension(:), intent(in) :: Ex
real(F64), dimension(:), intent(in) :: Ey
real(F64), dimension(:), intent(in) :: Ez
{Transform}end subroutine auto2e_BraTransform_{L}
""".format(L=LBra, Transform=BraTransform)
    return Code

def HermiteToCartesian_BraTransformModule(MaxL):
    Subroutines = ""
    for l in range(2, 2*MaxL+1):
        Subroutines += HermiteToCartesian_BraTransformSubroutine(l)
    Code = """!
! Code generated automatically.
!
module auto2e_BraTransform
use arithmetic
implicit none

contains
{sub}end module auto2e_BraTransform""".format(sub=Subroutines)
    return Code
        
def EtSubroutine(L):
    """ Compute the Hermite -> Cartesian transformation coefficients Et^(l0) for L=0,...,L and t=0...l. """
    s = ""
    for l in range(L+1):
        for t in range(l+1):
            if t == 0:
                if l == 0:
                    s += "E({i}) = E000\n".format(i=EtIdx[(0, 0)])
                elif l == 1:
                    s += "E({i}) = X * E({j})\n".format(i=EtIdx[(l, 0)], j=EtIdx[(l-1, 0)])
                else:
                    s += "E({i}) = X * E({j}) + E({k})\n".format(i=EtIdx[(l, 0)], j=EtIdx[(l-1, 0)], k=EtIdx[(l-1, 1)])
            else:
                s += "E({i}) = {l}.0_F64/({TwoT}.0_F64*p) * E({j})\n".format(i=EtIdx[(l, t)], j=EtIdx[(l-1, t-1)], l=l, TwoT=2*t)

    Code = """subroutine auto2e_EtCoeffs_{L}(E, p, X, E000)
!
! Compute the coefficients (Et) of the Hermite Gaussians -> Cartesian Gaussians transformation
! for the angular momentum pair ({L}, 0).
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: E
real(F64), intent(in) :: p
real(F64), intent(in) :: X
real(F64), intent(in) :: E000
{RecursiveFormulas}end subroutine auto2e_EtCoeffs_{L}
""".format(L=L, RecursiveFormulas=s)
    return Code
    
def RtuvSubroutine(n):
    """ Compute the Hermite Coulomb integrals matrix R(t, u, v) for 0<=t+u+v<=n. """
    Code = ""
    #
    # Interchange "RA" and "RB" so that the last write operation
    # stores the final numbers in "RA", that is, the output of
    # the subroutine
    #
    RecursionCode = ""
    if n % 2 == 0:
        R1, R2 = "RA", "RB"
    else:
        R1, R2 = "RB", "RA"
    for k in range(n+1):
        if k == 0:
            RecursionCode += "{ROut}(1) = (-2*p)**{n} * Fm({m})\n".format(n=n,  m=n+1, ROut=R1)
        else:
            RecursionCode += "call auto2e_RtuvK_{K}({ROut}, {RIn}, {N}, Fm, X, Y, Z, p)\n".format(K=k, ROut=R1, RIn=R2, N=n-k)
        R1, R2 = R2, R1
    if n == 0:
        RBDeclaration = ""
    else:
        RBDeclaration = "real(F64), dimension({RBDim}) :: RB\n".format(RBDim=RMatrixDim(n-1))
    Code = """subroutine auto2e_Rtuv_{n}(RA, Fm, p, X, Y, Z)
!
! Recursively compute the Hermite Coulomb integrals Rtuv for 0<=t+u+v<={n},
! Eqs. 9.9.18-20 p. 375 in Helgaker's textbook.
!
! Code generated automatically.
!
real(F64), dimension(:), intent(out) :: RA
real(F64), dimension(:), intent(in) :: Fm
real(F64), intent(in) :: p
real(F64), intent(in) :: X
real(F64), intent(in) :: Y
real(F64), intent(in) :: Z
{RB}{RecursionCode}end subroutine auto2e_Rtuv_{n}
""".format(n=n, RecursionCode=RecursionCode, RB=RBDeclaration)
    return Code


def HermiteModule(MaxL):
    subs = ""
    for n in range(1, 2*MaxL+1):
        subs += "\n" + EtSubroutine(n)
    for n in range(1, 4*MaxL+1):
        subs += "\n" + RtuvKSubroutine(n)
    for n in range(1, 4*MaxL+1):
        subs += "\n" + RtuvSubroutine(n)
        
    Module = """!
! Code generated automatically.
!
module auto2e_Hermite
use arithmetic
implicit none

contains
{Subroutines}end module auto2e_Hermite
""".format(Subroutines=subs)
    return Module

def BoysFunction(m):
    """ Boys function call. """
    MaxSpecialCases = 6
    if m <= MaxSpecialCases:
        s = "call f{}(Alpha * dot_product(Rpq, Rpq), F)".format(m)
    else:
        s = "call fm({}, Alpha * dot_product(Rpq, Rpq), F)".format(m)
    return s

def WMatrixSubroutineTransposed(LBra, LKet):
    Subroutine = """
subroutine auto2e_W_{la}_{lc}(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64), dimension({nc}, {na}) :: WT
integer :: i, j
WT = ZERO
call auto2e_W_{lc}_{la}(WT, Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, &
Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, Kappa, 1, 1)
do j = 1, {nc}
do i = 1, {na}
W(Row0+i-1, Col0+j-1) = WT(j, i)
end do
end do
end subroutine auto2e_W_{la}_{lc}
""".format(la=LBra, lc=LKet, na=NFuncCart(LBra), nc=NFuncCart(LKet))
    return Subroutine
    
def WMatrixSubroutine(LBra, LKet):
    """ Compute (LS|KS) integrals, L=LBra, K=LKet """
    RDim = RMatrixDim(LBra+LKet)
    SDim = RMatrixDim(LBra)

    if EnableRCopy(LBra, LKet):
        RCopyVars = """real(F64), dimension({SDim}) :: T
integer :: lx, ly, lz, k
""".format(SDim=SDim)
    else:
        RCopyVars = ""
    
    if LKet > 0:
        EtCD = """Rqc = Rq - Rc
call auto2e_EtCoeffs_{lc}(ExCD, Q, Rqc(1), ONE)
call auto2e_EtCoeffs_{lc}(EyCD, Q, Rqc(2), ONE)
call auto2e_EtCoeffs_{lc}(EzCD, Q, Rqc(3), ONE)""".format(lc=LKet)
        Rqc = "Rqc, "
    else:
        EtCD = """ExCD(1) = ONE
EyCD(1) = ONE
EzCD(1) = ONE"""
        Rqc = ""
    if LBra > 0:
        EtAB = """Rpa = Rp - Ra
call auto2e_EtCoeffs_{la}(ExAB, P, Rpa(1), E000)
call auto2e_EtCoeffs_{la}(EyAB, P, Rpa(2), ONE)
call auto2e_EtCoeffs_{la}(EzAB, P, Rpa(3), ONE)""".format(la=LBra)
        Rpa = "Rpa, "
    else:
        EtAB = """ExAB(1) = E000
EyAB(1) = ONE
EzAB(1) = ONE"""
        Rpa = ""
    if LBra+LKet > 0:
        Rtuv = "call auto2e_Rtuv_{}(R, F, Alpha, Rpq(1), Rpq(2), Rpq(3))".format(LBra+LKet)
    else:
        Rtuv = "R(1) = F(1)"
    #
    # Code for the Harmite->Cartesian Gaussian transformation
    #
    if EnableRCopy(LBra, LKet):
        HCTransf = HermiteToCartesian_RCopy_BraKetTransform(LBra, LKet)
    else:
        HCTransf = HermiteToCartesian_BraKetTransform(LBra, LKet)
    Subroutine = """
subroutine auto2e_W_{la}_{lc}(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, Row0, Col0)
real(F64), dimension(:, :), intent(inout) :: W
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
integer, intent(in) :: Row0, Col0
real(F64) :: Rab2, Rcd2, P, Q, MuAB, MuCD, Alpha, N, Kabcd, E000, PQ
real(F64), dimension(3) :: Rp, Rq, {Rqc}{Rpa}Rpq
integer :: Ga, Gb, Gc, Gd
real(F64), dimension({RDim}) :: R
real(F64), dimension({SDim}) :: S
real(F64), dimension({FDim}) :: F
real(F64), dimension({EtABdim}) :: ExAB, EyAB, EzAB
real(F64), dimension({EtCDdim}) :: ExCD, EyCD, EzCD
{RCopyVars}Rab2 = (Ra(1)-Rb(1))**2 + (Ra(2)-Rb(2))**2 + (Ra(3)-Rb(3))**2
Rcd2 = (Rc(1)-Rd(1))**2 + (Rc(2)-Rd(2))**2 + (Rc(3)-Rd(3))**2
do Gd = 1, NprimD
do Gc = 1, NprimC
Q = ExpC(Gc) + ExpD(Gd)
MuCD = ExpC(Gc) * ExpD(Gd) / Q
Rq = (ExpC(Gc) * Rc + ExpD(Gd) * Rd) / Q
{EtCD}
do Gb = 1, NprimB
do Ga = 1, NprimA
P = ExpA(Ga) + ExpB(Gb)
MuAB = ExpA(Ga) * ExpB(Gb) / P
Rp = (ExpA(Ga) * Ra + ExpB(Gb) * Rb) / P
Rpq = Rp - Rq
Kabcd = exp(-MuAB * Rab2 - MuCD * Rcd2)
PQ = P * Q
Alpha = PQ / (P + Q + Kappa*PQ)
N = TWO * PI**(FIVE/TWO) / (PQ * Sqrt(P + Q + Kappa*PQ))
E000 = N * CntrA(Ga) * CntrB(Gb) * CntrC(Gc) * CntrD(Gd) * Kabcd
{EtAB}
{Boys}
{Rtuv}
{HermiteToCartesian}end do
end do
end do
end do
end subroutine auto2e_W_{la}_{lc}
""".format(la=LBra, lc=LKet, RDim=RDim, SDim=SDim, FDim=LBra+LKet+1, HermiteToCartesian=HCTransf,
           EtABdim=EtIdx[(LBra, LBra)], EtCDdim=EtIdx[(LKet, LKet)], EtAB=EtAB, EtCD=EtCD, Rpa=Rpa, Rqc=Rqc,
           Boys=BoysFunction(LBra+LKet), Rtuv=Rtuv, RCopyVars=RCopyVars)
    return Subroutine

def WMatrix_Module(MaxL, ModuleNames):
    Comment = """!
! Compute the electron repulsion integrals (LS|KS), 0<=L,K<={MaxL}, and store them in W(Row0:, Col0:).
! Because this subroutine is used before the horizontal recurrence step, the contraction coefficients
! (CntrA, CntrB, ...) depend on the primitive Gaussian index, but do not depend on the angular function.
!
! The contracted Gaussian functions are defined as
!
! PhiA(lx,ly,lz) = (x-Ra(1))**lx * (y-Ra(2))**ly * (z-Ra(3))**lz * 
!     Sum(Ga) CntrA(Ga) Exp(-ExpA(Ga) * ((x-Ra(1))**2 + (y-Ra(2))**2 + (z-Ra(3))**2))
!
! The order of the angular-dependent prefactors is defined by the following loop:
!
! I = 1
! for lx in L..0
!     for ly in L-lx..0
!         lz = L-lx-ly
!         Index[(lx,ly,lz)] = I
!         I = I + 1
!
! Set Kappa=0 for conventional integrals with the Coulomb operator 1/r12.
! Set Kappa=1/Omega**2 for integrals with the screened interaction potential=Erf(Omega*r12)/r12.
!
! Integrals with the Erf(Omega*r12)/r12 operator
! ----------------------------------------------
! The formulas for Alpha(omega) and N(omega) (see the code) are drived from Eq. 12 in
! Gill, P.M.W, Adamson, R.D. Chem. Phys. Lett. 261, 105 (1996); doi: 10.1016/0009-2614(96)00931-1,
! where Alpha=1/P+1/Q+1/Omega**2. N(omega) has to cancel the 1/Sqrt(Alpha) factor in the Boys function,
! and that's why it depends on Omega.
!
! Code generated automatically.
!
""".format(MaxL=2*MaxL)

    ImportedModules = ""
    for s in sorted(ModuleNames):
        ImportedModules += "use {name}\n".format(name=s)

    MaxNLines = 2000
    Subroutines = [""]
    NModules = 1
    for L in range(2*MaxL+1):
        for K in range(L+1):
            Subroutines[NModules-1] += WMatrixSubroutine(L, K)
            if L != K:
                Subroutines[NModules-1] += WMatrixSubroutineTransposed(K, L)
            NLines = Subroutines[NModules-1].count("\n")
            if NLines > MaxNLines:
                Subroutines.append("")
                NModules += 1
    
    Code = []
    for k in range(NModules):
        if NModules > 1:
            ModuleName = "auto2e_WMatrix_part{k}".format(k=k+1)
        else:
            ModuleName = "auto2e_WMatrix"
        Code.append((ModuleName, """{Comment}module {ModuleName}
use arithmetic
use boys
use auto2e_Hermite
use auto2e_BraTransform
{Import}
implicit none
contains
{Subroutines}end module {ModuleName}
""".format(Subroutines=Subroutines[k], Import=ImportedModules, ModuleName=ModuleName, Comment=Comment if k==0 else "")))
    return Code

def NormalizeTranspose_Subroutine(lA, lB, lC, lD, TargetOrder, InputOrder, SpherTransf):
    label = "".join(TargetOrder).upper()
    if SpherTransf:
        N = {"a":NFuncSpher(lA), "b":NFuncSpher(lB), "c":NFuncSpher(lC), "d":NFuncSpher(lD)}
    else:
        N = {"a":NFuncCart(lA), "b":NFuncCart(lB), "c":NFuncCart(lC), "d":NFuncCart(lD)}
    p, q, r, s = InputOrder
    x0, x1, x2, x3 = TargetOrder
    Gin = "G{}{}{}{}".format(p, q, r, s)
    Gout = "G{}{}{}{}".format(x0, x1, x2, x3)
    Code = """
subroutine auto2e_Normalize{Spher}_{a}_{b}_{c}_{d}_{LABEL}({Gout}, {Gin}, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants and change the memory layout:
!
! Input memory layout: {InputOrder}
! Output memory layout: {OutputOrder}
!
real(F64), dimension(*), intent(out) :: {Gout}
real(F64), dimension(*), intent(in) :: {Gin}
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v, w
v = 1
do {s} = 1, {Ns}
do {r} = 1, {Nr}
do {q} = 1, {Nq}
do {p} = 1, {Np}
w = {x0} + ({x1}-1)*{N1} + ({x2}-1)*{N2} + ({x3}-1)*{N3}
{Gout}(w) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*{Gin}(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize{Spher}_{a}_{b}_{c}_{d}_{LABEL}
""".format(a=lA, b=lB, c=lC, d=lD,
           p=p, q=q, r=r, s=s,
           x0=x0, x1=x1, x2=x2, x3=x3,
           Ns=N[s], Nr=N[r], Nq=N[q], Np=N[p],
           N1=N[x0],
           N2=N[x0]*N[x1],
           N3=N[x0]*N[x1]*N[x2],
           LABEL=label, Gin=Gin, Gout=Gout,
           InputOrder="".join(InputOrder).upper(),
           OutputOrder="".join(TargetOrder).upper(),
           Spher="_Spher" if (NeedsSpherTransf(lA, lB, lC, lD) and SpherTransf) else "")
    return Code

def Normalize_Subroutine(lA, lB, lC, lD, InputOrder, SpherTransf):
    label = "".join(InputOrder).upper()
    if SpherTransf:
        N = {"a":NFuncSpher(lA), "b":NFuncSpher(lB), "c":NFuncSpher(lC), "d":NFuncSpher(lD)}
    else:
        N = {"a":NFuncCart(lA), "b":NFuncCart(lB), "c":NFuncCart(lC), "d":NFuncCart(lD)}
    p, q, r, s = InputOrder
    Code = """
subroutine auto2e_Normalize{Spher}_{a}_{b}_{c}_{d}_{LABEL}(G, NormA, NormB, NormC, NormD)
!
! Apply L-dependent normalization constants
!
real(F64), dimension(*), intent(inout) :: G
real(F64), dimension(*), intent(in) :: NormA, NormB, NormC, NormD
integer :: a, b, c, d, v
v = 1
do {s} = 1, {Ns}
do {r} = 1, {Nr}
do {q} = 1, {Nq}
do {p} = 1, {Np}
G(v) = NormA(a)*NormB(b)*NormC(c)*NormD(d)*G(v)
v = v + 1
end do
end do
end do
end do
end subroutine auto2e_Normalize{Spher}_{a}_{b}_{c}_{d}_{LABEL}
""".format(a=lA, b=lB, c=lC, d=lD,
           p=p, q=q, r=r, s=s,
           Ns=N[s], Nr=N[r], Nq=N[q], Np=N[p],
           LABEL=label,
           Spher="_Spher" if (NeedsSpherTransf(lA, lB, lC, lD) and SpherTransf) else "")
    return Code

def BraKetTransfer_Call(lA, lB, lC, lD, SpherTransf):
    if (lA, lB, lC, lD) == (0, 0, 0, 0):
        s = "G(1) = W(1, 1)"
    else:
        if lA == 0 or lB == 0:
            RAB = ""
        else:
            RAB = ", Ra, Rb"
        if lC == 0 or lD == 0:
            RCD = ""
        else:
            RCD = ", Rc, Rd" 
        s = "call auto2e_BraKetTransfer{Spher}_{la}_{lb}_{lc}_{ld}(G, W{RAB}{RCD})".format(
            RAB=RAB, RCD=RCD, la=lA, lb=lB, lc=lC, ld=lD,
            Spher="_Spher" if (SpherTransf and NeedsSpherTransf(lA, lB, lC, lD)) else "")
    return s

def ERI_Subroutine(lA, lB, lC, lD, MemoryLayout, SpherTransf):
    """ Computation of the electron repulsion integral (lA,lB|lC,lD).
    Angular-momentum dependent normalization constants are not applied yet at this level.
    """

    NeedsBraKetTransf = (lB > 0 or lD > 0) or (NeedsSpherTransf(lA, lB, lC, lD) and SpherTransf)
    NabCart = NFuncCart(lA) * NFuncCart(lB)
    NcdCart = NFuncCart(lC) * NFuncCart(lD)
    NabSpher = NFuncSpher(lA) * NFuncSpher(lB)
    NcdSpher = NFuncSpher(lC) * NFuncSpher(lD)
    if not NeedsBraKetTransf:
        Body = ("""G(1:{Nab}, 1:{Ncd}) = ZERO
call auto2e_W_{L}_{K}(G(1:{Nab}, 1:{Ncd}), Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, {Row0}, {Col0})"""
                 .format(L=lA, K=lC, Row0=1, Col0=1, Nab=NabCart, Ncd=NcdCart))
    else:
        DimA, DimB = WMatrixDims(lA, lB, lC, lD)
        Body = """real(F64), dimension({DimA}, {DimB}) :: W
!
! Compute (LS|KS) integrals, L={l0}..{l1} and K={k0}..{k1}
!
W = ZERO
""".format(DimA=DimA, DimB=DimB, l0=lA, l1=lA+lB, k0=lC, k1=lC+lD,)
        Col0 = 1
        for K in range(lC, lC+lD+1):
            Row0 = 1
            for L in range(lA, lA+lB+1):
                Body += """call auto2e_W_{L}_{K}(W, Ra, CntrA, ExpA, NprimA, Rb, CntrB, ExpB, NprimB, &
   Rc, CntrC, ExpC, NprimC, Rd, CntrD, ExpD, NprimD, Kappa, {Row0}, {Col0})\n""".format(L=L, K=K, Row0=Row0, Col0=Col0)
                Row0 += NFuncCart(L)
            Col0 += NFuncCart(K)
    Code = """
subroutine auto2e_eri{Spher}_{la}_{lb}_{lc}_{ld}(G, Ra, CntrA, ExpA, NprimA, &
Rb, CntrB, ExpB, NprimB, Rc, CntrC, ExpC, NprimC, &
Rd, CntrD, ExpD, NprimD, Kappa)
!
! Two-electron integrals class: ({a}{b}|{c}{d})
! L-dependent normalization constants are not applied.
! Output memory layout: {Layout}
!
real(F64), dimension({Nab}, *), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, ExpA, CntrB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, ExpC, CntrD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
{Body}{BraKetTransfer}
end subroutine auto2e_eri{Spher}_{la}_{lb}_{lc}_{ld}
""".format(la=lA, lb=lB, lc=lC, ld=lD, Body=Body,
           BraKetTransfer=BraKetTransfer_Call(lA, lB, lC, lD, SpherTransf) if NeedsBraKetTransf else "",
           a=LSymbols[lA], b=LSymbols[lB], c=LSymbols[lC], d=LSymbols[lD],
           Layout="".join(MemoryLayout).upper(),
           Nab=NabSpher if SpherTransf else NabCart,
           Spher="_Spher" if (NeedsSpherTransf(lA, lB, lC, lD) and SpherTransf) else "")
    return Code

def EquivalentPermutations(Order, lA, lB, lC, lD):
    """
    List all permutations of orbital indices that yield the same memory
    layout of integrals as the reference  permutation (Order).
    """
    L = {"a":lA, "b":lB, "c":lC, "d":lD}
    equiv = [Order]

    j = len(equiv)
    for i in range(j):
        p, q, r, s = equiv[i]
        if (L[p], L[q]) == (0, 0) or (L[r], L[s]) == (0, 0):
            equiv.append([r, s, p, q])

    j = len(equiv)
    for i in range(j):
        p, q, r, s = equiv[i]
        if L[p] == 0 or L[q] == 0:
            equiv.append([q, p, r, s])

    j = len(equiv)
    for i in range(j):
        p, q, r, s = equiv[i]
        if L[r] == 0 or L[s] == 0:
            equiv.append([p, q, s, r])
            
    return equiv

def UniquePermutations(lA, lB, lC, lD):
    All = []
    
    All.append(["a", "b", "c", "d"])
    All.append(["b", "a", "c", "d"])
    All.append(["a", "b", "d", "c"])
    All.append(["b", "a", "d", "c"])
    
    All.append(["c", "d", "a", "b"])
    All.append(["d", "c", "a", "b"])
    All.append(["c", "d", "b", "a"])
    All.append(["d", "c", "b", "a"])

    unique = []
    s = []
    for x in All:
        if x not in s:
            unique.append(sorted(EquivalentPermutations(x, lA, lB, lC, lD))[0])
            s += EquivalentPermutations(x, lA, lB, lC, lD)

    return unique

def MomentumQuadruples(A, B, C, D):
    M = {(A, B, C, D),
         (B, A, C, D),
         (A, B, D, C),
         (B, A, D, C),
         (C, D, A, B),
         (D, C, A, B),
         (C, D, B, A),
         (D, C, B, A)}
    return M

def ERI_Module(lA, lB, lC, lD, TargetLayout, WMatrix_Modules, SpherTransf):
    NeedsBraKetTransf = (lB > 0 or lD > 0) or (NeedsSpherTransf(lA, lB, lC, lD) and SpherTransf)
    if NeedsBraKetTransf:
        #
        # Generate the bra-ket momentum transfer subroutine for both Cartesian GTOs
        # and real spherical harmonics. The bra-ket transfer subroutine won't apply
        # the horizontal momentum transfer for LSLS integrals, but it will carry out
        # the Cartesian->spherical transformation.
        #
        if lB > 0 or lD > 0:
            BraKetTransfCart, OrderCart = BraKetTransfer_Subroutine(lA, lB, lC, lD, False)
        else:
            BraKetTransfCart = ""
            OrderCart = ["a", "b", "c", "d"]
        if NeedsSpherTransf(lA, lB, lC, lD) and SpherTransf:
            BraKetTransfSpher, OrderSpher = BraKetTransfer_Subroutine(lA, lB, lC, lD, True)
        else:
            BraKetTransfSpher = ""
            OrderSpher = OrderCart
        BraKetTransfer = BraKetTransfCart + BraKetTransfSpher
    else:
        BraKetTransfer = ""
        OrderCart = ["a", "b", "c", "d"]
        OrderSpher = ["a", "b", "c", "d"]
    ImportWMatrix = ""
    for ModuleName in WMatrix_Modules:
        ImportWMatrix += "use {}\n".format(ModuleName)
    InitialLayoutCart = sorted(EquivalentPermutations(OrderCart, lA, lB, lC, lD))[0]
    InitialLayoutSpher = sorted(EquivalentPermutations(OrderSpher, lA, lB, lC, lD))[0]

    Frontend = ""
    RequestedPermsCart = []
    RequestedPermsSpher = []
    for Momenta in MomentumQuadruples(lA, lB, lC, lD):
        p, s = Frontend_Subroutine(Momenta, TargetLayout, InitialLayoutCart, False)
        Frontend += s
        if p not in RequestedPermsCart:
            RequestedPermsCart.append(p)
        if NeedsSpherTransf(lA, lB, lC, lD) and SpherTransf:
            p, s = Frontend_Subroutine(Momenta, TargetLayout, InitialLayoutSpher, True)
            Frontend += s
            if p not in RequestedPermsSpher:
                RequestedPermsSpher.append(p)
            
    Normalize = ""
    for P in sorted(RequestedPermsCart):
        if P in EquivalentPermutations(InitialLayoutCart, lA, lB, lC, lD):
            Normalize += Normalize_Subroutine(lA, lB, lC, lD, InitialLayoutCart, False)
        else:
            Normalize += NormalizeTranspose_Subroutine(lA, lB, lC, lD, P, InitialLayoutCart, False)

    for P in sorted(RequestedPermsSpher):
        if P in EquivalentPermutations(InitialLayoutSpher, lA, lB, lC, lD):
            Normalize += Normalize_Subroutine(lA, lB, lC, lD, InitialLayoutSpher, True)
        else:
            Normalize += NormalizeTranspose_Subroutine(lA, lB, lC, lD, P, InitialLayoutSpher, True)
        
    DriverSubroutine = ERI_Subroutine(lA, lB, lC, lD, InitialLayoutCart, False)
    if NeedsSpherTransf(lA, lB, lC, lD) and SpherTransf:
            DriverSubroutine += ERI_Subroutine(lA, lB, lC, lD, InitialLayoutSpher, True)
    Module = """!
! Compute a batch of two-electron repulsion integrals, (AB|CD), in the Cartesian Gaussian
! and real spherical harmonic basis.
! Angular momenta: A={la}, B={lb}, C={lc}, and D={ld}
!
! The order of the angular-dependent prefactors (x-Ra(1))**lx*(y-Ra(2))**ly*(z-Ra(3))**lz
! is defined by the following loop:
! I = 1
! for lx in A..0
!     for ly in A-lx..0
!         lz = A-lx-ly
!         Index[(lx,ly,lz)] = I
!         I = I + 1
!
! The contracted Gaussian functions are defined as
!
! PhiA(lx,ly,lz) = NormA(Index[(lx,ly,lz)]) * (x-Ra(1))**lx * (y-Ra(2))**ly * (z-Ra(3))**lz * 
!     Sum(Ga=1..NprimA) CntrA(Ga) Exp(-ExpA(Ga) * ((x-Ra(1))**2 + (y-Ra(2))**2 + (z-Ra(3))**2))
!
! Kappa specifies the form of the interaction potential:
! Kappa=0 for the Coulomb interaction operator 1/r12
! Kappa=1/Omega**2 for the screened interaction potential Erf(Omega*r12)/r12
!
! Code generated automatically.
!
module auto2e_eri_{a}{b}{c}{d}
use arithmetic
use auto2e_KetTransfer
{SpherTransfModule}
{ImportWMatrix}
implicit none

contains
{Driver}
{BraKetTransfer}
{Normalize}
{Frontend}
end module auto2e_eri_{a}{b}{c}{d}
""".format(Driver=DriverSubroutine,
           BraKetTransfer=BraKetTransfer.strip(),
           Normalize=Normalize,
           a=LSymbols[lA], b=LSymbols[lB],
           c=LSymbols[lC], d=LSymbols[lD], 
           la=lA, lb=lB, lc=lC, ld=lD,
           Frontend=Frontend.strip(),
           ImportWMatrix=ImportWMatrix.strip(),
           SpherTransfModule="use auto2e_SpherTransf" if SpherTransf else "").strip()
    return Module

def Canonical4Tuple(Momenta, TargetLayout):
    #
    # The computational subroutines assume canonical order of angular momenta.
    # Determine the correct order of arguments, including the parameters of
    # PhiA, PhiB, PhiC, and PhiD, and the permutation that changes the memory
    # layout to the requested one.
    #
    # The canonical form of the angular momenta 4-tuple is A >= B, C >= D, AB >= CD
    #
    l1, l2, l3, l4 = Momenta
    m1, m2, m3, m4 = "a", "b", "c", "d"
    if l2 > l1:
        l1, l2 = l2, l1
        m1, m2 = m2, m1
    if l4 > l3:
        l3, l4 = l4, l3
        m3, m4 = m4, m3
    if l3 > l1 or (l3 == l1 and l4 > l2):
        l1, l2, l3, l4 = l3, l4, l1, l2
        m1, m2, m3, m4 = m3, m4, m1, m2
    #
    # Order in which the arguments (contraction coefficients etc.) should be passed
    # to subroutines which assume canonical ordering
    #
    CanonicalParams = [m1, m2, m3, m4]
    CanonicalMomenta = [l1, l2, l3, l4]
    d = {m1:"a", m2:"b", m3:"c", m4:"d"}
    #
    # CanonicalPermutation changes the memory layout to TargetLayout
    #
    Permutation = []
    for x in TargetLayout:
        Permutation.append(d[x])
    lA, lB, lC, lD = CanonicalMomenta
    CanonicalPermutation = sorted(EquivalentPermutations(Permutation, lA, lB, lC, lD))[0]
    return CanonicalParams, CanonicalMomenta, CanonicalPermutation

def Frontend_Subroutine(Momenta, TargetLayout, BaseLayout, SpherTransf):
    """
    Generate a front-end subroutine for a given quartet of angular momenta:
    1) compute the two-electron integrals using the memory alignment optimal for
    reducing the integral computation bottlenecks;
    2) normalize the resulting integrals using the angular-function dependent
    normalization factors and transpose the output to the requested memory alginment.
    """
    CanonicalParams, CanonicalMomenta, Permutation = Canonical4Tuple(Momenta, TargetLayout)
    a, b, c, d = CanonicalParams
    lA, lB, lC, lD = CanonicalMomenta
    IsTransposed = (BaseLayout != Permutation)
    if IsTransposed:
        if SpherTransf:
            Nabcd = NFuncSpher(lA) * NFuncSpher(lB) * NFuncSpher(lC) * NFuncSpher(lD)
        else:
            Nabcd = NFuncCart(lA) * NFuncCart(lB) * NFuncCart(lC) * NFuncCart(lD)
        Body = ("""real(F64), dimension({HDim}) :: H
call auto2e_eri{Spher}_{la}_{lb}_{lc}_{ld}(H, R{a}, Cntr{A}, Exp{A}, Nprim{A}, &
   R{b}, Cntr{B}, Exp{B}, Nprim{B}, &
   R{c}, Cntr{C}, Exp{C}, Nprim{C}, &
   R{d}, Cntr{D}, Exp{D}, Nprim{D}, Kappa)
call auto2e_Normalize{Spher}_{la}_{lb}_{lc}_{ld}_{LABEL}(G, H, Norm{A}, Norm{B}, Norm{C}, Norm{D})"""
                .format(HDim=Nabcd, la=lA, lb=lB, lc=lC, ld=lD, LABEL="".join(Permutation).upper(),
                        a=a, b=b, c=c, d=d, A=a.upper(), B=b.upper(), C=c.upper(), D=d.upper(),
                        Spher="_Spher" if (NeedsSpherTransf(lA, lB, lC, lD) and SpherTransf) else ""))
    else:
        Body = ("""call auto2e_eri{Spher}_{la}_{lb}_{lc}_{ld}(G, R{a}, Cntr{A}, Exp{A}, Nprim{A}, &
   R{b}, Cntr{B}, Exp{B}, Nprim{B}, &
   R{c}, Cntr{C}, Exp{C}, Nprim{C}, &
   R{d}, Cntr{D}, Exp{D}, Nprim{D}, Kappa)
call auto2e_Normalize{Spher}_{la}_{lb}_{lc}_{ld}_{LABEL}(G, Norm{A}, Norm{B}, Norm{C}, Norm{D})"""
                .format(la=lA, lb=lB, lc=lC, ld=lD, LABEL="".join(Permutation).upper(),
                        a=a, b=b, c=c, d=d, A=a.upper(), B=b.upper(), C=c.upper(), D=d.upper(),
                        Spher="_Spher" if (NeedsSpherTransf(lA, lB, lC, lD) and SpherTransf) else ""))
        
    Code = """
subroutine auto2e_frontend{Spher}_{p}_{q}_{r}_{s}(G, Ra, CntrA, NormA, ExpA, NprimA, &
Rb, CntrB, NormB, ExpB, NprimB, Rc, CntrC, NormC, ExpC, NprimC, &
Rd, CntrD, NormD, ExpD, NprimD, Kappa)
real(F64), dimension(*), intent(out) :: G
real(F64), dimension(3), intent(in) :: Ra, Rb, Rc, Rd
real(F64), dimension(*), intent(in) :: CntrA, NormA, ExpA, CntrB, NormB, ExpB
real(F64), dimension(*), intent(in) :: CntrC, NormC, ExpC, CntrD, NormD, ExpD
integer, intent(in) :: NprimA, NprimB, NprimC, NprimD
real(F64), intent(in) :: Kappa
{Body}
end subroutine auto2e_frontend{Spher}_{p}_{q}_{r}_{s}
""".format(p=Momenta[0], q=Momenta[1], r=Momenta[2], s=Momenta[3], Body=Body,
           Spher="_Spher" if (NeedsSpherTransf(lA, lB, lC, lD) and SpherTransf) else "")
    
    return Permutation, Code

def Init_Subroutine(MaxL, SpherTransf):
    PtrAssignments = []
    #
    # Cartesian integrals
    #
    for lD in range(MaxL+1):
        for lC in range(MaxL+1):
            for lB in range(MaxL+1):
                for lA in range(MaxL+1):
                    Index = 1 + lA + lB * (MaxL+1) + lC * (MaxL+1)**2 + lD * (MaxL+1)**3
                    PtrAssignments.append("Auto2eERI({})%ptr => auto2e_frontend_{}_{}_{}_{}".format(Index, lA, lB, lC, lD))
    #
    # Spherical integrals
    #
    offset = (MaxL+1)**4
    if SpherTransf:
        for lD in range(MaxL+1):
            for lC in range(MaxL+1):
                for lB in range(MaxL+1):
                    for lA in range(MaxL+1):
                        Index = 1 + lA + lB * (MaxL+1) + lC * (MaxL+1)**2 + lD * (MaxL+1)**3
                        if NeedsSpherTransf(lA, lB, lC, lD):
                            SpherPostfix = "_Spher"
                        else:
                            SpherPostfix = ""
                        PtrAssignments.append(
                            "Auto2eERI({})%ptr => auto2e_frontend{}_{}_{}_{}_{}"
                            .format(offset+Index, SpherPostfix, lA, lB, lC, lD))
    Subroutine = """
subroutine auto2e_init()
call auto2e_SpherTransf_init()
{}
end subroutine auto2e_init
""".format("\n".join(PtrAssignments))
    return Subroutine

def Driver_Module(MaxL, MemoryLayout, SpherTransf):
    Modules = "use auto2e_SpherTransf\n"
    for lD in range(MaxL+1):
        for lC in range(lD, MaxL+1):
            for lB in range(MaxL+1):
                for lA in range(lB, MaxL+1):
                    if lA > lC or (lA == lC and lB >= lD):
                        ModuleName = "auto2e_eri_{}{}{}{}".format(LSymbols[lA], LSymbols[lB], LSymbols[lC], LSymbols[lD])
                        Modules += "use {}\n".format(ModuleName)
    Code = ("""!
! -------------------------------------------------------------------------
! Automatically-generated code for two-electron Coulomb repulsion integrals
! -------------------------------------------------------------------------
! The main steps which dominate the computation time are (1) the Hermite ->
! Coulomb transformation of (LS|KS)-type integrals, and (2) the horizontal
! momentum transfer (A+BS|C+DS) -> (AB|CD).
!
! This module contains the array of pointers to the subroutines generated
! for each angular momenta 4-tuple. Example usage:
!
! Cartesian basis
! ---------------
! call AUTO2EERI(auto2e_idx(SHTYPE(ShA), SHTYPE(ShB), SHTYPE(ShC), SHTYPE(ShD)))%ptr(Gabcd, &
!      ATOMR(:, A), CNTR(:, ShA), CNTRNORM(:, ShA), EXPN(:, ShA), NPRM(ShA), &
!      ATOMR(:, B), CNTR(:, ShB), CNTRNORM(:, ShB), EXPN(:, ShB), NPRM(ShB), &
!      ATOMR(:, C), CNTR(:, ShC), CNTRNORM(:, ShC), EXPN(:, ShC), NPRM(ShC), &
!      ATOMR(:, D), CNTR(:, ShD), CNTRNORM(:, ShD), EXPN(:, ShD), NPRM(ShD), &
!      Kappa)
!
! Spherical basis
! ---------------
! call AUTO2EERI(AUTO2E_SPHER_OFFSET+auto2e_idx(SHTYPE(ShA), SHTYPE(ShB), &
!      SHTYPE(ShC), SHTYPE(ShD)))%ptr(Gabcd, &
!      ATOMR(:, A), CNTR(:, ShA), CNTRNORM(:, ShA), EXPN(:, ShA), NPRM(ShA), &
!      ATOMR(:, B), CNTR(:, ShB), CNTRNORM(:, ShB), EXPN(:, ShB), NPRM(ShB), &
!      ATOMR(:, C), CNTR(:, ShC), CNTRNORM(:, ShC), EXPN(:, ShC), NPRM(ShC), &
!      ATOMR(:, D), CNTR(:, ShD), CNTRNORM(:, ShD), EXPN(:, ShD), NPRM(ShD), &
!      Kappa)
!
! Maximum angular momentum: {MaxL}
! Memory layout of the output array: {Layout}
! Code generated automatically on {datestr}
!
module auto2e
{ModuleImport}
implicit none

type TAuto2ePtr
procedure(auto2e_frontend_0_0_0_0), pointer, nopass :: ptr
end type TAuto2ePtr
integer, parameter :: AUTO2E_SPHER_OFFSET = {M}
integer, parameter :: AUTO2E_MAXL = {MaxL}
type(TAuto2ePtr), dimension({N}), save :: Auto2eERI

contains
{Init}
function auto2e_idx(La, Lb, Lc, Ld)
integer :: auto2e_idx
integer, intent(in) :: La, Lb, Lc, Ld
auto2e_idx = 1 + La + Lb * {Wb} + Lc * {Wc} + Ld * {Wd}
end function auto2e_idx
end module auto2e"""
.format(ModuleImport=Modules, M=(MaxL+1)**4, N=2*(MaxL+1)**4 if SpherTransf else (MaxL+1)**4,
        Init=Init_Subroutine(MaxL, SpherTransf), Wb=MaxL+1, Wc=(MaxL+1)**2,
        Wd=(MaxL+1)**3, Layout="".join(MemoryLayout).upper(), MaxL=MaxL, datestr=DateStamp()))
    return Code

def MakeAll(MaxL, TargetLayout):
    #
    # Generate code for integrals in the real spherical harmonics basis
    #
    SpherTransf = True
    #
    # Set up directories
    #
    RootDir = path.dirname(path.realpath(__file__))
    SourceDir = path.relpath(path.join(RootDir, "src"), RootDir)
    if not path.exists(SourceDir):
        os.makedirs(SourceDir)
    FileList = []
    Batch = []
    #
    # Cartesian->spherical transformation
    #
    if SpherTransf:
        Spher = SpherTransf_Module(MaxL)
        FilePath = path.join(SourceDir, "auto2e_SpherTransf.f90")
        f = open(FilePath, "w")
        f.write(Spher)
        f.close()
        Batch.append(FilePath)
    #
    # Two-electron integrals in the Hermite Gaussian basis
    # Hermite -> Gaussian transformation coefficients
    #
    Hermite = HermiteModule(MaxL)
    FilePath = path.join(SourceDir, "auto2e_Hermite.f90")
    f = open(FilePath, "w")
    f.write(Hermite)
    f.close()
    Batch.append(FilePath)
    #
    # Hermite -> Gaussian transformation of bra shell pairs
    #
    BraTransform = HermiteToCartesian_BraTransformModule(MaxL)
    FilePath = path.join(SourceDir, "auto2e_BraTransform.f90")
    f = open(FilePath, "w")
    f.write(BraTransform)
    f.close()
    Batch.append(FilePath)
    #
    # Hermite -> Gaussian transformation of ket shell pairs
    #
    KetTransformModules = {}
    for LKet in range(1, 2*MaxL+1):
        LBraMin = max(2, LKet)
        LBraMax = 2*MaxL
        for LBra in range(LBraMin, LBraMax+1):
            if EnableRCopy(LBra, LKet):
                KetTransform = RCopy_KetTransform_Module(LBra, LKet)
            else:
                KetTransform = KetTransform_Module(LBra, LKet)
            KetTransformModules[(LBra, LKet)] = []
            for name, SourceCode in KetTransform:
                FilePath = path.join(SourceDir, "{}.f90".format(name))
                f = open(FilePath, "w")
                f.write(SourceCode)
                f.close()
                Batch.append(FilePath)
                KetTransformModules[(LBra, LKet)].append(name)
    #
    # Transfer ket shell pair angular momentum using the horizontal
    # recurrence relation (Head-Gordon and Pople)
    #
    KetTransfer = KetTransfer_Module(MaxL)
    FilePath = path.join(SourceDir, "auto2e_KetTransfer.f90")
    f = open(FilePath, "w")
    f.write(KetTransfer)
    f.close()
    Batch.append(FilePath)
    FileList.append(sorted(Batch))
    #
    # (LS|KS)
    #
    WMatrix_ModuleNames = []
    Batch = []
    WMatrix = WMatrix_Module(MaxL, reduce(lambda x, y: x+y, KetTransformModules.values(), []))
    for ModuleName, Code in WMatrix:
        FilePath = path.join(SourceDir, "{}.f90".format(ModuleName))
        f = open(FilePath, "w")
        f.write(Code)
        f.close()
        Batch.append(FilePath)
        WMatrix_ModuleNames.append(ModuleName)
    FileList.append(Batch)
    ERI_batch = []
    for lD in range(MaxL+1):
        for lC in range(lD, MaxL+1):
            for lB in range(MaxL+1):
                for lA in range(lB, MaxL+1):
                    if lA > lC or (lA == lC and lB >= lD):
                        ERIModule = ERI_Module(lA, lB, lC, lD, TargetLayout, WMatrix_ModuleNames, SpherTransf)
                        FilePath = path.join(SourceDir, "auto2e_eri_{a}{b}{c}{d}.f90".format(a=LSymbols[lA], b=LSymbols[lB], c=LSymbols[lC],
                                            d=LSymbols[lD]))
                        f = open(FilePath, "w")
                        f.write(ERIModule)
                        f.close()
                        ERI_batch.append(FilePath)
    FileList.append(sorted(ERI_batch))
    driver = Driver_Module(MaxL, TargetLayout, SpherTransf)
    FilePath = path.join(SourceDir, "auto2e.f90")
    f = open(FilePath, "w")
    f.write(driver)
    f.close()
    FileList.append([FilePath])
    ParallelBuild.BuildScript(FileList, path.join(RootDir, "Build.py"))

MakeAll(5, ["d", "c", "b", "a"])
