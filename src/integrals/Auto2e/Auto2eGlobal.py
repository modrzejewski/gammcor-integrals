#!/usr/bin/env python3
# ---------------------------------------------------------------------
#        Global parameters controlling automatic code generation
#        for molecular integrals. In addition, this module contains
#        small subroutines shared by all other code generation
#        modules.
# ---------------------------------------------------------------------
import datetime
#
# Real kind used in Fortran code
#
REAL_KIND_LABEL = "F64"
REAL_KIND_PREC = 17
#
# Maximum angular momentum for which indices are precomputed.
# The angular momenta requested by the user are less or equal MaxL.
#
MaxL = 6
LSymbols = ["s", "p", "d", "f", "g", "h", "i", "k"]
#
# Minimum value of the angular momentum for which the Cartesian->spherical
# transformation is applied. No rescaling or reordering is applied
# for L<SPHER_TRANSF_LMIN, thus, these functions are treated as
# Cartesian Gaussians.
#
SPHER_TRANSF_LMIN = 2
#
# Ordering of the Gaussian orbital angular (xyz) prefactors.
# Note that the ordering of xyz functions has to be
# the same as in the loops over lx,ly,lz in BraKetTransform
# subroutines in the Auto2e module.
#
AngularFunctions = {}
XYZIdx = {}
for L in range(2*MaxL+1):
    AngularFunctions[L] = []
    i = 1
    for x in range(L, -1, -1):
        for y in range(L-x, -1, -1):
            z = L - x - y
            XYZIdx[(x, y, z)] = i
            AngularFunctions[L].append((x, y, z))
            i += 1
#
# Indices of the Hermite->Cartesian transformation coefficients Et^{L,0}
#
EtIdx = {}
i = 1
for L in range(2*MaxL+1):
    for t in range(L+1):
        EtIdx[(L, t)] = i
        i += 1

        
def NFuncCart(l):
    """ Number of angular functions in a Cartesian orbital shell."""
    return ((l + 1) * (l + 2)) // 2


def NFuncSpher(l):
    """ Number of angular spherical functionas in an orbital shell. """
    return 2 * l + 1

def DateStamp():
    utc_datetime = datetime.datetime.utcnow()
    return utc_datetime.strftime("%Y-%m-%d %H:%M:%S")


def SplitLine(line, delims, indent=3, TargetSize=80, CommentLine=False):
    FormattedLine = ""
    CurrentBlock = ""
    for x in line:
        if x in delims:
            if len(CurrentBlock) >= TargetSize:
                FormattedLine += CurrentBlock + " &\n"
                CurrentBlock = " "*indent + x
                if CommentLine:
                    CurrentBlock = "!" + CurrentBlock
            else:
                CurrentBlock += x
        else:
            CurrentBlock += x
    FormattedLine += CurrentBlock
    return FormattedLine


def FormatSum(y, c, x):
    #
    # Generate a formatted string y = c[0] * x[0] + c[1] * x[1] + ...
    #
    s = "{} =".format(y)
    for k in range(len(c)):
        if k == 0:
            s += " {c:.{prec}g}_{fmt} * {x}".format(c=c[k], x=x[k], prec=REAL_KIND_PREC, fmt=REAL_KIND_LABEL)
        else:
            if c[k] > 0:
                sign = "+"
            else:
                sign = "-"
            s += " {sign} {c:.{prec}g}_{fmt} * {x}".format(c=abs(c[k]), x=x[k], prec=REAL_KIND_PREC, fmt=REAL_KIND_LABEL, sign=sign)
    return s


def GabIndex(A, B):
    """ Compound index of a matrix Gab, where a and b are Cartesian angular function indices. """
    lA = sum(A)
    nA = NFuncCart(lA)
    iA = XYZIdx[tuple(A)]
    iB = XYZIdx[tuple(B)]
    return iA + (iB - 1) * nA


def GabIndexSpher(mA, lA, mB, lB):
    """ Compound index of a matrix Gab, where a and b are spherical angular function indices. """
    iA = mA + lA + 1
    iB = mB + lB + 1
    nA = 2 * lA + 1
    return iA + (iB - 1) * nA
    
