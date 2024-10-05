module auto2e_KetTransform_10_2
use arithmetic
use math_constants
implicit none

integer, dimension(66, 6), parameter :: RCopyIdx_10_2 = reshape([1, 14, 26, 37, 47, &
    56, 64, 71, 77, 82, 86, 92, 104, 115, 125, 134, 142, 149, 155, 160, 164, 170, &
    181, 191, 200, 208, 215, 221, 226, 230, 236, 246, 255, 263, 270, 276, 281, 285, &
    291, 300, 308, 315, 321, 326, 330, 336, 344, 351, 357, 362, 366, 372, 379, 385, &
    390, 394, 400, 406, 411, 415, 421, 426, 430, 436, 440, 446, 14, 26, 37, 47, 56, &
    64, 71, 77, 82, 86, 89, 104, 115, 125, 134, 142, 149, 155, 160, 164, 167, 181, &
    191, 200, 208, 215, 221, 226, 230, 233, 246, 255, 263, 270, 276, 281, 285, 288, &
    300, 308, 315, 321, 326, 330, 333, 344, 351, 357, 362, 366, 369, 379, 385, 390, &
    394, 397, 406, 411, 415, 418, 426, 430, 433, 440, 443, 449, 26, 37, 47, 56, 64, &
    71, 77, 82, 86, 89, 91, 115, 125, 134, 142, 149, 155, 160, 164, 167, 169, 191, &
    200, 208, 215, 221, 226, 230, 233, 235, 255, 263, 270, 276, 281, 285, 288, 290, &
    308, 315, 321, 326, 330, 333, 335, 351, 357, 362, 366, 369, 371, 385, 390, 394, &
    397, 399, 411, 415, 418, 420, 430, 433, 435, 443, 445, 451, 92, 104, 115, 125, &
    134, 142, 149, 155, 160, 164, 167, 170, 181, 191, 200, 208, 215, 221, 226, 230, &
    233, 236, 246, 255, 263, 270, 276, 281, 285, 288, 291, 300, 308, 315, 321, 326, &
    330, 333, 336, 344, 351, 357, 362, 366, 369, 372, 379, 385, 390, 394, 397, 400, &
    406, 411, 415, 418, 421, 426, 430, 433, 436, 440, 443, 446, 449, 452, 104, 115, &
    125, 134, 142, 149, 155, 160, 164, 167, 169, 181, 191, 200, 208, 215, 221, 226, &
    230, 233, 235, 246, 255, 263, 270, 276, 281, 285, 288, 290, 300, 308, 315, 321, &
    326, 330, 333, 335, 344, 351, 357, 362, 366, 369, 371, 379, 385, 390, 394, 397, &
    399, 406, 411, 415, 418, 420, 426, 430, 433, 435, 440, 443, 445, 449, 451, 454, &
    170, 181, 191, 200, 208, 215, 221, 226, 230, 233, 235, 236, 246, 255, 263, 270, &
    276, 281, 285, 288, 290, 291, 300, 308, 315, 321, 326, 330, 333, 335, 336, 344, &
    351, 357, 362, 366, 369, 371, 372, 379, 385, 390, 394, 397, 399, 400, 406, 411, &
    415, 418, 420, 421, 426, 430, 433, 435, 436, 440, 443, 445, 446, 449, 451, 452, &
    454, 455], [66, 6])

contains

subroutine auto2e_RCopy_10_2(T, R, idx)
real(F64), dimension(:), intent(out) :: T
real(F64), dimension(:), intent(in) :: R
integer, dimension(:), intent(in) :: idx
T(1:11) = R(idx(1):idx(1)+10)
T(12:21) = R(idx(2):idx(2)+9)
T(22:30) = R(idx(3):idx(3)+8)
T(31:38) = R(idx(4):idx(4)+7)
T(39:45) = R(idx(5):idx(5)+6)
T(46:51) = R(idx(6):idx(6)+5)
T(52:56) = R(idx(7):idx(7)+4)
T(57:60) = R(idx(8):idx(8)+3)
T(61:63) = R(idx(9):idx(9)+2)
T(64:65) = R(idx(10):idx(10)+1)
T(66) = R(idx(11))
T(67:76) = R(idx(12):idx(12)+9)
T(77:85) = R(idx(13):idx(13)+8)
T(86:93) = R(idx(14):idx(14)+7)
T(94:100) = R(idx(15):idx(15)+6)
T(101:106) = R(idx(16):idx(16)+5)
T(107:111) = R(idx(17):idx(17)+4)
T(112:115) = R(idx(18):idx(18)+3)
T(116:118) = R(idx(19):idx(19)+2)
T(119:120) = R(idx(20):idx(20)+1)
T(121) = R(idx(21))
T(122:130) = R(idx(22):idx(22)+8)
T(131:138) = R(idx(23):idx(23)+7)
T(139:145) = R(idx(24):idx(24)+6)
T(146:151) = R(idx(25):idx(25)+5)
T(152:156) = R(idx(26):idx(26)+4)
T(157:160) = R(idx(27):idx(27)+3)
T(161:163) = R(idx(28):idx(28)+2)
T(164:165) = R(idx(29):idx(29)+1)
T(166) = R(idx(30))
T(167:174) = R(idx(31):idx(31)+7)
T(175:181) = R(idx(32):idx(32)+6)
T(182:187) = R(idx(33):idx(33)+5)
T(188:192) = R(idx(34):idx(34)+4)
T(193:196) = R(idx(35):idx(35)+3)
T(197:199) = R(idx(36):idx(36)+2)
T(200:201) = R(idx(37):idx(37)+1)
T(202) = R(idx(38))
T(203:209) = R(idx(39):idx(39)+6)
T(210:215) = R(idx(40):idx(40)+5)
T(216:220) = R(idx(41):idx(41)+4)
T(221:224) = R(idx(42):idx(42)+3)
T(225:227) = R(idx(43):idx(43)+2)
T(228:229) = R(idx(44):idx(44)+1)
T(230) = R(idx(45))
T(231:236) = R(idx(46):idx(46)+5)
T(237:241) = R(idx(47):idx(47)+4)
T(242:245) = R(idx(48):idx(48)+3)
T(246:248) = R(idx(49):idx(49)+2)
T(249:250) = R(idx(50):idx(50)+1)
T(251) = R(idx(51))
T(252:256) = R(idx(52):idx(52)+4)
T(257:260) = R(idx(53):idx(53)+3)
T(261:263) = R(idx(54):idx(54)+2)
T(264:265) = R(idx(55):idx(55)+1)
T(266) = R(idx(56))
T(267:270) = R(idx(57):idx(57)+3)
T(271:273) = R(idx(58):idx(58)+2)
T(274:275) = R(idx(59):idx(59)+1)
T(276) = R(idx(60))
T(277:279) = R(idx(61):idx(61)+2)
T(280:281) = R(idx(62):idx(62)+1)
T(282) = R(idx(63))
T(283:284) = R(idx(64):idx(64)+1)
T(285) = R(idx(65))
T(286) = R(idx(66))
end subroutine auto2e_RCopy_10_2

subroutine auto2e_RCopy_KetTransform_10_2(S, T, Ex, Ey, Ez, R, lx, ly, lz)
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
i = ((2*2+3-phi)*phi)/2+nu+1
call auto2e_Rcopy_10_2(T, R(tau+1:), RCopyIdx_10_2(:, i))
x = x0 + tau
y = y0 + nu
z = z0 + phi
c = (-1)**modulo(tau+nu+phi, 2)*Ex(x)*Ey(y)*Ez(z)
S = S + c * T
end do
end do
end do
end subroutine auto2e_RCopy_KetTransform_10_2
end module auto2e_KetTransform_10_2