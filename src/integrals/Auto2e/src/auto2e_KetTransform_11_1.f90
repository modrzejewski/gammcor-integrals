module auto2e_KetTransform_11_1
use arithmetic
use math_constants    
implicit none

integer, dimension(78, 3), parameter :: RCopyIdx_11_1 = reshape([1, 14, 26, 37, 47, &
    56, 64, 71, 77, 82, 86, 89, 92, 104, 115, 125, 134, 142, 149, 155, 160, 164, &
    167, 170, 181, 191, 200, 208, 215, 221, 226, 230, 233, 236, 246, 255, 263, 270, &
    276, 281, 285, 288, 291, 300, 308, 315, 321, 326, 330, 333, 336, 344, 351, 357, &
    362, 366, 369, 372, 379, 385, 390, 394, 397, 400, 406, 411, 415, 418, 421, 426, &
    430, 433, 436, 440, 443, 446, 449, 452, 14, 26, 37, 47, 56, 64, 71, 77, 82, 86, &
    89, 91, 104, 115, 125, 134, 142, 149, 155, 160, 164, 167, 169, 181, 191, 200, &
    208, 215, 221, 226, 230, 233, 235, 246, 255, 263, 270, 276, 281, 285, 288, 290, &
    300, 308, 315, 321, 326, 330, 333, 335, 344, 351, 357, 362, 366, 369, 371, 379, &
    385, 390, 394, 397, 399, 406, 411, 415, 418, 420, 426, 430, 433, 435, 440, 443, &
    445, 449, 451, 454, 92, 104, 115, 125, 134, 142, 149, 155, 160, 164, 167, 169, &
    170, 181, 191, 200, 208, 215, 221, 226, 230, 233, 235, 236, 246, 255, 263, 270, &
    276, 281, 285, 288, 290, 291, 300, 308, 315, 321, 326, 330, 333, 335, 336, 344, &
    351, 357, 362, 366, 369, 371, 372, 379, 385, 390, 394, 397, 399, 400, 406, 411, &
    415, 418, 420, 421, 426, 430, 433, 435, 436, 440, 443, 445, 446, 449, 451, 452, &
    454, 455], [78, 3])

contains

subroutine auto2e_RCopy_11_1(T, R, idx)
real(F64), dimension(:), intent(out) :: T
real(F64), dimension(:), intent(in) :: R
integer, dimension(:), intent(in) :: idx
T(1:12) = R(idx(1):idx(1)+11)
T(13:23) = R(idx(2):idx(2)+10)
T(24:33) = R(idx(3):idx(3)+9)
T(34:42) = R(idx(4):idx(4)+8)
T(43:50) = R(idx(5):idx(5)+7)
T(51:57) = R(idx(6):idx(6)+6)
T(58:63) = R(idx(7):idx(7)+5)
T(64:68) = R(idx(8):idx(8)+4)
T(69:72) = R(idx(9):idx(9)+3)
T(73:75) = R(idx(10):idx(10)+2)
T(76:77) = R(idx(11):idx(11)+1)
T(78) = R(idx(12))
T(79:89) = R(idx(13):idx(13)+10)
T(90:99) = R(idx(14):idx(14)+9)
T(100:108) = R(idx(15):idx(15)+8)
T(109:116) = R(idx(16):idx(16)+7)
T(117:123) = R(idx(17):idx(17)+6)
T(124:129) = R(idx(18):idx(18)+5)
T(130:134) = R(idx(19):idx(19)+4)
T(135:138) = R(idx(20):idx(20)+3)
T(139:141) = R(idx(21):idx(21)+2)
T(142:143) = R(idx(22):idx(22)+1)
T(144) = R(idx(23))
T(145:154) = R(idx(24):idx(24)+9)
T(155:163) = R(idx(25):idx(25)+8)
T(164:171) = R(idx(26):idx(26)+7)
T(172:178) = R(idx(27):idx(27)+6)
T(179:184) = R(idx(28):idx(28)+5)
T(185:189) = R(idx(29):idx(29)+4)
T(190:193) = R(idx(30):idx(30)+3)
T(194:196) = R(idx(31):idx(31)+2)
T(197:198) = R(idx(32):idx(32)+1)
T(199) = R(idx(33))
T(200:208) = R(idx(34):idx(34)+8)
T(209:216) = R(idx(35):idx(35)+7)
T(217:223) = R(idx(36):idx(36)+6)
T(224:229) = R(idx(37):idx(37)+5)
T(230:234) = R(idx(38):idx(38)+4)
T(235:238) = R(idx(39):idx(39)+3)
T(239:241) = R(idx(40):idx(40)+2)
T(242:243) = R(idx(41):idx(41)+1)
T(244) = R(idx(42))
T(245:252) = R(idx(43):idx(43)+7)
T(253:259) = R(idx(44):idx(44)+6)
T(260:265) = R(idx(45):idx(45)+5)
T(266:270) = R(idx(46):idx(46)+4)
T(271:274) = R(idx(47):idx(47)+3)
T(275:277) = R(idx(48):idx(48)+2)
T(278:279) = R(idx(49):idx(49)+1)
T(280) = R(idx(50))
T(281:287) = R(idx(51):idx(51)+6)
T(288:293) = R(idx(52):idx(52)+5)
T(294:298) = R(idx(53):idx(53)+4)
T(299:302) = R(idx(54):idx(54)+3)
T(303:305) = R(idx(55):idx(55)+2)
T(306:307) = R(idx(56):idx(56)+1)
T(308) = R(idx(57))
T(309:314) = R(idx(58):idx(58)+5)
T(315:319) = R(idx(59):idx(59)+4)
T(320:323) = R(idx(60):idx(60)+3)
T(324:326) = R(idx(61):idx(61)+2)
T(327:328) = R(idx(62):idx(62)+1)
T(329) = R(idx(63))
T(330:334) = R(idx(64):idx(64)+4)
T(335:338) = R(idx(65):idx(65)+3)
T(339:341) = R(idx(66):idx(66)+2)
T(342:343) = R(idx(67):idx(67)+1)
T(344) = R(idx(68))
T(345:348) = R(idx(69):idx(69)+3)
T(349:351) = R(idx(70):idx(70)+2)
T(352:353) = R(idx(71):idx(71)+1)
T(354) = R(idx(72))
T(355:357) = R(idx(73):idx(73)+2)
T(358:359) = R(idx(74):idx(74)+1)
T(360) = R(idx(75))
T(361:362) = R(idx(76):idx(76)+1)
T(363) = R(idx(77))
T(364) = R(idx(78))
end subroutine auto2e_RCopy_11_1

subroutine auto2e_RCopy_KetTransform_11_1(S, T, Ex, Ey, Ez, R, lx, ly, lz)
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
i = ((2*1+3-phi)*phi)/2+nu+1
call auto2e_Rcopy_11_1(T, R(tau+1:), RCopyIdx_11_1(:, i))
x = x0 + tau
y = y0 + nu
z = z0 + phi
c = (-1)**modulo(tau+nu+phi, 2)*Ex(x)*Ey(y)*Ez(z)
S = S + c * T
end do
end do
end do
end subroutine auto2e_RCopy_KetTransform_11_1
end module auto2e_KetTransform_11_1