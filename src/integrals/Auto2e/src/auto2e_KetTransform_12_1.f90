module auto2e_KetTransform_12_1
use arithmetic
use math_constants
implicit none

integer, dimension(91, 3), parameter :: RCopyIdx_12_1 = reshape([1, 15, 28, 40, 51, &
    61, 70, 78, 85, 91, 96, 100, 103, 106, 119, 131, 142, 152, 161, 169, 176, 182, &
    187, 191, 194, 197, 209, 220, 230, 239, 247, 254, 260, 265, 269, 272, 275, 286, &
    296, 305, 313, 320, 326, 331, 335, 338, 341, 351, 360, 368, 375, 381, 386, 390, &
    393, 396, 405, 413, 420, 426, 431, 435, 438, 441, 449, 456, 462, 467, 471, 474, &
    477, 484, 490, 495, 499, 502, 505, 511, 516, 520, 523, 526, 531, 535, 538, 541, &
    545, 548, 551, 554, 557, 15, 28, 40, 51, 61, 70, 78, 85, 91, 96, 100, 103, 105, &
    119, 131, 142, 152, 161, 169, 176, 182, 187, 191, 194, 196, 209, 220, 230, 239, &
    247, 254, 260, 265, 269, 272, 274, 286, 296, 305, 313, 320, 326, 331, 335, 338, &
    340, 351, 360, 368, 375, 381, 386, 390, 393, 395, 405, 413, 420, 426, 431, 435, &
    438, 440, 449, 456, 462, 467, 471, 474, 476, 484, 490, 495, 499, 502, 504, 511, &
    516, 520, 523, 525, 531, 535, 538, 540, 545, 548, 550, 554, 556, 559, 106, 119, &
    131, 142, 152, 161, 169, 176, 182, 187, 191, 194, 196, 197, 209, 220, 230, 239, &
    247, 254, 260, 265, 269, 272, 274, 275, 286, 296, 305, 313, 320, 326, 331, 335, &
    338, 340, 341, 351, 360, 368, 375, 381, 386, 390, 393, 395, 396, 405, 413, 420, &
    426, 431, 435, 438, 440, 441, 449, 456, 462, 467, 471, 474, 476, 477, 484, 490, &
    495, 499, 502, 504, 505, 511, 516, 520, 523, 525, 526, 531, 535, 538, 540, 541, &
    545, 548, 550, 551, 554, 556, 557, 559, 560], [91, 3])

contains

subroutine auto2e_RCopy_12_1(T, R, idx)
real(F64), dimension(:), intent(out) :: T
real(F64), dimension(:), intent(in) :: R
integer, dimension(:), intent(in) :: idx
T(1:13) = R(idx(1):idx(1)+12)
T(14:25) = R(idx(2):idx(2)+11)
T(26:36) = R(idx(3):idx(3)+10)
T(37:46) = R(idx(4):idx(4)+9)
T(47:55) = R(idx(5):idx(5)+8)
T(56:63) = R(idx(6):idx(6)+7)
T(64:70) = R(idx(7):idx(7)+6)
T(71:76) = R(idx(8):idx(8)+5)
T(77:81) = R(idx(9):idx(9)+4)
T(82:85) = R(idx(10):idx(10)+3)
T(86:88) = R(idx(11):idx(11)+2)
T(89:90) = R(idx(12):idx(12)+1)
T(91) = R(idx(13))
T(92:103) = R(idx(14):idx(14)+11)
T(104:114) = R(idx(15):idx(15)+10)
T(115:124) = R(idx(16):idx(16)+9)
T(125:133) = R(idx(17):idx(17)+8)
T(134:141) = R(idx(18):idx(18)+7)
T(142:148) = R(idx(19):idx(19)+6)
T(149:154) = R(idx(20):idx(20)+5)
T(155:159) = R(idx(21):idx(21)+4)
T(160:163) = R(idx(22):idx(22)+3)
T(164:166) = R(idx(23):idx(23)+2)
T(167:168) = R(idx(24):idx(24)+1)
T(169) = R(idx(25))
T(170:180) = R(idx(26):idx(26)+10)
T(181:190) = R(idx(27):idx(27)+9)
T(191:199) = R(idx(28):idx(28)+8)
T(200:207) = R(idx(29):idx(29)+7)
T(208:214) = R(idx(30):idx(30)+6)
T(215:220) = R(idx(31):idx(31)+5)
T(221:225) = R(idx(32):idx(32)+4)
T(226:229) = R(idx(33):idx(33)+3)
T(230:232) = R(idx(34):idx(34)+2)
T(233:234) = R(idx(35):idx(35)+1)
T(235) = R(idx(36))
T(236:245) = R(idx(37):idx(37)+9)
T(246:254) = R(idx(38):idx(38)+8)
T(255:262) = R(idx(39):idx(39)+7)
T(263:269) = R(idx(40):idx(40)+6)
T(270:275) = R(idx(41):idx(41)+5)
T(276:280) = R(idx(42):idx(42)+4)
T(281:284) = R(idx(43):idx(43)+3)
T(285:287) = R(idx(44):idx(44)+2)
T(288:289) = R(idx(45):idx(45)+1)
T(290) = R(idx(46))
T(291:299) = R(idx(47):idx(47)+8)
T(300:307) = R(idx(48):idx(48)+7)
T(308:314) = R(idx(49):idx(49)+6)
T(315:320) = R(idx(50):idx(50)+5)
T(321:325) = R(idx(51):idx(51)+4)
T(326:329) = R(idx(52):idx(52)+3)
T(330:332) = R(idx(53):idx(53)+2)
T(333:334) = R(idx(54):idx(54)+1)
T(335) = R(idx(55))
T(336:343) = R(idx(56):idx(56)+7)
T(344:350) = R(idx(57):idx(57)+6)
T(351:356) = R(idx(58):idx(58)+5)
T(357:361) = R(idx(59):idx(59)+4)
T(362:365) = R(idx(60):idx(60)+3)
T(366:368) = R(idx(61):idx(61)+2)
T(369:370) = R(idx(62):idx(62)+1)
T(371) = R(idx(63))
T(372:378) = R(idx(64):idx(64)+6)
T(379:384) = R(idx(65):idx(65)+5)
T(385:389) = R(idx(66):idx(66)+4)
T(390:393) = R(idx(67):idx(67)+3)
T(394:396) = R(idx(68):idx(68)+2)
T(397:398) = R(idx(69):idx(69)+1)
T(399) = R(idx(70))
T(400:405) = R(idx(71):idx(71)+5)
T(406:410) = R(idx(72):idx(72)+4)
T(411:414) = R(idx(73):idx(73)+3)
T(415:417) = R(idx(74):idx(74)+2)
T(418:419) = R(idx(75):idx(75)+1)
T(420) = R(idx(76))
T(421:425) = R(idx(77):idx(77)+4)
T(426:429) = R(idx(78):idx(78)+3)
T(430:432) = R(idx(79):idx(79)+2)
T(433:434) = R(idx(80):idx(80)+1)
T(435) = R(idx(81))
T(436:439) = R(idx(82):idx(82)+3)
T(440:442) = R(idx(83):idx(83)+2)
T(443:444) = R(idx(84):idx(84)+1)
T(445) = R(idx(85))
T(446:448) = R(idx(86):idx(86)+2)
T(449:450) = R(idx(87):idx(87)+1)
T(451) = R(idx(88))
T(452:453) = R(idx(89):idx(89)+1)
T(454) = R(idx(90))
T(455) = R(idx(91))
end subroutine auto2e_RCopy_12_1

subroutine auto2e_RCopy_KetTransform_12_1(S, T, Ex, Ey, Ez, R, lx, ly, lz)
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
call auto2e_Rcopy_12_1(T, R(tau+1:), RCopyIdx_12_1(:, i))
x = x0 + tau
y = y0 + nu
z = z0 + phi
c = (-1)**modulo(tau+nu+phi, 2)*Ex(x)*Ey(y)*Ez(z)
S = S + c * T
end do
end do
end do
end subroutine auto2e_RCopy_KetTransform_12_1
end module auto2e_KetTransform_12_1