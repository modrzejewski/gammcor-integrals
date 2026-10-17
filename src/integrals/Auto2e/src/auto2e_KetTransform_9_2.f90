module auto2e_KetTransform_9_2
use arithmetic
use math_constants
implicit none

integer, dimension(55, 6), parameter :: RCopyIdx_9_2 = reshape([1, 13, 24, 34, 43, &
    51, 58, 64, 69, 73, 79, 90, 100, 109, 117, 124, 130, 135, 139, 145, 155, 164, &
    172, 179, 185, 190, 194, 200, 209, 217, 224, 230, 235, 239, 245, 253, 260, 266, &
    271, 275, 281, 288, 294, 299, 303, 309, 315, 320, 324, 330, 335, 339, 345, 349, &
    355, 13, 24, 34, 43, 51, 58, 64, 69, 73, 76, 90, 100, 109, 117, 124, 130, 135, &
    139, 142, 155, 164, 172, 179, 185, 190, 194, 197, 209, 217, 224, 230, 235, 239, &
    242, 253, 260, 266, 271, 275, 278, 288, 294, 299, 303, 306, 315, 320, 324, 327, &
    335, 339, 342, 349, 352, 358, 24, 34, 43, 51, 58, 64, 69, 73, 76, 78, 100, 109, &
    117, 124, 130, 135, 139, 142, 144, 164, 172, 179, 185, 190, 194, 197, 199, 217, &
    224, 230, 235, 239, 242, 244, 260, 266, 271, 275, 278, 280, 294, 299, 303, 306, &
    308, 320, 324, 327, 329, 339, 342, 344, 352, 354, 360, 79, 90, 100, 109, 117, &
    124, 130, 135, 139, 142, 145, 155, 164, 172, 179, 185, 190, 194, 197, 200, 209, &
    217, 224, 230, 235, 239, 242, 245, 253, 260, 266, 271, 275, 278, 281, 288, 294, &
    299, 303, 306, 309, 315, 320, 324, 327, 330, 335, 339, 342, 345, 349, 352, 355, &
    358, 361, 90, 100, 109, 117, 124, 130, 135, 139, 142, 144, 155, 164, 172, 179, &
    185, 190, 194, 197, 199, 209, 217, 224, 230, 235, 239, 242, 244, 253, 260, 266, &
    271, 275, 278, 280, 288, 294, 299, 303, 306, 308, 315, 320, 324, 327, 329, 335, &
    339, 342, 344, 349, 352, 354, 358, 360, 363, 145, 155, 164, 172, 179, 185, 190, &
    194, 197, 199, 200, 209, 217, 224, 230, 235, 239, 242, 244, 245, 253, 260, 266, &
    271, 275, 278, 280, 281, 288, 294, 299, 303, 306, 308, 309, 315, 320, 324, 327, &
    329, 330, 335, 339, 342, 344, 345, 349, 352, 354, 355, 358, 360, 361, 363, 364], &
    [55, 6])

contains

subroutine auto2e_RCopy_9_2(T, R, idx)
real(F64), dimension(:), intent(out) :: T
real(F64), dimension(:), intent(in) :: R
integer, dimension(:), intent(in) :: idx
T(1:10) = R(idx(1):idx(1)+9)
T(11:19) = R(idx(2):idx(2)+8)
T(20:27) = R(idx(3):idx(3)+7)
T(28:34) = R(idx(4):idx(4)+6)
T(35:40) = R(idx(5):idx(5)+5)
T(41:45) = R(idx(6):idx(6)+4)
T(46:49) = R(idx(7):idx(7)+3)
T(50:52) = R(idx(8):idx(8)+2)
T(53:54) = R(idx(9):idx(9)+1)
T(55) = R(idx(10))
T(56:64) = R(idx(11):idx(11)+8)
T(65:72) = R(idx(12):idx(12)+7)
T(73:79) = R(idx(13):idx(13)+6)
T(80:85) = R(idx(14):idx(14)+5)
T(86:90) = R(idx(15):idx(15)+4)
T(91:94) = R(idx(16):idx(16)+3)
T(95:97) = R(idx(17):idx(17)+2)
T(98:99) = R(idx(18):idx(18)+1)
T(100) = R(idx(19))
T(101:108) = R(idx(20):idx(20)+7)
T(109:115) = R(idx(21):idx(21)+6)
T(116:121) = R(idx(22):idx(22)+5)
T(122:126) = R(idx(23):idx(23)+4)
T(127:130) = R(idx(24):idx(24)+3)
T(131:133) = R(idx(25):idx(25)+2)
T(134:135) = R(idx(26):idx(26)+1)
T(136) = R(idx(27))
T(137:143) = R(idx(28):idx(28)+6)
T(144:149) = R(idx(29):idx(29)+5)
T(150:154) = R(idx(30):idx(30)+4)
T(155:158) = R(idx(31):idx(31)+3)
T(159:161) = R(idx(32):idx(32)+2)
T(162:163) = R(idx(33):idx(33)+1)
T(164) = R(idx(34))
T(165:170) = R(idx(35):idx(35)+5)
T(171:175) = R(idx(36):idx(36)+4)
T(176:179) = R(idx(37):idx(37)+3)
T(180:182) = R(idx(38):idx(38)+2)
T(183:184) = R(idx(39):idx(39)+1)
T(185) = R(idx(40))
T(186:190) = R(idx(41):idx(41)+4)
T(191:194) = R(idx(42):idx(42)+3)
T(195:197) = R(idx(43):idx(43)+2)
T(198:199) = R(idx(44):idx(44)+1)
T(200) = R(idx(45))
T(201:204) = R(idx(46):idx(46)+3)
T(205:207) = R(idx(47):idx(47)+2)
T(208:209) = R(idx(48):idx(48)+1)
T(210) = R(idx(49))
T(211:213) = R(idx(50):idx(50)+2)
T(214:215) = R(idx(51):idx(51)+1)
T(216) = R(idx(52))
T(217:218) = R(idx(53):idx(53)+1)
T(219) = R(idx(54))
T(220) = R(idx(55))
end subroutine auto2e_RCopy_9_2

subroutine auto2e_RCopy_KetTransform_9_2(S, T, Ex, Ey, Ez, R, lx, ly, lz)
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
call auto2e_Rcopy_9_2(T, R(tau+1:), RCopyIdx_9_2(:, i))
x = x0 + tau
y = y0 + nu
z = z0 + phi
c = (-1)**modulo(tau+nu+phi, 2)*Ex(x)*Ey(y)*Ez(z)
S = S + c * T
end do
end do
end do
end subroutine auto2e_RCopy_KetTransform_9_2
end module auto2e_KetTransform_9_2