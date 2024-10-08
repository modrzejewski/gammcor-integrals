module auto2e_KetTransform_8_5
use arithmetic
use math_constants
implicit none

integer, dimension(45, 21), parameter :: RCopyIdx_8_5 = reshape([1, 15, 28, 40, 51, &
    61, 70, 78, 85, 106, 119, 131, 142, 152, 161, 169, 176, 197, 209, 220, 230, 239, &
    247, 254, 275, 286, 296, 305, 313, 320, 341, 351, 360, 368, 375, 396, 405, 413, &
    420, 441, 449, 456, 477, 484, 505, 15, 28, 40, 51, 61, 70, 78, 85, 91, 119, 131, &
    142, 152, 161, 169, 176, 182, 209, 220, 230, 239, 247, 254, 260, 286, 296, 305, &
    313, 320, 326, 351, 360, 368, 375, 381, 405, 413, 420, 426, 449, 456, 462, 484, &
    490, 511, 28, 40, 51, 61, 70, 78, 85, 91, 96, 131, 142, 152, 161, 169, 176, 182, &
    187, 220, 230, 239, 247, 254, 260, 265, 296, 305, 313, 320, 326, 331, 360, 368, &
    375, 381, 386, 413, 420, 426, 431, 456, 462, 467, 490, 495, 516, 40, 51, 61, &
    70, 78, 85, 91, 96, 100, 142, 152, 161, 169, 176, 182, 187, 191, 230, 239, 247, &
    254, 260, 265, 269, 305, 313, 320, 326, 331, 335, 368, 375, 381, 386, 390, 420, &
    426, 431, 435, 462, 467, 471, 495, 499, 520, 51, 61, 70, 78, 85, 91, 96, 100, &
    103, 152, 161, 169, 176, 182, 187, 191, 194, 239, 247, 254, 260, 265, 269, 272, &
    313, 320, 326, 331, 335, 338, 375, 381, 386, 390, 393, 426, 431, 435, 438, 467, &
    471, 474, 499, 502, 523, 61, 70, 78, 85, 91, 96, 100, 103, 105, 161, 169, 176, &
    182, 187, 191, 194, 196, 247, 254, 260, 265, 269, 272, 274, 320, 326, 331, 335, &
    338, 340, 381, 386, 390, 393, 395, 431, 435, 438, 440, 471, 474, 476, 502, 504, &
    525, 106, 119, 131, 142, 152, 161, 169, 176, 182, 197, 209, 220, 230, 239, 247, &
    254, 260, 275, 286, 296, 305, 313, 320, 326, 341, 351, 360, 368, 375, 381, 396, &
    405, 413, 420, 426, 441, 449, 456, 462, 477, 484, 490, 505, 511, 526, 119, 131, &
    142, 152, 161, 169, 176, 182, 187, 209, 220, 230, 239, 247, 254, 260, 265, 286, &
    296, 305, 313, 320, 326, 331, 351, 360, 368, 375, 381, 386, 405, 413, 420, 426, &
    431, 449, 456, 462, 467, 484, 490, 495, 511, 516, 531, 131, 142, 152, 161, 169, &
    176, 182, 187, 191, 220, 230, 239, 247, 254, 260, 265, 269, 296, 305, 313, 320, &
    326, 331, 335, 360, 368, 375, 381, 386, 390, 413, 420, 426, 431, 435, 456, 462, &
    467, 471, 490, 495, 499, 516, 520, 535, 142, 152, 161, 169, 176, 182, 187, 191, &
    194, 230, 239, 247, 254, 260, 265, 269, 272, 305, 313, 320, 326, 331, 335, 338, &
    368, 375, 381, 386, 390, 393, 420, 426, 431, 435, 438, 462, 467, 471, 474, 495, &
    499, 502, 520, 523, 538, 152, 161, 169, 176, 182, 187, 191, 194, 196, 239, 247, &
    254, 260, 265, 269, 272, 274, 313, 320, 326, 331, 335, 338, 340, 375, 381, 386, &
    390, 393, 395, 426, 431, 435, 438, 440, 467, 471, 474, 476, 499, 502, 504, 523, &
    525, 540, 197, 209, 220, 230, 239, 247, 254, 260, 265, 275, 286, 296, 305, 313, &
    320, 326, 331, 341, 351, 360, 368, 375, 381, 386, 396, 405, 413, 420, 426, 431, &
    441, 449, 456, 462, 467, 477, 484, 490, 495, 505, 511, 516, 526, 531, 541, 209, &
    220, 230, 239, 247, 254, 260, 265, 269, 286, 296, 305, 313, 320, 326, 331, 335, &
    351, 360, 368, 375, 381, 386, 390, 405, 413, 420, 426, 431, 435, 449, 456, 462, &
    467, 471, 484, 490, 495, 499, 511, 516, 520, 531, 535, 545, 220, 230, 239, 247, &
    254, 260, 265, 269, 272, 296, 305, 313, 320, 326, 331, 335, 338, 360, 368, 375, &
    381, 386, 390, 393, 413, 420, 426, 431, 435, 438, 456, 462, 467, 471, 474, 490, &
    495, 499, 502, 516, 520, 523, 535, 538, 548, 230, 239, 247, 254, 260, 265, 269, &
    272, 274, 305, 313, 320, 326, 331, 335, 338, 340, 368, 375, 381, 386, 390, 393, &
    395, 420, 426, 431, 435, 438, 440, 462, 467, 471, 474, 476, 495, 499, 502, 504, &
    520, 523, 525, 538, 540, 550, 275, 286, 296, 305, 313, 320, 326, 331, 335, 341, &
    351, 360, 368, 375, 381, 386, 390, 396, 405, 413, 420, 426, 431, 435, 441, 449, &
    456, 462, 467, 471, 477, 484, 490, 495, 499, 505, 511, 516, 520, 526, 531, 535, &
    541, 545, 551, 286, 296, 305, 313, 320, 326, 331, 335, 338, 351, 360, 368, 375, &
    381, 386, 390, 393, 405, 413, 420, 426, 431, 435, 438, 449, 456, 462, 467, 471, &
    474, 484, 490, 495, 499, 502, 511, 516, 520, 523, 531, 535, 538, 545, 548, 554, &
    296, 305, 313, 320, 326, 331, 335, 338, 340, 360, 368, 375, 381, 386, 390, 393, &
    395, 413, 420, 426, 431, 435, 438, 440, 456, 462, 467, 471, 474, 476, 490, 495, &
    499, 502, 504, 516, 520, 523, 525, 535, 538, 540, 548, 550, 556, 341, 351, 360, &
    368, 375, 381, 386, 390, 393, 396, 405, 413, 420, 426, 431, 435, 438, 441, 449, &
    456, 462, 467, 471, 474, 477, 484, 490, 495, 499, 502, 505, 511, 516, 520, 523, &
    526, 531, 535, 538, 541, 545, 548, 551, 554, 557, 351, 360, 368, 375, 381, 386, &
    390, 393, 395, 405, 413, 420, 426, 431, 435, 438, 440, 449, 456, 462, 467, 471, &
    474, 476, 484, 490, 495, 499, 502, 504, 511, 516, 520, 523, 525, 531, 535, 538, &
    540, 545, 548, 550, 554, 556, 559, 396, 405, 413, 420, 426, 431, 435, 438, 440, &
    441, 449, 456, 462, 467, 471, 474, 476, 477, 484, 490, 495, 499, 502, 504, 505, &
    511, 516, 520, 523, 525, 526, 531, 535, 538, 540, 541, 545, 548, 550, 551, 554, &
    556, 557, 559, 560], [45, 21])

contains

subroutine auto2e_RCopy_8_5(T, R, idx)
real(F64), dimension(:), intent(out) :: T
real(F64), dimension(:), intent(in) :: R
integer, dimension(:), intent(in) :: idx
T(1:9) = R(idx(1):idx(1)+8)
T(10:17) = R(idx(2):idx(2)+7)
T(18:24) = R(idx(3):idx(3)+6)
T(25:30) = R(idx(4):idx(4)+5)
T(31:35) = R(idx(5):idx(5)+4)
T(36:39) = R(idx(6):idx(6)+3)
T(40:42) = R(idx(7):idx(7)+2)
T(43:44) = R(idx(8):idx(8)+1)
T(45) = R(idx(9))
T(46:53) = R(idx(10):idx(10)+7)
T(54:60) = R(idx(11):idx(11)+6)
T(61:66) = R(idx(12):idx(12)+5)
T(67:71) = R(idx(13):idx(13)+4)
T(72:75) = R(idx(14):idx(14)+3)
T(76:78) = R(idx(15):idx(15)+2)
T(79:80) = R(idx(16):idx(16)+1)
T(81) = R(idx(17))
T(82:88) = R(idx(18):idx(18)+6)
T(89:94) = R(idx(19):idx(19)+5)
T(95:99) = R(idx(20):idx(20)+4)
T(100:103) = R(idx(21):idx(21)+3)
T(104:106) = R(idx(22):idx(22)+2)
T(107:108) = R(idx(23):idx(23)+1)
T(109) = R(idx(24))
T(110:115) = R(idx(25):idx(25)+5)
T(116:120) = R(idx(26):idx(26)+4)
T(121:124) = R(idx(27):idx(27)+3)
T(125:127) = R(idx(28):idx(28)+2)
T(128:129) = R(idx(29):idx(29)+1)
T(130) = R(idx(30))
T(131:135) = R(idx(31):idx(31)+4)
T(136:139) = R(idx(32):idx(32)+3)
T(140:142) = R(idx(33):idx(33)+2)
T(143:144) = R(idx(34):idx(34)+1)
T(145) = R(idx(35))
T(146:149) = R(idx(36):idx(36)+3)
T(150:152) = R(idx(37):idx(37)+2)
T(153:154) = R(idx(38):idx(38)+1)
T(155) = R(idx(39))
T(156:158) = R(idx(40):idx(40)+2)
T(159:160) = R(idx(41):idx(41)+1)
T(161) = R(idx(42))
T(162:163) = R(idx(43):idx(43)+1)
T(164) = R(idx(44))
T(165) = R(idx(45))
end subroutine auto2e_RCopy_8_5

subroutine auto2e_RCopy_KetTransform_8_5(S, T, Ex, Ey, Ez, R, lx, ly, lz)
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
i = ((2*5+3-phi)*phi)/2+nu+1
call auto2e_Rcopy_8_5(T, R(tau+1:), RCopyIdx_8_5(:, i))
x = x0 + tau
y = y0 + nu
z = z0 + phi
c = (-1)**modulo(tau+nu+phi, 2)*Ex(x)*Ey(y)*Ez(z)
S = S + c * T
end do
end do
end do
end subroutine auto2e_RCopy_KetTransform_8_5
end module auto2e_KetTransform_8_5