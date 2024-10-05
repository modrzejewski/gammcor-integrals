module auto2e_KetTransform_7_6
use arithmetic
use math_constants
implicit none

integer, dimension(36, 28), parameter :: RCopyIdx_7_6 = reshape([1, 15, 28, 40, 51, &
    61, 70, 78, 106, 119, 131, 142, 152, 161, 169, 197, 209, 220, 230, 239, 247, &
    275, 286, 296, 305, 313, 341, 351, 360, 368, 396, 405, 413, 441, 449, 477, 15, &
    28, 40, 51, 61, 70, 78, 85, 119, 131, 142, 152, 161, 169, 176, 209, 220, 230, &
    239, 247, 254, 286, 296, 305, 313, 320, 351, 360, 368, 375, 405, 413, 420, 449, &
    456, 484, 28, 40, 51, 61, 70, 78, 85, 91, 131, 142, 152, 161, 169, 176, 182, &
    220, 230, 239, 247, 254, 260, 296, 305, 313, 320, 326, 360, 368, 375, 381, 413, &
    420, 426, 456, 462, 490, 40, 51, 61, 70, 78, 85, 91, 96, 142, 152, 161, 169, &
    176, 182, 187, 230, 239, 247, 254, 260, 265, 305, 313, 320, 326, 331, 368, 375, &
    381, 386, 420, 426, 431, 462, 467, 495, 51, 61, 70, 78, 85, 91, 96, 100, 152, &
    161, 169, 176, 182, 187, 191, 239, 247, 254, 260, 265, 269, 313, 320, 326, 331, &
    335, 375, 381, 386, 390, 426, 431, 435, 467, 471, 499, 61, 70, 78, 85, 91, 96, &
    100, 103, 161, 169, 176, 182, 187, 191, 194, 247, 254, 260, 265, 269, 272, 320, &
    326, 331, 335, 338, 381, 386, 390, 393, 431, 435, 438, 471, 474, 502, 70, 78, &
    85, 91, 96, 100, 103, 105, 169, 176, 182, 187, 191, 194, 196, 254, 260, 265, &
    269, 272, 274, 326, 331, 335, 338, 340, 386, 390, 393, 395, 435, 438, 440, 474, &
    476, 504, 106, 119, 131, 142, 152, 161, 169, 176, 197, 209, 220, 230, 239, 247, &
    254, 275, 286, 296, 305, 313, 320, 341, 351, 360, 368, 375, 396, 405, 413, 420, &
    441, 449, 456, 477, 484, 505, 119, 131, 142, 152, 161, 169, 176, 182, 209, 220, &
    230, 239, 247, 254, 260, 286, 296, 305, 313, 320, 326, 351, 360, 368, 375, 381, &
    405, 413, 420, 426, 449, 456, 462, 484, 490, 511, 131, 142, 152, 161, 169, 176, &
    182, 187, 220, 230, 239, 247, 254, 260, 265, 296, 305, 313, 320, 326, 331, 360, &
    368, 375, 381, 386, 413, 420, 426, 431, 456, 462, 467, 490, 495, 516, 142, 152, &
    161, 169, 176, 182, 187, 191, 230, 239, 247, 254, 260, 265, 269, 305, 313, 320, &
    326, 331, 335, 368, 375, 381, 386, 390, 420, 426, 431, 435, 462, 467, 471, 495, &
    499, 520, 152, 161, 169, 176, 182, 187, 191, 194, 239, 247, 254, 260, 265, 269, &
    272, 313, 320, 326, 331, 335, 338, 375, 381, 386, 390, 393, 426, 431, 435, 438, &
    467, 471, 474, 499, 502, 523, 161, 169, 176, 182, 187, 191, 194, 196, 247, 254, &
    260, 265, 269, 272, 274, 320, 326, 331, 335, 338, 340, 381, 386, 390, 393, 395, &
    431, 435, 438, 440, 471, 474, 476, 502, 504, 525, 197, 209, 220, 230, 239, 247, &
    254, 260, 275, 286, 296, 305, 313, 320, 326, 341, 351, 360, 368, 375, 381, 396, &
    405, 413, 420, 426, 441, 449, 456, 462, 477, 484, 490, 505, 511, 526, 209, 220, &
    230, 239, 247, 254, 260, 265, 286, 296, 305, 313, 320, 326, 331, 351, 360, 368, &
    375, 381, 386, 405, 413, 420, 426, 431, 449, 456, 462, 467, 484, 490, 495, 511, &
    516, 531, 220, 230, 239, 247, 254, 260, 265, 269, 296, 305, 313, 320, 326, 331, &
    335, 360, 368, 375, 381, 386, 390, 413, 420, 426, 431, 435, 456, 462, 467, 471, &
    490, 495, 499, 516, 520, 535, 230, 239, 247, 254, 260, 265, 269, 272, 305, 313, &
    320, 326, 331, 335, 338, 368, 375, 381, 386, 390, 393, 420, 426, 431, 435, 438, &
    462, 467, 471, 474, 495, 499, 502, 520, 523, 538, 239, 247, 254, 260, 265, 269, &
    272, 274, 313, 320, 326, 331, 335, 338, 340, 375, 381, 386, 390, 393, 395, 426, &
    431, 435, 438, 440, 467, 471, 474, 476, 499, 502, 504, 523, 525, 540, 275, 286, &
    296, 305, 313, 320, 326, 331, 341, 351, 360, 368, 375, 381, 386, 396, 405, 413, &
    420, 426, 431, 441, 449, 456, 462, 467, 477, 484, 490, 495, 505, 511, 516, 526, &
    531, 541, 286, 296, 305, 313, 320, 326, 331, 335, 351, 360, 368, 375, 381, 386, &
    390, 405, 413, 420, 426, 431, 435, 449, 456, 462, 467, 471, 484, 490, 495, 499, &
    511, 516, 520, 531, 535, 545, 296, 305, 313, 320, 326, 331, 335, 338, 360, 368, &
    375, 381, 386, 390, 393, 413, 420, 426, 431, 435, 438, 456, 462, 467, 471, 474, &
    490, 495, 499, 502, 516, 520, 523, 535, 538, 548, 305, 313, 320, 326, 331, 335, &
    338, 340, 368, 375, 381, 386, 390, 393, 395, 420, 426, 431, 435, 438, 440, 462, &
    467, 471, 474, 476, 495, 499, 502, 504, 520, 523, 525, 538, 540, 550, 341, 351, &
    360, 368, 375, 381, 386, 390, 396, 405, 413, 420, 426, 431, 435, 441, 449, 456, &
    462, 467, 471, 477, 484, 490, 495, 499, 505, 511, 516, 520, 526, 531, 535, 541, &
    545, 551, 351, 360, 368, 375, 381, 386, 390, 393, 405, 413, 420, 426, 431, 435, &
    438, 449, 456, 462, 467, 471, 474, 484, 490, 495, 499, 502, 511, 516, 520, 523, &
    531, 535, 538, 545, 548, 554, 360, 368, 375, 381, 386, 390, 393, 395, 413, 420, &
    426, 431, 435, 438, 440, 456, 462, 467, 471, 474, 476, 490, 495, 499, 502, 504, &
    516, 520, 523, 525, 535, 538, 540, 548, 550, 556, 396, 405, 413, 420, 426, 431, &
    435, 438, 441, 449, 456, 462, 467, 471, 474, 477, 484, 490, 495, 499, 502, 505, &
    511, 516, 520, 523, 526, 531, 535, 538, 541, 545, 548, 551, 554, 557, 405, 413, &
    420, 426, 431, 435, 438, 440, 449, 456, 462, 467, 471, 474, 476, 484, 490, 495, &
    499, 502, 504, 511, 516, 520, 523, 525, 531, 535, 538, 540, 545, 548, 550, 554, &
    556, 559, 441, 449, 456, 462, 467, 471, 474, 476, 477, 484, 490, 495, 499, 502, &
    504, 505, 511, 516, 520, 523, 525, 526, 531, 535, 538, 540, 541, 545, 548, 550, &
    551, 554, 556, 557, 559, 560], [36, 28])

contains

subroutine auto2e_RCopy_7_6(T, R, idx)
real(F64), dimension(:), intent(out) :: T
real(F64), dimension(:), intent(in) :: R
integer, dimension(:), intent(in) :: idx
T(1:8) = R(idx(1):idx(1)+7)
T(9:15) = R(idx(2):idx(2)+6)
T(16:21) = R(idx(3):idx(3)+5)
T(22:26) = R(idx(4):idx(4)+4)
T(27:30) = R(idx(5):idx(5)+3)
T(31:33) = R(idx(6):idx(6)+2)
T(34:35) = R(idx(7):idx(7)+1)
T(36) = R(idx(8))
T(37:43) = R(idx(9):idx(9)+6)
T(44:49) = R(idx(10):idx(10)+5)
T(50:54) = R(idx(11):idx(11)+4)
T(55:58) = R(idx(12):idx(12)+3)
T(59:61) = R(idx(13):idx(13)+2)
T(62:63) = R(idx(14):idx(14)+1)
T(64) = R(idx(15))
T(65:70) = R(idx(16):idx(16)+5)
T(71:75) = R(idx(17):idx(17)+4)
T(76:79) = R(idx(18):idx(18)+3)
T(80:82) = R(idx(19):idx(19)+2)
T(83:84) = R(idx(20):idx(20)+1)
T(85) = R(idx(21))
T(86:90) = R(idx(22):idx(22)+4)
T(91:94) = R(idx(23):idx(23)+3)
T(95:97) = R(idx(24):idx(24)+2)
T(98:99) = R(idx(25):idx(25)+1)
T(100) = R(idx(26))
T(101:104) = R(idx(27):idx(27)+3)
T(105:107) = R(idx(28):idx(28)+2)
T(108:109) = R(idx(29):idx(29)+1)
T(110) = R(idx(30))
T(111:113) = R(idx(31):idx(31)+2)
T(114:115) = R(idx(32):idx(32)+1)
T(116) = R(idx(33))
T(117:118) = R(idx(34):idx(34)+1)
T(119) = R(idx(35))
T(120) = R(idx(36))
end subroutine auto2e_RCopy_7_6

subroutine auto2e_RCopy_KetTransform_7_6(S, T, Ex, Ey, Ez, R, lx, ly, lz)
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
i = ((2*6+3-phi)*phi)/2+nu+1
call auto2e_Rcopy_7_6(T, R(tau+1:), RCopyIdx_7_6(:, i))
x = x0 + tau
y = y0 + nu
z = z0 + phi
c = (-1)**modulo(tau+nu+phi, 2)*Ex(x)*Ey(y)*Ez(z)
S = S + c * T
end do
end do
end do
end subroutine auto2e_RCopy_KetTransform_7_6
end module auto2e_KetTransform_7_6