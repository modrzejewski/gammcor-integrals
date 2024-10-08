module auto2e_KetTransform_7_4
use arithmetic
use math_constants
implicit none

integer, dimension(36, 15), parameter :: RCopyIdx_7_4 = reshape([1, 13, 24, 34, 43, &
    51, 58, 64, 79, 90, 100, 109, 117, 124, 130, 145, 155, 164, 172, 179, 185, 200, &
    209, 217, 224, 230, 245, 253, 260, 266, 281, 288, 294, 309, 315, 330, 13, 24, &
    34, 43, 51, 58, 64, 69, 90, 100, 109, 117, 124, 130, 135, 155, 164, 172, 179, &
    185, 190, 209, 217, 224, 230, 235, 253, 260, 266, 271, 288, 294, 299, 315, 320, &
    335, 24, 34, 43, 51, 58, 64, 69, 73, 100, 109, 117, 124, 130, 135, 139, 164, &
    172, 179, 185, 190, 194, 217, 224, 230, 235, 239, 260, 266, 271, 275, 294, 299, &
    303, 320, 324, 339, 34, 43, 51, 58, 64, 69, 73, 76, 109, 117, 124, 130, 135, &
    139, 142, 172, 179, 185, 190, 194, 197, 224, 230, 235, 239, 242, 266, 271, 275, &
    278, 299, 303, 306, 324, 327, 342, 43, 51, 58, 64, 69, 73, 76, 78, 117, 124, &
    130, 135, 139, 142, 144, 179, 185, 190, 194, 197, 199, 230, 235, 239, 242, 244, &
    271, 275, 278, 280, 303, 306, 308, 327, 329, 344, 79, 90, 100, 109, 117, 124, &
    130, 135, 145, 155, 164, 172, 179, 185, 190, 200, 209, 217, 224, 230, 235, 245, &
    253, 260, 266, 271, 281, 288, 294, 299, 309, 315, 320, 330, 335, 345, 90, 100, &
    109, 117, 124, 130, 135, 139, 155, 164, 172, 179, 185, 190, 194, 209, 217, 224, &
    230, 235, 239, 253, 260, 266, 271, 275, 288, 294, 299, 303, 315, 320, 324, 335, &
    339, 349, 100, 109, 117, 124, 130, 135, 139, 142, 164, 172, 179, 185, 190, 194, &
    197, 217, 224, 230, 235, 239, 242, 260, 266, 271, 275, 278, 294, 299, 303, 306, &
    320, 324, 327, 339, 342, 352, 109, 117, 124, 130, 135, 139, 142, 144, 172, 179, &
    185, 190, 194, 197, 199, 224, 230, 235, 239, 242, 244, 266, 271, 275, 278, 280, &
    299, 303, 306, 308, 324, 327, 329, 342, 344, 354, 145, 155, 164, 172, 179, 185, &
    190, 194, 200, 209, 217, 224, 230, 235, 239, 245, 253, 260, 266, 271, 275, 281, &
    288, 294, 299, 303, 309, 315, 320, 324, 330, 335, 339, 345, 349, 355, 155, 164, &
    172, 179, 185, 190, 194, 197, 209, 217, 224, 230, 235, 239, 242, 253, 260, 266, &
    271, 275, 278, 288, 294, 299, 303, 306, 315, 320, 324, 327, 335, 339, 342, 349, &
    352, 358, 164, 172, 179, 185, 190, 194, 197, 199, 217, 224, 230, 235, 239, 242, &
    244, 260, 266, 271, 275, 278, 280, 294, 299, 303, 306, 308, 320, 324, 327, 329, &
    339, 342, 344, 352, 354, 360, 200, 209, 217, 224, 230, 235, 239, 242, 245, 253, &
    260, 266, 271, 275, 278, 281, 288, 294, 299, 303, 306, 309, 315, 320, 324, 327, &
    330, 335, 339, 342, 345, 349, 352, 355, 358, 361, 209, 217, 224, 230, 235, 239, &
    242, 244, 253, 260, 266, 271, 275, 278, 280, 288, 294, 299, 303, 306, 308, 315, &
    320, 324, 327, 329, 335, 339, 342, 344, 349, 352, 354, 358, 360, 363, 245, 253, &
    260, 266, 271, 275, 278, 280, 281, 288, 294, 299, 303, 306, 308, 309, 315, 320, &
    324, 327, 329, 330, 335, 339, 342, 344, 345, 349, 352, 354, 355, 358, 360, 361, &
    363, 364], [36, 15])

contains

subroutine auto2e_RCopy_7_4(T, R, idx)
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
end subroutine auto2e_RCopy_7_4

subroutine auto2e_RCopy_KetTransform_7_4(S, T, Ex, Ey, Ez, R, lx, ly, lz)
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
i = ((2*4+3-phi)*phi)/2+nu+1
call auto2e_Rcopy_7_4(T, R(tau+1:), RCopyIdx_7_4(:, i))
x = x0 + tau
y = y0 + nu
z = z0 + phi
c = (-1)**modulo(tau+nu+phi, 2)*Ex(x)*Ey(y)*Ez(z)
S = S + c * T
end do
end do
end do
end subroutine auto2e_RCopy_KetTransform_7_4
end module auto2e_KetTransform_7_4