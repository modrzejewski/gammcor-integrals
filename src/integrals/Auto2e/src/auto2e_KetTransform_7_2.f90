module auto2e_KetTransform_7_2
use arithmetic
use math_constants
implicit none

integer, dimension(36, 6), parameter :: RCopyIdx_7_2 = reshape([1, 11, 20, 28, 35, &
    41, 46, 50, 56, 65, 73, 80, 86, 91, 95, 101, 109, 116, 122, 127, 131, 137, 144, &
    150, 155, 159, 165, 171, 176, 180, 186, 191, 195, 201, 205, 211, 11, 20, 28, &
    35, 41, 46, 50, 53, 65, 73, 80, 86, 91, 95, 98, 109, 116, 122, 127, 131, 134, &
    144, 150, 155, 159, 162, 171, 176, 180, 183, 191, 195, 198, 205, 208, 214, 20, &
    28, 35, 41, 46, 50, 53, 55, 73, 80, 86, 91, 95, 98, 100, 116, 122, 127, 131, &
    134, 136, 150, 155, 159, 162, 164, 176, 180, 183, 185, 195, 198, 200, 208, 210, &
    216, 56, 65, 73, 80, 86, 91, 95, 98, 101, 109, 116, 122, 127, 131, 134, 137, &
    144, 150, 155, 159, 162, 165, 171, 176, 180, 183, 186, 191, 195, 198, 201, 205, &
    208, 211, 214, 217, 65, 73, 80, 86, 91, 95, 98, 100, 109, 116, 122, 127, 131, &
    134, 136, 144, 150, 155, 159, 162, 164, 171, 176, 180, 183, 185, 191, 195, 198, &
    200, 205, 208, 210, 214, 216, 219, 101, 109, 116, 122, 127, 131, 134, 136, 137, &
    144, 150, 155, 159, 162, 164, 165, 171, 176, 180, 183, 185, 186, 191, 195, 198, &
    200, 201, 205, 208, 210, 211, 214, 216, 217, 219, 220], [36, 6])

contains

subroutine auto2e_RCopy_7_2(T, R, idx)
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
end subroutine auto2e_RCopy_7_2

subroutine auto2e_RCopy_KetTransform_7_2(S, T, Ex, Ey, Ez, R, lx, ly, lz)
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
call auto2e_Rcopy_7_2(T, R(tau+1:), RCopyIdx_7_2(:, i))
x = x0 + tau
y = y0 + nu
z = z0 + phi
c = (-1)**modulo(tau+nu+phi, 2)*Ex(x)*Ey(y)*Ez(z)
S = S + c * T
end do
end do
end do
end subroutine auto2e_RCopy_KetTransform_7_2
end module auto2e_KetTransform_7_2