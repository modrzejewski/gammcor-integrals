module auto2e_KetTransform_7_1
use arithmetic
use math_constants
implicit none

integer, dimension(36, 3), parameter :: RCopyIdx_7_1 = reshape([1, 10, 18, 25, 31, &
    36, 40, 43, 46, 54, 61, 67, 72, 76, 79, 82, 89, 95, 100, 104, 107, 110, 116, &
    121, 125, 128, 131, 136, 140, 143, 146, 150, 153, 156, 159, 162, 10, 18, 25, &
    31, 36, 40, 43, 45, 54, 61, 67, 72, 76, 79, 81, 89, 95, 100, 104, 107, 109, 116, &
    121, 125, 128, 130, 136, 140, 143, 145, 150, 153, 155, 159, 161, 164, 46, 54, &
    61, 67, 72, 76, 79, 81, 82, 89, 95, 100, 104, 107, 109, 110, 116, 121, 125, 128, &
    130, 131, 136, 140, 143, 145, 146, 150, 153, 155, 156, 159, 161, 162, 164, 165], &
    [36, 3])

contains

subroutine auto2e_RCopy_7_1(T, R, idx)
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
end subroutine auto2e_RCopy_7_1

subroutine auto2e_RCopy_KetTransform_7_1(S, T, Ex, Ey, Ez, R, lx, ly, lz)
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
call auto2e_Rcopy_7_1(T, R(tau+1:), RCopyIdx_7_1(:, i))
x = x0 + tau
y = y0 + nu
z = z0 + phi
c = (-1)**modulo(tau+nu+phi, 2)*Ex(x)*Ey(y)*Ez(z)
S = S + c * T
end do
end do
end do
end subroutine auto2e_RCopy_KetTransform_7_1
end module auto2e_KetTransform_7_1