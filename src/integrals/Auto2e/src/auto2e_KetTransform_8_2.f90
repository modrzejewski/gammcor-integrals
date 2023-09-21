module auto2e_KetTransform_8_2
use arithmetic
implicit none

integer, dimension(45, 6), parameter :: RCopyIdx_8_2 = reshape([1, 12, 22, 31, 39, &
    46, 52, 57, 61, 67, 77, 86, 94, 101, 107, 112, 116, 122, 131, 139, 146, 152, &
    157, 161, 167, 175, 182, 188, 193, 197, 203, 210, 216, 221, 225, 231, 237, 242, &
    246, 252, 257, 261, 267, 271, 277, 12, 22, 31, 39, 46, 52, 57, 61, 64, 77, 86, &
    94, 101, 107, 112, 116, 119, 131, 139, 146, 152, 157, 161, 164, 175, 182, 188, &
    193, 197, 200, 210, 216, 221, 225, 228, 237, 242, 246, 249, 257, 261, 264, 271, &
    274, 280, 22, 31, 39, 46, 52, 57, 61, 64, 66, 86, 94, 101, 107, 112, 116, 119, &
    121, 139, 146, 152, 157, 161, 164, 166, 182, 188, 193, 197, 200, 202, 216, 221, &
    225, 228, 230, 242, 246, 249, 251, 261, 264, 266, 274, 276, 282, 67, 77, 86, &
    94, 101, 107, 112, 116, 119, 122, 131, 139, 146, 152, 157, 161, 164, 167, 175, &
    182, 188, 193, 197, 200, 203, 210, 216, 221, 225, 228, 231, 237, 242, 246, 249, &
    252, 257, 261, 264, 267, 271, 274, 277, 280, 283, 77, 86, 94, 101, 107, 112, &
    116, 119, 121, 131, 139, 146, 152, 157, 161, 164, 166, 175, 182, 188, 193, 197, &
    200, 202, 210, 216, 221, 225, 228, 230, 237, 242, 246, 249, 251, 257, 261, 264, &
    266, 271, 274, 276, 280, 282, 285, 122, 131, 139, 146, 152, 157, 161, 164, 166, &
    167, 175, 182, 188, 193, 197, 200, 202, 203, 210, 216, 221, 225, 228, 230, 231, &
    237, 242, 246, 249, 251, 252, 257, 261, 264, 266, 267, 271, 274, 276, 277, 280, &
    282, 283, 285, 286], [45, 6])

contains

subroutine auto2e_RCopy_8_2(T, R, idx)
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
end subroutine auto2e_RCopy_8_2

subroutine auto2e_RCopy_KetTransform_8_2(S, T, Ex, Ey, Ez, R, lx, ly, lz)
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
call auto2e_Rcopy_8_2(T, R(tau+1:), RCopyIdx_8_2(:, i))
x = x0 + tau
y = y0 + nu
z = z0 + phi
c = (-1)**modulo(tau+nu+phi, 2)*Ex(x)*Ey(y)*Ez(z)
S = S + c * T
end do
end do
end do
end subroutine auto2e_RCopy_KetTransform_8_2
end module auto2e_KetTransform_8_2