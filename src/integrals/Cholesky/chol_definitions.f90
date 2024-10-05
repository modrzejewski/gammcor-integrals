module chol_definitions
      use arithmetic
      use display

      implicit none

      integer, parameter :: CHOL_ACCURACY_DEFAULT = 1
      integer, parameter :: CHOL_ACCURACY_TIGHT = 2
      integer, parameter :: CHOL_ACCURACY_LUDICROUS = 3
      integer, parameter :: CHOL_ACCURACY_DEBUG = 4
      
contains

      subroutine chol_NumericalThresholds(TargetTraceError, TargetTraceErrorPrescreen, &
            TargetMaxError, MaxNAOMult, Accuracy)
            
            real(F64), intent(out) :: TargetTraceError
            real(F64), intent(out) :: TargetTraceErrorPrescreen
            real(F64), intent(out) :: TargetMaxError
            integer, intent(out)   :: MaxNAOMult
            integer, intent(in)    :: Accuracy
            
            select case (Accuracy)
            case (CHOL_ACCURACY_DEFAULT)
                  TargetTraceError = 1.0E-2_F64
                  TargetTraceErrorPrescreen = 1.0E-14_F64
                  TargetMaxError = 1.0E-11_F64
                  MaxNAOMult = 8
            case (CHOL_ACCURACY_TIGHT)
                  TargetTraceError = 1.0E-3_F64
                  TargetTraceErrorPrescreen = 1.0E-14_F64
                  TargetMaxError = 1.0E-11_F64
                  MaxNAOMult = 9
            case (CHOL_ACCURACY_LUDICROUS)
                  TargetTraceError = 1.0E-4_F64
                  TargetTraceErrorPrescreen = 1.0E-14_F64
                  TargetMaxError = 1.0E-11_F64
                  MaxNAOMult = 10
            case (CHOL_ACCURACY_DEBUG)
                  call msg("Using debug accuracy, can be numerically unstable")
                  TargetTraceError = 1.0E-10_F64
                  TargetTraceErrorPrescreen = 1.0E-18_F64
                  TargetMaxError = 1.0E-11_F64
                  MaxNAOMult = 20
            end select
      end subroutine chol_NumericalThresholds
end module chol_definitions
