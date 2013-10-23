module class_Planform
    implicit none

    public :: Planform

    real*8, parameter :: pi = acos(-1.0d0)

    ! Supported wing types
    enum, bind(C)
        enumerator :: Tapered = 1, Elliptic = 2
    end enum

    type Planform
        ! Wing Parameters
        integer :: WingType = Tapered ! Wing type
        integer :: NNodes = 7 ! Total number of nodes
        real*8 :: AspectRatio = 8.0d0 ! Aspect ratio
        real*8 :: TaperRatio = 1.0d0 ! Taper ratio (tapered wing only)
        real*8 :: LiftSlope = 2.0d0 * pi ! Section lift slope
        real*8 :: AileronRoot = 0.25d0 ! Location of aileron root (z/b)
        real*8 :: AileronTip = 0.45d0 ! Location of aileron tip (z/b)

        ! Output Options
        logical :: WriteCMatrix = .true.  ! Write C Matrix to output file?
        logical :: WriteCInverse = .true.  ! Write Inv(C) Matrix to output file?
        logical :: WriteFourier = .true.  ! Write Fourier coefficients?
        logical :: WriteOther = .true.  ! Write KL, KD, eps, and lift slope?
        character*80 :: FileName = "PrandtlsLiftingLine.out" ! Name of output file

        ! Operating Conditions
        real*8 :: DesiredAngleOfAttack = pi / 36.0d0 ! Desired root aerodynamic angle of Attack
                                                     ! (alpha - alpha_L0), in radians
                                                     ! When specified, a new LiftCoefficient is calculated
        real*8 :: AngleOfAttack = pi / 36.0d0 ! Root Aerodynamic Angle of Attack
                                              ! (alpha - alpha_L0), in radians
        real*8 :: DesiredLiftCoefficient = 0.4d0 ! Desired lift coefficient
                                                 ! When specified, a new AngleOfAttack is calculated
        real*8 :: LiftCoefficient = 0.4d0 ! Lift coefficient
        real*8 :: Omega = 0.0d0 ! Amount of linear twist, in radians
        real*8 :: AileronDeflection = 0.0d0 ! Aileron deflection, in radians
        real*8 :: RollingRate = 0.0d0 ! Dimensionless rolling rate (constant over wingspan)
        real*8 :: FlapEffectiveness = 0.445d0 ! Section flap effectiveness (constant over wingspan)
        logical :: SpecifyAlpha = .true. ! Was alpha specified?
                                         ! .true.  = Use desired alpha to calculate CL
                                         ! .false. = Use desired CL to calculate alpha

    end type Planform

    contains
        function GetWingType(this) result(name)
            type(Planform) :: this
            character*8 :: name

            if (this%WingType .eq. Tapered) then
                name = "Tapered"
            else if (this%WingType .eq. Elliptic) then
                name = "Elliptic"
            else
                name = "Unknown"
            end if
        end function GetWingType

end module class_Planform


