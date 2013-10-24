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
        real*8 :: FlapFractionRoot = 0.25d0 ! Flap fraction at aileron root (cf/c)
        real*8 :: FlapFractionTip = 0.25d0 ! Flap fraction at aileron tip (cf/c)
        logical :: ParallelHingeLine = .true. ! Is the hinge line parallel to the
                                              ! quarter-chord line? When true,
                                              ! FlapFractionTip will be calculated

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
        function GetWingType(pf) result(name)
            type(Planform), intent(in) :: pf
            character*8 :: name

            if (pf%WingType .eq. Tapered) then
                name = "Tapered"
            else if (pf%WingType .eq. Elliptic) then
                name = "Elliptic"
            else
                name = "Unknown"
            end if
        end function GetWingType

        real*8 function theta_i(i, nnodes) result(theta)
            integer, intent(in) :: i
            integer, intent(in) :: nnodes

            if (i < 1 .or. i > nnodes) then
                write(6, '(a, f7.4)') "ERROR: Function theta_i called with i = ", i
                if (i < 1) then
                    theta = 0.0d0
                else
                    theta = pi
                end if
            else
                theta = real(i-1, 8) / real(nnodes - 1, 8) * pi
            end if
        end function theta_i

        real*8 function theta_zb(zb) result(theta)
            real*8, intent(in) :: zb  ! z/b

            if (zb < -0.5d0 .or. zb > 0.5d0) then
                write(6, '(a, f7.4)') "ERROR: Function theta_d called with z/b = ", zb
                if (zb < -0.5d0) then
                    theta = 0.0d0
                else
                    theta = pi
                end if
            else
                theta = acos(-2.0d0 * zb)
            end if
        end function theta_zb

        real*8 function c_over_b_i(pf, i) result(cb)
            type(Planform), intent(in) :: pf
            integer, intent(in) :: i

            real*8 :: theta

            theta = theta_i(i, pf%NNodes)
            cb = c_over_b(pf, theta)
        end function c_over_b_i

        real*8 function c_over_b_zb(pf, zb) result(cb)
            type(Planform), intent(in) :: pf
            real*8, intent(in) :: zb  ! z/b

            real*8 :: theta

            theta = theta_zb(zb)
            cb = c_over_b(pf, theta)
        end function c_over_b_zb

        real*8 function c_over_b(pf, theta) result(cb)
            type(Planform), intent(in) :: pf
            real*8, intent(in) :: theta

            if (pf%WingType == Tapered) then
                ! Calculate c/b for tapered wing
                cb = (2.0d0 * (1.0d0 - (1.0d0 - pf%TaperRatio) * &
                    & dabs(cos(theta)))) / (pf%AspectRatio * (1.0d0 + pf%TaperRatio))
            else if (pf%WingType == Elliptic) then
                ! Calculate c/b for elliptic wing
                cb = (4.0d0 * sin(theta)) / &
                    & (pi * pf%AspectRatio)
            else
                ! Unknown wing type!
                stop "*** Unknown Wing Type ***"
            end if
        end function c_over_b

        real*8 function z_over_b_i(i, nnodes) result(zb)
            integer, intent(in) :: i
            integer, intent(in) :: nnodes

            zb = z_over_b(theta_i(i, nnodes))
        end function z_over_b_i

        real*8 function z_over_b(theta) result(zb)
            real*8, intent(in) :: theta

            zb = -0.5d0 * cos(theta)
        end function z_over_b

        subroutine CalculateAileronTipFlapFraction(pf)
            type(Planform), intent(inout) :: pf

            real*8 :: cb_root, cfc_root
            real*8 :: cb_tip, cfc_tip

            if (pf%ParallelHingeLine) then
                cb_root = c_over_b_zb(pf, pf%AileronRoot)
                cfc_root = pf%FlapFractionRoot
                cb_tip = c_over_b_zb(pf, pf%AileronTip)
                cfc_tip = 0.75d0 - cb_root / cb_tip * (0.75d0 - cfc_root)
                pf%FlapFractionTip = cfc_tip
            end if
        end subroutine CalculateAileronTipFlapFraction

end module class_Planform


