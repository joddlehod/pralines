module class_Planform
    use Utilities
    implicit none

    public :: Planform

    ! Supported wing types
    enum, bind(C)
        enumerator :: Tapered = 1, Elliptic = 2, Combination = 3
    end enum

    ! Supported washout distribution types
    enum, bind(C)
        enumerator :: Linear = 1, Optimum = 2
    end enum

    ! Supported low-aspect-ratio methods
    enum, bind(C)
        enumerator :: Classical = 1, Hodson = 2, ModifiedSlender = 3, Kuchemann = 4
    end enum

    type Planform
        ! Wing Parameters
        integer :: WingType = Tapered ! Wing type
        integer :: WashoutDistribution = Linear ! Washout distribution type
        integer :: NNodes = 99 ! Total number of nodes
        real*8 :: AspectRatio = 5.56d0 ! Aspect ratio
        real*8 :: TaperRatio = 1.0d0 ! Taper ratio (tapered wing only)
        real*8 :: TransitionPoint = 0.25d0 ! Transition point (Combination wing only)
        real*8 :: TransitionChord = 1.0d0 ! c/croot at transtion point (Combination wing only)
        real*8 :: SectionLiftSlope = 2.0d0 * pi ! Section lift slope
        real*8 :: AileronRoot = 0.253d0 ! Location of aileron root (z/b)
        real*8 :: AileronTip = 0.438d0 ! Location of aileron tip (z/b)
        logical :: ParallelHingeLine = .true. ! Is the hinge line parallel to the
                                              ! quarter-chord line? When true,
                                              ! FlapFractionTip will be calculated
        real*8 :: DesiredFlapFractionRoot = 0.28d0 ! Desired flap fraction at aileron root (cf/c)
        real*8 :: FlapFractionRoot = 0.28d0 ! Flap fraction at aileron root (cf/c)
        real*8 :: FlapFractionTip = 0.25d0 ! Flap fraction at aileron tip (cf/c)
        real*8 :: HingeEfficiency = 0.85d0 ! Aileron hinge efficiency
        real*8 :: DeflectionEfficiency = 1.0d0 ! Aileron deflection efficiency
        integer :: LowAspectRatioMethod = Classical ! Low-Aspect-Ratio correction method

        ! Coefficients for Tapered wing with elliptic tip
        real*8 :: C1 = 0.0d0 ! Represents transition point
        real*8 :: C2 = 0.0d0 ! Represents slope of tapered section
        real*8 :: C3 = 0.0d0 ! Represents secondary axis of ellipse
        real*8 :: C4 = 0.0d0 ! Represents ellipse center offset
        real*8 :: C5 = 0.0d0 ! Represents croot/b

        ! Output Options
        logical :: OutputMatrices = .true.  ! Write C Matrix and Fourier coefficients to output file?
        character*80 :: FileName = "planform.out" ! Name of output file

        ! Operating Conditions
        real*8 :: DesiredAngleOfAttack = pi / 36.0d0 ! Desired root aerodynamic angle of Attack
                                                     ! (alpha - alpha_L0), in radians
                                                     ! When specified, a new LiftCoefficient is calculated
        real*8 :: AngleOfAttack = pi / 36.0d0 ! Root Aerodynamic Angle of Attack
                                              ! (alpha - alpha_L0), in radians
        real*8 :: DesiredLiftCoefficient = 0.4d0 ! Desired lift coefficient
                                                 ! When specified, a new AngleOfAttack is calculated
        real*8 :: LiftCoefficient = 0.4d0 ! Lift coefficient (user input, ignored if SpecifyAlpha == .true.)
        real*8 :: DesiredWashout = 0.0d0 ! Desired total washout, in radians
        real*8 :: OptimumWashout1 = 0.0d0 ! Optimum total washout, in radians (Eq. 1.8.37)
        real*8 :: OptimumWashout2 = 0.0d0 ! Optimum total washout, in radians (Eq. 1.8.42)
        real*8 :: Washout = 0.0d0 ! Total washout to use
        logical :: UseOptimumWashout = .true. ! Use the optimum total washout?
        real*8 :: AileronDeflection = 0.0d0 ! Aileron deflection, in radians
        real*8 :: DesiredRollingRate = 0.0d0 ! Desired dimensionless rolling rate (constant over wingspan)
        real*8 :: RollingRate = 0.0d0 ! Dimensionless rolling rate (constant over wingspan)
        logical :: SpecifyAlpha = .true. ! Was alpha specified?
                                         ! .true.  = Use desired alpha to calculate CL
                                         ! .false. = Use desired CL to calculate alpha
        logical :: UseSteadyRollingRate = .true. ! Use the steady dimensionless rolling rate?

        ! Planform Calculations
        real*8, allocatable, dimension(:,:) :: BigC, BigC_Inv
        real*8, allocatable, dimension(:) :: a, b, c, d, BigA
        real*8, allocatable, dimension(:) :: Omega
        logical :: IsAllocated = .false.

        ! Lift Coefficient Calculations
        real*8 :: KL  ! Lift slope factor
        real*8 :: EW  ! Washout effectiveness (epsilon omega)
        real*8 :: CLa ! Wing lift slope (derivative of CL with respect to alpha)
        real*8 :: CL1 ! Lift Coefficient (Eq. 1.8.24)
        real*8 :: CL2 ! Lift Coefficient (Eq. 1.8.5)

        ! Drag Coefficient Calculations
        real*8 :: KD   ! Induced drag factor
        real*8 :: KDL  ! Lift-washout contribution to induced drag
        real*8 :: KDW  ! Washout contribution to induced drag
        real*8 :: ES   ! Span efficiency factor
        real*8 :: CDi1 ! Induced drag coefficient (Eq. 1.8.25)
        real*8 :: CDi2 ! Induced drag coefficient (Eq. 1.8.6)
        real*8 :: CDi3 ! Induced drag coefficient (Eq. 32, Wing Flapping paper)

        ! Roll/yaw calculations
        real*8 :: CRM_da   ! Change in rolling moment coefficient with respect to alpha
        real*8 :: CRM_pbar ! Change in rolling moment coefficient with respect to rolling rate
        real*8 :: CRM      ! Rolling moment coefficient
        real*8 :: CYM      ! Yawing moment coefficient

    end type Planform

    contains
        character*80 function GetWingType(pf) result(name)
            type(Planform), intent(in) :: pf

            if (pf%WingType .eq. Tapered) then
                name = "Tapered"
            else if (pf%WingType .eq. Elliptic) then
                name = "Elliptic"
            else if (pf%WingType .eq. Combination) then
                name = "Tapered with elliptic tip"
            else
                name = "Unknown"
            end if
        end function GetWingType

        character*80 function GetWashoutDistributionType(pf) result(name)
            type(Planform), intent(in) :: pf

            if (pf%WashoutDistribution .eq. Linear) then
                name = "Linear"
            else if (pf%WashoutDistribution .eq. Optimum) then
                name = "Optimum"
            else
                name = "Unknown"
            end if
        end function GetWashoutDistributionType

        character*80 function GetLowAspectRatioMethod(pf) result(name)
            type(Planform), intent(in) :: pf

            if (pf%LowAspectRatioMethod .eq. Classical) then
                name = "Classical"
            else if (pf%LowAspectRatioMethod .eq. Hodson) then
                name = "Hodson"
            else if (pf%LowAspectRatioMethod .eq. ModifiedSlender) then
                name = "Modified Slender Wing"
            else if (pf%LowAspectRatioMethod .eq. Kuchemann) then
                name = "Kuchemann"
            else
                name = "Unknown"
            end if
        end function GetLowAspectRatioMethod

        real*8 function theta_i(i, nnodes) result(theta)
            integer, intent(in) :: i
            integer, intent(in) :: nnodes

            if (i < 1 .or. i > nnodes) then
                write(6, '(a, i3)') "ERROR: Function theta_i called with i = ", i
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

            real*8 :: zb, u

            if (pf%WingType == Tapered) then
                ! Calculate c/b for tapered wing
                cb = (2.0d0 * (1.0d0 - (1.0d0 - pf%TaperRatio) * &
                    & dabs(cos(theta)))) / (pf%AspectRatio * (1.0d0 + pf%TaperRatio))
            else if (pf%WingType == Elliptic) then
                ! Calculate c/b for elliptic wing
                cb = (4.0d0 * sin(theta)) / &
                    & (pi * pf%AspectRatio)
            else if (pf%WingType == Combination) then
                ! Calculate c/b for combination wing
                zb = abs(z_over_b(theta))
                if (zb <= pf%TransitionPoint) then
                    cb = pf%C5 * (1.0d0 - pf%C2 * zb)
                else
                    u = (zb - pf%C4) / (0.5d0 - pf%C4)
                    cb = pf%C5 * pf%C3 * sqrt(1.0d0 - u**2)
                end if
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

        real*8 function cf_over_c_i(pf, i) result(cfc)
            type(Planform), intent(in) :: pf
            integer, intent(in) :: i

            real*8 :: zbi

            zbi = z_over_b_i(i, pf%NNodes)
            if (Compare(dabs(zbi), pf%AileronRoot, zero) == -1 .or. &
                & Compare(dabs(zbi), pf%AileronTip, zero) == 1) then
                cfc = 0.0d0
            else
                cfc = 0.75d0 - y_i(pf, i) / c_over_b_i(pf, i)
            end if
        end function cf_over_c_i

        real*8 function y_i(pf, i) result(y)
            type(Planform), intent(in) :: pf
            integer, intent(in) :: i

            real*8 :: zb_i, cb_i
            real*8 :: zb_root, cfc_root, theta_root, cb_root, y_root
            real*8 :: zb_tip, cfc_tip, theta_tip, cb_tip, y_tip
            real*8 :: slope, offst

            zb_root = pf%AileronRoot
            cfc_root = pf%FlapFractionRoot
            theta_root = theta_zb(zb_root)
            cb_root = c_over_b(pf, theta_root)
            y_root = (0.75d0 - cfc_root) * cb_root

            zb_tip = pf%AileronTip
            cfc_tip = pf%FlapFractionTip
            theta_tip = theta_zb(zb_tip)
            cb_tip = c_over_b(pf, theta_tip)
            y_tip = (0.75d0 - cfc_tip) * cb_tip

            slope = (y_tip - y_root) / (zb_tip - zb_root)
            offst = y_root - slope * zb_root

            zb_i = z_over_b_i(i, pf%NNodes)
            y = slope * dabs(zb_i) + offst
        end function y_i

        real*8 function FlapEffectiveness(pf, i) result(eps_f)
            type(Planform), intent(in) :: pf
            integer, intent(in) :: i

            real*8 :: theta_f, eps_fi

            theta_f = acos(2.0d0 * cf_over_c_i(pf, i) - 1.0d0)
            eps_fi = 1.0d0 - (theta_f - sin(theta_f)) / pi
            eps_f = eps_fi * pf%HingeEfficiency * pf%DeflectionEfficiency
        end function FlapEffectiveness

        subroutine DeallocateArrays(pf)
            type(Planform), intent(inout) :: pf

            if (pf%IsAllocated) then
                deallocate(pf%BigC)
                deallocate(pf%BigC_Inv)
                deallocate(pf%a)
                deallocate(pf%b)
                deallocate(pf%c)
                deallocate(pf%d)
                deallocate(pf%BigA)
                deallocate(pf%Omega)

                pf%IsAllocated = .false.
            end if
        end subroutine DeallocateArrays

        subroutine AllocateArrays(pf)
            type(Planform), intent(inout) :: pf

            if (pf%IsAllocated) call DeallocateArrays(pf)

            allocate(pf%BigC(pf%NNodes, pf%NNodes))
            allocate(pf%BigC_Inv(pf%NNodes, pf%NNodes))
            allocate(pf%a(pf%NNodes))
            allocate(pf%b(pf%NNodes))
            allocate(pf%c(pf%NNodes))
            allocate(pf%d(pf%NNodes))
            allocate(pf%BigA(pf%NNodes))
            allocate(pf%Omega(pf%NNodes))
            pf%IsAllocated = .true.
        end subroutine AllocateArrays

end module class_Planform


