module LiftingLineSolver
    use class_Planform
    use LiftingLineSetters
    use LiftingLineOutput
    use matrix
    implicit none

contains
    subroutine ComputeCMatrixAndCoefficients(pf)
        type(Planform), intent(inout) :: pf

        write(6, '(a)') "Calculating C matrix and Fourier coefficients, please wait..."
        write(6, '(a, a, a)') "Estimated calculation time: ", &
            & trim(FormatReal(pf%NNodes**2 * 1.0d-5, 3)), " seconds"
        write(6, *)

        if (.not. pf%IsAllocated) then
            if (pf%WingType == Combination) then
                call SetCombinationWingCoefficients(pf)
            end if

            call AllocateArrays(pf)

            call ComputeC(pf, pf%BigC)
            call ComputeCInverse(pf, pf%BigC_Inv)
            call ComputeFourierCoefficients_a(pf, pf%a)
            call ComputeFourierCoefficients_b(pf, pf%b, pf%Omega)
            call ComputeFourierCoefficients_c(pf, pf%c)
            call ComputeFourierCoefficients_d(pf, pf%d)

            call ComputeLiftCoefficientParameters(pf)
            call ComputeDragCoefficientParameters(pf)
            call ComputeRollCoefficientParameters(pf)
            call ComputeFlightConditions(pf)
        end if
    end subroutine ComputeCMatrixAndCoefficients

    subroutine ComputeFourierCoefficients_a(pf, a)
        type(Planform), intent(in) :: pf
        real*8, intent(out) :: a(pf%NNodes)

        real*8 :: ones(pf%NNodes)
        integer :: i
        integer :: nnodes

        nnodes = pf%NNodes
        ones = (/ (1.0d0, i=1, nnodes) /)
        if (pf%WingType == Tapered .and. Compare(pf%TaperRatio, 0.0d0, zero) == 0) then
            ones(1) = 0.0d0
            ones(nnodes) = 0.0d0
        end if

        a = matmul(pf%BigC_Inv, ones)
    end subroutine ComputeFourierCoefficients_a

    subroutine ComputeFourierCoefficients_b(pf, b, omega)
        type(Planform), intent(in) :: pf
        real*8, intent(out) :: b(pf%NNodes)
        real*8, intent(out) :: omega(pf%NNodes)

        real*8 :: croot_over_b, theta
        integer :: i
        integer :: nnodes

        nnodes = pf%NNodes
        if (pf%WashoutDistribution == Linear) then
            omega = (/ (dabs(cos(theta_i(i, nnodes))), i=1, nnodes) /)
            if (pf%WingType == Tapered .and. Compare(pf%TaperRatio, 0.0d0, zero) == 0) then
                omega(1) = 0.0d0
                omega(nnodes) = 0.0d0
            end if
        else if (pf%WashoutDistribution == Optimum) then
            croot_over_b = c_over_b(pf, pi / 2.0d0)
            do i = 1, nnodes
                theta = theta_i(i, nnodes)
                omega(i) = 1.0d0 - sin(theta) / (c_over_b(pf, theta) / croot_over_b)
            end do

            if (pf%WingType == Combination) then
                omega(1) = 1.0d0 - sqrt(1.0d0 - 2.0d0 * pf%C4) / pf%C3
                omega(nnodes) = omega(1)
            else if (pf%WingType == Tapered .and. Compare(pf%TaperRatio, 0.0d0, zero) == 0) then
                omega(1) = 2.0d0
                omega(nnodes) = 2.0d0
            end if
        else
            write(6, '(a)') "Unknown washout distribution type!"
            stop
        end if

        b = matmul(pf%BigC_Inv, omega)
    end subroutine ComputeFourierCoefficients_b

    subroutine ComputeFourierCoefficients_c(pf, c)
        type(Planform), intent(in) :: pf
        real*8, intent(out) :: c(pf%NNodes)

        real*8 :: chi(pf%NNodes)
        real*8 :: zbi
        integer :: i
        integer :: nnodes

        nnodes = pf%NNodes
        do i = 1, nnodes
            zbi = z_over_b_i(i, nnodes)
            if (Compare(dabs(zbi), pf%AileronRoot, zero) /= -1 .and. &
                & Compare(dabs(zbi), pf%AileronTip, zero) /= 1) then
                chi(i) = -sign(FlapEffectiveness(pf, i), zbi)
            else
                chi(i) = 0.0d0
            end if
        end do

        if (pf%WingType == Tapered .and. Compare(pf%TaperRatio, 0.0d0, zero) == 0) then
            chi(1) = 0.0d0
            chi(nnodes) = 0.0d0
        end if

        c = matmul(pf%BigC_Inv, chi)
    end subroutine ComputeFourierCoefficients_c

    subroutine ComputeFourierCoefficients_d(pf, d)
        type(Planform), intent(in) :: pf
        real*8, intent(out) :: d(pf%NNodes)

        real*8 :: cos_theta(pf%NNodes)
        integer :: i
        integer :: nnodes

        nnodes = pf%NNodes
        cos_theta = (/ (cos(theta_i(i, nnodes)), i=1, nnodes) /)
        if (pf%WingType == Tapered .and. Compare(pf%TaperRatio, 0.0d0, zero) == 0) then
            cos_theta(1) = 0.0d0
            cos_theta(nnodes) = 0.0d0
        end if

        d = matmul(pf%BigC_Inv, cos_theta)
    end subroutine ComputeFourierCoefficients_d

    subroutine ComputeBigACoefficients(pf, bigA)
        type(Planform), intent(in) :: pf
        real*8, intent(out) :: bigA(pf%NNodes)

        integer :: i

        do i = 1, pf%NNodes
            bigA(i) = pf%a(i) * pf%AngleOfAttack - pf%b(i) * pf%Washout + &
                & pf%c(i) * pf%AileronDeflection + pf%d(i) * pf%RollingRate
        end do
    end subroutine ComputeBigACoefficients

    subroutine ComputeLiftCoefficientParameters(pf)
        type(Planform), intent(inout) :: pf

        pf%KL = Kappa_L(pf%AspectRatio, pf%SectionLiftSlope, pf%a(1))
        pf%EW = Epsilon_Omega(pf%a(1), pf%b(1))
        pf%CLa = C_L_alpha(pf%AspectRatio, pf%a(1))
    end subroutine ComputeLiftCoefficientParameters

    subroutine ComputeDragCoefficientParameters(pf)
        type(Planform), intent(inout) :: pf

        pf%KD = Kappa_D(pf%NNodes, pf%a)
        pf%ES = SpanEfficiencyFactor(pf%KD)
        pf%KDL = Kappa_DL(pf%NNodes, pf%a, pf%b)
        pf%KDW = Kappa_DOmega(pf%NNodes, pf%a, pf%b)
    end subroutine ComputeDragCoefficientParameters

    subroutine ComputeRollCoefficientParameters(pf)
        type(Planform), intent(inout) :: pf

        pf%CRM_da = CRM_dAlpha(pf%AspectRatio, pf%c(2))
        pf%CRM_pbar = CRM_PBar(pf%AspectRatio, pf%d(2))
    end subroutine ComputeRollCoefficientParameters

    real*8 function Kappa_L(ra, cla_section, a1) result(kl)
        real*8, intent(in) :: ra
        real*8, intent(in) :: cla_section
        real*8, intent(in) :: a1

        kl = 1.0d0 / ((1.0d0 + pi * ra / cla_section) * a1) - 1.0d0
    end function Kappa_L

    real*8 function Epsilon_Omega(a1, b1) result(ew)
        real*8, intent(in) :: a1
        real*8, intent(in) :: b1

        ew = b1 / a1
    end function Epsilon_Omega

    real*8 function C_L_alpha(ra, a1) result(cla)
        real*8, intent(in) :: ra
        real*8, intent(in) :: a1

        cla = pi * ra * a1
    end function C_L_alpha

    real*8 function Kappa_D(nnodes, a) result(kd)
        integer, intent(in) :: nnodes
        real*8, intent(in) :: a(nnodes)

        integer :: i

        kd = 0.0d0
        do i = 2, nnodes
            kd = kd + real(i, 8) * (a(i) / a(1))**2
        end do
    end function Kappa_D

    real*8 function SpanEfficiencyFactor(kd) result(es)
        real*8, intent(in) :: kd
        es = 1.0d0 / (1.0d0 + kd)
    end function SpanEfficiencyFactor

    real*8 function Kappa_DL(nnodes, a, b) result (kdl)
        integer, intent(in) :: nnodes
        real*8, intent(in) :: a(nnodes)
        real*8, intent(in) :: b(nnodes)

        integer :: i

        kdl = 0.0d0
        do i = 2, nnodes
            kdl = kdl + real(i, 8) * a(i) / a(1) * &
                & (b(i) / b(1) - a(i) / a(1))
        end do
        kdl = kdl * 2.0d0 * b(1) / a(1)
    end function Kappa_DL

    real*8 function Kappa_DOmega(nnodes, a, b) result(kdw)
        integer, intent(in) :: nnodes
        real*8, intent(in) :: a(nnodes)
        real*8, intent(in) :: b(nnodes)

        integer :: i

        kdw = 0.0d0
        do i = 2, nnodes
            kdw = kdw + real(i, 8) * (b(i) / b(1) - a(i) / a(1))**2
        end do
        kdw = kdw * (b(1) / a(1))**2
    end function Kappa_DOmega

    real*8 function CRM_dAlpha(ra, c2) result(crmda)
        real*8, intent(in) :: ra
        real*8, intent(in) :: c2

        crmda = -pi * ra / 4.0d0 * c2
    end function CRM_dAlpha

    real*8 function CRM_PBar(ra, d2) result(crmpbar)
        real*8, intent(in) :: ra
        real*8, intent(in) :: d2

        crmpbar = -pi * ra / 4.0d0 * d2
    end function CRM_PBar

    subroutine ComputeFlightConditions(pf)
        type(Planform), intent(inout) :: pf

        ! Make sure planform characteristics have been computed
        if (.not. pf%IsAllocated) then
            call ComputeCMatrixAndCoefficients(pf)
        end if

        ! Compute root aerodynamic angle of attack, if necessary
        if (.not. pf%SpecifyAlpha) then
            pf%AngleOfAttack = RootAlpha(pf%CLa, pf%LiftCoefficient, pf%EW, pf%Washout)
        else
            pf%LiftCoefficient = CL1(pf%CLa, pf%AngleOfAttack, pf%EW, pf%Washout)
        end if

        ! Compute optimum total washout, if necessary
        if (pf%UseOptimumWashout) then
            call SetOptimumWashout(pf)
        else
            call SetWashout(pf, pf%DesiredWashout * 180.0d0 / pi)
        end if

        ! Compute steady rolling rate, if necessary
        if (pf%UseSteadyRollingRate) then
            call SetSteadyRollingRate(pf)
        end if

        ! Compute BigA Fourier Coefficients
        call ComputeBigACoefficients(pf, pf%BigA)

        ! Compute lift coefficients
        call ComputeLiftCoefficients(pf)

        ! Compute drag coefficient
        call ComputeDragCoefficients(pf)

        ! Compute roll coefficient
        pf%CRM = CRoll(pf%CRM_da, pf%CRM_pbar, pf%AileronDeflection, pf%RollingRate)

        ! Compute yaw coefficient
        pf%CYM = CYaw(pf, pf%CL1, pf%BigA)
    end subroutine ComputeFlightConditions

    subroutine ComputeLiftCoefficients(pf)
        type(Planform), intent(inout) :: pf

        pf%CL1 = CL1(pf%CLa, pf%AngleOfAttack, pf%EW, pf%Washout)
        pf%CL2 = CL2(pf%AspectRatio, pf%BigA(1))
    end subroutine ComputeLiftCoefficients

    subroutine ComputeDragCoefficients(pf)
        type(Planform), intent(inout) :: pf

        pf%CDi1 = CDi1(pf)
        pf%CDi2 = CDi2(pf)
        pf%CDi3 = CDi3(pf)
    end subroutine ComputeDragCoefficients

    real*8 function CRoll(crmda, crmpbar, da, pbar) result(crm)
        real*8, intent(in) :: crmda
        real*8, intent(in) :: crmpbar
        real*8, intent(in) :: da
        real*8, intent(in) :: pbar

        crm = crmda * da + crmpbar * pbar
    end function CRoll

    real*8 function CYaw(pf, cl, bigA) result(cym)
        type(Planform), intent(in) :: pf
        real*8, intent(in) :: cl
        real*8, intent(in) :: bigA(pf%NNodes)

        integer :: i
        integer :: nnodes
        nnodes = pf%NNodes

        cym = cl / 8.0d0 * (6.0d0 * bigA(2) - pf%RollingRate) + &
            & pi * pf%AspectRatio / 8.0d0 * (10.0d0 * bigA(2) - &
            & pf%RollingRate) * bigA(3)
        do i = 4, nnodes
            cym = cym + 0.25d0 * pi * pf%AspectRatio * &
                & (2.0d0 * real(i, 8) - 1.0d0) * bigA(i-1) * bigA(i)
        end do
    end function CYaw

    subroutine ComputeC(pf, c)
        type(Planform), intent(in) :: pf
        real*8, intent(inout) :: c(pf%NNodes, pf%NNodes)

        integer :: i
        integer :: nnodes

        nnodes = pf%NNodes

        ! Compute values for i=1, i=N
        call C1j_Nj(c, pf)

        ! Compute values for i=2 to i=N-1
        do i = 2, nnodes-1
            call Cij(c, i, pf)
        end do
    end subroutine ComputeC

    subroutine ComputeCInverse(pf, c_inv)
        type(Planform), intent(in) :: pf
        real*8, intent(inout) :: c_inv(pf%NNodes, pf%NNodes)

        call matinv_gauss(pf%NNodes, pf%BigC, c_inv)
    end subroutine ComputeCInverse

    subroutine C1j_Nj(c, pf)
        real*8, dimension(:,:), intent(inout) :: c
        type(Planform), intent(in) :: pf

        integer :: j
        integer :: jsq
        integer :: nnode
        real*8 :: cb0

        nnode = pf%NNodes
        do j = 1, nnode
            jsq = j**2
            c(1, j) = real(jsq, 8)
            c(nnode, j) = real((-1)**(j + 1) * jsq, 8)
        end do

        cb0 = c_over_b(pf, pi)
        if (dabs(cb0) < 1.0d-10) then
            call C1j_Nj_zero_chord(c, pf)
        end if

    end subroutine C1j_Nj

    subroutine Cij(c, i, pf)
        real*8, dimension(:,:), intent(inout) :: c
        integer, intent(in) :: i
        type(Planform), intent(in) :: pf

        integer :: j
        integer :: nnode
        real*8 :: theta
        real*8 :: cb
        real*8 :: sin_theta

        nnode = pf%NNodes
        theta = theta_i(i, nnode)
        cb = c_over_b_i(pf, i)
        sin_theta = sin(theta)

        do j = 1, nnode
            c(i, j) = (4.0d0 / (pf%SectionLiftSlope * cb) + &
                & real(j, 8) / sin_theta) * sin(real(j, 8) * theta)
        end do
    end subroutine Cij

    subroutine C1j_Nj_zero_chord(c, pf)
        real*8, dimension(:,:), intent(inout) :: c
        type(Planform), intent(in) :: pf

        integer :: j, n

        n = pf%NNodes

        if (pf%WingType == Tapered) then
            do j = 1, n
                c(1, j) = 2.0d0 * pf%AspectRatio * (1.0d0 + real(j, 8))
                c(n, j) = real((-1)**(j + 1), 8) * c(1, j)
            end do
        else if (pf%WingType == Elliptic) then
            do j = 1, n
                c(1, j) = c(1, j) + real(j, 8) * pi * &
                    & pf%AspectRatio / pf%SectionLiftSlope
                c(n, j) = c(n, j) + real((-1)**(j + 1) * j, 8) * pi * &
                    & pf%AspectRatio / pf%SectionLiftSlope
            end do
        else if (pf%WingType == Combination) then
            ! TODO: Add code for combination wing type
            do j = 1, n
                c(1, j) = c(1, j) + 4.0d0 * real(j, 8) * &
                    & sqrt(1.0d0 - 2.0d0 * pf%C4) / &
                    & (pf%C3 * pf%C5 * pf%SectionLiftSlope)
                c(n, j) = c(n, j) + 4.0d0 * real((-1)**(j + 1) * j, 8) * &
                    & sqrt(1.0d0 - 2.0d0 * pf%C4) / &
                    & (pf%C3 * pf%C5 * pf%SectionLiftSlope)
            end do
        else
            stop "*** Unknown Wing Type ***"
        end if

    end subroutine C1j_Nj_zero_chord

end module LiftingLineSolver
