module LiftingLineSolver
    use class_Planform
    use LiftingLineOutput
    use matrix
    implicit none

contains
    subroutine ComputeFourierCoefficients_a(pf, a)
        type(Planform), intent(in) :: pf
        real*8, intent(out) :: a(pf%NNodes)

        real*8 :: ones(pf%NNodes)
        integer :: i
        integer :: nnodes

        nnodes = pf%NNodes
        ones = (/ (1.0d0, i=1, nnodes) /)
        if (pf%WingType == Tapered .and. pf%TaperRatio < 1.0d-10) then
            ones(1) = 0.0d0
            ones(nnodes) = 0.0d0
        end if

        a = matmul(pf%BigC_Inv, ones)
    end subroutine ComputeFourierCoefficients_a

    subroutine ComputeFourierCoefficients_b(pf, b)
        type(Planform), intent(in) :: pf
        real*8, intent(out) :: b(pf%NNodes)

        real*8 :: omega(pf%NNodes)
        integer :: i
        integer :: nnodes

        nnodes = pf%NNodes
        omega = (/ (dabs(cos(theta_i(i, nnodes))), i=1, nnodes) /)
        if (pf%WingType == Tapered .and. pf%TaperRatio < 1.0d-10) then
            omega(1) = 0.0d0
            omega(nnodes) = 0.0d0
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

        if (pf%WingType == Tapered .and. pf%TaperRatio < zero) then
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
        if (pf%WingType == Tapered .and. pf%TaperRatio < 1.0d-10) then
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
            bigA(i) = pf%a(i) * pf%AngleOfAttack - pf%b(i) * pf%Omega + &
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

    real*8 function C_L(cla, alpha, ew, omega) result(cl)
        real*8, intent(in) :: cla
        real*8, intent(in) :: alpha
        real*8, intent(in) :: ew
        real*8, intent(in) :: omega

        cl = cla * (alpha - ew * omega)
    end function C_L

    real*8 function RootAlpha(cla, cl, ew, omega) result(alpha)
        real*8, intent(in) :: cla
        real*8, intent(in) :: cl
        real*8, intent(in) :: ew
        real*8, intent(in) :: omega

        alpha = cl / cla + ew * omega
    end function RootAlpha

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

    real*8 function C_Di(cl, kd, kdl, cla, omega, kdw, ra) result (cdi)
        real*8, intent(in) :: cl
        real*8, intent(in) :: kd
        real*8, intent(in) :: kdl
        real*8, intent(in) :: cla
        real*8, intent(in) :: omega
        real*8, intent(in) :: kdw
        real*8, intent(in) :: ra

        cdi = (cl * cl * (1.0d0 + kd) - kdl * cl * cla * omega + &
            & kdw * (cla * omega)**2) / (pi * ra)
    end function C_Di

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
            cym = cym + pi * pf%AspectRatio / 4.0d0 * (2.0d0 * real(i, 8) - &
                & 1.0d0) * bigA(i-1) * bigA(i)
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
        real*8 :: cb1

        nnode = pf%NNodes
        do j = 1, nnode
            jsq = j**2
            c(1, j) = real(j**2, 8)
            c(nnode, j) = real((-1)**(j + 1) * jsq, 8)
        end do

        cb1 = c_over_b_i(pf, 1)
        if (dabs(cb1) < 1.0d-10) then
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
            ! TODO: Add code for RT = 0.0 here!
        else if (pf%WingType == Elliptic) then
            do j = 1, pf%NNodes
                c(1, j) = c(1, j) + real(j, 8) * pi * &
                    & pf%AspectRatio / pf%SectionLiftSlope
                c(n, j) = c(n, j) + real((-1)**(j + 1) * j, 8) * pi * &
                    & pf%AspectRatio / pf%SectionLiftSlope
            end do
        else
            stop "*** Unknown Wing Type ***"
        end if

    end subroutine C1j_Nj_zero_chord

end module LiftingLineSolver
