module liftinglinesolver
    use class_Planform
    use matrix
    implicit none

contains
    subroutine RunSimulation(pf)
        type(Planform), intent(in) :: pf

        real*8 :: c(pf%NNodes, pf%NNodes)
        real*8 :: c_inv(pf%NNodes, pf%NNodes)
        real*8 :: a(pf%NNodes)
        logical :: print_to_file

        print_to_file = (pf%WriteCMatrix .or. pf%WriteCInverse &
            & .or. pf%WriteFourier .or. pf%WriteOther)
        if (print_to_file) then
            open(unit=10, file=pf%FileName)
        end if

        call PrintPlanformSummary(6, pf)
        if (print_to_file) then
            call PrintPlanformSummary(10, pf)
        end if

        call ComputeC(pf, c)
        call ComputeCInverse(pf, c, c_inv)
        call ComputeFourierCoefficients_a(pf, c_inv, a)
        call ComputeWingCoefficients(pf, a)

        if (print_to_file) then
            close(unit=10)
        end if

    end subroutine RunSimulation

    subroutine PrintPlanformSummary(u, pf)
        integer, intent(in) :: u
        type(Planform), intent(in) :: pf

        character*8 :: wtype
        character*16 :: fmt_str
        integer :: log_nnodes

        if (pf%WingType == Tapered) then
            wtype = "Tapered"
        else if (pf%WingType == Elliptic) then
            wtype = "Elliptic"
        else
            wtype = "Unknown"
        end if

        log_nnodes = log10(real(pf%NNodes)) + 1
        write(fmt_str, '(a,i1,a,i1,a)') "(2x,a,i", log_nnodes, ",a,i", &
            & log_nnodes, ",a)"

        write(u, *) "Planform Summary:"
        write(u, '(2x,a,a)') "Wing type = ", wtype
        write(u, '(2x,a,ES22.15)') "Airfoil section lift slope = ", pf%LiftSlope
        write(u, fmt_str) "N  = ", pf%NNodes, " (", &
            & (pf%NNodes + 1) / 2, " nodes per semispan)"
        write(u, '(2x,a,ES22.15)') "RA = ", pf%AspectRatio
        if (pf%WingType == Tapered) then
            write(u, '(2x,a,ES22.15)') "RT = ", pf%TaperRatio
        end if
        write(u, '(2x,a,ES22.15,a)') "alpha = ", pf%AngleOfAttack, " rad"
        write(u, *)
    end subroutine PrintPlanformSummary

    subroutine ComputeFourierCoefficients_a(pf, c_inv, a)
        type(Planform), intent(in) :: pf
        real*8, intent(in) :: c_inv(pf%NNodes, pf%NNodes)
        real*8, intent(out) :: a(pf%NNodes)

        real*8 :: ones(pf%NNodes)
        integer :: i
        integer :: nnodes

        nnodes = pf%NNodes
        ones = (/ (1.0d0, i=1, nnodes) /)
        a = matmul(c_inv, ones)

        call PrintFourierCoefficients_a(6, nnodes, a)
        if (pf%WriteFourier) then
            call PrintFourierCoefficients_a(10, nnodes, a)
        end if
    end subroutine ComputeFourierCoefficients_a

    subroutine PrintFourierCoefficients_a(u, nnodes, a)
        integer, intent(in) :: u
        integer, intent(in) :: nnodes
        real*8, intent(in) :: a(nnodes)

        integer :: i

        write(u, '(a)') "Fourier Coefficients:"
        do i = 1, nnodes
            write(u, '(a, i2, a, ES22.15)') "a", i, " = ", a(i)
        end do
        write(u, *)
    end subroutine PrintFourierCoefficients_a

    subroutine ComputeWingCoefficients(pf, a)
        type(Planform), intent(in) :: pf
        real*8, intent(in) :: a(pf%NNodes)

        real*8 :: kl   ! Lift slope factor
        real*8 :: kd   ! Induced drag factor
        real*8 :: es   ! Span efficiency factor
        real*8 :: cla  ! Wing lift slope
        real*8 :: cl   ! Wing lift coefficient
        real*8 :: cdi  ! Wing induced drag coefficient

        kl = Kappa_L(pf, a(1))
        kd = Kappa_D(pf, a)
        es = SpanEfficiencyFactor(kd)
        cla = C_L_alpha(pf, kl)
        cl = C_L(cla, pf%AngleOfAttack, 0.0d0)
        cdi = C_Di(cl, pf%AspectRatio, es)

        call PrintWingCoefficients(6, kl, kd, es, cla, cl, cdi)
        if (pf%WriteOther) then
            call PrintWingCoefficients(10, kl, kd, es, cla, cl, cdi)
        end if
    end subroutine ComputeWingCoefficients

    subroutine PrintWingCoefficients(u, kl, kd, es, cla, cl, cdi)
        integer, intent(in) :: u   ! Output unit (6 = standard output)
        real*8, intent(in) :: kl   ! Lift slope factor
        real*8, intent(in) :: kd   ! Induced drag factor
        real*8, intent(in) :: es   ! Span efficiency factor
        real*8, intent(in) :: cla  ! Wing lift slope
        real*8, intent(in) :: cl   ! Wing lift coefficient
        real*8, intent(in) :: cdi  ! Wing induced drag coefficient

        write(u, '(a, ES22.15)') "KL  = ", kl
        write(u, '(a, ES22.15)') "CLa = ", cla
        write(u, '(a, ES22.15)') "CL  = ", cl
        write(u, *)

        write(u, '(a, ES22.15)') "KD  = ", kd
        write(u, '(a, ES22.15)') "es  = ", es
        write(u, '(a, ES22.15)') "CDi = ", cdi
        write(u, *)
    end subroutine PrintWingCoefficients

    real*8 function Kappa_L(pf, a1) result(kl)
        type(Planform), intent(in) :: pf
        real*8, intent(in) :: a1

        kl = 1.0d0 / ((1.0d0 + pi * pf%AspectRatio / pf%LiftSlope) * a1) - 1.0d0
    end function Kappa_L

    real*8 function Kappa_D(pf, a) result(kd)
        type(Planform), intent(in) :: pf
        real*8, intent(in) :: a(pf%NNodes)

        integer :: i

        kd = 0.0d0
        do i = 2, pf%NNodes
            kd = kd + real(i, 8) * (a(i) / a(1))**2
        end do
    end function Kappa_D

    real*8 function SpanEfficiencyFactor(kd) result(es)
        real*8, intent(in) :: kd
        es = 1.0d0 / (1.0d0 + kd)
    end function SpanEfficiencyFactor

    real*8 function C_L_alpha(pf, kl) result(cla)
        type(Planform), intent(in) :: pf
        real*8, intent(in) :: kl

        cla = pf%LiftSlope / (1.0d0 + pf%LiftSlope / (pi * pf%AspectRatio)) / &
            & (1.0d0 + kl)
    end function C_L_alpha

    real*8 function C_L(cla, alpha, alpha_l0) result(cl)
        real*8, intent(in) :: cla
        real*8, intent(in) :: alpha
        real*8, intent(in) :: alpha_l0

        cl = cla * (alpha - alpha_l0)
    end function C_L

    real*8 function C_Di(cl, ra, es) result (cdi)
        real*8, intent(in) :: cl
        real*8, intent(in) :: ra
        real*8, intent(in) :: es

        cdi = (cl * cl) / (pi * ra * es)
    end function C_Di

    real*8 function theta(i, nnode) result(theta_i)
        integer, intent(in) :: i
        integer, intent(in) :: nnode

        theta_i = real(i - 1, 8) * pi / real(nnode - 1, 8)
    end function theta

    real*8 function c_over_b(i, pf) result(c_over_b_i)
        integer, intent(in) :: i
        type(Planform), intent(in) :: pf

        if (pf%WingType == Tapered) then
            ! Calculate c/b for tapered wing
            c_over_b_i = (2.0d0 * dabs(z_over_b(i, pf%NNodes)) * &
                & (pf%TaperRatio - 1.0d0) + 1.0d0) / pf%AspectRatio
        else if (pf%WingType == Elliptic) then
            ! Calculate c/b for elliptic wing
            c_over_b_i = (4.0d0 * sin(theta(i, pf%NNodes))) / &
                & (pi * pf%AspectRatio)
        else
            ! Unknown wing type!
            stop "*** Unknown Wing Type ***"
        end if

    end function c_over_b

    real*8 function z_over_b(i, nnodes) result(z_over_b_i)
        integer, intent(in) :: i
        integer, intent(in) :: nnodes

        z_over_b_i = -0.5d0 * cos(theta(i, nnodes))
    end function z_over_b

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

        call PrintC(6, nnodes, c)
        if (pf%WriteCMatrix) then
            call PrintC(10, nnodes, c)
        end if
    end subroutine ComputeC

    subroutine PrintC(u, nnodes, c)
        integer, intent(in) :: u
        integer, intent(in) :: nnodes
        real*8, intent(in) :: c(nnodes, nnodes)

        write(u, *) "[C] Matrix:"
        call printmat(u, nnodes, nnodes, c)
        write(u, *)

    end subroutine PrintC

    subroutine ComputeCInverse(pf, c, c_inv)
        type(Planform), intent(in) :: pf
        real*8, intent(inout) :: c(pf%NNodes, pf%NNodes)
        real*8, intent(inout) :: c_inv(pf%NNodes, pf%NNodes)

        call matinv_gauss(pf%NNodes, c, c_inv)

        call PrintCInverse(6, pf%NNodes, c_inv)
        if (pf%WriteCInverse) then
            call PrintCInverse(10, pf%NNodes, c_inv)
        end if
    end subroutine ComputeCInverse

    subroutine PrintCInverse(u, nnodes, c_inv)
        integer, intent(in) :: u
        integer, intent(in) :: nnodes
        real*8, intent(in) :: c_inv(nnodes, nnodes)

        write(u, *) "[C]^-1 Matrix:"
        call printmat(u, nnodes, nnodes, c_inv)
        write(u, *)

    end subroutine PrintCInverse

    subroutine C1j_Nj(c, pf)
        real*8, dimension(:,:), intent(inout) :: c
        type(Planform), intent(in) :: pf

        integer :: j
        integer :: jsq
        integer :: nnode
        real*8 :: c_over_b_1

        nnode = pf%NNodes
        do j = 1, nnode
            jsq = j**2
            c(1, j) = real(j**2, 8)
            c(nnode, j) = real((-1)**(j + 1) * jsq, 8)
        end do

        c_over_b_1 = c_over_b(1, pf)
        if (dabs(c_over_b_1) < 1.0d-10) then
            call C1j_Nj_zero_chord(c, pf)
        end if

    end subroutine C1j_Nj

    subroutine Cij(c, i, pf)
        real*8, dimension(:,:), intent(inout) :: c
        integer, intent(in) :: i
        type(Planform), intent(in) :: pf

        integer :: j
        integer :: nnode
        real*8 :: theta_i
        real*8 :: c_over_b_i
        real*8 :: sin_theta_i

        nnode = pf%NNodes
        theta_i = theta(i, nnode)
        c_over_b_i = c_over_b(i, pf)
        sin_theta_i = sin(theta_i)

        do j = 1, nnode
            c(i, j) = (4.0d0 / (pf%LiftSlope * c_over_b_i) + &
                & real(j, 8) / sin_theta_i) * sin(real(j, 8) * theta_i)
        end do
    end subroutine Cij

    subroutine C1j_Nj_zero_chord(c, pf)
        real*8, dimension(:,:), intent(inout) :: c
        type(Planform), intent(in) :: pf

        integer :: j, n

        n = pf%NNodes

        if (pf%WingType == Tapered) then
            ! limit = 0, so do nothing
        else if (pf%WingType == Elliptic) then
            do j = 1, pf%NNodes
                c(1, j) = c(1, j) + real(j, 8) * pi * &
                    & pf%AspectRatio / pf%LiftSlope
                c(n, j) = c(n, j) + real((-1)**(j + 1) * j, 8) * pi * &
                    & pf%AspectRatio / pf%LiftSlope
            end do
        else
            stop "*** Unknown Wing Type ***"
        end if

    end subroutine C1j_Nj_zero_chord


end module liftinglinesolver
