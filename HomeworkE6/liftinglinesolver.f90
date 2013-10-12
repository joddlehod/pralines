module liftinglinesolver
    use class_Planform
    use matrix
    implicit none

contains
    subroutine runsimulation(pf)
        type(Planform), intent(in) :: pf

        integer :: i
        real*8 :: c(pf%NNodes, pf%NNodes)
        real*8 :: c_inv(pf%NNodes, pf%NNodes)
        real*8 :: a(pf%NNodes)

        write(6, '(a2, 2x, a22, 2x, a22, 2x, a22)') "i", "Theta", "z/b", "c/b"
        do i = 1,pf%NNodes
            write(6, '(i2, 2x, ES22.15, 2x, ES22.15, 2x, ES22.15)') &
                & i, theta(i, pf%NNodes), z_over_b(i, pf%NNodes), c_over_b(i, pf)
        end do

        call ComputeC(pf, c)
        call ComputeCInverse(pf, c, c_inv)

        call ComputeFourierCoefficients_a(pf%NNodes, c_inv, a)
        write(6, *)
        write(6, '(a)') "Fourier Coefficients:"
        do i = 1, pf%NNodes
            write(6, '(a, i2, a, f22.15)') "a", i, " = ", a(i)
        end do

    end subroutine runsimulation

    subroutine ComputeFourierCoefficients_a(nnodes, c_inv, a)
        integer, intent(in) :: nnodes
        real*8, intent(in) :: c_inv(nnodes, nnodes)
        real*8, intent(out) :: a(nnodes)

        integer :: i

        real*8 :: ones(nnodes)
        ones = (/ (1.0d0, i=1, nnodes) /)

        a = matmul(c_inv, ones)
    end subroutine ComputeFourierCoefficients_a

    real*8 function Kappa_L(pf, a1) result(kl)
        type(Planform), intent(in) :: pf
        real*8, intent(in) :: a1

        kl = 1.0d0 / ((1.0d0 + pi * pf%AspectRatio / pf%LiftSlope) * a1) - 1.0d0
    end function Kappa_L

    real*8 function Kappa_D(pf, nnodes, a) result(kd)
        type(Planform), intent(in) :: pf
        integer, intent(in) :: nnodes
        real*8, intent(in) :: a(nnodes)

        integer :: i

        kd = 0.0d0
        do i = 2, nnodes
            kd = kd + nnodes * (a(i) / a(1))**2
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
            c_over_b_i = (2.0d0 * (pf%TaperRatio - 1.0d0) * &
                & dabs(z_over_b(i, pf%NNodes)) + 1.0d0) / pf%AspectRatio
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

        z_over_b_i = (real(i, 8) - 1.0d0) / (real(nnodes, 8) - 1.0d0) - 0.5d0
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

        write(6, *)
        write(6, *) "[C] Matrix:"
        call printmat(nnodes, nnodes, c)
    end subroutine ComputeC

    subroutine ComputeCInverse(pf, c, c_inv)
        type(Planform), intent(in) :: pf
        real*8, intent(inout) :: c(pf%NNodes, pf%NNodes)
        real*8, intent(inout) :: c_inv(pf%NNodes, pf%NNodes)

        call matinv(pf%NNodes, c, c_inv)

        write(6, *)
        write(6, *) "[C]^-1 Matrix:"
        call printmat(pf%NNodes, pf%NNodes, c_inv)
    end subroutine ComputeCInverse

    subroutine C1j_Nj(c, pf)
        real*8, dimension(:,:), intent(inout) :: c
        type(Planform), intent(in) :: pf

        integer :: j
        integer :: jsq
        integer :: nnode
        real*8 :: c_over_b_1
        real*8 :: c_over_b_N

        nnode = pf%NNodes
        do j = 1, nnode
            jsq = j**2
            c(1, j) = real(j**2, 8)
            c(nnode, j) = real((-1)**(j + 1) * jsq, 8)
        end do

        c_over_b_1 = c_over_b(1, pf)
        if (dabs(c_over_b_1) < 1.0d-10) then
            call C1j_Nj_zero_chord(c, 1, pf)
        end if

        c_over_b_N = c_over_b(nnode, pf)
        if (dabs(c_over_b_N) < 1.0d-10) then
            call C1j_Nj_zero_chord(c, nnode, pf)
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
            c(i, j) = (4.0d0 / pf%LiftSlope / c_over_b_i + &
                & real(j, 8) / sin_theta_i) * sin(real(j, 8) * theta_i)
        end do
    end subroutine Cij

    subroutine C1j_Nj_zero_chord(c, i, pf)
        real*8, dimension(:,:), intent(inout) :: c
        integer, intent(in) :: i
        type(Planform), intent(in) :: pf

        integer :: j

        if (pf%WingType == Tapered) then
            ! limit = 0, so do nothing
        else if (pf%WingType == Elliptic) then
            do j = 1, pf%NNodes
                c(i, j) = c(i, j) + pi * pf%AspectRatio / pf%LiftSlope
            end do
        else
            stop "*** Unknown Wing Type ***"
        end if

    end subroutine C1j_Nj_zero_chord


end module liftinglinesolver
