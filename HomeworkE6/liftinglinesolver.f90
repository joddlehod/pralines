module liftinglinesolver
    use class_Planform
    implicit none

contains
    subroutine runsimulation(pf)
        type(Planform), intent(in) :: pf

        integer :: i, j
        character*80 :: format_string
        real*8 :: C(pf%NNodes, pf%NNodes)

        write(6, '(a2, 2x, a22, 2x, a22, 2x, a22)') "i", "Theta", "z/b", "c/b"
        do i = 1,pf%NNodes
            write(6, '(i2, 2x, ES22.15, 2x, ES22.15, 2x, ES22.15)') &
                & i, theta(i, pf%NNodes), z_over_b(i, pf%NNodes), c_over_b(i, pf)
        end do

        call ComputeC(C, pf)

        write(format_string, '(a, i3, a)') "(", pf%NNodes, "(F9.5, 2x))"
        write(6, *)
        write(6, '(a)') "[C] Matrix:"
        do i = 1,pf%NNodes
            write(6, format_string) (C(i, j), j=1, pf%NNodes)
        end do

    end subroutine runsimulation

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
                & abs(z_over_b(i, pf%NNodes)) + 1.0d0) / pf%AspectRatio
        else if (pf%WingType == Elliptic) then
            ! Calculate c/b for elliptic wing
            c_over_b_i = (4.0d0 * sin(theta(i, pf%NNodes))) / &
                & (pi * pf%AspectRatio)
        else
            ! Unknown wing type!
            stop "*** Unknown Wing Type ***"
        end if

    end function c_over_b

    real*8 function z_over_b(i, nnode) result(z_over_b_i)
        integer, intent(in) :: i
        integer, intent(in) :: nnode

        z_over_b_i = (real(i, 8) - 1.0d0) / (real(nnode, 8) - 1.0d0) - 0.5d0
    end function z_over_b

    subroutine ComputeC(C, pf)
        real*8, dimension(:,:), intent(inout) :: C
        type(Planform), intent(in) :: pf

        integer :: i
        integer :: nnode

        nnode = pf%NNodes

        ! Compute values for i=1, i=N
        call C1j_Nj(C, pf)

        ! Compute values for i=2 to i=N-1
        do i = 2, nnode-1
            call Cij(C, i, pf)
        end do

    end subroutine ComputeC

    subroutine C1j_Nj(C, pf)
        real*8, dimension(:,:), intent(inout) :: C
        type(Planform), intent(in) :: pf

        integer :: j
        integer :: jsq
        integer :: nnode
        real*8 :: c_over_b_1
        real*8 :: c_over_b_N

        nnode = pf%NNodes
        do j = 1, nnode
            jsq = j**2
            C(1, j) = real(j**2, 8)
            C(nnode, j) = real((-1)**(j + 1) * jsq, 8)
        end do

        c_over_b_1 = c_over_b(1, pf)
        if (abs(c_over_b_1) < 1.0d-10) then
            call C1j_Nj_zero_chord(C, 1, pf)
        end if

        c_over_b_N = c_over_b(nnode, pf)
        if (abs(c_over_b_N) < 1.0d-10) then
            call C1j_Nj_zero_chord(C, nnode, pf)
        end if

    end subroutine C1j_Nj

    subroutine Cij(C, i, pf)
        real*8, dimension(:,:), intent(inout) :: C
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
            C(i, j) = (4.0d0 / pf%LiftSlope / c_over_b_i + &
                & real(j, 8) / sin_theta_i) * sin(real(j, 8) * theta_i)
        end do
    end subroutine Cij

    subroutine C1j_Nj_zero_chord(C, i, pf)
        real*8, dimension(:,:), intent(inout) :: C
        integer, intent(in) :: i
        type(Planform), intent(in) :: pf

        integer :: j

        if (pf%WingType == Tapered) then
            ! limit = 0, so do nothing
        else if (pf%WingType == Elliptic) then
            do j = 1, pf%NNodes
                C(i, j) = C(i, j) + pi * pf%AspectRatio / pf%LiftSlope
            end do
        else
            stop "*** Unknown Wing Type ***"
        end if

    end subroutine C1j_Nj_zero_chord


end module liftinglinesolver
