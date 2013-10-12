module liftinglinesolver
    use class_Planform
    implicit none

contains
    subroutine runsimulation(pf)
        type(Planform), intent(in) :: pf

        integer :: i

        write(6, '(a2, 2x, a22, 2x, a22, 2x, a22)') "i", "Theta", "z/b", "c/b"
        do i = 1,pf%NNodes
            write(6, '(i2, 2x, ES22.15, 2x, ES22.15, 2x, ES22.15)') &
                & i, theta(i, pf%NNodes), z_over_b(i, pf%NNodes), c_over_b(i, pf)
        end do
    end subroutine runsimulation

    real*8 function theta(i, nnode) result(theta_i)
        integer, intent(in) :: i
        integer, intent(in) :: nnode

        theta_i = real(i - 1) * pi / real(nnode - 1)
    end function theta

    real*8 function c_over_b(i, pf) result(c_over_b_i)
        integer, intent(in) :: i
        type(Planform), intent(in) :: pf

        if (pf%WingType == Tapered) then
            ! Calculate c/b for tapered wing
            c_over_b_i = 2.0d0 * (pf%TaperRatio - 1.0d0) * &
                & abs(z_over_b(i, pf%NNodes)) + 1.0d0
        else if (pf%WingType == Elliptic) then
            ! Calculate c/b for elliptic wing
            c_over_b_i = (4.0d0 * sin(theta(i, pf%NNodes))) / &
                & (pi * pf%AspectRatio)
        else
            ! Unknown wing type!
            c_over_b_i = 0.0d0
        end if

    end function c_over_b

    real*8 function z_over_b(i, nnode) result(z_over_b_i)
        integer, intent(in) :: i
        integer, intent(in) :: nnode

        z_over_b_i = (real(i) - 1.0d0) / (real(nnode) - 1.0d0) - 0.5d0
    end function z_over_b

end module liftinglinesolver
