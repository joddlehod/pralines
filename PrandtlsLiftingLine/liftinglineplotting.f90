module LiftingLinePlotting
    use class_Planform
    implicit none

contains
    subroutine PlotPlanform(pf)
        type(Planform), intent(in) :: pf

        integer :: i

        ! Generate temporary text file for plotting
        open(unit=11, file='planform.dat')
        write(11, '(a)') "$ Planform Geometry"

        ! Write data points for planform
        write(11, '(a)') "! Wing"
        do i=1, pf%NNodes
            write(11, '(f22.15, a, 2x, f22.15)') &
                & z_over_b_i(i, pf%NNodes), ";", 0.25d0 * c_over_b_i(pf, i)
            write(11, '(f22.15, a, 2x, f22.15)') &
                & z_over_b_i(i, pf%NNodes), ";", -0.75d0 * c_over_b_i(pf, i)
            write(11, '(f22.15, a, 2x, f22.15)') &
                & z_over_b_i(i, pf%NNodes), ";", 0.25d0 * c_over_b_i(pf, i)
        end do
        do i=pf%NNodes, 1, -1
            write(11, '(f22.15, a, 2x, f22.15)') &
                & z_over_b_i(i, pf%NNodes), ";", -0.75d0 * c_over_b_i(pf, i)
        end do
        write(11, '(f22.15, a, 2x, f22.15)') &
            & z_over_b_i(1, pf%NNodes), ";", 0.25d0 * c_over_b_i(pf, 1)

        ! Write data points for right aileron
        write(11, '(a)') "$"
        write(11, '(a)') "! Right Aileron"
        write(11, '(f22.15, a, 2x, f22.15)') pf%AileronRoot, ";", &
            & -0.75d0 * c_over_b_zb(pf, pf%AileronRoot)
        write(11, '(f22.15, a, 2x, f22.15)') pf%AileronRoot, ";", &
            & (-0.75d0 + pf%FlapFractionRoot) * c_over_b_zb(pf, pf%AileronRoot)
        write(11, '(f22.15, a, 2x, f22.15)') pf%AileronTip, ";", &
            & (-0.75d0 + pf%FlapFractionTip) * c_over_b_zb(pf, pf%AileronTip)
        write(11, '(f22.15, a, 2x, f22.15)') pf%AileronTip, ";", &
            & -0.75d0 * c_over_b_zb(pf, pf%AileronTip)

        ! Write data points for left aileron
        write(11, '(a)') "$"
        write(11, '(a)') "! Left Aileron"
        write(11, '(f22.15, a, 2x, f22.15)') -pf%AileronRoot, ";", &
            & -0.75d0 * c_over_b_zb(pf, pf%AileronRoot)
        write(11, '(f22.15, a, 2x, f22.15)') -pf%AileronRoot, ";", &
            & (-0.75d0 + pf%FlapFractionRoot) * c_over_b_zb(pf, pf%AileronRoot)
        write(11, '(f22.15, a, 2x, f22.15)') -pf%AileronTip, ";", &
            & (-0.75d0 + pf%FlapFractionTip) * c_over_b_zb(pf, pf%AileronTip)
        write(11, '(f22.15, a, 2x, f22.15)') -pf%AileronTip, ";", &
            & -0.75d0 * c_over_b_zb(pf, pf%AileronTip)

        ! Close the geometry file
        close(unit=11)

        ! System call to plot planform
        call system('"C:\Program Files (x86)\ESPlot v1.3c\esplot.exe" planform.dat planform.qtp')
    end subroutine PlotPlanform
end module LiftingLinePlotting
