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

    subroutine PlotWashout(pf)
        type(Planform), intent(in) :: pf

        integer :: i

        open(unit=11, file='washout.dat')
        write(11, '(a)') "$ Dimensionless Washout Distribution"

        ! Write washout distribution
        do i = 1, pf%NNodes
            write(11, '(f22.15, a, 2x, f22.15)') &
                & z_over_b_i(i, pf%NNodes), ";", pf%Omega(i)
        end do

        close(unit=11)

        call system('"C:\Program Files (x86)\ESPlot v1.3c\esplot.exe" washout.dat washout.qtp')
    end subroutine PlotWashout

    subroutine PlotSectionLiftDistribution(pf)
        type(Planform), intent(in) :: pf

        integer :: i
        real*8 :: zb, cl(pf%NNodes)

        call GetLiftDistribution(pf, cl)

        open(unit=11, file='liftdistribution.dat')
        write(11, '(a)') "$ Section Lift Distribution"

        do i = 1, pf%NNodes
            zb = z_over_b_i(i, pf%NNodes)
            write(11, '(f22.15, a, 2x, f22.15)') zb, ";", cl(i)
        end do

        close(unit=11)

        call system('"C:\Program Files (x86)\ESPlot v1.3c\esplot.exe" liftdistribution.dat liftdistribution.qtp')
    end subroutine PlotSectionLiftDistribution

    subroutine PlotNormalizedLiftCoefficient(pf)
        type(Planform), intent(in) :: pf

        integer :: i
        real*8 :: zb, cb, cl(pf%NNodes)

        call GetLiftDistribution(pf, cl)

        open(unit=11, file='liftcoefficient.dat')
        write(11, '(a)') "$ Normalized Section Lift Coefficient"

        do i = 1, pf%NNodes
            zb = z_over_b_i(i, pf%NNodes)
            cb = c_over_b_zb(pf, zb)
            write(11, '(f22.15, a, 2x, f22.15)') zb, ";", cl(i) / cb / pf%CL1
        end do

        close(unit=11)

        call system('"C:\Program Files (x86)\ESPlot v1.3c\esplot.exe" liftcoefficient.dat liftcoefficient.qtp')
    end subroutine PlotNormalizedLiftCoefficient

    subroutine GetLiftDistribution(pf, cl)
        type(Planform), intent(in) :: pf
        real*8, intent(out) :: cl(pf%NNodes)

        integer :: i, j
        real*8 :: zb, theta

        do i = 1, pf%NNodes
            zb = z_over_b_i(i, pf%NNodes)
            theta = theta_zb(zb)
            cl(i) = 0.0d0
            do j = 1, pf%NNodes
                cl(i) = cl(i) + pf%BigA(j) * sin(real(j, 8) * theta)
            end do
            cl(i) = cl(i) * 4.0d0
        end do
    end subroutine GetLiftDistribution

end module LiftingLinePlotting
