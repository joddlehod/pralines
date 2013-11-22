module LiftingLineOutput
    use Utilities
    use class_Planform
    use LiftingLineSetters
    use matrix
    implicit none

contains
    subroutine OutputHeader()
        integer :: i

        write(6, ('(80a)')) ("*", i=1, 80)
        write(6, '(34x, a)') "Pralines v1.0"
        write(6, *)
        write(6, '(28x, a)') "Author: Josh Hodson"
        write(6, '(28x, a)') "Release Date: 20 Nov 2013"
        write(6, *)
        write(6, ('(80a)')) ("*", i=1, 80)
    end subroutine OutputHeader

    subroutine OutputPlanform(pf)
        type(Planform), intent(in) :: pf

        ! Open a clean file for output
        open(unit=10, file=pf%FileName)

        ! Output the planform summary to output file
        call OutputPlanformSummary(10, pf)

        ! Output C matrix and fourier coefficients to output file
        if (pf%OutputMatrices) then
            call OutputC(10, pf%NNodes, pf%BigC)
            call OutputCInverse(10, pf%NNodes, pf%BigC_Inv)
            call OutputFourierCoefficients(10, pf)
        end if

        ! Close the output file
        close(unit=10)
    end subroutine OutputPlanform

    subroutine OutputLiftCoefficientParameters(u, pf)
        integer, intent(in) :: u  ! Output unit
        type(Planform), intent(in) :: pf

        write(u, '(a)') "Lift Coefficient Parameters:"
        write(u, '(2x, a, f20.15)') "KL    = ", pf%KL
        write(u, '(2x, a, f20.15)') "CL,a  = ", pf%CLa
        write(u, '(2x, a, f20.15)') "EW    = ", pf%EW
        write(u, *)
    end subroutine OutputLiftCoefficientParameters

    subroutine OutputDragCoefficientParameters(u, pf)
        integer, intent(in) :: u  ! Output unit
        type(Planform), intent(in) :: pf

        write(u, '(a)') "Drag Coefficient Parameters:"
        write(u, '(2x, a, f20.15)') "KD    = ", pf%KD
        write(u, '(2x, a, f20.15)') "KDL   = ", pf%KDL
        write(u, '(2x, a, f20.15)') "KDW   = ", pf%KDW
        write(u, '(2x, a, f20.15)') "es    = ", pf%ES
        write(u, *)
    end subroutine OutputDragCoefficientParameters

    subroutine OutputRollCoefficientParameters(u, pf)
        integer, intent(in) :: u  ! Output unit
        type(Planform), intent(in) :: pf

        write(u, '(a)') "Rolling Moment Coefficient Parameters:"
        write(u, '(2x, a, f20.15)') "Cl,da = ", pf%Crm_da
        write(u, '(2x, a, f20.15)') "Cl,pb = ", pf%Crm_pbar
        write(u, *)
    end subroutine OutputRollCoefficientParameters

    subroutine OutputFlightConditions(pf)
        type(Planform), intent(in) :: pf

        ! Open the file and append flight conditions to end
        open(unit=10, file=pf%FileName, access="append")

        ! Output flight conditions to output file
        call OutputOperatingConditions(10, pf)
        call OutputFlightCoefficients(10, pf)

        ! Close the output file
        close(unit=10)
    end subroutine OutputFlightConditions

    subroutine OutputOperatingConditions(u, pf)
        integer, intent(in) :: u  ! Output unit
        type(Planform), intent(in) :: pf

        write(u, '(a15, 19x, 1x, a1, f20.15, 1x, a)') "Optimum washout", &
            & "=", pf%OptimumWashout1 * 180.0d0 / pi, "degrees (Eq. 1.8.37)"
        if (pf%WashoutDistribution == Optimum) then
            write(u, '(a15, 19x, 1x, a1, f20.15, 1x, a)') "Optimum washout", &
                & "=", pf%OptimumWashout2 * 180.0d0 / pi, "degrees (Eq. 1.8.42)"
        end if
        write(u, '(a28, 6x, 1x, a1, f20.15, 1x, a)') "Washout used in calculations", &
            & "=", pf%Washout * 180.0d0 / pi, "degrees"
        write(u, '(a18, 16x, 1x, a1, f20.15, 1x, a)') &
            & "Aileron deflection", "=", &
            & pf%AileronDeflection * 180.0d0 / pi, "degrees"
        write(u, '(a33, 1x, 1x, a1, f20.15)') &
            & "Steady dimensionless rolling rate", "=", SteadyRollingRate(pf)
        write(u, '(a31, 3x, 1x, a1, f20.15)') &
            & "Dimensionless rolling rate used", "=", pf%RollingRate
        write(u, '(a32, 2x, 1x, a1, f20.15, 1x, a)') &
            & "Root aerodynamic angle of attack", "=", &
            & pf%AngleOfAttack * 180.0d0 / pi, "degrees"
        write(u, *)
    end subroutine OutputOperatingConditions

    subroutine OutputFlightCoefficients(u, pf)
        integer, intent(in) :: u  ! Output unit
        type(Planform), intent(in) :: pf

        write(u, '(a)') "Flight Coefficients:"
        write(u, '(2x, a, f20.15, a)') "CL    = ", pf%CL1, " (Eq. 1.8.24)"
        write(u, '(2x, a, f20.15, a)') "CL    = ", pf%CL2, " (Eq. 1.8.5)"
        write(u, '(2x, a, f20.15, a)') "CDi   = ", pf%CDi1, " (Eq. 1.8.25)"
        write(u, '(2x, a, f20.15, a)') "CDi   = ", pf%CDi2, " (Eq. 1.8.6)"
        write(u, '(2x, a, f20.15, a)') "CDi   = ", pf%CDi3, " (Exact)"
        write(u, '(2x, a, f20.15)') "Croll = ", pf%CRM
        write(u, '(2x, a, f20.15)') "Cyaw  = ", pf%CYM
        write(u, *)
    end subroutine OutputFlightCoefficients

    subroutine OutputPlanformSummary(u, pf)
        integer, intent(in) :: u
        type(Planform), intent(in) :: pf

        character*80 :: fmt_str
        integer :: len_nnodes

        len_nnodes = int(log10(real(pf%NNodes))) + 1
        write(fmt_str, '(a,i1,a,i1,a)') "(2x, a15, 11x, 1x, a1, 3x, i", &
            & len_nnodes, ", 1x, a, i", len_nnodes, ",a)"

        write(u, '(a)') "Planform Summary:"

        ! Wing type
        write(u, '(2x, a9, 17x, 1x, a1, 3x, a)') "Wing type", "=", &
            & trim(GetWingType(pf))

        ! Number of nodes
        write(u, fmt_str) "Number of nodes", "=", pf%NNodes, " (", &
            & (pf%NNodes + 1) / 2, " nodes per semispan)"

        ! Section Lift Slope
        write(u, '(2x, a26, 1x, a1, f20.15)') &
            & "Airfoil section lift slope", "=", pf%SectionLiftSlope

        ! Aspect Ratio
        write(u, '(2x, a12, 14x, 1x, a1, f20.15)') &
            & "Aspect Ratio", "=", pf%AspectRatio

        ! Taper Ratio
        if (pf%WingType == Tapered) then
            write(u, '(2x, a11, 15x, 1x, a1, f20.15)') &
                & "Taper Ratio", "=", pf%TaperRatio
        end if

        ! Transition from tapered to elliptic
        if (pf%WingType == Combination) then
            write(u, '(2x, a20, 6x, 1x, a1, f20.15)') "Transition Point z/b", &
                & "=", pf%TransitionPoint
            write(u, '(2x, a20, 6x, 1x, a1, f20.15)') "Transition Point c/croot", &
                & "=", pf%TransitionChord
        end if

        ! Washout distribution type
        if (pf%WingType /= Elliptic) then
            write(u, '(2x, a20, 6x, 1x, a1, 3x, a)') "Washout Distribution", &
                & "=", trim(GetWashoutDistributionType(pf))
        end if

        ! Location of aileron root, tip
        write(u, '(2x, a19, 7x, 1x, a1, f20.15)') &
            & "z/b at aileron root", "=", pf%AileronRoot
        write(u, '(2x, a18, 8x, 1x, a1, f20.15)') &
            & "z/b at aileron tip", "=", pf%AileronTip

        ! Flap fraction at aileron root, tip
        write(u, '(2x, a20, 6x, 1x, a1, f20.15)') &
            & "cf/c at aileron root", "=", pf%FlapFractionRoot
        write(u, '(2x, a19, 7x, 1x, a1, f20.15)') &
            & "cf/c at aileron tip", "=", pf%FlapFractionTip

        ! Hinge Efficiency Factor
        write(u, '(2x, a16, 10x, 1x, a1, f20.15)') &
            & "Hinge Efficiency", "=", pf%HingeEfficiency

        ! Deflection efficiency factor
        write(u, '(2x, a21, 5x, 1x, a1, f20.15)') &
            & "Deflection Efficiency", "=", pf%DeflectionEfficiency

        write(u, *)

        call OutputLiftCoefficientParameters(u, pf)
        call OutputDragCoefficientParameters(u, pf)
        call OutputRollCoefficientParameters(u, pf)
    end subroutine OutputPlanformSummary

    subroutine OutputFourierCoefficients(u, pf)
        integer, intent(in) :: u
        type(Planform), intent(in) :: pf

        integer :: i

        write(u, '(a)') "Fourier Coefficients:"
        write(u, '(a3, 4(2x, a20))') &
            & "i", "a(i)", "b(i)", "c(i)", "d(i)"
        do i = 1, pf%NNodes
            write(u, '(i3, 4(2x, f20.15))') &
                & i, pf%a(i), pf%b(i), pf%c(i), pf%d(i)
        end do
        write(u, *)
    end subroutine OutputFourierCoefficients

    subroutine OutputC(u, nnodes, c)
        integer, intent(in) :: u
        integer, intent(in) :: nnodes
        real*8, intent(in) :: c(nnodes, nnodes)

        write(u, *) "[C] Matrix:"
        call printmat(u, nnodes, nnodes, c)
        write(u, *)

    end subroutine OutputC

    subroutine OutputCInverse(u, nnodes, c_inv)
        integer, intent(in) :: u
        integer, intent(in) :: nnodes
        real*8, intent(in) :: c_inv(nnodes, nnodes)

        write(u, *) "[C]^-1 Matrix:"
        call printmat(u, nnodes, nnodes, c_inv)
        write(u, *)

    end subroutine OutputCInverse

    subroutine OutputHingeLine(u, pf)
        integer, intent(in) :: u
        type(Planform), intent(in) :: pf

        integer :: i

        do i = 1, pf%NNodes
            write(u, '(i3, 2x, f20.15, 2x, f20.15)') i, z_over_b_i(i, pf%NNodes), y_i(pf, i)
        end do
    end subroutine OutputHingeLine

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
        real*8 :: zb, cb, cl1, cl_over_cl, cl(pf%NNodes)

        ! Don't normalize if CL1 == 0
        if (Compare(pf%CL1, 0.0d0, zero) == 0) then
            cl1 = 1.0d0
        else
            cl1 = pf%CL1
        end if

        call GetLiftDistribution(pf, cl)

        open(unit=11, file='liftcoefficient.dat')
        write(11, '(a)') "$ Normalized Section Lift Coefficient"

        do i = 1, pf%NNodes
            zb = z_over_b_i(i, pf%NNodes)
            cb = c_over_b_zb(pf, zb)
            if (Compare(cb, 0.0d0, zero) == 0) then
                if (Compare(cl(i), 0.0d0, zero) == 0) then
                    if (pf%WingType == Elliptic) then
                        cl_over_cl = NLC_ZeroChord_Elliptic(pf, zb, cb, cl1)
                    else if (pf%WingType == Tapered) then
                        cl_over_cl = NLC_ZeroChord_Tapered(pf, zb, cb, cl1)
                    else if (pf%WingType == Combination) then
                        cl_over_cl = NLC_ZeroChord_Tapered(pf, zb, cb, cl1)
                    else
                        stop "***Unknown Wing Type***"
                    end if
                else
                    ! Finite lift from zero-chord section, should never happen
                    cl_over_cl = 1.0d0 / zero
                end if
            else
                cl_over_cl = cl(i) / cb / cl1
            end if
            write(11, '(f22.15, a, 2x, f22.15)') zb, ";", cl_over_cl
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

    real*8 function NLC_ZeroChord_Elliptic(pf, zb, cb, cl) result(cl_over_cl)
        type(Planform), intent(in) :: pf
        real*8, intent(in) :: zb
        real*8, intent(in) :: cb
        real*8, intent(in) :: cl

        integer :: i
        real*8 :: theta

        theta = theta_zb(zb)
        cl_over_cl = 0.0d0
        do i = 1, pf%NNodes
            cl_over_cl = cl_over_cl + real(i, 8) * pf%BigA(i) * &
                & cos(real(i, 8) * theta) / cos(theta)
        end do
        cl_over_cl = cl_over_cl * pi * pf%AspectRatio / cl
    end function NLC_ZeroChord_Elliptic

    real*8 function NLC_ZeroChord_Tapered(pf, zb, cb, cl) result(cl_over_cl)
        type(Planform), intent(in) :: pf
        real*8, intent(in) :: zb
        real*8, intent(in) :: cb
        real*8, intent(in) :: cl

        integer :: i
        real*8 :: theta, cb2

        if (zb < 0) then
            theta = 1.0d-5
        else
            theta = pi - 1.0d-5
        end if
        cb2 = c_over_b(pf, theta)
        cl_over_cl = 0.0d0
        do i = 1, pf%NNodes
            cl_over_cl = cl_over_cl + pf%BigA(i) * sin(real(i, 8) * theta)
        end do
        cl_over_cl = 4.0d0 * cl_over_cl / cb2 / cl
    end function NLC_ZeroChord_Tapered
end module LiftingLineOutput
