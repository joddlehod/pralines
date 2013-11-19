module LiftingLineOutput
    use Utilities
    use class_Planform
    use matrix
    implicit none

contains
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

        ! Output flight conditions to console window
        call OutputOperatingConditions(6, pf)
        call OutputFlightCoefficients(6, pf)

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

        write(u, '(a15, 17x, 1x, a1, f20.15, 1x, a)') "Washout (twist)", &
            & "=", pf%Washout * 180.0d0 / pi, "degrees"
        write(u, '(a15, 17x, 1x, a1, f20.15, 1x, a)') "Optimum washout", &
            & "=", pf%OptimumWashout1 * 180.0d0 / pi, "degrees (Eq. 1.8.37)"
        if (pf%WashoutDistribution == Optimum) then
            write(u, '(a15, 17x, 1x, a1, f20.15, 1x, a)') "Optimum washout", &
                & "=", pf%OptimumWashout2 * 180.0d0 / pi, "degrees (Eq. 1.8.42)"
        end if
        write(u, '(a18, 14x, 1x, a1, f20.15, 1x, a)') &
            & "Aileron deflection", "=", &
            & pf%AileronDeflection * 180.0d0 / pi, "degrees"
        write(u, '(a26, 6x, 1x, a1, f20.15)') &
            & "Dimensionless rolling rate", "=", pf%RollingRate
        write(u, '(a32, 1x, a1, f20.15, 1x, a)') &
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
        write(u, '(2x, a9, 17x, 1x, a1, 3x, a)') "Wing type", "=", GetWingType(pf)
        write(u, fmt_str) "Number of nodes", "=", pf%NNodes, " (", &
            & (pf%NNodes + 1) / 2, " nodes per semispan)"
        write(u, '(2x, a26, 1x, a1, f20.15)') &
            & "Airfoil section lift slope", "=", pf%SectionLiftSlope
        write(u, '(2x, a12, 14x, 1x, a1, f20.15)') &
            & "Aspect Ratio", "=", pf%AspectRatio
        if (pf%WingType == Tapered) then
            write(u, '(2x, a11, 15x, 1x, a1, f20.15)') &
                & "Taper Ratio", "=", pf%TaperRatio
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

end module LiftingLineOutput
