module LiftingLineInterface
    use class_Planform
    use LiftingLineSetters
    use LiftingLineSolver
    use LiftingLineOutput
    use LiftingLinePlotting
    use LiftingLineSolver_Test

    implicit none

contains
    subroutine BeginLiftingLineInterface()
        type(Planform) :: pf
        character*2 :: inp = 'A'

        call InitPlanform(pf)

        do while(inp /= 'Q')
            inp = PlanformParamters(pf)
            if (inp == 'A') then
                call ComputeCMatrixAndCoefficients(pf)
                call OutputPlanform(pf)
                do while(inp /= 'Q' .and. inp /= 'B')
                    inp = OperatingConditions(pf)
                    if (inp /= 'Q') then
                        call UpdateOperatingConditions(pf, inp)
                    end if
                end do
            else if (inp /= 'Q') then
                call UpdatePlanformParameters(pf, inp)
            end if
        end do
    end subroutine BeginLiftingLineInterface

    character*2 function PlanformParamters(pf) result(inp)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        ! Clear the screen
        call system('cls')

        ! Display options to user
        write(6, '(a)') "Select from the following menu options:"

        ! Wing parameters
        write(6, '(2x, a)') "Wing Parameters:"

        msg = "WT - Edit wing type"
        call DisplayMessageWithTextDefault(msg, GetWingType(pf), 4)

        if (pf%WingType == Combination) then
            msg = "TZ - Edit z/b at the transition from tapered to elliptic"
            call DisplayMessageWithRealDefault(msg, pf%TransitionPoint, 4)

            msg = "TC - Edit c/croot at the transition from tapered to elliptic"
            call DisplayMessageWithRealDefault(msg, pf%TransitionChord, 4)
        end if

        if (pf%WingType /= Elliptic) then
            msg = "WD - Toggle washout distribution type"
            call DisplayMessageWithTextDefault(msg, GetWashoutDistributionType(pf), 4)
        end if

        msg = "N  - Edit number of nodes per semispan"
        call DisplayMessageWithIntegerDefault(msg, (pf%NNodes + 1) / 2, 4)

        msg = "RA - Edit aspect ratio"
        call DisplayMessageWithRealDefault(msg, pf%AspectRatio, 4)

        if (pf%WingType == Tapered) then
            msg = "RT - Edit taper ratio"
            call DisplayMessageWithRealDefault(msg, pf%TaperRatio, 4)
        end if

        msg = "S  - Edit section lift slope"
        call DisplayMessageWithRealDefault(msg, pf%SectionLiftSlope, 4)

        ! Aileron parameters
        write(6, *)
        write(6, '(2x, a)') "Aileron Parameters:"

        msg = "ZR - Edit z/b of aileron root"
        call DisplayMessageWithRealDefault(msg, pf%AileronRoot, 4)

        msg = "ZT - Edit z/b of aileron tip"
        call DisplayMessageWithRealDefault(msg, pf%AileronTip, 4)

        msg = "PH - Make hinge line parallel with quarter-chord line?"
        call DisplayMessageWithLogicalDefault(msg, pf%ParallelHingeLine, 4)

        msg = "CR - Edit cf/c of aileron root"
        call DisplayMessageWithRealDefault(msg, pf%FlapFractionRoot, 4)

        msg = "CT - Edit cf/c of aileron tip"
        call DisplayMessageWithRealDefault(msg, pf%FlapFractionTip, 4)

        msg = "HE - Edit aileron hinge efficiency"
        call DisplayMessageWithRealDefault(msg, pf%HingeEfficiency, 4)

        msg = "DE - Edit aileron deflection efficiency"
        call DisplayMessageWithRealDefault(msg, pf%DeflectionEfficiency, 4)

        ! Output and Plotting options
        write(6, *)
        write(6, '(2x, a)') "Output and Plotting Options:"
        msg = "C  - Output C matrix and Fourier Coefficients?"
        call DisplayMessageWithLogicalDefault(msg, pf%OutputMatrices, 4)

        msg = "F  - Edit output file name"
        call DisplayMessageWithTextDefault(msg, pf%FileName, 4)

        write(6, '(4x, a)') "PP - Plot planform in ES-Plot"

        ! Main Execution commands
        write(6, *)
        write(6, '(2x, a)') "A - Advance to operating conditions menu"
        write(6, '(2x, a)') "Q - Quit"

        write(6, *)
        write(6, '(a)') "Your selection: "

        inp = GetCharacterInput("  ")
        write(6, *)
    end function PlanformParamters

    character*2 function OperatingConditions(pf) result(inp)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        ! Clear the screen
        call system('cls')

        ! Output the Planform summary
        call OutputPlanformSummary(6, pf)
        write(6, *)

        ! Display options to user
        write(6, '(a)') "Select from the following menu options:"

        ! Operating Conditions
        write(6, *)
        write(6, '(2x, a)') "Operating Conditions:"

        msg = "AA - Edit root aerodynamic angle of attack"
        call DisplayMessageWithAngleDefault(msg, pf%AngleOfAttack, 4)

        msg = "CL - Edit coefficient of lift"
        call DisplayMessageWithRealDefault(msg, pf%LiftCoefficient, 4)

        msg = "OW - Use optimum total washout"
        call DisplayMessageWithLogicalDefault(msg, pf%UseOptimumWashout, 4)

        msg = "W  - Edit total amount of washout"
        call DisplayMessageWithAngleDefault(msg, pf%Washout, 4)

        msg = "AD - Edit aileron deflection"
        call DisplayMessageWithAngleDefault(msg, pf%AileronDeflection, 4)

        msg = "SR - Use steady dimensionless rolling rate"
        call DisplayMessageWithLogicalDefault(msg, pf%UseSteadyRollingRate, 4)

        msg = "R  - Edit dimensionless rolling rate"
        call DisplayMessageWithRealDefault(msg, pf%RollingRate, 4)

        ! Output and Plotting options
        write(6, *)
        write(6, '(2x, a)') "Output and Plotting Options:"
        msg = "S  - Save Flight coefficients to output file"
        call DisplayMessageWithTextDefault(msg, pf%FileName, 4)

        write(6, '(4x, a)') "PP - Plot Planform in ES-Plot"
        write(6, '(4x, a)') "PW - Plot Dimensionless Washout Distribution in ES-Plot"
        write(6, '(4x, a)') "PL - Plot Section Lift Distribution in ES-Plot"
        write(6, '(4x, a)') "PN - Plot Normalized Section Lift Coefficient in ES-Plot"

        ! Main Execution commands
        write(6, *)
        write(6, '(2x, a)') "B - Back to Planform Parameters"
        write(6, '(2x, a)') "Q - Quit"

        write(6, *)
        write(6, '(a)') "Your selection: "

        inp = GetCharacterInput("  ")
        write(6, *)
    end function OperatingConditions

    subroutine UpdatePlanformParameters(pf, input)
        type(Planform), intent(inout) :: pf
        character*2, intent(in) :: input

        ! Process input command
        ! Wing parameters
        if (input == 'WT') then
            call EditWingType(pf)
        else if (input == 'TZ' .and. pf%WingType == Combination) then
            call EditTransitionPoint(pf)
        else if (input == 'TC' .and. pf%WingType == Combination) then
            call EditTransitionChord(pf)
        else if (input == 'WD' .and. pf%WingType /= Elliptic) then
            call EditWashoutDistribution(pf)
        else if (input == 'N') then
            call EditNNodes(pf)
        else if (input == 'RA') then
            call EditAspectRatio(pf)
        else if (input == 'RT' .and. pf%WingType == Tapered) then
            call EditTaperRatio(pf)
        else if (input == 'S') then
            call EditLiftSlope(pf)

        ! Aileron parameters
        else if (input == 'ZR') then
            call EditAileronRoot(pf)
        else if (input == 'ZT') then
            call EditAileronTip(pf)
        else if (input == 'PH') then
            call ToggleParallelHinge(pf)
        else if (input == 'CR') then
            call EditFlapFractionRoot(pf)
        else if (input == 'CT') then
            call EditFlapFractionTip(pf)
        else if (input == 'HE') then
            call EditHingeEfficiency(pf)
        else if (input == 'DE') then
            call EditDeflectionEfficiency(pf)

        ! Output options
        else if (input == 'C') then
            pf%OutputMatrices = .not. pf%OutputMatrices
        else if (input == 'F') then
            call EditFileName(pf)
        else if (input == 'PP') then
            call PlotPlanform(pf)

        ! Testing options
        else if (input == 'T') then
            call TestLiftingLineSolver()
        end if
    end subroutine UpdatePlanformParameters

    subroutine UpdateOperatingConditions(pf, input)
        type(Planform), intent(inout) :: pf
        character*2, intent(in) :: input

        ! Operating Conditions
        if (input == 'AA') then
            call EditAngleOfAttack(pf)
        else if (input == 'CL') then
            call EditLiftCoefficient(pf)
        else if (input == 'OW') then
            call ToggleUseOptimumWashout(pf)
        else if (input == 'W') then
            call EditWashout(pf)
        else if (input == 'AD') then
            call EditAileronDeflection(pf)
        else if (input == 'SR') then
            call ToggleUseSteadyRollingRate(pf)
        else if (input == 'R') then
            call EditRollingRate(pf)

        ! Output and Plotting options
        else if (input == 'S') then
            call OutputFlightConditions(pf)
            call system('pause')
        else if (input == 'PP') then
            call PlotPlanform(pf)
        else if (input == 'PW') then
            call PlotWashout(pf)
        else if (input == 'PL') then
            call PlotSectionLiftDistribution(pf)
        else if (input == 'PN') then
            call PlotNormalizedLiftCoefficient(pf)
        end if

        call ComputeFlightConditions(pf)
    end subroutine UpdateOperatingConditions

    subroutine EditWingType(pf)
        type(Planform), intent(inout) :: pf

        character*2 :: inp

        write(6, *)
        write(6, '(a)') "Select from the following wing type options:"
        write(6, '(2x, a)') "T - Tapered"
        write(6, '(2x, a)') "E - Elliptic"
        write(6, '(2x, a)') "C - Combination (Tapered with elliptic tip)"
        write(6, *)
        write(6, '(a)') "Your selection: "

        inp = GetCharacterInput("  ")
        write(6, *)

        if (inp == "T") then
            call SetWingType(pf, Tapered)
        else if (inp == "E") then
            call SetWingType(pf, Elliptic)
        else if (inp == "C") then
            call SetWingType(pf, Combination)
        end if
    end subroutine EditWingType

    subroutine EditTransitionPoint(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg
        real*8 :: tp_old
        logical :: isValid

        tp_old = pf%TransitionPoint

        msg = "Enter z/b at the transition point from tapered to elliptic"
        call DisplayMessageWithRealDefault(msg, pf%TransitionPoint, 0)
        call SetTransitionPoint(pf, GetRealInput(0.0d0, 0.5d0, pf%TransitionPoint))

        isValid = AreCombinationWingCoefficientsValid(pf)
        do while(.not. isValid)
            call SetTransitionPoint(pf, tp_old)
            write(6, *)
            write(6, '(a)') "The input provided results in invalid ellipse coefficients."
            write(6, '(a)') "Try a new value or press <ENTER> to accept default."

            call SetTransitionPoint(pf, GetRealInput(0.0d0, 0.5d0, pf%TransitionPoint))
            isValid = AreCombinationWingCoefficientsValid(pf)
        end do
    end subroutine EditTransitionPoint

    subroutine EditTransitionChord(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg
        real*8 :: tc_old
        logical :: isValid

        tc_old = pf%TransitionChord

        msg = "Enter c/croot at the transition point from tapered to elliptic"
        call DisplayMessageWithRealDefault(msg, pf%TransitionChord, 0)
        call SetTransitionChord(pf, GetRealInput(0.0d0, 10.0d0, pf%TransitionChord))

        isValid = AreCombinationWingCoefficientsValid(pf)
        do while(.not. isValid)
            call SetTransitionChord(pf, tc_old)
            write(6, *)
            write(6, '(a)') "The input provided results in invalid ellipse coefficients."
            write(6, '(a)') "Try a new value or press <ENTER> to accept default."

            call SetTransitionChord(pf, GetRealInput(0.0d0, 2.0d0, pf%TransitionChord))
            isValid = AreCombinationWingCoefficientsValid(pf)
        end do
    end subroutine EditTransitionChord

    subroutine EditWashoutDistribution(pf)
        type(Planform), intent(inout) :: pf

        if (pf%WashoutDistribution == Linear) then
            call SetWashoutDistribution(pf, Optimum)
        else
            call SetWashoutDistribution(pf, Linear)
        end if
    end subroutine EditWashoutDistribution

    subroutine EditNNodes(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg
        integer :: npss
        character*80 :: int_str

        npss = (pf%NNodes + 1) / 2

        msg = "Enter number of nodes per semispan or press <ENTER> to accept default"
        call DisplayMessageWithIntegerDefault(msg, npss, 0)
        call SetNNodes(pf, GetIntInput(4, 1000, npss))
    end subroutine EditNNodes

    subroutine EditAspectRatio(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter new aspect ratio or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%AspectRatio, 0)
        call SetAspectRatio(pf, GetRealInput(1.0d0, 100.0d0, pf%AspectRatio))
    end subroutine EditAspectRatio

    subroutine EditTaperRatio(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter new taper ratio or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%TaperRatio, 0)
        call SetTaperRatio(pf, GetRealInput(0.0d0, 100.0d0, pf%TaperRatio))
    end subroutine EditTaperRatio

    subroutine EditLiftSlope(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter new section lift slope or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%SectionLiftSlope, 0)
        call SetSectionLiftSlope(pf, GetRealInput(-100.0d0 * pi, 100.0d0 * pi, &
            & pf%SectionLiftSlope))
    end subroutine EditLiftSlope

    subroutine EditAileronRoot(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter new z/b for aileron root or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%AileronRoot, 0)
        call SetAileronRoot(pf, GetRealInput(0.0d0, pf%AileronTip, pf%AileronRoot))
    end subroutine EditAileronRoot

    subroutine EditAileronTip(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter new z/b for aileron tip or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%AileronTip, 0)
        call SetAileronTip(pf, GetRealInput(pf%AileronRoot, 0.5d0, pf%AileronTip))
    end subroutine EditAileronTip

    subroutine EditFlapFractionRoot(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        if (pf%ParallelHingeLine) then
            write(6, '(a)') "NOTE: Aileron hinge is no longer constrained to be parallel with quarter-chord line."
        end if

        msg = "Enter new cf/c at aileron root or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%DesiredFlapFractionRoot, 0)
        call SetFlapFractionRoot(pf, GetRealInput(0.0d0, 1.0d0, pf%DesiredFlapFractionRoot))
    end subroutine EditFlapFractionRoot

    subroutine EditFlapFractionTip(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter new cf/c at aileron tip or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%FlapFractionTip, 0)
        call SetFlapFractionTip(pf, GetRealInput(0.0d0, 1.0d0, pf%FlapFractionTip))
    end subroutine EditFlapFractionTip

    subroutine ToggleParallelHinge(pf)
        type(Planform), intent(inout) :: pf

        if (pf%ParallelHingeLine) then
            call SetFlapFractionRoot(pf, pf%DesiredFlapFractionRoot)
        else
            call SetParallelHingeLine(pf)
        end if
    end subroutine ToggleParallelHinge

    subroutine EditHingeEfficiency(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter aileron hinge efficiency or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%HingeEfficiency, 0)
        call SetHingeEfficiency(pf, GetRealInput(0.0d0, 1.0d0, pf%HingeEfficiency))
    end subroutine EditHingeEfficiency

    subroutine EditDeflectionEfficiency(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter deflection efficiency or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%DeflectionEfficiency, 0)
        call SetDeflectionEfficiency(pf, GetRealInput(0.0d0, 1.0d0, &
            & pf%DeflectionEfficiency))
    end subroutine EditDeflectionEfficiency

    subroutine EditFileName(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter output file name or press <ENTER> to accept default"
        call DisplayMessageWithTextDefault(msg, pf%FileName, 0)
        call SetFileName(pf, GetStringInput(pf%FileName))
    end subroutine EditFileName

    subroutine EditAngleOfAttack(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        write(6, '(a, a)') "NOTE: This operation will calculate a new lift ", &
            & "coefficient and optimum washout."

        msg = "Enter angle of attack or press <ENTER> to accept default"
        call DisplayMessageWithAngleDefault(msg, pf%DesiredAngleOfAttack, 0)
        call SetAngleOfAttack(pf, GetRealInput(-12.0d0, 12.0d0, &
            & pf%DesiredAngleOfAttack * 180.0d0 / pi))
    end subroutine EditAngleOfAttack

    subroutine EditLiftCoefficient(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg
        real*8 :: mn, mx, dflt

        mn = CL1(pf%CLa, -12.0d0 * pi / 180.0d0, pf%EW, pf%Washout)
        mx = CL1(pf%CLa,  12.0d0 * pi / 180.0d0, pf%EW, pf%Washout)
        if (pf%DesiredLiftCoefficient < mn) then
            dflt = mn
        else if (pf%DesiredLiftCoefficient > mx) then
            dflt = mx
        else
            dflt = pf%DesiredLiftCoefficient
        end if

        write(6, *)
        write(6, '(a, a)') "NOTE: This operation will calculate a new alpha ", &
            & "and optimum washout"
        msg = "Enter lift coefficient or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, dflt, 0)
        call SetLiftCoefficient(pf, GetRealInput(mn, mx, dflt))
    end subroutine EditLiftCoefficient

    subroutine EditWashout(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        if (pf%UseOptimumWashout) then
            write(6, '(a)') "NOTE: Use of optimum total washout has been disabled."
        end if
        msg = "Enter total washout or press <ENTER> to accept default"
        call DisplayMessageWithAngleDefault(msg, pf%DesiredWashout, 0)
        call SetWashout(pf, GetRealInput(-12.0d0, 12.0d0, pf%DesiredWashout * 180.0d0 / pi))
    end subroutine EditWashout

    subroutine ToggleUseOptimumWashout(pf)
        type(Planform), intent(inout) :: pf

        if (pf%UseOptimumWashout) then
            call SetWashout(pf, pf%DesiredWashout * 180.0d0 / pi)
        else
            call SetOptimumWashout(pf)
        end if
    end subroutine ToggleUseOptimumWashout

    subroutine EditAileronDeflection(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter aileron deflection or press <ENTER> to accept default"
        call DisplayMessageWithAngleDefault(msg, pf%AileronDeflection, 0)
        call SetAileronDeflection(pf, GetRealInput(-12.0d0, 12.0d0, &
            & pf%AileronDeflection * 180.0d0 / pi))
    end subroutine EditAileronDeflection

    subroutine ToggleUseSteadyRollingRate(pf)
        type(Planform), intent(inout) :: pf

        pf%UseSteadyRollingRate = .not. pf%UseSteadyRollingRate
        if (pf%UseSteadyRollingRate) then
            call SetSteadyRollingRate(pf)
        else
            call SetRollingRate(pf, pf%DesiredRollingRate)
        end if
    end subroutine ToggleUseSteadyRollingRate

    subroutine EditRollingRate(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)

        if (pf%UseSteadyRollingRate) then
            write(6, '(a)') "NOTE: Use of steady rolling rate has been disabled."
        end if

        msg = "Enter dimensionless rolling rate or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%DesiredRollingRate, 0)
        call SetRollingRate(pf, GetRealInput(-100.0d0, 100.0d0, pf%DesiredRollingRate))
    end subroutine EditRollingRate

    subroutine DisplayMessageWithRealDefault(msg, dflt, tab)
        character*80, intent(in) :: msg ! Message to be displayed
        real*8, intent(in) :: dflt ! Default value to show in parenthesis
        integer, intent(in) :: tab ! Size of indentation to use

        call DisplayMessageWithTextDefault(msg, FormatReal(dflt, 5), tab)
    end subroutine DisplayMessageWithRealDefault

    subroutine DisplayMessageWithAngleDefault(msg, dflt, tab)
        character*80, intent(in) :: msg ! Message to be displayed
        real*8, intent(in) :: dflt ! Default value to show in parenthesis
        integer, intent(in) :: tab ! Size of indentation to use

        character*80 :: dflt_deg

        write(dflt_deg, '(a, a)') trim(FormatReal(dflt * 180.0d0 / pi, 5)), " degrees"
        call DisplayMessageWithTextDefault(msg, dflt_deg, tab)
    end subroutine DisplayMessageWithAngleDefault

    subroutine DisplayMessageWithIntegerDefault(msg, dflt, tab)
        character*80, intent(in) :: msg ! Message to be displayed
        integer, intent(in) :: dflt ! Default value to show in parenthesis
        integer, intent(in) :: tab ! Size of indentation to use

        call DisplayMessageWithTextDefault(msg, FormatInteger(dflt), tab)
    end subroutine DisplayMessageWithIntegerDefault

    subroutine DisplayMessageWithTextDefault(msg, dflt, tab)
        character*80, intent(in) :: msg ! Message to be displayed
        character*80, intent(in) :: dflt ! Default value to show in parenthesis
        integer, intent(in) :: tab ! Size of indentation to use

        character*80 :: msg_fmt

        if (tab == 0) then
            msg_fmt = "(a, a, a, a)"
        else
            write(msg_fmt, '(a, i1, a)') "(", tab, "x, a, a, a, a)"
        end if

        write(6, msg_fmt) trim(msg), " ( ", trim(dflt), " )"
    end subroutine DisplayMessageWithTextDefault

    subroutine DisplayMessageWithLogicalDefault(msg, dflt, tab)
        character*80, intent(in) :: msg ! Message to be displayed
        logical, intent(in) :: dflt ! Default value to show in parenthesis
        integer, intent(in) :: tab ! Size of indentation to use

        character*80 :: tf

        if (dflt) then
            tf = "True"
        else
            tf = "False"
        end if
        call DisplayMessageWithTextDefault(msg, tf, tab)

    end subroutine DisplayMessageWithLogicalDefault

end module LiftingLineInterface
