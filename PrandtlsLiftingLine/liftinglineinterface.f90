module LiftingLineInterface
    use class_Planform
    use LiftingLineSolver
    use LiftingLineOutput
    implicit none

contains
    subroutine InitApp()
        type(Planform) :: pf
        character*2 :: inp = 'A'

        call ComputeCMatrixAndCoefficients(pf)
        call ComputeFlightConditions(pf)

        do while(inp /= 'Q')
            inp = PlanformParamters(pf)
            if (inp == 'A') then
                call OutputPlanform(pf)
                do while(inp /= 'Q' .and. inp /= 'B')
                    inp = OperatingConditions(pf)
                    if (inp == 'E') then
                        call ComputeFlightConditions(pf)
                        call OutputFlightConditions(pf)
                        call system('pause')
                    else if (inp /= 'Q') then
                        call UpdateOperatingConditions(pf, inp)
                    end if
                end do
            else if (inp /= 'Q') then
                call UpdatePlanformParameters(pf, inp)
            end if
        end do
    end subroutine InitApp

    character*2 function PlanformParamters(pf) result(inp)
        type(Planform), intent(inout) :: pf

        ! Clear the screen
        call system('cls')

        ! Display options to user
        write(6, '(a)') "Select from the following menu options:"

        ! Wing parameters
        write(6, '(2x, a)') "Wing Parameters:"
        write(6, '(4x, a, a, a)') "W  - Toggle wing type ( ", &
            & trim(GetWingType(pf)), " )"
        write(6, '(4x, a, i3, a)') "N  - Edit number of nodes per semispan (", &
            & (pf%NNodes + 1) / 2, " )"
        write(6, '(4x, a, f7.4, a)') "RA - Edit aspect ratio (", &
            & pf%AspectRatio, " )"
        if (pf%WingType == Tapered) then
            write(6, '(4x, a, f7.4, a)') "RT - Edit taper ratio (", &
                & pf%TaperRatio, " )"
        end if
        write(6, '(4x, a, f11.7, a)') "S  - Edit section lift slope (", &
            & pf%LiftSlope, " )"
        write(6, *)

        ! Aileron parameters
        write(6, '(2x, a)') "Aileron Parameters:"
        write(6, '(4x, a, f7.4, a)') "AR - Edit location of aileron root (z/b = ", &
            & pf%AileronRoot, " )"
        write(6, '(4x, a, f7.4, a)') "AT - Edit location of aileron tip (z/b = ", &
            & pf%AileronTip, " )"
        write(6, '(4x, a, f7.4, a)') "FR - Edit the flap fraction at the aileron root (cf/c = ", &
            & pf%FlapFractionRoot, " )"
        write(6, '(4x, a, f7.4, a)') "FT - Edit the flap fraction at the aileron tip (cf/c = ", &
            & pf%FlapFractionTip, " )"
        write(6, '(4x, a)') "P  - Toggle calculation of flap fraction at aileron tip to make"
        write(6, '(4x, 4x, a, l1, a)') "hinge line parallel with quarter-chord line ( ", &
            & pf%ParallelHingeLine, " )"
        write(6, '(4x, a, f7.4, a)') "H  - Edit aileron hinge efficiency ( ", &
            & pf%HingeEfficiency, " )"
        write(6, '(4x, a, f7.4, a)') "D  - Edit aileron deflection efficiency ( ", &
            & pf%DeflectionEfficiency, " )"
        write(6, *)

        ! Output options
        write(6, '(2x, a)') "Output Options:"
        write(6, '(4x, a, l1, a)') "C  - Toggle output of C matrix and Fourier Coefficients ( ", &
            & pf%OutputMatrices, " )"
        write(6, '(4x, a, a, a)') "F  - Edit output file name ( ", &
            & trim(pf%FileName), " )"

        ! Main Execution commands
        write(6, *)
        write(6, '(2x, a)') "A - Advance to operating conditions menu"
        write(6, '(2x, a)') "Q - Quit"

        write(6, *)
        write(6, '(a)') "Your selection: "

        inp = GetCharacterInput()
        write(6, *)
    end function PlanformParamters

    character*2 function OperatingConditions(pf) result(inp)
        type(Planform), intent(inout) :: pf

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
        write(6, '(4x, a, f7.4, a)') "A - Edit root aerodynamic angle of attack (", &
            & pf%AngleOfAttack * 180.0d0 / pi, " deg )"
        write(6, '(4x, a, f7.4, a)') "L - Edit lift coefficient (", &
            & pf%LiftCoefficient, " )"
        write(6, '(4x, a, f7.4, a)') "W - Edit amount of linear washout (", &
            & pf%Omega * 180.0d0 / pi, " deg )"
        write(6, '(4x, a, f7.4, a)') "D - Edit the aileron deflection (", &
            & pf%AileronDeflection * 180.0d0 / pi, " deg )"
        write(6, '(4x, a, f7.4, a)') "R - Edit dimensionless rolling rate (", &
            & pf%RollingRate, " )"

        ! Main Execution commands
        write(6, *)
        write(6, '(2x, a)') "B - Back to Planform Parameters"
        write(6, '(2x, a)') "E - Execute Solver for specified operating conditions"
        write(6, '(2x, a)') "Q - Quit"

        write(6, *)
        write(6, '(a)') "Your selection: "

        inp = GetCharacterInput()
        write(6, *)
    end function OperatingConditions

    subroutine UpdatePlanformParameters(pf, input)
        type(Planform), intent(inout) :: pf
        character*2, intent(in) :: input

        ! Process input command
        ! Wing parameters
        if (input == 'W') then
            call ToggleWingType(pf)
        else if (input == 'N') then
            call EditNodes(pf)
        else if (input == 'RA') then
            call EditAspectRatio(pf)
        else if (input == 'RT' .and. pf%WingType == Tapered) then
            call EditTaperRatio(pf)
        else if (input == 'S') then
            call EditLiftSlope(pf)

        ! Aileron parameters
        else if (input == 'AR') then
            call EditAileronRoot(pf)
        else if (input == 'AT') then
            call EditAileronTip(pf)
        else if (input == 'FR') then
            call EditFlapFractionRoot(pf)
        else if (input == 'FT') then
            call EditFlapFractionTip(pf)
        else if (input == 'P') then
            call ToggleParallelHinge(pf)
        else if (input == 'H') then
            call EditHingeEfficiency(pf)
        else if (input == 'D') then
            call EditDeflectionEfficiency(pf)

        ! Output options
        else if (input == 'C') then
            pf%OutputMatrices = .not. pf%OutputMatrices
        else if (input == 'F') then
            call EditFileName(pf)
        end if

        call ComputeCMatrixAndCoefficients(pf)
    end subroutine UpdatePlanformParameters

    subroutine UpdateOperatingConditions(pf, input)
        type(Planform), intent(inout) :: pf
        character*2, intent(in) :: input

        ! Operating Conditions
        if (input == 'A') then
            call EditAngleOfAttack(pf)
        else if (input == 'L') then
            call EditLiftCoefficient(pf)
        else if (input == 'W') then
            call EditOmega(pf)
        else if (input == 'D') then
            call EditAileronDeflection(pf)
        else if (input == 'R') then
            call EditRollingRate(pf)
        end if

        call ComputeFlightConditions(pf)
    end subroutine UpdateOperatingConditions

    subroutine ToggleWingType(pf)
        type(Planform), intent(inout) :: pf
        if (pf%WingType == Tapered) then
            pf%WingType = Elliptic
        else if (pf%WingType == Elliptic) then
            pf%WingType = Tapered
        end if
    end subroutine ToggleWingType

    subroutine EditNodes(pf)
        type(Planform), intent(inout) :: pf

        integer :: nnodes

        write(6, *)
        write(6, '(a)') "Enter number of nodes per semispan:"

        nnodes = GetIntInput()
        if (nnodes > 0) then
            pf%NNodes = nnodes * 2 - 1
        end if
    end subroutine EditNodes

    subroutine EditAspectRatio(pf)
        type(Planform), intent(inout) :: pf

        real*8 :: aratio

        write(6, *)
        write(6, '(a)') "Enter aspect ratio:"

        aratio = GetRealInput()
        if (aratio > 0.0) then
            pf%AspectRatio = aratio
        end if
    end subroutine EditAspectRatio

    subroutine EditTaperRatio(pf)
        type(Planform), intent(inout) :: pf

        real*8 :: tratio

        write(6, *)
        write(6, '(a)') "Enter taper ratio:"

        tratio = GetRealInput()
        if (tratio >= 0.0) then
            pf%TaperRatio = tratio
        end if
    end subroutine EditTaperRatio

    subroutine EditLiftSlope(pf)
        type(Planform), intent(inout) :: pf

        write(6, *)
        write(6, '(a)') "Enter lift slope:"

        pf%LiftSlope = GetRealInput()
    end subroutine EditLiftSlope

    subroutine EditAileronRoot(pf)
        type(Planform), intent(inout) :: pf

        write(6, *)
        write(6, '(a)') "Enter z/b for aileron root location:"

        pf%AileronRoot = GetRealInput()
    end subroutine EditAileronRoot

    subroutine EditAileronTip(pf)
        type(Planform), intent(inout) :: pf

        write(6, *)
        write(6, '(a)') "Enter z/b for aileron tip location:"

        pf%AileronTip = GetRealInput()
    end subroutine EditAileronTip

    subroutine EditFlapFractionRoot(pf)
        type(Planform), intent(inout) :: pf

        write(6, *)
        write(6, '(a)') "Enter cf/c at aileron root:"

        pf%FlapFractionRoot = GetRealInput()
    end subroutine EditFlapFractionRoot

    subroutine EditFlapFractionTip(pf)
        type(Planform), intent(inout) :: pf

        write(6, *)
        if (pf%ParallelHingeLine) then
            pf%ParallelHingeLine = .false.
            write(6, '(a)') "NOTE: Calculation of flap fraction at aileron tip has been turned off."
        end if
        write(6, '(a)') "Enter cf/c at aileron tip:"

        pf%FlapFractionTip = GetRealInput()
    end subroutine EditFlapFractionTip

    subroutine ToggleParallelHinge(pf)
        type(Planform), intent(inout) :: pf

        pf%ParallelHingeLine = .not. pf%ParallelHingeLine
    end subroutine ToggleParallelHinge

    subroutine EditHingeEfficiency(pf)
        type(Planform), intent(inout) :: pf

        write(6, *)
        write(6, '(a)') "Enter aileron hinge efficiency:"

        pf%HingeEfficiency = GetRealInput()
    end subroutine EditHingeEfficiency

    subroutine EditDeflectionEfficiency(pf)
        type(Planform), intent(inout) :: pf

        write(6, *)
        write(6, '(a)') "Enter aileron deflection efficiency:"

        pf%DeflectionEfficiency = GetRealInput()
    end subroutine EditDeflectionEfficiency

    subroutine EditFileName(pf)
        type(Planform), intent(inout) :: pf

        write(6, *)
        write(6, '(a)') "Enter name of output file:"

        read(5, '(a)') pf%FileName
    end subroutine EditFileName

    subroutine EditAngleOfAttack(pf)
        type(Planform), intent(inout) :: pf

        write(6, *)
        write(6, '(a)') "NOTE: This operation will calculate a new lift coefficient."
        write(6, '(a,a,f7.4,a)') "Enter desired root aerodynamic angle of ", &
            " attack in degrees (", pf%AngleOfAttack * 180.0d0 / pi, " ):"

        pf%DesiredAngleOfAttack = GetRealInput() * pi / 180.0d0
        pf%AngleOfAttack = pf%DesiredAngleOfAttack
        pf%SpecifyAlpha = .true.
        call ComputeFlightConditions(pf)
    end subroutine EditAngleOfAttack

    subroutine EditLiftCoefficient(pf)
        type(Planform), intent(inout) :: pf

        write(6, *)
        write(6, '(a)') "NOTE: This operation will calculate a new root aerodynamic angle of attack."
        write(6, '(a, f7.4, a)') "Enter desired lift coefficient (", &
            pf%DesiredLiftCoefficient, " ):"

        pf%DesiredLiftCoefficient = GetRealInput()
        pf%LiftCoefficient = pf%DesiredLiftCoefficient
        pf%SpecifyAlpha = .false.
        call ComputeFlightConditions(pf)
    end subroutine EditLiftCoefficient

    subroutine EditOmega(pf)
        type(Planform), intent(inout) :: pf

        write(6, *)
        write(6, '(a)') "Enter amount of linear twist (degrees):"

        pf%Omega = GetRealInput() * pi / 180.0d0
    end subroutine EditOmega

    subroutine EditAileronDeflection(pf)
        type(Planform), intent(inout) :: pf

        write(6, *)
        write(6, '(a)') "Enter the aileron deflection (degrees):"

        pf%AileronDeflection = GetRealInput() * pi / 180.0d0
    end subroutine EditAileronDeflection

    subroutine EditRollingRate(pf)
        type(Planform), intent(inout) :: pf

        write(6, *)
        write(6, '(a)') "Enter the dimensionless rolling rate:"

        pf%RollingRate = GetRealInput()
    end subroutine EditRollingRate

    character*2 function GetCharacterInput() result(inp)
        integer :: i
        character :: a

        read(5, '(a)') inp
        do i = 1, 2
            a = inp(i:i)
            if(iachar(a) >= iachar('a') .and. iachar(a) <= iachar('z')) then
                inp(i:i) = char(iachar(a) - 32)
            end if
        end do
    end function GetCharacterInput

    integer function GetIntInput() result(input)
        read(5, *) input
    end function GetIntInput

    real*8 function GetRealInput() result(input)
        read(5, *) input
    end function GetRealInput

end module LiftingLineInterface
