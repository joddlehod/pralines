module LiftingLineInterface
    use class_Planform
    use LiftingLineSolver
    use LiftingLineOutput
    implicit none

contains
    subroutine BeginLiftingLineInterface()
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
        write(6, '(4x, a, a, a)') "W  - Toggle wing type ( ", &
            & trim(GetWingType(pf)), " )"
        write(6, '(4x, a, i3, a)') "N  - Edit number of nodes per semispan (", &
            & (pf%NNodes + 1) / 2, " )"

        msg = "RA - Edit aspect ratio"
        call DisplayMessageWithRealDefault(msg, pf%AspectRatio, 4)

        if (pf%WingType == Tapered) then
            msg = "RT - Edit taper ratio"
            call DisplayMessageWithRealDefault(msg, pf%TaperRatio, 4)
        end if

        msg = "S  - Edit section lift slope"
        call DisplayMessageWithRealDefault(msg, pf%LiftSlope, 4)
        write(6, *)

        ! Aileron parameters
        write(6, '(2x, a)') "Aileron Parameters:"

        msg = "AR - Edit z/b of aileron root"
        call DisplayMessageWithRealDefault(msg, pf%AileronRoot, 4)

        msg = "AT - Edit z/b of aileron tip"
        call DisplayMessageWithRealDefault(msg, pf%AileronTip, 4)

        msg = "FR - Edit cf/c of aileron root"
        call DisplayMessageWithRealDefault(msg, pf%FlapFractionRoot, 4)

        msg = "FT - Edit cf/c of aileron tip"
        call DisplayMessageWithRealDefault(msg, pf%FlapFractionTip, 4)

        write(6, '(4x, a)') "P  - Toggle calculation of flap fraction at aileron root to make"
        write(6, '(4x, 5x, a, l1, a)') "hinge line parallel with quarter-chord line ( ", &
            & pf%ParallelHingeLine, " )"

        msg = "H  - Edit aileron hinge efficiency"
        call DisplayMessageWithRealDefault(msg, pf%HingeEfficiency, 4)

        msg = "D  - Edit aileron deflection efficiency"
        call DisplayMessageWithRealDefault(msg, pf%DeflectionEfficiency, 4)

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

        inp = GetStringInput("  ")
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

        msg = "A - Edit root aerodynamic angle of attack"
        call DisplayMessageWithAngleDefault(msg, pf%AngleOfAttack, 4)

        msg = "L - Edit lift coefficient"
        call DisplayMessageWithRealDefault(msg, pf%LiftCoefficient, 4)

        msg = "W - Edit amount of linear washout"
        call DisplayMessageWithAngleDefault(msg, pf%Omega, 4)

        msg = "D - Edit the aileron deflection"
        call DisplayMessageWithAngleDefault(msg, pf%AileronDeflection, 4)

        msg = "R - Edit dimensionless rolling rate"
        call DisplayMessageWithRealDefault(msg, pf%RollingRate, 4)

        ! Main Execution commands
        write(6, *)
        write(6, '(2x, a)') "B - Back to Planform Parameters"
        write(6, '(2x, a)') "E - Execute Solver for specified operating conditions"
        write(6, '(2x, a)') "Q - Quit"

        write(6, *)
        write(6, '(a)') "Your selection: "

        inp = GetStringInput("  ")
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

        integer :: npss
        character*80 :: int_str

        npss = (pf%NNodes + 1) / 2

        write(6, *)
        write(6, '(a, a, a, a)') "Enter number of nodes per semispan ", &
            & "( ", trim(FormatInteger(npss)), " ) :"

        npss = GetIntInput(4, 1000, npss)
        pf%NNodes = npss * 2 - 1
    end subroutine EditNodes

    subroutine EditAspectRatio(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter new aspect ratio or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%AspectRatio, 0)
        pf%AspectRatio = GetRealInput(1.0d0, 100.0d0, pf%AspectRatio)
    end subroutine EditAspectRatio

    subroutine EditTaperRatio(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter new taper ratio or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%TaperRatio, 0)
        pf%TaperRatio = GetRealInput(0.0d0, 100.0d0, pf%TaperRatio)
    end subroutine EditTaperRatio

    subroutine EditLiftSlope(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter new section lift slope or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%LiftSlope, 0)
        pf%LiftSlope = GetRealInput(-100.0d0 * pi, 100.0d0 * pi, pf%LiftSlope)
    end subroutine EditLiftSlope

    subroutine EditAileronRoot(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter new z/b for aileron root or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%AileronRoot, 0)
        pf%AileronRoot = GetRealInput(0.0d0, pf%AileronTip, pf%AileronRoot)
    end subroutine EditAileronRoot

    subroutine EditAileronTip(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter new z/b for aileron tip or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%AileronTip, 0)
        pf%AileronTip = GetRealInput(pf%AileronRoot, 0.5d0, pf%AileronTip)
    end subroutine EditAileronTip

    subroutine EditFlapFractionRoot(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        if (pf%ParallelHingeLine) then
            pf%ParallelHingeLine = .false.
            write(6, '(a)') "NOTE: Calculation of flap fraction at aileron root has been turned off."
        end if

        msg = "Enter new cf/c at aileron root or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%FlapFractionRoot, 0)
        pf%FlapFractionRoot = GetRealInput(0.0d0, 1.0d0, pf%FlapFractionRoot)
    end subroutine EditFlapFractionRoot

    subroutine EditFlapFractionTip(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter new cf/c at aileron tip or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%FlapFractionTip, 0)
        pf%FlapFractionTip = GetRealInput(0.0d0, 1.0d0, pf%FlapFractionTip)
    end subroutine EditFlapFractionTip

    subroutine ToggleParallelHinge(pf)
        type(Planform), intent(inout) :: pf

        pf%ParallelHingeLine = .not. pf%ParallelHingeLine
    end subroutine ToggleParallelHinge

    subroutine EditHingeEfficiency(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter aileron hinge efficiency or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%HingeEfficiency, 0)
        pf%HingeEfficiency = GetRealInput(0.0d0, 1.0d0, pf%HingeEfficiency)
    end subroutine EditHingeEfficiency

    subroutine EditDeflectionEfficiency(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter new aileron deflection efficiency or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%DeflectionEfficiency, 0)
        pf%DeflectionEfficiency = GetRealInput(0.0d0, 1.0d0, pf%DeflectionEfficiency)
    end subroutine EditDeflectionEfficiency

    subroutine EditFileName(pf)
        type(Planform), intent(inout) :: pf

        write(6, *)
        write(6, '(a)') "Enter name of output file:"

        read(5, '(a)') pf%FileName
    end subroutine EditFileName

    subroutine EditAngleOfAttack(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        write(6, '(a)') "NOTE: This operation will calculate a new lift coefficient."

        msg = "Enter new angle of attack or press <ENTER> to accept default"
        call DisplayMessageWithAngleDefault(msg, pf%AngleOfAttack, 0)
        pf%DesiredAngleOfAttack = GetRealInput(-12.0d0, 12.0d0, &
            & pf%DesiredAngleOfAttack * 180.0d0 / pi) * pi / 180.0d0
        pf%AngleOfAttack = pf%DesiredAngleOfAttack
        pf%SpecifyAlpha = .true.
    end subroutine EditAngleOfAttack

    subroutine EditLiftCoefficient(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        write(6, '(a)') "NOTE: This operation will calculate a new root aerodynamic angle of attack."
        msg = "Enter a new lift coefficient or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%DesiredLiftCoefficient, 0)
        pf%DesiredLiftCoefficient = GetRealInput(-100.0d0, 100.0d0, &
            & pf%DesiredLiftCoefficient)
        pf%LiftCoefficient = pf%DesiredLiftCoefficient
        pf%SpecifyAlpha = .false.
    end subroutine EditLiftCoefficient

    subroutine EditOmega(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter amount of linear washout or press <ENTER> to accept default"
        call DisplayMessageWithAngleDefault(msg, pf%Omega, 0)
        pf%Omega = GetRealInput(-100.0d0, 100.0d0, pf%Omega * 180.0d0 / pi) * &
            & pi / 180.0d0
    end subroutine EditOmega

    subroutine EditAileronDeflection(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter the aileron deflection or press <ENTER> to accept default"
        call DisplayMessageWithAngleDefault(msg, pf%AileronDeflection, 0)
        pf%AileronDeflection = GetRealInput(-12.0d0, 12.0d0, &
            & pf%AileronDeflection * 180.0d0 / pi) * pi / 180.0d0
    end subroutine EditAileronDeflection

    subroutine EditRollingRate(pf)
        type(Planform), intent(inout) :: pf

        character*80 :: msg

        write(6, *)
        msg = "Enter the dimensionless rolling rate or press <ENTER> to accept default"
        call DisplayMessageWithRealDefault(msg, pf%RollingRate, 0)
        pf%RollingRate = GetRealInput(-100.0d0, 100.0d0, pf%RollingRate)
    end subroutine EditRollingRate

    character*2 function GetStringInput(def) result(inp)
        character*2, intent(in) :: def  ! Default value if invalid input

        integer :: i
        character :: a

        read(5, '(a)') inp
        if (len(inp) > 2 .or. len(inp) < 1) then
            inp = def
        else
            do i = 1, 2
                a = inp(i:i)
                if(iachar(a) >= iachar('a') .and. iachar(a) <= iachar('z')) then
                    inp(i:i) = char(iachar(a) - 32)
                end if
            end do
        end if
    end function GetStringInput

    integer function GetIntInput(mn, mx, def) result(inp)
        integer, intent(in) :: mn   ! Minimum accepted value
        integer, intent(in) :: mx   ! Maximum accepted value
        integer, intent(in) :: def  ! Default value if invalid input

        logical :: cont
        character*80 :: inp_str
        integer :: ios
        integer :: len_mn, len_mx
        character*80 :: msg_fmt

        cont = .true.
        do while (cont)
            read(5, '(a)', iostat=ios) inp_str
            if (ios == 0 .and. trim(inp_str) /= "") then
                read(inp_str, *, iostat=ios) inp
                if (ios /= 0 .or. inp < mn .or. inp > mx) then
                    len_mn = int(log10(real(abs(mn)))) + 1
                    if (mn < 0) len_mn = len_mn + 1

                    len_mx = int(log10(real(abs(mx)))) + 1
                    if (mx < 0) len_mx = len_mx + 1

                    write(msg_fmt, '(a, i1, a, i1, a)') "(a, a, i", len_mn, &
                        & ", a, i", len_mx, ", a)"

                    write(6, *)
                    write(6, msg_fmt) "Invalid input. Please ", &
                        & "specify an integer between ", mn, " and ", mx, ","
                    write(6, '(a)') "or press <ENTER> to accept the default value."
                else
                    cont = .false.
                end if
            else
                inp = def
                cont = .false.
            end if
        end do
    end function GetIntInput

    subroutine DisplayMessageWithRealDefault(msg, dflt, tab)
        character*80, intent(in) :: msg ! Message to be displayed
        real*8, intent(in) :: dflt ! Default value to show in parenthesis
        integer, intent(in) :: tab ! Size of indentation to use

        integer :: len_inp
        character*80 :: msg_fmt

        if (tab == 0) then
            msg_fmt = "(a, a, a, a)"
        else
            write(msg_fmt, '(a, i1, a)') "(", tab, "x, a, a, a, a)"
        end if

        write(6, msg_fmt) trim(msg), " ( ", &
            & trim(FormatReal(dflt, 5)), " )"
    end subroutine DisplayMessageWithRealDefault

    subroutine DisplayMessageWithAngleDefault(msg, dflt, tab)
        character*80, intent(in) :: msg ! Message to be displayed
        real*8, intent(in) :: dflt ! Default value to show in parenthesis
        integer, intent(in) :: tab ! Size of indentation to use

        integer :: len_inp
        character*80 :: msg_fmt
        real*8 :: dflt_deg

        dflt_deg = dflt * 180.0d0 / pi

        if (tab == 0) then
            msg_fmt = "(a, a, a, a)"
        else
            write(msg_fmt, '(a, i1, a)') "(", tab, "x, a, a, a, a)"
        end if

        write(6, msg_fmt) trim(msg), " ( ", &
            & trim(FormatReal(dflt_deg, 5)), " degrees )"
    end subroutine DisplayMessageWithAngleDefault

    real*8 function GetRealInput(mn_orig, mx_orig, dflt_orig) result(inp)
        real*8, intent(in) :: mn_orig   ! Minimum accepted value
        real*8, intent(in) :: mx_orig   ! Maximum accepted value
        real*8, intent(in) :: dflt_orig ! Default value for input

        logical :: cont
        character*80 :: inp_str
        integer :: ios
        integer :: len_mn, ndec_mn
        integer :: len_mx, ndec_mx
        character*80 :: msg_fmt
        real*8 :: mn, mx, dflt

        if (Compare(mn_orig, 0.0d0, zero) == 0) then
            mn = 0.0d0
        else
            mn = mn_orig
        end if

        if (Compare(mx_orig, 0.0d0, zero) == 0) then
            mx = 0.0d0
        else
            mx = mx_orig
        end if

        if (Compare(dflt_orig, 0.0d0, zero) == 0) then
            dflt = 0.0d0
        else
            dflt = dflt_orig
        end if

        cont = .true.
        do while (cont)
            read(5, '(a)', iostat=ios) inp_str
            if (ios == 0 .and. trim(inp_str) /= "") then
                read(inp_str, *, iostat=ios) inp
                if (ios /= 0 .or. inp < mn .or. inp > mx) then
                    write(6, *)
                    write(6, '(a, a, a, a, a, a)') "Invalid input. Please ", &
                        & "enter a number between ", trim(FormatReal(mn, 5)), &
                        & " and ", trim(FormatReal(mx, 5)), ","
                    write(6, '(a)') "or press <ENTER> to accept the default value."
                else
                    cont = .false.
                end if
            else
                inp = dflt
                cont = .false.
            end if
        end do
    end function GetRealInput

    character*80 function FormatReal(r, ndigits) result(real_str)
        real*8, intent(in) :: r
        integer, intent(in) :: ndigits

        integer :: order, width, ndecimal
        character*80 :: real_fmt
        real*8 :: r_div_pi
        integer :: num, denom

        if (Compare(r, 0.0d0, zero) /= 0 .and. (IsFactorOfPi(r, ndigits) &
            & .or. IsFractionOfPi(r, num, denom))) then
            if (Compare(r, pi, zero) == 0) then
                write(real_str, '(a)') "PI"
            else if (IsFractionOfPi(r, num, denom)) then
                if (denom == 1) then
                    write(real_str, '(a, a)') trim(FormatInteger(num)), "*PI"
                else
                    write(real_str, '(a, a, a, a)') trim(FormatInteger(num)), &
                        & "/", trim(FormatInteger(denom)), "*PI"
                end if
            else
                r_div_pi = r / pi
                write(real_str, '(a, a)') trim(FormatReal(r_div_pi, ndigits)), &
                    & "*PI )"
            end if
        else
            if (Compare(r, 0.0d0, zero) == 0) then
                order = 1
            else
                ! Determine the location of the first non-zero digit in the number
                order = int(log10(real(abs(r), 8))) + 1
            end if

            ! Check for sizes that should use exponential format
            if (order <= -4 .or. order >= ndigits) then
                if (r < 0.0d0) then
                    width = ndigits + 6  ! e.g. -1.2345E+67 - 5 digits + 6 other
                else
                    width = ndigits + 5  ! e.g. 1.2345E-67 - 5 digits + 5 other
                end if

                write(real_fmt, '(a, i2, a, i2, a)') "(ES", width, ".", &
                    & ndigits - 1, ")"
            else
                if (r < 0.0d0) then
                    width = ndigits + 2  ! e.g. -12.345 - 5 digits + 2 other
                else
                    width = ndigits + 1  ! e.g. 123.45 - 5 digits + 1 other
                end if

                if (order < 0) then
                    width = width + 1  ! e.g. -0.012345 - additional for leading 0
                end if

                write(real_fmt, '(a, i2, a, i2, a)') "(F", width, ".", &
                    & ndigits - order, ")"
            end if
            write(real_str, real_fmt) r
        end if
    end function FormatReal

    character*80 function FormatInteger(i) result(int_str)
        integer, intent(in) :: i

        integer :: len_i
        character*80 :: int_fmt

        if (i == 0) then
            len_i = 1
        else
            len_i = int(log10(real(abs(i)))) + 1
            if (i < 0) then
                len_i = len_i + 1
            end if
        end if

        write(int_fmt, '(a, i2, a)') "(i", len_i, ")"
        write(int_str, int_fmt) i
    end function FormatInteger

    logical function IsFactorOfPi(r, ndigits) result(isFactor)
        real*8, intent(in) :: r
        integer, intent(in) :: ndigits

        real*8 :: rx, rx_trunc

        rx = r / pi * 10**ndigits
        rx_trunc = real(int(rx), 8)

        if (Compare(rx, rx_trunc, zero) == 0) then
            isFactor = .true.
        else
            isFactor = .false.
        end if
    end function IsFactorOfPi

    logical function IsFractionOfPi(r, num, denom) result(isFraction)
        real*8, intent(in) :: r
        integer, intent(out) :: num
        integer, intent(out) :: denom

        integer :: i, num2, denom2
        real*8 :: r_div_pi, r_div_pi_i

        isFraction = .false.
        r_div_pi = r / pi
        do i = 1, 360
            r_div_pi_i = r_div_pi * real(i, 8)
            if (Compare(r_div_pi_i , real(int(r_div_pi_i), 8), zero) == 0) then
                num = int(r_div_pi_i)
                denom = i
                isFraction = .true.
                return
            end if
        end do
    end function IsFractionOfPi

end module LiftingLineInterface
