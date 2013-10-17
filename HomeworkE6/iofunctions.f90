module iofunctions
    use class_Planform
    use liftinglinesolver
    implicit none

contains
    subroutine MainPage(pf)
        type(Planform), intent(inout) :: pf

        character :: input
        logical :: cont = .true.

        do while (cont)
            ! Clear the screen
            call system('cls')

            ! Display options to user
            write(6, '(a)') "Select from the following options:"

            ! Wing parameters
            write(6, '(2x, a)') "Wing Parameters:"
            write(6, '(4x, a, a, a)') "W - Toggle wing type ( ", &
                & trim(GetWingType(pf)), " )"
            write(6, '(4x, a, i3, a)') "N - Edit number of nodes per semispan (", &
                & (pf%NNodes + 1) / 2, " )"
            write(6, '(4x, a, f7.4, a)') "A - Edit aspect ratio (", &
                & pf%AspectRatio, " )"
            if (pf%WingType == Tapered) then
                write(6, '(4x, a, f7.4, a)') "T - Edit taper ratio (", &
                    & pf%TaperRatio, " )"
            end if
            write(6, '(4x, a, f11.7, a)') "S - Edit section lift slope (", &
                & pf%LiftSlope, " )"

            ! Output options
            write(6, *)
            write(6, '(2x, a)') "Output Options:"
            write(6, '(4x, a, l1, a)') "C - Toggle output of C matrix ( ", &
                & pf%WriteCMatrix, " )"
            write(6, '(4x, a, l1, a)') "I - Toggle output of C-inverse matrix ( ", &
                & pf%WriteCInverse, " )"
            write(6, '(4x, a, l1, a)') "F - Toggle output of Fourier Coefficients ( ", &
                & pf%WriteFourier, " )"
            write(6, '(4x, a, a, l1, a)') "K - Toggle output of KL, KD, es, ", &
                & "and section lift slope ( ", pf%WriteOther, " )"
            write(6, '(4x, a, a, a)') "O - Edit output file name ( ", &
                & trim(pf%FileName), " )"

            ! Operating Conditions
            write(6, *)
            write(6, '(2x, a)') "Operating Conditions:"
            write(6, '(4x, a, f7.4, a)') "G - Edit root aerodynamic angle of attack (", &
                & pf%AngleOfAttack * 180.0d0 / pi, " deg )"
            write(6, '(4x, a, f7.4, a)') "B - Edit lift coefficient (", &
                & pf%LiftCoefficient, " )"
            write(6, '(4x, a, f7.4, a)') "L - Edit amount of linear twist (", &
                & pf%Omega * 180.0d0 / pi, " deg )"

            ! Main Execution commands
            write(6, *)
            write(6, '(2x, a)') "R - Run Simulation"
            write(6, '(2x, a)') "Q - Quit"

            write(6, *)
            write(6, '(a)') "Your selection: "

            input = GetCharInput()
            write(6, *)

            cont = MainPageResponse(pf, input)
        end do
    end subroutine MainPage

    function MainPageResponse(pf, input) result(cont)
        type(Planform), intent(inout) :: pf
        character, intent(in) :: input
        logical :: cont

        cont = .true.

        ! Process input command
        ! Wing parameters
        if (input == 'W') then
            call ToggleWingType(pf)
        else if (input == 'N') then
            call EditNodes(pf)
        else if (input == 'A') then
            call EditAspectRatio(pf)
        else if (input == 'T' .and. pf%WingType == Tapered) then
            call EditTaperRatio(pf)
        else if (input == 'S') then
            call EditLiftSlope(pf)

        ! Output options
        else if (input == 'C') then
            pf%WriteCMatrix = .not. pf%WriteCMatrix
        else if (input == 'I') then
            pf%WriteCInverse = .not. pf%WriteCInverse
        else if (input == 'F') then
            pf%WriteFourier = .not. pf%WriteFourier
        else if (input == 'K') then
            pf%WriteOther = .not. pf%WriteOther
        else if (input == 'O') then
            call EditFileName(pf)

        ! Operating Conditions
        else if (input == 'G') then
            call EditAngleOfAttack(pf)
        else if (input == 'B') then
            call EditLiftCoefficient(pf)
        else if (input == 'L') then
            call EditOmega(pf)

        ! Main Execution Commands
        else if (input == 'R') then
            call RunSimulation(pf)
            call system('pause')
        else if (input == 'Q') then
            cont = .false.
        end if
    end function MainPageResponse

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
        call RunSimulation(pf)
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
        call RunSimulation(pf)
    end subroutine EditLiftCoefficient

    subroutine EditOmega(pf)
        type(Planform), intent(inout) :: pf

        write(6, *)
        write(6, '(a)') "Enter amount of linear twist (degrees):"

        pf%Omega = GetRealInput() * pi / 180.0d0
    end subroutine EditOmega

    character function GetCharInput() result(input)
        read(5, '(a)') input
        if(iachar(input) >= iachar('a') .AND. iachar(input) <= iachar('z')) then
            input = char(iachar(input) - 32)
        end if
    end function GetCharInput

    integer function GetIntInput() result(input)
        read(5, *) input
    end function GetIntInput

    real*8 function GetRealInput() result(input)
        read(5, *) input
    end function GetRealInput

end module iofunctions
