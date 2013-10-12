module iofunctions
    use class_Planform
    implicit none

contains
    subroutine MainPage(pf)
        type(Planform), intent(inout) :: pf

        character :: input
        logical :: cont = .true.

        do while (cont)
            ! Clear the screan
            call system('cls')

            ! Display options to user
            write(6, '(a)') "Select from the following options:"

            ! Input parameters
            write(6, '(2x, a)') "Input Parameters:"
            write(6, '(4x, a, a, a)') "W - Edit wing type ( ", &
                & trim(GetWingType(pf)), " )"
            write(6, '(4x, a, i3, a)') "N - Edit number of nodes per semispan (", &
                & pf%NodesPerSemispan, " )"
            write(6, '(4x, a, f7.4, a)') "A - Edit aspect ratio (", &
                & pf%AspectRatio, " )"
            if (pf%WingType == Tapered) then
                write(6, '(4x, a, f7.4, a)') "T - Edit taper ratio (", &
                    & pf%TaperRatio, " )"
            end if
            write(6, '(4x, a, f10.7, a)') "S - Edit section lift slope (", &
                & pf%LiftSlope, " )"
            write(6, '(4x, a, a, a)') "O - Edit output file name ( ", &
                & trim(pf%FileName), " )"

            ! Output options
            write(6, *)
            write(6, '(2x, a)') "Output Options:"
            write(6, '(4x, a, l1, a)') "C - Toggle output of C matrix ( ", &
                & pf%WriteCMatrix, " )"
            write(6, '(4x, a, l1, a)') "I - Toggle output of C-inverse matrix ( ", &
                & pf%WriteCInverse, " )"
            write(6, '(4x, a, l1, a)') "F - Toggle output of Fourier Coefficients ( ", &
                & pf%WriteFourier, " )"
            write(6, '(4x, a, a, l1, a)') "K - Toggle output of KL, KD, and ", &
                & " section lift slope ( ", pf%WriteOther, " )"

            ! Main Execution commands
            write(6, *)
            write(6, '(2x, a)') "R - Run Simulation"
            write(6, '(2x, a)') "Q - Quit"

            read(5, '(a)') input
            if(iachar(input) >= iachar('a') .AND. iachar(input) <= iachar('z')) then
                input = char(iachar(input) - 32)
            end if

            cont = MainPageResponse(pf, input)
        end do
    end subroutine MainPage

    function MainPageResponse(pf, input) result(cont)
        type(Planform), intent(inout) :: pf
        character, intent(in) :: input
        logical :: cont

        cont = .true.

        ! Process input command
        ! Input parameters
        if (input == 'W') then
!            call EditWingType(pf)
        else if (input == 'N') then
!            call EditNodes(pf)
        else if (input == 'A') then
!            call EditAspectRatio(pf)
        else if (input == 'T' .and. pf%WingType == Tapered) then
!            call EditTaperRatio(pf)
        else if (input == 'S') then
!            call EditLiftSlope(pf)
        else if (input == 'O') then
!            call EditFileName(pf)

            ! Output options
        else if (input == 'C') then
            pf%WriteCMatrix = .not. pf%WriteCMatrix
        else if (input == 'I') then
            pf%WriteCInverse = .not. pf%WriteCInverse
        else if (input == 'F') then
            pf%WriteFourier = .not. pf%WriteFourier
        else if (input == 'K') then
            pf%WriteOther = .not. pf%WriteOther

            ! Main Execution Commands
        else if (input == 'R') then
!            call RunSimulation(pf)
        else if (input == 'Q') then
            cont = .false.
        end if
    end function MainPageResponse

end module iofunctions
