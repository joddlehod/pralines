module LiftingLineSolver_Test
    use class_Planform
    use LiftingLineSolver
    use LiftingLineOutput

    implicit none

contains
    subroutine TestLiftingLineSolver()
        integer :: nerror

        nerror = 0
        nerror = nerror + TestProblem1p34b()

        write(6, *)

        if (nerror == 0) then
            write(6, '(a)') "SUCCESS - All tests passed."
        else
            write(6, '(a, i3, a)') "FAIL - ", nerror, " tests failed."
        end if

        call system('pause')
    end subroutine TestLiftingLineSolver

    integer function TestProblem1p34b() result(fail)
        type(Planform) :: pf

        integer :: i
        integer :: badline, ios1, ios2
        character*5000 :: results_line, work_line

        pf%NNodes = 199
        pf%LiftSlope = 2.0d0 * pi
        pf%AspectRatio = 5.56d0
        pf%AileronRoot = 0.253d0
        pf%AileronTip = 0.438d0
        pf%FlapFractionRoot = 0.28d0
        pf%FlapFractionTip = 0.25d0
        pf%ParallelHingeLine = .false.
        pf%HingeEfficiency = 0.85d0
        pf%FileName = "Problem1p34b_work.txt"

        call ComputeCMatrixAndCoefficients(pf)
        call OutputPlanform(pf)

        pf%Omega = 2.0d0 * pi / 180.0d0
        pf%AileronDeflection = 5.0d0 * pi / 180.0d0
        pf%RollingRate = -0.04361952421640
        pf%DesiredLiftCoefficient = 0.5d0
        pf%LiftCoefficient = 0.5d0
        pf%SpecifyAlpha = .false.

        call ComputeFlightConditions(pf)
        call OutputFlightConditions(pf)

        pf%Omega = 2.0d0 * pi / 180.0d0
        pf%AileronDeflection = 5.0d0 * pi / 180.0d0
        pf%RollingRate = -0.02d0
        pf%DesiredLiftCoefficient = 0.5d0
        pf%LiftCoefficient = 0.5d0
        pf%SpecifyAlpha = .false.

        call ComputeFlightConditions(pf)
        call OutputFlightConditions(pf)

        pf%Omega = 2.0d0 * pi / 180.0d0
        pf%AileronDeflection = 5.0d0 * pi / 180.0d0
        pf%RollingRate = 0.0d0
        pf%DesiredLiftCoefficient = 0.5d0
        pf%LiftCoefficient = 0.5d0
        pf%SpecifyAlpha = .false.

        call ComputeFlightConditions(pf)
        call OutputFlightConditions(pf)

        open(unit=11, file="Problem1p34b_results.txt")
        open(unit=12, file=pf%FileName)

        badline = 0
        ios1 = 0
        i = 0
        do while (badline == 0 .and. ios1 == 0)
            i = i + 1

            results_line(1:5000) = " "
            read(11, '(A)', iostat=ios1, end=99) results_line

            work_line(1:5000) = " "
            read(12, '(A)', iostat=ios2, end=99) work_line
99  continue

            if (ios1 == 0) then
                if (work_line /= results_line) then
                    badline = i
                end if
            else if (len(trim(work_line)) /= 0) then
                badline = i
            end if
        end do

        close(unit=11)
        close(unit=12)

        if (badline /= 0) then
            write(6, '(a)') "Test Failed - Problem1p34b"
            write(6, '(2x, a, i3, a)') "Comparison of line ", badline, " failed!"
            fail = 1
        else
            write(6, '(a)') "Test Successful - Problem1p34b"
            fail = 0
        end if
    end function TestProblem1p34b

end module LiftingLineSolver_Test
