module LiftingLineSolver_Test
    use Utilities
    use class_Planform
    use LiftingLineSetters
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

        integer :: badline
        character*80 :: fname_work = "Problem1p34b_work.txt"
        character*80 :: fname_results = "Problem1p34b_results.txt"

		call InitPlanform(pf)

        call SetWingType(pf, Elliptic)
        call SetNNodes(pf, 100)
        call SetSectionLiftSlope(pf, 2.0d0 * pi)
        call SetAspectRatio(pf, 5.56d0)
        call SetAileronRoot(pf, 0.253d0)
        call SetAileronTip(pf, 0.438d0)
        call SetFlapFractionRoot(pf, 0.28d0)
        call SetFlapFractionTip(pf, 0.25d0)
        call SetHingeEfficiency(pf, 0.85d0)
        call SetFileName(pf, fname_work)

        call ComputeCMatrixAndCoefficients(pf)
        call OutputPlanform(pf)

        call SetWashout(pf, 2.0d0)
        call SetAileronDeflection(pf, 5.0d0)
        call SetSteadyRollingRate(pf)
        call SetLiftCoefficient(pf, 0.5d0)

        call ComputeFlightConditions(pf)
        call OutputFlightConditions(pf)

        call SetWashout(pf, 2.0d0)
        call SetAileronDeflection(pf, 5.0d0)
        call SetRollingRate(pf, -0.02d0)
        call SetLiftCoefficient(pf, 0.5d0)

        call ComputeFlightConditions(pf)
        call OutputFlightConditions(pf)

        call SetWashout(pf, 2.0d0)
        call SetAileronDeflection(pf, 5.0d0)
        call SetRollingRate(pf, 0.0d0)
        call SetLiftCoefficient(pf, 0.5d0)

        call ComputeFlightConditions(pf)
        call OutputFlightConditions(pf)

        badline = CompareFiles(fname_work, fname_results)
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
