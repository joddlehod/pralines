module LiftingLineSetters
    use class_Planform
    implicit none

contains
    subroutine InitPlanform(pf)
        type(Planform), intent(inout) :: pf

        call SetParallelHingeLine(pf)
    end subroutine InitPlanform

    ! Planform Parameters
    subroutine SetWingType(pf, wingType)
        type(Planform), intent(inout) :: pf
        integer, intent(in) :: wingType

        if (pf%WingType /= wingType) then
            pf%WingType = wingType
            call DeallocateArrays(pf)

            if (pf%WingType == Combination) then
                call SetCombinationWingCoefficients(pf)
            else if (pf%ParallelHingeLine) then
                call SetParallelHingeLine(pf)
            end if
            if (pf%WingType == Elliptic) then
                pf%WashoutDistribution = Linear
            end if
        end if
    end subroutine SetWingType

    subroutine SetWashoutDistribution(pf, washoutDist)
        type(Planform), intent(inout) :: pf
        integer, intent(in) :: washoutDist

        if (pf%WingType /= Elliptic .and. pf%WashoutDistribution /= washoutDist) then
            pf%WashoutDistribution = washoutDist
            call DeallocateArrays(pf)
        end if
    end subroutine SetWashoutDistribution

    subroutine SetTransitionPoint(pf, tp)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: tp

        if (Compare(pf%TransitionPoint, tp, zero) /= 0) then
            pf%TransitionPoint = tp
            call DeallocateArrays(pf)
            call SetCombinationWingCoefficients(pf)
        end if
    end subroutine SetTransitionPoint

    subroutine SetTransitionChord(pf, tc)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: tc

        if (Compare(pf%TransitionChord, tc, zero) /= 0) then
            pf%TransitionChord = tc
            call DeallocateArrays(pf)
            call SetCombinationWingCoefficients(pf)
        end if
    end subroutine SetTransitionChord

    subroutine SetCombinationWingCoefficients(pf)
        type(Planform), intent(inout) :: pf

        real*8 :: u, asin_u

        pf%C1 = pf%TransitionPoint
        pf%C2 = (1.0d0 - pf%TransitionChord) / pf%C1
        pf%C4 = (pf%C1 - 2.0d0 * pf%C1**2 * pf%C2 + 0.25d0 * pf%C2) / &
            & (pf%C1 * pf%C2 - pf%C2 + 1.0d0)
        u = (pf%C1 - pf%C4) / (0.5d0 - pf%C4)
        asin_u = asin(u)
        pf%C3 = (1.0d0 - pf%C1 * pf%C2) / sqrt(1.0d0 - u**2)
        pf%C5 = 1.0d0 / (pf%AspectRatio * (2.0d0 * pf%C1 - pf%C1**2 * pf%C2 + &
            & 0.5d0 * pf%C3 * (0.5d0 - pf%C4) * (pi - 2.0d0 * asin_u - &
            & sin(2.0d0 * asin_u))))
        !pf%C5 = 1.0d0 / (2.0d0 * pf%AspectRatio * &
        !    & (pf%C1 - 0.5d0 * pf%C1**2 * pf%C2 + &
        !    & 0.25d0 * pf%C3 * (0.5d0 - pf%C4) * &
        !    & (0.5d0 * pi - sin(2.0d0 * asin_u) - asin_u)))

        if (pf%ParallelHingeLine) then
            call SetParallelHingeLine(pf)
        end if
    end subroutine SetCombinationWingCoefficients

    logical function AreCombinationWingCoefficientsValid(pf) result(isValid)
        type(Planform), intent(in) :: pf

        real*8 :: u

        u = (pf%C1 - pf%C4) / (0.5d0 - pf%C4)

        if (Compare(pf%C1, 0.0d0, zero) == 0) then
            isValid = .false.
        else if (Compare(pf%C1 * pf%C2 - pf%C2 + 1.0d0, 0.0d0, zero) == 0) then
            isValid = .false.
        else if (Compare(pf%C4, 0.5d0, zero) == 0) then
            isValid = .false.
        else if (Compare(u, 0.0d0, zero) /= 1) then
            isValid = .false.
        else if (Compare(u, 1.0d0, zero) /= -1) then
            isValid = .false.
        else
            isValid = .true.
        end if
    end function AreCombinationWingCoefficientsValid

    subroutine SetNNodes(pf, npss)
        type(Planform), intent(inout) :: pf
        integer, intent(in) :: npss

        integer :: nnodes

        nnodes = npss * 2 - 1
        if (pf%NNodes /= nnodes) then
            pf%NNodes = nnodes
            call DeallocateArrays(pf)
        end if
    end subroutine SetNNodes

    subroutine SetAspectRatio(pf, ra)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: ra

        if (Compare(pf%AspectRatio, ra, zero) /= 0) then
            pf%AspectRatio = ra
            call DeallocateArrays(pf)
        end if
    end subroutine SetAspectRatio

    subroutine SetTaperRatio(pf, rt)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: rt

        if (Compare(pf%TaperRatio, rt, zero) /= 0) then
            pf%TaperRatio = rt
            call DeallocateArrays(pf)
        end if
    end subroutine SetTaperRatio

    subroutine SetSectionLiftSlope(pf, cla_sec)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: cla_sec

        if (Compare(pf%SectionLiftSlope, cla_sec, zero) /= 0) then
            pf%SectionLiftSlope = cla_sec
            call DeallocateArrays(pf)
        end if
    end subroutine SetSectionLiftSlope

    subroutine SetAileronRoot(pf, ar)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: ar

        if (Compare(pf%AileronRoot, ar, zero) /= 0) then
            pf%AileronRoot = ar
            call DeallocateArrays(pf)
        end if
    end subroutine SetAileronRoot

    subroutine SetAileronTip(pf, at)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: at

        if (Compare(pf%AileronTip, at, zero) /= 0) then
            pf%AileronTip = at
            call DeallocateArrays(pf)
        end if
    end subroutine SetAileronTip

    subroutine SetParallelHingeLine(pf)
        type(Planform), intent(inout) :: pf

        real*8 :: cfc_root_par

        pf%ParallelHingeLine = .true.
        cfc_root_par = ParallelRootFlapFraction(pf)
        if (Compare(pf%FlapFractionRoot, cfc_root_par, zero) /= 0) then
            pf%FlapFractionRoot = cfc_root_par
            call DeallocateArrays(pf)
        end if
    end subroutine SetParallelHingeLine

    subroutine SetFlapFractionRoot(pf, cfc_root)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: cfc_root

        pf%ParallelHingeLine = .false.
        pf%DesiredFlapFractionRoot = cfc_root
        if (Compare(pf%FlapFractionRoot, cfc_root, zero) /= 0) then
            pf%FlapFractionRoot = cfc_root
            call DeallocateArrays(pf)
        end if
    end subroutine SetFlapFractionRoot

    subroutine SetFlapFractionTip(pf, cfc_tip)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: cfc_tip

        if (Compare(pf%FlapFractionTip, cfc_tip, zero) /= 0) then
            pf%FlapFractionTip = cfc_tip
            call DeallocateArrays(pf)

            if (pf%ParallelHingeLine) then
                call SetParallelHingeLine(pf)
            end if
        end if
    end subroutine SetFlapFractionTip

    subroutine SetHingeEfficiency(pf, eff_hinge)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: eff_hinge

        if (Compare(pf%HingeEfficiency, eff_hinge, zero) /= 0) then
            pf%HingeEfficiency = eff_hinge
            call DeallocateArrays(pf)
        end if
    end subroutine SetHingeEfficiency

    subroutine SetDeflectionEfficiency(pf, eff_def)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: eff_def

        if (Compare(pf%DeflectionEfficiency, eff_def, zero) /= 0) then
            pf%DeflectionEfficiency = eff_def
            call DeallocateArrays(pf)
        end if
    end subroutine SetDeflectionEfficiency

    subroutine ToggleOutputMatricies(pf)
        type(Planform), intent(inout) :: pf

        pf%OutputMatrices = .not. pf%OutputMatrices
    end subroutine ToggleOutputMatricies

    subroutine SetFileName(pf, filename)
        type(Planform), intent(inout) :: pf
        character*80, intent(in) :: filename

        pf%FileName = trim(filename)
    end subroutine SetFileName


    ! Operating Conditions
    subroutine SetAngleOfAttack(pf, alpha)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: alpha

        pf%SpecifyAlpha = .true.
        pf%DesiredAngleOfAttack = alpha * pi / 180.0d0
        pf%AngleOfAttack = pf%DesiredAngleOfAttack
        pf%LiftCoefficient = CL1(pf%CLa, pf%AngleOfAttack, pf%EW, pf%Washout)
    end subroutine SetAngleOfAttack

    subroutine SetLiftCoefficient(pf, cl)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: cl

        pf%SpecifyAlpha = .false.
        pf%DesiredLiftCoefficient = cl
        pf%LiftCoefficient = pf%DesiredLiftCoefficient
        pf%AngleOfAttack = RootAlpha(pf%CLa, pf%LiftCoefficient, pf%EW, pf%Washout)
    end subroutine SetLiftCoefficient

    subroutine SetAileronDeflection(pf, da)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: da

        pf%AileronDeflection = da * pi / 180.0d0
    end subroutine SetAileronDeflection

    subroutine SetWashout(pf, washout)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: washout

        pf%DesiredWashout = washout * pi / 180.0d0
        pf%Washout = pf%DesiredWashout
        pf%UseOptimumWashout = .false.
    end subroutine SetWashout

    subroutine SetOptimumWashout(pf)
        type(Planform), intent(inout) :: pf

        pf%OptimumWashout1 = (pf%KDL * pf%LiftCoefficient) / &
            & (2.0d0 * pf%KDW * pf%CLa)
        if (pf%WashoutDistribution == Optimum) then
            pf%OptimumWashout2 = (4.0d0 * pf%CL1) / (pi * pf%AspectRatio &
                & * pf%SectionLiftSlope * c_over_b(pf, pi / 2.0d0))
        end if
        pf%Washout = pf%OptimumWashout1
        pf%UseOptimumWashout = .true.
    end subroutine SetOptimumWashout

    subroutine SetRollingRate(pf, rollingrate)
        type(Planform), intent(inout) :: pf
        real*8, intent(in) :: rollingrate

        pf%DesiredRollingRate = rollingrate
        pf%RollingRate = rollingrate
        pf%UseSteadyRollingRate = .false.
    end subroutine

    subroutine SetSteadyRollingRate(pf)
        type(Planform), intent(inout) :: pf

        real*8 :: steady_pbar

        pf%RollingRate = SteadyRollingRate(pf)
        pf%UseSteadyRollingRate = .true.
    end subroutine SetSteadyRollingRate

    real*8 function ParallelRootFlapFraction(pf) result(cfc_root_par)
        type(Planform), intent(in) :: pf

        real*8 :: cb_root, cfc_root
        real*8 :: cb_tip, cfc_tip

        cb_tip = c_over_b_zb(pf, pf%AileronTip)
        cfc_tip = pf%FlapFractionTip
        cb_root = c_over_b_zb(pf, pf%AileronRoot)
        cfc_root_par = 0.75d0 - cb_tip / cb_root * (0.75d0 - cfc_tip)
    end function ParallelRootFlapFraction

    real*8 function RootAlpha(cla, cl, ew, omega) result(alpha)
        real*8, intent(in) :: cla
        real*8, intent(in) :: cl
        real*8, intent(in) :: ew
        real*8, intent(in) :: omega

        alpha = cl / cla + ew * omega
    end function RootAlpha

    real*8 function SteadyRollingRate(pf) result(pbar_steady)
        type(Planform), intent(inout) :: pf

        ! Calculate steady dimensionless rolling rate
        pbar_steady = -pf%CRM_da / pf%CRM_pbar * pf%AileronDeflection
    end function SteadyRollingRate

    real*8 function CL1(cla, alpha, ew, w) result(cl)
        real*8, intent(in) :: cla
        real*8, intent(in) :: alpha
        real*8, intent(in) :: ew
        real*8, intent(in) :: w

        cl = cla * (alpha - ew * w)  ! Eq. 1.8.24
    end function CL1

    real*8 function CL2(ra, bigA1) result(cl)
        real*8, intent(in) :: ra
        real*8, intent(in) :: bigA1

        cl = pi * ra * bigA1  ! Eq. 1.8.5
    end function CL2

    real*8 function CDi1(pf) result(cdi)
        type(Planform), intent(in) :: pf

        cdi = (pf%CL1**2 * (1.0d0 + pf%KD) - pf%KDL * pf%CL1 * pf%CLa * pf%Washout + &
            & pf%KDW * (pf%CLa * pf%Washout)**2) / (pi * pf%AspectRatio)
    end function CDi1

    real*8 function CDi2(pf) result(cdi)
        type(Planform), intent(in) :: pf

        integer :: i

        cdi = 0.0d0
        do i = 1, pf%NNodes
            cdi = cdi + real(i, 8) * pf%BigA(i)**2
        end do
        cdi = cdi * pi * pf%AspectRatio
    end function CDi2

    real*8 function CDi3(pf) result(cdi)
        type(Planform), intent(in) :: pf

        integer :: i
        real*8 :: ri

        cdi = 0.0d0
        do i = 1, pf%NNodes
            ri = real(i, 8)
            cdi = cdi + ri * pf%BigA(i)**2
        end do
        cdi = (cdi - 0.5d0 * pf%RollingRate * pf%BigA(2)) * pi * pf%AspectRatio
    end function CDi3

end module LiftingLineSetters
