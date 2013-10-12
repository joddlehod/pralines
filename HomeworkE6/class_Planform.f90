module class_Planform
    implicit none

    public :: Planform

    real*8, parameter :: pi = acos(-1.0d0)

    ! Supported wing types
    enum, bind(C)
        enumerator :: Tapered = 1, Elliptic = 2
    end enum

    type Planform
        integer :: WingType = Tapered ! Wing type
        integer :: NNodes = 7 ! Total number of nodes
        real*8 :: AspectRatio = 8.0d0 ! Aspect ratio
        real*8 :: TaperRatio = 1.0d0 ! Taper ratio (tapered wing only)
        real*8 :: LiftSlope = 2.0d0 * pi ! Section lift slope

        character*80 :: FileName = ".\HomeworkE6.out" ! Name of output file
        logical :: WriteCMatrix = .true.  ! Write C Matrix to output file?
        logical :: WriteCInverse = .true.  ! Write Inv(C) Matrix to output file?
        logical :: WriteFourier = .true.  ! Write Fourier coefficients?
        logical :: WriteOther = .true.  ! Write KL, KD, eps, and lift slope?

    end type Planform

    contains
        function GetWingType(this) result(name)
            type(Planform) :: this
            character*8 :: name

            if (this%WingType .eq. Tapered) then
                name = "Tapered"
            else if (this%WingType .eq. Elliptic) then
                name = "Elliptic"
            else
                name = "Unknown"
            end if
        end function GetWingType

end module class_Planform


