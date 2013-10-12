program hello
    use class_Planform
    use iofunctions

    implicit none

    logical :: cont = .true.

    type(Planform) :: planform1
    planform1%WingType = Elliptic

    !Begin execution loop
    do while (cont)
        call MainPage(planform1)
        cont = .false.
    end do

end program
