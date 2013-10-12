program hello
    use class_Planform
    use iofunctions

    implicit none

    type(Planform) :: planform1
    planform1%WingType = Elliptic

    !Begin execution loop
    call MainPage(planform1)

end program
