program hello
    use class_Planform
    use liftinglinesolver
    use iofunctions

    implicit none

    type(Planform) :: pf
    call RunSimulation(pf)

    !Begin execution
    call MainPage(pf)

end program
