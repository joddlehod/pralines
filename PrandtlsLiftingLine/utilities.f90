module Utilities
    implicit none

    real*8, parameter :: pi = acos(-1.0d0)
    real*8, parameter :: zero = 1.0d-10

contains
    integer function Compare(a, b, tol) result(eq)
    ! Comparison function
    ! Inputs:
    !   a = First argument to compare
    !   b = Second argument to compare
    !   tol = Relative tolerance for comparison
    ! Return Value (eq):
    !   -1 = a < (b - tol)
    !    0 = a == b (within tolerance)
    !    1 = a > (b + tol)
        real*8, intent(in) :: a, b, tol

        if (abs(a) < tol .and. abs(b) < tol) then
            eq = 0
        else if (abs(a - b) / max(abs(a), abs(b)) < tol) then
            eq = 0
        else if (a < b) then
            eq = -1
        else
            eq = 1
        end if
    end function Compare

    real*8 function Residual(oldVal, newVal) result(res)
        real*8, intent(in) :: oldVal
        real*8, intent(in) :: newVal

        res = dabs(oldVal - newVal) / max(dabs(oldVal), dabs(newVal))
    end function Residual

    integer function CompareFiles(a, b) result(badline)
        character*80, intent(in) :: a, b  ! Filenames of files to compare

        integer :: i, ios1, ios2
        character*5000 :: results_line, work_line

        open(unit=11, file=a)
        open(unit=12, file=b)

        badline = 0
        ios1 = 0
        i = 0
        do while (badline == 0 .and. ios1 == 0)
            i = i + 1

            results_line(1:5000) = " "
            read(11, '(A)', iostat=ios1, end=99) results_line

            work_line(1:5000) = " "
            read(12, '(A)', iostat=ios2, end=99) work_line

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
99  continue
    end function CompareFiles

    character*2 function GetCharacterInput(def) result(inp)
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
    end function GetCharacterInput

    character*80 function GetStringInput(def) result(inp)
        character*80, intent(in) :: def  ! Default value if invalid input

        integer :: i, ios
        character :: a

        read(5, '(a)', iostat=ios) inp
        if (len(trim(inp)) < 1 .or. ios /= 0) then
            inp = def
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
                ios = ParseFormula(trim(inp_str), inp)
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

    integer function ParseFormula(inp_str, num) result(estat)
        character*80, intent(in) :: inp_str
        real*8, intent(out) :: num

        integer :: i, j, n_oper, last_ind, strlen, ios
        character*40 :: operators, temp_num
        real*8, Dimension(41) :: numbers

        estat = 0
        strlen = len(trim(inp_str))
        n_oper = 0
        last_ind = 0
        do i = 2, strlen
            if (inp_str(i:i) == '*' .or. inp_str(i:i) == '/') then
                n_oper = n_oper + 1
                operators(n_oper:n_oper) = inp_str(i:i)
                temp_num = "                                        "
                temp_num(1:i-last_ind-1) = inp_str(last_ind+1:i-1)
                if ((temp_num(1:1) == 'P' .or. temp_num(1:1) == 'p') .and. &
                    & (temp_num(2:2) == 'I' .or. temp_num(2:2) == 'i')) then
                    numbers(n_oper) = pi
                else
                    read(temp_num, *, iostat=ios) numbers(n_oper)
                    if (ios /= 0) then
                        estat = 1
                    end if
                end if
                last_ind = i
            end if
        end do

        temp_num = "                                        "
        temp_num(1:strlen-last_ind) = inp_str(last_ind+1:strlen)
        if ((temp_num(1:1) == 'P' .or. temp_num(1:1) == 'p') .and. &
            & (temp_num(2:2) == 'I' .or. temp_num(2:2) == 'i')) then
            numbers(n_oper + 1) = pi
        else
            read(temp_num, *, iostat = ios) numbers(n_oper + 1)
            if (ios /= 0) then
                estat = 1
            end if
        end if

        num = numbers(1)
        do i = 1, n_oper
            if (operators(i:i) == '*') then
                num = num * numbers(i + 1)
            else if (operators(i:i) == '/') then
                num = num / numbers(i + 1)
            else
                estat = 2
            end if
        end do
    end function ParseFormula

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

                if (order <= 0) then
                    width = width - order + 1  ! e.g. -0.012345 - additional for leading 0
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

end module Utilities
