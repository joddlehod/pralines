module matrix
    implicit none

contains
    subroutine matinv_gauss(n, mat, mat_inv)
        integer, intent(in) :: n
        real*8, intent(in) :: mat(n,n)
        real*8, intent(out) :: mat_inv(n,n)

        real*8 :: b(n, n), c, d, temp(n)
        integer :: i, j, k, m, imax(1), ipvt(n)

        b = mat
        ipvt = (/ (i, i=1, n) /)

        do k = 1, n
            imax = maxloc(abs(b(k:n, k)))
            m = k - 1 + imax(1)

            if (m /= k) then
                ipvt( (/m, k/) ) = ipvt( (/k, m/) )
                b( (/m, k/), :) = b( (/k, m/), :)
            end if

            d = 1.0d0 / b(k, k)

            temp = b(:, k)
            do j = 1, n
                c = b(k, j) * d
                b(:, j) = b(:, j) - temp * c
                b(k, j) = c
            end do
            b(:, k) = temp * (-d)
            b(k, k) = d
        end do

        mat_inv(:, ipvt) = b
    end subroutine matinv_gauss


    subroutine printmat(u, m, n, mat)
        integer, intent(in) :: u
        integer, intent(in) :: m
        integer, intent(in) :: n
        real*8, intent(in) :: mat(m, n)

        integer :: i, j
        character*80 :: format_string

        if (u == 6) then
            write(format_string, '(a, i10, a)') "(", n, "(F9.5, 2x))"
        else
            write(format_string, '(a, i10, a)') "(", n, "(ES22.15, 2x))"
        end if

        do i = 1, m
            write(u, format_string) (mat(i, j), j=1, n)
        end do
    end subroutine printmat

end module matrix
