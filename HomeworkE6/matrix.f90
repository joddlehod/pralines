module matrix
    implicit none

contains
    subroutine matinv(n, mat, mat_inv)
        integer, intent(in) :: n
        real*8, intent(in) :: mat(n,n)
        real*8, intent(out) :: mat_inv(n,n)

        integer :: i, j, k
        real*8 :: mat_aug(n, 2*n)
        real*8 :: pivot, factor

        mat_aug = reshape((/ ((0.0d0, i=1, n), j=1, 2*n) /), shape(mat_aug))

        do i = 1, n
            mat_aug(i, i + n) = 1.0d0
            do j = 1, n
                mat_aug(i, j) = mat(i, j)
            end do
        end do

        do i = 1, n
            do j = 1, n
                if (dabs(mat_aug(i, j)) > 0.0d0) then
                    pivot = mat_aug(i, j)
                    exit
                end if
            end do

            do j = 1, 2*n
                mat_aug(i, j) = mat_aug(i, j) / pivot
            end do
            pivot = 1.0d0

            do j = 1, n
                if (j /= i) then
                    factor = mat_aug(j, i) / pivot
                    do k = 1, 2*n
                        mat_aug(j, k) = mat_aug(j, k) - factor * mat_aug(i, k)
                    end do
                end if
            end do
        end do

        mat_inv = reshape((/ ((mat_aug(i, j), i=1, n), j=n+1, 2*n) /), shape(mat_inv))

    end subroutine matinv

    subroutine printmat(m, n, mat)
        integer, intent(in) :: m
        integer, intent(in) :: n
        real*8, intent(in) :: mat(m, n)

        integer :: i, j
        character*80 :: format_string

        write(format_string, '(a, i3, a)') "(", n, "(F9.5, 2x))"
        do i = 1, m
            write(6, format_string) (mat(i, j), j=1, n)
        end do
    end subroutine printmat

end module matrix
