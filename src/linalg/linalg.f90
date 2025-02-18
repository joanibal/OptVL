!!!
! This module aims to privide an efficent yet portable replacement of a LAPACK library
! 
! The following subroutines are reimplemented:
!       1. LUDCMP
!       CALL SGETRF(N,NSIZ,A,NSIZ,INDX,INFO)  
!       2. BAKSUB
!       CALL SGETRS('N',N,MRHS,A,NSIZ,INDX,B,NSIZ,INFO)
!       3. BACKSUB with transpose
!       CALL SGETRS('T',N,MRHS,A,NSIZ,INDX,B,NSIZ,INFO)
! 
!!!
module linalg
    implicit none
    integer, parameter :: real_kind = 8
    
    character(7) :: fmt_real = '(f10.4)'
    character(5) :: fmt_int = '(I10)'
    
contains


    subroutine ludcmp(nsiz, n, a, indx, work)
        ! implementation of the court decomposition algorithm
        ! references 
        ! lapack documentation
        ! https://www.netlib.org/lapack/explore-html-3.6.1/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html
        
        ! openblas algorithm overview
        ! dgetf2
        
        
        ! https://github.com/openmathlib/openblas/blob/33bb4b98a4725efef26517d6fbec7a81260ee531/lapack-netlib/src/dgetf2.f
        ! numerical recipes page 2.3.1
        
        ! explanation
        ! https://www.andreinc.net/2021/01/20/writing-your-own-linear-algebra-matrix-library-in-c#the-lup-algorithm-as-an-example
        implicit none
        integer nsiz, n
        real(kind=real_kind) a(nsiz, nsiz), work(nsiz)
        integer indx(nsiz)
        real(kind=real_kind) aamax, sum, dum
        integer i, j, k, imax

        do i = 1, n
            aamax = 0.0
            ! find the largest element in column j
            do j = 1, n
                aamax = max(abs(a(i, j)), aamax)
            end do
            work(i) = 1.0 / aamax
        end do

        do j = 1, n
            do i = 1, j - 1
                sum = a(i, j)
                do k = 1, i - 1
                    sum = sum - a(i, k) * a(k, j)
                end do
                a(i, j) = sum
            end do

            aamax = 0.0
            do i = j, n
                sum = a(i, j)
                do k = 1, j - 1
                    sum = sum - a(i, k) * a(k, j)
                end do
                a(i, j) = sum

                dum = work(i) * abs(sum)
                if (dum .ge. aamax) then
                    imax = i
                    aamax = dum
                end if
            end do

            if (j .ne. imax) then
                do k = 1, n
                    dum = a(imax, k)
                    a(imax, k) = a(j, k)
                    a(j, k) = dum
                end do
                work(imax) = work(j)
            end if

            indx(j) = imax
            if (j .ne. n) then
                dum = 1.0 / a(j, j)
                do i = j + 1, n
                    a(i, j) = a(i, j) * dum
                end do
            end if
        end do
    end subroutine ludcmp

    subroutine baksub(nsiz, n, a, indx, b)
        implicit none
        integer nsiz, n
        real(kind=real_kind) a(nsiz, nsiz), b(nsiz)
        integer indx(nsiz)
        integer ii, ll, i, j
        real(kind=real_kind) sum

        
        ! Solve A * X = B.
                 
        ! Apply row interchanges to the right hand sides.
        ! Solve L*X = B, overwriting B with X.
        ! forward sub
        ii = 0
        do i = 1, n
            ll = indx(i)
            sum = b(ll)
            b(ll) = b(i)
            if (ii .ne. 0) then
                do j = ii, i - 1
                    sum = sum - a(i, j) * b(j)
                end do
            else if (sum .ne. 0.0) then
                ii = i
            end if
            b(i) = sum
        end do
        ! Solve U*X = B, overwriting B with X.
        do i = n, 1, -1
            sum = b(i)
            if (i .lt. n) then
                do j = i + 1, n
                    sum = sum - a(i, j) * b(j)
                end do
            end if
            b(i) = sum / a(i, i)
        end do
    end subroutine baksub
    
    subroutine baksubtrans(nsiz, n, a, indx, b)
        implicit none
        integer nsiz, n
        real(kind=real_kind) a(nsiz, nsiz), b(nsiz), tmp
        integer indx(nsiz)
        integer ii, ll, i, j
        real(kind=real_kind) sum
        
        ! Solve A**T * X = B.
        ! the A matrix passed has L and U factors stored in place
        ! notes diagonal of L is implied to be all ones
        
        ! Solve U**T *X = B, overwriting B with X.
        ! forward substitution but indexing gets weird because we won'T
        ! actually transpose the matrix 
        do i = 1, n
            sum = b(i)
            if (i > 1) then ! skip the first
                do j = 1, i-1
                    sum = sum - a(j, i) * b(j)
                end do
            end if
            b(i) = sum /a(i,i)
        end do
        ! Solve L**T *X = B, overwriting B with X.
        do i = n, 1, -1
            sum = b(i)
            if (i .lt. n) then
                do j = i + 1, n
                    sum = sum - a(j, i) * b(j)
                end do
            end if
            b(i) = sum 
        end do
        
        ! Apply row interchanges to the solution vectors.
        do i = N,1, -1
            tmp = b(indx(i))
            b(indx(i)) = b(i)
            b(i) = tmp
        enddo
        
        
    end subroutine baksubtrans

! utilities
    subroutine print_matrix(A)
        real(kind=real_kind), intent(in) :: A(:, :)
        integer :: i, j
        integer :: m, n
        
        ! get the dimensions of the matrix
        m = size(A, 1)
        n = size(A, 2)
        
        
        do j = 1, m
            do i = 1, n
                write(*, fmt_real, advance='no') a(j, i)
                if (i < m) write(*, '(a)', advance='no') ' '  ! add a space between elements
            end do
            write(*,*) ''  ! new line after each row is printed
        end do
    end subroutine print_matrix

    subroutine print_vector(b)
        real(kind=real_kind), intent(in) :: b(:)
        integer :: n, i
        
        ! get the size of the vector
        n = size(b)
        
        do i = 1, n
            write(*,fmt_real) b(i)
        end do
    end subroutine print_vector

    subroutine print_ivector(b)
        integer, intent(in) :: b(:)
        integer :: n, i
        
        ! get the size of the vector
        n = size(b)
        
        do i = 1, n
            write(*,fmt_int) b(i)
        end do
    end subroutine print_ivector

end module linalg
