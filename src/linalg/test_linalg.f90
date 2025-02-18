program main
    use linalg
    implicit none

    
    integer, parameter :: nsiz = 2
    integer :: n
    real :: a(nsiz, nsiz)
    integer :: indx(nsiz)
    real :: work(nsiz)
    real :: b(nsiz)
    
    

    
    ! call test_baksub(4)
    call test_baksubtrans(10)

    ! n = 2

    ! ! initialize matrix a
    ! a(1, 1) = 4.0
    ! a(1, 2) = 3.0
    ! a(2, 1) = 6.0
    ! a(2, 2) = 3.0

    ! ! print the original matrix
    ! print *, 'original matrix a:'
    ! call print_matrix(nsiz, n, a)

    ! ! perform lu decomposition
    ! call ludcmp(nsiz, n, a, indx, work)

    ! ! print the lu decomposed matrix
    ! print *, 'lu decomposed matrix a:'
    ! call print_matrix(nsiz, n, a)

    ! ! initialize right-hand side vector b
    ! b(1) = 10.0
    ! b(2) = 15.0

    ! ! perform back-substitution
    ! call baksub(nsiz, n, a, indx, b)

    ! ! print the solution vector
    ! print *, 'solution vector b:'
    ! call printvector(nsiz, n, b)
    
    ! if (b(1) .eq. 2.50 .and. b(2) .eq. 0.00) then 
    !     print *, 'test passed'
    ! else
    !     print *, 'test failed!!!'
    ! end if 

contains

    subroutine test_baksub(N)
        integer, intent(in) :: N ! the size of the test matrix
        real(kind=real_kind) :: A(N, N)
        real(kind=real_kind) :: X(N), X_lp(N)
        
        ! working 
        integer :: rand_size
        integer :: idxs_pivots(N), info, MRHS
        integer, allocatable, dimension(:) :: seed
        
        
        
        ! generate a random matrix A
        ! call random_seed(size=size)
        
        ! look at this block of code just to set the random seed, smh fortran.
        call random_seed(size = rand_size)
        allocate(seed(rand_size))
        call random_seed(get=seed)
        seed = 0.0
        call random_seed(put=seed)
        
        ! fill A with random numbers
        call random_number(A)
        call random_number(X)
        X_lp = X
        
     
        ! call print_matrix(A)
        
        ! first use lapack to get the factorization of A
        CALL DGETRF(N,N,A,N,idxs_pivots,info)
        
        ! now use the factorization to solve Ax = b
        ! use lapack to get the reference value
        MRHS = 1
        CALL DGETRS('N',N,MRHS,A,N,idxs_pivots,X_lp,N,info)
        
        ! print the solution vector
        call print_vector(X_lp) 
        write(*,*) ''
        
        
        ! now use the baksub function to get the solution
        call baksub(N, N, A, idxs_pivots, X)
        
        call print_vector(X) 
        
    end subroutine test_baksub
    
    subroutine test_baksubtrans(N)
        integer, intent(in) :: N ! the size of the test matrix
        real(kind=real_kind) :: A(N, N)
        real(kind=real_kind) :: X(N), X_lp(N)
        
        ! working 
        integer :: rand_size
        integer :: idxs_pivots(N), info, MRHS
        integer, allocatable, dimension(:) :: seed
        
        
        
        
        ! look at this block of code just to set the random seed, smh fortran.
        call random_seed(size = rand_size)
        allocate(seed(rand_size))
        call random_seed(get=seed)
        seed = 0.0
        call random_seed(put=seed)
        
        ! fill A with random numbers
        call random_number(A)
        call random_number(X)
        ! a(1, 1) = 1.0
        ! a(1, 2) = 2.0
        ! a(2, 1) = 2.0
        ! a(2, 2) = 3.0
        ! X(1) = 1.0
        ! X(2) = 1.0
        X_lp = X
        
        write(*,*) 'input matrix'
        call print_matrix(A)
        print*, ''


        
        
     
        ! call print_matrix(A)
        
        ! first use lapack to get the factorization of A
        CALL DGETRF(N,N,A,N,idxs_pivots,info)
        print*, 'factorized'
        call print_matrix(A)
        call print_ivector(idxs_pivots)
        print*, ''
        
        ! now use the factorization to solve Ax = b
        ! use lapack to get the reference value
        MRHS = 1
        CALL DGETRS('T',N,MRHS,A,N,idxs_pivots,X_lp,N,info)
        
        ! print the solution vector
        write(*,*) 'right answer'
        call print_vector(X_lp) 
        write(*,*) ''
        write(*,*) 'my answer'
        
        ! now use the baksub function to get the solution
        call baksubtrans(N, N, A, idxs_pivots, X)
        
        call print_vector(X) 
        
    end subroutine test_baksubtrans
end program main

