!=============================================================================80
! Allocates AIC arrays that may be too big for COMMONS
!=============================================================================80
subroutine avlheap_init(n)

  use avl_heap_inc
  
  integer :: n

! Allocate AIC variable storage
! Use heap_allocated flag instead of allocated() to avoid false positives
! when Fortran allocatable descriptors contain garbage with -fno-init-global-zero

  if (.not. heap_allocated) then
    allocate(AICN(n,n))
    allocate(AICN_LU(n,n))
    allocate(WC_GAM(3,n,n))
    allocate(WV_GAM(3,n,n))
    heap_allocated = .TRUE.
  endif

  NAIC = n
end subroutine avlheap_init

!=============================================================================80
! Free AIC array storage
!=============================================================================80
subroutine avlheap_clean()

  use avl_heap_inc

! Deallocate heap storage for AIC's 

  if (heap_allocated) then
    deallocate(AICN)
    deallocate(AICN_LU)
    deallocate(WC_GAM)
    deallocate(WV_GAM)
    heap_allocated = .FALSE.
  endif

  NAIC = -1

end subroutine avlheap_clean

