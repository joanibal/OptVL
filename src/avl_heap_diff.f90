!=============================================================================80
! Allocates AIC arrays that may be too big for COMMONS
!=============================================================================80
subroutine avlheap_diff_init(n)

  use avl_heap_diff_inc
  
  integer :: n

! Allocate AIC variable storage
 
! Use heap_diff_allocated flag instead of allocated() to avoid false positives
! when Fortran allocatable descriptors contain garbage with -fno-init-global-zero

  if (.not. heap_diff_allocated) then
    allocate(AICN_DIFF(n,n))
    allocate(AICN_LU_DIFF(n,n))
    allocate(WC_GAM_DIFF(3,n,n))
    allocate(WV_GAM_DIFF(3,n,n))
    heap_diff_allocated = .TRUE.
  endif
end subroutine avlheap_diff_init

!=============================================================================80
! Free AIC array storage
!=============================================================================80
subroutine avlheap_diff_clean()

  use avl_heap_diff_inc

! Deallocate heap storage for AIC's 

  if (heap_diff_allocated) then
    deallocate(AICN_DIFF)
    deallocate(AICN_LU_DIFF)
    deallocate(WC_GAM_DIFF)
    deallocate(WV_GAM_DIFF)
    heap_diff_allocated = .FALSE.
  endif

end subroutine avlheap_diff_clean

