!=============================================================================80
! Allocates AIC arrays that may be too big for COMMONS
!=============================================================================80
subroutine avlheap_diff_init(n)

  use avl_heap_diff_inc
  
  integer :: n

! Allocate AIC variable storage
 
  if (.not. allocated(AICN_DIFF)) then
    allocate(AICN_DIFF(n,n))
    allocate(AICN_LU_DIFF(n,n))
    allocate(WC_GAM_DIFF(3,n,n))
    allocate(WV_GAM_DIFF(3,n,n))
  endif
end subroutine avlheap_diff_init

!=============================================================================80
! Free AIC array storage
!=============================================================================80
subroutine avlheap_diff_clean()

  use avl_heap_diff_inc

! Deallocate heap storage for AIC's 

  if (allocated(AICN_DIFF)) then
    deallocate(AICN_DIFF)
    deallocate(AICN_LU_DIFF)
    deallocate(WC_GAM_DIFF)
    deallocate(WV_GAM_DIFF)
  endif

end subroutine avlheap_diff_clean

