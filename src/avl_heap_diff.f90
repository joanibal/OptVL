!=============================================================================80
! Allocates AIC arrays that may be too big for COMMONS
!=============================================================================80
subroutine avlheap_diff_init()

  use avl_heap_diff_inc

! Allocate AIC variable storage
 
  write(*,*) 'allocating ', 8*NVX**2/1024**2
  if (.not. allocated(AICN_DIFF)) then
    allocate(AICN_DIFF(NVX,NVX))
    allocate(AICN_LU_DIFF(NVX,NVX))
    allocate(WC_GAM_DIFF(3,NVX,NVX))
    allocate(WV_GAM_DIFF(3,NVX,NVX))
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

