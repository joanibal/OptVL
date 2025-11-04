!=============================================================================80
! Allocates AIC arrays that may be too big for COMMONS
!=============================================================================80
subroutine avlheap_init(n)

  use avl_heap_inc
  
  integer :: n

! Allocate AIC variable storage

  if (.not. allocated(AICN)) then
    allocate(AICN(n,n))
    allocate(AICN_LU(n,n))
    allocate(WC_GAM(3,n,n))
    allocate(WV_GAM(3,n,n))
  endif
  
  NAIC = n
end subroutine avlheap_init

!=============================================================================80
! Free AIC array storage
!=============================================================================80
subroutine avlheap_clean()

  use avl_heap_inc

! Deallocate heap storage for AIC's 

  if (allocated(AICN)) then
    deallocate(AICN)
    deallocate(AICN_LU)
    deallocate(WC_GAM)
    deallocate(WV_GAM)
  endif
  
  NAIC = -1

end subroutine avlheap_clean

