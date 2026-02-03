! Declarations for heap memory storage for AVL
! Harold Youngren 6/21/2023

module avl_heap_inc

! Heap array include file for large AVL AIC arrays

  INTEGER, PARAMETER :: NVX=5000   ! must match nvmax in ADIMEN.INC)
  INTEGER :: NAIC 
  

! All non-constant variables are declared as threadprivate for OpenMP
 ! warning! hardcoding precision!
  REAL(8), DIMENSION(:,:), ALLOCATABLE :: AICN
!!  !$omp threadprivate(AICN)

  REAL(8), DIMENSION(:,:), ALLOCATABLE :: AICN_LU
!!  !$omp threadprivate(AICN_LU)

  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: WC_GAM
!!  !$omp threadprivate(WC_GAM)

  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: WV_GAM
!!  !$omp threadprivate(WV_GAM)


end module avl_heap_inc
