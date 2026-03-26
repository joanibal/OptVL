cc#ifdef USE_CPOML
      REAL(kind=avl_real) XYN1
      REAL(kind=avl_real) XYN2
      REAL(kind=avl_real) ZLON1
      REAL(kind=avl_real) ZLON2
      REAL(kind=avl_real) ZUPN1
      REAL(kind=avl_real) ZUPN2
      REAL(kind=avl_real) XYZSURF
      REAL(kind=avl_real) CPSURF
      COMMON /VRTX_S/
     & XYN1(2,NVMAX), XYN2(2,NVMAX),  ! left/right aft-node of element
     & ZLON1(NVMAX), ZLON2(NVMAX),    ! left/right lower z-coord of aft-node
     & ZUPN1(NVMAX), ZUPN2(NVMAX),     ! left/right upper z-coord of aft-node
     & XYZSURF(3, NVMAX, NFMAX),
     & CPSURF(NVMAX,NFMAX)
cc#endif