c     This file is to trick tapenade into treating the allocated module variables as normal common block variables
      REAL(kind=avl_real) AMACH
      REAL(kind=avl_real) AICN
      REAL(kind=avl_real) AICN_LU
      REAL(kind=avl_real) WC_GAM
      REAL(kind=avl_real) WV_GAM
      COMMON /AFAKE/
     & AICN(NVMAX,NVMAX),     ! normalwash AIC matrix (and VL system matrix)
     & AICN_LU(NVMAX,NVMAX),     ! LU facotrization of AICN
     & WC_GAM(3,NVMAX,NVMAX), ! c.p. velocity/Gamma influence matrix
     & WV_GAM(3,NVMAX,NVMAX) ! h.v. velocity/Gamma influence matrix
C