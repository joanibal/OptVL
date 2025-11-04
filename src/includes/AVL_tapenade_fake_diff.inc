c     This file is to trick tapenade into treating the allocated module variables as normal common block variables
      REAL(kind=avl_real) AICN_diff
      REAL(kind=avl_real) AICN_LU_diff
      REAL(kind=avl_real) WC_GAM_diff
      REAL(kind=avl_real) WV_GAM_diff
      COMMON /AFAKE_R/
     & AICN_diff(NVMAX,NVMAX),     ! normalwash AIC matrix _diff(and VL system matrix)
     & AICN_LU_diff(NVMAX,NVMAX),     ! LU facotrization of AICN
     & WC_GAM_diff(3,NVMAX,NVMAX), ! c.p. velocity/Gamma influence matrix
     & WV_GAM_diff(3,NVMAX,NVMAX) ! h.v. velocity/Gamma influence matrix
C