      SUBROUTINE solve_adjoint():
        IMPLICIT NONE
        INCLUDE 'AVL_d.inc'
        REAL df_dr(NVOR)
        
        CALL BAKSUB(NVOR,NVOR,AICN_LU,IAPIV,df_dr)
        
        
        
      END ! ADJOINT
        
        
        