      subroutine write_tecplot(file_name, add_time, solution_time)
C
      INCLUDE 'AVL.INC'
      INCLUDE 'AVL_oml.INC'
C
      LOGICAL LRANGEALL
      INTEGER LU, nChords, nStrips, nPts
      CHARACTER*(*) file_name
      logical :: add_time
      real(kind=8) :: solution_time
C
C...  check that all surfaces use full range of airfoil definitions
      LRANGEALL = .TRUE.
      DO ISURF = 1, NSURF
        LRANGEALL = LRANGEALL .AND. LRANGE(ISURF)
      ENDDO
      IF (.NOT.LRANGEALL) THEN
      WRITE(*,*) 'ERROR in write_tecplot: implemented only for ',
     &    'surfaces defined using full range of input airfoils'
        WRITE(*,*) '  returning without writing OML output'
        RETURN
      ENDIF
      
c      write out the header
c      we have to use "block packing since the data is all cell cenered"

c      Open the file for writing. Use a unit number (out_unit) that is not in use.
      LU = 13
      OPEN(UNIT=LU, FILE=file_name, STATUS='REPLACE', ACTION='WRITE')
      
      WRITE(LU, '(A,A,A)') 'TITLE = "', trim(TITLE),'"'
      WRITE(LU, '(A)') 'VARIABLES = "X", "Y", "Z" "CP"'
      WRITE(LU, '(A,f4.2,A)') 'DATASETAUXDATA AVL="',VERSION,'"'
      
      
      DO ISURF = 1, NSURF
            nChords = NK(ISURF)
            nStrips = NJ(ISURF)
            
            nPts = (nStrips+1)*(nChords*2+1)
            ! nPts = (nStrips+1)*(nChords+1)
            nCCPts = (nStrips)*(nChords*2)
        WRITE(LU, '(A,A,A,A,2(A,I0))') 'ZONE T="',
     &   trim(STITLE(ISURF)),'", DATAPACKING=BLOCK, ',
     &   'VARLOCATION=(4=CELLCENTERED), ',
     &   'I=',(2*nChords+1), ', J=', (nStrips+1)
        if (add_time) then 
            WRITE(LU, '(A,f20.16)') 'SOLUTIONTIME=',solution_time
        endif 
            
        DO idim = 1,3
            DO iPt = 1,nPts
                WRITE(LU, '(F15.8)') XYZSURF(idim, iPt, ISURF)
            ENDDO 
        ENDDO
        
      DO iPt = 1,nCCPts
            WRITE(LU, '(F15.8)') CPSURF(iPt, ISURF)
      ENDDO 

        
      ENDDO  
      CLOSE(LU)
      
      end subroutine write_tecplot
      
      
