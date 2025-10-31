C***********************************************************************
C    Module:  aoper.f
C 
C    Copyright (C) 2002 Mark Drela, Harold Youngren
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************

      SUBROUTINE OPER
C---------------------------------------
C     Main driver routine for AVL
C---------------------------------------
      INCLUDE 'AVL.INC'
C       INCLUDE 'AVLPLT.INC'
      LOGICAL ERROR, LCERR, LCERRI, LWRIT, LMATCH, LMATCH_TEST
      LOGICAL USEMRF
C
      CHARACTER*1 ANS
      CHARACTER*2 OPT
      CHARACTER*4 COMAND, COMAND_TEST, ITEMC
      CHARACTER*256 FNOUT, FNDER, FNNEW, FNVB
      CHARACTER*120 LINE, COMARG, COMARG_TEST, CRUN, PROMPT, RTNEW
      CHARACTER*50 SATYPE, ROTTYPE
C
      LOGICAL LOWRIT
C
      REAL    RINPUT(20), RINP(20)
      INTEGER IINPUT(20), IINP(20), IRUN_TEST
C
      IF(.NOT.LGEO) THEN
       WRITE(*,*)
       WRITE(*,*) '* Configuration not defined'
       RETURN
      ENDIF
C
C
      FNVB = ' '
C
      CALL GETSA(LNASA_SA,SATYPE,DIR)
C
      IF(LSA_RATES) THEN
        ROTTYPE = 'Rates,moments about Stability axes'
       ELSE
        ROTTYPE = 'Rates,moments about Body axes'
      ENDIF
C
 1000 FORMAT (A)
C
C       LPLOT = .FALSE.
      LWRIT = .FALSE.
      LSOL  = .FALSE.
C
      USEMRF = .FALSE.
C
      FNOUT = ' '
C
C=================================================================
C---- start of user interaction loop
 800  CONTINUE
C
      LCERR = .FALSE.
C
C
      CALL CFRAC(IRUN,NRUN,CRUN,NPR)
C
      if(lverbose)then
            WRITE(*,1050) CRUN(1:NPR), RTITLE(IRUN)
            CALL CONLST(IRUN)
      end if
C       WRITE(*,1052)
C
 1050 FORMAT(
     &  /' Operation of run case ',A,':  ', A
     &  /' ==========================================================')
C
 1052 FORMAT(
     &  /'  C1  set level or banked  horizontal flight constraints'
     &  /'  C2  set steady pitch rate (looping) flight constraints'
     &  /'  M odify parameters                                    '
     & //' "#" select  run case          L ist defined run cases   '
     &  /'  +  add new run case          S ave run cases to file   '
     &  /'  -  delete  run case          F etch run cases from file'
     &  /'  N ame current run case       W rite forces to file     '
     & //' eX ecute run case             I nitialize variables     '
     & //'  G eometry plot               T refftz Plane plot       '
     & //'  ST  stability derivatives    FT  total   forces        '
     &  /'  SB  body-axis derivatives    FN  surface forces        '
     &  /'  RE  reference quantities     FS  strip   forces',
     &                                         ' (FSB bodyaxes)'
     &  /'  DE  design changes           FE  element forces        '
     &  /'  O ptions                     FB  body forces           '
     &  /'                               HM  hinge moments         '
     &  /'  OB offbody flow survey       VM  strip shear,moment    '
     &  /'  MRF  machine-readable format CPOM OML surface pressures')
C 
C   A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
C   x x x x   x x   x     x x x x x   x x x   x x x x

 810  CONTINUE
C       PROMPT = ' .OPER (case ' // CRUN(1:NPR) // ')^'
C       CALL ASKC(PROMPT,COMAND,COMARG)
C C
C------------------------------------------------------
C C       IF    (COMAND.EQ.'    ') THEN
C C        IF(LPLOT) CALL PLEND
C C C        LPLOT = .FALSE.
C C C     CALL CLRZOOM
C C        RETURN
C C
C       IF(COMAND.EQ.'?   ') THEN
C        GO TO 800
C C
C       ENDIF
C C
C C ------------------------------------------------------
C C---- check for run case commands
C       IF(COMAND.EQ.'+   ') THEN
C C----- add new case after current one
C C
C        IF(NRUN.EQ.NRMAX) THEN
C         WRITE(*,*)
C         WRITE(*,*) '* Run case array limit NRMAX reached'
C        ELSE
C         NRUN = NRUN + 1
C C
C         DO JR = NRUN, IRUN+1, -1
C           CALL RCOPY(JR,JR-1)
C         ENDDO
C         WRITE(*,*) 'Initializing new run case from current one'
C C        
C         IRUN = IRUN + 1
C        ENDIF
C C
C        GO TO 800
C C
C       ELSEIF(COMAND.EQ.'-   ') THEN
C C----- delete current case
C C
C        IF(NRUN.LE.1) THEN
C         WRITE(*,*)
C         WRITE(*,*) '* Cannot delete one remaining run case'
C        ELSE
C         DO JR = IRUN, NRUN-1
C           CALL RCOPY(JR,JR+1)
C         ENDDO
C         NRUN  = NRUN - 1
C         IRUN  = MAX( 1 , MIN( IRUN  , NRUN ) )
C         IRUNE = MAX( 1 , MIN( IRUNE , NRUN ) )
C        ENDIF
C C
C        GO TO 800
C C
C       ENDIF
C C
C------------------------------------------------------
C---- first check if point index was input
C       IF(INDEX('0123456789',COMAND(1:1)) .NE. 0) THEN
C         READ(COMAND,*,ERR=8) IRUN
C         IF(IRUN.LT.1 .OR. IRUN.GT.NRUN) THEN
C          IRUN = MAX( 1 , MIN(NRUN,IRUN) )
C          WRITE(*,*) '* Run case index was limited to available range'
C         ENDIF
C         GO TO 800
C C
C  8      WRITE(*,*) '* Invalid command: ', COMAND
C         GO TO 800
C       ENDIF
C C
C C------------------------------------------------------
C C---- extract command line numeric arguments
C       DO I=1, 20
C         IINPUT(I) = 0
C         RINPUT(I) = 0.0
C       ENDDO
C       NINPUT = 20

C       CALL GETINT(COMARG,IINPUT,NINPUT,ERROR)

C       NINPUT = 20

C       CALL GETFLT(COMARG,RINPUT,NINPUT,ERROR)
C C
C  14   CONTINUE
C C------------------------------------------------------
C C---- check for parameter toggle/set command
      
C       WRITE(*,*) '------------------COM 1-----------------------'
C       WRITE(*,*) COMAND, 'COMAND'
C       WRITE(*,*) COMARG, 'COMARG'
C       WRITE(*,*) LMATCH,  'LMATCH'
C       WRITE(*,*) IRUN, 'IRUN'
C       WRITE(*,*) 
C       WRITE(*,*) ALFA, 'ALFA'
C       WRITE(*,*) '-----------------------------------------------'





C      CALL CONSET(COMAND,COMARG)! ,LMATCH,IRUN)
C C
C       IF(LMATCH) THEN
C C----- match found... go back to OPER menu
C        GO TO 800
C       ENDIF
C C
C C------------------------------------------------------
C---- check for trim set command

C      CALL TRMSET(COMAND,COMARG,LMATCH,IRUN)
C
C       IF(LMATCH) THEN
C C----- match found... go back to OPER menu
C        GO TO 800
C       ENDIF

C------------------------------------------------------
C---- pick up here to try decoding for remaining commands
C
C       IF(COMAND .EQ. 'X   ') THEN
C------ execute calculation



        IF(LCERR) THEN
         WRITE(*,*) '** Flow solution is not possible.'
         WRITE(*,*) '** Cannot impose a constraint more than once.'
C          GO TO 800
        ENDIF
C
        IRUN0 = IRUN
C
        XYZREF(1) = PARVAL(IPXCG,IRUN)
        XYZREF(2) = PARVAL(IPYCG,IRUN)
        XYZREF(3) = PARVAL(IPZCG,IRUN)
        CDREF     = PARVAL(IPCD0,IRUN)
C
        INFO = 1
        CALL EXEC(NITMAX,INFO,IRUN)
C        IF(.NOT.LSOL) GO TO 810
C C
        if(lverbose) then
            IF(LPTOT)   CALL OUTTOT(6)
            IF(LPSURF)  CALL OUTSURF(6)
            IF(LPSTRP)  CALL OUTSTRP(6)
            IF(LPELE)   CALL OUTELE(6)
            IF(LPHINGE) CALL OUTHINGE(6)
        end if


C
C       ENDIF


C C------------------------------------------------------
C       ELSEIF(COMAND .EQ. 'XX  ') THEN
C C------ execute calculation for all run cases
C         DO 24 IR = 1, NRUN
C C-------- check for well-posedness
C           LCERRI = .FALSE.
C           DO IV = 1, NVTOT
C             IC = ICON(IV,IR)
C             DO JV = 1, NVTOT
C               IF(IV.NE.JV .AND. ICON(IV,IR).EQ.ICON(JV,IR)) THEN
C                LCERRI = .TRUE.
C               ENDIF
C             ENDDO
C           ENDDO
C           IF(LCERRI) THEN
C            WRITE(*,*) '** Run case', IR,' ...'
C            WRITE(*,*) '** Flow solution is not possible.'
C            WRITE(*,*) '** Cannot impose a constraint more than once.'
C            GO TO 24          
C           ENDIF
C C
C           INFO = 1
C           CALL EXEC(NITMAX,INFO,IR)
C           IF(.NOT.LSOL) GO TO 24
C C
C           IF(LPTOT)  CALL OUTTOT(6)
C           IF(LPSURF) CALL OUTSURF(6)
C           IF(LPSTRP) CALL OUTSTRP(6)
C           IF(LPELE)  CALL OUTELE(6)
C           IF(LPHINGE) CALL OUTHINGE(6)
C  24     CONTINUE
C C
C C------------------------------------------------------
C       ELSEIF(COMAND .EQ. 'M   ') THEN
C         CALL PARMOD(IRUN)
C C
C C------------------------------------------------------
C       ELSEIF(COMAND .EQ. 'N   ') THEN
C C------ change name of run case
C         IF(COMARG.NE.' ') THEN
C          RTITLE(IRUN) = COMARG
C         ELSE
C          WRITE(*,830) RTITLE(IRUN)
C  830     FORMAT(/' Enter run case name:  ', A)
C          READ(*,1000) RTNEW
C          IF(RTNEW.NE.' ') RTITLE(IRUN) = RTNEW
C         ENDIF
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'DE  ') THEN
C C------ design changes
C         IF(NDESIGN.EQ.0) THEN
C          WRITE(*,*) '* No design parameters are declared'
C          GO TO 810
C         ENDIF
C C
C  30     CONTINUE
C         WRITE(*,1036)
C         DO K = 1, NDESIGN
C           WRITE(*,1037) K, GNAME(K), DELDES(K)
C         ENDDO
C  1036   FORMAT(/' ================================================'
C      &         /' Current design parameter changes:' /
C      &         /'    k   Parameter      change')
C CCC              x1234xxx123456789012345612345678901234
C  1037   FORMAT(1X, I4,3X, A, G14.5) 
C C
C         WRITE(*,*)
C  35     WRITE(*,*) 'Enter  k, design changes (<return> if done) ...'
C  37     READ (*,1000) LINE
C         CALL STRIP(LINE,NLIN)
C         IF(LINE(1:1).EQ.'?') GO TO 30
C         IF(LINE(1:1).EQ.' ') GO TO 800
C C
C         NINP = 40
C         CALL GETFLT(LINE,RINP,NINP,ERROR)
C         IF(ERROR) THEN
C          WRITE(*,*) '* Bad input'
C          GO TO 35
C         ENDIF
C C
C         DO I = 1, NINP, 2
C           N = INT(RINP(I))
C           IF(N.LT.1 .OR. N.GT.NDESIGN) THEN
C            WRITE(*,*) 'Index k out of bounds. Input ignored.'
C           ENDIF
C           DELDES(N) = RINP(I+1)
C           LSOL = .FALSE.
C           GO TO 37
C         ENDDO
C         GO TO 30
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'I   ') THEN
C C------ clear operating parameters
C         ALFA = 0.
C         BETA = 0.
C         WROT(1) = 0.
C         WROT(2) = 0.
C         WROT(3) = 0.
C C
C         DO N = 1, NCONTROL
C           DELCON(N) = 0.
C         ENDDO
C C
C         DO N = 1, NDESIGN
C           DELDES(N) = 0.
C         ENDDO
C C
C         LSOL = .FALSE.
C C
C C------------------------------------------------------
C C       ELSE IF(COMAND.EQ.'G   ') THEN
C C------ plot geometry
C C       CALL PLOTVL(AZIMOB, ELEVOB, TILTOB, ROBINV)
C C         LPLOT = .TRUE.
C C
C C------------------------------------------------------
C C       ELSE IF(COMAND.EQ.'T   ') THEN
C C C------ plot spanloadings in Trefftz plane
C C C         CALL PLOTTP
C C         LPLOT = .TRUE.
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'FT  ') THEN
C C------ print total forces
C         IF(LSOL) THEN
C          CALL GETFILE(LU,COMARG)
C C
C          IF(LU.LE.-1) THEN
C           WRITE(*,*) '* Filename error *'
C          ELSEIF(LU.EQ.0) THEN
C           WRITE(*,*) '* Data not written'
C          ELSE
C           IF(USEMRF) THEN
C            CALL MRFTOT(LU, 'TOT')
C           ELSE
C            CALL OUTTOT(LU)
C           ENDIF
C           IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
C          ENDIF
C C
C         ELSE
C          WRITE(*,*) '* Execute flow calculation first!'
C         ENDIF
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'FN  ') THEN
C C------ print surface forces
C         IF(LSOL) THEN
C          CALL GETFILE(LU,COMARG)
C C
C          IF(LU.LE.-1) THEN
C           WRITE(*,*) '* Filename error *'
C          ELSEIF(LU.EQ.0) THEN
C           WRITE(*,*) '* Data not written'
C          ELSE
C           IF(USEMRF) THEN
C            CALL MRFSURF(LU)
C          ELSE
C           CALL OUTSURF(LU)
C          ENDIF
C           IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
C          ENDIF
C C
C         ELSE
C          WRITE(*,*) '* Execute flow calculation first!'
C         ENDIF
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'FS  ') THEN
C C------ print strip forces
C         IF(LSOL) THEN
C          CALL GETFILE(LU,COMARG)
C C
C          IF(LU.LE.-1) THEN
C           WRITE(*,*) '* Filename error *'
C          ELSEIF(LU.EQ.0) THEN
C           WRITE(*,*) '* Data not written'
C          ELSE
C           IF(USEMRF) THEN
C            CALL MRFSTRP(LU)
C          ELSE
C           CALL OUTSTRP(LU)
C           ENDIF
C           IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
C          ENDIF
C C
C         ELSE
C          WRITE(*,*) '* Execute flow calculation first!'
C         ENDIF
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'FSB ') THEN
C C------ print strip forces in body axes
C         IF(LSOL) THEN
C          CALL GETFILE(LU,COMARG)
C C
C          IF(LU.LE.-1) THEN
C           WRITE(*,*) '* Filename error *'
C          ELSEIF(LU.EQ.0) THEN
C           WRITE(*,*) '* Data not written'
C          ELSE
C           CALL OUTSTRPB(LU)
C           IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
C          ENDIF
C C
C         ELSE
C          WRITE(*,*) '* Execute flow calculation first!'
C         ENDIF
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'FE  ') THEN
C C------ print vortex element forces
C         IF(LSOL) THEN
C          CALL GETFILE(LU,COMARG)
C C
C          IF(LU.LE.-1) THEN
C           WRITE(*,*) '* Filename error *'
C          ELSEIF(LU.EQ.0) THEN
C           WRITE(*,*) '* Data not written'
C          ELSE
C           IF(USEMRF) THEN
C            CALL MRFELE(LU)
C          ELSE
C           CALL OUTELE(LU)
C          ENDIF
C           IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
C          ENDIF
C C
C         ELSE
C          WRITE(*,*) '* Execute flow calculation first!'
C         ENDIF
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'FB  ') THEN
C C------ print body forces
C         IF(LSOL) THEN
C          CALL GETFILE(LU,COMARG)
C C
C          IF(LU.LE.-1) THEN
C           WRITE(*,*) '* Filename error *'
C          ELSEIF(LU.EQ.0) THEN
C           WRITE(*,*) '* Data not written'
C         ELSE
C          IF(USEMRF) THEN
C           CALL MRFBODY(LU)
C          ELSE
C           CALL OUTBODY(LU)
C          ENDIF
C           IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
C          ENDIF
C C
C         ELSE
C          WRITE(*,*) '* Execute flow calculation first!'
C         ENDIF
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'HM  ') THEN
C C------ print hinge moments
C         IF(LSOL) THEN
C          CALL GETFILE(LU,COMARG)
C C
C          IF(LU.LE.-1) THEN
C           WRITE(*,*) '* Filename error *'
C          ELSEIF(LU.EQ.0) THEN
C           WRITE(*,*) '* Data not written'
C          ELSE
C           IF(USEMRF) THEN
C            CALL MRFHINGE(LU)
C          ELSE
C           CALL OUTHINGE(LU)
C           ENDIF
C           IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
C          ENDIF
C C
C         ELSE
C          WRITE(*,*) '* Execute flow calculation first!'
C         ENDIF
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'VM  ') THEN
C C------ calculate and print shear and bending on surfaces
C         IF(LSOL) THEN
C          CALL GETFILE(LU,COMARG)
C C
C          IF(LU.LE.-1) THEN
C           WRITE(*,*) '* Filename error *'
C          ELSEIF(LU.EQ.0) THEN
C           WRITE(*,*) '* Data not written'
C          ELSE
C           IF(USEMRF) THEN
C            CALL MRFVM(LU)
C           ELSE
C            CALL OUTVM(LU)
C           ENDIF
C           IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
C          ENDIF
C C
C         ELSE
C          WRITE(*,*) '* Execute flow calculation first!'
C         ENDIF
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'CN  ') THEN
C C------ print a spanloading file
C         IF(LSOL) THEN
C          CALL GETFILE(LU,COMARG)
C C
C          IF(LU.LE.-1) THEN
C           WRITE(*,*) '* Filename error *'
C          ELSEIF(LU.EQ.0) THEN
C           WRITE(*,*) '* Data not written'
C          ELSE
C           IF(USEMRF) THEN
C            CALL MRFCNC(LU)
C           ELSE
C            CALL OUTCNC(LU)
C           ENDIF
C           IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
C          ENDIF
C C
C         ELSE
C          WRITE(*,*) '* Execute flow calculation first!'
C         ENDIF
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'L   ') THEN
C C------ list run cases
C         LU = 6
C         CALL RUNSAV(LU)
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'OB  ') THEN
C C------ off-body flow survey
C         CALL OBOPER
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'S   ') THEN
C C------ save run case file
C         CALL STRIP(FRNDEF,NFR)
C         WRITE(*,2040) FRNDEF(1:NFR)
C  2040   FORMAT(' Enter run case filename: ', A)
C         READ (*,1000) FNNEW
C C
C         IF(FNNEW.NE.' ') FRNDEF = FNNEW
C         OPEN(LURUN,FILE=FRNDEF,STATUS='OLD',ERR=42)
C C
C         IF(LOWRIT(FRNDEF)) THEN
C          REWIND(LURUN)
C          GO TO 45
C         ELSE
C          WRITE(*,*) 'Run cases not saved'
C          GO TO 810
C         ENDIF
C C
C  42     OPEN(LURUN,FILE=FRNDEF,STATUS='NEW')
C C
C  45     CALL RUNSAV(LURUN)
C         CLOSE(LURUN)
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'F   ') THEN
C C------ fetch run case file
C         CALL STRIP(FRNDEF,NFR)
C         WRITE(*,2050) FRNDEF(1:NFR)
C  2050   FORMAT(' Enter run case filename: ', A)
C         READ (*,1000) FNNEW
C         IF(FNNEW.NE.' ') FRNDEF = FNNEW
C C
C         CALL RUNGET(LURUN,FRNDEF,ERROR)
C         IF(ERROR) THEN
C          GO TO 810
C         ELSE
C          WRITE(*,2055) (IR, RTITLE(IR), IR=1, NRUN)
C  2055    FORMAT(' Run cases read in ...',
C      &          100(/1X,I4,': ',A))
C         ENDIF
C C
C         IRUN = MIN( IRUN, NRUN )
C C
C C------------------------------------------------------
C       IF(COMAND.EQ.'ST  ') THEN
C C------ create stability derivatives
C         IF(LSOL) THEN
C          CALL GETFILE(LU,COMARG)
C C
C          IF(LU.LE.-1) THEN
C           WRITE(*,*) '* Filename error *'
C          ELSEIF(LU.EQ.0) THEN
C           WRITE(*,*) '* Data not written'
C          ELSE
C           CALL DERMATS(LU)


C           IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
C          ENDIF
C C
C         ELSE
C          WRITE(*,*) '* Execute flow calculation first!'
C         ENDIF
C       ENDIF
C C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'SM  ') THEN
C C------ create stability derivatives
C         IF(LSOL) THEN
C          CALL GETFILE(LU,COMARG)
C C
C          IF(LU.LE.-1) THEN
C           WRITE(*,*) '* Filename error *'
C          ELSEIF(LU.EQ.0) THEN
C           WRITE(*,*) '* Data not written'
C          ELSE
C           CALL DERMATM(LU)
C           IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
C          ENDIF
C C
C         ELSE
C          WRITE(*,*) '* Execute flow calculation first!'
C         ENDIF
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'SB  ') THEN
C C------ create stability derivatives
C         IF(LSOL) THEN
C          CALL GETFILE(LU,COMARG)
C C
C          IF(LU.LE.-1) THEN
C           WRITE(*,*) '* Filename error *'
C          ELSEIF(LU.EQ.0) THEN
C           WRITE(*,*) '* Data not written'
C          ELSE
C           CALL DERMATB(LU)
C           IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
C          ENDIF
C C
C         ELSE
C          WRITE(*,*) '* Execute flow calculation first!'
C         ENDIF
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'W   ') THEN
C C------ write force  data to a file
C         CALL STRIP(FNOUT,NFN)
C         WRITE(*,1080) FNOUT(1:NFN)
C  1080   FORMAT('Enter forces output file: ', A)
C         READ (*,1000) FNNEW
C C
C         IF(FNNEW.NE.' ') THEN
C C-------- new filename was entered...
C C-------- if previous file is open, close it
C           IF(LWRIT) CLOSE(LUOUT)
C           FNOUT = FNNEW
C C-------- open new file and write header
C           OPEN(LUOUT,FILE=FNOUT,STATUS='UNKNOWN')
C           LWRIT = .TRUE.
C C
C         ELSE
C C-------- just a <return> was entered...
C           IF(.NOT.LWRIT) THEN
C             WRITE(*,*) 'No action taken.'
C             GO TO 800
C           ENDIF
C C
C         ENDIF
C C
C         IF(USEMRF) THEN
C          IF(LPTOT)   CALL MRFTOT(LUOUT, 'TOT')
C          IF(LPSURF)  CALL MRFSURF(LUOUT)
C          IF(LPSTRP)  CALL MRFSTRP(LUOUT)
C          IF(LPELE)   CALL MRFELE(LUOUT)
C         ELSE
C         IF(LPTOT)   CALL OUTTOT(LUOUT)
C         IF(LPSURF)  CALL OUTSURF(LUOUT)
C         IF(LPSTRP)  CALL OUTSTRP(LUOUT)
C         IF(LPELE)   CALL OUTELE(LUOUT)
C ccc     IF(LPDERIV) CALL DERMAT(LUOUT)
C         ENDIF
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'MRF ') THEN
C C------ all subsequent data file writes are full-precision machine readable format
C         USEMRF = .TRUE.
C         WRITE(*,'(A,A)')
C      &    'MRF: all subsequent output in full-precision ',
C      &    'machine readable format'
C cc#ifdef USE_CPOML
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'CPOM') THEN
C C------ write OML surface pressures to a file
C         CALL CPOML()
C cc#endif
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'RE  ') THEN
C C------ Change reference data 
C  89     XYZREF(1) = PARVAL(IPXCG,IRUN)
C         XYZREF(2) = PARVAL(IPYCG,IRUN)
C         XYZREF(3) = PARVAL(IPZCG,IRUN)
C         CDREF     = PARVAL(IPCD0,IRUN)
C C
C         WRITE(*,2090) SREF,CREF,BREF,XYZREF,CDREF
C  2090   FORMAT(/' ==========================='
C      &         /'  S ref: ', G11.5,
C      &         /'  C ref: ', G11.5,
C      &         /'  B ref: ', G11.5,
C      &         /'  X ref: ', G11.5,
C      &         /'  Y ref: ', G11.5,
C      &         /'  Z ref: ', G11.5,
C      &         /' CD ref: ', G11.5 )
C C
C  90     CALL ASKC(' Select item,value^',ITEMC,COMARG)
C  2100   FORMAT(' Enter new ',A,': ', $)
C C
C         IF(ITEMC.EQ.'    ') THEN
C           GO TO 800
C         ENDIF
C C
C         NINP = 1
C         CALL GETFLT(COMARG,RINP,NINP,ERROR)
C C
C         IF    (INDEX('Cc',ITEMC(1:1)).NE.0 .AND.
C      &         INDEX('Dd',ITEMC(2:2)).NE.0) THEN
C           IF(NINP.EQ.0) THEN
C  97          WRITE(*,2100) 'Added CD0 profile drag'
C            READ (*,*,ERR=97) CDREF
C           ELSE
C            CDREF = RINP(1)
C           ENDIF
C           PARVAL(IPCD0,IRUN) = CDREF
C C
C         ELSEIF(INDEX('Ss',ITEMC(1:1)).NE.0) THEN
C           IF(NINP.EQ.0) THEN
C  91        WRITE(*,2100) 'reference area Sref'
C            READ (*,*,ERR=91) SREF
C           ELSE
C            SREF = RINP(1)
C           ENDIF
C           LSOL = .FALSE.
C C
C         ELSEIF(INDEX('Cc',ITEMC(1:1)).NE.0) THEN
C           IF(NINP.EQ.0) THEN
C  92        WRITE(*,2100) 'reference chord Cref'
C            READ (*,*,ERR=92) CREF
C           ELSE
C            CREF = RINP(1)
C           ENDIF
C           LSOL = .FALSE.
C C
C         ELSEIF(INDEX('Bb',ITEMC(1:1)).NE.0) THEN
C           IF(NINP.EQ.0) THEN
C  93        WRITE(*,2100) 'reference span Bref'
C            READ (*,*,ERR=93) BREF
C           ELSE
C            BREF = RINP(1)
C           ENDIF
C           LSOL = .FALSE.
C C
C         ELSEIF(INDEX('Xx',ITEMC(1:1)).NE.0) THEN
C           IF(NINP.EQ.0) THEN
C  94        WRITE(*,2100) 'reference Xref'
C            READ (*,*,ERR=94) XREF
C           ELSE
C            XREF = RINP(1)
C           ENDIF
C           PARVAL(IPXCG,IRUN) = XREF
C           LSOL = .FALSE.
C           LSRD = .FALSE.
C C
C         ELSEIF(INDEX('Yy',ITEMC(1:1)).NE.0) THEN
C           IF(NINP.EQ.0) THEN
C  95        WRITE(*,2100) 'reference Yref'
C            READ (*,*,ERR=95) YREF
C           ELSE
C            ZREF = RINP(1)
C           ENDIF
C           PARVAL(IPYCG,IRUN) = YREF
C           LSOL = .FALSE.
C           LSRD = .FALSE.
C C
C         ELSEIF(INDEX('Zz',ITEMC(1:1)).NE.0) THEN
C           IF(NINP.EQ.0) THEN
C  96        WRITE(*,2100) 'reference Zref'
C            READ (*,*,ERR=96) ZREF
C           ELSE
C            ZREF = RINP(1)
C           ENDIF
C           PARVAL(IPZCG,IRUN) = ZREF
C           LSOL = .FALSE.
C           LSRD = .FALSE.
C C
C         ELSE
C           WRITE(*,*) 'Item not recognized'
C           GO TO 89
C C
C         ENDIF
C         GO TO 89
C C
C C------------------------------------------------------
C       ELSE IF(COMAND.EQ.'O   ') THEN
C         CALL OPTGET


C C------------------------------------------------------
C       ELSE
C         WRITE(*,*)
C         WRITE(*,*) '* Option not recognized'
C C
C       ENDIF
C       GO TO 800
C
      RETURN
      END ! OPER




      SUBROUTINE CONLST(IR)
      INCLUDE 'AVL.INC'
C
      CHARACTER*4 CHSS
C
      WRITE(*,1010)
C
      DO IV = 1, NVTOT
        IC = ICON(IV,IR)
        CHSS = '  '
        DO JV = 1, NVTOT
          IF(IV.NE.JV .AND. ICON(IV,IR).EQ.ICON(JV,IR)) THEN
           CHSS = '**'
          ENDIF
        ENDDO
        WRITE(*,1020) VARKEY(IV), CONNAM(IC), CONVAL(IC,IR), CHSS
      ENDDO
C
      WRITE(*,1030)
      RETURN
C      
 1010 FORMAT(
     &  /'  variable          constraint              '
     &  /'  ------------      ------------------------')
 10200 FORMAT(
     &   '  ',A,'  ->  ', A, '=', G12.4, 1X, A)
 1030 FORMAT(
     &   '  ------------      ------------------------')
      END ! CONLST



      SUBROUTINE CONSET(COMAND,COMARG) !,LMATCH,IR)
      INCLUDE 'AVL.INC'
      CHARACTER*(*) COMAND, COMARG
      LOGICAL LMATCH
      INTEGER IR
      CHARACTER*120 PROMPT
      CHARACTER*4 ARROW
      REAL    RINP(20)
      INTEGER IINP(20)
      LOGICAL ERROR
      LMATCH=.True.
      IR=1





C       WRITE(*,*) '------------------CONSET 1-----------------------'
C       WRITE(*,*) COMAND, 'COMAND'
C       WRITE(*,*) COMARG, 'COMARG'
C       WRITE(*,*) '-----------------------------------------------'


C
C---- for control variable, first number in arg string should be part of command
      IF(COMAND(1:2) .EQ. 'D ') THEN
       COMAND(2:3) = COMARG(1:2)
       COMARG(1:2) = '  '
       CALL STRIP(COMARG,NARG)
      ENDIF
C
C---- length of non-blank part of command, if any
      KCLEN = INDEX(COMAND,' ') - 1
      IF(KCLEN.LE.0) KCLEN = LEN(COMAND)
C
C---- test command against variable keys, using only non-blank part of command
      DO IV = 1, NVTOT
        KVLEN = INDEX(VARKEY(IV),' ') - 1
        IF(KCLEN .EQ. KVLEN .AND.
     &     COMAND(1:KCLEN) .EQ. VARKEY(IV)(1:KVLEN)) GO TO 16
      ENDDO
C
C---- no variable key matched... go test for regular commands
      LMATCH = .FALSE.
      RETURN
C
C------------------------------------------------------
C---- found a variable-key match!
 16   CONTINUE
      LMATCH = .TRUE.
      CALL TOUPER(COMARG) 
C
C---- see if constraint was already specified as second command argument
      KCLEN = INDEX(COMARG,' ') - 1
      IF(KCLEN.LE.0) KCLEN = LEN(COMARG)
      KCLEN = MIN( KCLEN , LEN(CONKEY(1)) )
      DO IC = 1, NCTOT
        IF(COMARG(1:KCLEN) .EQ. CONKEY(IC)(1:KCLEN)) GO TO 18
      ENDDO
C
C---- constraint not given... get it from constraint-selection menu
      WRITE(*,1081)
      DO IC = 1, NCTOT
        IF(IC.EQ.ICON(IV,IR)) THEN
         ARROW = '->  '
        ELSE
         ARROW = '    '
        ENDIF
        WRITE(*,1082) ARROW, CONKEY(IC), CONNAM(IC), CONVAL(IC,IR)
      ENDDO
 1081 FORMAT(/'       constraint            value     '
     &       /'      - - - - - - - - - - - - - - - - -')
 1082 FORMAT( '   ', A, A, 2X, A, '=', G12.4)
C
      PROMPT= '      Select new  constraint,value  for '
     &        // VARNAM(IV) // '^'
      CALL ASKC(PROMPT,COMAND,COMARG)
      IF(COMAND.EQ.' ') RETURN
C
      IF(COMAND(1:2) .EQ. 'D ') THEN
       COMAND(2:3) = COMARG(1:2)
       COMARG(1:2) = '  '
       CALL STRIP(COMARG,NARG)
      ENDIF
C
C---- try to parse command again
      COMARG = COMAND(1:3) // ' ' // COMARG
      GO TO 16
C
C----------------------------------
C---- pick up here to set new constraint
 18   CONTINUE
C
C---- set new constraint index for selected variable IV
      ICON(IV,IR) = IC
C
C---- see if constraint value was already specified in command argument
      NINP = 1
      CALL GETFLT(COMARG(KCLEN+1:LEN(COMARG)),RINP,NINP,ERROR)
      IF(ERROR) NINP = 0
C
      IF(NINP.GE.1) THEN
C----- yep...  set constraint value to command argument
       CONVAL(IC,IR) = RINP(1)
      ELSE
C----- nope... get constraint value from user (current value is the default)
 19    WRITE(*,1090) CONNAM(IC), CONVAL(IC,IR)
 1090  FORMAT(/' Enter specified ', A,':', G12.4)
       CALL READR(1,CONVAL(IC,IR),ERROR)
       IF(ERROR) GO TO 19
      ENDIF
C
C---- go back to OPER menu
      RETURN
      END ! CONSET





      SUBROUTINE EXEC(NITER,INFO,IR)
C---------------------------------------------------
C     Solves for the flow condition specified by 
C     the global operating parameters:
C
C       CONVAL(ICALFA)     alpha (deg)
C       CONVAL(ICBETA)     beta  (deg)
C       CONVAL(ICROTX)     roll_rate * Bref / 2V
C       CONVAL(ICROTY)    pitch_rate * Cref / 2V
C       CONVAL(ICROTZ)      yaw_rate * Bref / 2V
C        .
C        .
C
C---------------------------------------------------
      INCLUDE 'AVL.INC'
      REAL VSYS(IVMAX,IVMAX), VRES(IVMAX), DDC(NDMAX), WORK(IVMAX)
      INTEGER IVSYS(IVMAX)
      real t0, t1, t2, t3, t4, t5, t6, t7, t8, t9
      real t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t99

C
C---- convergence epsilon, max angle limit (90 deg in radians)
      ! EPS replaced by EXEC_TOL
      ! DATA EPS, DMAX / 0.00002, 1.5708 /
      DATA DMAX / 1.5708 /

C
      IF(LNASA_SA) THEN
C----- NASA Std. Stability axes, X fwd, Z down
       DIR = -1.0
      ELSE
C----- Geometric Stability axes, X aft, Z up
       DIR =  1.0
      ENDIF
C===================================================================
C---- start a new solution
      LSOL = .FALSE.
C     
      call set_par_and_cons(NITER, IR)
      
      call cpu_time(t0)
C
C---- set, factor AIC matrix and induced-velocity matrix (if they don't exist)
      CALL SETUP
      if (ltiming) then 
            call cpu_time(t1)
            write(*,*) ' SETUP time: ', t1 - t0
      end if
      IF(.NOT.LAIC) THEN
            call factor_AIC
      ENDIF
      
      if (ltiming) then 
            call cpu_time(t2)
            write(*,*) ' factorize time: ', t2 - t1
      end if
C

C
C----- set GAM_U
      if (lverbose) then
      WRITE(*,*) ' Solving for unit-freestream vortex circulations...'
      endif
      CALL GUCALC
      if (ltiming) then 
            call cpu_time(t3)
            write(*,*) ' GUCALC time: ',  t3 - t2
      end if
C
C-------------------------------------------------------------
C---- calculate initial operating state
C
C---- set VINF() vector from initial ALFA,BETA
      CALL VINFAB
      if (ltiming) then 
            call cpu_time(t4)
            write(*,*) ' VINFAB time: ', t4 - t3
      end if
C
      IF(NCONTROL.GT.0) THEN
C----- set GAM_D
       if (lverbose) then
           WRITE(*,*) ' Solving for vortex control-var sensitivities...'
       end if 
       CALL GDCALC(NCONTROL,LCONDEF,ENC_D,GAM_D)
      ENDIF
      if (ltiming) then 
            call cpu_time(t5)  
            write(*,*) ' GDCALC time: ', t5 - t4
      end if
      
C
      IF(NDESIGN.GT.0) THEN
C----- set GAM_G
      if (lverbose) then
          WRITE(*,*) ' Solving for vortex  design-var sensitivities...'
      end if
       CALL GDCALC(NDESIGN ,LDESDEF,ENC_G,GAM_G)
      ENDIF
      if (ltiming) then 
            call cpu_time(t6)  
            write(*,*) ' GDCALC time: ', t6 - t5
      end if
C
C---- sum AIC matrices to get GAM,SRC,DBL
      CALL GAMSUM
      if (ltiming) then 
            call cpu_time(t7)  
            write(*,*) ' GAMSUM time: ', t7 - t6   
      end if
C
C---- sum AIC matrices to get WC,WV
      CALL VELSUM
      if (ltiming) then 
            call cpu_time(t8)  
            write(*,*) ' VELSUM time: ', t8 - t7
      end if
C
C---- compute forces
      CALL AERO
      if (ltiming) then 
            call cpu_time(t9)  
            write(*,*) ' AERO time: ', t9 - t8
      end if
C
C---- Newton loop for operating variables
      DO 190 ITER = 1, NITER
C
        IF(LSA_RATES) THEN
C-------- rates specified in NASA stability-axes, transform to body axes
          CA = COS(ALFA)
          SA = SIN(ALFA)
          CA_A = -SA
          SA_A =  CA
         ELSE
C-------- rates specified in body-axes, no transformation
          CA = 1.0
          SA = 0.0
          CA_A = 0.
          SA_A = 0.
        ENDIF
C
        DO K=1, IVMAX
          DO L=1, IVMAX
            VSYS(K,L) = 0.
          ENDDO
        ENDDO
C
C------ set up Newton system:  set constraints for all parameters
        DO 100 IV = 1, NVTOT
C
C-------- set index and value of constraint for this parameter
          IC = ICON(IV,IR)
C
C------------------------------------
          IF    (IC.EQ.ICALFA) THEN
           VRES(IV) = ALFA - CONVAL(IC,IR)*DTR
           VSYS(IV,IVALFA) = 1.0
C
C------------------------------------
          ELSEIF(IC.EQ.ICBETA) THEN
           VRES(IV) = BETA - CONVAL(IC,IR)*DTR
           VSYS(IV,IVBETA) = 1.0
C
C------------------------------------
          ELSEIF(IC.EQ.ICROTX) THEN
           VRES(IV) = (WROT(1)*CA + WROT(3)*SA)*DIR
     &              - CONVAL(IC,IR)*2.0/BREF
           VSYS(IV,IVROTX) = CA*DIR
           VSYS(IV,IVROTZ) = SA*DIR
           VSYS(IV,IVALFA) = (WROT(1)*CA_A + WROT(3)*SA_A)*DIR
C
C------------------------------------
          ELSEIF(IC.EQ.ICROTY) THEN
           VRES(IV) = WROT(2)
     &              - CONVAL(IC,IR)*2.0/CREF
           VSYS(IV,IVROTY) = 1.0
C
C------------------------------------
          ELSEIF(IC.EQ.ICROTZ) THEN
           VRES(IV) = (WROT(3)*CA - WROT(1)*SA)*DIR
     &              - CONVAL(IC,IR)*2.0/BREF
           VSYS(IV,IVROTX) = -SA*DIR
           VSYS(IV,IVROTZ) =  CA*DIR
           VSYS(IV,IVALFA) = (WROT(3)*CA_A - WROT(1)*SA_A)*DIR
C
C------------------------------------
          ELSEIF(IC.EQ.ICCL  ) THEN
           VRES(IV) = CLTOT - CONVAL(IC,IR)
           VSYS(IV,IVALFA) = CLTOT_U(1)*VINF_A(1)
     &                     + CLTOT_U(2)*VINF_A(2)
     &                     + CLTOT_U(3)*VINF_A(3) + CLTOT_A
           VSYS(IV,IVBETA) = CLTOT_U(1)*VINF_B(1)
     &                     + CLTOT_U(2)*VINF_B(2)
     &                     + CLTOT_U(3)*VINF_B(3)
           VSYS(IV,IVROTX) = CLTOT_U(4)
           VSYS(IV,IVROTY) = CLTOT_U(5)
           VSYS(IV,IVROTZ) = CLTOT_U(6)

C
           DO N = 1, NCONTROL
             NV = IVTOT + N
             VSYS(IV,NV) = CLTOT_D(N)
           ENDDO
C
C------------------------------------
          ELSEIF(IC.EQ.ICCY  ) THEN
           VRES(IV) = CYTOT - CONVAL(IC,IR)
           VSYS(IV,IVALFA) = CYTOT_U(1)*VINF_A(1)
     &                     + CYTOT_U(2)*VINF_A(2)
     &                     + CYTOT_U(3)*VINF_A(3)
           VSYS(IV,IVBETA) = CYTOT_U(1)*VINF_B(1)
     &                     + CYTOT_U(2)*VINF_B(2)
     &                     + CYTOT_U(3)*VINF_B(3)
           VSYS(IV,IVROTX) = CYTOT_U(4)
           VSYS(IV,IVROTY) = CYTOT_U(5)
           VSYS(IV,IVROTZ) = CYTOT_U(6)
C
           DO N = 1, NCONTROL
             NV = IVTOT + N
             VSYS(IV,NV) = CYTOT_D(N)
           ENDDO
C
C------------------------------------
          ELSEIF(IC.EQ.ICMOMX) THEN
           VRES(IV) = (CMTOT(1)*CA + CMTOT(3)*SA)*DIR - CONVAL(IC,IR)
           VSYS(IV,IVALFA) = ( CMTOT_U(1,1)*VINF_A(1)
     &                        +CMTOT_U(1,2)*VINF_A(2)
     &                        +CMTOT_U(1,3)*VINF_A(3))*CA*DIR
     &                     + ( CMTOT_U(3,1)*VINF_A(1)
     &                        +CMTOT_U(3,2)*VINF_A(2)
     &                        +CMTOT_U(3,3)*VINF_A(3))*SA*DIR
     &                     + (CMTOT(1)*CA_A + CMTOT(3)*SA_A)*DIR
           VSYS(IV,IVBETA) = ( CMTOT_U(1,1)*VINF_B(1)
     &                        +CMTOT_U(1,2)*VINF_B(2)
     &                        +CMTOT_U(1,3)*VINF_B(3))*CA*DIR
     &                     + ( CMTOT_U(3,1)*VINF_B(1)
     &                        +CMTOT_U(3,2)*VINF_B(2)
     &                        +CMTOT_U(3,3)*VINF_B(3))*SA*DIR
           VSYS(IV,IVROTX) = (CMTOT_U(1,4)*CA + CMTOT_U(3,4)*SA)*DIR
           VSYS(IV,IVROTY) = (CMTOT_U(1,5)*CA + CMTOT_U(3,5)*SA)*DIR
           VSYS(IV,IVROTZ) = (CMTOT_U(1,6)*CA + CMTOT_U(3,6)*SA)*DIR
C
           DO N = 1, NCONTROL
             NV = IVTOT + N
             VSYS(IV,NV) = (CMTOT_D(1,N)*CA + CMTOT_D(3,N)*SA)*DIR
           ENDDO
C
C------------------------------------
          ELSEIF(IC.EQ.ICMOMY) THEN
           VRES(IV) = CMTOT(2) - CONVAL(IC,IR)
           VSYS(IV,IVALFA) = CMTOT_U(2,1)*VINF_A(1)
     &                     + CMTOT_U(2,2)*VINF_A(2)
     &                     + CMTOT_U(2,3)*VINF_A(3)
           VSYS(IV,IVBETA) = CMTOT_U(2,1)*VINF_B(1)
     &                     + CMTOT_U(2,2)*VINF_B(2)
     &                     + CMTOT_U(2,3)*VINF_B(3)
           VSYS(IV,IVROTX) = CMTOT_U(2,4)
           VSYS(IV,IVROTY) = CMTOT_U(2,5)
           VSYS(IV,IVROTZ) = CMTOT_U(2,6)
C
           DO N = 1, NCONTROL
             NV = IVTOT + N
             VSYS(IV,NV) = CMTOT_D(2,N)
           ENDDO
C
C------------------------------------
          ELSEIF(IC.EQ.ICMOMZ) THEN
           VRES(IV) = (CMTOT(3)*CA - CMTOT(1)*SA)*DIR - CONVAL(IC,IR)
           VSYS(IV,IVALFA) = ( CMTOT_U(3,1)*VINF_A(1)
     &                        +CMTOT_U(3,2)*VINF_A(2)
     &                        +CMTOT_U(3,3)*VINF_A(3))*CA*DIR
     &                     - ( CMTOT_U(1,1)*VINF_A(1)
     &                        +CMTOT_U(1,2)*VINF_A(2)
     &                        +CMTOT_U(1,3)*VINF_A(3))*SA*DIR
     &                     + (CMTOT(3)*CA_A - CMTOT(1)*SA_A)*DIR
           VSYS(IV,IVBETA) = ( CMTOT_U(3,1)*VINF_B(1)
     &                        +CMTOT_U(3,2)*VINF_B(2)
     &                        +CMTOT_U(3,3)*VINF_B(3))*CA*DIR
     &                     - ( CMTOT_U(1,1)*VINF_B(1)
     &                        +CMTOT_U(1,2)*VINF_B(2)
     &                        +CMTOT_U(1,3)*VINF_B(3))*SA*DIR
           VSYS(IV,IVROTX) = (CMTOT_U(3,4)*CA - CMTOT_U(1,4)*SA)*DIR
           VSYS(IV,IVROTY) = (CMTOT_U(3,5)*CA - CMTOT_U(1,5)*SA)*DIR
           VSYS(IV,IVROTZ) = (CMTOT_U(3,6)*CA - CMTOT_U(1,6)*SA)*DIR
C
           DO N = 1, NCONTROL
             NV = IVTOT + N
             VSYS(IV,NV) = (CMTOT_D(3,N)*CA - CMTOT_D(1,N)*SA)*DIR
           ENDDO
C
C------------------------------------
          ELSE
           DO N = 1, NCONTROL
             ICCON = ICTOT + N
             IVCON = IVTOT + N
             IF(IC.EQ.ICCON) THEN
              VRES(IV) = DELCON(N) - CONVAL(ICCON,IR)
              VSYS(IV,IVCON) = 1.0
              GO TO 100
             ENDIF
           ENDDO
C
           WRITE(*,*) '? Illegal constraint index: ', IC
          ENDIF
C
 100    CONTINUE
C

C        write(*,*)
C        do k = 1, nvtot
C          write(*,'(1x,40f9.4)', advance='no') (vsys(k,l), l=1, nvtot)
C          write(*,'(1x,a,1x,40f9.4)') '|', vres(k)
C        enddo
C        write(*,*)
C
C------ LU-factor,  and back-substitute RHS
        CALL LUDCMP(IVMAX,NVTOT,VSYS,IVSYS,WORK)
        if (ltiming) then 
            call cpu_time(t10)  
            write(*,*) ITER, ' LUDCMP time: ', t10 - t9
        end if

        CALL BAKSUB(IVMAX,NVTOT,VSYS,IVSYS,VRES)
        if (ltiming) then 
            call cpu_time(t11)  
            write(*,*) ITER, ' BAKSUB time: ', t11 - t10
        end if


C
C------ set Newton deltas
        DAL = -VRES(IVALFA)
        DBE = -VRES(IVBETA)
        DWX = -VRES(IVROTX)
        DWY = -VRES(IVROTY)
        DWZ = -VRES(IVROTZ)
        DO N = 1, NCONTROL
          IV = IVTOT + N
          DDC(N) = -VRES(IV)
        ENDDO
C
        IF(INFO .GE. 1) THEN
C------- display Newton deltas
         IF(ITER.EQ.1) THEN
            if (lverbose) then
           WRITE(*,*)
           WRITE(*,1902) 'iter',
     &            ' d(alpha)  ',
     &            ' d(beta)   ',
     &            ' d(pb/2V)  ',
     &            ' d(qc/2V)  ',
     &            ' d(rb/2V)  ',
     &            (DNAME(K), K=1, NCONTROL)
            endif
 1902     FORMAT(1X,A4,5A11,1X,30A11)
         ENDIF
         if (lverbose) then
       WRITE(*,1905) ITER, 
     &                 DAL/DTR, DBE/DTR, 
     &                 DWX*BREF/2.0, DWY*CREF/2.0, DWZ*BREF/2.0,
     &                 (DDC(K), K=1, NCONTROL)
         end if 
 1905    FORMAT(1X,I3,40E11.3)
        ENDIF
C
C------ limits on angles and rates
        DMAXA = DMAX
        DMAXR = 5.0*DMAX/BREF
C
C------ if changes are too big, configuration is probably untrimmable
        IF(ABS(ALFA+DAL).GT.DMAXA) THEN
         WRITE(*,*) 'Cannot trim.  Alpha too large.  a =',(ALFA+DAL)/DTR
         RETURN
        ENDIF
C
        IF(ABS(BETA+DBE).GT.DMAXA) THEN
         WRITE(*,*) 'Cannot trim.  Beta too large.  b =',(BETA+DBE)/DTR
         RETURN
        ENDIF
C
        IF(ABS(WROT(1)+DWX).GT.DMAXR) THEN
         WRITE(*,*) 'Cannot trim.  Roll rate too large.  pb/2V =', 
     &               (WROT(1)+DWX)*BREF*0.5
         RETURN
        ENDIF
C
        IF(ABS(WROT(2)+DWY).GT.DMAXR) THEN
         WRITE(*,*) 'Cannot trim.  Pitch rate too large.  qc/2V =',
     &               (WROT(2)+DWY)*CREF*0.5
         RETURN
        ENDIF
C
        IF(ABS(WROT(3)+DWZ).GT.DMAXR) THEN
         WRITE(*,*) 'Cannot trim.  Yaw rate too large.  rb/2V =',
     &               (WROT(3)+DWZ)*BREF*0.5
         RETURN
        ENDIF
C
C------ update
        ALFA  = ALFA  + DAL
        BETA  = BETA  + DBE
        WROT(1) = WROT(1) + DWX
        WROT(2) = WROT(2) + DWY
        WROT(3) = WROT(3) + DWZ
        DO K = 1, NCONTROL
          DELCON(K) = DELCON(K) + DDC(K)
        ENDDO
C
C
C------ set VINF() vector from new ALFA,BETA
        CALL VINFAB
C
        IF(NCONTROL.GT.0) THEN
C------- set new GAM_D
         CALL GDCALC(NCONTROL,LCONDEF,ENC_D,GAM_D)
        ENDIF
C
        IF(NDESIGN.GT.0) THEN
C------- set new GAM_G
         CALL GDCALC(NDESIGN ,LDESDEF,ENC_G,GAM_G)
        ENDIF
C
C------ sum AIC matrices to get GAM,SRC,DBL
        CALL GAMSUM
        if (ltiming) then 
            call cpu_time(t12)  
            write(*,*) ITER, ' other1 time: ', t12 - t11
        end if

C
C------ sum AIC matrices to get WC,WV
        CALL VELSUM
        if (ltiming) then 
            call cpu_time(t13)  
            write(*,*) ITER, ' VELSUM time: ', t13 - t12
        end if
C
C
C------ compute forces
        CALL AERO
        if (ltiming) then 
            call cpu_time(t14)  
            write(*,*) ITER, ' other2 time: ', t14 - t13
        end if
C
C
C------ convergence check
        DELMAX = MAX( ABS(DAL), 
     &                ABS(DBE),
     &                ABS(DWX*BREF/2.0),
     &                ABS(DWY*CREF/2.0),
     &                ABS(DWZ*BREF/2.0) )
        DO K = 1, NCONTROL
          DELMAX = MAX( DELMAX , ABS(DDC(K)) )
        ENDDO
C
        IF(DELMAX.LT.EXEC_TOL) THEN
         LSOL = .TRUE.
         LOBVEL = .FALSE.
C------- mark trim case as being converged
         ITRIM(IR) = IABS(ITRIM(IR))
         GO TO 191
        ENDIF
C
 190  CONTINUE
      IF(NITER.GT.0) THEN
       WRITE(*,*) 'Trim convergence failed'
       LSOL = .FALSE.
       LOBVEL = .FALSE.
       RETURN
      ENDIF
C
 191  CONTINUE
      PARVAL(IPALFA,IR) = ALFA/DTR
      PARVAL(IPBETA,IR) = BETA/DTR
      PARVAL(IPROTX,IR) = WROT(1)*0.5*BREF
      PARVAL(IPROTY,IR) = WROT(2)*0.5*CREF
      PARVAL(IPROTZ,IR) = WROT(3)*0.5*BREF
      PARVAL(IPCL  ,IR) = CLTOT
C
      LSEN = .TRUE.
      if (ltiming) then 
            call cpu_time(t99)
            write(*,*) ' TOTAL time: ', t99 - t0
      end if

      RETURN
C
      END ! EXEC



      SUBROUTINE OPTGET
C-------------------------------------------------
C     Allows toggling and setting of various 
C     printing and plotting stuff.
C-------------------------------------------------
      INCLUDE 'AVL.INC'
      CHARACTER*4 ITEMC
      CHARACTER*120 COMARG
      CHARACTER*50 SATYPE, ROTTYPE
      LOGICAL ERROR
C
      REAL    RINPUT(20)
      INTEGER IINPUT(20)
      LOGICAL LINPUT(20)
C
 1000 FORMAT(A)
C
      CALL GETSA(LNASA_SA,SATYPE,DIR)
C
 100  CONTINUE
      IF(LSA_RATES) THEN
        ROTTYPE =
     &         'Rates,moments about Stability Axes, X along Vinf'
      ELSE
        ROTTYPE =
     &         'Rates,moments about Body Axes, X along geometric X axis'
      ENDIF
C  
      WRITE(*,1110) LPTOT,LPSURF,LPSTRP,LPELE,
     &              LPHINGE,LPDERIV,
     &              LTRFORCE,LVISC,LBFORCE,LNFLD_WV,
     &              SATYPE,ROTTYPE,IZSYM,ZSYM,SAXFR,
     &              VRCOREC,VRCOREW,SRCORE
 1110   FORMAT(/'   ======================================'
     &         /'    P rint default output for...'
     &         /'        total     :  ',L2,
     &         /'        surfaces  :  ',L2,
     &         /'        strips    :  ',L2,
     &         /'        elements  :  ',L2,
     &        //'    H inge mom. output:  ',L2,
     &         /'    D erivative output:  ',L2,
     &        //'    T rail.leg forces :  ',L2,
     &         /'    V iscous forces   :  ',L2,
     &         /'    B ody forces      :  ',L2,
     &         /'    N ear-field forces:  ',L2,
     &        //'    A xis orient. :  ', A, 
     &         /'    R ate,mom axes:  ', A,
     &         /'    Z  symmetry   :  ',I2,' @ Z =',F10.4
     &         /'    S pan axis x/c:  ',F10.4
     &         /'    CC core vortex radius/chord ratio:  ',F10.5
     &         /'    CW core radius/stripwidth ratio:  ',F10.5
     &         /'    CS source-doublet core/radius ratio:  ',F10.5)
C
C
      CALL ASKC(' ..Select item to change^',ITEMC,COMARG)
C
C------------------------------------------------------
      IF    (ITEMC.EQ.'    ') THEN
        RETURN
C
C---------------------------------
      ELSEIF(ITEMC.EQ.'A   ') THEN
        LNASA_SA = .NOT.LNASA_SA
        CALL GETSA(LNASA_SA,SATYPE,DIR)
C
C---------------------------------
      ELSEIF(ITEMC.EQ.'R   ') THEN
        LSA_RATES = .NOT.LSA_RATES
C
C---------------------------------
      ELSEIF(ITEMC.EQ.'V   ') THEN
        LVISC = .NOT.LVISC
        IF(LVISC) THEN
          WRITE(*,*) 'Forces will include profile drag'
         ELSE
          WRITE(*,*) 'Forces will not include profile drag'
        ENDIF
C
C---------------------------------
      ELSEIF(ITEMC.EQ.'B   ') THEN
        LBFORCE = .NOT.LBFORCE
        IF(LBFORCE) THEN
          WRITE(*,*) 'Forces will include body forces'
         ELSE
          WRITE(*,*) 'Forces will not include body forces'
        ENDIF
C
C---------------------------------
      ELSEIF(ITEMC.EQ.'N   ') THEN
        LNFLD_WV = .NOT.LNFLD_WV
        IF(LNFLD_WV) THEN
          WRITE(*,*) 'Near-field forces will include body ind. vel.'
         ELSE
          WRITE(*,*) 'Near-field forces only use h.v. ind. vel.'
        ENDIF
C
C---------------------------------
      ELSEIF(ITEMC.EQ.'T   ') THEN
        LTRFORCE = .NOT.LTRFORCE
        IF(LTRFORCE) THEN
          WRITE(*,*) 'Forces on trailing legs will be included'
         ELSE
          WRITE(*,*) 'Forces on trailing legs will not be included'
        ENDIF
C
C---------------------------------
      ELSEIF(ITEMC.EQ.'P   ') THEN
 128   IF(COMARG(1:1).NE.' ') THEN
         NINP = 4
         CALL GETLGV(COMARG,LINPUT,NINP,ERROR)
         IF(ERROR) GO TO 130
C
         IF(NINP.GE.1) LPTOT   = LINPUT(1)
         IF(NINP.GE.2) LPSURF  = LINPUT(2)
         IF(NINP.GE.3) LPSTRP  = LINPUT(3)
         IF(NINP.GE.4) LPELE   = LINPUT(4)
C
         GO TO 100
        ENDIF
C
 130    WRITE(*,2100)
     &  'Enter print flags T/F (total,surf,strip,elem)'
 2100   FORMAT(1X,A,': ', $)
C
        READ(*,1000) COMARG
        IF(COMARG.EQ.' ') THEN
         GO TO 100
        ELSE
         GO TO 128
        ENDIF
C
C---------------------------------
      ELSEIF(ITEMC.EQ.'H   ') THEN
        LPHINGE = .NOT.LPHINGE
C
C---------------------------------
      ELSEIF(ITEMC.EQ.'D   ') THEN
        LPDERIV = .NOT.LPDERIV
C
C---------------------------------
      ELSEIF(ITEMC.EQ.'Z   ') THEN
       WRITE(*,*) ' '
       WRITE(*,*) 'Currently:'
       IF(IZSYM.EQ.0) THEN
         WRITE (*,1015)
        ELSEIF(IZSYM.GT.0) THEN
         WRITE (*,1016) ZSYM
        ELSEIF(IZSYM.LT.0) THEN
         WRITE (*,1017) ZSYM
       ENDIF
       WRITE(*,*) 'Enter symmetry flag: -1 Free surface'
       WRITE(*,*) '                      0 no Z symmetry'
       WRITE(*,*) '                      1 Ground plane'
       ZSYMIN = 0.0
       READ(*,*,ERR=100) IZSYMIN
       IF(IZSYMIN.NE.0.0) THEN
         WRITE(*,2100) 'Enter Z for symmetry plane'
         READ(*,*,ERR=100) ZSYMIN
       ENDIF
       IZSYM = IZSYMIN
       ZSYM  = ZSYMIN
       LAIC = .FALSE.
       LSRD = .FALSE.
       LVEL = .FALSE.
       LSOL = .FALSE.
       LSEN = .FALSE.
       LOBAIC = .FALSE.
       LOBVEL = .FALSE.
C
 1015  FORMAT(' Z Symmetry: No symmetry assumed')
 1016  FORMAT(' Z Symmetry: Ground plane at Zsym =',F10.4)
 1017  FORMAT(' Z Symmetry: Free surface at Zsym =',F10.4)
C
C---------------------------------
      ELSEIF(ITEMC.EQ.'S   ') THEN
        NINP = 1
        CALL GETFLT(COMARG,RINPUT,NINP,ERROR)
C
        IF(ERROR .OR. NINP.LE.0) THEN
         RINPUT(1) = SAXFR
         WRITE(*,1030) RINPUT(1)
 1030    FORMAT(/' Enter x/c location of spanwise ref. axis:', F10.4)
         CALL READR(1,RINPUT(1),ERROR)
         IF(ERROR) GO TO 100
        ENDIF
C
        SAXFR = MAX( 0.0 , MIN(1.0,RINPUT(1)) )
        CALL ENCALC
        CALL AERO
C
C---------------------------------
      ELSEIF(ITEMC.EQ.'CC  ') THEN
        NINP = 1
        CALL GETFLT(COMARG,RINPUT,NINP,ERROR)
C
        IF(ERROR .OR. NINP.LE.0) THEN
         RINPUT(1) = VRCOREC
         WRITE(*,1040) RINPUT(1)
 1040    FORMAT(/' Enter core vortex radius/strip chord:', F10.5)
         CALL READR(1,RINPUT,ERROR)
         IF(ERROR) GO TO 100
        ENDIF
C
        VRCOREC = MAX( 0.0 , MIN(5.0,RINPUT(1)) )
        CALL ENCALC
        LAIC = .FALSE.
        LSRD = .FALSE.
        LVEL = .FALSE.
        LSOL = .FALSE.
        LSEN = .FALSE.
        LOBAIC = .FALSE.
        LOBVEL = .FALSE.
C
C---------------------------------
      ELSEIF(ITEMC.EQ.'CW  ') THEN
        NINP = 1
        CALL GETFLT(COMARG,RINPUT,NINP,ERROR)
C
        IF(ERROR .OR. NINP.LE.0) THEN
         RINPUT(1) = VRCOREW
         WRITE(*,1042) RINPUT(1)
 1042    FORMAT(/' Enter core vortex radius/strip width:', F10.5)
         CALL READR(1,RINPUT,ERROR)
         IF(ERROR) GO TO 100
        ENDIF
C
        VRCOREW = MAX( 0.0 , MIN(10.0,RINPUT(1)) )
        CALL ENCALC
        LAIC = .FALSE.
        LSRD = .FALSE.
        LVEL = .FALSE.
        LSOL = .FALSE.
        LSEN = .FALSE.

C---------------------------------
      ELSEIF(ITEMC.EQ.'CS  ') THEN
        NINP = 1
        CALL GETFLT(COMARG,RINPUT,NINP,ERROR)
C
        IF(ERROR .OR. NINP.LE.0) THEN
         RINPUT(1) = SRCORE
         WRITE(*,1044) RINPUT(1)
 1044    FORMAT(/' Enter core source radius/strip width:', F10.5)
         CALL READR(1,RINPUT,ERROR)
         IF(ERROR) GO TO 100
        ENDIF
C
C---- SRCCORE>0 uses body radius, SRCORE<0 uses source line spacing
        SRCORE = RINPUT(1)
        SRCORE = SIGN( MIN(5.0,ABS(SRCORE)), SRCORE )
        CALL ENCALC
        LAIC = .FALSE.
        LSRD = .FALSE.
        LVEL = .FALSE.
        LSOL = .FALSE.
        LSEN = .FALSE.
        LOBAIC = .FALSE.
        LOBVEL = .FALSE.
C
C---------------------------------
      ELSE
        WRITE(*,*) 'Item not recognized'
        GO TO 100
      ENDIF
      GO TO 100
C
      END ! OPTGET



      SUBROUTINE CFRAC(IRUN,NRUN,CPR,NPR)
      CHARACTER*(*) CPR
C
      IZERO = ICHAR('0')
      ITEN = IRUN/10
      IONE = IRUN - 10*(IRUN/10)
      IF(ITEN.LE.0) THEN
       CPR = CHAR(IZERO+IONE) // '/'
      ELSE
       CPR = CHAR(IZERO+ITEN) // CHAR(IZERO+IONE) // '/'
      ENDIF
C
      NPR = INDEX(CPR,'/')
      ITEN = NRUN/10
      IONE = NRUN - 10*(NRUN/10)
      IF(ITEN.LE.0) THEN
       CPR = CPR(1:NPR)
     &    // CHAR(IZERO+IONE) // '^'
      ELSE
       CPR = CPR(1:NPR)
     &    // CHAR(IZERO+ITEN) // CHAR(IZERO+IONE) // '^'
      ENDIF
C
      NPR = INDEX(CPR,'^') - 1
C
      RETURN
      END ! CFRAC



      SUBROUTINE RCOPY(IRSET,IR)
      INCLUDE 'AVL.INC'
C
      DO IV = 1, NVTOT
        ICON(IV,IRSET) = ICON(IV,IR)
      ENDDO
      DO IC = 1, NCTOT
        CONVAL(IC,IRSET) = CONVAL(IC,IR)
      ENDDO
      DO IP = 1, NPTOT
        PARVAL(IP,IRSET) = PARVAL(IP,IR)
      ENDDO
C
      RTITLE(IRSET) = RTITLE(IR)
      ITRIM(IRSET) = ITRIM(IR)
      NEIGEN(IRSET) = NEIGEN(IR)
C
      DO KE = 1, JEMAX
        EVAL(KE,IRSET) = EVAL(KE,IR)
        DO JE = 1, JEMAX
          EVEC(KE,JE,IRSET) = EVEC(KE,JE,IR)
        ENDDO
      ENDDO
C
      RETURN
      END ! RCOPY





      SUBROUTINE GETFILE(LU,FNAME)
      CHARACTER(*) FNAME
C
      CHARACTER*1 ANS, DUMMY
C
 1000 FORMAT(A)
C
      IF(FNAME.EQ.' ') THEN
       CALL ASKS('Enter filename, or <return> for screen output^',FNAME)
      ENDIF
C
      IF(FNAME.EQ.' ') THEN
       LU = 6
       RETURN
C
      ELSE
       LU = 11
       OPEN(LU,FILE=FNAME,STATUS='OLD',FORM='FORMATTED',ERR=44)
       WRITE(*,*) 
       WRITE(*,*) 'File exists.  Append/Overwrite/Cancel  (A/O/C)?  C'
       READ(*,1000) ANS
       IF (INDEX('Aa',ANS).NE.0) THEN
C---- reopen file with append status HHY 4/17/18
         CLOSE(LU)
         OPEN(LU,FILE=FNAME,STATUS='OLD',FORM='FORMATTED',
     &        ACCESS='APPEND',ERR=44)
C---- old append code (re-reads to EOF)
cc 40     CONTINUE
cc        READ(LU,1000,END=42) DUMMY
cc        GO TO 40
cc 42     CONTINUE
       ELSEIF(INDEX('Oo',ANS).NE.0) THEN
         REWIND(LU) 
       ELSE
         CLOSE(LU)
         LU = 0
       ENDIF
       RETURN
C
 44    OPEN(LU,FILE=FNAME,STATUS='UNKNOWN',FORM='FORMATTED',ERR=48)
       REWIND(LU)
       RETURN
C
 48    CONTINUE
       LU = -1
       RETURN
      ENDIF
      END ! GETFILE



C ==================== addition wrapper helper programs ==============

      SUBROUTINE calc_stab_derivs
C C---------------------------------------------------------
C C     Calculates and outputs stability derivative matrix
C C     for current ALFA, BETA for both the stability axis and the body axis.
C C---------------------------------------------------------
      INCLUDE 'AVL.INC'
      CHARACTER*50 SATYPE
      REAL WROT_RX(3), WROT_RZ(3), WROT_A(3)
      REAL
     & CRSAX_U(NUMAX),CMSAX_U(NUMAX),CNSAX_U(NUMAX),
     & CRSAX_G(NGMAX),CMSAX_G(NGMAX),CNSAX_G(NGMAX)
C
      CALL GETSA(LNASA_SA,SATYPE,DIR)
C
C---- set freestream velocity components from alpha, beta
      ! CALL VINFAB
C
C---- calculate forces and sensitivities
      ! CALL AERO
C
C---- set stability-axes rates (RX,RY,RZ) in terms of body-axes rates
      CA = COS(ALFA)
      SA = SIN(ALFA)
C
      RX = (WROT(1)*CA + WROT(3)*SA) * DIR
      RY =  WROT(2)
      RZ = (WROT(3)*CA - WROT(1)*SA) * DIR
C
C---- now vice-versa, and set sensitivities (which is what's really needed)
cc    WROT(1)    =  RX*CA - RZ*SA
cc    WROT(2)    =  RY
cc    WROT(3)    =  RZ*CA + RX*SA
C
      WROT_RX(1) = CA     * DIR
      WROT_RX(2) = 0.
      WROT_RX(3) =     SA * DIR
C
      WROT_RZ(1) =    -SA * DIR
      WROT_RZ(2) = 0.
      WROT_RZ(3) = CA     * DIR
C
      WROT_A(1)  = -RX*SA - RZ*CA   !!! = -WROT(3)
      WROT_A(2)  =  0.
      WROT_A(3)  = -RZ*SA + RX*CA   !!! =  WROT(1)
C
C
      CRSAX_A = -CMTOT(1)*SA + CMTOT(3)*CA
      CNSAX_A = -CMTOT(3)*SA - CMTOT(1)*CA
C
      DO K = 1, 6
        CRSAX_U(K) = CMTOT_U(1,K)*CA + CMTOT_U(3,K)*SA
        CMSAX_U(K) = CMTOT_U(2,K)              
        CNSAX_U(K) = CMTOT_U(3,K)*CA - CMTOT_U(1,K)*SA  
      ENDDO
      DO K = 1, NCONTROL
        CRSAX_D(K) = DIR*(CMTOT_D(1,K)*CA + CMTOT_D(3,K)*SA)
        CMSAX_D(K) = CMTOT_D(2,K)              
        CNSAX_D(K) = DIR*(CMTOT_D(3,K)*CA - CMTOT_D(1,K)*SA)
      ENDDO

      DO K = 1, NDESIGN
        CRSAX_G(K) = CMTOT_G(1,K)*CA + CMTOT_G(3,K)*SA
        CMSAX_G(K) = CMTOT_G(2,K)              
        CNSAX_G(K) = CMTOT_G(3,K)*CA - CMTOT_G(1,K)*SA  
      ENDDO

C
C---- set force derivatives in stability axes
      CLTOT_AL = CLTOT_U(1)*VINF_A(1) + CLTOT_U(4)*WROT_A(1)
     &         + CLTOT_U(2)*VINF_A(2) + CLTOT_U(5)*WROT_A(2)
     &         + CLTOT_U(3)*VINF_A(3) + CLTOT_U(6)*WROT_A(3) + CLTOT_A
      CLTOT_BE = CLTOT_U(1)*VINF_B(1)
     &         + CLTOT_U(2)*VINF_B(2)
     &         + CLTOT_U(3)*VINF_B(3)
      CLTOT_RX = CLTOT_U(4)*WROT_RX(1) + CLTOT_U(6)*WROT_RX(3)
      CLTOT_RY = CLTOT_U(5)
      CLTOT_RZ = CLTOT_U(6)*WROT_RZ(3) + CLTOT_U(4)*WROT_RZ(1)
C
      CYTOT_AL = CYTOT_U(1)*VINF_A(1) + CYTOT_U(4)*WROT_A(1)
     &         + CYTOT_U(2)*VINF_A(2) + CYTOT_U(5)*WROT_A(2)
     &         + CYTOT_U(3)*VINF_A(3) + CYTOT_U(6)*WROT_A(3)
      CYTOT_BE = CYTOT_U(1)*VINF_B(1)
     &         + CYTOT_U(2)*VINF_B(2)
     &         + CYTOT_U(3)*VINF_B(3)
      CYTOT_RX = CYTOT_U(4)*WROT_RX(1) + CYTOT_U(6)*WROT_RX(3)
      CYTOT_RY = CYTOT_U(5)
      CYTOT_RZ = CYTOT_U(6)*WROT_RZ(3) + CYTOT_U(4)*WROT_RZ(1)
C
      CDTOT_AL = CDTOT_U(1)*VINF_A(1) + CDTOT_U(4)*WROT_A(1)
     &         + CDTOT_U(2)*VINF_A(2) + CDTOT_U(5)*WROT_A(2)
     &         + CDTOT_U(3)*VINF_A(3) + CDTOT_U(6)*WROT_A(3) + CDTOT_A
      CDTOT_BE = CDTOT_U(1)*VINF_B(1)
     &         + CDTOT_U(2)*VINF_B(2)
     &         + CDTOT_U(3)*VINF_B(3)
      CDTOT_RX = CDTOT_U(4)*WROT_RX(1) + CDTOT_U(6)*WROT_RX(3)
      CDTOT_RY = CDTOT_U(5)
      CDTOT_RZ = CDTOT_U(6)*WROT_RZ(3) + CDTOT_U(4)*WROT_RZ(1)
C
      CRTOT_AL = CRSAX_U(1)*VINF_A(1) + CRSAX_U(4)*WROT_A(1)
     &         + CRSAX_U(2)*VINF_A(2) + CRSAX_U(5)*WROT_A(2)
     &         + CRSAX_U(3)*VINF_A(3) + CRSAX_U(6)*WROT_A(3) + CRSAX_A
      CRTOT_BE = CRSAX_U(1)*VINF_B(1)
     &         + CRSAX_U(2)*VINF_B(2)
     &         + CRSAX_U(3)*VINF_B(3)
      CRTOT_RX = CRSAX_U(4)*WROT_RX(1) + CRSAX_U(6)*WROT_RX(3)
      CRTOT_RY = CRSAX_U(5)
      CRTOT_RZ = CRSAX_U(6)*WROT_RZ(3) + CRSAX_U(4)*WROT_RZ(1)
C
      CMTOT_AL = CMSAX_U(1)*VINF_A(1) + CMSAX_U(4)*WROT_A(1)
     &         + CMSAX_U(2)*VINF_A(2) + CMSAX_U(5)*WROT_A(2)
     &         + CMSAX_U(3)*VINF_A(3) + CMSAX_U(6)*WROT_A(3)
      CMTOT_BE = CMSAX_U(1)*VINF_B(1)
     &         + CMSAX_U(2)*VINF_B(2)
     &         + CMSAX_U(3)*VINF_B(3)
      CMTOT_RX = CMSAX_U(4)*WROT_RX(1) + CMSAX_U(6)*WROT_RX(3)
      CMTOT_RY = CMSAX_U(5)
      CMTOT_RZ = CMSAX_U(6)*WROT_RZ(3) + CMSAX_U(4)*WROT_RZ(1)
C
      CNTOT_AL = CNSAX_U(1)*VINF_A(1) + CNSAX_U(4)*WROT_A(1)
     &         + CNSAX_U(2)*VINF_A(2) + CNSAX_U(5)*WROT_A(2)
     &         + CNSAX_U(3)*VINF_A(3) + CNSAX_U(6)*WROT_A(3) + CNSAX_A
      CNTOT_BE = CNSAX_U(1)*VINF_B(1)
     &         + CNSAX_U(2)*VINF_B(2)
     &         + CNSAX_U(3)*VINF_B(3)
      CNTOT_RX = CNSAX_U(4)*WROT_RX(1) + CNSAX_U(6)*WROT_RX(3)
      CNTOT_RY = CNSAX_U(5)
      CNTOT_RZ = CNSAX_U(6)*WROT_RZ(3) + CNSAX_U(4)*WROT_RZ(1)
C
C        

      
      ! apply the facotors to the outputs as done in the print statements of DERMATS
      CRTOT_AL = DIR*CRTOT_AL
      CRTOT_BE = DIR*CRTOT_BE
      CNTOT_AL = DIR*CNTOT_AL
      CNTOT_BE = DIR*CNTOT_BE
      
      CLTOT_RX = CLTOT_RX*2.0/BREF 
      CLTOT_RY = CLTOT_RY*2.0/CREF 
      CLTOT_RZ = CLTOT_RZ*2.0/BREF

      CYTOT_RX = CYTOT_RX*2.0/BREF
      CYTOT_RY = CYTOT_RY*2.0/CREF
      CYTOT_RZ = CYTOT_RZ*2.0/BREF      
      
      CDTOT_RX = CDTOT_RX*2.0/BREF 
      CDTOT_RY = CDTOT_RY*2.0/CREF 
      CDTOT_RZ = CDTOT_RZ*2.0/BREF

      CRTOT_RX = DIR*CRTOT_RX*2.0/BREF 
      CRTOT_RY = DIR*CRTOT_RY*2.0/CREF
      CRTOT_RZ = DIR*CRTOT_RZ*2.0/BREF
      CMTOT_RX = CMTOT_RX*2.0/BREF
      CMTOT_RY = CMTOT_RY*2.0/CREF
      CMTOT_RZ = CMTOT_RZ*2.0/BREF
      CNTOT_RX = DIR*CNTOT_RX*2.0/BREF
      CNTOT_RY = DIR*CNTOT_RY*2.0/CREF
      CNTOT_RZ = DIR*CNTOT_RZ*2.0/BREF

      
      IF(CLTOT_AL .NE. 0.0) THEN
      !  XNP = XYZREF(1) - CREF*CMTOT_AL/CLTOT_AL
      !  SM = (XNP - XYZREF(1))/CREF This is the same as below
       SM = -CMTOT_AL/CLTOT_AL
       XNP = XYZREF(1) + CREF*SM


C        WRITE(*,8401) XNP
 8401  FORMAT(/'Neutral point  Xnp =', F11.6)
      ENDIF
C
      IF(ABS(CRTOT_RZ*CNTOT_BE) .GT. 0.0) THEN
       BB = CRTOT_BE*CNTOT_RZ / (CRTOT_RZ*CNTOT_BE)
C        WRITE(LU,8402) BB 
C  8402  FORMAT(/' Clb Cnr / Clr Cnb  =', F11.6,
C      &    '    (  > 1 if spirally stable )')
      ENDIF
      IF(ABS(CRTOT_BE) .GT. 0.0) THEN
       ! this is the lateral stablity parameter suggested by John Yost
       RR = CNTOT_BE/abs(CRTOT_BE)
      ENDIF
      
      ! ---- body axis derivatives
      ! these expressions are taken from DERMATB
      
      !'  axial   vel. u', ' sideslip vel. v', '  normal  vel. w'
      
      CXTOT_U_BA(1) = -    CFTOT_U(1,1)
      CXTOT_U_BA(2) = -DIR*CFTOT_U(1,2)
      CXTOT_U_BA(3) = -    CFTOT_U(1,3)
      
      CYTOT_U_BA(1) = -DIR*CFTOT_U(2,1)
      CYTOT_U_BA(2) = -    CFTOT_U(2,2)
      CYTOT_U_BA(3) = -DIR*CFTOT_U(2,3)
      
      CZTOT_U_BA(1) = -    CFTOT_U(3,1)
      CZTOT_U_BA(2) = -DIR*CFTOT_U(3,2)
      CZTOT_U_BA(3) = -    CFTOT_U(3,3)
      
      CRTOT_U_BA(1) = -    CMTOT_U(1,1)
      CRTOT_U_BA(2) = -DIR*CMTOT_U(1,2)
      CRTOT_U_BA(3) = -    CMTOT_U(1,3)
      
      CMTOT_U_BA(1) = -DIR*CMTOT_U(2,1)
      CMTOT_U_BA(2) = -    CMTOT_U(2,2)
      CMTOT_U_BA(3) = -DIR*CMTOT_U(2,3)
      
      CNTOT_U_BA(1) = -    CMTOT_U(3,1)
      CNTOT_U_BA(2) = -DIR*CMTOT_U(3,2)
      CNTOT_U_BA(3) = -    CMTOT_U(3,3)
      
      !'    roll rate  p','   pitch rate  q','     yaw rate  r'
      CXTOT_U_BA(4) =     CFTOT_U(1,4)*2.0/BREF
      CXTOT_U_BA(5) = DIR*CFTOT_U(1,5)*2.0/CREF
      CXTOT_U_BA(6) =     CFTOT_U(1,6)*2.0/BREF
      
      CYTOT_U_BA(4) = DIR*CFTOT_U(2,4)*2.0/BREF
      CYTOT_U_BA(5) =     CFTOT_U(2,5)*2.0/CREF
      CYTOT_U_BA(6) = DIR*CFTOT_U(2,6)*2.0/BREF
      
      CZTOT_U_BA(4) =      CFTOT_U(3,4)*2.0/BREF
      CZTOT_U_BA(5) =  DIR*CFTOT_U(3,5)*2.0/CREF
      CZTOT_U_BA(6) =      CFTOT_U(3,6)*2.0/BREF
      
      CRTOT_U_BA(4) =     CMTOT_U(1,4)*2.0/BREF
      CRTOT_U_BA(5) = DIR*CMTOT_U(1,5)*2.0/CREF
      CRTOT_U_BA(6) =     CMTOT_U(1,6)*2.0/BREF
      
      CMTOT_U_BA(4) = DIR*CMTOT_U(2,4)*2.0/BREF
      CMTOT_U_BA(5) =     CMTOT_U(2,5)*2.0/CREF
      CMTOT_U_BA(6) = DIR*CMTOT_U(2,6)*2.0/BREF
      
      CNTOT_U_BA(4) =     CMTOT_U(3,4)*2.0/BREF
      CNTOT_U_BA(5) = DIR*CMTOT_U(3,5)*2.0/CREF
      CNTOT_U_BA(6) =     CMTOT_U(3,6)*2.0/BREF
      
      ! control surfaces
      
      do K=1, NCONTROL
            CXTOT_D_BA(K) = DIR*CFTOT_D(1,K)
            CYTOT_D_BA(K) =     CFTOT_D(2,K)
            CZTOT_D_BA(K) = DIR*CFTOT_D(3,K)
            CRTOT_D_BA(K) = DIR*CMTOT_D(1,K)
            CMTOT_D_BA(K) =     CMTOT_D(2,K)
            CNTOT_D_BA(K) = DIR*CMTOT_D(3,K)
      enddo
      
      ! design variables
      
      do K=1, NDESIGN
            CXTOT_G_BA(K) = DIR*CFTOT_G(1,K)
            CYTOT_G_BA(K) =     CFTOT_G(2,K)
            CZTOT_G_BA(K) = DIR*CFTOT_G(3,K)
            CRTOT_G_BA(K) = DIR*CMTOT_G(1,K)
            CMTOT_G_BA(K) =     CMTOT_G(2,K)
            CNTOT_G_BA(K) = DIR*CMTOT_G(3,K)
      enddo
      
C C
      RETURN
      
      
      
      END ! calc_stab_derivs

C ============= Added to AVL ===============

      subroutine exec_rhs
      use avl_heap_inc
      include "AVL.INC"
      integer i, IU
      ! CALL EXEC(10,1,1)
      call set_par_and_cons(NITMAX, IRUN)
      
      CALL SETUP
      IF(.NOT.LAIC) THEN
            call factor_AIC
      ENDIF
      
      CALL VINFAB
      
      DO IU = 1,6
            call set_vel_rhs_u(IU)
            do i = 1,NVOR
                  GAM_U(i,IU) = RHS_U(i, IU)
            enddo
            CALL BAKSUB(NVMAX,NVOR,AICN_LU,IAPIV,GAM_U(:,IU))
      enddo
      
      call set_vel_rhs
      
C---- copy RHS vector into GAM that will be used for soluiton
      do i = 1,NVOR
            GAM(i) = RHS(i)
      enddo

      CALL BAKSUB(NVMAX,NVOR,AICN_LU,IAPIV,GAM)
      
      
      IF(NCONTROL.GT.0) THEN
C------- set new GAM_D
            CALL GDCALC(NCONTROL,LCONDEF,ENC_D,GAM_D)
      ENDIF

      end !subroutine exec_rhs


C ======================== res and Adjoint for GAM ========      
      
      subroutine get_res
      use avl_heap_inc
      INCLUDE "AVL.INC"
      integer I, IC
      real RHS_D(NVOR)
      call set_par_and_cons(NITMAX, IRUN)
C---  
      ! Do not use this routine in the sovler
      ! IF(.NOT.LAIC) THEN
      !      CALL build_AIC
      ! end if
      CALL build_AIC
      AMACH = MACH
      BETM = SQRT(1.0 - AMACH**2)
       CALL VVOR(BETM,IYSYM,YSYM,IZSYM,ZSYM,
     &           VRCOREC,VRCOREW,
     &           NVOR,RV1,RV2,LVCOMP,CHORDV,
     &           NVOR,RV ,    LVCOMP,.TRUE.,
     &           WV_GAM,NVMAX)

C---- set VINF() vector from initial ALFA,BETA
      CALL VINFAB
      
      DO IU = 1,6
            
            call mat_prod(AICN, GAM_U(:,IU), NVOR, res_U(:,IU))
            
            call set_vel_rhs_u(IU)
            do I = 1, NVOR
                  res_U(I,IU) = res_U(I,IU) - RHS_U(I,IU)
            enddo 
            
      enddo
      
            
      call set_vel_rhs

      
      call mat_prod(AICN, GAM, NVOR, res)
C---- add the RHS vector to the residual
      !$AD II-LOOP
      do I = 1, NVOR
            res(I) = res(I) - RHS(I)
      enddo 
      
      !$AD II-LOOP
      DO IC = 1, NCONTROL
C------ don't bother if this control variable is undefined
            IF (lcondef(ic)) THEN
            
                  call mat_prod(AICN, GAM_D(:,IC), NVOR, res_D(:, IC))
                  !  RHS_D(:) = 0.D0
                  call set_gam_d_rhs(IC, ENC_D, RHS_D)
                  !$AD II-LOOP
                  do I = 1, NVOR
                        res_D(I, IC) = res_D(I, IC) - RHS_D(I)
                  enddo      
            endif 
      ENDDO

      end !get_res

      
      
      subroutine solve_adjoint(solve_stab_deriv_adj, solve_con_surf_adj)
      use avl_heap_inc
      use avl_heap_diff_inc
      include "AVL.INC"
      include "AVL_ad_seeds.inc"
      integer i 
      logical :: solve_stab_deriv_adj, solve_con_surf_adj
      
      CALL SETUP
      IF(.NOT.LAIC) THEN
            call factor_AIC
      ENDIF
      

      do i =1,NVOR
            RES_diff(i) = GAM_diff(i)
      enddo

      CALL BAKSUBTRANS(NVMAX,NVOR,AICN_LU,IAPIV,RES_diff)
      
      if (solve_con_surf_adj) then 
      DO IC = 1, NCONTROL
            do i =1,NVOR
                  RES_D_diff(i,IC) = GAM_D_diff(i,IC)
            enddo
            CALL BAKSUBTRANS(NVMAX,NVOR,AICN_LU,IAPIV,RES_D_diff(:,IC))
      enddo
      endif 

      if (solve_stab_deriv_adj) then 
      DO IU = 1,6
            do i =1,NVOR
                  RES_U_diff(i,IU) = GAM_U_diff(i,IU)
            enddo
            
            CALL BAKSUBTRANS(NVMAX,NVOR,AICN_LU,IAPIV,RES_U_diff(:,IU))
      enddo
      endif
      
      end !subroutine solve_adjoint