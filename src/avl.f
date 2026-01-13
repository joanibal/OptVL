C***********************************************************************
C    Module:  avl.f
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

C       PROGRAM AVL
      SUBROUTINE AVL
C C=======================================================================
C C     3-D Vortex Lattice code.
C C     See file avl_doc.txt for user guide.
C C     See file version_notes.txt for most recent changes.
C C=======================================================================
      INCLUDE 'AVL.INC'
C      INCLUDE 'AVLPLT.INC'
      LOGICAL ERROR, LINPFILE
C
      CHARACTER*4 COMAND
      CHARACTER*128 COMARG
      CHARACTER*256 FNNEW
C
      REAL    RINPUT(20)
      INTEGER IINPUT(20)
C

      VERSION = 3.40

C  1000 FORMAT(A)
C
C       WRITE (*,1005) VERSION
C  1005 FORMAT (
C      &  /' ==================================================='
C      &  /'  Athena Vortex Lattice  Program      Version ',F5.2
C      &  /'  Copyright (C) 2002   Mark Drela, Harold Youngren'
C      & //'  This software comes with ABSOLUTELY NO WARRANTY,' 
C      &  /'    subject to the GNU General Public License.'
C      & //'  Caveat computor'
C      &  /' ===================================================')
C C
C
      PI = 4.0*ATAN(1.0)
      DTR = PI/180.0
C
C---- Flag for having valid input data from file
      LINPFILE = .FALSE.
C
C---- logical units
      LUINP = 4   ! configuration file
      LURUN = 7   ! run case file
      LUMAS = 8   ! mass file
      LUPRM = 9   ! parameter file
      LUOUT = 19  ! output dump file
      LUSTD = 20  ! stability derivative dump file
      LUSYS = 22  ! dynamic system matrix dump file
C
C---- set basic defaults
      CALL DEFINI
      CALL MASINI
C
C---- initialize Xplot, and AVL plot stuff
C       CALL PLINIT
      CPFAC = MIN(0.4*CREF,0.1*BREF)  / CREF

C
C
C---- Read a new input geometry from input file
C       CALL GETARG0(1,FILDEF)
C C
C       IF(FILDEF.NE.' ') THEN
C        CALL INPUT(LUINP,FILDEF,ERROR)
C C
C C----- no valid geometry... skip reading run and mass files
C        IF(ERROR) GO TO 100
C C
C C----- This is moved to the loadgeo function
CC----- we have an input dataset to process
C       LINPFILE = .TRUE.
CC
C C----- set up all parameters
C        CALL PARSET
C
C----- process geometry to define strip and vortex data
       LPLTNEW = .TRUE.
       ! TODO remove?
       CALL ENCALC
C
C----- initialize state
       CALL VARINI
C C
C       ELSE
C C----- no geometry... skip reading run and mass files
C        GO TO 100
C C
C       ENDIF
C C
C-------------------------------------------------------------------
C C---- try to read mass file
C       CALL GETARG0(3,FMSDEF)
C       IF(FMSDEF.EQ.' ') THEN
C        KDOT = INDEX(FILDEF,'.')
C        IF(KDOT.EQ.0) THEN
C         CALL STRIP(FILDEF,LENF)
C         FMSDEF = FILDEF(1:LENF) // '.mass'
C        ELSE
C         FMSDEF = FILDEF(1:KDOT) // 'mass'
C        ENDIF
C       ENDIF
C C
C       CALL STRIP(FMSDEF,NMS)
C       WRITE(*,*) 
C       WRITE(*,*)
C      & '---------------------------------------------------------------'
C       WRITE(*,*) 'Trying to read file: ', FMSDEF(1:NMS), '  ...'
C       CALL MASGET(LUMAS,FMSDEF,ERROR)
C C
C       IF(ERROR) THEN
C        WRITE(*,*) 'Internal mass defaults used'
C        CALL MASINI
C C
C       ELSE
C        WRITE(*,*)
C        WRITE(*,*) 'Mass distribution read ...'
C        CALL MASSHO(6)
C C
C        CALL APPGET
C        WRITE(*,*) 
C      & '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
C        CALL APPSHO(6,RHO0)
C C
C       ENDIF
C C
C C-------------------------------------------------------------------
C C---- try to read run case file
C       CALL GETARG0(2,FRNDEF)
C       IF(FRNDEF.EQ.' ') THEN
C        KDOT = INDEX(FILDEF,'.')
C        IF(KDOT.EQ.0) THEN
C         CALL STRIP(FILDEF,LENF)
C         FRNDEF = FILDEF(1:LENF) // '.run'
C        ELSE
C         FRNDEF = FILDEF(1:KDOT) // 'run'
C        ENDIF
C       ENDIF
C C
C       CALL STRIP(FRNDEF,NFR)
C       WRITE(*,*)
C       WRITE(*,*)
C      & '---------------------------------------------------------------'
C       WRITE(*,*) 'Trying to read file: ', FRNDEF(1:NFR), '  ...'
C       CALL RUNGET(LURUN,FRNDEF,ERROR)
C C
C       IF(ERROR) THEN
C        WRITE(*,*) 'Internal run case defaults used'
C        CALL RUNINI
C C
C       ELSE
C        WRITE(*,1025) (IR, RTITLE(IR), IR=1, NRUN)
C  1025  FORMAT(//' Run cases read  ...',
C      &        100(/1X,I4,': ',A))
C C
C       ENDIF
C C
C C-------------------------------------------------------------------
C  100  CONTINUE
C C
C---- set up plotting parameters for geometry (if any)
C      CALL PLPARS
C C
C       WRITE(*,2000) 
C  2000 FORMAT(
C      &  /' =========================================================='
C      &  /'   Quit    Exit program'
C      & //'  .OPER    Compute operating-point run cases'
C      &  /'  .MODE    Eigenvalue analysis of run cases'
C      &  /'  .TIME    Time-domain calculations'
C      & //'   LOAD f  Read configuration input file'
C      &  /'   MASS f  Read mass distribution file'
C      &  /'   CASE f  Read run case file'
C      & //'   CINI    Clear and initialize run cases'
C      &  /'   MSET i  Apply mass file data to stored run case(s)'
C      & //'  .PLOP    Plotting options'
C      &  /'   NAME s  Specify new configuration name')
C C
C======================================================================
C C---- start of menu loop
C   500 CONTINUE
C       CALL ASKC(' AVL^',COMAND,COMARG)
C C
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
C C===============================================
C       IF(COMAND.EQ.'    ') THEN
C        GO TO 500
C C
C C===============================================
C       ELSEIF(COMAND.EQ.'?   ') THEN
C        WRITE(*,2000)
C C
C C===============================================
C       ELSEIF(COMAND.EQ.'QUIT' .OR.
C      &       COMAND.EQ.'Q   '      ) THEN
C C        CALL PLCLOSE
C        STOP
C C
C C===============================================
C       ELSEIF(COMAND.EQ.'OPER') THEN
C        CALL OPER
C C
C C===============================================
C       ELSEIF(COMAND.EQ.'MODE') THEN
C        CALL MODE
C C
C C===============================================
C       ELSEIF(COMAND.EQ.'TIME') THEN
C ccc       CALL TIME
C C
C C===============================================
C       ELSE IF(COMAND.EQ.'LOAD') THEN
C C----- Read a new input geometry from input file
C        IF(COMARG.NE.' ') THEN
C         FILDEF = COMARG
C C
C        ELSE
C         CALL STRIP(FILDEF,LENF)
C         LENF1 = MAX(LENF,1)
C C
C         WRITE(*,2010) FILDEF(1:LENF1)
C  2010   FORMAT(' Enter input filename: ', A)
C         READ (*,1000)  FNNEW
C C
C         IF(FNNEW.EQ.' ') THEN
C          IF(LENF.EQ.0) GO TO 500
C         ELSE
C          FILDEF = FNNEW
C         ENDIF
C C
C        ENDIF
C C
C        CALL INPUT(LUINP,FILDEF,ERROR)
C        IF(ERROR) THEN
C         WRITE(*,*) 
C      &    '** File not processed. Current geometry may be corrupted.'
C        GO TO 500
C        ENDIF
C C
C        CALL PARSET
C C
C        IF(NRUN.EQ.0) THEN
C         CALL RUNINI
C        ELSE
C         WRITE(*,*)
C         WRITE(*,*) 'Existing run cases will be used.'
C         WRITE(*,*) 'Issue CASE or CINI command if necessary.'
C        ENDIF
C C
C C----- process geometry to define strip and vortex data
C        LPLTNEW = .TRUE.
C        CALL ENCALC
C C
C C----- initialize state
C       CALL VARINI
C C
C        LAIC = .FALSE.
C        LSRD = .FALSE.
C        LVEL = .FALSE.
C        LSOL = .FALSE.
C        LSEN = .FALSE.
C C
C C----- set up plotting parameters for new geometry 
C C      CALL PLPARS
C C
C C===============================================
C       ELSE IF(COMAND.EQ.'MASS') THEN
C C----- Read a new mass distribution file
C        IF(COMARG.NE.' ') THEN
C         FMSDEF = COMARG
C C
C        ELSE
C         CALL STRIP(FMSDEF,LENF)
C         LENF1 = MAX(LENF,1)
C C
C         WRITE(*,3010) FMSDEF(1:LENF1)
C  3010   FORMAT(' Enter mass filename: ', A)
C         READ (*,1000)  FNNEW
C C
C         IF(FNNEW.EQ.' ') THEN
C          IF(LENF.EQ.0) GO TO 500
C         ELSE
C          FMSDEF = FNNEW
C         ENDIF
C        ENDIF
C C
C        CALL STRIP(FMSDEF,NMS)
C        CALL MASGET(LUMAS,FMSDEF,ERROR)
C        IF(ERROR) THEN
C        ELSE
C         WRITE(*,*)
C         WRITE(*,*) 'Mass distribution read ...'
C         CALL MASSHO(6)
C C
C         CALL APPGET
C         WRITE(*,*) 
C      & '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
C         CALL APPSHO(6,RHO0)
C C
C         WRITE(*,*)
C         WRITE(*,*) 
C      &    'Use MSET to apply these mass,inertias to run cases'
C ccc        CALL MASPUT(1,NRMAX)
C        ENDIF
C C
C C===============================================
C       ELSE IF(COMAND.EQ.'CASE') THEN
C C----- Read a new run case file
C        IF(COMARG.NE.' ') THEN
C         FRNDEF = COMARG
C C
C        ELSE
C         CALL STRIP(FRNDEF,LENF)
C         LENF1 = MAX(LENF,1)
C C
C         WRITE(*,3020) FRNDEF(1:LENF1)
C  3020   FORMAT(' Enter run case filename: ', A)
C         READ (*,1000)  FNNEW
C C
C         IF(FNNEW.EQ.' ') THEN
C          IF(LENF.EQ.0) GO TO 500
C         ELSE
C          FRNDEF = FNNEW
C         ENDIF
C        ENDIF
C C
C        CALL STRIP(FRNDEF,NFR)
C        CALL RUNGET(LURUN,FRNDEF,ERROR)
C        IF(ERROR) THEN
C C        ELSE
C C         WRITE(*,1025) (IR, RTITLE(IR), IR=1, NRUN)
C        ENDIF
C C
C C----- initialize state
C       CALL VARINI
C C
C        LSOL = .FALSE.
C        LSEN = .FALSE.
C C
C C===============================================
C       ELSE IF(COMAND.EQ.'CINI') THEN
C        IF(LGEO) THEN
C         CALL RUNINI
C        ELSE
C         WRITE(*,*) 'No configuration available.'
C         NRUN = 0
C        ENDIF
C C
C===============================================
C       ELSE IF(COMAND.EQ.'MSET') THEN
C C----- set input mass,inertias
C        IF(NINPUT.GE.1) THEN
C         IR1 = IINPUT(1)
C        ELSE
C  60     WRITE(*,3060) 
C  3060   FORMAT(/
C      &     ' Enter index of target run case (0=all, -1=abort):  0')
C         IR1 = 0
C         CALL READI(1,IR1,ERROR)
C         IF(ERROR) GO TO 60
C        ENDIF
C C
C        IF(IR1.EQ.0) THEN
C         IR1 = 1
C         IR2 = NRUN
C        ELSE
C         IR2 = IR1
C        ENDIF
C C
C        IF(IR1.LT.1 .OR. IR1.GT.NRUN) GO TO 500
C C
C        CALL MASPUT(IR1,IR2)
C C
C        LSOL = .FALSE.
C        LSEN = .FALSE.
C C
C===============================================
C       ELSEIF(COMAND.EQ.'PLOP') THEN
C         CALL OPLSET(IDEV,IDEVH,IPSLU,LSVMOV,
C      &             SIZE,PLOTAR,
C      &             XMARG,YMARG,XPAGE,YPAGE,
C      &             CH,SCRNFRAC,LCURS,LCREV)
C C
C C===============================================
C       ELSEIF(COMAND.EQ.'NAME') THEN
C        IF(COMARG.EQ.' ') THEN
C         CALL ASKS('Enter new name^',TITLE)
C        ELSE
C         TITLE = COMARG
C        ENDIF
C C
C===============================================
C       ELSE
C        WRITE(*,1050) COMAND
C  1050  FORMAT(1X,A4,' command not recognized.  Type a "?" for list')
C C
C       ENDIF
C C
C       GO TO 500

      END 

      SUBROUTINE loadGEO(geom_file)
      INCLUDE 'AVL.INC'
      CHARACTER*256 geom_file
      LOGICAL ERROR
C     
      FILDEF = geom_file
      if (FILDEF .NE. '') then
       CALL INPUT(LUINP,FILDEF,ERROR)
       IF(ERROR) THEN
        WRITE(*,*) 
     &    '** File not processed. Current geometry may be corrupted.'
C        GO TO 500
       ENDIF
      endif
C
       CALL PARSET
C
       IF(NRUN.EQ.0) THEN
        CALL RUNINI
       ELSE
        WRITE(*,*)
        WRITE(*,*) 'Existing run cases will be used.'
        WRITE(*,*) 'Issue CASE or CINI command if necessary.'
       ENDIF
C
C----- process geometry to define strip and vortex data
       LPLTNEW = .TRUE.
       CALL ENCALC
C
C----- initialize state
C        CALL VARINI
C
       LAIC = .FALSE. ! Tell AVL that the AIC is no longer valid and to regenerate it
       LSRD = .FALSE. ! Tell AVL that unit source+doublet strengths are no longer valid and to regenerate them
       LVEL = .FALSE. ! Tell AVL that the induced velocity matrix is no longer valid and to regenerate it
       LSOL = .FALSE. ! Tell AVL that a valid solution no longer exists
       LSEN = .FALSE. ! Tell AVL that valid sensitives no longer exists
       
      END 

      SUBROUTINE loadMASS(mass_file)
      INCLUDE 'AVL.INC'
      CHARACTER*256 mass_file
      LOGICAL ERROR
C       LOGICAL LSOL      
C       LOGICAL LSEN
C C
C        CALL MASGET(LUMAS,FMSDEF,ERROR)
C        IF(ERROR) THEN
C        ELSE
C         WRITE(*,*)
C         WRITE(*,*) 'Mass distribution read ...'
C C        CALL MASSHO(6)
C C
C         CALL APPGET
C C         WRITE(*,*) 
C C      & '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
C         CALL APPSHO(6,RHO0)
C C C
C C         apply to all run cases
C         IR1 = 1
C         IR2 = NRUN

C       ENDIF
C        CALL MASPUT(IR1,IR2)
C C
C        LSOL = .FALSE.
C        LSEN = .FALSE.
C C         
       if (LVERBOSE) then
            WRITE(*,*) 'Trying to read file: ', FMSDEF(1:NMS), '  ...'
       endif 

C        CALL STRIP(FMSDEF,NMS)
       FMSDEF = mass_file
       CALL MASGET(LUMAS,FMSDEF,ERROR)
       IF(ERROR) THEN
       ELSE

        CALL APPGET
        if (LVERBOSE) then
         WRITE(*,*)
         WRITE(*,*) 'Mass distribution read ...'
         CALL MASSHO(6)

        WRITE(*,*) 
     & '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
        CALL APPSHO(6,RHO0)
C
        WRITE(*,*)
        WRITE(*,*) 
     &    'Use MSET to apply these mass,inertias to run cases'
ccc        CALL MASPUT(1,NRMAX)
        endif 
       ENDIF

        IR1 = 0
        NRUN = 1
        
       IF(IR1.EQ.0) THEN
        IR1 = 1
        IR2 = NRUN
       ELSE
        IR2 = IR1
       ENDIF
C
       IF(IR1.LT.1 .OR. IR1.GT.NRUN) WRITE(*,*) "ERROR in avl.f"
C
       CALL MASPUT(IR1,IR2)
C
       LSOL = .FALSE.
       LSEN = .FALSE.
C

      END 

C       SUBROUTINE PLINIT
C C---- Initialize plotting variables
C C
C       INCLUDE 'AVL.INC'
C C       INCLUDE 'AVLPLT.INC'
C C
C       REAL ORG(3)
C C
C C---- Plotting flag
C       IDEV = 1   ! X11 window only
C c     IDEV = 2   ! B&W PostScript output file only (no color)
C c     IDEV = 3   ! both X11 and B&W PostScript file
C c     IDEV = 4   ! Color PostScript output file only 
C c     IDEV = 5   ! both X11 and Color PostScript file 
C C
C C---- Re-plotting flag (for hardcopy)
C c     IDEVH = 2    ! B&W PostScript
C       IDEVH = 4    ! Color PostScript
C C
C C---- Movie-plotting flag
C cc    IDEVH = 3    ! B&W PostScript
C       IDEVM = 5   ! both X11 and Color PostScript file 
C C
C       LSVMOV = .FALSE.   ! no movie PS output yet
C C
C C---- PostScript output logical unit and file specification
C ccc   IPSLU = -1  ! output to files plotNNN.ps on LU 80, with NNN = 001, 002, ...
C       IPSLU = 0   ! output to file  plot.ps    on LU 80   (default case)
C ccc   IPSLU = nnn ! output to file  plotNNN.ps on LU NNN
C C
C C---- screen fraction taken up by plot window upon opening
C       SCRNFRAC = 0.70    ! Landscape
C C     SCRNFRAC = -0.85   ! Portrait  specified if < 0
C C
C C---- Default plot size in inches
C C-    (Default plot window is 11.0 x 8.5)
C       SIZE = 9.0
C C
C C---- plot aspect ratio
C       PLOTAR = 0.75
C C
C C---- character width/SIZE
C       CH = 0.017
C C
C C      CALL PLINITIALIZE
C C
C       NCOLORS = 0
C C---- set up color spectrum
C ccc      NCOLORS = 32
C ccc      CALL COLORSPECTRUMHUES(NCOLORS,'RYGCBM')
C C
C C---- plot-window dimensions in inches for plot blowup calculations
C C-    currently,  11.0 x 8.5  default window is hard-wired in libPlt
C       XPAGE = 11.0
C       YPAGE = 8.5
C C
C       XWIND = 11.0
C       YWIND = 8.5
C C
C C---- page margins in inches
C       XMARG = 0.0
C       YMARG = 0.0
C C
C C---- bottom,left plot margin from edge
C       PMARG = 0.15
C C
C       IF(IDEV.EQ.0) THEN 
C         LPLOT = .FALSE.
C       ENDIF
C C

C C---- set colors for run cases
C C       DO IR = 1, NRMAX
C C         IRCOLOR(IR) = MOD(IR-1,8) + 3
C C       ENDDO
C C
C C---- set vectors for little axes
C       SLEN = 0.5
C       HLEN = 0.5
C C
C       RHEAD = HLEN * 0.25
C       NHEAD = NHAXIS
C C
C       ORG(1) = 0.
C       ORG(2) = 0.
C       ORG(3) = 0.
C       DO IAX = 1, 3
C         UAXDIR(1,IAX) = 0.
C         UAXDIR(2,IAX) = 0.
C         UAXDIR(3,IAX) = 0.
C         UAXDIR(IAX,IAX) = 1.0
C         CALL ARWSET(ORG,UAXDIR(1,IAX),SLEN,HLEN,RHEAD,NHEAD,
C      &                  UAXARW(1,1,1,IAX),NLINAX)
C       ENDDO
C C
C C---- initial phase, eigenvector scale, slo-mo scale (for mode plots)
C       EPHASE = 0.0
C       EIGENF = 1.0
C       SLOMOF = 1.0
C       TMOFAC = 1.0

C       RETURN
C       END ! PLINIT



C       SUBROUTINE PLPARS
C       INCLUDE 'AVL.INC'
C       INCLUDE 'AVLPLT.INC'
C C
C       IMARKSURF = 0
C       DO N = 1, NSURF
C         LPLTSURF(N) = .TRUE. 
C       END DO
C       DO N = 1, NBODY
C         LPLTBODY(N) = .TRUE. 
C       END DO
C C
C C---- Scaling factors for velocity and pressure
C       CPFAC = MIN(0.4*CREF,0.1*BREF)  / CREF
C       ENFAC = MIN(0.3*CREF,0.06*BREF) / CREF
C       HNFAC = MIN(CREF,0.5*BREF)      / CREF
C C
C C---- initialize observer position angles and perspective 1/distance
C       AZIMOB = -45.0
C       ELEVOB =  20.0
C       TILTOB =   0.
C       ROBINV = 0.
C C
C C---- slo-mo factor
C       SLOMOF = 1.0
C C
C C---- eigenmode animation integration time step
C       DTIMED = 0.025
C C
C C---- movie-dump frame time step
C       DTMOVIE = 0.05
C C
C C---- max length of movie
C       TMOVIE = 10.0
C C
C C...Flags 
C       LABEL_BODY = .FALSE.
C       LABEL_SURF = .FALSE.
C       LABEL_STRP = .FALSE.
C       LABEL_VRTX = .FALSE.
C       LWAKEPLT   = .FALSE.
C       LHINGEPLT  = .FALSE.
C       LLOADPLT   = .FALSE.
C       LCNTLPTS   = .FALSE.
C       LNRMLPLT   = .FALSE.
C       LAXESPLT   = .TRUE.
C       LRREFPLT   = .TRUE.
C       LCLPERPLT  = .TRUE.
C       LDWASHPLT  = .TRUE.
C       LLABSURF   = .FALSE.
C       LCAMBER    = .FALSE.
C       LCHORDLINE = .TRUE.
C       LBOUNDLEG  = .TRUE.
C C
C C---- Initially assume nothing hidden
C       LHID = .TRUE.
C C
C C---- Initially assume no reverse color output
C       LCREV = .FALSE.
C C
C C---- flags to plot parameter values above eigenmode map
C       DO IP = 1, IPTOT
C         LPPAR(IP) = .FALSE.
C       ENDDO

C       LPPAR(IPALFA) = .TRUE.
C       LPPAR(IPBETA) = .TRUE.
C c      LPPAR(IPROTX) = .TRUE.
C c      LPPAR(IPROTY) = .TRUE.
C c      LPPAR(IPROTZ) = .TRUE.
C       LPPAR(IPCL  ) = .TRUE.
C       LPPAR(IPCD0 ) = .TRUE.

C       LPPAR(IPPHI ) = .TRUE.
C c      LPPAR(IPTHE ) = .TRUE.
C c      LPPAR(IPPSI ) = .TRUE.

C c      LPPAR(IPMACH) = .TRUE.
C       LPPAR(IPVEE ) = .TRUE.
C       LPPAR(IPRHO ) = .TRUE.
C c      LPPAR(IPGEE ) = .TRUE.

C       LPPAR(IPRAD ) = .TRUE.
C c      LPPAR(IPFAC ) = .TRUE.

C       LPPAR(IPXCG ) = .TRUE.
C c      LPPAR(IPYCG ) = .TRUE.
C       LPPAR(IPZCG ) = .TRUE.

C       LPPAR(IPMASS) = .TRUE.
C c      LPPAR(IPIXX ) = .TRUE.
C c      LPPAR(IPIYY ) = .TRUE.
C c      LPPAR(IPIZZ ) = .TRUE.
C c      LPPAR(IPIXY ) = .TRUE.
C c      LPPAR(IPIYZ ) = .TRUE.
C c      LPPAR(IPIZX ) = .TRUE.

C c      LPPAR(IPCLA ) = .TRUE.
C c      LPPAR(IPCLU ) = .TRUE.
C c      LPPAR(IPCMA ) = .TRUE.
C c      LPPAR(IPCMU ) = .TRUE.

C       RETURN
C       END ! PLPARS



      SUBROUTINE VARINI
      INCLUDE 'AVL.INC'
C
C---- initialize state
      ALFA = 0.
      BETA = 0.
      WROT(1) = 0.
      WROT(2) = 0.
      WROT(3) = 0.
C
      DO N = 1, NCONTROL
        DELCON(N) = 0.0
      ENDDO
C
      DO N = 1, NDESIGN
        DELDES(N) = 0.0
      ENDDO
      LSOL = .FALSE.
C
      RETURN
      END ! VARINI


      SUBROUTINE DEFINI
      INCLUDE 'AVL.INC'
      
      LVERBOSE = .FALSE.
      LTIMING = .FALSE.
C
C---- flag for forces in standard NASA stability axes (as in Etkin)
      LNASA_SA  = .TRUE.
C
C---- flag for rotations defined in stability axes or body axes
      LSA_RATES = .TRUE.
C
      LPTOT   = .TRUE.
      LPSURF  = .FALSE.
      LPSTRP  = .FALSE.
      LPELE   = .FALSE.
      LPHINGE = .FALSE.
      LPDERIV = .FALSE.
C
      LGEO  = .FALSE.
      LENC  = .FALSE.
C
      LAIC  = .FALSE.
      LSRD  = .FALSE.
      LVEL  = .FALSE.
      LSOL  = .FALSE.
      LSEN  = .FALSE.
C
      LVISC    = .TRUE.
      LBFORCE  = .TRUE.
      LTRFORCE = .TRUE.
C
      LMWAIT = .FALSE.
C
      MATSYM = 0
      NITMAX = 20
C
      SAXFR = 0.25  ! x/c location of spanwise axis for Vperp definition
C
      VRCORE = 0.25   ! vortex core radius / vortex span
      SRCORE = 0.75   ! source core radius / body radius
C
C---- dafault basic units
      UNITL = 1.
      UNITM = 1.
      UNITT = 1.
      UNCHL = 'Lunit'
      UNCHM = 'Munit'
      UNCHT = 'Tunit'
      NUL = 5
      NUM = 5
      NUT = 5
C
C---- set corresponding derived units
      CALL UNITSET
C
C---- default air density and grav. accel.
      RHO0 = 1.0
      GEE0 = 1.0
C
C---- no eigenvalue reference data yet
      FEVDEF = ' '
      DO IR = 1, NRMAX
        NEIGENDAT(IR) = 0
      ENDDO
C
C---- no run cases defined yet
      NRUN = 0
      IRUN = 1
C
C---- number of valid time levels stored
      NTLEV = 0
C
C---- default time step, and number of time steps to take
      DELTAT = 0.0
      NTSTEPS = 0
C
      RETURN
      END ! DEFINI



      SUBROUTINE PARSET
      INCLUDE 'AVL.INC'
C
C---- variable names
      VARNAM(IVALFA) = 'alpha '
      VARNAM(IVBETA) = 'beta  '
      VARNAM(IVROTX) = 'pb/2V '
      VARNAM(IVROTY) = 'qc/2V '
      VARNAM(IVROTZ) = 'rb/2V '
C
C---- variable selection keys
      VARKEY(IVALFA) = 'A lpha'
      VARKEY(IVBETA) = 'B eta'
      VARKEY(IVROTX) = 'R oll  rate'
      VARKEY(IVROTY) = 'P itch rate'
      VARKEY(IVROTZ) = 'Y aw   rate'
C
C---- constraint names
CCC                     123456789012
      CONNAM(ICALFA) = 'alpha '
      CONNAM(ICBETA) = 'beta  '
      CONNAM(ICROTX) = 'pb/2V '
      CONNAM(ICROTY) = 'qc/2V '
      CONNAM(ICROTZ) = 'rb/2V '
      CONNAM(ICCL  ) = 'CL    '
      CONNAM(ICCY  ) = 'CY    '
      CONNAM(ICMOMX) = 'Cl roll mom'
      CONNAM(ICMOMY) = 'Cm pitchmom'
      CONNAM(ICMOMZ) = 'Cn yaw  mom'
C
C---- constraint selection keys
      CONKEY(ICALFA) = 'A '
      CONKEY(ICBETA) = 'B '
      CONKEY(ICROTX) = 'R '
      CONKEY(ICROTY) = 'P '
      CONKEY(ICROTZ) = 'Y '
      CONKEY(ICCL  ) = 'C '
      CONKEY(ICCY  ) = 'S '
      CONKEY(ICMOMX) = 'RM'
      CONKEY(ICMOMY) = 'PM'
      CONKEY(ICMOMZ) = 'YM'
C
C------------------------------------------------------------------------
      IZERO = ICHAR('0')
C
C---- add control variables, direct constraints
      DO N = 1, NCONTROL
        ITEN = N/10
        IONE = N - 10*ITEN
C
C------ assign slots in variable ond constraint lists
        IV = IVTOT + N
        IC = ICTOT + N
        VARNAM(IV) = DNAME(N)
        CONNAM(IC) = DNAME(N)
        IF(ITEN.EQ.0) THEN
         VARKEY(IV) = 'D' // CHAR(IZERO+IONE) // ' '
     &             // ' ' // DNAME(N)(1:8)
         CONKEY(IC) = 'D' // CHAR(IZERO+IONE)
        ELSE
         VARKEY(IV) = 'D' // CHAR(IZERO+ITEN) // CHAR(IZERO+IONE)
     &             // ' ' // DNAME(N)(1:8)
         CONKEY(IC) = 'D' // CHAR(IZERO+ITEN) // CHAR(IZERO+IONE)
        ENDIF
C
        LCONDEF(N) = .TRUE.
      ENDDO
C
C---- default design-variable flags, names
      DO N = 1, NDESIGN
        LDESDEF(N) = .TRUE.
      ENDDO
C
C---- total number of variables, constraints
      NVTOT = IVTOT + NCONTROL
      NCTOT = ICTOT + NCONTROL
C
C---- run-case parameter names
      PARNAM(IPALFA) = 'alpha    '
      PARNAM(IPBETA) = 'beta     '
      PARNAM(IPROTX) = 'pb/2V    '
      PARNAM(IPROTY) = 'qc/2V    '
      PARNAM(IPROTZ) = 'rb/2V    '
      PARNAM(IPCL )  = 'CL       '
      PARNAM(IPCD0)  = 'CDo      '
      PARNAM(IPPHI)  = 'bank     '
      PARNAM(IPTHE)  = 'elevation'
      PARNAM(IPPSI)  = 'heading  '
      PARNAM(IPMACH) = 'Mach     '
      PARNAM(IPVEE)  = 'velocity '
      PARNAM(IPRHO)  = 'density  '
      PARNAM(IPGEE)  = 'grav.acc.'
      PARNAM(IPRAD)  = 'turn_rad.'
      PARNAM(IPFAC)  = 'load_fac.'
      PARNAM(IPXCG)  = 'X_cg     '
      PARNAM(IPYCG)  = 'Y_cg     '
      PARNAM(IPZCG)  = 'Z_cg     '
      PARNAM(IPMASS) = 'mass     '
      PARNAM(IPIXX)  = 'Ixx      '
      PARNAM(IPIYY)  = 'Iyy      '
      PARNAM(IPIZZ)  = 'Izz      '
      PARNAM(IPIXY)  = 'Ixy      '
      PARNAM(IPIYZ)  = 'Iyz      '
      PARNAM(IPIZX)  = 'Izx      '
      PARNAM(IPCLA)  = 'visc CL_a'
      PARNAM(IPCLU)  = 'visc CL_u'
      PARNAM(IPCMA)  = 'visc CM_a'
      PARNAM(IPCMU)  = 'visc CM_u'
C
C---- total number of parameters
      NPTOT = IPTOT
C
C---- set default parameter unit names
      CALL PARNSET
C
      RETURN
      END ! PARSET



      SUBROUTINE PARNSET
      INCLUDE 'AVL.INC'
C
C---- set parameter unit name
      DO IP = 1, IPTOT
        PARUNCH(IP) = ' '
      ENDDO
C
      PARUNCH(IPALFA) = 'deg'
      PARUNCH(IPBETA) = 'deg'
      PARUNCH(IPPHI)  = 'deg'
      PARUNCH(IPTHE)  = 'deg'
      PARUNCH(IPPSI)  = 'deg'
      PARUNCH(IPVEE)  = UNCHV
      PARUNCH(IPRHO)  = UNCHD
      PARUNCH(IPGEE)  = UNCHA
      PARUNCH(IPRAD)  = UNCHL

C---- bug  21 Feb 13   MD
c      PARUNCH(IPXCG)  = UNCHL
c      PARUNCH(IPYCG)  = UNCHL
c      PARUNCH(IPZCG)  = UNCHL
      PARUNCH(IPXCG)  = 'Lunit'
      PARUNCH(IPYCG)  = 'Lunit'
      PARUNCH(IPZCG)  = 'Lunit'

      PARUNCH(IPMASS) = UNCHM
      PARUNCH(IPIXX)  = UNCHI
      PARUNCH(IPIYY)  = UNCHI
      PARUNCH(IPIZZ)  = UNCHI
      PARUNCH(IPIXY)  = UNCHI
      PARUNCH(IPIYZ)  = UNCHI
      PARUNCH(IPIZX)  = UNCHI
C
      RETURN
      END ! PARNSET




      SUBROUTINE RUNINI
      INCLUDE 'AVL.INC'
C
      if (lverbose) then 
            WRITE(*,*)
            WRITE(*,*) 'Initializing run cases...'
      end if
C
C---- go over all run cases
      DO IR = 1, NRMAX
C------ index of default constraint for each variable
        ICON(IVALFA,IR) = ICALFA
        ICON(IVBETA,IR) = ICBETA
        ICON(IVROTX,IR) = ICROTX
        ICON(IVROTY,IR) = ICROTY
        ICON(IVROTZ,IR) = ICROTZ
C
C------ default constraint values
        DO IC = 1, ICTOT
          CONVAL(IC,IR) = 0.
        ENDDO
C
C------ default run case titles
        RTITLE(IR) = ' -unnamed- '
C
C------ default dimensional run case parameters
        DO IP = 1, NPTOT
          PARVAL(IP,IR) = 0.
        ENDDO
        PARVAL(IPGEE,IR) = GEE0
        PARVAL(IPRHO,IR) = RHO0
        PARVAL(IPMACH,IR) = MACH0
C
C------ default CG location is the input reference location
        PARVAL(IPXCG,IR) = XYZREF0(1)
        PARVAL(IPYCG,IR) = XYZREF0(2)
        PARVAL(IPZCG,IR) = XYZREF0(3)
C       just set the initial values         
        XYZREF(1) = XYZREF0(1)
        XYZREF(2) = XYZREF0(2)
        XYZREF(3) = XYZREF0(3)
C
        PARVAL(IPMASS,IR) = RMASS0
        PARVAL(IPIXX,IR) = RINER0(1,1)
        PARVAL(IPIYY,IR) = RINER0(2,2)
        PARVAL(IPIZZ,IR) = RINER0(3,3)
        PARVAL(IPIXY,IR) = RINER0(1,2)
        PARVAL(IPIYZ,IR) = RINER0(2,3)
        PARVAL(IPIZX,IR) = RINER0(3,1)
C
        PARVAL(IPCD0,IR) = CDREF0
C
        PARVAL(IPCLA,IR) = DCL_A0
        PARVAL(IPCLU,IR) = DCL_U0
        PARVAL(IPCMA,IR) = DCM_A0
        PARVAL(IPCMU,IR) = DCM_U0
C
        ITRIM(IR) = 0
        NEIGEN(IR) = 0
      ENDDO
C
C---- add control variables, direct constraints
      DO N = 1, NDMAX
        IV = IVTOT + N
        IC = ICTOT + N
        DO IR = 1, NRMAX
          ICON(IV,IR) = IC
          CONVAL(IC,IR) = 0.
        ENDDO
      ENDDO
C
C---- default number of run cases
      IRUN = 1
      NRUN = 1
C
C---- all run cases are targets for eigenmode calculation
      IRUNE = 0
C
C---- first run case is default for time march initial state
      IRUNT = 1
C
      RETURN
      END ! RUNINI



      SUBROUTINE RUNGET(LU,FNAME,ERROR)
C-------------------------------------------------
C     Reads run case file into run case arrays
C-------------------------------------------------
      INCLUDE 'AVL.INC'
      CHARACTER*(*) FNAME
      LOGICAL ERROR
C
      CHARACTER*120 LINE, REST
      CHARACTER*12 VARN, CONN
      CHARACTER*8  PARN
C
      OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=90)
      ILINE = 0
C
      IR = 0
C
C==============================================================
C---- start line-reading loop
 10   CONTINUE
C
      READ(LU,1000,END=50) LINE
 1000 FORMAT(A)
      ILINE = ILINE + 1
C
      KCOL = INDEX(LINE,':' )
      KARR = INDEX(LINE,'->')
      KEQU = INDEX(LINE,'=' )
      IF(KCOL.NE.0) THEN
C----- start of new run case
       READ(LINE(KCOL-3:KCOL-1),*,ERR=80) IR
C
       IF(IR.LT.1 .OR. IR.GT.NRMAX) THEN
        WRITE(*,*) 'RUNGET:  Run case array limit NRMAX exceeded:', IR
        IR = 0
        GO TO 10
       ENDIF
C
       NRUN = MAX(NRUN,IR)
C
       RTITLE(IR) = LINE(KCOL+1:80)
       CALL STRIP(RTITLE(IR),NRT)
C
      ELSEIF(IR.EQ.0) THEN
C----- keep ignoring lines if valid run case index is not set
       GO TO 10
C
      ELSEIF(KARR.NE.0 .AND. KEQU.NE.0) THEN
C----- variable/constraint declaration line
       VARN = LINE(1:KARR-1)
       CONN = LINE(KARR+2:KEQU-1)
       CALL STRIP(VARN,NVARN)
       CALL STRIP(CONN,NCONN)
C
       DO IV = 1, NVTOT
         IF(INDEX(VARNAM(IV),VARN(1:NVARN)).NE.0) GO TO 20
       ENDDO
       WRITE(*,*) 'Ignoring unrecognized variable: ', VARN(1:NVARN)
       GO TO 10
C
 20    CONTINUE
       DO IC = 1, NCTOT
         IF(INDEX(CONNAM(IC),CONN(1:NCONN)).NE.0) GO TO 25
       ENDDO
       WRITE(*,*) 'Ignoring unrecognized constraint: ', CONN(1:NCONN)
       GO TO 10
C
 25    CONTINUE
       READ(LINE(KEQU+1:80),*,ERR=80) CONV
C
       ICON(IV,IR) = IC
       CONVAL(IC,IR) = CONV
C
      ELSEIF(KARR.EQ.0 .AND. KEQU.NE.0) THEN
C----- run case parameter data line
       PARN = LINE(1:KEQU-1)
       CALL STRIP(PARN,NPARN)
       DO IP = 1, NPTOT
         IF(INDEX(PARNAM(IP),PARN(1:NPARN)).NE.0) GO TO 30
       ENDDO
       WRITE(*,*) 'Ignoring unrecognized parameter: ', PARN(1:NPARN)
       GO TO 10
C
 30    CONTINUE
       REST = LINE(KEQU+1:80)
       READ(REST,*,ERR=80) PARV
       PARVAL(IP,IR) = PARV


       if(.false.) then
       CALL STRIP(REST,NREST)
       KBLK = INDEX(REST,' ')
       IF(KBLK .NE. 0) THEN
        REST = REST(KBLK+1:80)
        CALL STRIP(REST,NREST)
        IF(NREST.GT.0) THEN
         PARUNCH(IP) = REST
        ENDIF
       ENDIF
       endif

      ENDIF
C
C---- keep reading lines
      GO TO 10
C
C==============================================================
C
 50   CONTINUE
      CLOSE(LU)
      ERROR = .FALSE.
      RETURN
C
 80   CONTINUE
      CALL STRIP(FNAME,NFN)
      CALL STRIP(LINE ,NLI)
      WRITE(*,8000) FNAME(1:NFN), ILINE, LINE(1:NLI)
 8000 FORMAT(/' Run case file  ',A,'  read error on line', I4,':',A)
      CLOSE(LU)
      ERROR = .TRUE.
      NRUN = 0
      RETURN
C
 90   CONTINUE
      CALL STRIP(FNAME,NFN)
      WRITE(*,9000) FNAME(1:NFN)
 9000 FORMAT(/' Run case file  ',A,'  open error')
      ERROR = .TRUE.
      RETURN
      END ! RUNGET



      SUBROUTINE RUNSAV(LU)
      INCLUDE 'AVL.INC'
C
      DO IR = 1, NRUN
        WRITE(LU,1010) IR, RTITLE(IR)
        DO IV = 1, NVTOT
          IC = ICON(IV,IR)
          WRITE(LU,1050) VARNAM(IV), CONNAM(IC), CONVAL(IC,IR)
        ENDDO
C
        WRITE(LU,*)
C
        DO IP = 1, NPTOT
          WRITE(LU,1080) PARNAM(IP), PARVAL(IP,IR), PARUNCH(IP)
        ENDDO
      ENDDO
C
 1010 FORMAT(/' ---------------------------------------------'
     &       /' Run case', I3,':  ', A /)
 1050 FORMAT(1X,A,' ->  ', A, '=', G14.6, 1X, A)
 1080 FORMAT(1X,A,'=', G14.6, 1X, A)
C
      RETURN
      END ! RUNSAV



      LOGICAL FUNCTION LOWRIT(FNAME)
      CHARACTER*(*) FNAME
C
      CHARACTER*1 ANS
 1000 FORMAT(A)
C
      K = INDEX(FNAME,' ')
C
      WRITE(*,*) 'File  ', FNAME(1:K), ' exists.  Overwrite?  Y'
      READ (*,1000) ANS
      LOWRIT = INDEX('Nn',ANS) .EQ. 0
C
      RETURN
      END


      SUBROUTINE AOCFIL(FNAME,IFILE)
      CHARACTER*(*) FNAME
C
      CHARACTER*1 ANS
 1000 FORMAT(A)
C
      K = INDEX(FNAME,' ')
C
      WRITE(*,*) 'File  ', FNAME(1:K), 
     &     ' exists.  Append, Overwrite, or Cancel?  A'
      READ (*,1000) ANS
      IFILE = INDEX('AOC',ANS) + INDEX('aoc',ANS)
C
      IF(IFILE.EQ.0) IFILE = 1
C
      RETURN
      END

