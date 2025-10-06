C***********************************************************************
C    Module:  amake.f
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

      SUBROUTINE MAKESURF(ISURF)
C--------------------------------------------------------------
C     Sets up all stuff for surface ISURF, 
C     using info from configuration input file.
C--------------------------------------------------------------
      INCLUDE 'AVL.INC'
C
C
      REAL XYZLEL(3), XYZLER(3)
C
      PARAMETER (KCMAX=50,
     &           KSMAX=500)
      REAL XPT0(KCMAX), XCP0(KCMAX), XVR0(KCMAX), XSR0(KCMAX),
     &     XPT1(KCMAX), XCP1(KCMAX), XVR1(KCMAX), XSR1(KCMAX),
     &     XPT2(KCMAX), XCP2(KCMAX), XVR2(KCMAX), XSR2(KCMAX)
      REAL XPT(KCMAX), XCP(KCMAX), XVR(KCMAX), XSR(KCMAX),
     &     YPT(KSMAX), YCP(KSMAX)
      REAL YZLEN(KSMAX)
      INTEGER IPTLOC(KSMAX)
C
      PARAMETER (KPMAX=2*KCMAX+2*KSMAX)
      REAL FSPACE(KPMAX)
C
      REAL CHSINL_G(NGMAX),CHCOSL_G(NGMAX),
     &     CHSINR_G(NGMAX),CHCOSR_G(NGMAX)
      INTEGER ISCONL(NDMAX), ISCONR(NDMAX)
      REAL XLED(NDMAX), XTED(NDMAX), GAINDA(NDMAX)
      integer idx_vor, idx_strip
C
C
      IF(NSEC(ISURF).LT.2) THEN
       WRITE(*,*) '*** Need at least 2 sections per surface.'
       STOP
      ENDIF
C
C
      IF(NVC(ISURF).GT.KCMAX) THEN
       WRITE(*,*) '* MAKESURF: Array overflow.  Increase KCMAX to',
     &                                                   NVC(ISURF)
       NVC(ISURF) = KCMAX
      ENDIF
C
      IF(NVS(ISURF).GT.KSMAX) THEN
       WRITE(*,*) '* MAKESURF: Array overflow.  Increase KSMAX to', 
     &                                                   NVS(ISURF)
       NVS(ISURF) = KSMAX
      ENDIF
C
C--- Image flag set to indicate section definition direction
C    IMAGS= 1  defines edge 1 located at surface root edge 
C    IMAGS=-1  defines edge 2 located at surface root edge (reflected surfaces)
      IMAGS(ISURF) = 1
      
      if (ISURF == 1) then
            IFRST(ISURF) = 1
      else
            IFRST(ISURF) =  IFRST(ISURF-1) +  NK(ISURF-1)*NJ(ISURF-1)        
      endif
      ! write(*,*) 'IFRST(ISURF)', IFRST(ISURF)
      ! IFRST(ISURF) = NVOR   + 1 
      ! write(*,*) 'IFRST(ISURF) 2', IFRST(ISURF)
      
      
      ! JFRST(ISURF) = NSTRIP + 1
      if (ISURF == 1) then
            JFRST(ISURF) = 1
      else
            JFRST(ISURF) =  JFRST(ISURF-1) +  NJ(ISURF-1)
      endif
      
      NK(ISURF) = NVC(ISURF)
      idx_strip = JFRST(ISURF)
C
C-----------------------------------------------------------------
C---- Arc length positions of sections in wing trace in y-z plane
      YZLEN(1) = 0.
      DO ISEC = 2, NSEC(ISURF)
        DY = XYZLES(2,ISEC, ISURF) - XYZLES(2,ISEC-1, ISURF)
        DZ = XYZLES(3,ISEC, ISURF) - XYZLES(3,ISEC-1, ISURF)
        YZLEN(ISEC) = YZLEN(ISEC-1) + SQRT(DY*DY + DZ*DZ)
      ENDDO
C
      ! we can not rely on the original condition becuase NVS(ISURF) is filled 
      ! and we may want to rebuild the surface later
      IF((NVS(ISURF).EQ.0) .or.
     &   (LSURFSPACING(ISURF) .EQV. .FALSE.)) THEN
C----- set spanwise spacing using spacing parameters for each section interval
       DO ISEC = 1, NSEC(ISURF)-1
         NVS(ISURF) = NVS(ISURF) + NSPANS(ISEC, ISURF)
       ENDDO
       IF(NVS(ISURF).GT.KSMAX) THEN
           WRITE(*,*) '*** MAKESURF: Array overflow. Increase',
     &      'KSMAX to',NVS(ISURF)
        STOP
       ENDIF
C
       NVS(ISURF) = 0
       YPT(1) = YZLEN(1)
       IPTLOC(1) = 1
C
       DO ISEC = 1, NSEC(ISURF)-1
         DYZLEN = YZLEN(ISEC+1) - YZLEN(ISEC)
C
         NVINT = NSPANS(ISEC, ISURF)
C
C------- set spanwise spacing array
         NSPACE = 2*NVINT + 1
         IF(NSPACE.GT.KPMAX) THEN
          WRITE(*,*) '*** MAKESURF: Array overflow. Increase KPMAX to', 
     &                 NSPACE
          STOP
         ENDIF
         CALL SPACER(NSPACE,SSPACES(ISEC, ISURF),FSPACE)
C
         DO N = 1, NVINT
           IVS = NVS(ISURF) + N
           YCP(IVS)   = YPT(NVS(ISURF)+1) + DYZLEN*FSPACE(2*N)
           YPT(IVS+1) = YPT(NVS(ISURF)+1) + DYZLEN*FSPACE(2*N+1)
         ENDDO
         IPTLOC(ISEC+1) = NVS(ISURF) + NVINT + 1
C
         NVS(ISURF) = NVS(ISURF) + NVINT
       ENDDO
C
      ELSE
C
C----- Otherwise, set spanwise spacing using the SURFACE spanwise
C      parameters NVS, SSPACE
C
C      This spanwise spacing is modified (fudged) to align vortex edges
C      with SECTIONs as defined.  This allows CONTROLs to be defined
C      without bridging vortex strips
C
       NSPACE = 2*NVS(ISURF) + 1
       IF(NSPACE.GT.KPMAX) THEN
        WRITE(*,*) '*** MAKESURF: Array overflow. Increase KPMAX to', 
     &              NSPACE
        STOP
       ENDIF
       CALL SPACER(NSPACE,SSPACE(ISURF),FSPACE)
C
       YPT(1) = YZLEN(1)
       DO IVS = 1, NVS(ISURF)
         YCP(IVS)   = YZLEN(1) + (YZLEN(NSEC(ISURF))
     &                         -  YZLEN(1))*FSPACE(2*IVS)
         YPT(IVS+1) = YZLEN(1) + (YZLEN(NSEC(ISURF))
     &                         -  YZLEN(1))*FSPACE(2*IVS+1)
       ENDDO
C
       NPT = NVS(ISURF) + 1
C
C----- find node nearest each section
       DO ISEC = 2, NSEC(ISURF)-1
         YPTLOC = 1.0E9
         IPTLOC(ISEC) = 1
         DO IPT = 1, NPT
           YPTDEL = ABS(YZLEN(ISEC) - YPT(IPT))
           IF(YPTDEL .LT. YPTLOC) THEN
            YPTLOC = YPTDEL
            IPTLOC(ISEC) = IPT
           ENDIF
         ENDDO
       ENDDO
       IPTLOC(1)    = 1
       IPTLOC(NSEC(ISURF)) = NPT
C
C----- fudge spacing array to make nodes match up exactly with interior sections
       DO ISEC = 2, NSEC(ISURF)-1
         IPT1 = IPTLOC(ISEC-1)
         IPT2 = IPTLOC(ISEC  )
         IF(IPT1.EQ.IPT2) THEN
          CALL STRIP(STITLE(ISURF),NST)
          WRITE(*,7000) ISEC, STITLE(ISURF)(1:NST)
          STOP
         ENDIF
C
C----- fudge spacing to this section so that nodes match up exactly with section
         YPT1 = YPT(IPT1)
         YSCALE = (YZLEN(ISEC)-YZLEN(ISEC-1)) / (YPT(IPT2)-YPT(IPT1))
         DO IPT = IPT1, IPT2-1
           YPT(IPT) = YZLEN(ISEC-1) + YSCALE*(YPT(IPT)-YPT1)
         ENDDO
         DO IVS = IPT1, IPT2-1
           YCP(IVS) = YZLEN(ISEC-1) + YSCALE*(YCP(IVS)-YPT1)
         ENDDO
C
C----- check for unique spacing node for next section, if not we need more nodes
         IPT1 = IPTLOC(ISEC  )
         IPT2 = IPTLOC(ISEC+1)
         IF(IPT1.EQ.IPT2) THEN
          CALL STRIP(STITLE(ISURF),NST)
          WRITE(*,7000) ISEC, STITLE(ISURF)(1:NST)
          STOP
         ENDIF
C
C----- fudge spacing to this section so that nodes match up exactly with section
         YPT1 = YPT(IPT1)
         YSCALE = (YPT(IPT2)-YZLEN(ISEC)) / (YPT(IPT2)-YPT(IPT1))
         DO IPT = IPT1, IPT2-1
           YPT(IPT) = YZLEN(ISEC) + YSCALE*(YPT(IPT)-YPT1)
         ENDDO
         DO IVS = IPT1, IPT2-1
           YCP(IVS) = YZLEN(ISEC) + YSCALE*(YCP(IVS)-YPT1)
         ENDDO
C
 7000    FORMAT(
     &   /' *** Cannot adjust spanwise spacing at section', I3, 
     &    ', on surface ', A
     &   /' *** Insufficient number of spanwise vortices to work with')
       ENDDO
C
      ENDIF
cc#ifdef USE_CPOML
C...  store section counters
      IF (ISURF .EQ. 1) THEN
        ICNTFRST(ISURF) = 1
      ELSE
        ICNTFRST(ISURF) = ICNTFRST(ISURF-1) + NCNTSEC(ISURF-1)
      ENDIF
      NCNTSEC(ISURF) = NSEC(ISURF)
      DO ISEC = 1, NSEC(ISURF)
        II = ICNTFRST(ISURF) + (ISEC-1)
        ICNTSEC(II) = IPTLOC(ISEC)
      ENDDO
cc#endif
C
C
C====================================================
C---- define strips between input sections
C
      NJ(ISURF) = 0
C
      IF(NCONTROL.GT.NDMAX) THEN
       WRITE(*,*) '*** Too many control variables.  Increase NDMAX to',
     &            NCONTROL
       STOP
      ENDIF
C
      IF(NDESIGN.GT.NGMAX) THEN
       WRITE(*,*) '*** Too many design variables.  Increase NGMAX to',
     &            NDESIGN
       STOP
      ENDIF
C
C---- go over section intervals
      DO 200 ISEC = 1, NSEC(ISURF)-1
        XYZLEL(1) = XYZSCAL(1,ISURF)*XYZLES(1,ISEC,ISURF)   
     &              + XYZTRAN(1,ISURF)
        XYZLEL(2) = XYZSCAL(2,ISURF)*XYZLES(2,ISEC,ISURF)   
     &              + XYZTRAN(2,ISURF)
        XYZLEL(3) = XYZSCAL(3,ISURF)*XYZLES(3,ISEC,ISURF)   
     &              + XYZTRAN(3,ISURF)
        XYZLER(1) = XYZSCAL(1,ISURF)*XYZLES(1,ISEC+1,ISURF) 
     &              + XYZTRAN(1,ISURF)
        XYZLER(2) = XYZSCAL(2,ISURF)*XYZLES(2,ISEC+1,ISURF) 
     &              + XYZTRAN(2,ISURF)
        XYZLER(3) = XYZSCAL(3,ISURF)*XYZLES(3,ISEC+1,ISURF) 
     &              + XYZTRAN(3,ISURF)
C
        WIDTH = SQRT(  (XYZLER(2)-XYZLEL(2))**2
     &               + (XYZLER(3)-XYZLEL(3))**2 )
C
        CHORDL = XYZSCAL(1,ISURF)*CHORDS(ISEC, ISURF)
        CHORDR = XYZSCAL(1,ISURF)*CHORDS(ISEC+1, ISURF)
C
        CLAFL = CLAF(ISEC,  ISURF)
        CLAFR = CLAF(ISEC+1,ISURF)
C
C------ removed CLAF influence on zero-lift angle  (MD  21 Mar 08)
        AINCL = AINCS(ISEC  , ISURF)*DTR + ADDINC(ISURF)*DTR
        AINCR = AINCS(ISEC+1, ISURF)*DTR + ADDINC(ISURF)*DTR
cc      AINCL = AINCS(ISEC)   + ADDINC(ISURF) - 4.0*DTR*(CLAFL-1.0)
cc      AINCR = AINCS(ISEC+1) + ADDINC(ISURF) - 4.0*DTR*(CLAFR-1.0)
C
        CHSINL = CHORDL*SIN(AINCL)
        CHSINR = CHORDR*SIN(AINCR)
        CHCOSL = CHORDL*COS(AINCL)
        CHCOSR = CHORDR*COS(AINCR)
C
C------ set control-declaration lines for each control variable
        DO N = 1, NCONTROL
          ISCONL(N) = 0
          ISCONR(N) = 0
          DO ISCON = 1, NSCON(ISEC,ISURF)
            IF(ICONTD(ISCON,ISEC,ISURF)  .EQ.N) ISCONL(N) = ISCON
          ENDDO
          DO ISCON = 1, NSCON(ISEC+1,ISURF)
            IF(ICONTD(ISCON,ISEC+1,ISURF).EQ.N) ISCONR(N) = ISCON
          ENDDO
        ENDDO
C
C------ set design-variable sensitivities of CHSIN and CHCOS
        DO N = 1, NDESIGN
          CHSINL_G(N) = 0.
          CHSINR_G(N) = 0.
          CHCOSL_G(N) = 0.
          CHCOSR_G(N) = 0.
C
          DO ISDES = 1, NSDES(ISEC,ISURF)
            IF(IDESTD(ISDES,ISEC,ISURF).EQ.N) THEN
             CHSINL_G(N) =  CHCOSL * GAING(ISDES,ISEC,ISURF)*DTR
             CHCOSL_G(N) = -CHSINL * GAING(ISDES,ISEC,ISURF)*DTR
            ENDIF
          ENDDO
C
          DO ISDES = 1, NSDES(ISEC+1,ISURF)
            IF(IDESTD(ISDES,ISEC+1,ISURF).EQ.N) THEN
             CHSINR_G(N) =  CHCOSR * GAING(ISDES,ISEC+1,ISURF)*DTR
             CHCOSR_G(N) = -CHSINR * GAING(ISDES,ISEC+1,ISURF)*DTR
            ENDIF
          ENDDO
        ENDDO
C
C
C------ go over chord strips
        IPTL = IPTLOC(ISEC)
        IPTR = IPTLOC(ISEC+1)
        NSPAN = IPTR - IPTL       
        NJ(ISURF) = NJ(ISURF) +  NSPAN

        DO 150 ISPAN = 1, NSPAN
C-------- define left and right edges of vortex strip
C-          note that incidence angle is set by ATAN of chord projections,
C-          not by linear interpolation of AINC
          IPT1 = IPTL + ISPAN - 1
          IPT2 = IPTL + ISPAN
          IVS  = IPTL + ISPAN - 1
          F1 = (YPT(IPT1)-YPT(IPTL))/(YPT(IPTR)-YPT(IPTL))
          F2 = (YPT(IPT2)-YPT(IPTL))/(YPT(IPTR)-YPT(IPTL))
          FC = (YCP(IVS) -YPT(IPTL))/(YPT(IPTR)-YPT(IPTL))
C
C-------- store strip in global data arrays
      !     NSTRIP = NSTRIP + 1
      !     NJ(ISURF) = NJ(ISURF) + 1
C
          RLE1(1,idx_strip) = (1.0-F1)*XYZLEL(1) + F1*XYZLER(1)
          RLE1(2,idx_strip) = (1.0-F1)*XYZLEL(2) + F1*XYZLER(2)
          RLE1(3,idx_strip) = (1.0-F1)*XYZLEL(3) + F1*XYZLER(3)
          CHORD1(idx_strip) = (1.0-F1)*CHORDL    + F1*CHORDR
C
          RLE2(1,idx_strip) = (1.0-F2)*XYZLEL(1) + F2*XYZLER(1)
          RLE2(2,idx_strip) = (1.0-F2)*XYZLEL(2) + F2*XYZLER(2)
          RLE2(3,idx_strip) = (1.0-F2)*XYZLEL(3) + F2*XYZLER(3)
          CHORD2(idx_strip) = (1.0-F2)*CHORDL    + F2*CHORDR
C
          RLE(1,idx_strip)  = (1.0-FC)*XYZLEL(1) + FC*XYZLER(1)
          RLE(2,idx_strip)  = (1.0-FC)*XYZLEL(2) + FC*XYZLER(2)
          RLE(3,idx_strip)  = (1.0-FC)*XYZLEL(3) + FC*XYZLER(3)
          CHORD(idx_strip)  = (1.0-FC)*CHORDL    + FC*CHORDR
C
          WSTRIP(idx_strip) = ABS(F2-F1)*WIDTH
          TANLE(idx_strip)  = (XYZLER(1)-XYZLEL(1))/WIDTH
          TANTE(idx_strip)  = (XYZLER(1)+CHORDR - 
     &                         XYZLEL(1)-CHORDL)/WIDTH
C
cc#ifdef USE_CPOML
          CHSIN = CHSINL + F1*(CHSINR-CHSINL)
          CHCOS = CHCOSL + F1*(CHCOSR-CHCOSL)
          AINC1(idx_strip) = ATAN2(CHSIN,CHCOS)
          CHSIN = CHSINL + F2*(CHSINR-CHSINL)
          CHCOS = CHCOSL + F2*(CHCOSR-CHCOSL)
          AINC2(idx_strip) = ATAN2(CHSIN,CHCOS)
C
cc#endif
          CHSIN = CHSINL + FC*(CHSINR-CHSINL)
          CHCOS = CHCOSL + FC*(CHCOSR-CHCOSL)
          AINC(idx_strip) = ATAN2(CHSIN,CHCOS)
C
          DO N = 1, NDESIGN
            CHSIN_G = (1.0-FC)*CHSINL_G(N) + FC*CHSINR_G(N)
            CHCOS_G = (1.0-FC)*CHCOSL_G(N) + FC*CHCOSR_G(N)
            AINC_G(idx_strip,N) = (CHCOS*CHSIN_G - CHSIN*CHCOS_G)
     &                       / (CHSIN**2 + CHCOS**2)
          ENDDO
C
          DO N = 1, NCONTROL
            ICL = ISCONL(N)
            ICR = ISCONR(N)
C
            IF(ICL.EQ.0 .OR. ICR.EQ.0) THEN
C----------- no control effect here
             GAINDA(N) = 0.
             XLED(N) = 0.
             XTED(N) = 0.
C
             VHINGE(1,idx_strip,N) = 0.
             VHINGE(2,idx_strip,N) = 0.
             VHINGE(3,idx_strip,N) = 0.
C
             VREFL(idx_strip,N) = 0.
C
             PHINGE(1,idx_strip,N) = 0.
             PHINGE(2,idx_strip,N) = 0.
             PHINGE(3,idx_strip,N) = 0.
C
            ELSE
C----------- control variable # N is active here
             GAINDA(N) = GAIND(ICL,ISEC  ,ISURF)*(1.0-FC)
     &                 + GAIND(ICR,ISEC+1,ISURF)*     FC
C
             XHD = CHORDL*XHINGED(ICL,ISEC  ,ISURF)*(1.0-FC)
     &           + CHORDR*XHINGED(ICR,ISEC+1,ISURF)*     FC
             IF(XHD.GE.0.0) THEN
C------------ TE control surface, with hinge at XHD
              XLED(N) = XHD
              XTED(N) = CHORD(idx_strip)
             ELSE
C------------ LE control surface, with hinge at -XHD
              XLED(N) =  0.0
              XTED(N) = -XHD
             ENDIF
C
             VHX = VHINGED(1,ICL,ISEC,ISURF)*XYZSCAL(1,ISURF)
             VHY = VHINGED(2,ICL,ISEC,ISURF)*XYZSCAL(2,ISURF)
             VHZ = VHINGED(3,ICL,ISEC,ISURF)*XYZSCAL(3,ISURF)
             VSQ = VHX**2 + VHY**2 + VHZ**2
             IF(VSQ.EQ.0.0) THEN
C------------ default: set hinge vector along hingeline
              VHX = XYZLES(1,ISEC+1,ISURF)
     &              + ABS(CHORDR*XHINGED(ICR,ISEC+1,ISURF))
     &              - XYZLES(1,ISEC  ,ISURF)
     &              - ABS(CHORDL*XHINGED(ICL,ISEC,ISURF))
              VHY = XYZLES(2,ISEC+1,ISURF)
     &            - XYZLES(2,ISEC  ,ISURF)
              VHZ = XYZLES(3,ISEC+1,ISURF)
     &            - XYZLES(3,ISEC  ,ISURF)
              VHX = VHX*XYZSCAL(1,ISURF)
              VHY = VHY*XYZSCAL(2,ISURF)
              VHZ = VHZ*XYZSCAL(3,ISURF)
              VSQ = VHX**2 + VHY**2 + VHZ**2
             ENDIF
C
             VMOD = SQRT(VSQ)
             VHINGE(1,idx_strip,N) = VHX/VMOD
             VHINGE(2,idx_strip,N) = VHY/VMOD
             VHINGE(3,idx_strip,N) = VHZ/VMOD
C
             VREFL(idx_strip,N) = REFLD(ICL,ISEC, ISURF)
C
             IF(XHD .GE. 0.0) THEN
              PHINGE(1,idx_strip,N) = RLE(1,idx_strip) + XHD
              PHINGE(2,idx_strip,N) = RLE(2,idx_strip)
              PHINGE(3,idx_strip,N) = RLE(3,idx_strip)
             ELSE
              PHINGE(1,idx_strip,N) = RLE(1,idx_strip) - XHD
              PHINGE(2,idx_strip,N) = RLE(2,idx_strip)
              PHINGE(3,idx_strip,N) = RLE(3,idx_strip)
             ENDIF
C
            ENDIF
          ENDDO
C
C--- Interpolate CD-CL polar defining data from input sections to strips
          DO L = 1, 6
            CLCD(L,idx_strip) = (1.0-FC)*CLCDSEC(L,ISEC  ,ISURF) 
     &                      +     FC *CLCDSEC(L,ISEC+1,ISURF)
          END DO
C--- If the min drag is zero flag the strip as no-viscous data
          LVISCSTRP(idx_strip) = (CLCD(4,idx_strip).NE.0.0)
C
C
      !     IJFRST(idx_strip) = NVOR + 1
          if (idx_strip ==1) then 
            IJFRST(idx_strip) = 1
          ELSE
            IJFRST(idx_strip) = IJFRST(idx_strip - 1) + 
     &                          NVSTRP(idx_strip - 1)
          endif
          
          NVSTRP(idx_strip) = NVC(ISURF)
!           write(*,*) 'IJFRST(idx_strip)', IJFRST(idx_strip),
!      &               'NVSTRP(idx_strip)', IJFRST(idx_strip - 1) + NVC(ISURF)
C
          NSURFS(idx_strip) = ISURF
C
          NSL = NASEC(ISEC  , ISURF)
          NSR = NASEC(ISEC+1, ISURF)
C
          CHORDC = CHORD(idx_strip)
C
          CLAFC =  (1.-FC)*(CHORDL/CHORDC)*CLAFL
     &           +     FC *(CHORDR/CHORDC)*CLAFR
C
C-------- set chordwise spacing fraction arrays
          CALL CSPACER(NVC(ISURF),CSPACE(ISURF),CLAFC, XPT,XVR,XSR,XCP)
c
C-------- go over vortices in this strip
          idx_vor = IJFRST(idx_strip)
          DO 1505 IVC = 1, NVC(ISURF)
            ! NVOR = NVOR + 1
            ! change all NVOR indices into idx_vor
            ! change all NSTRIP indices into idx_strip
C
            RV1(1,idx_vor) = RLE1(1,idx_strip)
     &                        + XVR(IVC)*CHORD1(idx_strip)
            RV1(2,idx_vor) = RLE1(2,idx_strip)
            RV1(3,idx_vor) = RLE1(3,idx_strip)
C
            RV2(1,idx_vor) = RLE2(1,idx_strip) 
     &                       + XVR(IVC)*CHORD2(idx_strip)
            RV2(2,idx_vor) = RLE2(2,idx_strip)
            RV2(3,idx_vor) = RLE2(3,idx_strip)
C
            RV(1,idx_vor) = RLE(1,idx_strip) + XVR(IVC)*CHORDC
            RV(2,idx_vor) = RLE(2,idx_strip)
            RV(3,idx_vor) = RLE(3,idx_strip)
C
            RC(1,idx_vor) = RLE(1,idx_strip) + XCP(IVC)*CHORDC
            RC(2,idx_vor) = RLE(2,idx_strip)
            RC(3,idx_vor) = RLE(3,idx_strip)
C
            RS(1,idx_vor) = RLE(1,idx_strip) + XSR(IVC)*CHORDC
            RS(2,idx_vor) = RLE(2,idx_strip)
            RS(3,idx_vor) = RLE(3,idx_strip)
C
            CALL AKIMA(XASEC(1,ISEC,  ISURF),SASEC(1,ISEC,  ISURF),NSL,
     &                 XCP(IVC),SLOPEL, DSDX)
            CALL AKIMA(XASEC(1,ISEC+1,ISURF),SASEC(1,ISEC+1,ISURF),NSR,
     &                 XCP(IVC),SLOPER, DSDX)
            SLOPEC(idx_vor) =  (1.-FC)*(CHORDL/CHORDC)*SLOPEL 
     &                    +     FC *(CHORDR/CHORDC)*SLOPER
C
            CALL AKIMA(XASEC(1,ISEC  ,ISURF),SASEC(1,ISEC  ,ISURF),NSL,
     &                 XVR(IVC),SLOPEL, DSDX)
            CALL AKIMA(XASEC(1,ISEC+1,ISURF),SASEC(1,ISEC+1,ISURF),NSR,
     &                 XVR(IVC),SLOPER, DSDX)
            SLOPEV(idx_vor) =  (1.-FC)*(CHORDL/CHORDC)*SLOPEL 
     &                    +     FC *(CHORDR/CHORDC)*SLOPER
C
            DXOC = XPT(IVC+1) - XPT(IVC)
            DXV(idx_vor) = DXOC*CHORDC
            CHORDV(idx_vor) = CHORDC
            NSURFV(idx_vor) = LSCOMP(ISURF)

            LVNC(idx_vor) = .TRUE.
C
C---------- element inherits alpha,beta flag from surface
            LVALBE(idx_vor) = LFALBE(ISURF)
C
            DO N = 1, NCONTROL
C------------ scale control gain by factor 0..1, (fraction of element on control surface)
              FRACLE = (XLED(N)/CHORDC-XPT(IVC)) / DXOC
              FRACTE = (XTED(N)/CHORDC-XPT(IVC)) / DXOC
C
              FRACLE = MIN( 1.0 , MAX( 0.0 , FRACLE ) )
              FRACTE = MIN( 1.0 , MAX( 0.0 , FRACTE ) )
C
              DCONTROL(idx_vor,N) = GAINDA(N)*(FRACTE-FRACLE)
            ENDDO
C
C---------- TE control point used only if surface sheds a wake
            LVNC(idx_vor) = LFWAKE(ISURF)
C
cc#ifdef USE_CPOML
C...        nodal grid associated with vortex strip (aft-panel nodes)
C...        NOTE: airfoil in plane of wing, but not rotated perpendicular to dihedral;
C...        retained in (x,z) plane at this point
            CALL AKIMA( XLASEC(1,ISEC,ISURF), ZLASEC(1,ISEC,ISURF), NSL,
     &                  XPT(IVC+1), ZL_L, DSDX )
            CALL AKIMA( XUASEC(1,ISEC,ISURF), ZUASEC(1,ISEC,ISURF), NSL,
     &                  XPT(IVC+1), ZU_L, DSDX )
C
            CALL AKIMA( XLASEC(1,ISEC+1,ISURF), ZLASEC(1,ISEC+1,ISURF),
     &                  NSR, XPT(IVC+1), ZL_R, DSDX )
            CALL AKIMA( XUASEC(1,ISEC+1,ISURF), ZUASEC(1,ISEC+1,ISURF),
     &                  NSR, XPT(IVC+1), ZU_R, DSDX )
C
            XYN1(1,idx_vor) = RLE1(1,idx_strip) + 
     &                        XPT(IVC+1)*CHORD1(idx_strip)
            XYN1(2,idx_vor) = RLE1(2,idx_strip)
            
            ZL =  (1.-F1)*ZL_L + F1 *ZL_R
            ZU =  (1.-F1)*ZU_L + F1 *ZU_R
            
            ZLON1(idx_vor)  = RLE1(3,idx_strip) + ZL*CHORD1(idx_strip)
            ZUPN1(idx_vor)  = RLE1(3,idx_strip) + ZU*CHORD1(idx_strip)
C
            XYN2(1,idx_vor) = RLE2(1,idx_strip) + 
     &                        XPT(IVC+1)*CHORD2(idx_strip)
            XYN2(2,idx_vor) = RLE2(2,idx_strip)
            
            ZL =  (1.-F2)*ZL_L + F2 *ZL_R
            ZU =  (1.-F2)*ZU_L + F2 *ZU_R
          
            ZLON2(idx_vor)  = RLE2(3,idx_strip) + ZL*CHORD2(idx_strip)
            ZUPN2(idx_vor)  = RLE2(3,idx_strip) + ZU*CHORD2(idx_strip)
C
cc#endif
            idx_vor = idx_vor + 1
 1505     CONTINUE
C           
        idx_strip = idx_strip + 1
 150    CONTINUE
C
 200  CONTINUE
C
C---- Find wetted surface area (one side)
      SUM  = 0.0
      WTOT = 0.0
      DO JJ = 1, NJ(ISURF)
        J = JFRST(ISURF) + JJ-1 
        ASTRP = WSTRIP(J)*CHORD(J)
        SUM  = SUM + ASTRP
        WTOT = WTOT + WSTRIP(J)
      ENDDO
      SSURF(ISURF) = SUM
C
      IF(WTOT .EQ. 0.0) THEN
       CAVESURF(ISURF) = 0.0
      ELSE
       CAVESURF(ISURF) = SUM/WTOT
      ENDIF
C     add number of strips to the global count
      NSTRIP = NSTRIP + NJ(ISURF)
C     add number of of votrices
      NVOR = NVOR + NK(ISURF)*NJ(ISURF) 
C
      RETURN
      END ! MAKESURF
      
      subroutine update_surfaces()
c--------------------------------------------------------------
c     Updates all surfaces, using the stored data.
c--------------------------------------------------------------
      
      include 'AVL.INC'
      integer ISURF
      
      NSTRIP = 0
      NVOR = 0
      
      do ISURF=1,NSURF
            if (lverbose) then 
                  write(*,*) 'Updating surface ',ISURF
            end if
            if (ISURF.ne.1) then
                  if(ldupl(isurf-1)) then 
                        ! this surface has already been created
                        ! it was probably duplicated from the previous one
                        cycle
                  end if
                  call makesurf(ISURF)
            else
                  call makesurf(ISURF)
            endif
            
            if(ldupl(isurf)) then
                  call sdupl(isurf,ydupl(isurf),'ydup')
            endif
      end do 
      
      CALL ENCALC
      
      LAIC = .FALSE.
      LSRD = .FALSE.
      LVEL = .FALSE.
      LSOL = .FALSE.
      LSEN = .FALSE.
      
      end subroutine update_surfaces
            


      SUBROUTINE MAKEBODY(IBODY)
C     &       XBOD,YBOD,TBOD,NBOD)
C--------------------------------------------------------------
C     Sets up all stuff for body IBODY,
C     using info from configuration input file.
C--------------------------------------------------------------
      INCLUDE 'AVL.INC'
C
C      REAL XBOD(IBX), YBOD(IBX), TBOD(IBX)
C
      PARAMETER (KLMAX=101)
      REAL XPT(KLMAX), FSPACE(KLMAX)
C
C
c      IF(NSEC(IBODY).LT.2) THEN
c       WRITE(*,*) '*** Need at least 2 sections per body.'
c       STOP
c      ENDIF
C
C
      IF(NVB(IBODY).GT.KLMAX) THEN
       WRITE(*,*) '* MAKEBODY: Array overflow.  Increase KLMAX to',
     & NVB(IBODY)
       NVB(IBODY) = KLMAX
      ENDIF
C
C
      LFRST(IBODY) = NLNODE + 1 
      NL(IBODY) = NVB(IBODY)
C
      IF(NLNODE+NVB(IBODY).GT.NLMAX) THEN
       WRITE(*,*) '*** MAKEBODY: Array overflow. Increase NLMAX to',
     &             NLNODE+NVB(IBODY)
       STOP
      ENDIF
C
C-----------------------------------------------------------------
C---- set lengthwise spacing fraction arrays
      NSPACE = NVB(IBODY) + 1
      IF(NSPACE.GT.KLMAX) THEN
       WRITE(*,*) '*** MAKEBODY: Array overflow. Increase KLMAX to', 
     &             NSPACE
       STOP
      ENDIF
      CALL SPACER(NSPACE,BSPACE(IBODY),FSPACE)
C
      DO IVB = 1, NVB(IBODY)
        XPT(IVB) = FSPACE(IVB)
      ENDDO
      XPT(1) = 0.0
      XPT(NVB(IBODY)+1) = 1.0
C
C---- set body nodes and radii
      VOLB = 0.0
      SRFB = 0.0
      DO IVB = 1, NVB(IBODY)+1
        NLNODE = NLNODE + 1
C
        XVB = XBOD(1, IBODY) + (XBOD(NBOD(IBODY), IBODY)-XBOD(1,IBODY))
     &  *XPT(IVB)
        CALL AKIMA(XBOD(1,IBODY),YBOD(1,IBODY),NBOD(IBODY),XVB,YVB,DYDX)
        RL(1,NLNODE) = XYZTRAN_B(1,IBODY) + XYZSCAL_B(1,IBODY)*XVB
        RL(2,NLNODE) = XYZTRAN_B(2,IBODY)
        RL(3,NLNODE) = XYZTRAN_B(3,IBODY) + XYZSCAL_B(3,IBODY)*YVB
C
        CALL AKIMA(XBOD(1,IBODY),TBOD(1,IBODY),NBOD(IBODY),XVB,TVB,DRDX)
        RADL(NLNODE) = SQRT(XYZSCAL_B(2,IBODY)*XYZSCAL_B(3,IBODY)) 
     & * 0.5*TVB
      ENDDO
C---- get surface length, area and volume
      VOLB = 0.0
      SRFB = 0.0
      XBMN = RL(1,LFRST(IBODY))
      XBMX = XBMN
      DO IVB = 1, NVB(IBODY)
        NL0 = LFRST(IBODY) + IVB-1
        NL1 = NL0 + 1
        X0 = RL(1,NL0)
        X1 = RL(1,NL1)
        DX = ABS(X1 - X0)
        R0 = RADL(NL0)
        R1 = RADL(NL1)
        DVOL = PI*DX * (R0**2 + R0*R1 + R1**2) / 3.0
        DS = SQRT((R0-R1)**2 + DX**2)
        DSRF = PI*DS * (R0+R1)
C
        SRFB = SRFB + DSRF
        VOLB = VOLB + DVOL
        XBMN = MIN(XBMN,X0,X1)
        XBMX = MAX(XBMX,X0,X1)
      ENDDO
      VOLBDY(IBODY) = VOLB
      SRFBDY(IBODY) = SRFB
      ELBDY(IBODY) = XBMX - XBMN
C
      RETURN
      END ! MAKEBODY

      subroutine update_bodies()
c--------------------------------------------------------------
c     Updates all bodies, using the stored data.
c--------------------------------------------------------------
      include 'AVL.INC'
      integer IBODY, NBOD, NBLDS
      character*120 upname

      NLNODE = 0

      do IBODY=1,NBODY
      if (lverbose) then 
        write(*,*) 'Updating body ',IBODY
      end if
      if (IBODY.ne.1) then
       if(ldupl_b(ibody-1)) then 
        ! this body has already been created
        ! it was probably duplicated from the previous one
        cycle
       end if
       call makebody(IBODY)
      else
       call makebody(IBODY)
      endif
      
      if(ldupl_b(ibody)) then
       call bdupl(ibody,ydupl_b(ibody),'ydup')
      endif
      end do 
      
      CALL ENCALC
      
      LAIC = .FALSE.
      LSRD = .FALSE.
      LVEL = .FALSE.
      LSOL = .FALSE.
      LSEN = .FALSE.
      
      end subroutine update_bodies

      subroutine  set_section_coordinates(isec,isurf,x,y,n,nin,xfmin,
     &       xfmax, storecoords)
c--------------------------------------------------------------
c     Sets the airfoil coodinate data for the given section and surface
c--------------------------------------------------------------
      include 'AVL.INC'
      integer isec, isurf, n, nin
      real x(n), y(n)
      real xin(IBX), yin(IBX), tin(IBX)
      real xfmin, xfmax
      logical storecoords

        if (storecoords) then
c--------------------------------------------------------------
c     Store the raw input data into the common block for general purposes
c--------------------------------------------------------------
        do i = 1,n
          XSEC(i,isec, isurf) = x(i)
          YSEC(i,isec, isurf) = y(i)
        end do
        XFMIN_R(isec,isurf) = xfmin
        XFMAX_R(isec,isurf) = xfmax
        end if

        if((xfmin .gt. 0.01) .or. (xfmax .lt. 0.99)) then
          write(*,*) 'aifoil Lrange false', isurf, isec
          write(*,*) (xfmin .gt. 0.01)
          write(*,*) (xfmax .lt. 0.99)
          LRANGE(isurf) = .false.
        else
          LRANGE(isurf) = .true.
        end if

        call GETCAM(x,y,n,xin,yin,tin,nin,.true.)

        NASEC(isec,isurf) = nin

        do i = 1, nin
          xf = xfmin + 
     &         (xfmax-xfmin)*float(i-1)/float(NASEC(isec,isurf)-1)
          XASEC(i,isec,isurf) = xin(1) + xf*(xin(nin)-xin(1))
          call AKIMA(xin,yin,nin,XASEC(i,isec,isurf),zc,
     &               SASEC(i,isec,isurf))
          call AKIMA(xin,tin,nin,XASEC(i,isec,isurf),
     &               TASEC(i,isec,isurf),dummy)
          XLASEC(i,isec,isurf) = XASEC(i,isec,isurf)
          XUASEC(i,isec,isurf) = XASEC(i,isec,isurf)
          ZLASEC(i,isec,isurf) = zc - 0.5*TASEC(i,isec,isurf)
          ZUASEC(i,isec,isurf) = zc + 0.5*TASEC(i,isec,isurf)
          CASEC(i,isec,isurf) = zc

        end do
        call NRMLIZ(NASEC(isec,isurf),XASEC(1,isec,isurf))
      
      end subroutine set_section_coordinates

      subroutine set_body_coordinates(ibod,xb,yb,nb,nin,storecoords)
c--------------------------------------------------------------
c     Sets the body oml coodinate data for the given section and surface
c--------------------------------------------------------------
      include 'AVL.INC'
      integer ibod, nb, nin
      real xb(nb), yb(nb)
      logical storecoords
C      real xin(ibx), yin(ibx), tin(ibx)
C------ xfmin and xfmax don't appear to be supported for bodies?

      if (storecoords) then
c--------------------------------------------------------------
c     Store the raw input data into the common block for general purposes
c--------------------------------------------------------------
      do i = 1,nb
      XBOD_R(i,ibod) = xb(i)
      YBOD_R(i,ibod) = yb(i)
      end do
      end if

C------ set thread line y, and thickness t ( = 2r)
      nbod = MIN( 50 , IBX )
      call GETCAM(xb,yb,nb,xbod(1,ibod),ybod(1,ibod),tbod(1,ibod),nin,
     & .false.)
      end subroutine set_body_coordinates

      SUBROUTINE SDUPL(NN, Ypt,MSG)
C-----------------------------------
C     Adds image of surface NN,
C     reflected about y=Ypt.
C-----------------------------------
      INCLUDE 'AVL.INC'
      CHARACTER*(*) MSG
      integer idx_vor
C
C     
      NNI = NN + 1
      IF(NNI.GT.NFMAX) THEN
        WRITE(*,*) 'SDUPL: Surface array overflow. Increase NFMAX',
     &             ' currently ',NFMAX
        STOP
      ENDIF
C
      KLEN = LEN(STITLE(NN))
      DO K = KLEN, 1, -1
        IF(STITLE(NN)(K:K) .NE. ' ') GO TO 6
      ENDDO
 6    STITLE(NNI) = STITLE(NN)(1:K) // ' (' // MSG // ')'
      if(lverbose)then
       WRITE(*,*) ' '
       WRITE(*,*) '  Building duplicate image-surface: ',STITLE(NNI)
      endif
C
C---- duplicate surface is assumed to be the same logical component surface
      LSCOMP(NNI) = LSCOMP(NN)
C
C---- same various logical flags
      LFWAKE(NNI) = LFWAKE(NN)
      LFALBE(NNI) = LFALBE(NN)
      LFLOAD(NNI) = LFLOAD(NN)
      LRANGE(NNI) = LRANGE(NN)
      LSURFSPACING(NNI) = LSURFSPACING(NN)

C---- accumulate stuff for new image surface 
      ! IFRST(NNI) = NVOR   + 1
      IFRST(NNI) =  IFRST(NNI-1) +  NK(NNI-1)*NJ(NNI-1)
      JFRST(NNI) =  JFRST(NNI-1) +  NJ(NNI-1)
      ! JFRST(NNI) = NSTRIP + 1
      NJ(NNI) = NJ(NN)
      NK(NNI) = NK(NN)
C
      NVC(NNI) = NK(NNI)
      NVS(NNi) = NJ(NNI)
C
      SSURF(NNI)    = SSURF(NN)
      CAVESURF(NNI) = CAVESURF(NN)
C--- Note hinge axis is flipped to reverse the Y component of the hinge
C    vector.   This means that deflections need to be reversed for image
C    surfaces.
C
C--- Image flag reversed (set to -IMAGS) for imaged surfaces
      IMAGS(NNI) = -IMAGS(NN)
C
cc#ifdef USE_CPOML
      ICNTFRST(NNI) = ICNTFRST(NN) + NCNTSEC(NN)
      NCNTSEC(NNI)  = NCNTSEC(NN)
      DO ISEC = 1, NCNTSEC(NNI)
            IDUP = ICNTFRST(NNI) + (ISEC-1)
            IORG = ICNTFRST(NN ) + (ISEC-1)
            ICNTSEC(IDUP) = ICNTSEC(IORG)
      ENDDO
cc#endif
C
      YOFF = 2.0*Ypt
C
C--- Create image strips, to maintain the same sense of positive GAMMA
C    these have the 1 and 2 strip edges reversed (i.e. root is edge 2, 
C    not edge 1 as for a strip with IMAGS=1
      idx_strip = JFRST(NNI)
      DO 100 IVS = 1, NVS(NNI)
      !   NSTRIP = NSTRIP + 1
        IF(idx_strip.GT.NSMAX) THEN
          WRITE(*,*) 'SDUPL: Strip array overflow. Increase NSMAX',
     &               ' currently ',NSMAX
          STOP
        ENDIF
C
        JJI = JFRST(NNI) + IVS-1
        JJ  = JFRST(NN)  + IVS-1
        RLE1(1,JJI)   =  RLE2(1,JJ)
        RLE1(2,JJI)   = -RLE2(2,JJ) + YOFF
        RLE1(3,JJI)   =  RLE2(3,JJ)
        CHORD1(JJI) =  CHORD2(JJ)
        RLE2(1,JJI)   =  RLE1(1,JJ)
        RLE2(2,JJI)   = -RLE1(2,JJ) + YOFF
        RLE2(3,JJI)   =  RLE1(3,JJ)
        CHORD2(JJI) =  CHORD1(JJ)
        RLE(1,JJI)    =  RLE(1,JJ)
        RLE(2,JJI)    = -RLE(2,JJ) + YOFF
        RLE(3,JJI)    =  RLE(3,JJ)
        CHORD(JJI)  =  CHORD(JJ)
        WSTRIP(JJI) =  WSTRIP(JJ)
        TANLE(JJI)  = -TANLE(JJ)
        AINC (JJI)  =  AINC(JJ)
C
cc#ifdef USE_CPOML
        AINC1(JJI) = AINC2(JJ)
        AINC2(JJI) = AINC1(JJ)
C
cc#endif
        NSURFS(idx_strip) = NNI
C
        DO N = 1, NDESIGN
          AINC_G(JJI,N) = AINC_G(JJ,N)
        ENDDO
C
        DO N = 1, NCONTROL
          VREFL(JJI,N) = VREFL(JJ,N)
C
          VHINGE(1,JJI,N) =  VHINGE(1,JJ,N)
          VHINGE(2,JJI,N) = -VHINGE(2,JJ,N)
          VHINGE(3,JJI,N) =  VHINGE(3,JJ,N)
C
          PHINGE(1,JJI,N) =  PHINGE(1,JJ,N)
          PHINGE(2,JJI,N) = -PHINGE(2,JJ,N) + YOFF
          PHINGE(3,JJI,N) =  PHINGE(3,JJ,N)
        ENDDO
C
C--- The defined section for image strip is flagged with (-)
      !   IJFRST(JJI)  = NVOR + 1
      !   IJFRST(JJI) = IJFRST(NSTRIP - 1) + NVC(NNI)
        IJFRST(JJI) = IJFRST(JJI - 1) + NVSTRP(JJI - 1)

        NVSTRP(JJI)  = NVC(NNI)
        DO L = 1, 6
          CLCD(L,JJI) = CLCD(L,JJ) 
        END DO
        LVISCSTRP(JJI) = LVISCSTRP(JJ)
C
        idx_vor = IJFRST(JJI)

        DO 80 IVC = 1, NVC(NNI)
      !     NVOR = NVOR + 1
          IF(idx_vor.GT.NVMAX) THEN
            WRITE(*,*) 'SDUPL: Vortex array overflow. Increase NVMAX',
     &                 ' currently ',NVMAX
            STOP
          ENDIF
C
          III = IJFRST(JJI) + IVC-1
          II  = IJFRST(JJ)  + IVC-1
          RV1(1,III)     =  RV2(1,II)
          RV1(2,III)     = -RV2(2,II) + YOFF
          RV1(3,III)     =  RV2(3,II)
          RV2(1,III)     =  RV1(1,II)
          RV2(2,III)     = -RV1(2,II) + YOFF
          RV2(3,III)     =  RV1(3,II)
          RV(1,III)     =  RV(1,II)
          RV(2,III)     = -RV(2,II) + YOFF
          RV(3,III)     =  RV(3,II)
          RC(1,III)     =  RC(1,II)
          RC(2,III)     = -RC(2,II) + YOFF
          RC(3,III)     =  RC(3,II)
          SLOPEC(III) = SLOPEC(II)
          SLOPEV(III) = SLOPEV(II)
          DXV(III)     = DXV(II)
          CHORDV(III) = CHORDV(II)
          NSURFV(III) = LSCOMP(NNI)
          LVALBE(III) = LVALBE(II)
          LVNC(III) = LVNC(II)
C
          DO N = 1, NCONTROL
ccc         RSGN = SIGN( 1.0 , VREFL(JJ,N) )
            RSGN = VREFL(JJ,N)
            DCONTROL(III,N) = -DCONTROL(II,N)*RSGN
          ENDDO
C          
cc#ifdef USE_CPOML
C...      nodal grid associated with vortex strip
          XYN1(1,III) =  XYN2(1,II)
          XYN1(2,III) = -XYN2(2,II) + YOFF
          XYN2(1,III) =  XYN1(1,II)
          XYN2(2,III) = -XYN1(2,II) + YOFF
C
          ZLON1(III)  = ZLON2(II)
          ZUPN1(III)  = ZUPN2(II)
          ZLON2(III)  = ZLON1(II)
          ZUPN2(III)  = ZUPN1(II)
cc#endif
          idx_vor = idx_vor + 1
C
   80   CONTINUE
        idx_strip = idx_strip + 1 
C
  100 CONTINUE
C
      
C
      NSTRIP = NSTRIP + NJ(NNI)
      NVOR = NVOR + NK(NNI)*NJ(NNI) 

      RETURN
      END ! SDUPL




      SUBROUTINE BDUPL(NN,Ypt,MSG)
C-----------------------------------
C     Adds image of surface NN,
C     reflected about y=Ypt.
C-----------------------------------
      INCLUDE 'AVL.INC'
      CHARACTER*(*) MSG
C
      NNI = NBODY + 1
      IF(NNI.GT.NBMAX) THEN
        WRITE(*,*) 'BDUPL: Body array overflow. Increase NBMAX',
     &             ' currently ',NBMAX
        STOP
      ENDIF
C
      KLEN = LEN(BTITLE(NN))
      DO K = KLEN, 1, -1
        IF(BTITLE(NN)(K:K) .NE. ' ') GO TO 6
      ENDDO
 6    BTITLE(NNI) = BTITLE(NN)(1:K) // ' (' // MSG // ')'
      if (lverbose) then
      WRITE(*,*) ' '
      WRITE(*,*) '  Building duplicate image-body: ',BTITLE(NNI)
      endif
C
      LFRST(NNI) = NLNODE + 1
      NL(NNI) = NL(NN)
C
      NVB(NNI) = NL(NNI)
C
      IF(NLNODE+NVB(NNI).GT.NLMAX) THEN
       WRITE(*,*) '*** MAKEBODY: Array overflow. Increase NLMAX to',
     &             NLNODE+NVB(NNI)
       STOP
      ENDIF
C
C
      ELBDY(NNI)  = ELBDY(NN)
      SRFBDY(NNI) = SRFBDY(NN)
      VOLBDY(NNI) = VOLBDY(NN)
C
      YOFF = 2.0*Ypt
C
C---- set body nodes and radii
      DO IVB = 1, NVB(NNI)+1
        NLNODE = NLNODE + 1
C
        LLI = LFRST(NNI) + IVB-1
        LL  = LFRST(NN)  + IVB-1
C
        RL(1,LLI) =  RL(1,LL)
        RL(2,LLI) = -RL(2,LL) + YOFF
        RL(3,LLI) =  RL(3,LL)
C
        RADL(LLI) =  RADL(LL)
      ENDDO
C
      NBODY = NBODY + 1
C
      RETURN
      END ! BDUPL




      SUBROUTINE ENCALC
C
C...PURPOSE  To calculate the normal vectors for the strips, 
C            the horseshoe vortices, and the control points.
C            Incorporates surface deflections.
C
C...INPUT    NVOR      Number of vortices
C            X1        Coordinates of endpoint #1 of the vortices
C            X2        Coordinates of endpoint #2 of the vortices
C            SLOPEV    Slope at bound vortices
C            SLOPEC    Slope at control points
C            NSTRIP    Number of strips
C            IJFRST    Index of first element in strip
C            NVSTRP    No. of vortices in strip
C            AINC      Angle of incidence of strip
C            LDES      include design-variable deflections if TRUE
C
C...OUTPUT   ENC(3)        Normal vector at control point
C            ENV(3)        Normal vector at bound vortices
C            ENSY, ENSZ    Strip normal vector (ENSX=0)
C            LSTRIPOFF     Non-used strip (T) (below z=ZSYM)
C
C...COMMENTS   
C
      INCLUDE 'AVL.INC'
C
      REAL EP(3), EQ(3), ES(3), EB(3), EC(3), ECXB(3)
      REAL EC_G(3,NDMAX), ECXB_G(3)
C
C...Calculate the normal vector at control points and bound vortex midpoints
C
      DO 10 J = 1, NSTRIP
C
C...Calculate normal vector for the strip (normal to X axis)
        I = IJFRST(J)
        DXLE =  RV2(1,I)-RV1(1,I)
        DYLE =  RV2(2,I)-RV1(2,I)
        DZLE =  RV2(3,I)-RV1(3,I)
c       AXLE = (RV2(1,I)+RV1(1,I))*0.5
c       AYLE = (RV2(2,I)+RV1(2,I))*0.5
c       AZLE = (RV2(3,I)+RV1(3,I))*0.5
        AXLE = RV(1,I)
        AYLE = RV(2,I)
        AZLE = RV(3,I)
C
        I = IJFRST(J) + (NVSTRP(J)-1)
        DXTE =  RV2(1,I)-RV1(1,I)
        DYTE =  RV2(2,I)-RV1(2,I)
        DZTE =  RV2(3,I)-RV1(3,I)
c       AXTE = (RV2(1,I)+RV1(1,I))*0.5
c       AYTE = (RV2(2,I)+RV1(2,I))*0.5
c       AZTE = (RV2(3,I)+RV1(3,I))*0.5
        AXTE = RV(1,I)
        AYTE = RV(2,I)
        AZTE = RV(3,I)
C
        DXT = (1.0-SAXFR)*DXLE + SAXFR*DXTE
        DYT = (1.0-SAXFR)*DYLE + SAXFR*DYTE
        DZT = (1.0-SAXFR)*DZLE + SAXFR*DZTE
C
        ESS(1,J) =  DXT/SQRT(DXT*DXT + DYT*DYT + DZT*DZT)
        ESS(2,J) =  DYT/SQRT(DXT*DXT + DYT*DYT + DZT*DZT)
        ESS(3,J) =  DZT/SQRT(DXT*DXT + DYT*DYT + DZT*DZT)
C
        ENSY(J) = -DZT/SQRT(DYT*DYT + DZT*DZT)
        ENSZ(J) =  DYT/SQRT(DYT*DYT + DZT*DZT)
C
        XSREF(J) = (1.0-SAXFR)*AXLE + SAXFR*AXTE
        YSREF(J) = (1.0-SAXFR)*AYLE + SAXFR*AYTE
        ZSREF(J) = (1.0-SAXFR)*AZLE + SAXFR*AZTE
C
C
        ES(1) = 0.
        ES(2) = ENSY(J)
        ES(3) = ENSZ(J)
C
        LSTRIPOFF(J) = .FALSE.
C
        NV = NVSTRP(J)
        DO 105 II = 1, NV
C
          I = IJFRST(J) + (II-1)
C
          DO N = 1, NCONTROL
            ENV_D(1,I,N) = 0.
            ENV_D(2,I,N) = 0.
            ENV_D(3,I,N) = 0.
            ENC_D(1,I,N) = 0.
            ENC_D(2,I,N) = 0.
            ENC_D(3,I,N) = 0.
          ENDDO
C
          DO N = 1, NDESIGN
            ENV_G(1,I,N) = 0.
            ENV_G(2,I,N) = 0.
            ENV_G(3,I,N) = 0.
            ENC_G(1,I,N) = 0.
            ENC_G(2,I,N) = 0.
            ENC_G(3,I,N) = 0.
          ENDDO
C
C...Define unit vector along bound leg
          DXB = RV2(1,I)-RV1(1,I) ! right h.v. pt - left h.v. pt 
          DYB = RV2(2,I)-RV1(2,I)
          DZB = RV2(3,I)-RV1(3,I)
          EMAG = SQRT(DXB**2 + DYB**2 + DZB**2)
          EB(1) = DXB/EMAG
          EB(2) = DYB/EMAG
          EB(3) = DZB/EMAG
C
C...Define direction of normal vector at control point 
C   The YZ projection of the normal vector matches the camber slope
C   + section local incidence in the YZ defining plane for the section
          ANG = AINC(J) - ATAN(SLOPEC(I))
cc          IF(LDES) THEN
C--------- add design-variable contribution to angle
           DO N = 1, NDESIGN
             ANG = ANG + AINC_G(J,N)*DELDES(N)
           ENDDO
cc          ENDIF
C
          SINC = SIN(ANG)
          COSC = COS(ANG)
          EC(1) =  COSC
          EC(2) = -SINC*ES(2)
          EC(3) = -SINC*ES(3)
      !     EC  = rotation of strip normal vector? or along chord?
          DO N = 1, NDESIGN
            EC_G(1,N) = -SINC      *AINC_G(J,N)
            EC_G(2,N) = -COSC*ES(2)*AINC_G(J,N)
            EC_G(3,N) = -COSC*ES(3)*AINC_G(J,N)
          ENDDO
C
C...Normal vector is perpendicular to camberline vector and to the bound leg
          CALL CROSS(EC,EB,ECXB)
          EMAG = SQRT(ECXB(1)**2 + ECXB(2)**2 + ECXB(3)**2)
          IF(EMAG.NE.0.0) THEN
            ENC(1,I) = ECXB(1)/EMAG
            ENC(2,I) = ECXB(2)/EMAG
            ENC(3,I) = ECXB(3)/EMAG
            DO N = 1, NDESIGN
              CALL CROSS(EC_G(1,N),EB,ECXB_G)
              EMAG_G = ENC(1,I)*ECXB_G(1)
     &               + ENC(2,I)*ECXB_G(2)
     &               + ENC(3,I)*ECXB_G(3)
              ENC_G(1,I,N) = (ECXB_G(1) - ENC(1,I)*EMAG_G)/EMAG
              ENC_G(2,I,N) = (ECXB_G(2) - ENC(2,I)*EMAG_G)/EMAG
              ENC_G(3,I,N) = (ECXB_G(3) - ENC(3,I)*EMAG_G)/EMAG
            ENDDO
          ELSE
            ENC(1,I) = ES(1)
            ENC(2,I) = ES(2)
            ENC(3,I) = ES(3)
          ENDIF
C
C
C...Define direction of normal vector at vortex mid-point. 
C   The YZ projection of the normal vector matches the camber slope
C   + section local incidence in the YZ defining plane for the section
          ANG = AINC(J) - ATAN(SLOPEV(I)) 
cc          IF(LDES) THEN
C--------- add design-variable contribution to angle
           DO N = 1, NDESIGN
             ANG = ANG + AINC_G(J,N)*DELDES(N)
           ENDDO
cc          ENDIF
C
          SINC = SIN(ANG)
          COSC = COS(ANG)
          EC(1) =  COSC
          EC(2) = -SINC*ES(2)
          EC(3) = -SINC*ES(3)
          DO N = 1, NDESIGN
            EC_G(1,N) = -SINC      *AINC_G(J,N)
            EC_G(2,N) = -COSC*ES(2)*AINC_G(J,N)
            EC_G(3,N) = -COSC*ES(3)*AINC_G(J,N)
          ENDDO
C
C...Normal vector is perpendicular to camberline vector and to the bound leg
          CALL CROSS(EC,EB,ECXB)
          EMAG = SQRT(ECXB(1)**2 + ECXB(2)**2 + ECXB(3)**2)
          IF(EMAG.NE.0.0) THEN
            ENV(1,I) = ECXB(1)/EMAG
            ENV(2,I) = ECXB(2)/EMAG
            ENV(3,I) = ECXB(3)/EMAG
            DO N = 1, NDESIGN
              CALL CROSS(EC_G(1,N),EB,ECXB_G)
              EMAG_G = ENC(1,I)*ECXB_G(1)
     &               + ENC(2,I)*ECXB_G(2)
     &               + ENC(3,I)*ECXB_G(3)
              ENV_G(1,I,N) = (ECXB_G(1) - ENV(1,I)*EMAG_G)/EMAG
              ENV_G(2,I,N) = (ECXB_G(2) - ENV(2,I)*EMAG_G)/EMAG
              ENV_G(3,I,N) = (ECXB_G(3) - ENV(3,I)*EMAG_G)/EMAG
            ENDDO
          ELSE
            ENV(1,I) = ES(1)
            ENV(2,I) = ES(2)
            ENV(3,I) = ES(3)
          ENDIF
C
C
ccc       write(*,*) i, dcontrol(i,1), dcontrol(i,2)
C
C=======================================================
C-------- rotate normal vectors for control surface
          DO 100 N = 1, NCONTROL
C
C---------- skip everything if this element is unaffected by control variable N
            IF(DCONTROL(I,N).EQ.0.0) GO TO 100
C
            ANG     = DTR*DCONTROL(I,N)*DELCON(N)
            ANG_DDC = DTR*DCONTROL(I,N)
C
            COSD = COS(ANG)
            SIND = SIN(ANG)
C
C---------- EP = normal-vector component perpendicular to hinge line
            ENDOT = DOT(ENC(1,I),VHINGE(1,J,N))
            EP(1) = ENC(1,I) - ENDOT*VHINGE(1,J,N)
            EP(2) = ENC(2,I) - ENDOT*VHINGE(2,J,N)
            EP(3) = ENC(3,I) - ENDOT*VHINGE(3,J,N)
C---------- EQ = unit vector perpendicular to both EP and hinge line
            CALL CROSS(VHINGE(1,J,N),EP,EQ)
C
C---------- rotated vector would consist of sin,cos parts from EP and EQ,
C-          with hinge-parallel component ENDOT restored 
cc          ENC(1,I) = EP(1)*COSD + EQ(1)*SIND + ENDOT*VHINGE(1,J,N)
cc          ENC(2,I) = EP(2)*COSD + EQ(2)*SIND + ENDOT*VHINGE(2,J,N)
cc          ENC(3,I) = EP(3)*COSD + EQ(3)*SIND + ENDOT*VHINGE(3,J,N)
C
C---------- linearize about zero deflection (COSD=1, SIND=0)
            ENC_D(1,I,N) = ENC_D(1,I,N) + EQ(1)*ANG_DDC
            ENC_D(2,I,N) = ENC_D(2,I,N) + EQ(2)*ANG_DDC
            ENC_D(3,I,N) = ENC_D(3,I,N) + EQ(3)*ANG_DDC
C
C
C---------- repeat for ENV vector
C
C---------- EP = normal-vector component perpendicular to hinge line
            ENDOT = DOT(ENV(1,I),VHINGE(1,J,N))
            EP(1) = ENV(1,I) - ENDOT*VHINGE(1,J,N)
            EP(2) = ENV(2,I) - ENDOT*VHINGE(2,J,N)
            EP(3) = ENV(3,I) - ENDOT*VHINGE(3,J,N)
C---------- EQ = unit vector perpendicular to both EP and hinge line
            CALL CROSS(VHINGE(1,J,N),EP,EQ)
C
C---------- rotated vector would consist of sin,cos parts from EP and EQ,
C-          with hinge-parallel component ENDOT restored 
cc          ENV(1,I) = EP(1)*COSD + EQ(1)*SIND + ENDOT*VHINGE(1,J,N)
cc          ENV(2,I) = EP(2)*COSD + EQ(2)*SIND + ENDOT*VHINGE(2,J,N)
cc          ENV(3,I) = EP(3)*COSD + EQ(3)*SIND + ENDOT*VHINGE(3,J,N)
C
C---------- linearize about zero deflection (COSD=1, SIND=0)
            ENV_D(1,I,N) = ENV_D(1,I,N) + EQ(1)*ANG_DDC
            ENV_D(2,I,N) = ENV_D(2,I,N) + EQ(2)*ANG_DDC
            ENV_D(3,I,N) = ENV_D(3,I,N) + EQ(3)*ANG_DDC
 100      CONTINUE
 101      CONTINUE
C
 105    CONTINUE
  10  CONTINUE
C
      LENC = .TRUE.
C
      RETURN
      END ! ENCALC

