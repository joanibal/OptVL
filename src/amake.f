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
         ! Throws an error in the case where the same node is the closest node 
         ! to two consecutive sections
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
          LSSURF(idx_strip) = ISURF
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
            LVCOMP(idx_vor) = LNCOMP(ISURF)

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

      integer function flatidx(idx_x, idx_y, idx_surf)
      include 'AVL.INC'
      ! store MFRST and  NVC in the common block
      integer idx_x, idx_y, idx_surf
      flatidx = idx_x + (idx_y - 1) * (NVC(idx_surf)+1)
      return
      end function flatidx

      subroutine makesurf_mesh(isurf)
c--------------------------------------------------------------
c     Sets up all stuff for surface ISURF, 
C     using info from configuration input 
C     and the given mesh coordinate array.
c--------------------------------------------------------------
      INCLUDE 'AVL.INC'
      ! input/output
      integer isurf

      ! working variables (AVL original)
      PARAMETER (KCMAX=50,
     &           KSMAX=500)
      REAL CHSIN, CHCOS, CHSINL, CHSINR, CHCOSL, CHCOSR, AINCL, AINCR, 
     &     CHORDL, CHORDR, CLAFL, CLAFR, SLOPEL, SLOPER, DXDX, ZU_L, 
     &     ZL_L, ZU_R, ZL_R, ZL, ZR, SUM, WTOT, ASTRP
      REAL CHSINL_G(NGMAX),CHCOSL_G(NGMAX),
     &     CHSINR_G(NGMAX),CHCOSR_G(NGMAX),
     &     XLED(NDMAX), XTED(NDMAX), GAINDA(NDMAX)
      INTEGER ISCONL(NDMAX), ISCONR(NDMAX)
      
      ! working variables (OptVL additions)
      real m1, m2, m3, f1, f2, fc, dc1, dc2, dc, a1, a2, a3, xptxind1,
     & xptxind2
      real mesh_surf(3,(NVC(isurf)+1)*(NVS(isurf)+1))
      integer idx_vor, idx_strip, idx_sec, idx_dim, idx_coef, idx_x, 
     & idx_node, idx_nodel, idx_noder, idx_node_yp1, idx_node_nx, 
     & idx_node_nx_yp1, idx_y, nx, ny

      ! functions
      integer flatidx

      ! Get data from common block
      nx = NVC(isurf) + 1
      ny = NVS(isurf) + 1

      ! Check MFRST
      if (MFRST(isurf) .eq. 0) then
      print *, "* Provide the index where the mesh begins for surface",
     & isurf 
      end if

      ! Get the mesh for this surface from the the common block
      mesh_surf = MSHBLK(:,MFRST(isurf):MFRST(isurf)+(nx*ny)-1)

      ! Perform input checks from makesurf (section check removed)

      IF(NVC(ISURF).GT.KCMAX) THEN
       WRITE(*,*) '* makesurf_mesh: Array overflow.  Increase KCMAX to',
     &                                                   NVC(ISURF)
       NVC(ISURF) = KCMAX
      ENDIF

      IF(NVS(ISURF).GT.KSMAX) THEN
       WRITE(*,*) '* makesurf_mesh: Array overflow.  Increase KSMAX to', 
     &                                                   NVS(ISURF)
       NVS(ISURF) = KSMAX
      ENDIF

      ! Image flag set to indicate section definition direction
      ! IMAGS= 1  defines edge 1 located at surface root edge 
      ! IMAGS=-1  defines edge 2 located at surface root edge (reflected surfaces)
      IMAGS(ISURF) = 1
      
      ! Start accumulating the element and strip index references
      ! Accumulate the first element in surface
      if (ISURF == 1) then
            IFRST(ISURF) = 1
      else
            IFRST(ISURF) =  IFRST(ISURF-1) +  NK(ISURF-1)*NJ(ISURF-1)        
      endif

      ! Accumulate the first strip in surface
      if (ISURF == 1) then
            JFRST(ISURF) = 1
      else
            JFRST(ISURF) =  JFRST(ISURF-1) +  NJ(ISURF-1)
      endif
      
      ! Set NK from input data (python layer will ensure this is consistent)
      NK(ISURF) = NVC(ISURF)

      ! We need to start counting strips now since it is a global count
      idx_strip = JFRST(ISURF)

      ! Bypass the entire spanwise node generation routine and go straight to store counters
      ! Index of first strip in surface
      ! This is normally used to store the index of each section in AVL
      ! but since we use strips now each is effectively just a section
      ! We assign this variable accordingly so as not to break anything else
      IF (ISURF .EQ. 1) THEN
        ICNTFRST(ISURF) = 1
      ELSE
        ICNTFRST(ISURF) = ICNTFRST(ISURF-1) + NCNTSEC(ISURF-1)
      ENDIF
      ! Number of strips/sections in surface
      NCNTSEC(ISURF) = NSEC(ISURF)
      ! Store the spanwise index of each strip in each surface
      DO ISEC = 1, NSEC(ISURF)
        II = ICNTFRST(ISURF) + (ISEC-1)
        ICNTSEC(II) = idx_strip
      ENDDO


      ! Apply the scaling and translations to the mesh as a whole
      do idx_y = 1,ny
       do idx_x = 1,nx
        do idx_dim = 1,3
            idx_node = flatidx(idx_x, idx_y, isurf)
            mesh_surf(idx_dim,idx_node) = XYZSCAL(idx_dim,isurf)
     &      *mesh_surf(idx_dim,idx_node) + XYZTRAN(idx_dim,isurf)
        end do
       end do
      end do


      ! Setup the strips

      ! Set spanwise elements to 0
      NJ(ISURF) = 0

      ! Check control and design vars
      IF(NCONTROL.GT.NDMAX) THEN
       WRITE(*,*) '*** Too many control variables.  Increase NDMAX to',
     &            NCONTROL
       STOP
      ENDIF

      IF(NDESIGN.GT.NGMAX) THEN
       WRITE(*,*) '*** Too many design variables.  Increase NGMAX to',
     &            NDESIGN
       STOP
      ENDIF

      ! Instead of looping over sections just loop over all strips in the surface
      do ispan = 1,ny-1


      ! Set reference information for the strip
      ! This code was used in the original to loop over strips in a section. 
      ! We will just reuse the variables here 
      idx_y = idx_strip - JFRST(isurf) + 1
      iptl = idx_y
      iptr = idx_y + 1    
      NJ(isurf) = NJ(isurf) + 1


      ! We need to compute the chord and claf values at the left and right edge of the strip
      ! This code was used in the original to interpolate over sections. 
      ! We will just reuse here to interpolate over a strip which is trivial but avoids pointless code rewrites.
      idx_node = flatidx(1,iptl,isurf)
      idx_node_nx = flatidx(nx,iptl,isurf)
      CHORDL = sqrt((mesh_surf(1,idx_node_nx)-mesh_surf(1,idx_node))**2 
     & + (mesh_surf(3,idx_node_nx)-mesh_surf(3,idx_node))**2)
      idx_node = flatidx(1,iptr,isurf)
      idx_node_nx = flatidx(nx,iptr,isurf)
      CHORDR = sqrt((mesh_surf(1,idx_node_nx)-mesh_surf(1,idx_node))**2
     & + (mesh_surf(3,idx_node_nx)-mesh_surf(3,idx_node))**2)
      CLAFL = CLAF(iptl,  isurf)
      CLAFR = CLAF(iptr,isurf)

      ! Linearly interpolate the incidence projections over the STRIP
      AINCL = AINCS(iptl,isurf)*DTR + ADDINC(isurf)*DTR
      AINCR = AINCS(iptr,isurf)*DTR + ADDINC(isurf)*DTR
      CHSINL = CHORDL*SIN(AINCL)
      CHSINR = CHORDR*SIN(AINCR)
      CHCOSL = CHORDL*COS(AINCL)
      CHCOSR = CHORDR*COS(AINCR)

      ! We need to determine which controls belong to this section 
      ! Bring over the routine for this from makesurf but do it for each strip now
      DO N = 1, NCONTROL
      ISCONL(N) = 0
      ISCONR(N) = 0
      DO ISCON = 1, NSCON(iptl,isurf)
      IF(ICONTD(ISCON,iptl,isurf)  .EQ.N) ISCONL(N) = ISCON
      ENDDO
      DO ISCON = 1, NSCON(iptr,isurf)
      IF(ICONTD(ISCON,iptr,isurf).EQ.N) ISCONR(N) = ISCON
      ENDDO
      ENDDO

      ! We need to determine which dvs belong to this strip 
      ! and setup the chord projection gains
      ! Bring over the routine for this from makesurf but setup for strips
      DO N = 1, NDESIGN
      CHSINL_G(N) = 0.
      CHSINR_G(N) = 0.
      CHCOSL_G(N) = 0.
      CHCOSR_G(N) = 0.

      DO ISDES = 1, NSDES(iptl,isurf)
      IF(IDESTD(ISDES,iptl,isurf).EQ.N) THEN
            CHSINL_G(N) =  CHCOSL * GAING(ISDES,iptl,isurf)*DTR
            CHCOSL_G(N) = -CHSINL * GAING(ISDES,iptl,isurf)*DTR
      ENDIF
      ENDDO

      DO ISDES = 1, NSDES(iptr,isurf)
      IF(IDESTD(ISDES,iptr,isurf).EQ.N) THEN
            CHSINR_G(N) =  CHCOSR * GAING(ISDES,iptr,isurf)*DTR
            CHCOSR_G(N) = -CHSINR * GAING(ISDES,iptr,isurf)*DTR
      ENDIF
      ENDDO
      ENDDO


      ! Set the strip geometry data
      ! Note these computations assume the mesh is not necessarily planar
      ! ultimately if/when we flatten the mesh into a planar one we will want
      ! to use the leading edge positions and chords from the original input mesh
      

      ! Strip left side
      idx_node = flatidx(1,idx_y,isurf)
      idx_node_nx = flatidx(nx,idx_y,isurf)
      do idx_dim = 1,3
       RLE1(idx_dim,idx_strip) = mesh_surf(idx_dim,idx_node)
      end do 
      CHORD1(idx_strip) = sqrt((mesh_surf(1,idx_node_nx)
     & -mesh_surf(1,idx_node))**2 + (mesh_surf(3,idx_node_nx)
     & -mesh_surf(3,idx_node))**2)

      ! Strip right side
      idx_node_yp1 = flatidx(1,idx_y+1,isurf)
      idx_node_nx_yp1 = flatidx(nx,idx_y+1,isurf)
      do idx_dim = 1,3
       RLE2(idx_dim,idx_strip) = mesh_surf(idx_dim,idx_node_yp1)
      end do 
      CHORD2(idx_strip) = sqrt((mesh_surf(1,idx_node_nx_yp1)
     & -mesh_surf(1,idx_node_yp1))**2 + (mesh_surf(3,idx_node_nx_yp1)
     & -mesh_surf(3,idx_node_yp1))**2)

      ! Strip mid-point 
      do idx_dim = 1,3
       ! Since the strips are linear SPANWISE we can just interpolate
       RLE(idx_dim,idx_strip) = (RLE1(idx_dim,idx_strip)
     &   + RLE2(idx_dim,idx_strip))/2.
      end do 
       ! The strips are not necessarily linear chord wise but by definition the chord value is
       ! so we can interpolate
       CHORD(idx_strip) = (CHORD1(idx_strip)+CHORD2(idx_strip))/2.

      ! Strip geometric incidence angle at the mid-point
      ! This is strip incidence angle is computed from the LE and TE points
      ! of the given geometry and is completely independent of AINC
      ! This quantity is needed to correctly handle nonplanar meshes and is only needed if the mesh isnt flattened
      GINCSTRIP(idx_strip) = atan2(((mesh_surf(3,idx_node_nx) 
     &  + mesh_surf(3,idx_node_nx_yp1))/2.- (mesh_surf(3,idx_node) + 
     &  mesh_surf(3,idx_node_yp1))/2.),
     & ((mesh_surf(1,idx_node_nx) + mesh_surf(1,idx_node_nx_yp1))/2. 
     &  - (mesh_surf(1,idx_node) + mesh_surf(1,idx_node_yp1))/2.))

      ! Strip width
      m2 = mesh_surf(2,idx_node_yp1)-mesh_surf(2,idx_node)
      m3 = mesh_surf(3,idx_node_yp1)-mesh_surf(3,idx_node)
      WSTRIP(idx_strip) = sqrt(m2**2 + m3**2)

      ! Strip LE and TE sweep slopes
      tanle(idx_strip) = (mesh_surf(1,idx_node_yp1)
     &  -mesh_surf(1,idx_node))/WSTRIP(idx_strip)
      idx_node = flatidx(nx,idx_y,isurf)
      idx_node_yp1 = flatidx(nx,idx_y+1,isurf)
      tante(idx_strip) = (mesh_surf(1,idx_node_yp1)
     &  -mesh_surf(1,idx_node))/WSTRIP(idx_strip)

      ! Compute chord projections and strip twists
      ! In AVL the AINCS are not interpolated. The chord projections are
      ! So we have to replicate this effect.
 
      ! LINEAR interpolation over the strip: left, right, and midpoint
      idx_nodel = flatidx(1,iptl,isurf)
      idx_noder = flatidx(1,iptr,isurf)

!       f1 = (mesh_surf(2,idx_node)-mesh_surf(2,idx_nodel))/
!      & (mesh_surf(2,idx_noder)-mesh_surf(2,idx_nodel))
!       f2 = (mesh_surf(2,idx_node_yp1)-mesh_surf(2,idx_nodel))/
!      & (mesh_surf(2,idx_noder)-mesh_surf(2,idx_nodel))
!       fc = (((mesh_surf(2,idx_node_yp1)+mesh_surf(2,idx_node))/2.) 
!      & -mesh_surf(2,idx_nodel))/(mesh_surf(2,idx_noder)
!      & -mesh_surf(2,idx_nodel))

      ! the above expressions will always evaluate to the following for individual strips
      f1 = 0.0
      f2 = 1.0
      fc = 0.5

      ! Strip left side incidence
      ! CHSIN = CHSINL + f1*(CHSINR-CHSINL)
      ! CHCOS = CHCOSL + f1*(CHCOSR-CHCOSL)
      AINC1(idx_strip) = ATAN2(CHSINL,CHCOSL)

      ! Strip right side incidence
      ! CHSIN = CHSINL + f2*(CHSINR-CHSINL)
      ! CHCOS = CHCOSL + f2*(CHCOSR-CHCOSL)
      AINC2(idx_strip) = ATAN2(CHSINR,CHCOSR)

      ! Strip mid-point incidence
      CHSIN = CHSINL + fc*(CHSINR-CHSINL)
      CHCOS = CHCOSL + fc*(CHCOSR-CHCOSL)
      AINC(idx_strip) = ATAN2(CHSIN,CHCOS)

      ! Set dv gains for incidence angles
      ! Bring over the routine for this from make surf
      DO N = 1, NDESIGN
         CHSIN_G = (1.0-FC)*CHSINL_G(N) + FC*CHSINR_G(N)
         CHCOS_G = (1.0-FC)*CHCOSL_G(N) + FC*CHCOSR_G(N)
         AINC_G(idx_strip,N) = (CHCOS*CHSIN_G - CHSIN*CHCOS_G)
     &                       / (CHSIN**2 + CHCOS**2)
      ENDDO

      ! We have to now setup any control surfaces we defined for this strip
      ! Bring over the routine for this from makesurf but modified for a strip
      DO N = 1, NCONTROL
      ICL = ISCONL(N)
      ICR = ISCONR(N)

      IF(ICL.EQ.0 .OR. ICR.EQ.0) THEN
      ! no control effect here
            GAINDA(N) = 0.
            XLED(N) = 0.
            XTED(N) = 0.

            VHINGE(1,idx_strip,N) = 0.
            VHINGE(2,idx_strip,N) = 0.
            VHINGE(3,idx_strip,N) = 0.

            VREFL(idx_strip,N) = 0.

            PHINGE(1,idx_strip,N) = 0.
            PHINGE(2,idx_strip,N) = 0.
            PHINGE(3,idx_strip,N) = 0.

      ELSE
      ! control variable # N is active here
            GAINDA(N) = GAIND(ICL,iptl  ,isurf)*(1.0-FC)
     &                 + GAIND(ICR,iptr,isurf)*     FC

            ! SAB Note: This interpolation ensures that the hinge line is 
            ! is linear which I think it is an ok assumption for arbitrary wings as long as the user is aware
            ! A curve hinge line could work if needed if we just interpolate XHINGED and scaled by local chord
            XHD = CHORDL*XHINGED(ICL,iptl  ,isurf)*(1.0-FC)
     &           + CHORDR*XHINGED(ICR,iptr,isurf)*     FC
            IF(XHD.GE.0.0) THEN
            ! TE control surface, with hinge at XHD
            XLED(N) = XHD
            XTED(N) = CHORD(idx_strip)
            ELSE
            ! LE control surface, with hinge at -XHD
            XLED(N) =  0.0
            XTED(N) = -XHD
            ENDIF

            VHX = VHINGED(1,ICL,iptl,isurf)*XYZSCAL(1,isurf)
            VHY = VHINGED(2,ICL,iptl,isurf)*XYZSCAL(2,isurf)
            VHZ = VHINGED(3,ICL,iptl,isurf)*XYZSCAL(3,isurf)
            VSQ = VHX**2 + VHY**2 + VHZ**2
            IF(VSQ.EQ.0.0) THEN
            ! default: set hinge vector along hingeline
            ! We are just setting the hinge line across the section
            ! this assumes the hinge is linear even for a nonlinear wing
            VHX = mesh_surf(1,idx_noder)
     &              + ABS(CHORDR*XHINGED(ICR,iptr,isurf))
     &              - mesh_surf(1,idx_nodel)
     &              - ABS(CHORDL*XHINGED(ICL,iptl,isurf))
            VHY = mesh_surf(2,idx_noder)
     &            - mesh_surf(2,idx_nodel)
            VHZ = mesh_surf(3,idx_noder)
     &            - mesh_surf(3,idx_nodel)
            VHX = VHX*XYZSCAL(1,isurf)
            VHY = VHY*XYZSCAL(2,isurf)
            VHZ = VHZ*XYZSCAL(3,isurf)
            VSQ = VHX**2 + VHY**2 + VHZ**2
            ENDIF

            VMOD = SQRT(VSQ)
            VHINGE(1,idx_strip,N) = VHX/VMOD
            VHINGE(2,idx_strip,N) = VHY/VMOD
            VHINGE(3,idx_strip,N) = VHZ/VMOD

            VREFL(idx_strip,N) = REFLD(ICL,iptl, isurf)

            IF(XHD .GE. 0.0) THEN
            PHINGE(1,idx_strip,N) = RLE(1,idx_strip) + XHD
            PHINGE(2,idx_strip,N) = RLE(2,idx_strip)
            PHINGE(3,idx_strip,N) = RLE(3,idx_strip)
            ELSE
            PHINGE(1,idx_strip,N) = RLE(1,idx_strip) - XHD
            PHINGE(2,idx_strip,N) = RLE(2,idx_strip)
            PHINGE(3,idx_strip,N) = RLE(3,idx_strip)
            ENDIF
      ENDIF
      ENDDO      

      ! Interpolate CD-CL polar defining data from input to strips
      DO idx_coef = 1, 6
      CLCD(idx_coef,idx_strip) = (1.0-fc)* 
     & CLCDSEC(idx_coef,iptl,isurf) + 
     & fc*CLCDSEC(idx_coef,iptr,isurf)
      END DO
      ! If the min drag is zero flag the strip as no-viscous data
      LVISCSTRP(idx_strip) = (CLCD(4,idx_strip).NE.0.0)  
      
      
      ! Set the panel (vortex) geometry data

      ! Accumulate the strip element indicies and start counting vorticies
      if (idx_strip .eq. 1) then 
            IJFRST(idx_strip) = 1
      else
            IJFRST(idx_strip) = IJFRST(idx_strip - 1) + 
     &                          NVSTRP(idx_strip - 1)
      endif
      idx_vor = IJFRST(idx_strip)
      NVSTRP(idx_strip) = NVC(isurf)

      ! Associate the strip with the surface
      LSSURF(idx_strip) = isurf

      ! Prepare for cross section interpolation
      NSL = NASEC(iptl  , isurf)
      NSR = NASEC(iptr, isurf)

      ! CHORDC = CHORD(idx_strip)      


      ! Funny story. this original line is now valid now that we interpolate over the strip
      clafc =  (1.-FC)*(CHORDL/CHORD(idx_strip))*CLAFL
     &           +     FC *(CHORDR/CHORD(idx_strip))*CLAFR
      ! Suggestion from Hal Yougren for non linear sections:
      ! clafc =  (1.-fc)*clafl + fc*clafr

      ! loop over vorticies for the strip
      do idx_x = 1, nvc(isurf)
       
       ! Left bound vortex points 
       idx_node = flatidx(idx_x,idx_y,isurf)
      ! Compute the panel left side chord
       dc1 = sqrt((mesh_surf(1,idx_node+1) - mesh_surf(1,idx_node))**2
     &          + (mesh_surf(3,idx_node+1) - mesh_surf(3,idx_node))**2)

       if (LMESHFLAT(isurf)) then
       ! Place vortex at panel quarter chord of the flat mesh
       dx1 = sqrt((mesh_surf(1,idx_node) - RLE1(1,idx_strip))**2
     &          + (mesh_surf(3,idx_node) - RLE1(3,idx_strip))**2)
       RV1(2,idx_vor) = RLE1(2,idx_strip)
       RV1(3,idx_vor) = RLE1(3,idx_strip)
       RV1(1,idx_vor) = RLE1(1,idx_strip) + dx1 + (dc1/4.)

       ! Compute the panel left side angle
       a1 = atan2((mesh_surf(3,idx_node+1) - mesh_surf(3,idx_node)),
     &            (mesh_surf(1,idx_node+1) - mesh_surf(1,idx_node)))
      ! Place vortex at panel quarter chord of the true mesh
       RV1MSH(2,idx_vor) = mesh_surf(2,idx_node)  
       RV1MSH(1,idx_vor) = mesh_surf(1,idx_node) + (dc1/4.)*cos(a1)
       RV1MSH(3,idx_vor) = mesh_surf(3,idx_node) + (dc1/4.)*sin(a1) 
       else
      ! Compute the panel left side angle
       a1 = atan2((mesh_surf(3,idx_node+1) - mesh_surf(3,idx_node)),
     &            (mesh_surf(1,idx_node+1) - mesh_surf(1,idx_node)))
      ! Place vortex at panel quarter chord
       RV1(2,idx_vor) = mesh_surf(2,idx_node)  
       RV1(1,idx_vor) = mesh_surf(1,idx_node) + (dc1/4.)*cos(a1)
       RV1(3,idx_vor) = mesh_surf(3,idx_node) + (dc1/4.)*sin(a1)  
       
       ! Make a copy in the true mesh array for post processing
       RV1MSH(2,idx_vor) = RV1(2,idx_vor)
       RV1MSH(1,idx_vor) = RV1(1,idx_vor)
       RV1MSH(3,idx_vor) = RV1(3,idx_vor)
       end if

       ! Right bound vortex points 
       idx_node_yp1 = flatidx(idx_x,idx_y+1,isurf)
      ! Compute the panel right side chord
       dc2 = sqrt((mesh_surf(1,idx_node_yp1+1) 
     &  - mesh_surf(1,idx_node_yp1))**2 + (mesh_surf(3,idx_node_yp1+1) 
     &  - mesh_surf(3,idx_node_yp1))**2)

      if (LMESHFLAT(isurf)) then
      ! Place vortex at panel quarter chord of the flat mesh
       dx2 = sqrt((mesh_surf(1,idx_node_yp1) - RLE2(1,idx_strip))**2
     & + (mesh_surf(3,idx_node_yp1) - RLE2(3,idx_strip))**2)
       RV2(2,idx_vor) = RLE2(2,idx_strip)
       RV2(3,idx_vor) = RLE2(3,idx_strip)
       RV2(1,idx_vor) = RLE2(1,idx_strip) + dx2 + (dc2/4.)

      ! Compute the panel right side angle
       a2 = atan2((mesh_surf(3,idx_node_yp1+1) -
     & mesh_surf(3,idx_node_yp1)), (mesh_surf(1,idx_node_yp1+1) - 
     & mesh_surf(1,idx_node_yp1)))
      ! Place vortex at panel quarter chord of the true mesh
       RV2MSH(2,idx_vor) = mesh_surf(2,idx_node_yp1)  
       RV2MSH(1,idx_vor) = mesh_surf(1,idx_node_yp1) + (dc2/4.)*cos(a2)
       RV2MSH(3,idx_vor) = mesh_surf(3,idx_node_yp1) + (dc2/4.)*sin(a2)
      else
      ! Compute the panel right side angle
       a2 = atan2((mesh_surf(3,idx_node_yp1+1) -
     & mesh_surf(3,idx_node_yp1)), (mesh_surf(1,idx_node_yp1+1) - 
     & mesh_surf(1,idx_node_yp1)))
      ! Place vortex at panel quarter chord
       RV2(2,idx_vor) = mesh_surf(2,idx_node_yp1)  
       RV2(1,idx_vor) = mesh_surf(1,idx_node_yp1) + (dc2/4.)*cos(a2)
       RV2(3,idx_vor) = mesh_surf(3,idx_node_yp1) + (dc2/4.)*sin(a2)

      ! Make a copy in the true mesh array for post processing
       RV2MSH(2,idx_vor) = RV2(2,idx_vor)
       RV2MSH(1,idx_vor) = RV2(1,idx_vor)
       RV2MSH(3,idx_vor) = RV2(3,idx_vor)    
      end if

      ! Mid-point bound vortex points 
      ! Compute the panel mid-point chord
      ! Panels themselves can never be curved so just interpolate the chord
      ! store as the panel chord in common block
       DXV(idx_vor) = (dc1+dc2)/2.
      ! We need to compute the midpoint angle and panel strip chord projection 
      ! as we need them to compute normals based on the real mesh
       a3 = atan2(((mesh_surf(3,idx_node_yp1+1) 
     &  + mesh_surf(3,idx_node+1))/2.- (mesh_surf(3,idx_node_yp1) + 
     &  mesh_surf(3,idx_node))/2.),
     & ((mesh_surf(1,idx_node_yp1+1) + mesh_surf(1,idx_node+1))/2. 
     &  - (mesh_surf(1,idx_node_yp1) + mesh_surf(1,idx_node))/2.))
      ! project the panel chord onto the strip chord
       DXSTRPV(idx_vor) = DXV(idx_vor)*cos(a3-GINCSTRIP(idx_strip))

       if (LMESHFLAT(isurf)) then
      ! Place vortex at panel quarter chord of the flat mesh
       dx3 = sqrt(((mesh_surf(1,idx_node_yp1)+mesh_surf(1,idx_node))/2 
     &            - RLE(1,idx_strip))**2
     &            + ((mesh_surf(3,idx_node_yp1)+mesh_surf(3,idx_node))/2
     &            - RLE(3,idx_strip))**2)
       RV(2,idx_vor) = RLE(2,idx_strip)
       RV(3,idx_vor) = RLE(3,idx_strip)
       RV(1,idx_vor) = RLE(1,idx_strip) + dx3 + (DXV(idx_vor)/4.)

       ! Place vortex at panel quarter chord of the true mesh
       RVMSH(2,idx_vor) = (mesh_surf(2,idx_node_yp1) 
     &  + mesh_surf(2,idx_node))/2.
       RVMSH(1,idx_vor) = (mesh_surf(1,idx_node_yp1)
     & +mesh_surf(1,idx_node))/2.+ (DXV(idx_vor)/4.)*cos(a3)
       RVMSH(3,idx_vor) = (mesh_surf(3,idx_node_yp1) 
     & +mesh_surf(3,idx_node))/2. + (DXV(idx_vor)/4.)*sin(a3)
      else
      ! Place vortex at panel quarter chord
       RV(2,idx_vor) = (mesh_surf(2,idx_node_yp1) 
     &  + mesh_surf(2,idx_node))/2.
       RV(1,idx_vor) = (mesh_surf(1,idx_node_yp1)
     & +mesh_surf(1,idx_node))/2.+ (DXV(idx_vor)/4.)*cos(a3)
       RV(3,idx_vor) = (mesh_surf(3,idx_node_yp1) 
     & +mesh_surf(3,idx_node))/2. + (DXV(idx_vor)/4.)*sin(a3)

      ! Make a copy in the true mesh array for post processing
       RVMSH(2,idx_vor) = RV(2,idx_vor)
       RVMSH(1,idx_vor) = RV(1,idx_vor)
       RVMSH(3,idx_vor) = RV(3,idx_vor)
      end if


      ! Panel Control points
      ! Y- point 
      ! is just the panel midpoint
       RC(2,idx_vor) = RV(2,idx_vor)
      ! Place the control point at the quarter chord + half chord*clafc
      ! note that clafc is a scaler so is 1. is for 2pi
      ! use data from vortex mid-point computation
       if (LMESHFLAT(isurf)) then
       RC(1,idx_vor) = RV(1,idx_vor) + clafc*(DXV(idx_vor)/2.)
       RC(3,idx_vor) = RV(3,idx_vor)

       RCMSH(1,idx_vor) = RVMSH(1,idx_vor) 
     &  + clafc*(DXV(idx_vor)/2.)*cos(a3)
       RCMSH(3,idx_vor) = RVMSH(3,idx_vor) 
     &  + clafc*(DXV(idx_vor)/2.)*sin(a3)
       RCMSH(2,idx_vor) = RVMSH(2,idx_vor)
       else
       RC(1,idx_vor) = RV(1,idx_vor) + clafc*(DXV(idx_vor)/2.)*cos(a3)
       RC(3,idx_vor) = RV(3,idx_vor) + clafc*(DXV(idx_vor)/2.)*sin(a3)

      ! Make a copy in the true mesh array for post processing
       RCMSH(1,idx_vor) = RC(1,idx_vor)
       RCMSH(3,idx_vor) = RC(3,idx_vor)
       RCMSH(2,idx_vor) = RC(2,idx_vor)
       end if

      ! Source points
      ! Y- point
       RS(2,idx_vor) = RV(2,idx_vor)
      ! Place the source point at the half chord
      ! use data from vortex mid-point computation
      ! add another quarter chord to the quarter chord
       if (LMESHFLAT(isurf)) then
       RS(1,idx_vor) = RV(1,idx_vor) + (DXV(idx_vor)/4.)
       RS(3,idx_vor) = RV(3,idx_vor) + (DXV(idx_vor)/4.)
       else
       RS(1,idx_vor) = RV(1,idx_vor) + (DXV(idx_vor)/4.)*cos(a3) 
       RS(3,idx_vor) = RV(3,idx_vor) + (DXV(idx_vor)/4.)*sin(a3)
       end if


       ! Set the camber slopes for the panel
       
       ! Camber slope at control point
       CALL AKIMA(XASEC(1,iptl,  isurf),SASEC(1,iptl,  isurf),
     &               NSL,(RC(1,idx_vor)-RLE(1,idx_strip))
     &                /CHORD(idx_strip),SLOPEL, DSDX)
       CALL AKIMA(XASEC(1,iptr,isurf),SASEC(1,iptr,isurf),
     &               NSR,(RC(1,idx_vor)-RLE(1,idx_strip))
     &                /CHORD(idx_strip),SLOPER, DSDX)

       ! Alternative for nonlinear sections per Hal Youngren
      ! SLOPEC(idx_vor) =  (1.-fc)*SLOPEL + fc*SLOPER
      ! The original line is valid for interpolation over a strip
       SLOPEC(idx_vor) =  (1.-fc)*(CHORDL/CHORD(idx_strip))*SLOPEL 
     &                    +     fc *(CHORDR/CHORD(idx_strip))*SLOPER

       ! Camber slope at vortex mid-point
       CALL AKIMA(XASEC(1,iptl,  isurf),SASEC(1,iptl,  isurf),
     &               NSL,(RV(1,idx_vor)-RLE(1,idx_strip))
     &               /CHORD(idx_strip),SLOPEL, DSDX)
       CALL AKIMA(XASEC(1,iptr,isurf),SASEC(1,iptr,isurf),
     &               NSR,(RV(1,idx_vor)-RLE(1,idx_strip))
     &               /CHORD(idx_strip),SLOPER, DSDX)

      ! Alternative for nonlinear sections per Hal Youngren
      ! SLOPEV(idx_vor) =  (1.-fc)*SLOPEL + fc*SLOPER
      ! The original line is valid for interpolation over a strip
       SLOPEV(idx_vor) =  (1.-fc)*(CHORDL/CHORD(idx_strip))*SLOPEL 
     &                    + fc *(CHORDR/CHORD(idx_strip))*SLOPER

      ! Associate the panel with strip chord and component
      CHORDV(idx_vor) = CHORD(idx_strip)
      LVCOMP(idx_vor) = LNCOMP(isurf)
      
      ! Enforce no penetration at the control point
      LVNC(idx_vor) = .true.

      ! element inherits alpha,beta flag from surface
      LVALBE(idx_vor) = LFALBE(isurf)

      ! We need to scale the control surface gains by the fraction
      ! of the element on the control surface
      do N = 1, NCONTROL
      !scale control gain by factor 0..1, (fraction of element on control surface)
       xpt = ((mesh_surf(1,idx_node)+mesh_surf(1,idx_node_yp1))
     &         /2 - RLE(1,idx_strip))/CHORD(idx_strip)

       FRACLE = (XLED(N)/CHORD(idx_strip)-xpt) / 
     & (DXV(idx_vor)/CHORD(idx_strip))

       FRACTE = (XTED(N)/CHORD(idx_strip)-xpt) / 
     & (DXV(idx_vor)/CHORD(idx_strip))

       FRACLE = MIN( 1.0 , MAX( 0.0 , FRACLE ) )
       FRACTE = MIN( 1.0 , MAX( 0.0 , FRACTE ) )

       DCONTROL(idx_vor,N) = GAINDA(N)*(FRACTE-FRACLE)
      end do

      ! TE control point used only if surface sheds a wake
      LVNC(idx_vor) = LFWAKE(isurf)

      ! Use the cross sections to generate the OML
      ! nodal grid associated with vortex strip (aft-panel nodes)
      ! NOTE: airfoil in plane of wing, but not rotated perpendicular to dihedral;
      ! retained in (x,z) plane at this point

      ! Store the panel LE mid point for the next panel in the strip
      ! This gets used a lot here 
      ! We use the original input mesh (true mesh) to compute points for the OML
      xptxind1 = ((mesh_surf(1,idx_node+1)+mesh_surf(1,idx_node_yp1+1))
     &         /2 - RLE(1,idx_strip))/CHORD(idx_strip)

!       xptxind2 = (mesh_surf(1,idx_node_yp1+1)
!      &           - RLE2(1,idx_strip))/CHORD2(idx_strip) 

      ! Interpolate cross section on left side
      CALL AKIMA( XLASEC(1,iptl,isurf), ZLASEC(1,iptl,isurf),
     &               NSL,xptxind1, ZL_L, DSDX )
      CALL AKIMA( XUASEC(1,iptl,isurf), ZUASEC(1,iptl,isurf),
     &               NSL,xptxind1, ZU_L, DSDX )

      ! Interpolate cross section on right side
      CALL AKIMA( XLASEC(1,iptr,isurf),
     & ZLASEC(1,iptr,isurf),NSR, xptxind1, ZL_R, DSDX)
                      
      CALL AKIMA( XUASEC(1,iptr,isurf),
     & ZUASEC(1,iptr,isurf),NSR, xptxind1, ZU_R, DSDX)


      ! Compute the left aft node of panel 
      ! X-point
      XYN1(1,idx_vor) = RLE1(1,idx_strip) + 
     &                        xptxind1*CHORD1(idx_strip)

      ! Y-point
      XYN1(2,idx_vor) = RLE1(2,idx_strip)

      ! Interpolate z from sections to left aft node of panel
      ZL =  (1.-f1)*ZL_L + f1 *ZL_R
      ZU =  (1.-f1)*ZU_L + f1 *ZU_R

      ! Store left aft z-point
      ZLON1(idx_vor)  = RLE1(3,idx_strip) + ZL*CHORD1(idx_strip)
      ZUPN1(idx_vor)  = RLE1(3,idx_strip) + ZU*CHORD1(idx_strip)

      ! Compute the right aft node of panel 
      ! X-point
      XYN2(1,idx_vor) = RLE2(1,idx_strip) + 
     &                        xptxind1*CHORD2(idx_strip)

      ! Y-point
      XYN2(2,idx_vor) = RLE2(2,idx_strip)
            
      ! Interpolate z from sections to right aft node of panel
      ZL =  (1.-f2)*ZL_L + f2 *ZL_R
      ZU =  (1.-f2)*ZU_L + f2 *ZU_R

      ! Store right aft z-point
      ZLON2(idx_vor)  = RLE2(3,idx_strip) + ZL*CHORD2(idx_strip)
      ZUPN2(idx_vor)  = RLE2(3,idx_strip) + ZU*CHORD2(idx_strip)
      

      idx_vor = idx_vor + 1
      end do ! End vortex loop
      idx_strip = idx_strip + 1
      end do ! End strip loop

      ! Compute the wetted area and cave from the true mesh
      sum = 0.0
      wtot = 0.0
      DO JJ = 1, NJ(isurf)
        J = JFRST(isurf) + JJ-1 
        ASTRP = WSTRIP(J)*CHORD(J)
        SUM  = SUM + ASTRP
        WTOT = WTOT + WSTRIP(J)
      ENDDO
      SSURF(isurf) = SUM

      IF(WTOT .EQ. 0.0) THEN
       CAVESURF(isurf) = 0.0
      ELSE
       CAVESURF(isurf) = sum/wtot
      ENDIF
      ! add number of strips to the global count
      NSTRIP = NSTRIP + NJ(isurf)
      ! add number of of votrices to the global count
      NVOR = NVOR + NK(isurf)*NJ(isurf)
      
      end subroutine makesurf_mesh
      
      subroutine update_surfaces()
c--------------------------------------------------------------
c     Updates all surfaces, using the stored data.
c     Resets the strips and vorticies so that AVL rebuilds them
c     from the updated geometry. Recomputes panels normals and
c     tells AVL to rebuild the AIC and other aero related data 
c     arrays on the next execution. 
c--------------------------------------------------------------
      use avl_heap_inc
      include 'AVL.INC'
      integer ii
      
      NSTRIP = 0
      NVOR = 0
      
      ISURF = 1
      
      NSURFDUPL = 0
      do ii=1,(NSURF)
            if (ldupl(ii)) NSURFDUPL = NSURFDUPL + 1
      enddo 
      
c     the iterations of this loop are not independent because we count
c     up the size information as we make each surface
      do ii=1,(NSURF-NSURFDUPL)
            if (lverbose) write(*,*) 'Updating surface ',ISURF
            if (lsurfmsh(isurf)) then
                  call makesurf_mesh(ISURF) 
            else
                  call makesurf(ISURF)
            end if
            
            if(ldupl(isurf)) then
                  if (lverbose) write(*,*) ' reduplicating ',ISURF
                  call sdupl(isurf,ydupl(isurf),'ydup')
                  ISURF = ISURF + 1
            endif
            
            ISURF = ISURF + 1
            
      end do 
      
      CALL ENCALC
      
c     reset all the flags related to the analysis pipline      
      LAIC = .FALSE.
      LSRD = .FALSE.
      LVEL = .FALSE.
      LSOL = .FALSE.
      LSEN = .FALSE.

      if (NAIC /= NVOR) then 
            call avlheap_clean()
            call avlheap_diff_clean()
            call avlheap_init(NVOR)
            call avlheap_diff_init(NVOR)
      endif 
      
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
      NSPACE = NVB(IBODY)
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
      XPT(NVB(IBODY)) = 1.0
C
C---- set body nodes and radii
      VOLB = 0.0
      SRFB = 0.0
      DO IVB = 1, NVB(IBODY)
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
ccc        write(43,*) 'NLNODE,RL,RADL ',NLNODE,(RL(K,NLNODE),K=1,3),
ccc     &               RADL(NLNODE)
      ENDDO
C---- get surface length, area and volume
      VOLB = 0.0
      SRFB = 0.0
      XBMN = RL(1,LFRST(IBODY))
      XBMX = XBMN
      DO IVB = 1, NVB(IBODY)-1
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
c      input
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
            if (lverbose) then 
                  write(*,*) 'aifoil Lrange false', isurf, isec
                  write(*,*) (xfmin .gt. 0.01)
                  write(*,*) (xfmax .lt. 0.99)
            endif
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
      LNCOMP(NNI) = LNCOMP(NN)
C
C---- same various logical flags
      LFWAKE(NNI) = LFWAKE(NN)
      LFALBE(NNI) = LFALBE(NN)
      LFLOAD(NNI) = LFLOAD(NN)
      LRANGE(NNI) = LRANGE(NN)
      LSURFSPACING(NNI) = LSURFSPACING(NN)
      LMESHFLAT(NNI) = LMESHFLAT(NN)
      LSURFMSH(NNI) = LSURFMSH(NN)

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
        GINCSTRIP(JJI) = GINCSTRIP(JJ)
        WSTRIP(JJI) =  WSTRIP(JJ)
        TANLE(JJI)  = -TANLE(JJ)
        AINC (JJI)  =  AINC(JJ)
C
cc#ifdef USE_CPOML
        AINC1(JJI) = AINC2(JJ)
        AINC2(JJI) = AINC1(JJ)
C
cc#endif
        LSSURF(idx_strip) = NNI
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
          DXSTRPV(III) = DXSTRPV(II)
          CHORDV(III) = CHORDV(II)
          LVCOMP(III) = LNCOMP(NNI)
          LVALBE(III) = LVALBE(II)
          LVNC(III) = LVNC(II)
          ! Duplicate mesh data if we are using a mesh
          if (lsurfmsh(NN)) then
            RV1MSH(1,III) = RV2MSH(1,II)
            RV1MSH(2,III) = -RV2MSH(2,II) + YOFF
            RV1MSH(3,III) = RV2MSH(3,II)
            RV2MSH(1,III) = RV1MSH(1,II)
            RV2MSH(2,III) = -RV1MSH(2,II) + YOFF
            RV2MSH(3,III) = RV1MSH(3,II)
            RVMSH(1,III) = RVMSH(1,II)
            RVMSH(2,III) = -RVMSH(2,II) + YOFF
            RVMSH(3,III) = RVMSH(3,II)
            RCMSH(1,III) = RCMSH(1,II)
            RCMSH(2,III) = -RCMSH(2,II) + YOFF
            RCMSH(3,III) = RCMSH(3,II)
          end if
          
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
             ! Also checks if surface has been assigned a point cloud mesh
             ! and uses the real mesh to compute normals if it is
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
      real(kind=avl_real) :: dchstrip, DXT, DYT, DZT, ec_msh(3)
C
C...Calculate the normal vector at control points and bound vortex midpoints
C
      DO 10 J = 1, NSTRIP

        ! Since we cannot seperate the encalc routine for direct mesh assignment we have to make it a branch here
        if (lsurfmsh(lssurf(J))) then

        ! Calculate normal vector for the strip (normal to X axis)
        ! we can't just interpolate this anymore given that 
        ! the strip is no longer necessarily linear chordwise

        ! We want the spanwise unit vector for the strip at the 
        ! chordwise location specified by SAXFR (usually set to 0.25)
        ! Loop over all panels in the strip until we find the one that contains
        ! the SAXFR position in it's projected chord. Since the panels themselves are still linear
        ! we can just use the bound vortex unit vector of that panel as 
        ! the spanwise unit vector of the strip at SAXFR

        ! SAB: This is slow, find a better way to do this
        dchstrip = 0.0
        searchSAXFR: do i = IJFRST(J),IJFRST(J) + (NVSTRP(J)-1)
            dchstrip = dchstrip+DXSTRPV(i)
            if (dchstrip .ge. CHORD(J)*SAXFR) then
                  exit searchSAXFR
            end if
        end do searchSAXFR


        ! compute the spanwise unit vector for Vperp def
        DXT =  RV2MSH(1,I)-RV1MSH(1,I)
        DYT =  RV2MSH(2,I)-RV1MSH(2,I)
        DZT =  RV2MSH(3,I)-RV1MSH(3,I)
        XSREF(J) = RVMSH(1,I)
        YSREF(J) = RVMSH(2,I)
        ZSREF(J) = RVMSH(3,I)

        else
        ! original encalc routine for standard AVL geometry
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
        XSREF(J) = (1.0-SAXFR)*AXLE + SAXFR*AXTE
        YSREF(J) = (1.0-SAXFR)*AYLE + SAXFR*AYTE
        ZSREF(J) = (1.0-SAXFR)*AZLE + SAXFR*AZTE
        end if
C
C
        ESS(1,J) =  DXT/SQRT(DXT*DXT + DYT*DYT + DZT*DZT)
        ESS(2,J) =  DYT/SQRT(DXT*DXT + DYT*DYT + DZT*DZT)
        ESS(3,J) =  DZT/SQRT(DXT*DXT + DYT*DYT + DZT*DZT)

        ! Treffz plane normals
        ENSY(J) = -DZT/SQRT(DYT*DYT + DZT*DZT)
        ENSZ(J) =  DYT/SQRT(DYT*DYT + DZT*DZT)

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

          if (lsurfmsh(lssurf(J))) then
          ! Define unit vector along bound leg
          DXB = RV2MSH(1,I)-RV1MSH(1,I) ! right h.v. pt - left h.v. pt 
          DYB = RV2MSH(2,I)-RV1MSH(2,I)
          DZB = RV2MSH(3,I)-RV1MSH(3,I)
          else
C
C...Define unit vector along bound leg
          DXB = RV2(1,I)-RV1(1,I) ! right h.v. pt - left h.v. pt 
          DYB = RV2(2,I)-RV1(2,I)
          DZB = RV2(3,I)-RV1(3,I)
          end if
          EMAG = SQRT(DXB**2 + DYB**2 + DZB**2)
          EB(1) = DXB/EMAG
          EB(2) = DYB/EMAG
          EB(3) = DZB/EMAG
C
C...Define direction of normal vector at control point 
C   The YZ projection of the normal vector matches the camber slope
C   + section local incidence in the YZ defining plane for the section
      ! First start by combining the contributions to the panel 
      ! incidence from AVL incidence and camberline slope variables
      ! these are not actual geometric transformations of the mesh
      ! but rather further modifications to the chordwise vector that 
      ! will get used to compute normals
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

          if (lsurfmsh(lssurf(J))) then
          ! direct mesh assignemnt branch
          ! now we compute the chordwise panel vector
          ! note that panel's chordwise vector has contributions
          ! from both the geometry itself and the incidence modification
          ! from the AVL AINC and camber slope variables

          ! Get the geometric chordwise vector using RVMSH and RCMSH which should
          ! be located in the same plane given that each individual panel is a 
          ! plane

          EMAG = SQRT((RCMSH(1,I)-RVMSH(1,I))**2 
     &  + (RCMSH(2,I)-RVMSH(2,I))**2 
     &  + (RCMSH(3,I)-RVMSH(3,I))**2)
          ec_msh(1) = (RCMSH(1,I)-RVMSH(1,I))/EMAG
          ec_msh(2) = (RCMSH(2,I)-RVMSH(2,I))/EMAG
          ec_msh(3) = (RCMSH(3,I)-RVMSH(3,I))/EMAG

          ! Now we have to rotate this vector by the incidence contribution from AINC and CAMBER
          ! However, this rotation needs to be done about the local y-axis of the wing 
          ! Earlier we computed ES the normal vector of the strip projected to the Trefftz plane
          ! The axis we need to rotate about is the one purpendicular to this ES.
          ! As a result all panel normals in a given strip will be rotated about the same axis defined by the that strip
          ! The components of the rotation axis are obtained from ES as follows
          ! rot_axis(1) = 0
          ! rot_axis(2) = -ES(3)
          ! rot_axis(3) = ES(2)
          ! We can then multiply ec_msh by the rotation matrix for a rotation about an arbitrary axis
          ! see https://pubs.aip.org/aapt/ajp/article/44/1/63/1050167/Formalism-for-the-rotation-matrix-of-rotations
          ! Note that standard AVL also does this exact same thing but since they always rotate the vector [1,0,0]
          ! the result collapses into the ridiculously simple expression for EC that you see in the other branch

          EC(1) =  COSC*ec_msh(1) + ES(2)*SINC*ec_msh(2) 
     &     + ES(3)*SINC*ec_msh(3)
          EC(2) = -ES(2)*SINC + ((ES(3)**2)*(1-COSC)+COSC)*ec_msh(2)
     &     - (ES(2)*ES(3)*(1-COSC))*ec_msh(3)  
          EC(3) = -ES(3)*SINC*ec_msh(1) - 
     &     (ES(2)*ES(3)*(1-COSC))*ec_msh(2) + 
     &    ((ES(2)**2)*(1-COSC) + COSC)*ec_msh(3)

          else
          EC(1) =  COSC
          EC(2) = -SINC*ES(2)
          EC(3) = -SINC*ES(3)
          end if

          DO N = 1, NDESIGN
            ! The derivative here also changes if we use a custom mesh
            ! Note the derivative is only wrt to AVL incidence vars
            ! as those are the vars AVL DVs can support
            if (lsurfmsh(lssurf(J))) then
            EC(1) =  -SINC*ec_msh(1) + ES(2)*COSC*ec_msh(2) 
     &     + ES(3)*COSC*ec_msh(3)
          EC(2) = -ES(2)*COSC + ((ES(3)**2)*(1+SINC)-SINC)*ec_msh(2)
     &     - (ES(2)*ES(3)*(1+SINC))*ec_msh(3)  
          EC(3) = -ES(3)*COSC*ec_msh(1) - 
     &     (ES(2)*ES(3)*(1+SINC))*ec_msh(2) + 
     &    ((ES(2)**2)*(1+SINC) - SINC)*ec_msh(3)

            else
            EC_G(1,N) = -SINC      *AINC_G(J,N)
            EC_G(2,N) = -COSC*ES(2)*AINC_G(J,N)
            EC_G(3,N) = -COSC*ES(3)*AINC_G(J,N)
            end if
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
      ! This section is identical to the normal vector at the control
      ! point. The only different is that the AVL camberline slope 
      ! is taken at the bound vortex point rather than the control point
      ! the geometric contributions to the normal vector at both of these
      ! point is identical as the lie in the plane of the same panel.
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
          if (lsurfmsh(lssurf(J))) then
          ! direct mesh assignment branch
          ! see explanation in section above for control point normals
          ! ec_msh was already computed in that section
          EC(1) =  COSC*ec_msh(1) + ES(2)*SINC*ec_msh(2) 
     &     + ES(3)*SINC*ec_msh(3)
          EC(2) = -ES(2)*SINC + ((ES(3)**2)*(1-COSC)+COSC)*ec_msh(2)
     &     - (ES(2)*ES(3)*(1-COSC))*ec_msh(3)  
          EC(3) = -ES(3)*SINC*ec_msh(1) - 
     &     (ES(2)*ES(3)*(1-COSC))*ec_msh(2) + 
     &    ((ES(2)**2)*(1-COSC) + COSC)*ec_msh(3)

          else
          EC(1) =  COSC
          EC(2) = -SINC*ES(2)
          EC(3) = -SINC*ES(3)
          end if

          DO N = 1, NDESIGN
            if (lsurfmsh(lssurf(J))) then
            ! Direct mesh assignment branch
            EC(1) =  -SINC*ec_msh(1) + ES(2)*COSC*ec_msh(2) 
     &     + ES(3)*COSC*ec_msh(3)
          EC(2) = -ES(2)*COSC + ((ES(3)**2)*(1+SINC)-SINC)*ec_msh(2)
     &     - (ES(2)*ES(3)*(1+SINC))*ec_msh(3)  
          EC(3) = -ES(3)*COSC*ec_msh(1) - 
     &     (ES(2)*ES(3)*(1+SINC))*ec_msh(2) + 
     &    ((ES(2)**2)*(1+SINC) - SINC)*ec_msh(3)

            else
            EC_G(1,N) = -SINC      *AINC_G(J,N)
            EC_G(2,N) = -COSC*ES(2)*AINC_G(J,N)
            EC_G(3,N) = -COSC*ES(3)*AINC_G(J,N)
            end if
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
          ! this is a pure rotation of the normal vector
          ! the geometric contribution from the mesh is already accounted for
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
