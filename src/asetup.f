C***********************************************************************
C    Module:  asetup.f
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

      SUBROUTINE SETUP
C
C...PURPOSE  To set up the vortex lattice calculation.
C
C            Additional geometry data is calculated.
C            An AIC matrix is assembled and solutions
C            obtained for 0 deg and 90 deg angle of attack 
C            (for later superposition in AERO). The matrix
C            defining normalwash at the bound vortex midpoints
C            is calculated. This is used to define the bound vortex
C            normalwash for 0 and 90 deg angles of attack (used 
C            for superposition in AERO).
C            
C...INPUT    Global Data in labelled commons, defining configuration
C          
C...OUTPUT   GAM(i,1..3)  Vortex strengths for unit VINF1,2,3
C            GAM(i,4..6)  Vortex strengths for unit WROT1,2,3
C            W            Normalwash at vortex midpoints for GAM
C
C...COMMENTS   
C
      use avl_heap_inc  
      INCLUDE 'AVL.INC'
      REAL WORK(NVMAX)
      real t0, t1, t2, t3, t4, t5
C
      AMACH = MACH
      BETM = SQRT(1.0 - AMACH**2)
C
C
      if (ltiming) then 
        call cpu_time(t0)
      end if
      CALL build_AIC
!       IF(.NOT.LAIC) THEN
!         CALL build_AIC
!         if (ltiming) then 
!           call cpu_time(t2)
!           write(*,*) '  nowake time: ', t2 - t1
!         end if 

CC...Holdover from HPV hydro project for forces near free surface
CC...Eliminates excluded vortices from eqns which are below z=Zsym 
C      CALL MUNGEA
C
!       ENDIF
      if (ltiming) then 
        call cpu_time(t3)
        write(*,*) '  LUDCMP time: ', t3 - t0
      end if
C
C
      IF(.NOT.LSRD) THEN
        if(lverbose) then
          WRITE(*,*) ' Building source+doublet strength AIC matrix...'
        end if
        CALL SRDSET(BETM,XYZREF,IYSYM,
     &              NBODY,LFRST,NLMAX,
     &              NL,RL,RADL,
     &              SRC_U,DBL_U)
        if(lverbose) then
          WRITE(*,*) ' Building source+doublet velocity AIC matrix...'
        end if 
        NU = 6
        CALL VSRD(BETM,IYSYM,YSYM,IZSYM,ZSYM,SRCORE,
     &            NBODY,LFRST,NLMAX,
     &            NL,RL,RADL,
     &            NU,SRC_U,DBL_U,
     &            NVOR,RC,
     &            WCSRD_U,NVMAX)
        LSRD = .TRUE.
      ENDIF
      if (ltiming) then 
        call cpu_time(t4)
        write(*,*) '  s+doub time: ', t4 - t3
      end if

C
      IF(.NOT.LVEL) THEN
      if(lverbose) then
        WRITE(*,*) ' Building bound-vortex velocity matrix...'
       end if 
       CALL VVOR(BETM,IYSYM,YSYM,IZSYM,ZSYM,
     &           VRCOREC,VRCOREW,
     &           NVOR,RV1,RV2,LVCOMP,CHORDV,
     &           NVOR,RV ,    LVCOMP,.TRUE.,
     &           WV_GAM,NVMAX)
C
       NU = 6
       CALL VSRD(BETM,IYSYM,YSYM,IZSYM,ZSYM,SRCORE,
     &           NBODY,LFRST,NLMAX,
     &           NL,RL,RADL,
     &           NU,SRC_U,DBL_U,
     &           NVOR,RV ,
     &           WVSRD_U,NVOR)
C
       LVEL = .TRUE.
      ENDIF
C
      if (ltiming) then 
        call cpu_time(t5)
        write(*,*) '  bo vel mat time: ', t5 - t4
      end if
      RETURN
      END ! SETUP


      SUBROUTINE build_AIC
      use avl_heap_inc  
      INCLUDE 'AVL.INC'

      AMACH = MACH
      BETM = SQRT(1.0 - AMACH**2)
      
      
        if (lverbose) then
          WRITE(*,*) ' Building normalwash AIC matrix...'
        end if
        
       CALL VVOR(BETM,IYSYM,YSYM,IZSYM,ZSYM,
     &           VRCOREC,VRCOREW,
     &           NVOR,RV1,RV2,LVCOMP,CHORDV,
     &           NVOR,RC ,    LVCOMP,.FALSE.,
     &           WC_GAM,NVMAX)
C
C$AD II-LOOP
      DO J = 1, NVOR
C$AD II-LOOP
        DO I = 1, NVOR
           AICN(I,J) = WC_GAM(1,I,J)*ENC(1,I)
     &               + WC_GAM(2,I,J)*ENC(2,I)
     &               + WC_GAM(3,I,J)*ENC(3,I)
           LVNC(I) = .TRUE.
         ENDDO
       ENDDO

C----- process each surface which does not shed a wake
C$AD II-LOOP
       DO 10 N = 1, NSURF
         IF(LFWAKE(N)) GO TO 10

C------- go over TE control points on this surface
         J1 = JFRST(N)
         JN = JFRST(N) + NJ(N)-1
C 
C$AD II-LOOP
         DO J = J1, JN
           I1 = IJFRST(J)
           IV = IJFRST(J) + NVSTRP(J) - 1

C--------- clear system row for TE control point
C$AD II-LOOP
           DO JV = 1, NVOR
             AICN(IV,JV) = 0.
           ENDDO
           LVNC(IV) = .FALSE.

C--------- set  sum_strip(Gamma) = 0  for this strip
C$AD II-LOOP
           DO JV = I1, IV
             AICN(IV,JV) = 1.0
           ENDDO
         ENDDO
 10    CONTINUE
       end ! build_AIC
       
      SUBROUTINE factor_AIC
       use avl_heap_inc  
       INCLUDE 'AVL.INC'
       REAL WORK(NVMAX)
       integer i, j
       
       if(lverbose) then
        WRITE(*,*) ' Factoring normalwash AIC matrix...'
       end if 
       do j = 1, NVOR
         do i = 1, NVOR
            AICN_LU(i,j) = AICN(i,j)
            ! write(*,*) i,j, 'AICN(i,j)', AICN(i,j)
         enddo 
        ENDDO
       CALL LUDCMP(NVMAX,NVOR,AICN_LU,IAPIV,WORK)
C
       LAIC = .TRUE.
      END ! factor_AIC
 
      SUBROUTINE GUCALC
      use avl_heap_inc  
      INCLUDE 'AVL.INC'
      REAL RROT(3), VUNIT(3), WUNIT(3)
C
C---- setup BC's at control points for unit freestream,rotation vectors,
C      and back-substitute to obtain corresponding vortex circulations
C
C---- go over freestream velocity components u,v,w = f(V,beta,alpha)
      DO 10 IU = 1, 3
C------ go over all control points
        DO I = 1, NVOR
          IF(LVNC(I)) THEN
C--------- this c.p. has V.n equation...
           VUNIT(1) = 0.
           VUNIT(2) = 0.
           VUNIT(3) = 0.
           IF(LVALBE(I)) THEN
C---------- direct freestream influence is enabled for this c.p.
            VUNIT(IU) = VUNIT(IU) + 1.0
           ENDIF

C--------- always add on indirect freestream influence via BODY sources and doublets
           VUNIT(1) = VUNIT(1) + WCSRD_U(1,I,IU)
           VUNIT(2) = VUNIT(2) + WCSRD_U(2,I,IU)
           VUNIT(3) = VUNIT(3) + WCSRD_U(3,I,IU)

C--------- set r.h.s. for V.n equation
           GAM_U_0(I,IU) = -DOT(ENC(1,I),VUNIT)
           DO N = 1, NCONTROL
             GAM_U_D(I,IU,N) = -DOT(ENC_D(1,I,N),VUNIT)
           ENDDO
           DO N = 1, NDESIGN
             GAM_U_G(I,IU,N) = -DOT(ENC_G(1,I,N),VUNIT)
           ENDDO

          ELSE
C--------- this c.p. does not satisfy V.n equation...
C--------- just clear r.h.s.
           GAM_U_0(I,IU) = 0.
           DO N = 1, NCONTROL
             GAM_U_D(I,IU,N) = 0.
           ENDDO
           DO N = 1, NDESIGN
             GAM_U_G(I,IU,N) = 0.
           ENDDO
          ENDIF
        ENDDO

        CALL BAKSUB(NVMAX,NVOR,AICN_LU,IAPIV,GAM_U_0(1,IU))
        DO N = 1, NCONTROL
          CALL BAKSUB(NVMAX,NVOR,AICN_LU,IAPIV,GAM_U_D(1,IU,N))
        ENDDO
        DO N = 1, NDESIGN
          CALL BAKSUB(NVMAX,NVOR,AICN_LU,IAPIV,GAM_U_G(1,IU,N))
        ENDDO
 10   CONTINUE
C
C---- go over freestream rotation components p,q,r
      DO 20 IU = 4, 6
C------ go over all control points
        DO I = 1, NVOR
          IF(LVNC(I)) THEN
C--------- this c.p. has V.n equation
           WUNIT(1) = 0.
           WUNIT(2) = 0.
           WUNIT(3) = 0.
           IF(LVALBE(I)) THEN
C---------- direct freestream influence is enabled for this c.p.
            WUNIT(IU-3) = WUNIT(IU-3) + 1.0
           ENDIF
           RROT(1) = RC(1,I) - XYZREF(1)
           RROT(2) = RC(2,I) - XYZREF(2)
           RROT(3) = RC(3,I) - XYZREF(3)
           CALL CROSS(RROT,WUNIT,VUNIT)

C--------- always add on indirect freestream influence via BODY sources and doublets
           VUNIT(1) = VUNIT(1) + WCSRD_U(1,I,IU)
           VUNIT(2) = VUNIT(2) + WCSRD_U(2,I,IU)
           VUNIT(3) = VUNIT(3) + WCSRD_U(3,I,IU)

C--------- set r.h.s. for V.n equation
           GAM_U_0(I,IU) = -DOT(ENC(1,I),VUNIT)
           DO N = 1, NCONTROL
             GAM_U_D(I,IU,N) = -DOT(ENC_D(1,I,N),VUNIT)
           ENDDO
           DO N = 1, NDESIGN
             GAM_U_G(I,IU,N) = -DOT(ENC_G(1,I,N),VUNIT)
           ENDDO
          ELSE
C--------- just clear r.h.s.
           GAM_U_0(I,IU) = 0.
           DO N = 1, NCONTROL
             GAM_U_D(I,IU,N) = 0.
           ENDDO
           DO N = 1, NDESIGN
             GAM_U_G(I,IU,N) = 0.
           ENDDO
          ENDIF
        ENDDO
        CALL BAKSUB(NVMAX,NVOR,AICN_LU,IAPIV,GAM_U_0(1,IU))
        DO N = 1, NCONTROL
          CALL BAKSUB(NVMAX,NVOR,AICN_LU,IAPIV,GAM_U_D(1,IU,N))
        ENDDO
        DO N = 1, NDESIGN
          CALL BAKSUB(NVMAX,NVOR,AICN_LU,IAPIV,GAM_U_G(1,IU,N))
        ENDDO
   20 CONTINUE
C
      RETURN
      END ! GUCALC


 
      SUBROUTINE GDCALC(NQDEF,LQDEF,ENC_Q,GAM_Q)
      use avl_heap_inc  
      INCLUDE 'AVL.INC'
      LOGICAL LQDEF(*)
      REAL ENC_Q(3,NVMAX,*), GAM_Q(NVMAX,*)
C
C
      REAL RROT(3), VROT(3), VQ(3)
C
      IF(NQDEF.EQ.0) RETURN
C
C---- Setup variational BC's at the control points
      DO 100 IQ = 1, NQDEF
C------ don't bother if this control variable is undefined
        IF(.NOT.LQDEF(IQ)) GO TO 100
C
        DO I = 1, NVOR
          IF(LVNC(I)) THEN
           IF(LVALBE(I)) THEN
            RROT(1) = RC(1,I) - XYZREF(1)
            RROT(2) = RC(2,I) - XYZREF(2)
            RROT(3) = RC(3,I) - XYZREF(3)
            CALL CROSS(RROT,WROT,VROT)
            DO K = 1, 3
              VQ(K) = VINF(K)
     &              + VROT(K)
            ENDDO
           ELSE
            VQ(1) = 0.
            VQ(2) = 0.
            VQ(3) = 0.
           ENDIF

           DO K = 1, 3
             VQ(K) = VQ(K)
     &             + WCSRD_U(K,I,1)*VINF(1)
     &             + WCSRD_U(K,I,2)*VINF(2)
     &             + WCSRD_U(K,I,3)*VINF(3)
     &             + WCSRD_U(K,I,4)*WROT(1)
     &             + WCSRD_U(K,I,5)*WROT(2)
     &             + WCSRD_U(K,I,6)*WROT(3)
           ENDDO
C
           GAM_Q(I,IQ) = -DOT(ENC_Q(1,I,IQ),VQ)
          ELSE
           GAM_Q(I,IQ) = 0.
          ENDIF
        ENDDO

C********************************************************************
C...Holdover from HPV hydro project for forces near free surface
C...Eliminates excluded vortex equations for strips with z<Zsym 
ccc      CALL MUNGEB(GAM_Q(1,IQ))
C********************************************************************
C
        CALL BAKSUB(NVMAX,NVOR,AICN_LU,IAPIV,GAM_Q(1,IQ))
 100  CONTINUE
C
      RETURN
      END ! GDCALC



      SUBROUTINE MUNGEA
C
C...PURPOSE  To remove hidden vortex equations in AIC matrix
C          
C...OUTPUT   A(.,.)  AIC matrix with affected rows replaced with 1 on diagonal
C                
C
      use avl_heap_inc  
      INCLUDE 'AVL.INC'
C
      DO 30 J = 1, NSTRIP
        IF (.NOT. LSTRIPOFF(J)) GO TO 30
        I1 = IJFRST(J)
        DO 20 K = 1, NVSTRP(J)
          II = I1+K-1
          DO 10 I = 1, NVOR
            AICN(II,I) = 0.0
   10     CONTINUE
          AICN(II,II) = 1.0
   20   CONTINUE
   30 CONTINUE
C
      RETURN
      END


      SUBROUTINE MUNGEB(B)
C
C...PURPOSE  To remove hidden vortex equations in RHS's
C          
C...OUTPUT   B(.)  RHS vector with affected rows replaced with 0
C
      use avl_heap_inc  
      INCLUDE 'AVL.INC'
      REAL B(NVMAX)
C
      DO 30 J = 1, NSTRIP
        IF (.NOT. LSTRIPOFF(J)) GO TO 30
        I1 = IJFRST(J)
        DO 20 K = 1, NVSTRP(J)
          II = I1+K-1
          B(II) = 0.
   20   CONTINUE
   30 CONTINUE
C
      RETURN
      END


      SUBROUTINE GAMSUM
      INCLUDE 'AVL.INC'
C--------------------------------------------------
C     Sums AIC components to get GAM, SRC, DBL
C--------------------------------------------------
C
C---- Set vortex strengths
      DO I = 1, NVOR
        DO IU = 1, 6
          GAM_U(I,IU) = GAM_U_0(I,IU)

          DO N = 1, NCONTROL
            GAM_U(I,IU) = GAM_U(I,IU) + GAM_U_D(I,IU,N)*DELCON(N)
          ENDDO
          DO N = 1, NDESIGN
            GAM_U(I,IU) = GAM_U(I,IU) + GAM_U_G(I,IU,N)*DELDES(N)
          ENDDO
        ENDDO

        DO N = 1, NCONTROL
          GAM_D(I,N) = GAM_U_D(I,1,N)*VINF(1)
     &               + GAM_U_D(I,2,N)*VINF(2)
     &               + GAM_U_D(I,3,N)*VINF(3)
     &               + GAM_U_D(I,4,N)*WROT(1)
     &               + GAM_U_D(I,5,N)*WROT(2)
     &               + GAM_U_D(I,6,N)*WROT(3)
        ENDDO

        DO N = 1, NDESIGN
          GAM_G(I,N) = GAM_U_G(I,1,N)*VINF(1)
     &               + GAM_U_G(I,2,N)*VINF(2)
     &               + GAM_U_G(I,3,N)*VINF(3)
     &               + GAM_U_G(I,4,N)*WROT(1)
     &               + GAM_U_G(I,5,N)*WROT(2)
     &               + GAM_U_G(I,6,N)*WROT(3)
        ENDDO

        GAM(I) = GAM_U(I,1)*VINF(1)
     &         + GAM_U(I,2)*VINF(2)
     &         + GAM_U(I,3)*VINF(3)
     &         + GAM_U(I,4)*WROT(1)
     &         + GAM_U(I,5)*WROT(2)
     &         + GAM_U(I,6)*WROT(3)

c        DO N = 1, NCONTROL
c          GAM(I) = GAM(I) + GAM_D(I,N)*DELCON(N)
c        ENDDO
c        DO N = 1, NDESIGN
c          GAM(I) = GAM(I) + GAM_G(I,N)*DELDES(N)
c        ENDDO
      END DO
C
C---- Set source and doublet strengths
      DO L = 1, NLNODE
        SRC(L) = SRC_U(L,1)*VINF(1)
     &         + SRC_U(L,2)*VINF(2)
     &         + SRC_U(L,3)*VINF(3)
     &         + SRC_U(L,4)*WROT(1)
     &         + SRC_U(L,5)*WROT(2)
     &         + SRC_U(L,6)*WROT(3)
        DO K = 1, 3
          DBL(K,L) = DBL_U(K,L,1)*VINF(1)
     &             + DBL_U(K,L,2)*VINF(2)
     &             + DBL_U(K,L,3)*VINF(3)
     &             + DBL_U(K,L,4)*WROT(1)
     &             + DBL_U(K,L,5)*WROT(2)
     &             + DBL_U(K,L,6)*WROT(3)
        ENDDO
      ENDDO
C
      RETURN
      END ! GAMSUM


      SUBROUTINE VELSUM()
      use avl_heap_inc  
      INCLUDE 'AVL.INC'
      
C---------------------------------------------------------
C     Sums AIC components to get induced velocities and
C     velocities per unit freestream VXX_U,WXX_U
C
C     VC, VV    .h.v. induced at control pts and vortex pts
C     WC, WV    total induced at control pts and vortex pts
C     WCSRD, WVSRD  body source/doublet induced velocities
C---------------------------------------------------------

C
      DO I = 1, NVOR
        DO K = 1, 3
          VC(K,I) = 0.0
          VV(K,I) = 0.0
C--- h.v. velocity at control points and vortex midpoints
          DO J = 1, NVOR
            VC(K,I) = VC(K,I) + WC_GAM(K,I,J)*GAM(J)
            VV(K,I) = VV(K,I) + WV_GAM(K,I,J)*GAM(J)
          ENDDO
          DO N = 1, NUMAX
            VC_U(K,I,N) = 0.0
            VV_U(K,I,N) = 0.0
            DO J = 1, NVOR
              VC_U(K,I,N) = VC_U(K,I,N) + WC_GAM(K,I,J)*GAM_U(J,N)
              VV_U(K,I,N) = VV_U(K,I,N) + WV_GAM(K,I,J)*GAM_U(J,N)
            ENDDO
          ENDDO
          DO N = 1, NCONTROL
            VC_D(K,I,N) = 0.0
            VV_D(K,I,N) = 0.0
            DO J = 1, NVOR
              VC_D(K,I,N) = VC_D(K,I,N) + WC_GAM(K,I,J)*GAM_D(J,N)
              VV_D(K,I,N) = VV_D(K,I,N) + WV_GAM(K,I,J)*GAM_D(J,N)
            ENDDO
          ENDDO
          DO N = 1, NDESIGN
            VC_G(K,I,N) = 0.0
            VV_G(K,I,N) = 0.0
            DO J = 1, NVOR
              VC_G(K,I,N) = VC_G(K,I,N) + WC_GAM(K,I,J)*GAM_G(J,N)
              VV_G(K,I,N) = VV_G(K,I,N) + WV_GAM(K,I,J)*GAM_G(J,N)
            ENDDO
          ENDDO
C--- velocity contribution from body sources and doublets
          WCSRD(K,I) = WCSRD_U(K,I,1)*VINF(1)
     &               + WCSRD_U(K,I,2)*VINF(2)
     &               + WCSRD_U(K,I,3)*VINF(3)
     &               + WCSRD_U(K,I,4)*WROT(1)
     &               + WCSRD_U(K,I,5)*WROT(2)
     &               + WCSRD_U(K,I,6)*WROT(3)
          WVSRD(K,I) = WVSRD_U(K,I,1)*VINF(1)
     &               + WVSRD_U(K,I,2)*VINF(2)
     &               + WVSRD_U(K,I,3)*VINF(3)
     &               + WVSRD_U(K,I,4)*WROT(1)
     &               + WVSRD_U(K,I,5)*WROT(2)
     &               + WVSRD_U(K,I,6)*WROT(3)
C--- total velocity at control points and vortex midpoints
          WC(K,I) = VC(K,I) + WCSRD(K,I)
          WV(K,I) = VV(K,I) + WVSRD(K,I)
          DO N = 1, NUMAX
            WC_U(K,I,N) = VC_U(K,I,N) + WCSRD_U(K,I,N)
            WV_U(K,I,N) = VV_U(K,I,N) + WVSRD_U(K,I,N)
          ENDDO
          DO N = 1, NDMAX
            WC_D(K,I,N) = VC_D(K,I,N)
            WV_D(K,I,N) = VV_D(K,I,N)
          ENDDO
          DO N = 1, NGMAX
            WC_G(K,I,N) = VC_G(K,I,N)
            WV_G(K,I,N) = VV_G(K,I,N)
          ENDDO
C
        ENDDO
      ENDDO
C
      RETURN
      END ! VELSUM
      
      subroutine set_par_and_cons(NITER, IR)
      include 'AVL.INC'
      integer NITER, IR
      
      ! mach setting is done in here. Should it be?
      call set_params(IR)
      ! Additionally set the reference point to be at the cg
      XYZREF(1) = PARVAL(IPXCG,IR)
      XYZREF(2) = PARVAL(IPYCG,IR)
      XYZREF(3) = PARVAL(IPZCG,IR)
C
      CDREF = PARVAL(IPCD0,IR)
C
C
      IF(MACH.NE.AMACH) THEN
C----- new Mach number invalidates close to everything that's stored
       LAIC = .FALSE.
       LSRD = .FALSE.
       LVEL = .FALSE.
       LSOL = .FALSE.
       LSEN = .FALSE.
       LOBAIC = .FALSE.
       LOBVEL = .FALSE.
      ENDIF
      
      IF(NITER.GT.0) THEN
C----- might as well directly set operating variables if they are known
       IF(ICON(IVALFA,IR).EQ.ICALFA) ALFA    = CONVAL(ICALFA,IR)*DTR
       IF(ICON(IVBETA,IR).EQ.ICBETA) BETA    = CONVAL(ICBETA,IR)*DTR
       IF(ICON(IVROTX,IR).EQ.ICROTX) WROT(1) = CONVAL(ICROTX,IR)*2./BREF
       IF(ICON(IVROTY,IR).EQ.ICROTY) WROT(2) = CONVAL(ICROTY,IR)*2./CREF
       IF(ICON(IVROTZ,IR).EQ.ICROTZ) WROT(3) = CONVAL(ICROTZ,IR)*2./BREF
      
C$AD II-LOOP
      DO N = 1, NDMAX
        IV = IVTOT + N
        IC = ICTOT + N
        IF(ICON(IV,IR) .eq. IC) DELCON(N) = CONVAL(IC,IR)
      ENDDO
      
      ENDIF
      end subroutine set_par_and_cons
        
      subroutine set_vel_rhs
      include 'AVL.INC'
      real rrot(3), vunit(3), VUNIT_W_term(3),  wunit(3)
      
C$AD II-LOOP
      DO I = 1, NVOR
        IF(LVNC(I)) THEN
          VUNIT(1) = 0.
          VUNIT(2) = 0.
          VUNIT(3) = 0.
          
          WUNIT(1) = 0.
          WUNIT(2) = 0.
          WUNIT(3) = 0.
           IF(LVALBE(I)) THEN
            VUNIT(1) = VINF(1)
            VUNIT(2) = VINF(2)
            VUNIT(3) = VINF(3)
            
            WUNIT(1) = WROT(1)
            WUNIT(2) = WROT(2)
            WUNIT(3) = WROT(3)
            
           ENDIF
           
           RROT(1) = RC(1,I) - XYZREF(1)
           RROT(2) = RC(2,I) - XYZREF(2)
           RROT(3) = RC(3,I) - XYZREF(3)
           CALL CROSS(RROT,WUNIT,VUNIT_W_term)
           
           VUNIT = VUNIT + VUNIT_W_term
           
           RHS(I) = -DOT(ENC(1,I),VUNIT)
            
           ! Add contribution from control surfaces
C$AD II-LOOP
          DO N = 1, NCONTROL
            RHS(I) = RHS(I) -DOT(ENC_D(1,I,N),VUNIT)*DELCON(N)
          ENDDO
        ELSE
          RHS(I) = 0
        endif
        
      enddo

        
      end !set_vel_rhs

      subroutine set_vel_rhs_u(IU)
        include 'AVL.INC'
        real rrot(3), vunit(3), VUNIT_W_term(3),  wunit(3)
        
C$AD II-LOOP
        DO I = 1, NVOR
          IF(LVNC(I)) THEN
            VUNIT(1) = 0.
            VUNIT(2) = 0.
            VUNIT(3) = 0.
            
            WUNIT(1) = 0.
            WUNIT(2) = 0.
            WUNIT(3) = 0.
             IF(LVALBE(I)) THEN
              if (IU .le. 3) then
                VUNIT(IU) = VUNIT(IU) + 1.0
              else if (IU .le. 6) then
                WUNIT(IU-3) = WUNIT(IU-3) + 1.0
              endif
             ENDIF
C--------- always add on indirect freestream influence via BODY sources and doublets
             VUNIT(1) = VUNIT(1) + WCSRD_U(1,I,IU)
             VUNIT(2) = VUNIT(2) + WCSRD_U(2,I,IU)
             VUNIT(3) = VUNIT(3) + WCSRD_U(3,I,IU)

             
             RROT(1) = RC(1,I) - XYZREF(1)
             RROT(2) = RC(2,I) - XYZREF(2)
             RROT(3) = RC(3,I) - XYZREF(3)
             CALL CROSS(RROT,WUNIT,VUNIT_W_term)
             
             VUNIT = VUNIT + VUNIT_W_term
             
             RHS_U(I,IU) = -DOT(ENC(1,I),VUNIT)
              
             ! Add contribution from control surfaces
C$AD II-LOOP
             DO N = 1, NCONTROL
              RHS_U(I,IU) = RHS_U(I,IU) 
     &                      - DOT(ENC_D(1,I,N),VUNIT)*DELCON(N)
            ENDDO
          ELSE
            RHS_U(I,IU) = 0
          endif
          
        enddo
        end !set_vel_rhs_u
      
      SUBROUTINE set_gam_d_rhs(IQ,ENC_Q, rhs_vec)
      INCLUDE 'AVL.INC'
      REAL ENC_Q(3,NVMAX,*), rhs_vec(NVMAX)
      REAL RROT(3), VROT(3), VQ(3)
C
C$AD II-LOOP
      DO I = 1, NVOR
        IF(LVNC(I)) THEN
          IF(LVALBE(I)) THEN
          RROT(1) = RC(1,I) - XYZREF(1)
          RROT(2) = RC(2,I) - XYZREF(2)
          RROT(3) = RC(3,I) - XYZREF(3)
          CALL CROSS(RROT,WROT,VROT)
C$AD II-LOOP
          DO K = 1, 3
            VQ(K) = VINF(K)
     &              + VROT(K)
          ENDDO
          ELSE
          VQ(1) = 0.
          VQ(2) = 0.
          VQ(3) = 0.
          ENDIF
C$AD II-LOOP
          DO K = 1, 3
            VQ(K) = VQ(K)
     &             + WCSRD_U(K,I,1)*VINF(1)
     &             + WCSRD_U(K,I,2)*VINF(2)
     &             + WCSRD_U(K,I,3)*VINF(3)
     &             + WCSRD_U(K,I,4)*WROT(1)
     &             + WCSRD_U(K,I,5)*WROT(2)
     &             + WCSRD_U(K,I,6)*WROT(3)
          ENDDO
C 
          rhs_vec(I) = -DOT(ENC_Q(1,I,IQ),VQ)
        ELSE
          rhs_vec(I) = 0.
        ENDIF
      ENDDO
C
      RETURN
      END ! set_gam_d_rhs

      
      subroutine mat_prod(mat, vec, n, out_vec)
        include 'AVL.INC'
        real mat(NVMAX,NVMAX), vec(NVMAX), out_vec(NVMAX)
        
        out_vec = 0.
C$AD II-LOOP
        DO J = 1, n
C$AD II-LOOP
          DO I = 1, n
            out_vec(I) = out_vec(I) +  mat(I,J)*vec(J)
          enddo 
        enddo
        
      end !mat_prod
        
                 
        
