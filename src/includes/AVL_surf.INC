C      !start added variables for python geometry minipulation

      LOGICAL LDUPL_B
      COMMON /BODY_GEOM_L/
     & LDUPL_B(NBMAX)      ! T if body is a duplicated
c
      COMMON /BODY_GEOM_I/
     & NSEC_B(NBMAX),  ! number of sections in body
     & NVB(NFMAX),     ! number of body elements
     & NBOD(NBMAX)    ! number of body oml points

      REAL(kind=avl_real) XYZSCAL_B
      REAL(kind=avl_real) XYZTRAN_B
      REAL(kind=avl_real) YDUPL_B
      REAL(kind=avl_real) XYZLES_B 
      REAL(kind=avl_real) BSPACE
      REAL(kind=avl_real) XBOD
      REAL(kind=avl_real) YBOD
      REAL(kind=avl_real) TBOD
      REAL(kind=avl_real) XBOD_R
      REAL(kind=avl_real) YBOD_R
      COMMON /BODY_GEOM_R/
     & XYZSCAL_B(3, NBMAX), ! scaling factors for XYZ coordinates
     & XYZTRAN_B(3, NBMAX), ! translation factors for XYZ coordinates
     & YDUPL_B(NBMAX), ! y duplicate
     & XYZLES_B(3, NSMAX, NBMAX), ! leading edge cordinate vector
     & BSPACE(NBMAX), ! body spacing
     & XBOD(IBX, NBMAX), ! body oml x-coordinates
     & YBOD(IBX, NBMAX), ! body oml y-coordinates
     & TBOD(IBX, NBMAX), ! body oml thickness
     & XBOD_R(IBX, NBMAX), ! raw input oml x-coordinates
     & YBOD_R(IBX, NBMAX) ! raw input oml y-coordinates
c
      LOGICAL LDUPL
      LOGICAL LSURFSPACING
      COMMON /SURF_GEOM_L/
     & LDUPL(NFMAX),      ! T if surface is a duplicated
     & LSURFSPACING(NFMAX)   ! surface spacing set under the surface heeading

      COMMON /SURF_GEOM_I/
     & NSEC(NFMAX),  ! number of sections in surface
     & NVC(NFMAX), ! number of chordwise elements
     & NVS(NFMAX), ! number of spanwise elements
     & NSPANS(NSECMAX, NFMAX), ! number of spanwise elements vector
     & NASEC(NSECMAX, NFMAX), ! number of points used to represent section geometry
     & NAPTSSEC(NSECMAX, IBX), ! number of points of input airfoil geometry
     & ICONTD(ICONX, NSECMAX, NFMAX), ! control variable index
     & NSCON(NSECMAX, NFMAX), ! number of control variables
     & IDESTD(ICONX, NSECMAX, NFMAX), ! design variable index
     & NSDES(NSECMAX, NFMAX) ! number of design variables

      REAL(kind=avl_real) XYZSCAL
      REAL(kind=avl_real) XYZTRAN
      REAL(kind=avl_real) YDUPL
      REAL(kind=avl_real) ADDINC
      REAL(kind=avl_real) CSPACE
      REAL(kind=avl_real) SSPACE
      REAL(kind=avl_real) SSPACES
      REAL(kind=avl_real) XYZLES
      REAL(kind=avl_real) XSEC
      REAL(kind=avl_real) YSEC
      REAL(kind=avl_real) XFMIN_R
      REAL(kind=avl_real) XFMAX_R
      REAL(kind=avl_real) XLASEC
      REAL(kind=avl_real) ZLASEC
      REAL(kind=avl_real) XUASEC
      REAL(kind=avl_real) ZUASEC
      REAL(kind=avl_real) CHORDS
      REAL(kind=avl_real) AINCS
      REAL(kind=avl_real) XASEC
      REAL(kind=avl_real) CASEC
      REAL(kind=avl_real) SASEC
      REAL(kind=avl_real) TASEC
      REAL(kind=avl_real) CLCDSEC
      REAL(kind=avl_real) CLCDSRF
      REAL(kind=avl_real) CLAF
      REAL(kind=avl_real) XHINGED
      REAL(kind=avl_real) VHINGED
      REAL(kind=avl_real) GAIND
      REAL(kind=avl_real) REFLD
      REAL(kind=avl_real) GAING  
      COMMON /SURF_GEOM_R/
     & XYZSCAL(3, NFMAX), ! scaling factors for XYZ coordinates
     & XYZTRAN(3, NFMAX), ! translation factors for XYZ coordinates
     & YDUPL(NFMAX), ! y position of duplicate plane
     & ADDINC(NFMAX), ! additional incidence angle
     & CSPACE(NFMAX), ! chordwise spacing
     & SSPACE(NFMAX), ! spanwise spacing
     & SSPACES(NSECMAX, NFMAX), ! spanwise spacing vector
     & XYZLES(3, NSECMAX, NFMAX), ! leading edge cordinate vector
     & XSEC(IBX, NSECMAX, NFMAX), ! raw section input x coordinates
     & YSEC(IBX, NSECMAX, NFMAX), ! raw section input y coordinates
     & XFMIN_R(NSECMAX, NFMAX), ! section input lower x bound 
     & XFMAX_R(NSECMAX, NFMAX), ! section input upper x bound 
     & XLASEC(IBX, NSECMAX, NFMAX), ! airfoil lower surface x-coords
     & ZLASEC(IBX, NSECMAX, NFMAX), ! airfoil lower surface z-coords
     & XUASEC(IBX, NSECMAX, NFMAX), ! airfoil upper surface x-coords
     & ZUASEC(IBX, NSECMAX, NFMAX), ! airfoil upper surface z-coords
     & CHORDS(NSECMAX, NFMAX), ! chord length vectorm
     & AINCS(NSECMAX, NFMAX), ! incidence angle vector
     & XASEC(IBX, NSECMAX, NFMAX),  ! the x coordinate aifoil section
     & CASEC(IBX, NSECMAX, NFMAX),  ! camber line at xasec
     & SASEC(IBX, NSECMAX, NFMAX),  ! slope of camber line at xasec
     & TASEC(IBX, NSECMAX, NFMAX),   ! thickness at xasec
     & CLCDSEC(6, NSECMAX, NFMAX),
     & CLCDSRF(6, NFMAX),
     & CLAF(NSECMAX, NFMAX),   ! CL alpha (dCL/da) scaling factor
     & XHINGED(ICONX, NSECMAX, NFMAX), ! hinge location
     & VHINGED(3, ICONX, NSECMAX, NFMAX), ! hinge vector
     & GAIND(ICONX, NSECMAX, NFMAX), ! control surface gain
     & REFLD(ICONX, NSECMAX, NFMAX), ! control surface reflection
     & GAING(ICONX, NSECMAX, NFMAX)  ! desgin variable gain

      COMMON /SURF_MESH_I/
     & MFRST(NFMAX), ! stores the index in the MSHBLK where each surface's mesh begins
     & IPTSEC(NSMAX,NFMAX) ! stores the iptloc vector for each surface

      REAL(kind=avl_real) MSHBLK
      REAL(kind=avl_real) RV1MSH
      REAL(kind=avl_real) RV2MSH
      REAL(kind=avl_real) RVMSH
      REAL(kind=avl_real) RCMSH
      COMMON /SURF_MESH_R/
     & MSHBLK(3, 4*NVMAX), ! block to store all surface meshes
     & RV1MSH(3,NVMAX),    !  mesh h.v. vortex left  points
     & RV2MSH(3,NVMAX),    !  mesh h.v. vortex right  points
     & RVMSH(3,NVMAX),    ! mesh h.v. vortex center points
     & RCMSH(3,NVMAX)    ! mesh h.v. control points

      LOGICAL LSURFMSH
      LOGICAL LMESHFLAT
      COMMON /SURF_MESH_L/
     & LSURFMSH(NFMAX), ! T if surface uses a mesh cordinates for its geometry
     & LMESHFLAT(NFMAX) ! T if the surface's mesh should be flattened when computing vorticies and control points


C      !!--- end added variables for python geometry manipulation ---
