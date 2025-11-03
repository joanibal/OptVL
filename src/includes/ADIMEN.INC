!   Parameters for Athena 3-D vortex lattice code
!
!   Primary array limits
!    NVMAX   number of horseshoe vortices
!
!    NSMAX   number of chord strips
!    NFMAX   number of surfaces
!
!    NLMAX   number of source/doublet line nodes
!    NBMAX   number of bodies
!
!    NUMAX   number of freestream parameters (V,Omega)
!    NDMAX   number of control deflection parameters
!    NGMAX   number of design variables
!
!    NRMAX   number of stored run cases
!    NTMAX   number of stored time levels
!
!    NOBMAX  number of off-body points
!
!    ICONX   max number of control or design variable declaration lines per section (20)
!    IBX     max number of airfoil section points
      INTEGER NVMAX
      INTEGER NSMAX,NFMAX,NLMAX,NBMAX
      INTEGER NUMAX,NDMAX,NGMAX
      INTEGER NRMAX,NTMAX
      INTEGER NOBMAX
      INTEGER ICONX
      INTEGER IBX

!     WARNING: This values have been tweaked compared to AVL 
!              to make it easier to termine the true size in clear_ad_seeds_fast()
!              I increased some values to avoid matching

      PARAMETER (NVMAX=2500)

      PARAMETER (NSMAX=500) ! max number of strips
      PARAMETER (NSECMAX=501) ! max number of geometry sections
      PARAMETER (NFMAX=100) ! max number of surfaces
      PARAMETER (NLMAX=502) ! max number of body nodes
      PARAMETER (NBMAX=20)  ! max number of bodies

      PARAMETER (NUMAX=6)  ! number of velocity components
      PARAMETER (NDMAX=30) ! max number of control surfaces
      PARAMETER (NGMAX=21) ! max number of design variables

      PARAMETER (NRMAX=25)  ! max number of run cases
      PARAMETER (NTMAX=503) ! max number of time values

      PARAMETER (NOBMAX=100) ! max number of off body points
      
      PARAMETER (ICONX=20) ! 
      PARAMETER (IBX=300)  ! max number of airfoil coordinates

