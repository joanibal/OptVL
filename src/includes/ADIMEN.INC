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

      PARAMETER (NVMAX=5000)

      PARAMETER (NSMAX=500)
      PARAMETER (NFMAX=100)
      PARAMETER (NLMAX=500)
      PARAMETER (NBMAX=20)

      PARAMETER (NUMAX=6)
      PARAMETER (NDMAX=30)
      PARAMETER (NGMAX=20)

      PARAMETER (NRMAX=25)
      PARAMETER (NTMAX=500)

      PARAMETER (NOBMAX=1000)
      
      PARAMETER (ICONX=20)
      PARAMETER (IBX=300)

