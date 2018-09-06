!----------------------------------------
!     Number of sites along t-direction
integer nsite
parameter(nsite=6)
!----------------------------------------
!     Size of matrices
integer nmat
parameter(nmat=4)
!----------------------------------------
!     Number of pseudo fermions
integer npf !(1,2 at this moment)
parameter (npf=2)
!----------------------------------------
!     number of remez coefficients
integer nremez_md,nremez_pf
parameter(nremez_md=15)
parameter(nremez_pf=15)
!--------------------------------------
integer ndim ! fix to 9
parameter(ndim=9)
integer nspin ! fix to 16
parameter(nspin=16)
!----------------------------------------
!   nimprove=0 -> no improvement
!   nimprove=1 -> O(a^2) improvement. Note that communication cost increases.
integer nimprove
parameter(nimprove=1)
integer nmargin!size of the margin
parameter(nmargin=nimprove+1)
!----------------------------------------
!  ngauge=1 -> gauged
!  ngauge=0 -> ungauged
integer ngauge
parameter(ngauge=1)
double precision, parameter :: pi1=2d0*dasin(1d0)
!----------------------------------------
! debugging switches
! solver: 0,1,2,3
integer,parameter :: solver_verbose=0
integer,parameter :: rhmc_verbose=0
integer,parameter :: check_host_metropolis=0
integer,parameter :: cublasmult=0
integer,parameter :: havecublas=1
integer,parameter :: usetimer=1

