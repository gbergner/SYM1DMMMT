!----------------------------------------
!     Number of blocks along t-direction and sites in each block
integer nsite_local,nsublat
parameter(nsublat=2)
parameter(nsite_local=4)
!----------------------------------------
!     Size of matrices
!     nmat=nmat_block*nblock
integer nmat_block,nblock
parameter(nmat_block=4)
parameter(nblock=2)
!----------------------------------------
!     number of remez coefficients
integer nremez_md,nremez_pf
parameter(nremez_md=15)
parameter(nremez_pf=15)
!----------------------------------------
!     Number of scalars (fix to 9)
integer ndim
parameter(ndim=9)
!----------------------------------------
!     Number of spinor indices (fix to 16)
integer nspin
parameter(nspin=16)
!----------------------------------------
!   nimprove=0 -> no improvement
!   nimprove=1 -> O(a^2) improvement. Note that communication cost increases.
integer nimprove
parameter(nimprove=1)
integer nmargin!size of the margin
parameter(nmargin=nimprove+1)
