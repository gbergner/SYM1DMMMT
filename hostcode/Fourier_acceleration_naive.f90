!aho=1 -> xmat to xmat_mom
!aho=2 -> xmat_mom to xmat
subroutine Fourier_acceleration_naive(acceleration)

  implicit none
  include '../staticparameters.f90'

  double precision acceleration(1:nsite)
  integer imom,k,p,Lambda

  do imom=1,nsite
     ! k = 1,2,....,nsite = 1,2,...,Lambda
     k=imom
     Lambda=nsite
     ! When Lambda is even, 
     ! p = -Lambda/2+1,...,Lambda/2. 
     ! k > Lambda/2 corresponds to p=k-Lambda. 
     
     ! When Lambda is odd, 
     ! p = -Lambda/2+1/2,...,Lambda/2-1/2. 
     ! k > Lambda/2-1/2 corresponds to p=k-Lambda.
     
     ! In any case, k > int(Lamdbda/2+0.01) corresponds to p=k-Lambda. 
     
     if(k.GT.int(dble(Lambda)*0.5d0+0.01d0))then
        p=k-Lambda
     else
        p=k
     end if
     if(p.GT.0)then     
        acceleration(imom)=1d0/dble(p)
     else if(p.LT.0)then
        acceleration(imom)=1d0/dble(-p)
     else
        acceleration(imom)=2d0
     end if
  end do
!  if(myrank.EQ.nprocs-1)then
!     acceleration(nsite_local)=2d0
!  end if

  return

END subroutine Fourier_acceleration_naive
