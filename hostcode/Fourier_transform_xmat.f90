!nxp=1 -> xmat to xmat_mom
!nxp=2 -> xmat_mom to xmat
subroutine Fourier_transform_xmat(xmat,xmat_mom,nxp)

  implicit none

  include '../staticparameters.f90'

  integer nxp
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double complex xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)

  
  double complex phase
  double precision temp,pi

  integer isite,idim,imom
  integer imat,jmat

  pi=2d0*dasin(1d0)

  if (nxp.EQ.1)then
     

     xmat_mom=(0d0,0d0)
     do imom=1,nsite
        do isite=1,nsite
           temp=2d0*pi*dble(imom*isite)/dble(nsite)
           phase=dcmplx(dcos(temp))-(0d0,1d0)*dcmplx(dsin(temp))
           phase=phase/dcmplx(dsqrt(dble(nsite)))
           do idim=1,ndim
              do imat=1,nmat
                 do jmat=1,nmat
                    xmat_mom(imat,jmat,idim,imom)=&
                         &xmat_mom(imat,jmat,idim,imom)&
                         &+xmat(imat,jmat,idim,isite)*phase
                 end do
              end do
           end do
        end do
     end do

  else if(nxp.EQ.2)then


     xmat=(0d0,0d0)
     do isite=1,nsite
        do imom=1,nsite
           temp=2d0*pi*dble(isite*imom)/dble(nsite)
           phase=dcmplx(dcos(temp))+(0d0,1d0)*dcmplx(dsin(temp))
           phase=phase/dcmplx(dsqrt(dble(nsite)))
           do idim=1,ndim
              do imat=1,nmat
                 do jmat=1,nmat
                    xmat(imat,jmat,idim,isite)=&
                         &xmat(imat,jmat,idim,isite)&
                         &+xmat_mom(imat,jmat,idim,imom)*phase
                 end do
              end do
           end do
        end do
     end do
 
     
  end if

  return

END subroutine Fourier_transform_xmat
