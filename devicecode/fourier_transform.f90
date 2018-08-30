!nxp=1 -> P_xmat to P_xmat_mom
!nxp=2 -> P_xmat_mom to P_xmat
module fourier_transform
  implicit none
  contains
subroutine FT_P_xmat_device(P_xmat,P_xmat_mom)

  use compiletimeconstants
  implicit none
  

  double complex,intent(in) :: P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double complex, intent(out) :: P_xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
   !$acc declare device_resident(P_xmat,P_xmat_mom)
  
  double complex phase
  double precision temp,pi

  integer isite,idim,imom
  integer imat,jmat

  pi=2d0*dasin(1d0)
     !$acc kernels
     P_xmat_mom=(0d0,0d0)
     do imom=1,nsite
        do isite=1,nsite
           temp=2d0*pi*dble(imom*isite)/dble(nsite)
           phase=dcmplx(dcos(temp))-(0d0,1d0)*dcmplx(dsin(temp))
           phase=phase/dcmplx(dsqrt(dble(nsite)))
           do idim=1,ndim
              do imat=1,nmat
                 do jmat=1,nmat
                    P_xmat_mom(imat,jmat,idim,imom)=&
                         &P_xmat_mom(imat,jmat,idim,imom)&
                         &+P_xmat(imat,jmat,idim,isite)*phase
                 end do
              end do
           end do
        end do
     end do
     !$acc end kernels
  return

END subroutine FT_P_xmat_device

subroutine FTinv_P_xmat_device(P_xmat,P_xmat_mom)

  implicit none
  
  include '../staticparameters.f90'

  double complex,intent(out) :: P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double complex, intent(in) :: P_xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
   !$acc declare device_resident(P_xmat,P_xmat_mom)
  
  double complex phase
  double precision temp,pi

  integer isite,idim,imom
  integer imat,jmat

  pi=2d0*dasin(1d0)

      !$acc kernels
     P_xmat=(0d0,0d0)
     do imom=1,nsite
        do isite=1,nsite
           temp=2d0*pi*dble(isite*imom)/dble(nsite)
           phase=dcmplx(dcos(temp))+(0d0,1d0)*dcmplx(dsin(temp))
           phase=phase/dcmplx(dsqrt(dble(nsite)))
           do idim=1,ndim
              do imat=1,nmat
                 do jmat=1,nmat
                    P_xmat(imat,jmat,idim,isite)=&
                         &P_xmat(imat,jmat,idim,isite)&
                         &+P_xmat_mom(imat,jmat,idim,imom)*phase
                 end do
              end do
           end do
        end do
     end do
     !$acc end kernels
  return

END subroutine FTinv_P_xmat_device


subroutine FT_xmat_device(xmat,xmat_mom)

  implicit none

  include '../staticparameters.f90'

  double complex, intent(in):: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double complex,intent(out):: xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
  !$acc declare device_resident(xmat_mom)
  !$acc declare present(xmat)
  
  double complex phase
  double precision temp,pi

  integer isite,idim,imom
  integer imat,jmat

  pi=2d0*dasin(1d0)
   !$acc kernels
   xmat_mom=(0d0,0d0)
   !$acc end kernels
   !$acc kernels
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
  !$acc end kernels

END subroutine FT_xmat_device

subroutine FTinv_xmat_device(xmat,xmat_mom)

  implicit none

  include '../staticparameters.f90'

  integer nxp
  double complex, intent(out):: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double complex,intent(in):: xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
  !$acc declare device_resident(xmat_mom)
  !$acc declare present(xmat)
  
  double complex phase
  double precision temp,pi

  integer isite,idim,imom
  integer imat,jmat

  pi=2d0*dasin(1d0)
  !$acc kernels
  xmat=(0d0,0d0)
  !$acc end kernels
  !$acc kernels
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
 !$acc end kernels
     
  return

END subroutine FTinv_xmat_device

end module fourier_transform
