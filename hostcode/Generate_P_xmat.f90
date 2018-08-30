! Generate P_xmat with Gaussian weight.
! We do not have to care about the traceless condition. 
SUBROUTINE Generate_P_xmat(P_xmat)
  
  implicit none
  include '../staticparameters.f90'
  integer imat,jmat,idim,isite
  double precision r1,r2

  double complex P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  
    
  do isite=1,nsite
     do idim=1,ndim
        do imat=1,nmat-1
           do jmat=imat+1,nmat
              call BoxMuller(r1,r2)
              P_xmat(imat,jmat,idim,isite)=&
                   dcmplx(r1/dsqrt(2d0))+dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
              P_xmat(jmat,imat,idim,isite)=&
                   dcmplx(r1/dsqrt(2d0))-dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
           end do
        end do
        do imat=1,nmat
           call BoxMuller(r1,r2)
           P_xmat(imat,imat,idim,isite)=dcmplx(r1)
        end do
     end do
  end do
  
  return
  
END SUBROUTINE Generate_P_xmat
