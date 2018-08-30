SUBROUTINE hermitian_projection_device(xmat)
  
  use compiletimeconstants
  implicit none
 
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  !$acc declare present(xmat)
  integer idim,imat,jmat,isite
  
  !$acc kernels
  do idim=1,ndim
     do isite=-(nmargin-1),nsite+nmargin
        do imat=1,nmat-1
           do jmat=imat+1,nmat
              xmat(jmat,imat,idim,isite)=dconjg(xmat(imat,jmat,idim,isite))
           end do
        end do
        do imat=1,nmat
           xmat(imat,imat,idim,isite)=dcmplx(dble(xmat(imat,imat,idim,isite)))
        end do
     end do
  end do
  !$acc end kernels
  return
  
END SUBROUTINE hermitian_projection_device
