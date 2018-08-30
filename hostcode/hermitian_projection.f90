SUBROUTINE hermitian_projection(xmat)
  
  implicit none
  include '../staticparameters.f90'
 
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  integer idim,imat,jmat,isite
  
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
  
  return
  
END SUBROUTINE 
