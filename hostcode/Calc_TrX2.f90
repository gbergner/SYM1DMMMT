!************************************************
!*** Calculate trx2 = (1/N)*¥int dt Tr(X_I^2) ***
!************************************************
subroutine Calc_TrX2(xmat,trx2)
  
  implicit none
  
  include '../staticparameters.f90'

  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double precision trx2
  
  integer isite,idim
  integer imat,jmat

  trx2=0d0
  do isite=1,nsite
     do idim=1,ndim
        do jmat=1,nmat
           do imat=1,nmat
              trx2=trx2&
                   +dble(xmat(imat,jmat,idim,isite)&
                   *dconjg(xmat(imat,jmat,idim,isite)))
           end do
        end do
     end do
  end do
  
  trx2=trx2/dble(nmat*nsite)
  
  return
  
END subroutine Calc_TrX2
