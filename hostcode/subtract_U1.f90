!trace part of (¥int dt X) and alpha are removed. 
SUBROUTINE subtract_U1(xmat,alpha)
  
  implicit none
  include '../staticparameters.f90'
  !***** input & output *****
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double precision alpha(1:nmat)
  !****************************
  integer idim,imat,isite
  double complex trace(1:ndim)
  double precision sum_alpha
  
  !***********************************************
  !**** trace part of (¥int dt X) is removed. ****
  !***********************************************
  trace=(0d0,0d0)!sum of contributions from all processes.
  do idim=1,ndim
     do isite=1,nsite
        do imat=1,nmat
           trace(idim)=trace(idim)+xmat(imat,imat,idim,isite)
        end do
     end do
  end do
  trace=trace/dcmplx(nsite*nmat)
  do idim=1,ndim
     !take care of the margin too. 
     do isite=-(nmargin-1),nsite+nmargin
        do imat=1,nmat
           xmat(imat,imat,idim,isite)=xmat(imat,imat,idim,isite)-trace(idim)
        end do
     end do
  end do

  !*********************************************************
  !*** Trace part of alpha is removed. *********************
  !*** Do the same at all nodes, to avoid communication. ***
  !*********************************************************  
  sum_alpha=0d0
  do imat=1,nmat
     sum_alpha=sum_alpha+alpha(imat)
  end do
  sum_alpha=sum_alpha/dble(nmat)
  do imat=1,nmat
     alpha(imat)=alpha(imat)-sum_alpha
  end do
  
  return
  
END SUBROUTINE subtract_U1
