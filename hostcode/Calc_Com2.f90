!*******************************************************
!*** Calculate com2 = (-1/N)*¥int dt Tr([X_I,X_J]^2) ***
!*******************************************************
subroutine Calc_Com2(xmat,com2)

  implicit none

  include '../staticparameters.f90'
  !********** input ********** 
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  !********** outout **********
  double precision com2
  !****************************
  double complex com(1:nmat,1:nmat)
  integer isite,idim,jdim
  integer imat,jmat,kmat

  com2=0d0
  do isite=1,nsite
     do idim=1,ndim-1
        do jdim=idim+1,ndim
           com=(0d0,0d0)
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    com(imat,jmat)=com(imat,jmat)+&
                         &xmat(imat,kmat,idim,isite)*xmat(kmat,jmat,jdim,isite)&
                         &-xmat(imat,kmat,jdim,isite)*xmat(kmat,jmat,idim,isite)
                 end do
              end do
           end do
           do jmat=1,nmat
              do imat=1,nmat
                 com2=com2+&
                      &dble(com(imat,jmat)*dconjg(com(imat,jmat)))
              end do
           end do
        end do
     end do
  end do

  com2=com2/dble(nmat*nsite)*2d0

  return

END subroutine Calc_Com2
