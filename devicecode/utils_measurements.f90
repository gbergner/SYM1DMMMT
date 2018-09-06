!************************************************
!*** Calculate trx2 = (1/N)*¥int dt Tr(X_I^2) ***
!************************************************
module utils_measurements
    implicit none
contains
    subroutine Calc_TrX2_device(xmat,trx2)
  
        use compiletimeconstants
        implicit none

        double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        !$acc declare present(xmat)
        double precision trx2
  
        integer isite,idim
        integer imat,jmat

        trx2=0d0
        !$acc kernels
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
        !$acc end kernels
        trx2=trx2/dble(nmat*nsite)
  
        return
  
    END subroutine Calc_TrX2_device

    !trace part of (¥int dt X) and alpha are removed.
    SUBROUTINE subtract_U1_device(xmat,alpha)

        use compiletimeconstants
        implicit none

        !***** input & output *****
        double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double precision alpha(1:nmat)
        !$acc declare present(xmat,alpha)
        !****************************
        integer idim,imat,isite
        double complex trace(1:ndim),tmp
        double precision sum_alpha

        !***********************************************
        !**** trace part of (¥int dt X) is removed. ****
        !***********************************************
        do idim=1,ndim
            tmp=(0d0,0d0)
            !$acc kernels
            do isite=1,nsite
                do imat=1,nmat
                    tmp=tmp+xmat(imat,imat,idim,isite)
                end do
            end do
            !$acc end kernels
            tmp=tmp/dcmplx(nsite*nmat)
            !take care of the margin too.
            !$acc kernels
            do isite=-(nmargin-1),nsite+nmargin
                do imat=1,nmat
                    xmat(imat,imat,idim,isite)=xmat(imat,imat,idim,isite)-tmp
                end do
            end do
            !$acc end kernels
            trace(idim)=tmp/dcmplx(nsite*nmat)
        end do
        !*********************************************************
        !*** Trace part of alpha is removed. *********************
        !*** Do the same at all nodes, to avoid communication. ***
        !*********************************************************
        !$acc kernels
        sum_alpha=Sum(alpha)
        !$acc end kernels
        sum_alpha=sum_alpha/dble(nmat)
        !$acc kernels
        do imat=1,nmat
            alpha(imat)=alpha(imat)-sum_alpha
        end do
        !$acc end kernels
        return

    END SUBROUTINE subtract_U1_device
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

end module utils_measurements
