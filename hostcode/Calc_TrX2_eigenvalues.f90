!************************************************
!*** Calculate trx2 = (1/N)*¥int dt Tr(X_I^2) ***
!************************************************
subroutine DiagMat(x2mat,ev,msize)
    implicit none
    double complex x2mat(1:msize,1:msize)
    double complex ev(1:msize)
    character jobvl,jobvr
    integer lda,ldvl,ldvr,info,lwork,msize
    double complex VL(1:msize),VR(1:msize,1:msize)
    ! work size = lwork determined without workspace query.
    double complex  work(1:2*msize)
    double precision rwork(1:2*msize)
    ! no left and right eigenvectors needed. (VL,VR not needed,but I still leave it.)
    jobvl='N'
    ldvl=1
    jobvr='N'
    ldvr=1
!    ldvr=msize
    !leading dimension
    lda=msize
    lwork=2*msize
    call ZGEEV(jobvl,jobvr,msize,x2mat,lda,ev,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
    if(info.ne.0) then
      print *, "ERROR IN EIGENVALUE DIAGONALIZATION"
    end if
end subroutine DiagMat


subroutine Calc_TrX2_eigenvalues(xmat,eigxmat)
  
    implicit none
  
  include '../staticparameters.f90'

    double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    double complex x2mat(1:ndim,1:ndim)
    double complex eigxmat(1:ndim)
  
    integer isite,idim1,idim2
    integer imat,jmat

    x2mat=0d0
 !   print *,'\n t=['
    do idim1=1,ndim
        do idim2=1,ndim
            do isite=1,nsite
                do jmat=1,nmat
                    do imat=1,nmat
                        x2mat(idim1,idim2)=x2mat(idim1,idim2)&
                            +(xmat(imat,jmat,idim1,isite)&
                            *dconjg(xmat(imat,jmat,idim2,isite)))
                    end do
                end do
            end do
            x2mat(idim1,idim2)=x2mat(idim1,idim2)/dble(nsite)
 !           print *,dble(x2mat(idim1,idim2)),','
        end do
 !       print *,';'
    end do
!    print *,']'
    call DiagMat(x2mat,eigxmat,ndim)
  
    return
  
END subroutine Calc_TrX2_eigenvalues
