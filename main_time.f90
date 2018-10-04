!##############################################################################
!######              BFSS matrix model on lattice                     #########
!######              Benchmark for CPU and GPU multiplications        #########
!######              written by Georg Bergner                         #########
!##############################################################################

SUBROUTINE Test_Multiply_Dirac_dagger_device(xmat,pf1,pf2,var)
    !********************************(needed for cublas version)
    use cublasinterface
    use compiletimeconstants
    use gammamatrix
    use cublas
    implicit none

    integer, intent(in) :: var
    double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    double complex, intent(in) :: pf1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex, intent(out) :: pf2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex :: gamtmp,tmp

    double complex :: xmat2(1:nmat,1:nmat,-(nmargin-1):nsite+nmargin,1:ndim)
    double complex :: pf12(1:nmat,1:nmat,-(nmargin-1):nsite+nmargin,1:nspin)
    double complex :: pf22(1:nmat,1:nmat,-(nmargin-1):nsite+nmargin,1:nspin)
    !******************
    integer :: imat,jmat,kmat
    integer :: idim
    integer :: ispin,jspin,kspin
    integer :: isite
    integer :: countr
    !$acc declare present(xmat,pf1,pf2)

    !$acc kernels
    pf2=(0d0,0d0)
    !$acc end kernels

    if(var==1) then
        do countr=1,144
            ispin=gamispin(countr)
            jspin=gamjspin(countr)
            idim=gamidim(countr)
            gamtmp=dconjg(gamgam(countr))
            !$acc kernels
            do isite=1,nsite
                do imat=1,nmat
                    do jmat=1,nmat
                        do kmat=1,nmat
                            pf2(imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)-gamtmp&
                                *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                                -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
                        end do
                    end do
                end do
            end do
           !$acc end kernels
        end do
    elseif(var==2) then
        do countr=1,144
            ispin=gamispin(countr)
            jspin=gamjspin(countr)
            idim=gamidim(countr)
            gamtmp=dconjg(gamgam(countr))
            !$acc parallel &
            !$acc loop independent &
            !$acc collapse(3)
            do isite=1,nsite
                do imat=1,nmat
                    do jmat=1,nmat
                        !$acc loop seq
                        do kmat=1,nmat
                            pf2(imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)-gamtmp&
                                *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                                -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
                        end do
                    end do
                end do
            end do
           !$acc end parallel loop
        end do
    elseif(var==3) then
        do countr=1,144
            ispin=gamispin(countr)
            jspin=gamjspin(countr)
            idim=gamidim(countr)
            gamtmp=dconjg(gamgam(countr))
            !$acc parallel &
            !$acc loop independent &
            !$acc collapse(3)&
            !$acc private(tmp)
            do isite=1,nsite
                do imat=1,nmat
                    do jmat=1,nmat
                        tmp=(0d0,0d0)
                        !$acc loop reduction(+:tmp)
                        do kmat=1,nmat
                            tmp=tmp-gamtmp&
                                *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                                -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
                        end do
                        pf2(imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)+tmp
                    end do
                end do
            end do
           !$acc end parallel loop
        end do
    elseif(var==4) then
        do countr=1,144
            ispin=gamispin(countr)
            jspin=gamjspin(countr)
            idim=gamidim(countr)
            gamtmp=dconjg(gamgam(countr))
            !$acc parallel &
            !$acc loop independent&
            !$acc async(ispin+5)&
            !$acc collapse(3)
            do isite=1,nsite
                do imat=1,nmat
                    do jmat=1,nmat
                        !$acc loop seq
                        do kmat=1,nmat
                            pf2(imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)-gamtmp&
                                *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                                -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
                        end do
                    end do
                end do
            end do
           !$acc end parallel loop
        end do
        do ispin=1,nspin
          !$acc wait(ispin+5)
        end do
    elseif(var==5) then
        do countr=1,144
            ispin=gamispin(countr)
            jspin=gamjspin(countr)
            idim=gamidim(countr)
            gamtmp=dconjg(gamgam(countr))
            !$acc kernels
            do isite=1,nsite
                do kmat=1,nmat
                    do imat=1,nmat
                        do jmat=1,nmat
                            pf2(imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)-gamtmp&
                                *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                                -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
                        end do
                    end do
                end do
            end do
           !$acc end kernels
        end do
    elseif(var==6) then
        do countr=1,144
            ispin=gamispin(countr)
            jspin=gamjspin(countr)
            idim=gamidim(countr)
            gamtmp=dconjg(gamgam(countr))
            !$acc kernels&
            !$acc async(ispin+5)
            do isite=1,nsite
                do kmat=1,nmat
                    do imat=1,nmat
                        do jmat=1,nmat
                            pf2(imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)-gamtmp&
                                *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                                -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
                        end do
                    end do
                end do
            end do
           !$acc end kernels
        end do
        do ispin=1,nspin
          !$acc wait(ispin+5)
        end do
    elseif(var==7) then
        do countr=1,144
            ispin=gamispin(countr)
            jspin=gamjspin(countr)
            idim=gamidim(countr)
            gamtmp=dconjg(gamgam(countr))
            !$acc kernels
            do isite=1,nsite
                do kmat=1,nmat
                    do imat=1,nmat
                        do jmat=1,nmat
                            pf2(imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)-gamtmp*xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)
                        end do
                    end do
                end do
                do kmat=1,nmat
                    do imat=1,nmat
                        do jmat=1,nmat
                            pf2(imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)+gamtmp&
                                *xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite)
                        end do
                    end do
                end do
            end do
            !$acc end kernels
        end do
    elseif(var==8) then
        do countr=1,144
            ispin=gamispin(countr)
            jspin=gamjspin(countr)
            idim=gamidim(countr)
            gamtmp=dconjg(gamgam(countr))
                call callZgemmBatched(pf2,xmat,pf1,gamtmp,dcmplx(1.d0),ispin,jspin,idim)
        end do
    elseif(var==9) then
        do countr=1,144
            ispin=gamispin(countr)
            jspin=gamjspin(countr)
            idim=gamidim(countr)
            gamtmp=dconjg(gamgam(countr))
                do isite=1,nsite
                    call  zgemm('n', 'n',nmat,nmat,nmat,-gamtmp,xmat(:,:,idim,isite),nmat,pf1(:,:,jspin,isite),nmat,dcmplx(1.d0),pf2(:,:,ispin,isite),nmat)
                    call  zgemm('n', 'n',nmat,nmat,nmat,gamtmp,pf1(:,:,jspin,isite),nmat,xmat(:,:,idim,isite),nmat,dcmplx(1.d0),pf2(:,:,ispin,isite),nmat)
                end do
        end do
    end if
END SUBROUTINE Test_Multiply_Dirac_dagger_device

SUBROUTINE reorder_fields_xmat1(xmatin,xmatout)
    use compiletimeconstants
    use gammamatrix
    implicit none
    double complex,intent(in) :: xmatin(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    double complex,intent(out) :: xmatout(1:nmat,1:nmat,-(nmargin-1):nsite+nmargin,1:ndim)
    !$acc declare present(xmatin,xmatout)
    integer :: imat,jmat,idim, isite
    !$acc kernels
    do isite=-(nmargin-1),nsite+nmargin
        do idim=1,ndim
            do imat=1,nmat
                do jmat=1,nmat
                    xmatout(imat,jmat,isite,idim)=xmatin(imat,jmat,idim,isite)
                end do
            end do
        end do
    end do
!$acc end kernels
END SUBROUTINE reorder_fields_xmat1

SUBROUTINE reorder_fields_pf1(pfin,pfout)
    use compiletimeconstants
    use gammamatrix
    implicit none
    double complex,intent(in):: pfin(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex,intent(out) :: pfout(1:nmat,1:nmat,-(nmargin-1):nsite+nmargin,1:nspin)
        !$acc declare present(pfin,pfout)
    integer :: imat,jmat,idim, isite
    !$acc kernels
    do isite=-(nmargin-1),nsite+nmargin
        do idim=1,nspin
            do imat=1,nmat
                do jmat=1,nmat
                    pfout(imat,jmat,isite,idim)=pfin(imat,jmat,idim,isite)
                end do
            end do
        end do
    end do
!$acc end kernels
END SUBROUTINE reorder_fields_pf1

SUBROUTINE inv_reorder_fields_pf1(pfin,pfout)
    use compiletimeconstants
    use gammamatrix
    implicit none
    double complex,intent(out):: pfout(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex,intent(in) :: pfin(1:nmat,1:nmat,-(nmargin-1):nsite+nmargin,1:nspin)
        !$acc declare present(pfin,pfout)
    integer :: imat,jmat,idim, isite
    !$acc kernels
    do isite=-(nmargin-1),nsite+nmargin
        do idim=1,nspin
            do imat=1,nmat
                do jmat=1,nmat
                    pfout(imat,jmat,idim,isite)=pfin(imat,jmat,isite,idim)
                end do
            end do
        end do
    end do
!$acc end kernels
END SUBROUTINE inv_reorder_fields_pf1

SUBROUTINE reorder_fields_xmat2(xmatin,xmatout)
    use compiletimeconstants
    use gammamatrix
    implicit none
    double complex,intent(in) :: xmatin(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    double complex,intent(out) :: xmatout(-(nmargin-1):nsite+nmargin,1:nmat,1:nmat,1:ndim)
    !$acc declare present(xmatin,xmatout)
    integer :: imat,jmat,idim, isite
    !$acc kernels
    do isite=-(nmargin-1),nsite+nmargin
        do idim=1,ndim
            do imat=1,nmat
                do jmat=1,nmat
                    xmatout(isite,imat,jmat,idim)=xmatin(imat,jmat,idim,isite)
                end do
            end do
        end do
    end do
!$acc end kernels
END SUBROUTINE reorder_fields_xmat2

SUBROUTINE reorder_fields_pf2(pfin,pfout)
    use compiletimeconstants
    use gammamatrix
    implicit none
    double complex,intent(in):: pfin(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex,intent(out) :: pfout(-(nmargin-1):nsite+nmargin,1:nmat,1:nmat,1:nspin)
        !$acc declare present(pfin,pfout)
    integer :: imat,jmat,idim, isite
    !$acc kernels
    do isite=-(nmargin-1),nsite+nmargin
        do idim=1,nspin
            do imat=1,nmat
                do jmat=1,nmat
                    pfout(isite,imat,jmat,idim)=pfin(imat,jmat,idim,isite)
                end do
            end do
        end do
    end do
!$acc end kernels
END SUBROUTINE reorder_fields_pf2

SUBROUTINE inv_reorder_fields_pf2(pfin,pfout)
    use compiletimeconstants
    use gammamatrix
    implicit none
    double complex,intent(out):: pfout(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex,intent(in) :: pfin(-(nmargin-1):nsite+nmargin,1:nmat,1:nmat,1:nspin)
        !$acc declare present(pfin,pfout)
    integer :: imat,jmat,idim, isite
    !$acc kernels
    do isite=-(nmargin-1),nsite+nmargin
        do idim=1,nspin
            do imat=1,nmat
                do jmat=1,nmat
                    pfout(imat,jmat,idim,isite)=pfin(isite,imat,jmat,idim)
                end do
            end do
        end do
    end do
!$acc end kernels
END SUBROUTINE inv_reorder_fields_pf2

SUBROUTINE Test_Multiply_Dirac_dagger_device_reorg1(xmat,pf1,pf2)
    use compiletimeconstants
    use gammamatrix
    implicit none

    double complex :: gamtmp,tmp

    double complex,intent(in) :: xmat(1:nmat,1:nmat,-(nmargin-1):nsite+nmargin,1:ndim)
    double complex,intent(in) :: pf1(1:nmat,1:nmat,-(nmargin-1):nsite+nmargin,1:nspin)
    double complex,intent(out) :: pf2(1:nmat,1:nmat,-(nmargin-1):nsite+nmargin,1:nspin)
    !******************
    integer :: imat,jmat,kmat
    integer :: idim
    integer :: ispin,jspin,kspin
    integer :: isite
    integer :: countr
    !$acc declare present(xmat,pf1,pf2)

    !$acc kernels
    pf2=(0d0,0d0)
    !$acc end kernels

    do countr=1,144
        ispin=gamispin(countr)
        jspin=gamjspin(countr)
        idim=gamidim(countr)
        gamtmp=dconjg(gamgam(countr))
        !$acc kernels
        do isite=1,nsite
            do imat=1,nmat
                do jmat=1,nmat
                    do kmat=1,nmat
                        pf2(imat,jmat,isite,ispin)=pf2(imat,jmat,isite,ispin)-gamtmp&
                            *(xmat(imat,kmat,isite,idim)*pf1(kmat,jmat,isite,jspin)&
                            -xmat(kmat,jmat,isite,idim)*pf1(imat,kmat,isite,jspin))
                    end do
                end do
            end do
        end do
       !$acc end kernels
    end do

END SUBROUTINE Test_Multiply_Dirac_dagger_device_reorg1

SUBROUTINE Test_Multiply_Dirac_dagger_device_reorg2(xmat,pf1,pf2)
    use compiletimeconstants
    use gammamatrix
    implicit none

    double complex :: gamtmp,tmp

    double complex,intent(in) :: xmat(-(nmargin-1):nsite+nmargin,1:nmat,1:nmat,1:ndim)
    double complex,intent(in) :: pf1(-(nmargin-1):nsite+nmargin,1:nmat,1:nmat,1:nspin)
    double complex,intent(out) :: pf2(-(nmargin-1):nsite+nmargin,1:nmat,1:nmat,1:nspin)
    !******************
    integer :: imat,jmat,kmat
    integer :: idim
    integer :: ispin,jspin,kspin
    integer :: isite
    integer :: countr
    !$acc declare present(xmat,pf1,pf2)

    !$acc kernels
    pf2=(0d0,0d0)
    !$acc end kernels

    do countr=1,144
        ispin=gamispin(countr)
        jspin=gamjspin(countr)
        idim=gamidim(countr)
        gamtmp=dconjg(gamgam(countr))
        do kmat=1,nmat
            !$acc parallel loop independent collapse(2)
            do imat=1,nmat
                do jmat=1,nmat
                    !$acc loop vector
                    do isite=1,nsite
                        pf2(isite,imat,jmat,ispin)=pf2(isite,imat,jmat,ispin)-gamtmp&
                            *(xmat(isite,imat,kmat,idim)*pf1(isite,kmat,jmat,jspin)&
                            -xmat(isite,kmat,jmat,idim)*pf1(isite,imat,kmat,jspin))
                    end do
                end do
            end do
          !$acc end parallel
        end do

    end do

END SUBROUTINE Test_Multiply_Dirac_dagger_device_reorg2


! These memory reduced parts are currently in a test status.
! They are working, but not more efficient.
SUBROUTINE reduce_mem_xmat(xmatin,xmatout)
    use compiletimeconstants
    implicit none
    double complex, intent(in) :: xmatin(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    double complex, intent(out) :: xmatout(1:(nmat*(nmat+1)/2),1:ndim,-(nmargin-1):nsite+nmargin)
    !$acc declare present(xmatin,xmatout)
    integer :: isite,imat,jmat,idim,vmat

    !$acc kernels
    do isite=-(nmargin-1),nsite+nmargin
        do idim=1,ndim
            ! diagonal part first
            do vmat=1,nmat
                xmatout(vmat,idim,isite)=xmatin(vmat,vmat,idim,isite)
            end do
            jmat=1
            imat=2
            do vmat=(nmat+1),(nmat*(nmat+1)/2)
                xmatout(vmat,idim,isite)=xmatin(imat,jmat,idim,isite)
                jmat=jmat+1
                if(imat==jmat) then
                    jmat=1
                    imat=imat+1
                end if
            end do
        end do
    end do
    !$acc end kernels
END SUBROUTINE reduce_mem_xmat

SUBROUTINE expand_mem_xmat(xmatin,xmatout)
    use compiletimeconstants
    implicit none
    double complex, intent(out) :: xmatout(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    double complex, intent(in) :: xmatin(1:(nmat*(nmat+1)/2),1:ndim,-(nmargin-1):nsite+nmargin)
    !$acc declare present(xmatin,xmatout)
    integer :: isite,imat,jmat,idim,vmat

    !$acc kernels
    do isite=-(nmargin-1),nsite+nmargin
        do idim=1,ndim
            ! diagonal part first
            vmat=1
            do vmat=1,nmat
                xmatout(vmat,vmat,idim,isite)=xmatin(vmat,idim,isite)
            end do
            jmat=1
            imat=2
            do vmat=(nmat+1),(nmat*(nmat+1)/2)
                xmatout(imat,jmat,idim,isite)=xmatin(vmat,idim,isite)
                xmatout(jmat,imat,idim,isite)=conjg(xmatin(vmat,idim,isite))
                jmat=jmat+1
                if(imat==jmat) then
                    jmat=1
                    imat=imat+1
                end if
            end do
        end do
    end do
    !$acc end kernels
END SUBROUTINE expand_mem_xmat

SUBROUTINE Multiply_Dirac_device_memred(xmat,pf1,pf2,var)
    use compiletimeconstants
    use gammamatrix
    implicit none
    integer var
    double complex, intent(in) :: xmat(1:(nmat*(nmat+1)/2),1:ndim,-(nmargin-1):nsite+nmargin)
    double complex, intent(in) :: pf1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex, intent(out) :: pf2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double precision :: lattice_spacing
    double complex :: gamtmp
    !******************
    integer :: imat,jmat,kmat,vmat
    integer :: idim
    integer :: ispin,jspin,kspin
    integer :: isite
    integer :: countr
    !double complex :: tmpsum ! private

    !$acc declare present(xmat,pf1,pf2)

    !**************************
    !**************************
    !***  interaction part  ***
    !**************************
    !**************************
    lattice_spacing=1d0
    !$acc kernels
    pf2=(0d0,0d0)
    !$acc end kernels

    do countr=1,144
        ispin=gamispin(countr)
        jspin=gamjspin(countr)
        idim=gamidim(countr)
        gamtmp=dconjg(gamgam(countr))
        if(var==1) then
            !$acc kernels
            do isite=1,nsite
                do kmat=1,nmat
                    do vmat=1,nmat
                        pf2(vmat,kmat,ispin,isite)=pf2(vmat,kmat,ispin,isite)-dcmplx(lattice_spacing)*gamtmp&
                            *xmat(vmat,idim,isite)*pf1(vmat,kmat,jspin,isite)
                        pf2(kmat,vmat,ispin,isite)=pf2(kmat,vmat,ispin,isite)+dcmplx(lattice_spacing)*gamtmp&
                            *pf1(kmat,vmat,jspin,isite)*xmat(vmat,idim,isite)
                    end do
                    jmat=1
                    imat=2
                    do vmat=(nmat+1),(nmat*(nmat+1)/2)
                        ! imat=vmat-(jmat-1)
                        pf2(imat,kmat,ispin,isite)=pf2(imat,kmat,ispin,isite)-dcmplx(lattice_spacing)*gamtmp&
                            *xmat(vmat,idim,isite)*pf1(jmat,kmat,jspin,isite)
                        pf2(jmat,kmat,ispin,isite)=pf2(jmat,kmat,ispin,isite)-dcmplx(lattice_spacing)*gamtmp&
                            *conjg(xmat(vmat,idim,isite))*pf1(imat,kmat,jspin,isite)
                        pf2(kmat,jmat,ispin,isite)=pf2(kmat,jmat,ispin,isite)+dcmplx(lattice_spacing)*gamtmp&
                            *pf1(kmat,imat,jspin,isite)*xmat(vmat,idim,isite)
                        pf2(kmat,imat,ispin,isite)=pf2(kmat,imat,ispin,isite)+dcmplx(lattice_spacing)*gamtmp&
                            *pf1(kmat,jmat,jspin,isite)*conjg(xmat(vmat,idim,isite))
                        jmat=jmat+1
                        if(imat==jmat) then
                            jmat=1
                            imat=imat+1
                        end if
                    end do
                end do
            end do
        !$acc end kernels
        end if
        if(var==2) then
            !$acc kernels
            do isite=1,nsite
                do kmat=1,nmat
                    do vmat=1,nmat
                        pf2(vmat,kmat,ispin,isite)=pf2(vmat,kmat,ispin,isite)-dcmplx(lattice_spacing)*gamtmp&
                            *xmat(vmat,idim,isite)*pf1(vmat,kmat,jspin,isite)
                    end do
                end do
            end do
            !$acc end kernels
            !$acc kernels
            do isite=1,nsite
                do kmat=1,nmat
                    do vmat=1,nmat
                        pf2(kmat,vmat,ispin,isite)=pf2(kmat,vmat,ispin,isite)+dcmplx(lattice_spacing)*gamtmp&
                            *pf1(kmat,vmat,jspin,isite)*xmat(vmat,idim,isite)
                    end do
                end do
            end do
            !$acc end kernels
            !$acc kernels
            do isite=1,nsite
                do kmat=1,nmat
                    do imat=2,nmat
                        do jmat=1,(imat-1)
                            vmat=((imat-1)*(imat-2))/2+jmat+nmat
                            pf2(imat,kmat,ispin,isite)=pf2(imat,kmat,ispin,isite)-dcmplx(lattice_spacing)*gamtmp&
                                *xmat(vmat,idim,isite)*pf1(jmat,kmat,jspin,isite)
                        end do
                    end do
                end do
            end do
            !$acc end kernels
            !$acc kernels
            do isite=1,nsite
                do kmat=1,nmat
                    do jmat=1,(nmat-1)
                        do imat=jmat+1,nmat
                            vmat=((imat-1)*(imat-2))/2+jmat+nmat
                            pf2(jmat,kmat,ispin,isite)=pf2(jmat,kmat,ispin,isite)-dcmplx(lattice_spacing)*gamtmp&
                                *conjg(xmat(vmat,idim,isite))*pf1(imat,kmat,jspin,isite)
                        end do
                    end do
                end do
            end do
            !$acc end kernels
            !$acc kernels
            do isite=1,nsite
                do kmat=1,nmat
                    do jmat=1,(nmat-1)
                        do imat=jmat+1,nmat
                            vmat=((imat-1)*(imat-2))/2+jmat+nmat
                            pf2(kmat,jmat,ispin,isite)=pf2(kmat,jmat,ispin,isite)+dcmplx(lattice_spacing)*gamtmp&
                                *pf1(kmat,imat,jspin,isite)*xmat(vmat,idim,isite)
                        end do
                    end do
                end do
            end do
            !$acc end kernels
            !$acc kernels
            do isite=1,nsite
                do kmat=1,nmat
                    do imat=2,nmat
                        do jmat=1,(imat-1)
                            vmat=((imat-1)*(imat-2))/2+jmat+nmat
                            pf2(kmat,imat,ispin,isite)=pf2(kmat,imat,ispin,isite)+dcmplx(lattice_spacing)*gamtmp&
                                *pf1(kmat,jmat,jspin,isite)*conjg(xmat(vmat,idim,isite))
                        end do
                    end do
                end do
            end do
        !$acc end kernels
        end if
        if(var==3) then
            !$acc kernels
            do isite=1,nsite
                do kmat=1,nmat
                    do vmat=1,nmat
                        pf2(vmat,kmat,ispin,isite)=pf2(vmat,kmat,ispin,isite)-dcmplx(lattice_spacing)*gamtmp&
                            *xmat(vmat,idim,isite)*pf1(vmat,kmat,jspin,isite)
                        pf2(kmat,vmat,ispin,isite)=pf2(kmat,vmat,ispin,isite)+dcmplx(lattice_spacing)*gamtmp&
                            *pf1(kmat,vmat,jspin,isite)*xmat(vmat,idim,isite)
                    end do
                end do
                do kmat=1,nmat
                    do jmat=1,(nmat-1)
                        do imat=jmat+1,nmat
                            vmat=((imat-1)*(imat-2))/2+jmat+nmat
                            pf2(imat,kmat,ispin,isite)=pf2(imat,kmat,ispin,isite)-dcmplx(lattice_spacing)*gamtmp&
                                *xmat(vmat,idim,isite)*pf1(jmat,kmat,jspin,isite)
                            pf2(jmat,kmat,ispin,isite)=pf2(jmat,kmat,ispin,isite)-dcmplx(lattice_spacing)*gamtmp&
                                *conjg(xmat(vmat,idim,isite))*pf1(imat,kmat,jspin,isite)
                            pf2(kmat,jmat,ispin,isite)=pf2(kmat,jmat,ispin,isite)+dcmplx(lattice_spacing)*gamtmp&
                                *pf1(kmat,imat,jspin,isite)*xmat(vmat,idim,isite)
                            pf2(kmat,imat,ispin,isite)=pf2(kmat,imat,ispin,isite)+dcmplx(lattice_spacing)*gamtmp&
                                *pf1(kmat,jmat,jspin,isite)*conjg(xmat(vmat,idim,isite))
                        end do
                    end do
                end do
            end do
        !$acc end kernels
        end if
    end do

END SUBROUTINE Multiply_Dirac_device_memred

SUBROUTINE Test_Multiply_Dirac_dagger_device_cublas(xmat,pf1,pf2,xptr_d,pf1ptr_d,pf2ptr_d)
    !********************************(needed for cublas version)
    use cublasinterface
    use compiletimeconstants
    use gammamatrix
    implicit none

    double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    double complex, intent(in) :: pf1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex, intent(out) :: pf2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex :: gamtmp,tmp

    double complex :: xmat2(1:nmat,1:nmat,-(nmargin-1):nsite+nmargin,1:ndim)
    double complex :: pf12(1:nmat,1:nmat,-(nmargin-1):nsite+nmargin,1:nspin)
    double complex :: pf22(1:nmat,1:nmat,-(nmargin-1):nsite+nmargin,1:nspin)
    !******************
    integer :: imat,jmat,kmat
    integer :: idim
    integer :: ispin,jspin,kspin
    integer :: isite
    integer :: countr
    !$acc declare present(xmat,pf1,pf2)
    type(c_devptr), device :: xptr_d(nsite,ndim)
    type(c_devptr), device :: pf1ptr_d(nsite,nspin)
    type(c_devptr), device :: pf2ptr_d(nsite,nspin)

    !$acc kernels
    pf2=(0d0,0d0)
    !$acc end kernels

    do countr=1,144
        ispin=gamispin(countr)
        jspin=gamjspin(countr)
        idim=gamidim(countr)
        gamtmp=dconjg(gamgam(countr))
            !callZgemmBatched(pf2,xmat,pf1,gamtmp,dcmplx(1.d0),ispin,jspin,idim)
            call multiply_cublas_pointer(xptr_d,pf1ptr_d,pf2ptr_d,gamtmp,dcmplx(1.d0),ispin,jspin,idim)
       end do

END SUBROUTINE Test_Multiply_Dirac_dagger_device_cublas


SUBROUTINE Test_Multiply_Dirac_dagger_device_cublas_stream(pf2,xptr_d,pf1ptr_d,pf2ptr_d)
    !********************************(needed for cublas version)
    use cublasinterface
    use compiletimeconstants
    use gammamatrix
    implicit none

    double complex, intent(out) :: pf2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    !$acc declare present(pf2)
    type(c_devptr), device :: xptr_d(nsite,ndim)
    type(c_devptr), device :: pf1ptr_d(nsite,nspin)
    type(c_devptr), device :: pf2ptr_d(nsite,nspin)

    !$acc kernels
    pf2=(0d0,0d0)
    !$acc end kernels
    call multiply_cublas_pointer_streams(xptr_d,pf1ptr_d,pf2ptr_d,dcmplx(1.d0),.TRUE.)

END SUBROUTINE Test_Multiply_Dirac_dagger_device_cublas_stream

SUBROUTINE diff_vector_test(diff,vect1,vect2)
    use compiletimeconstants
    implicit none

    double complex, intent(in) :: vect1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex, intent(in) :: vect2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double precision, intent(out) :: diff
    integer isite,ispin,imat,jmat
    !$acc declare device_resident(vect1,vect2)
    diff=0d0
    !$acc kernels
    do isite=1,nsite
        do ispin=1,nspin
            do imat=1,nmat
                do jmat=1,nmat
                    diff=diff+abs(vect1(jmat,imat,ispin,isite)-vect2(jmat,imat,ispin,isite))
                end do
            end do
        end do
    end do
    !$acc end kernels
end SUBROUTINE

subroutine check_mintime(time)
  implicit none
  real time
  real,save :: mtime=100000000
  if(time<mtime) then
    mtime=time
    print*, "=================> fastest so far"
  end if
end subroutine check_mintime

program timing_mult
    use compiletimeconstants
    use dirac_operator
    use outputstreamnumbers
    use mtmod !Mersenne twistor
    use cublasinterface
    use cublas
    use iso_c_binding
    implicit none
  
    double complex :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    ! double complex, dimension(:,:,:,:), allocatable :: xmat
    double precision alpha(1:nmat)

    integer nfuzzy
    !Gamma matrices
    double complex Gamma10d(1:ndim,1:nspin,1:nspin)
    double precision temperature
    double precision flux

    ! Run parameters
    integer, parameter :: testhost=1 !switch of on slow systems
    integer, parameter :: runnumber=10
    integer, parameter :: numvar=9
    ! This estimate has to be checked.
    double precision :: flopspermult

    integer :: nbc !boundary condition for fermions; 0 -> pbc, 1 -> apbc
    integer:: nbmn ! 0 -> BFSS, 1 -> BMN
  
    double complex :: phase(1:nmat,1:nmat,1:2)
    double complex :: Gam123(1:nspin,1:nspin)
    !$acc declare device_resident(phase,Gam123)
  
    double complex testvect1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex testvect2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex testvect0(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
  
    double complex :: testvect_d0(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex :: testvect_d1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex :: testvect_d2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex :: testvect_d0_reorg1(1:nmat,1:nmat,-(nmargin-1):nsite+nmargin,1:nspin)
    double complex :: testvect_d2_reorg1(1:nmat,1:nmat,-(nmargin-1):nsite+nmargin,1:nspin)
    double complex :: testvect_d0_reorg2(-(nmargin-1):nsite+nmargin,1:nmat,1:nmat,1:nspin)
    double complex :: testvect_d2_reorg2(-(nmargin-1):nsite+nmargin,1:nmat,1:nmat,1:nspin)
    !$acc declare device_resident(testvect_d0,testvect_d2,testvect_d1)
    !$acc declare device_resident(testvect_d0_reorg1,testvect_d2_reorg1)
    !$acc declare device_resident(testvect_d0_reorg2,testvect_d2_reorg2)

    double complex :: xmatreduced(1:(nmat*(nmat+1)/2),1:ndim,-(nmargin-1):nsite+nmargin)
    double complex :: xmatexpand(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    double complex :: xmatreorg1(1:nmat,1:nmat,-(nmargin-1):nsite+nmargin,1:ndim)
    double complex :: xmatreorg2(-(nmargin-1):nsite+nmargin,1:nmat,1:nmat,1:ndim)

    real :: start_time,stop_time,time1,time2,time1_cud,time2_cud
    integer :: counter,isite,idim,imat,jmat,ivar
    double complex :: tmp
    real mintime

    type(c_devptr), device :: xptr_d(nsite,ndim)
    type(c_devptr), device :: pf1ptr_d(nsite,nspin)
    type(c_devptr), device :: pf2ptr_d(nsite,nspin)

    nbc=1 !boundary condition for fermions; 0 -> pbc, 1 -> apbc
    nbmn=0 ! 0 -> BFSS, 1 -> BMN
    temperature=1d0*nbc
    flux=0.5d0*temperature
    nfuzzy=1
        ! This estimate has to be checked.
    flopspermult=dble((6+1)*2)*dble(144)*dble(nsite)*dble(nmat)*dble(nmat)*dble(nmat)*dble(runnumber)

     !allocate(xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin))

     
    call MakeGamma(Gamma10d)
 

    testvect0=(1d0,1d0)
    ! more or less random field xmat
    xmat=(0d0,0d0)
    do isite=-(nmargin-1),nsite+nmargin
        do idim=1,ndim
            do imat=1,nmat
                do jmat=1,nmat
                    xmat(imat,jmat,idim,isite)=CMPLX(imat+2.943*jmat,idim+0.02*isite)
                    xmat(jmat,imat,idim,isite)=conjg(xmat(imat,jmat,idim,isite))
                end do
            end do
        end do
    end do
    ! This is in case you want to have some diagonal field only
    !    do isite=-(nmargin-1),nsite+nmargin
    !        do idim=1,ndim
    !            do imat=1,nmat
    !                  xmat(imat,imat,idim,isite)=CMPLX(imat+2.943*imat,idim+0.02*isite)
    !            end do
    !        end do
    !    end do

    !$acc data &
    !$acc copyin(testvect0,temperature,xmat,alpha,GAMMA10d,nbmn,flux,nbc,xmatreduced,xmatexpand,xmatreorg1,xmatreorg2)

    ! do setup of all the data.
    call setup_cublas()
    call setup_cublas_pointers_xmat(xmat,xptr_d)
    call setup_cublas_pointers_pf(testvect_d0,pf1ptr_d)
    call setup_cublas_pointers_pf(testvect_d2,pf2ptr_d)
    call setup_data_device(alpha,flux,GAMMA10d,phase,Gam123,temperature)
    call cudaSetupTimer()

    testvect1=0d0
    if(testhost==1) then
        call cudaTimerStart()
        call cpu_time(start_time)
        do counter=1,runnumber
            call Multiply_Dirac_dagger(temperature,xmat,alpha,&
                testvect0,testvect1,GAMMA10d,nbmn,flux)
        end do
        call cpu_time(stop_time)
        call cudaTimerStop(time1_cud)
    endif
    time1=stop_time - start_time

    print *, "Host time:", &
        time1, "seconds ", dble(flopspermult)/time1, "Flops"
    call check_mintime(time1)
    !$acc kernels
    testvect_d0=testvect0
    !$acc end kernels
    !$acc kernels
    testvect_d1=testvect1
    !$acc end kernels
    call cudaTimerStart()
    call cpu_time(start_time)
    do counter=1,runnumber
        call Multiply_Dirac_dagger_device(temperature,xmat,phase,Gam123,nbmn,&
            testvect_d0,testvect_d2)
    end do
    call cpu_time(stop_time)
    call cudaTimerStop(time2_cud)
    time2=stop_time - start_time
    print *, "Device time:", &
        time2, "seconds ", flopspermult/time2, "Flops"
    print *,"relative timing",time1/time2," ",time1_cud/time2_cud
    call check_mintime(time2)
    !$acc kernels
    testvect2=testvect_d2
    !$acc end kernels
    !$acc kernels
    tmp=Sum(abs(testvect_d2))
    !$acc end kernels
    print *,"vect ",tmp," ", Sum(abs(testvect1)), " err ",Sum(abs(testvect2-testvect1)),&
        abs(tmp-Sum(abs(testvect1)))

    call cudaTimerStart()
    call cpu_time(start_time)
    do counter=1,runnumber
        call Multiply_Dirac_dagger_device_cuda(temperature,xmat,phase,Gam123,nbmn,&
            testvect_d0,testvect_d2,xptr_d,pf1ptr_d,pf2ptr_d)
    end do
    call cpu_time(stop_time)
    call cudaTimerStop(time2_cud)
    time2=stop_time - start_time
    call check_cublas()
    print *, "Device time cublas:", &
        time2, "seconds ", flopspermult/time2, "Flops"
    print *,"relative timing",time1/time2," ",time1_cud/time2_cud
    !$acc kernels
    testvect2=testvect_d2
    tmp=Sum(abs(testvect_d2))
    !$acc end kernels
    print *,"vect ",tmp," ", Sum(abs(testvect1)), " err ",Sum(abs(testvect2-testvect1)),&
        abs(tmp-Sum(abs(testvect1)))
    call check_mintime(time2)

    testvect1=0d0
    if(testhost==1) then
        call cudaTimerStart()
        call cpu_time(start_time)
        do counter=1,runnumber
            call Multiply_Dirac(temperature,xmat,alpha,&
                testvect0,testvect1,GAMMA10d,nbmn,flux)
        end do
        call cpu_time(stop_time)
        call cudaTimerStop(time1_cud)
    endif
    time1=stop_time - start_time
    print *, "Host time:", &
        time1, "seconds", flopspermult/time1, "Flops"

    call cudaTimerStart()
    call cpu_time(start_time)
    do counter=1,runnumber
        call Multiply_Dirac_device(temperature,xmat,phase,Gam123,nbmn,&
            testvect_d0,testvect_d2)
    end do
    call cpu_time(stop_time)
    call cudaTimerStop(time2_cud)
    time2=stop_time - start_time
    print *, "Device time:", &
        time2, "seconds", flopspermult/time2, "Flops"
    print *,"relative timing",time1/time2," ",time1_cud/time2_cud
    !$acc kernels
    testvect2=testvect_d2
    tmp=Sum(abs(testvect_d2))
    !$acc end kernels
    print *,"vect ",tmp," ", Sum(abs(testvect1)), " err ",Sum(abs(testvect2-testvect1))
    call check_mintime(time2)

    call cudaTimerStart()
    call cpu_time(start_time)
    do counter=1,runnumber
        call Multiply_Dirac_device_cuda(temperature,xmat,phase,Gam123,nbmn,&
            testvect_d0,testvect_d2,xptr_d,pf1ptr_d,pf2ptr_d)
    end do
    call cpu_time(stop_time)
    call cudaTimerStop(time2_cud)
    time2=stop_time - start_time
    call check_cublas()
    call finish_cublas()
    print *, "Device time cublas:", &
        time2, "seconds", flopspermult/time2, "Flops"
    print *,"relative timing",time1/time2," ",time1_cud/time2_cud
    !$acc kernels
    testvect2=testvect_d2
    tmp=Sum(abs(testvect_d2))
    !$acc end kernels
    print *,"vect ",tmp," ", Sum(abs(testvect1)), " err ",Sum(abs(testvect2-testvect1)),&
        abs(tmp-Sum(abs(testvect1)))
    call check_mintime(time2)
  
    do ivar=1,numvar
        print *, "Variant ", ivar
        call cpu_time(start_time)
        do counter=1,runnumber
            call Test_Multiply_Dirac_dagger_device(xmat,testvect_d0,testvect_d2,ivar)
        end do
        call cpu_time(stop_time)
        time2=stop_time - start_time
        print *, "Device time:", &
            time2, "seconds", flopspermult/time2, "Flops"
        print *,"relative timing",time1/time2
        if(ivar==1) then
            !$acc kernels
            testvect1=testvect_d2
            !$acc end kernels
        else
            !$acc kernels
            testvect2=testvect_d2
            !$acc end kernels
            print *,"diff to var1 ",Sum(abs(testvect2-testvect1))
        endif
        call check_mintime(time2)
    end do

    call setup_cublas()
    call setup_cublas_pointers_xmat(xmat,xptr_d)
    call setup_cublas_pointers_pf(testvect_d0,pf1ptr_d)
    call setup_cublas_pointers_pf(testvect_d2,pf2ptr_d)
    call cpu_time(start_time)
    do counter=1,runnumber
        call Test_Multiply_Dirac_dagger_device_cublas(xmat,testvect_d0,testvect_d2,xptr_d,pf1ptr_d,pf2ptr_d)
    end do
    call cpu_time(stop_time)
    time2=stop_time - start_time
    print *, "Device time cublas variant:", &
        time2, "seconds", flopspermult/time2, "Flops"
    print *,"relative timing",time1/time2
        !$acc kernels
    testvect2=testvect_d2
    tmp=Sum(abs(testvect_d2))
    !$acc end kernels
    print *,"vect ",tmp," ", Sum(abs(testvect1)), " diff var 1 ",Sum(abs(testvect2-testvect1))

    call cpu_time(start_time)
    do counter=1,runnumber
        call Test_Multiply_Dirac_dagger_device_cublas_stream(testvect_d2,xptr_d,pf1ptr_d,pf2ptr_d)
    end do
    call cpu_time(stop_time)
    time2=stop_time - start_time
    print *, "Device time cublas variant2:", &
        time2, "seconds", flopspermult/time2, "Flops"
    print *,"relative timing",time1/time2
        !$acc kernels
    testvect2=testvect_d2
    tmp=Sum(abs(testvect_d2))
    !$acc end kernels
    print *,"vect ",tmp," ", Sum(abs(testvect1)), " diff var 1 ",Sum(abs(testvect2-testvect1))

    call check_mintime(time2)
    call finish_cublas()



    print *, "Memory reduced variant "
    call reduce_mem_xmat(xmat,xmatreduced)
    call expand_mem_xmat(xmatreduced,xmatexpand)
    !$acc update host(xmatexpand)
    print *,"xmatred err ",Sum(abs(xmatexpand))," ", Sum(abs(xmat))," ",Sum(abs(xmatexpand-xmat))
    do ivar=2,2 ! all other variants are too slow.
        print *, "Variant ", ivar
        call cpu_time(start_time)
        do counter=1,runnumber
            call Multiply_Dirac_device_memred(xmatreduced,testvect_d0,testvect_d2,ivar)
        end do
        call cpu_time(stop_time)
        time2=stop_time - start_time
        print *, "Device time:", &
            time2, "seconds", flopspermult/time2, "Flops"
        print *,"relative timing",time1/time2
        !$acc kernels
        testvect2=testvect_d2
        tmp=Sum(abs(testvect_d2))
        !$acc end kernels
        print *,"vect ",tmp," ", Sum(abs(testvect1)), " err ",Sum(abs(testvect2-testvect1))
        call check_mintime(time2)
    end do

    print *, "Memory reorganized variant 1"
    call reorder_fields_xmat1(xmat,xmatreorg1)
    !$acc update host(xmatreorg1)
    call reorder_fields_pf1(testvect_d0,testvect_d0_reorg1)
    call cpu_time(start_time)
    do counter=1,runnumber
        call Test_Multiply_Dirac_dagger_device_reorg1(xmatreorg1,testvect_d0_reorg1,testvect_d2_reorg1)
    end do
    call cpu_time(stop_time)
    time2=stop_time - start_time
    print *, "Device time:", &
        time2, "seconds", flopspermult/time2, "Flops"
    print *,"relative timing",time1/time2
    call inv_reorder_fields_pf1(testvect_d2_reorg1,testvect_d2)
    !$acc kernels
    testvect2=testvect_d2
    tmp=Sum(abs(testvect_d2_reorg1))
    !$acc end kernels
    print *,"vect ",tmp," ", Sum(abs(testvect1)), " err ",Sum(abs(testvect2-testvect1))
    call check_mintime(time2)

    print *, "Memory reorganized variant 2"
    call reorder_fields_xmat2(xmat,xmatreorg2)
    !$acc update host(xmatreorg2)
    call reorder_fields_pf2(testvect_d0,testvect_d0_reorg2)
    call cpu_time(start_time)
    do counter=1,runnumber
        call Test_Multiply_Dirac_dagger_device_reorg2(xmatreorg2,testvect_d0_reorg2,testvect_d2_reorg2)
    end do
    call cpu_time(stop_time)
    time2=stop_time - start_time
    print *, "Device time:", &
        time2, "seconds", flopspermult/time2, "Flops"
    print *,"relative timing",time1/time2
    call inv_reorder_fields_pf2(testvect_d2_reorg2,testvect_d2)
    !$acc kernels
    testvect2=testvect_d2
    tmp=Sum(abs(testvect_d2_reorg2))
    !$acc end kernels
    print *,"vect ",tmp," ", Sum(abs(testvect1)), " err ",Sum(abs(testvect2-testvect1))
    call check_mintime(time2)


   call cudaFinishTimer()
  !End test part
  !deallocate(xmat)
      !$acc end data
end program timing_mult
