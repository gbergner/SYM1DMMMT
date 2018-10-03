! ###################################################
!  Dirac operator of the theory used in CGM solver
!  In order to keep compatibility with Masamoris code: extra factor (-i) included!
!  pf2=(-i)*M*pf1, M: Dirac op
!  In order to simplify the operation: variables phase, Gam123 precomputed
! ###################################################


module dirac_operator
    implicit none
contains
    ! Boundary setup for up down lookup.
    SUBROUTINE set_boundary_device(nbc,pf)
        use compiletimeconstants
        implicit none
        double complex,intent(inout) :: pf(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        integer, intent(in) :: nbc
        !$acc declare device_resident(pf)
        integer :: ispin,imat,jmat,isite
        if(nbc.EQ.0)then!pbc
            !$acc kernels
            do ispin=1,nspin
                do jmat=1,nmat
                    do imat=1,nmat
                        do isite=-(nmargin-1),0
                            pf(imat,jmat,ispin,isite)=&
                                pf(imat,jmat,ispin,nsite+isite)
                        end do
                    end do
                end do
            end do
            do ispin=1,nspin
                do jmat=1,nmat
                    do imat=1,nmat
                        do isite=nsite+1,nsite+nmargin
                            pf(imat,jmat,ispin,isite)=&
                                pf(imat,jmat,ispin,isite-nsite)
                        end do
                    end do
                end do
            end do
           !$acc end kernels
           
        else if(nbc.EQ.1)then!apbc
             !$acc kernels
            do ispin=1,nspin
                do jmat=1,nmat
                    do imat=1,nmat
                        do isite=-(nmargin-1),0
                            pf(imat,jmat,ispin,isite)=&
                                pf(imat,jmat,ispin,nsite+isite)*(-1d0,0d0)
                        end do
                    end do
                end do
            end do
            do ispin=1,nspin
                do jmat=1,nmat
                    do imat=1,nmat
                        do isite=nsite+1,nsite+nmargin
                            pf(imat,jmat,ispin,isite)=&
                                pf(imat,jmat,ispin,isite-nsite)*(-1d0,0d0)
                        end do
                    end do
                end do
            end do
           !$acc end kernels
        end if
    END SUBROUTINE set_boundary_device

    ! Precompute variables for simpler Dirac operator application.
    ! phase(:,:,1) are just the exponential phase factors
    ! phase(:,:,2) next neighbours interaction in improved operator.
    SUBROUTINE setup_data_device(alpha,flux,GAMMA10d,phase,Gam123,temperature)
        use compiletimeconstants
        implicit none
        double precision, intent(in) :: alpha(1:nmat)
        double precision, intent(in) :: flux
        double complex, intent(in):: GAMMA10d(1:ndim,1:nspin,1:nspin)
        double precision, intent(in) :: temperature
          !$acc declare present (alpha,flux,GAMMA10d, temperature)
        double complex, intent(out) :: phase(1:nmat,1:nmat,1:2)
        double complex, intent(out) :: Gam123(1:nspin,1:nspin)
        double precision :: aij
        double complex :: Gam12(1:nspin,1:nspin)
        !$acc declare device_resident(phase,Gam123,aij,Gam12)
        integer jmat,imat,jspin,ispin,kspin
        double precision :: lattice_spacing
        lattice_spacing=1d0/temperature/dble(nsite)
        !$acc kernels
        do jmat=1,nmat
            do imat=1,nmat
                aij=alpha(imat)-alpha(jmat)
                aij=aij/dble(nsite)
                phase(imat,jmat,1)=dcmplx(dcos(aij))+(0d0,1d0)*dcmplx(dsin(aij))
                phase(imat,jmat,2)=phase(imat,jmat,1)*phase(imat,jmat,1);
            end do
        end do
        Gam12=(0d0,0d0)
        Gam123=(0d0,0d0)
        do jspin=1,nspin
            do ispin=1,nspin
                do kspin=1,nspin
                    Gam12(ispin,jspin)=Gam12(ispin,jspin)&
                        +Gamma10d(1,ispin,kspin)*Gamma10d(2,kspin,jspin)
                end do
            end do
        end do
        do jspin=1,nspin
            do ispin=1,nspin
                do kspin=1,nspin
                    Gam123(ispin,jspin)=Gam123(ispin,jspin)&
                        +Gam12(ispin,kspin)*Gamma10d(3,kspin,jspin)
                end do
            end do
        end do
        Gam123=lattice_spacing*Gam123*dcmplx(flux)*(0d0,-0.75d0)
         !$acc end kernels
    END SUBROUTINE setup_data_device

    ! Of the precomputed variables only the phase needs to be updated.
    SUBROUTINE update_data_device(alpha,phase)
        use compiletimeconstants
        implicit none
        double precision, intent(in) :: alpha(1:nmat)
        !$acc declare present (alpha)
        double complex, intent(out) :: phase(1:nmat,1:nmat,1:2)
        double precision :: aij
        !$acc declare device_resident(phase)
        integer jmat,imat,jspin,ispin,kspin
        !$acc kernels
        do jmat=1,nmat
            do imat=1,nmat
                aij=alpha(imat)-alpha(jmat)
                aij=aij/dble(nsite)
                phase(imat,jmat,1)=dcmplx(dcos(aij))+(0d0,1d0)*dcmplx(dsin(aij))
                phase(imat,jmat,2)=phase(imat,jmat,1)*phase(imat,jmat,1);
            end do
        end do
        !$acc end kernels
    END SUBROUTINE update_data_device

    SUBROUTINE Multiply_Dirac_dagger_device(temperature,xmat,phase,Gam123,nbmn,&
        pf1,pf2)
        use cublasinterface
        use compiletimeconstants
        use gammamatrix
        implicit none

        integer, intent(in) :: nbmn
        double precision, intent(in) :: temperature
        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double complex, intent(in) :: phase(1:nmat,1:nmat,1:2)
        double complex, intent(in) :: Gam123(1:nspin,1:nspin)
        double complex, intent(in) :: pf1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex, intent(out) :: pf2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double precision :: lattice_spacing
        double complex :: gamtmp
        !******************
        integer :: imat,jmat,kmat
        integer :: idim
        integer :: ispin,jspin,kspin
        integer :: isite
        integer :: countr
        !$acc declare present(temperature,xmat,phase,Gam123,pf1,pf2)

        !**********************
        !**********************
        !***  kinetic part  ***
        !**********************
        !**********************
        !$acc kernels
        do isite=1,nsite
            do ispin=1,nspin
                do jmat=1,nmat
                    do imat=1,nmat
                        pf2(imat,jmat,ispin,isite)=(0d0,0d0)
                    end do
                end do
            end do
        end do
        !$acc end kernels
        !********************
        !*** Naive Action ***
        !********************
        if(nimprove.EQ.0)then
            !$acc parallel loop &
            !$acc collapse(4) &
            !$acc gang vector &
            !$acc independent
            do isite=1,nsite
                do ispin=1,8
                    do jmat=1,nmat
                        do imat=1,nmat
                            pf2(imat,jmat,ispin+8,isite)=&
                                pf2(imat,jmat,ispin+8,isite)&
                                -(1d0,0d0)*phase(imat,jmat,1)*pf1(imat,jmat,ispin,isite+1)&
                                +(1d0,0d0)*pf1(imat,jmat,ispin,isite)
                            pf2(imat,jmat,ispin,isite)=&
                                pf2(imat,jmat,ispin,isite)&
                                +(1d0,0d0)*dconjg(phase(imat,jmat,1))*pf1(imat,jmat,ispin+8,isite-1)&
                                -(1d0,0d0)*pf1(imat,jmat,ispin+8,isite)
                        end do
                    end do
                end do
            end do
           !$acc end parallel loop
           !***********************
           !*** Improved Action ***
           !***********************
        else if(nimprove.EQ.1)then
            !$acc parallel loop &
            !$acc collapse(4) &
            !$acc gang vector &
            !$acc independent
            do isite=1,nsite
                do ispin=1,8
                    do jmat=1,nmat
                        do imat=1,nmat
                            pf2(imat,jmat,ispin+8,isite)=&
                                pf2(imat,jmat,ispin+8,isite)&
                                +(0.5d0,0d0)*phase(imat,jmat,2)*pf1(imat,jmat,ispin,isite+2)&
                                -(2d0,0d0)*phase(imat,jmat,1)*pf1(imat,jmat,ispin,isite+1)&
                                +(1.5d0,0d0)*pf1(imat,jmat,ispin,isite)

                            pf2(imat,jmat,ispin,isite)=&
                                pf2(imat,jmat,ispin,isite)&
                                -(0.5d0,0d0)*dconjg(phase(imat,jmat,2))*pf1(imat,jmat,ispin+8,isite-2)&
                                +(2d0,0d0)*dconjg(phase(imat,jmat,1))*pf1(imat,jmat,ispin+8,isite-1)&
                                -(1.5d0,0d0)*pf1(imat,jmat,ispin+8,isite)
                        end do
                    end do
                end do
            end do
           !$acc end parallel loop
        end if
        !**************************
        !**************************
        !***  interaction part  ***
        !**************************
        !**************************
        lattice_spacing=1d0/temperature/dble(nsite)

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
                            pf2(imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)-dcmplx(lattice_spacing)*gamtmp&
                                *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                                -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
                        end do
                    end do
                end do
            end do
           !$acc end kernels
        end do

        !******************************
        !******************************
        !*** Plane wave deformation ***
        !******************************
        !******************************
        if(nbmn.EQ.1)then
            !$acc parallel loop &
            !$acc collapse(4) &
            !$acc gang vector &
            !$acc independent
            do isite=1,nsite
                do ispin=1,nspin
                    do imat=1,nmat
                        do jmat=1,nmat
                            !$acc loop seq
                            do jspin=1,nspin
                                pf2(imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)&
                                    +dconjg(Gam123(jspin,ispin))*pf1(imat,jmat,jspin,isite)
                            end do
                        end do
                    end do
                end do
            end do
          !$acc end parallel loop
        end if
        return
  
    END SUBROUTINE Multiply_Dirac_dagger_device

    SUBROUTINE Multiply_Dirac_device(temperature,xmat,phase,Gam123,nbmn,&
        pf1,pf2)
        use cublasinterface
        use compiletimeconstants
        use gammamatrix
        implicit none
        integer, intent(in) :: nbmn
        double precision, intent(in) :: temperature
        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double complex, intent(in) :: phase(1:nmat,1:nmat,1:2)
        double complex, intent(in) :: Gam123(1:nspin,1:nspin)
        double complex, intent(in) :: pf1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex, intent(out) :: pf2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double precision :: lattice_spacing
        double complex :: gamtmp
        !******************
        integer :: imat,jmat,kmat
        integer :: idim
        integer :: ispin,jspin,kspin
        integer :: isite
        integer :: countr

        !$acc declare present(temperature,xmat,phase,Gam123,pf1,pf2)


        !**********************
        !**********************
        !***  kinetic part  ***
        !**********************
        !**********************
        !$acc parallel loop &
        !$acc independent collapse(4) gang vector
        do isite=1,nsite
            do ispin=1,nspin
                do jmat=1,nmat
                    do imat=1,nmat
                        pf2(imat,jmat,ispin,isite)=(0d0,0d0)
                    end do
                end do
            end do
        end do
        !$acc end parallel loop
        !********************
        !*** Naive Action ***
        !********************
        if(nimprove.EQ.0)then
            !$acc parallel loop &
            !$acc collapse(4) &
            !$acc gang vector &
            !$acc independent
            do isite=1,nsite
                do ispin=1,8
                    do jmat=1,nmat
                        do imat=1,nmat
                            pf2(imat,jmat,ispin+8,isite)=&
                                pf2(imat,jmat,ispin+8,isite)&
                                +(1d0,0d0)*phase(imat,jmat,1)*pf1(imat,jmat,ispin,isite+1)&
                                -(1d0,0d0)*pf1(imat,jmat,ispin,isite)
                            pf2(imat,jmat,ispin,isite)=&
                                pf2(imat,jmat,ispin,isite)&
                                -(1d0,0d0)*dconjg(phase(imat,jmat,1))*pf1(imat,jmat,ispin+8,isite-1)&
                                +(1d0,0d0)*pf1(imat,jmat,ispin+8,isite)
                        end do
                    end do
                end do
            end do
           !$acc end parallel loop
           !***********************
           !*** Improved Action ***
           !***********************
        else if(nimprove.EQ.1)then
            !$acc parallel loop &
            !$acc collapse(4) &
            !$acc gang vector &
            !$acc independent
            do isite=1,nsite
                do ispin=1,8
                    do jmat=1,nmat
                        do imat=1,nmat
                            pf2(imat,jmat,ispin+8,isite)=&
                                pf2(imat,jmat,ispin+8,isite)&
                                -(0.5d0,0d0)*phase(imat,jmat,2)*pf1(imat,jmat,ispin,isite+2)&
                                +(2d0,0d0)*phase(imat,jmat,1)*pf1(imat,jmat,ispin,isite+1)&
                                -(1.5d0,0d0)*pf1(imat,jmat,ispin,isite)

                            pf2(imat,jmat,ispin,isite)=&
                                pf2(imat,jmat,ispin,isite)&
                                +(0.5d0,0d0)*dconjg(phase(imat,jmat,2))*pf1(imat,jmat,ispin+8,isite-2)&
                                -(2d0,0d0)*dconjg(phase(imat,jmat,1))*pf1(imat,jmat,ispin+8,isite-1)&
                                +(1.5d0,0d0)*pf1(imat,jmat,ispin+8,isite)
                        end do
                    end do
                end do
            end do
           !$acc end parallel loop
        end if
        !**************************
        !**************************
        !***  interaction part  ***
        !**************************
        !**************************
        lattice_spacing=1d0/temperature/dble(nsite)
        ! this is not the most efficient variant (see main_time)
        do countr=1,144
            ispin=gamispin(countr)
            jspin=gamjspin(countr)
            idim=gamidim(countr)
            gamtmp=gamgam(countr)
            !$acc kernels
            do isite=1,nsite
                do imat=1,nmat
                    do jmat=1,nmat
                        do kmat=1,nmat
                            pf2(imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)-dcmplx(lattice_spacing)*gamtmp&
                                *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                                -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
                        end do
                    end do
                end do
            end do
           !$acc end kernels
        end do
        !******************************
        !******************************
        !*** Plane wave deformation ***
        !******************************
        !******************************
        if(nbmn.EQ.1)then
            !$acc parallel loop &
            !$acc collapse(4) &
            !$acc gang vector &
            !$acc independent
            do isite=1,nsite
                do ispin=1,nspin
                    do imat=1,nmat
                        do jmat=1,nmat
                            !$acc loop seq
                            do jspin=1,nspin
                                pf2(imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)&
                                    +Gam123(ispin,jspin)*pf1(imat,jmat,jspin,isite)
                            end do
                        end do
                    end do
                end do
            end do
          !$acc end parallel loop
        end if
        return

    END SUBROUTINE Multiply_Dirac_device

    SUBROUTINE Multiply_Dirac_dagger_device_cuda(temperature,xmat,phase,Gam123,nbmn,&
        pf1,pf2,xptr_d,pf1ptr_d,pf2ptr_d)
        use cublasinterface
        use compiletimeconstants
        use gammamatrix
        implicit none

        integer, intent(in) :: nbmn
        double precision, intent(in) :: temperature
        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double complex, intent(in) :: phase(1:nmat,1:nmat,1:2)
        double complex, intent(in) :: Gam123(1:nspin,1:nspin)
        double complex, intent(in) :: pf1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex, intent(out) :: pf2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double precision :: lattice_spacing
        double complex :: gamtmp
        !******************
        integer :: imat,jmat,kmat
        integer :: idim
        integer :: ispin,jspin,kspin
        integer :: isite
        integer :: countr
        !$acc declare present(temperature,xmat,phase,Gam123,pf1,pf2)

        type(c_devptr), device :: xptr_d(nsite,ndim)
        type(c_devptr), device :: pf1ptr_d(nsite,nspin)
        type(c_devptr), device :: pf2ptr_d(nsite,nspin)


        !**********************
        !**********************
        !***  kinetic part  ***
        !**********************
        !**********************
        !$acc kernels
        do isite=1,nsite
            do ispin=1,nspin
                do jmat=1,nmat
                    do imat=1,nmat
                        pf2(imat,jmat,ispin,isite)=(0d0,0d0)
                    end do
                end do
            end do
        end do
        !$acc end kernels
        !********************
        !*** Naive Action ***
        !********************
        if(nimprove.EQ.0)then
            !$acc parallel loop &
            !$acc collapse(4) &
            !$acc gang vector &
            !$acc independent
            do isite=1,nsite
                do ispin=1,8
                    do jmat=1,nmat
                        do imat=1,nmat
                            pf2(imat,jmat,ispin+8,isite)=&
                                pf2(imat,jmat,ispin+8,isite)&
                                -(1d0,0d0)*phase(imat,jmat,1)*pf1(imat,jmat,ispin,isite+1)&
                                +(1d0,0d0)*pf1(imat,jmat,ispin,isite)
                            pf2(imat,jmat,ispin,isite)=&
                                pf2(imat,jmat,ispin,isite)&
                                +(1d0,0d0)*dconjg(phase(imat,jmat,1))*pf1(imat,jmat,ispin+8,isite-1)&
                                -(1d0,0d0)*pf1(imat,jmat,ispin+8,isite)
                        end do
                    end do
                end do
            end do
           !$acc end parallel loop
           !***********************
           !*** Improved Action ***
           !***********************
        else if(nimprove.EQ.1)then
            !$acc parallel loop &
            !$acc collapse(4) &
            !$acc gang vector &
            !$acc independent
            do isite=1,nsite
                do ispin=1,8
                    do jmat=1,nmat
                        do imat=1,nmat
                            pf2(imat,jmat,ispin+8,isite)=&
                                pf2(imat,jmat,ispin+8,isite)&
                                +(0.5d0,0d0)*phase(imat,jmat,2)*pf1(imat,jmat,ispin,isite+2)&
                                -(2d0,0d0)*phase(imat,jmat,1)*pf1(imat,jmat,ispin,isite+1)&
                                +(1.5d0,0d0)*pf1(imat,jmat,ispin,isite)

                            pf2(imat,jmat,ispin,isite)=&
                                pf2(imat,jmat,ispin,isite)&
                                -(0.5d0,0d0)*dconjg(phase(imat,jmat,2))*pf1(imat,jmat,ispin+8,isite-2)&
                                +(2d0,0d0)*dconjg(phase(imat,jmat,1))*pf1(imat,jmat,ispin+8,isite-1)&
                                -(1.5d0,0d0)*pf1(imat,jmat,ispin+8,isite)
                        end do
                    end do
                end do
            end do
           !$acc end parallel loop
        end if
        !**************************
        !**************************
        !***  interaction part  ***
        !**************************
        !**************************
        lattice_spacing=1d0/temperature/dble(nsite)
        if(cublasstream.EQ.1) then
            call multiply_cublas_pointer_streams(xptr_d,pf1ptr_d,pf2ptr_d,dcmplx(lattice_spacing),.TRUE.)
        else
            do countr=1,144
                ispin=gamispin(countr)
                jspin=gamjspin(countr)
                idim=gamidim(countr)
                gamtmp=dconjg(gamgam(countr))
                call multiply_cublas_pointer(xptr_d,pf1ptr_d,pf2ptr_d,dcmplx(lattice_spacing)*gamtmp,dcmplx(1.d0),ispin,jspin,idim)
                !call callZgemmBatched(pf2,xmat,pf1,dcmplx(lattice_spacing)*gamtmp,dcmplx(1.d0),ispin,jspin,idim)
                !call callZgemmBatched(pf2,xmat,pf1,dcmplx(lattice_spacing)*gamtmp,dcmplx(1.d0),ispin,jspin,idim)
            end do
        end if
        !******************************
        !******************************
        !*** Plane wave deformation ***
        !******************************
        !******************************
        if(nbmn.EQ.1)then
            !$acc parallel loop &
            !$acc collapse(4) &
            !$acc gang vector &
            !$acc independent
            do isite=1,nsite
                do ispin=1,nspin
                    do imat=1,nmat
                        do jmat=1,nmat
                            !$acc loop seq
                            do jspin=1,nspin
                                pf2(imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)&
                                    +dconjg(Gam123(jspin,ispin))*pf1(imat,jmat,jspin,isite)
                            end do
                        end do
                    end do
                end do
            end do
          !$acc end parallel loop
        end if
        return

    END SUBROUTINE Multiply_Dirac_dagger_device_cuda

    SUBROUTINE Multiply_Dirac_device_cuda(temperature,xmat,phase,Gam123,nbmn,&
        pf1,pf2,xptr_d,pf1ptr_d,pf2ptr_d)
        use cublasinterface
        use compiletimeconstants
        use gammamatrix
        implicit none
        integer, intent(in) :: nbmn
        double precision, intent(in) :: temperature
        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double complex, intent(in) :: phase(1:nmat,1:nmat,1:2)
        double complex, intent(in) :: Gam123(1:nspin,1:nspin)
        double complex, intent(in) :: pf1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex, intent(out) :: pf2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double precision :: lattice_spacing
        double complex :: gamtmp
        !******************
        integer :: imat,jmat,kmat
        integer :: idim
        integer :: ispin,jspin,kspin
        integer :: isite
        integer :: countr

        !$acc declare present(temperature,xmat,phase,Gam123,pf1,pf2)

        type(c_devptr), device :: xptr_d(nsite,ndim)
        type(c_devptr), device :: pf1ptr_d(nsite,nspin)
        type(c_devptr), device :: pf2ptr_d(nsite,nspin)



        !**********************
        !**********************
        !***  kinetic part  ***
        !**********************
        !**********************
        !$acc parallel loop &
        !$acc independent collapse(4) gang vector
        do isite=1,nsite
            do ispin=1,nspin
                do jmat=1,nmat
                    do imat=1,nmat
                        pf2(imat,jmat,ispin,isite)=(0d0,0d0)
                    end do
                end do
            end do
        end do
        !$acc end parallel loop
        !********************
        !*** Naive Action ***
        !********************
        if(nimprove.EQ.0)then
            !$acc parallel loop &
            !$acc collapse(4) &
            !$acc gang vector &
            !$acc independent
            do isite=1,nsite
                do ispin=1,8
                    do jmat=1,nmat
                        do imat=1,nmat
                            pf2(imat,jmat,ispin+8,isite)=&
                                pf2(imat,jmat,ispin+8,isite)&
                                +(1d0,0d0)*phase(imat,jmat,1)*pf1(imat,jmat,ispin,isite+1)&
                                -(1d0,0d0)*pf1(imat,jmat,ispin,isite)
                            pf2(imat,jmat,ispin,isite)=&
                                pf2(imat,jmat,ispin,isite)&
                                -(1d0,0d0)*dconjg(phase(imat,jmat,1))*pf1(imat,jmat,ispin+8,isite-1)&
                                +(1d0,0d0)*pf1(imat,jmat,ispin+8,isite)
                        end do
                    end do
                end do
            end do
           !$acc end parallel loop
           !***********************
           !*** Improved Action ***
           !***********************
        else if(nimprove.EQ.1)then
            !$acc parallel loop &
            !$acc collapse(4) &
            !$acc gang vector &
            !$acc independent
            do isite=1,nsite
                do ispin=1,8
                    do jmat=1,nmat
                        do imat=1,nmat
                            pf2(imat,jmat,ispin+8,isite)=&
                                pf2(imat,jmat,ispin+8,isite)&
                                -(0.5d0,0d0)*phase(imat,jmat,2)*pf1(imat,jmat,ispin,isite+2)&
                                +(2d0,0d0)*phase(imat,jmat,1)*pf1(imat,jmat,ispin,isite+1)&
                                -(1.5d0,0d0)*pf1(imat,jmat,ispin,isite)

                            pf2(imat,jmat,ispin,isite)=&
                                pf2(imat,jmat,ispin,isite)&
                                +(0.5d0,0d0)*dconjg(phase(imat,jmat,2))*pf1(imat,jmat,ispin+8,isite-2)&
                                -(2d0,0d0)*dconjg(phase(imat,jmat,1))*pf1(imat,jmat,ispin+8,isite-1)&
                                +(1.5d0,0d0)*pf1(imat,jmat,ispin+8,isite)
                        end do
                    end do
                end do
            end do
           !$acc end parallel loop
        end if
        !**************************
        !**************************
        !***  interaction part  ***
        !**************************
        !**************************
        lattice_spacing=1d0/temperature/dble(nsite)
        if(cublasstream.EQ.1) then
            call multiply_cublas_pointer_streams(xptr_d,pf1ptr_d,pf2ptr_d,dcmplx(lattice_spacing),.FALSE.)
        else
            do countr=1,144
                ispin=gamispin(countr)
                jspin=gamjspin(countr)
                idim=gamidim(countr)
                gamtmp=gamgam(countr)
                call multiply_cublas_pointer(xptr_d,pf1ptr_d,pf2ptr_d,dcmplx(lattice_spacing)*gamtmp,dcmplx(1.d0),ispin,jspin,idim)
                !call callZgemmBatched(pf2,xmat,pf1,dcmplx(lattice_spacing)*gamtmp,dcmplx(1.d0),ispin,jspin,idim)
            end do
        end if
        !******************************
        !******************************
        !*** Plane wave deformation ***
        !******************************
        !******************************
        if(nbmn.EQ.1)then
            !$acc parallel loop &
            !$acc collapse(4) &
            !$acc gang vector &
            !$acc independent
            do isite=1,nsite
                do ispin=1,nspin
                    do imat=1,nmat
                        do jmat=1,nmat
                            !$acc loop seq
                            do jspin=1,nspin
                                pf2(imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)&
                                    +Gam123(ispin,jspin)*pf1(imat,jmat,jspin,isite)
                            end do
                        end do
                    end do
                end do
            end do
          !$acc end parallel loop
        end if
        return

    END SUBROUTINE Multiply_Dirac_device_cuda

end module dirac_operator

