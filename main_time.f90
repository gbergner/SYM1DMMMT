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
    implicit none

    integer, intent(in) :: var
    double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    double complex, intent(in) :: pf1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex, intent(out) :: pf2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    double complex :: gamtmp,tmp
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
        end if
    END SUBROUTINE Test_Multiply_Dirac_dagger_device


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

    SUBROUTINE Multiply_Dirac_device_memred(temperature,xmat,phase,Gam123,nbmn,&
        pf1,pf2)
        use cublasinterface
        use compiletimeconstants
        use gammamatrix
        implicit none

        integer, intent(in) :: nbmn
        double precision, intent(in) :: temperature
        double complex, intent(in) :: xmat(1:(nmat*(nmat+1)/2),1:ndim,-(nmargin-1):nsite+nmargin)
        double complex, intent(in) :: phase(1:nmat,1:nmat,1:2)
        double complex, intent(in) :: Gam123(1:nspin,1:nspin)
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

        !$acc declare present(temperature,xmat,nbmn,phase,Gam123,pf1,pf2)


        !**********************
        !**********************
        !***  kinetic part  ***
        !**********************
        !**********************
        !$acc parallel loop &
        !$acc independent collapse(4) gang, vector
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
        !return
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
            gamtmp=gamgam(countr)
            if(.FALSE.) then
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
            if(.TRUE.) then
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
            if(.FALSE.) then
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

    END SUBROUTINE Multiply_Dirac_device_memred


    program timing_mult
        use compiletimeconstants
        use dirac_operator
        use outputstreamnumbers
        use mtmod !Mersenne twistor
        implicit none
  

        integer nbc !boundary condition for fermions; 0 -> pbc, 1 -> apbc
        integer nbmn ! 0 -> BFSS, 1 -> BMN
        integer init !initial condition; 0 -> continue, 1 -> new
        integer isave !0 -> save intermediate config, 1 -> do not save
        integer nsave!saved every nsave trajectories
        !matrices
        double complex :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        ! double complex, dimension(:,:,:,:), allocatable :: xmat
        !gauge field
        double precision alpha(1:nmat)
        !remez coefficients
        double precision acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)!molecular evolution
        double precision acoeff_pf(0:nremez_pf),bcoeff_pf(1:nremez_pf)!pseudo fermion
        double precision upper_approx!the largest eigenvalue of (M^¥dagger M)must be smaller than this.

        !CG solver
        integer max_iteration, iteration,n_bad_CG
        double precision max_err
        !For Mersenne Twister
        integer mersenne_seed
        !Fourier acceleration
        double precision acceleration(1:nsite)
        double precision fluctuation(1:nsite)
        integer iaccelerate,imeasure
        integer imetropolis

        !number of CG iteration for calculating the largest and smallest eigenvalues of D=(M^dag*M)
        integer neig_max,neig_min
        !number of fuzzy sphere, when init=2.
        integer nfuzzy
        !Gamma matrices
        double complex Gamma10d(1:ndim,1:nspin,1:nspin)

        integer iremez
        integer ncv!ncv=number of constraint[max(alpha_i-alpha_j)<2*pi] violation

        double precision temperature
        double precision flux
        !parameters for molecular evolution
        integer ntau,nratio
        double precision dtau_xmat,dtau_alpha

        !number of trajectories
        integer ntraj     !total number of trajectories at the end of the run
        integer itraj

        double precision ham_init,ham_fin

        !measurements
        integer nskip !measurement is performed every nskiptrajectories
        integer nacceptance !number of acceptance
        integer ntrial !number of trial
        character(1000) input_config,data_output,output_config,acc_input,acc_output,intermediate_config,CG_log


        ! For MPI
        !integer IERR,NPROCS,MYRANK




        !coefficient for potential for alpha
        double precision g_alpha
        !coefficient for potential for Tr(x^2)
        double precision g_R,RCUT

        !    include 'include.h'

  
        double complex :: phase(1:nmat,1:nmat,1:2)
        double complex :: Gam123(1:nspin,1:nspin)
        !$acc declare device_resident(phase,Gam123)
  
        double complex testvect1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex testvect2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex testvect0(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
  
        double complex :: testvect_d0(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex :: testvect_d2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        !$acc declare device_resident(testvect_d0,testvect_d2)

        double complex :: xmatreduced(1:(nmat*(nmat+1)/2),1:ndim,-(nmargin-1):nsite+nmargin)
        double complex :: xmatexpand(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)

        real :: start_time,stop_time,time1,time2
        integer :: counter,isite,idim,imat,jmat,ivar
        double complex :: tmp
  
        if(npf.EQ.1)then
     include 'remez_md_1.dat'
     include 'remez_pf_1.dat'
        else if(npf.EQ.2)then
     include 'remez_md_2.dat'
     include 'remez_pf_2.dat'
        end if
     
         !allocate(xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin))
     
        !---------------------------------
  
        !*************************
        !**** read parameters ****
        !*************************
        open(unit=unit_input_para,status='OLD',file='input_5.dat',&
            &action='READ')
        read(unit_input_para,*) input_config
        read(unit_input_para,*) output_config
        read(unit_input_para,*) data_output
        read(unit_input_para,*) intermediate_config
        read(unit_input_para,*) acc_input
        read(unit_input_para,*) acc_output
        read(unit_input_para,*) CG_log
        read(unit_input_para,*) nbc
        read(unit_input_para,*) nbmn
        read(unit_input_para,*) init
        read(unit_input_para,*) iaccelerate
        read(unit_input_para,*) isave
        read(unit_input_para,*) temperature
        read(unit_input_para,*) flux
        read(unit_input_para,*) ntraj
        read(unit_input_para,*) nskip
        read(unit_input_para,*) nsave
        read(unit_input_para,*) ntau
        read(unit_input_para,*) nratio
        read(unit_input_para,*) dtau_xmat
        read(unit_input_para,*) dtau_alpha
        read(unit_input_para,*) upper_approx
        read(unit_input_para,*) max_err
        read(unit_input_para,*) max_iteration
        read(unit_input_para,*) g_alpha
        read(unit_input_para,*) g_R
        read(unit_input_para,*) RCUT
        read(unit_input_para,*) neig_max
        read(unit_input_para,*) neig_min
        read(unit_input_para,*) nfuzzy
        read(unit_input_para,*) mersenne_seed
        read(unit_input_para,*) imetropolis
        close(unit_input_para)
        !Construc Gamma matrices.
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
        !$acc copyin(testvect0,temperature,xmat,alpha,GAMMA10d,nbmn,flux,nbc,xmatreduced,xmatexpand)
        !$acc kernels
        testvect_d0=testvect0
        !$acc end kernels
        call  setup_data_device(alpha,flux,GAMMA10d,phase,Gam123,temperature)
        call cpu_time(start_time)
        do counter=1,20
            call Multiply_Dirac_dagger(temperature,xmat,alpha,&
                testvect0,testvect1,GAMMA10d,nbmn,flux)
        end do
        call cpu_time(stop_time)
        time1=stop_time - start_time
        print *, "Host time:", &
            time1, "seconds "

        call cpu_time(start_time)
        do counter=1,20
            call Multiply_Dirac_dagger_device(temperature,xmat,phase,Gam123,nbmn,&
                testvect_d0,testvect_d2)
        end do
        call cpu_time(stop_time)
        time2=stop_time - start_time
        print *, "Device time:", &
            time2, "seconds "
        print *,"relative timing",time1/time2
        !$acc kernels
        testvect2=testvect_d2
        tmp=Sum(abs(testvect_d2))
        !$acc end kernels
        print *,"vect ",tmp," ", Sum(abs(testvect1)), " err ",Sum(abs(testvect2-testvect1))
        call cpu_time(start_time)
        do counter=1,20
            call Multiply_Dirac(temperature,xmat,alpha,&
                testvect0,testvect1,GAMMA10d,nbmn,flux)
        end do
        call cpu_time(stop_time)
        time1=stop_time - start_time
        print *, "Host time:", &
            time1, "seconds"

        call cpu_time(start_time)
        do counter=1,20
            call Multiply_Dirac_device(temperature,xmat,phase,Gam123,nbmn,&
                testvect_d0,testvect_d2)
        end do
        call cpu_time(stop_time)
        time2=stop_time - start_time
        print *, "Device time:", &
            time2, "seconds"
        print *,"relative timing",time1/time2
        !$acc kernels
        testvect2=testvect_d2
        tmp=Sum(abs(testvect_d2))
        !$acc end kernels
        print *,"vect ",tmp," ", Sum(abs(testvect1)), " err ",Sum(abs(testvect2-testvect1))
  
        print *, "Memory reduced variant "
        call reduce_mem_xmat(xmat,xmatreduced)
        call expand_mem_xmat(xmatreduced,xmatexpand)
        !$acc update host(xmatexpand)
        print *,"xmatred err ",Sum(abs(xmatexpand))," ", Sum(abs(xmat))," ",Sum(abs(xmatexpand-xmat))

        call cpu_time(start_time)
        do counter=1,20
            call Multiply_Dirac_device_memred(temperature,xmatreduced,phase,Gam123,nbmn,&
                testvect_d0,testvect_d2)
        end do
        call cpu_time(stop_time)
        time2=stop_time - start_time
        print *, "Device time:", &
            time2, "seconds"
        print *,"relative timing",time1/time2
        !$acc kernels
        testvect2=testvect_d2
        tmp=Sum(abs(testvect_d2))
        !$acc end kernels
        print *,"vect ",tmp," ", Sum(abs(testvect1)), " err ",Sum(abs(testvect2-testvect1))


        do ivar=1,4
            print *, "Variant ", ivar
            call cpu_time(start_time)
            do counter=1,20
                call Test_Multiply_Dirac_dagger_device(xmat,testvect_d0,testvect_d2,ivar)
            end do
            call cpu_time(stop_time)
            time2=stop_time - start_time
            print *, "Device time:", &
                time2, "seconds"
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
        end do
        !$acc end data
      !End test part
      !deallocate(xmat)
    end program timing_mult
