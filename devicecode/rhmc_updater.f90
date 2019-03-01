module RHMC_Updater
    implicit none
contains
    SUBROUTINE check_ham_host(temperature,xmat,alpha,P_xmat,P_alpha,&
        &pf,chi,acoeff_md,g_R,RCUT,nbmn,flux,nbc,nremez_md,gamma10d,bcoeff_md,max_err,max_iteration,iteration,&
        &hamold,gam123,phase,ngauge,purebosonic)
        use compiletimeconstants
        use hmc_molecular_dynamics
        implicit none
        integer nbc,nbmn,nremez_md,ngauge,purebosonic
        double precision temperature,flux
        double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
        double precision max_err
        integer max_iteration,iteration
        double precision acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)
        double precision g_alpha,g_R,RCUT
        !input & output
        double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double precision alpha(1:nmat)
        double precision ham, hamold
        integer info
        double complex P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision P_alpha(1:nmat)
        double complex pf(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
        !$acc declare device_resident(P_xmat,P_alpha,pf)
        double complex P_xmat_h(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision P_alpha_h(1:nmat)
        double complex pf_h(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
        double complex Chi(1:nremez_md,1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
        integer IERR,myrank,nprocs

        double complex, intent(in) :: phase(1:nmat,1:nmat,1:2)
        double complex, intent(in) :: Gam123(1:nspin,1:nspin)
        double complex :: testvect_d0(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex :: testvect_d1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex :: testvect0(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex :: testvect1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        !$acc declare device_resident(phase,Gam123,testvect_d0,testvect_d1)
        !$acc update host(alpha)
        !$acc update host(xmat)
        !$acc kernels
        P_xmat_h=P_xmat
        P_alpha_h=P_alpha
        pf_h=pf
        !$acc end kernels
        if(purebosonic.eq.0) then
            call solver_biCGm(nbc,nbmn,nremez_md,&
                xmat,alpha,pf_h,chi,GAMMA10d,&
                bcoeff_md,max_err,max_iteration,iteration,&
                temperature,flux,info)
        end if
        call Calc_Ham(temperature,xmat,alpha,P_xmat_h,P_alpha_h,ham,pf_h,chi,&
            &acoeff_md,g_R,RCUT,nbmn,flux,ngauge,purebosonic)
        print*,"host test ham ",ham
        if (abs(ham-hamold).GE.(0.010d0)) then
            print*, "WARNING: large miss device to host in Hamiltonian: ",ham, " ",hamold
            if (abs(ham-hamold).GE.(0.010d0)) then
                print*, "STOP: ERROR: large miss device host Hamiltonian."
                call exit(-1)
            end if
        end if
    end subroutine

    SUBROUTINE hamilton_calculation(temperature,xmat,alpha,P_xmat,P_alpha,ham,pf,&
        &acoeff_md,g_R,RCUT,nbmn,flux,phase,bcoeff_md,info_CG,max_err,max_iteration,&
        &iteration,gamma10d,gam123,nbc,ngauge,purebosonic)
        use compiletimeconstants
        use dirac_operator
        use cgm_solver
        use outputstreamnumbers
        use hmc_hamiltonean
        use hmc_molecular_dynamics
        use timer
        implicit none
        !***** input *****
        integer, intent(in) :: nbmn,nbc,ngauge,purebosonic
        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double precision, intent(in) :: alpha(1:nmat)
        double complex, intent(in) :: P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision, intent(in) :: P_alpha(1:nmat)
        double complex, intent(in) :: pf(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
        double precision, intent(in) :: acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)
        double complex, intent(inout) :: phase(1:nmat,1:nmat,1:2)
        double precision, intent(in) :: temperature,flux
        double precision, intent(in) :: g_R,RCUT
        double complex :: Chi(1:nmat,1:nmat,1:nspin,&
            &-(nmargin-1):nsite+nmargin,1:nremez_md,1:npf)
        double complex, intent(in):: GAMMA10d(1:ndim,1:nspin,1:nspin)
        double complex, intent(in) :: Gam123(1:nspin,1:nspin)
        !$acc declare present(nbmn,xmat,alpha,pf,temperature,flux,g_R,RCUT,acoeff_md,bcoeff_md,GAMMA10d)
        !$acc declare device_resident(Chi,P_xmat,P_alpha,phase,gam123)
        !***** output *****
        double precision, intent(out) :: ham
        double precision max_err
        integer max_iteration,iteration
        integer info_CG
        double complex tmp
        call print_time_step("hamilton calculation start")
        call update_data_device(alpha,phase)
        call print_time_step("hamilton cgm inverter start")
        if(purebosonic.eq.0) then
            call cgm_solver_device(nremez_md,bcoeff_md,nbmn,nbc,temperature,&
                max_err,max_iteration,xmat,phase,Gam123,pf,chi,info_CG,iteration)
        end if
        call print_time_step("hamilton cgm inverter end")
        !info_CG_init=0 -> OK (CG solver converged)
        !info_CG_init=1 -> error (CG solver did not converge)

        !Take CG_log
        !write(unit_CG_log,*)"ham_init",iteration
        ! requires summation and distribution
        call Calc_Ham_device(temperature,xmat,alpha,P_xmat,P_alpha,ham,pf,chi,&
            &acoeff_md,g_R,RCUT,nbmn,flux,phase,ngauge,purebosonic)
        if(rhmc_verbose.EQ.1) then
            !$acc kernels
            tmp=Sum(chi)
            !$acc end kernels
            print *,"check Ham ",ham," inverter result ",tmp
        end if
        if(check_host_metropolis.EQ.1) then
            call  check_ham_host(temperature,xmat,alpha,P_xmat,P_alpha,&
                &pf,chi,acoeff_md,g_R,RCUT,nbmn,flux,nbc,nremez_md,gamma10d,bcoeff_md,max_err,&
                &max_iteration,iteration,ham,gam123,phase,ngauge,purebosonic)
        end if
        call print_time_step("hamilton calculation end")
    end SUBROUTINE hamilton_calculation

    SUBROUTINE RHMC_evolution_device(xmat,alpha,phase,ncv,n_bad_CG,nacceptance,nbc,nbmn,&
        &temperature,flux,GAMMA10d,Gam123,ntau,nratio,dtau_xmat,dtau_alpha,&
        &acceleration,g_alpha,g_R,RCUT,&
        &acoeff_md,bcoeff_md,acoeff_pf,bcoeff_pf,max_err,max_iteration,iteration,&
        &ham_init,ham_fin,ntrial,imetropolis,ngauge,purebosonic)
        use compiletimeconstants
        use outputstreamnumbers
        use hmc_molecular_dynamics
        use rand_fermion_fields
        use lattice_action
        use mtmod !Mersenne twistor
        use timer
        use  utils_measurements
        implicit none

        !input
        integer, intent(in):: nbc,nbmn,ngauge,purebosonic
        double precision, intent(in):: temperature,flux
        double complex, intent(in):: GAMMA10d(1:ndim,1:nspin,1:nspin)
        !$acc declare present(nbc,nbmn,temperature,flux,GAMMA10d)
        double complex, intent(inout) :: phase(1:nmat,1:nmat,1:2)
        double complex, intent(in) :: Gam123(1:nspin,1:nspin)
        !$acc declare device_resident(phase,Gam123)
        double precision, intent(in) :: max_err
        integer, intent(in) :: max_iteration
        double precision, intent(in) :: acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)
        double precision, intent(in) :: acoeff_pf(0:nremez_pf),bcoeff_pf(1:nremez_pf)
        double precision, intent(in) :: g_alpha,g_R,RCUT
        !$acc declare present(acoeff_md,acoeff_pf,bcoeff_md,bcoeff_pf,g_alpha,g_R,RCUT)
        integer, intent(in) :: ntau,nratio
        integer, intent(in) :: imetropolis!1-> no Metropolis test
        double precision, intent(in) :: acceleration(1:nsite),dtau_alpha,dtau_xmat
        !$acc declare present(acceleration)
        !input & output
        double complex, intent(inout) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double precision, intent(inout) :: alpha(1:nmat)
        !$acc declare present(xmat,alpha)
        integer ncv,n_bad_CG,nacceptance,ntrial
        double precision ham_init,ham_fin
        !output
        integer iteration


        double precision metropolis
        integer info_pf,info_mol,info_CG_init,info_CG_fin,info_alpha,info_accept,info
        double complex backup_xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double precision backup_alpha(1:nmat)
        !$acc declare device_resident(backup_xmat,backup_alpha)
        double complex P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision P_alpha(1:nmat)
        double complex pf(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
        !double complex Chi(1:nmat,1:nmat,&
        !      1:nspin,-(nmargin-1):nsite+nmargin,1:nremez_md,1:npf)
        !!$acc declare device_resident(P_xmat,P_alpha,pf,Chi)
        !$acc declare device_resident(P_xmat,P_alpha,pf)
        integer IERR,myrank,nprocs
        double complex tmp
        !initialize info flugs
        info_pf=1
        info_mol=1
        info_CG_init=1
        info_CG_fin=1
        if(purebosonic.eq.1) then
            info_pf=0
            info_mol=0
            info_CG_init=0
            info_CG_fin=0
        end if
        !**********************************
        !**** Generate pseudo fermion. ****
        !**********************************
        call print_time_step("pf generation start")
        if(purebosonic.eq.0) then
            call generate_pseudo_fermion_SUN_device(pf,xmat,alpha,phase,&
                &Gam123,acoeff_pf,bcoeff_pf,max_err,max_iteration,iteration,&
                &nbc,nbmn,temperature,flux,info_pf)
        else
            pf=0d0
        end if
        call print_time_step("pf generation end")
        !info_pf=0 -> OK (CG solver converged)
        !info_pf=1 -> error (CG solver did not converge)
  
        !Take CG_log -> not device
        write(unit_CG_log,*)"pseudo-fermion generation",iteration
        !******************************
        !**** Molecular Evolution. ****
        !******************************
        !save the old config.
        ! on device
        !$acc kernels
        backup_xmat=xmat
        backup_alpha=alpha
        !$acc end kernels
        !generate auxiliary momenta.
        !all input on device
        call print_time_step("momentum generation start")
        call Generate_Momenta_device(P_xmat,P_alpha)
        call print_time_step("momentum generation end")
        if(rhmc_verbose.EQ.1) then
            !$acc kernels
            tmp=Sum(pf)
            !$acc end kernels
            print*, "check  device pseudofermion initial ", tmp
            !$acc kernels
            tmp=Sum(xmat)
            !$acc end kernels
            print*, "check device xmat initial ", tmp
            !$acc kernels
            tmp=Sum(P_xmat)
            !$acc end kernels
            print*, "check device momenta xmat initial ", tmp
            !$acc kernels
            tmp=Maxval(abs(P_alpha))
            !$acc end kernels
            print*, "check device momenta alpha initial ", tmp
        end if
        !Calculate ham_init
        !call update_data_device(alpha,phase)
        !call cgm_solver_device(nremez_md,bcoeff_md,nbmn,nbc,temperature,&
        !  max_err,max_iteration,xmat,phase,Gam123,pf,chi,info_CG_init,iteration)
        !info_CG_init=0 -> OK (CG solver converged)
        !info_CG_init=1 -> error (CG solver did not converge)
  
        !Take CG_log
        write(unit_CG_log,*)"ham_init",iteration
        ! requires summation and distribution
        call hamilton_calculation(temperature,xmat,alpha,P_xmat,P_alpha,ham_init,pf,&
            &acoeff_md,g_R,RCUT,nbmn,flux,phase,bcoeff_md,info_CG_init,max_err,max_iteration,&
            &iteration,gamma10d,gam123,nbc,ngauge,purebosonic)
        !call Calc_Ham_device(temperature,xmat,alpha,P_xmat,P_alpha,ham_init,pf,chi,&
         !    &acoeff_md,g_R,RCUT,nbmn,flux,phase)
        if(rhmc_verbose.EQ.1) then
            print *,"check initial Ham ",ham_init
        end if
        ! if(check_host_metropolis.EQ.1) then
        !   call  check_ham_host(temperature,xmat,alpha,P_xmat,P_alpha,&
        !        &pf,chi,acoeff_md,g_R,RCUT,nbmn,flux,nbc,nremez_md,gamma10d,bcoeff_md,max_err,max_iteration,iteration,ham_init,gam123,phase)
        ! end if
        !Molecular Evolution
        if(imetropolis.EQ.2) then
            call Molecular_Dynamics_device_SW(nbc,temperature,&
                &ntau,nratio,dtau_xmat,dtau_alpha,xmat,alpha,phase,P_xmat,P_alpha,&
                &acoeff_md,bcoeff_md,pf,max_iteration,max_err,iteration,&
                &gamma10d,Gam123,g_alpha,g_R,RCUT,acceleration,nbmn,flux,info_mol,ngauge,purebosonic)

        else
            call Molecular_Dynamics_device(nbc,temperature,&
                &ntau,nratio,dtau_xmat,dtau_alpha,xmat,alpha,phase,P_xmat,P_alpha,&
                &acoeff_md,bcoeff_md,pf,max_iteration,max_err,iteration,&
                &gamma10d,Gam123,g_alpha,g_R,RCUT,acceleration,nbmn,flux,info_mol,ngauge,purebosonic)
        end if
        !info_mol=0 -> OK (CG solver converged)
        !info_mol=1 -> error (CG solver did not converge)
        if(info_mol.EQ.0)then
            call hamilton_calculation(temperature,xmat,alpha,P_xmat,P_alpha,ham_fin,pf,&
                &acoeff_md,g_R,RCUT,nbmn,flux,phase,bcoeff_md,info_CG_fin,max_err,&
                &max_iteration,iteration,gamma10d,gam123,nbc,ngauge,purebosonic)
            !calculate ham_fin
            !call update_data_device(alpha,phase)
            !call cgm_solver_device(nremez_md,bcoeff_md,nbmn,nbc,temperature,&
            !      max_err,max_iteration,xmat,phase,Gam123,pf,chi,info_CG_fin,iteration)
            !info_CG_fin=0 -> OK (CG solver converged)
            !info_CG_fin=1 -> error (CG solver did not converge)
     
            !Take CG_log
     
           !call Calc_Ham_device(temperature,xmat,alpha,P_xmat,P_alpha,ham_fin,&
           !     &pf,chi,acoeff_md,g_R,RCUT,nbmn,flux,phase)
     
        end if
        if(rhmc_verbose.EQ.1) then
            !$acc kernels
            tmp=Sum(pf)
            !$acc end kernels
            print*, "check  device pseudofermion final ", tmp
            !$acc kernels
            tmp=Sum(xmat)
            !$acc end kernels
            print*, "check device xmat final ", tmp
            !$acc kernels
            tmp=Sum(P_xmat)
            !$acc end kernels
            print*, "check device momenta xmat final ", tmp
            !$acc kernels
            tmp=Maxval(abs(P_alpha))
            !$acc end kernels
            print*, "check device momenta alpha final ", tmp
        end if
         !if(check_host_metropolis.EQ.1) then
         !     call  check_ham_host(temperature,xmat,alpha,P_xmat,P_alpha,&
          !      &pf,chi,acoeff_md,g_R,RCUT,nbmn,flux,nbc,nremez_md,gamma10d,bcoeff_md,max_err,max_iteration,iteration,ham_fin,gam123,phase)
         !end if
        !#######################################
        !#### Did CG solver work correctly? ####
        !#######################################
        if((info_pf.EQ.0).AND.(info_mol.EQ.0).AND.(info_CG_init.EQ.0).AND.(info_CG_fin.EQ.0))then
            info=0!CG solver worked fine everywhere.
        else
            info=1!CG solver did not converge at least one place.
        end if
        !Count how many times the CG solver failed to converge.
        if(info.EQ.1)then
            n_bad_CG=n_bad_CG+1
        end if
        !#########################
        !#### metropolis test ####
        !#########################
        ntrial=ntrial+1
        metropolis=grnd()
        !"metropolis=grnd()" must be put outside "if(myrank.EQ.0)then"
        !so that random numbers synchronize.
        if(imetropolis.EQ.1)then
            !No Metropolis test
            metropolis=-1d0 !Then, it is always accepted.
            print*, "no metropolis test"
        end if

        !####################################
        !#### check constraint for alpha ####
        !####################################
        call check_alpha_constraint_device(alpha,info_alpha)
        !Count how many times the constraint is violated.
        if(info_alpha.EQ.1)then
            ncv=ncv+1
            info=1
            print*, "bad constraint"
        end if
        if(info.EQ.0)then
            !alpha satisfies the constraint.
            !CG converged everywhere during Molecular evolution.
        
            if(dexp(ham_init-ham_fin) > metropolis)THEN
                !accept
                info_accept=0
            else
                !reject
                info_accept=1
            end if
        else!automatic reject, before the Metropolis test.
            info_accept=1
            print*, "auto reject"
        end if

        if(info_accept.EQ.0)then
            !accept
            nacceptance=nacceptance+1
            ! Project out the trace parts of X and alpha.
            ! For X, we remove only the zero mode of X,
            !        i.e. (¥int dt X) is made traceless.
            call subtract_U1_device(xmat,alpha)
        else 
            !reject
            !$acc kernels
            xmat=backup_xmat
            alpha=backup_alpha
           !$acc end kernels
        end if
        if(rhmc_verbose.EQ.1) then
            print*,"RHMC ev finished accepted ",info_accept," ",info_alpha," ",metropolis, " ", info," ", dexp(ham_init-ham_fin),"\n"
        end if
        write(unit_CG_log,*)"ham_fin",iteration, " ", ham_init-ham_fin," ",ham_init," ",ham_fin," ",info_accept," ",info
        return

    END SUBROUTINE RHMC_evolution_device

end module RHMC_Updater
