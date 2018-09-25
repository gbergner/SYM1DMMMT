!##############################################################################
!######              BFSS matrix model on lattice                     #########
!######              The code of Masanori Hanada reorganized          #########
!######              and moved to the GPU by Georg Bergner.           #########
!##############################################################################

program BFSS_DEVICE
    use compiletimeconstants
    use dirac_operator
    use mtmod !Mersenne twistor
    use outputstreamnumbers
    use RHMC_Updater
    use timer
    use utils_measurements
    implicit none
  

    integer nbc !boundary condition for fermions; 0 -> pbc, 1 -> apbc
    integer nbmn ! 0 -> BFSS, 1 -> BMN
    integer init !initial condition; 0 -> continue, 1 -> new
    integer isave !0 -> save intermediate config, 1 -> do not save
    integer nsave!saved every nsave trajectories
    !matrices
    !double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    double complex, dimension(:,:,:,:),allocatable :: xmat
    !gauge field
    double precision alpha(1:nmat)
    !copyin
    !remez coefficients
    !double precision acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)!molecular evolution
    double precision,dimension(:),allocatable :: acoeff_md,bcoeff_md!molecular evolution
    !copyin
    !double precision acoeff_pf(0:nremez_pf),bcoeff_pf(1:nremez_pf)!pseudo fermion
    double precision,dimension(:),allocatable :: acoeff_pf,bcoeff_pf!pseudo fermion
    !copyin
    double precision upper_approx!the largest eigenvalue of (M^¥dagger M)must be smaller than this.
    !copyin

     !CG solver
    integer max_iteration, iteration,n_bad_CG
    double precision max_err
    !For Mersenne Twister
    integer mersenne_seed
    !Fourier acceleration
    double precision acceleration(1:nsite)
    !copyin
    double precision fluctuation(1:nsite)
    integer iaccelerate,imeasure
    integer imetropolis

    !number of CG iteration for calculating the largest and smallest eigenvalues of D=(M^dag*M)
    integer neig_max,neig_min
    !number of fuzzy sphere, when init=2.
    integer nfuzzy
    !Gamma matrices
    double complex Gamma10d(1:ndim,1:nspin,1:nspin)
    !copyin
    double complex :: phase(1:nmat,1:nmat,1:2)
    !$acc declare device_resident(phase)
    double complex :: Gam123(1:nspin,1:nspin)
    !$acc declare device_resident(Gam123)

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
    character(1000) input_config,data_output,output_config,acc_input,acc_output,intermediate_config,CG_log,checkpointn
    !character,dimension(:),allocatable :: input_config,data_output,output_config,acc_input,acc_output,intermediate_config,CG_log
    ! For MPI
    !integer IERR,NPROCS,MYRANK
    !coefficient for potential for alpha
    double precision g_alpha
    !coefficient for potential for Tr(x^2)
    double precision g_R,RCUT

    call print_timestamp("starting program")
    call start_timer()

    allocate(acoeff_pf(0:nremez_pf))
    allocate(bcoeff_pf(1:nremez_pf))
    allocate(acoeff_md(0:nremez_md))
    allocate(bcoeff_md(1:nremez_md))

    allocate(xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin))

    if(npf.EQ.1)then
     include 'remez_md_1.dat'
     include 'remez_pf_1.dat'
    else if(npf.EQ.2)then
     include 'remez_md_2.dat'
     include 'remez_pf_2.dat'
    end if
     



     checkpointn="chp.dat"
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
    if(ngauge.EQ.0)then
        dtau_alpha=0d0
    end if
    !***************************************
    !*** Rescaling of Remez coefficients ***
    !***************************************
    !acoeff_md(0)= acoeff_md(0)*upper_approx**(-0.25d0)
    !do iremez=1,nremez_md
    !   acoeff_md(iremez)=acoeff_md(iremez)*upper_approx**(0.75d0)
    !   bcoeff_md(iremez)=bcoeff_md(iremez)*upper_approx
    !end do
    !acoeff_pf(0)= acoeff_pf(0)*upper_approx**(0.125d0)
    !do iremez=1,nremez_pf
    !   acoeff_pf(iremez)=acoeff_pf(iremez)*upper_approx**(1.125d0)
    !   bcoeff_pf(iremez)=bcoeff_pf(iremez)*upper_approx
    !end do
  
    acoeff_md(0)= acoeff_md(0)*upper_approx**(-0.25d0/dble(npf))
    do iremez=1,nremez_md
        acoeff_md(iremez)=acoeff_md(iremez)*upper_approx**(1d0-0.25d0/dble(npf))
        bcoeff_md(iremez)=bcoeff_md(iremez)*upper_approx
    end do
    acoeff_pf(0)= acoeff_pf(0)*upper_approx**(0.125d0/dble(npf))
    do iremez=1,nremez_pf
        acoeff_pf(iremez)=acoeff_pf(iremez)*upper_approx**(1d0+0.125d0/dble(npf))
        bcoeff_pf(iremez)=bcoeff_pf(iremez)*upper_approx
    end do
     
    !*************************************
    !*** Set the initial configuration ***
    !*************************************
    call initial_configuration(xmat,alpha,acceleration,itraj,init,&
        &iaccelerate,nfuzzy,input_config,acc_input,flux,mersenne_seed)
     if(init.EQ.4) then
      call read_checkpoint(xmat,alpha,itraj)
     end if
    !***********************************
    !******  Make the output file ******
    !***********************************
    call output_header(data_output,temperature,flux,&
        &ntau,nratio,dtau_xmat,dtaU_alpha,neig_max,neig_min,nbc,nbmn,&
        &init,input_config,output_config,iaccelerate,acc_input,acc_output,&
        &g_alpha,g_R,RCUT,upper_approx,max_err,max_iteration,CG_log,&
        &isave,nsave,intermediate_config,imetropolis)

    !*******************************************************
    !******  Make the intermediate configuration file ******
    !*******************************************************
    if(isave.EQ.0)then
        open(unit=unit_intermediate_config,status='REPLACE',&
            &file=intermediate_config,action='WRITE')
    end if
    !**********************************
    !**** initialize counters *********
    !**********************************
    nacceptance=0 !number of acceptance
    ntrial=0 !number of trial
    ncv=0 !number of the constraint violation.
    n_bad_CG=0 !count how many times CG solver failed to converge.
  
    !************************************************************************
    !**** initialize the variables for optimizeing Fourier acceleration. ****
    !************************************************************************
    imeasure=0!Number of measurements of fluctuation.
    fluctuation=0d0 !fluctuation of scalar fields X.
    call print_time_step("finished init")
  
    !$acc data &
    !$acc copyin(temperature,xmat,alpha,GAMMA10d,nbmn,flux,nbc)&
    !$acc copyin(g_R,RCUT,acoeff_md,bcoeff_md,acoeff_pf,bcoeff_pf,g_alpha,acceleration)

    call print_time_step("finished copy to device")
    call  setup_data_device(alpha,flux,GAMMA10d,phase,Gam123,temperature)
    !*************************************
    !*************************************
    !***  Loop of molecular evolutions ***
    !*************************************
    !*************************************
    do while (itraj.LE.ntraj)
        !******************************
        !**** Molecular Evolution. ****
        !******************************
        !Take CG_log
        write(unit_CG_log,*)"itraj=",itraj
     
        if(rhmc_verbose.EQ.1) then
            print*, "calling rhmc"
        end if
        call print_time_step("start rhmc")
        call RHMC_evolution_device(xmat,alpha,phase,ncv,n_bad_CG,nacceptance,nbc,nbmn,&
            &temperature,flux,GAMMA10d,Gam123,ntau,nratio,dtau_xmat,dtau_alpha,&
            &acceleration,g_alpha,g_R,RCUT,&
            &acoeff_md,bcoeff_md,acoeff_pf,bcoeff_pf,max_err,max_iteration,iteration,&
            &ham_init,ham_fin,ntrial,imetropolis)
        if(rhmc_verbose.EQ.1) then
            print*, "RHMC finished traj=",itraj
        end if
        call print_time_step("finished rhmc")
        !**************************
        !**** measurements etc ****
        !**************************
        if(MOD(itraj,nskip).EQ.0)then
            call print_time_step("start measurements")
            !keep the hermiticity with a good precision.
            !Not really needed but just in case.
            if(rhmc_verbose.EQ.1) then
                print*, "calling Hermitian projection"
            end if
            call hermitian_projection_device(xmat)
            !$acc update host(xmat)
            !$acc update host(alpha)
            !measurements.
            if(rhmc_verbose.EQ.1) then
                print*, "calling measurements"
            end if
            ! These are the measurements on hostcode and device code
            call measure_host_device(xmat,alpha,nbc,nbmn,temperature,flux,&
                &GAMMA10d,neig_max,neig_min,ham_init,ham_fin,itraj,ntrial,&
                iteration,max_err,max_iteration,ncv,n_bad_CG,nacceptance,phase,Gam123)
           !call measurements(xmat,alpha,nbc,nbmn,temperature,flux,&
            !    &GAMMA10d,neig_max,neig_min,ham_init,ham_fin,itraj,ntrial,&
             !   iteration,max_err,max_iteration,ncv,n_bad_CG,nacceptance)
            !optimization of the Fourier acceleration parameters.
            !call Fourier_transform_xmat(xmat,&
            !     xmat_mom,x2p)
            !What does this mean?
            !call Fourier_acceleration_optimize(xmat_mom,fluctuation,&
            !     &imeasure,1)
            if(rhmc_verbose.EQ.1) then
                print*, "measurements done"
            end if
            call print_time_step("finished measurements")
        end if
        if((isave.EQ.0).AND.(nsave.GT.0))then
            if(MOD(itraj,nsave).EQ.0)then
                call print_time_step("start config writeout")
                !$acc update self(xmat)
                !$acc update self(alpha)
                if(rhmc_verbose.EQ.1) then
                    print*, "calling save config"
                end if
                call Save_Intermediate_Config(xmat,alpha,itraj)
                ! I have added this line in order to ensure that there are checkpoints.
                !call Save_Final_Config(xmat,alpha,itraj+1,checkpointn)
                call save_checkpoint(xmat,alpha,itraj+1)
                call print_time_step("finished configwriteout")
            end if
        end if
        call print_timestamp("trajectory finished")
        itraj=itraj+1
    end do

    !************************************************
    !************************************************
    !***  End of the loop of molecular evolutions ***
    !************************************************
    !************************************************
  
    !**********************************************************
    !*** close the output file for the measurement results. ***
    !**********************************************************
    close(unit_measurement)
    if(isave.EQ.0)then
        close(unit_intermediate_config)
    end if
    !*************************************
    !*** Save the final configuration ****
    !*************************************
    call print_time_step("start config writeout")
    !$acc update self(xmat)
    !$acc update self(alpha)
    call Save_Final_Config(xmat,alpha,itraj,output_config)
    call print_time_step("finished configwriteout")
    !$acc end data

    deallocate(acoeff_pf)
    deallocate(bcoeff_pf)
    deallocate(acoeff_md)
    deallocate(bcoeff_md)

    deallocate(xmat)
    call print_timestamp("end program")
end program BFSS_DEVICE
