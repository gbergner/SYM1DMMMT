SUBROUTINE RHMC_evolution(xmat,alpha,ncv,n_bad_CG,nacceptance,nbc,nbmn,&
    &temperature,flux,GAMMA10d,ntau,nratio,dtau_xmat,dtau_alpha,&
    &acceleration,g_alpha,g_R,RCUT,&
    &acoeff_md,bcoeff_md,acoeff_pf,bcoeff_pf,max_err,max_iteration,iteration,&
    &ham_init,ham_fin,ntrial,imetropolis,ngauge,purebosonic)

    use mtmod !Mersenne twistor
    implicit none

    include '../staticparameters.f90'
    include '../Fourier.inc'
    include '../unit_number.inc'

    !input
    integer nbc,nbmn,ngauge,purebosonic
    double precision temperature,flux
    double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
    double precision max_err
    integer max_iteration
    double precision acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)
    double precision acoeff_pf(0:nremez_pf),bcoeff_pf(1:nremez_pf)
    double precision g_alpha,g_R,RCUT
    integer ntau,nratio
    integer imetropolis!1-> no Metropolis test
    double precision acceleration(1:nsite),dtau_alpha,dtau_xmat
    !input & output
    double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    double precision alpha(1:nmat)
    integer ncv,n_bad_CG,nacceptance,ntrial
    double precision ham_init,ham_fin
    !output
    integer iteration


    double precision metropolis
    integer info_pf,info_mol,info_CG_init,info_CG_fin,info_alpha,info_accept,info
    double complex backup_xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    double precision backup_alpha(1:nmat)
    double complex P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
    double precision P_alpha(1:nmat)
    double complex pf(1:nmat,1:nmat,&
        1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
    double complex Chi(1:nremez_md,1:nmat,1:nmat,&
        1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
    integer IERR,myrank,nprocs

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
    ! all input on device execute on device
    if(purebosonic.eq.0) then
        call generate_pseudo_fermion_SUN(pf,xmat,alpha,&
            &GAMMA10d,acoeff_pf,bcoeff_pf,max_err,max_iteration,iteration,&
            &nbc,nbmn,temperature,flux,info_pf)
    end if
    !info_pf=0 -> OK (CG solver converged)
    !info_pf=1 -> error (CG solver did not converge)
  
    !Take CG_log -> not device
    write(unit_CG_log,*)"pseudo-fermion generation",iteration
    !******************************
    !**** Molecular Evolution. ****
    !******************************
    !save the old config.
    ! on device
    backup_xmat=xmat
    backup_alpha=alpha
    !generate auxiliary momenta.
    !all input on device
    call Generate_P_xmat(P_xmat)
    ! device
    call Generate_P_alpha(P_alpha)

    if(rhmc_verbose.EQ.1) then
        print*, "check  host pseudofermion initial ", Sum(pf)
        print*, "check host xmat initial ", Sum(xmat)
        print*, "check host momenta xmat initial ", Sum(P_xmat)
        print*, "check host momenta alpha initial ",Maxval(abs(P_alpha))
    end if

    !Calculate ham_init
    if(purebosonic.eq.0) then
        call solver_biCGm(nbc,nbmn,nremez_md,&
            xmat,alpha,pf,chi,GAMMA10d,&
            bcoeff_md,max_err,max_iteration,iteration,&
            temperature,flux,info_CG_init)
    end if
    !info_CG_init=0 -> OK (CG solver converged)
    !info_CG_init=1 -> error (CG solver did not converge)
  
    !Take CG_log
    write(unit_CG_log,*)"ham_init",iteration
    ! requires summation and distribution
    call Calc_Ham(temperature,xmat,alpha,P_xmat,P_alpha,ham_init,pf,chi,&
        &acoeff_md,g_R,RCUT,nbmn,flux,ngauge,purebosonic)

    if(rhmc_verbose.EQ.1) then
        print *,"check initial Ham ",ham_init," inverter res ",Sum(chi)
    end if

    !Molecular Evolution
    call Molecular_Dynamics(nbc,temperature,&
        &ntau,nratio,dtau_xmat,dtau_alpha,xmat,alpha,P_xmat,P_alpha,&
        &acoeff_md,bcoeff_md,pf,max_iteration,max_err,iteration,&
        &gamma10d,g_alpha,g_R,RCUT,acceleration,nbmn,flux,info_mol,ngauge,purebosonic)
    !info_mol=0 -> OK (CG solver converged)
    !info_mol=1 -> error (CG solver did not converge)
    if(info_mol.EQ.0)then
        !calculate ham_fin
        if(purebosonic.eq.0) then
            call solver_biCGm(nbc,nbmn,nremez_md,&
                xmat,alpha,pf,chi,GAMMA10d,&
                bcoeff_md,max_err,max_iteration,iteration,&
                temperature,flux,info_CG_fin)
        end if
        !info_CG_fin=0 -> OK (CG solver converged)
        !info_CG_fin=1 -> error (CG solver did not converge)
     
        !Take CG_log
        write(unit_CG_log,*)"ham_fin",iteration
     
        call Calc_Ham(temperature,xmat,alpha,P_xmat,P_alpha,ham_fin,&
            &pf,chi,acoeff_md,g_R,RCUT,nbmn,flux,ngauge,purebosonic)
     
    end if
  
    if(rhmc_verbose.EQ.1) then
        print*, "check  device pseudofermion final ", Sum(pf)
        print*, "check device xmat final ", Sum(xmat)
        print*, "check device momenta xmat final ", Sum(P_xmat)
        print*, "check device momenta alpha final ", Maxval(abs(P_alpha))
        print *,"check final Ham ",ham_fin," inverter res ",Sum(chi)
    end if

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
    end if
     
    !####################################
    !#### check constraint for alpha ####
    !####################################
    call check_alpha_constraint(alpha,info_alpha)
    !Count how many times the constraint is violated.
    if(info_alpha.EQ.1)then
        ncv=ncv+1
        info=1
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
    end if

    if(info_accept.EQ.0)then
        !accept
        nacceptance=nacceptance+1
        ! Project out the trace parts of X and alpha.
        ! For X, we remove only the zero mode of X, 
        !        i.e. (¥int dt X) is made traceless.   
        call subtract_U1(xmat,alpha)
    else
        !reject
        xmat=backup_xmat
        alpha=backup_alpha
    end if
    if(rhmc_verbose.EQ.1) then
        print*,"RHMC ev finished accepted ",info_accept," ",info_alpha," ",metropolis, " ", info
    end if
    return

END SUBROUTINE RHMC_evolution
