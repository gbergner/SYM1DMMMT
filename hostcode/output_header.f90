  !*************************************************************************
  !******  Open the output file (for measurements) and set the header ******
  !*************************************************************************

SUBROUTINE output_header(data_output,temperature,flux,&
    &ntau,nratio,dtau_xmat,dtaU_alpha,neig_max,neig_min,nbc,nbmn,&
    &init,input_config,output_config,iaccelerate,acc_input,acc_output,&
    &g_alpha,g_R,RCUT,upper_approx,max_err,max_iteration,CG_log,&
    &isave,nsave,intermediate_config,imetropolis)

    implicit none
  include '../staticparameters.f90'
  include '../unit_number.inc'
    integer isave,nsave
    character(1000) data_output,input_config,output_config,acc_input,acc_output,&
        &intermediate_config,CG_log
    character(100)input_config_short,output_config_short,&
        &acc_input_short,acc_output_short,intermediate_config_short,CG_log_short
    double precision temperature
    double precision flux
    double precision g_alpha,g_R,RCUT,upper_approx,max_err
    integer max_iteration
    integer neig_max,neig_min,myrank,nprocs,nbmn,nbc,ntau,nratio,init,iaccelerate,imetropolis
    double precision dtau_xmat,dtau_alpha
    logical exist
  
    !************************************************************
    !*** They are too long to write in the output file... :' ****
    !************************************************************
    input_config_short=input_config
    output_config_short=output_config
    acc_input_short=acc_input
    acc_output_short=acc_output
    intermediate_config_short=intermediate_config
    CG_log_short=CG_log
    !*******************
    !*** data_output ***
    !*******************
    inquire(file=data_output,exist=exist)
    if(exist.AND.(init.EQ.4)) then
        open(unit=unit_measurement,status='OLD',file=data_output,&
            &action='WRITE',position='APPEND')
    else
        open(unit=unit_measurement,status='REPLACE',file=data_output,&
            &action='WRITE')
        if(nbc.EQ.0)then
            write(unit_measurement,*)"#pbc"
        else if(nbc.EQ.1)then
            write(unit_measurement,*)"#apbc"
        end if
        if(nbmn.EQ.0)then
            write(unit_measurement,*)"#BFSS matrix model"
        else if(nbmn.EQ.1)then
            write(unit_measurement,*)"#BMN matrix model"
        end if
        if(nimprove.EQ.0)then
            write(unit_measurement,*) "#naive lattice action"
        else if(nimprove.EQ.1)then
            write(unit_measurement,*) "#improved lattice action"
        end if
        write(unit_measurement,*) "#size of the gauge group: nmat=",nmat
        write(unit_measurement,*) "#size of the lattice: Nt=",nsite
        write(unit_measurement,*) "#Temperature=",temperature
        if(nbmn.EQ.1)then
            write(unit_measurement,*)"#flux parameret: mu=", flux
        end if
        write(unit_measurement,*) "#g_alpha=",g_alpha
        write(unit_measurement,*) "#g_R=",g_R
        write(unit_measurement,*) "#RCUT=",RCUT
        write(unit_measurement,*) "#(TrX^2 is constrained so that TrX^2 < RCUT)"
        write(unit_measurement,*) "#ntau=",ntau
        write(unit_measurement,*) "#nratio=",nratio
        write(unit_measurement,*) "#dtau for xmat=",Dtau_xmat
        write(unit_measurement,*) "#dtau for alpha=",Dtau_alpha
        write(unit_measurement,*) "#CG stopping condition: relative error=",max_err
        write(unit_measurement,*) "#CG stopping condition: maximum CG iteration=",&
            max_iteration
        write(unit_measurement,*) "#upper limit of remez approximation &
       &(Largest eigenvalue of D^dag*D must be smaller than this):"        ,&
            &upper_approx
        write(unit_measurement,*) "#init=",init
        write(unit_measurement,*) "#(init; 0 -> old config, 1 -> new config, &
       &2-> fuzzy sphere)"        
        write(unit_measurement,*) "#imetropolis=",imetropolis
        write(unit_measurement,*) "#(imetropolis=1 -> no Metropolis test,&
       &just for thermalization)"        
        if(init.EQ.0)then
            write(unit_measurement,*) "#input configuration:",input_config_short
        end if
        write(unit_measurement,*) "#output configuration:",output_config_short
        write(unit_measurement,*) "#iaccelerate=",iaccelerate
        write(unit_measurement,*) "#(iaccelerate; 0 -> read from acc_input, &
       &1-> naive)"        
        if(iaccelerate.EQ.0)then
            write(unit_measurement,*) "#Fourier acceleration file &
          &used for this run:"            ,&
                &acc_input_short
        end if
        write(unit_measurement,*) "#Fourier acceleration file obtained &
       &from this run:"        ,acc_output_short
        if(isave.EQ.0)then
            write(unit_measurement,*) "#intermediate configurations are saved in:",&
                &intermediate_config_short
            write(unit_measurement,*) "#saved every nsave configurations, nsave=",&
                nsave
        end if
  
        write(unit_measurement,*) "#log of the CG solver is saved in:",&
            &CG_log_short
        !********************************************************
        !**** Explain the ordering of the measurement output ****
        !********************************************************
        if((neig_max.EQ.0).AND.(neig_min.EQ.0))then
            !The largest and smallest eigenvalues of (D^dag*D) are NOT calculated.
            write(unit_measurement,*) "# traj,ham_fin-ham_init, &
          &number of constraint violation,&
          &number of CG-err,CG-iteration,energy, &
          &|Pol.|, trX^2, tr[X,Y]^2,acceptance"            
     
        else  if((neig_max.GT.0).AND.(neig_min.EQ.0))then
            !The largest eigenvalue of (D^dag*D) is calculated.
            write(unit_measurement,*)"# iteration for evaluating largest eig=",&
                &neig_max
            write(unit_measurement,*) "# traj,ham_fin-ham_init, &
          &number of constraint violation,&
          &number of CG-err,CG-iteration,energy, &
          &|Pol.|, trX^2, tr[X,Y]^2, largest eig, acceptance"            
     
        else  if((neig_max.EQ.0).AND.(neig_min.GT.0))then
            !The smallest eigenvalue of (D^dag*D) is calculated.
            write(unit_measurement,*)"# iteration for evaluating smallest eig=",&
                &neig_min
            write(unit_measurement,*) "# traj,ham_fin-ham_init, &
          &number of constraint violation,&
          &number of CG-err,CG-iteration,energy, &
          &|Pol.|, trX^2, tr[X,Y]^2, smallest eig, acceptance"            
     
        else  if((neig_max.GT.0).AND.(neig_min.GT.0))then
            !Both the largest and smallest eigenvalues of (D^dag*D) are calculated.
            write(unit_measurement,*)"# iteration for evaluating largest eig=",&
                &neig_max
            write(unit_measurement,*)"# iteration for evaluating smallest eig=",&
                &neig_min
            write(unit_measurement,*) "# traj,ham_fin-ham_init, &
          &number of constraint violation,&
          &number of CG-err,CG-iteration,energy, &
          &|Pol.|, trX^2, tr[X,Y]^2, largest eig,&
          &smallest eig, acceptance"            
        end if
        write(unit_measurement,*)'#-----------------------------------------------'
  
    end if
  
    !**************
    !*** CG_log ***
    !**************
    inquire(file=CG_log,exist=exist)
    if(exist) then
        open(unit=unit_CG_log,status='OLD',file=CG_log,&
            &action='WRITE',position='APPEND')
    else
        open(unit=unit_CG_log,status='REPLACE',file=CG_log,&
            &action='WRITE')
    end if

    call flush(unit_measurement)
    return
END SUBROUTINE output_header
