!##############################################################################
!######              BFSS matrix model on lattice                     #########
!######                                                               #########
!######              This is a version of the providing the same
!######              environment as the main serial update code
!######              for the tests comparing host and device code.
!##############################################################################

program TEST_DEVICE_HOST
  
  use mtmod !Mersenne twistor
  use compare_host_device
  implicit none
  
  include 'staticparameters.f90'
  include 'Fourier.inc'
  include 'include.h' 
  include 'unit_number.inc'
  if(npf.EQ.1)then
     include 'remez_md_1.dat'
     include 'remez_pf_1.dat'
  else if(npf.EQ.2)then
     include 'remez_md_2.dat'
     include 'remez_pf_2.dat'
  end if
     
     
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

   call compare_host_device(xmat,alpha,ncv,n_bad_CG,nacceptance,nbc,nbmn,&
     &temperature,flux,GAMMA10d,ntau,nratio,dtau_xmat,dtau_alpha,&
     &acceleration,g_alpha,g_R,RCUT,&
     &acoeff_md,bcoeff_md,acoeff_pf,bcoeff_pf,max_err,max_iteration,iteration,&
     &ham_init,ham_fin,ntrial,imetropolis)

  
end program TEST_DEVICE_HOST
