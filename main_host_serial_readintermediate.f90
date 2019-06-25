!##############################################################################
!######              BFSS matrix model on lattice                     #########
!######                                                               #########
!######                -original written by Masanori Hanada           #########
!######                -reorganized by Georg Bergner                  #########
!######                -most parts of this CPU version kept           #########
!######                 as in the original to have a reference        #########
!######                 for comparisons                               #########
!######                                                               #########
!######            ver.1  bosonic, HMC, no parallelization.           #########
!######            ver.2  gauge fixed version is available.           #########
!######            ver.3  parallel code.                              #########
!######            ver.4  fermionic part is introduced.               #########
!######                   trace part of the fermion is left as it is. #########
!######            ver.5  U(1) part of pseudo-fermion is removed.     #########
!######                   potential term for alpha, which enforce     #########
!######                   them to satisfy the constraint, is introduced. ######
!######            ver.6  Constraint for alpha is improved.           #########
!######            ver.7  Constraint for TrX^2 is introduced.         #########
!######            ver.8  Fourier acceleration is introduced.         #########
!######                   Does not use Fast Fourier Transform yet.    #########
!######            ver.9  Flux introduced. nbmn=0 -> BFSS, nbmn=1 -> BMN ######
!######            ver.10 Mersenne twistor is introduced.             #########
!######                   (Taken from M.Matsumoto's webpage @ Hiroshima-U.) ###
!######            ver.11 "Fourier_acceleration_optimize" was at wrong place.##
!######                   Fast Fourier Transform is introduced. (A.Tsuchiya) ##
!######                   Debug (flux part of Dirac operator)         #########
!######            ver.12 Dirac operator is modified.                 #########
!######            ver.13 O(a^2) improvement is introduced.           #########
!######            ver.14 Increase the readability. (with I.Kanamori) #########
!######            ver.15 Make it easier to read. (with I.Kanamori)   #########
!######            ver.16 Minor improvement in output format.         #########
!######                   Intermediate configurations can be saved.   #########
!######            ver.17 Now you can take a log of CG solver.        #########
!##############################################################################
!######            ver.1  Serial version, rewritten from ver.17.      #########
!######            ver.2  Dirac multiplication made explicit.         #########
!######            ver.3  Multiple pseudo fermion.                    #########
!######            ver.4  Multi time scales.                          #########
!######            ver.5  minor improvements.                         #########
!##############################################################################
!Mersenne twister.
!include 'mt19937.f90'

program BFSS
  
  use mtmod !Mersenne twistor
  implicit none
  
  include 'staticparameters.f90'
  include 'Fourier.inc'
  include 'include.h' 
  include 'unit_number.inc'
  integer status
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
  read(unit_input_para,*) ngauge
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
  read(unit_input_para,*) purebosonic
  read(unit_input_para,*) Pol_phase
  read(unit_input_para,*) Eigenval
  close(unit_input_para)
  !Construc Gamma matrices. 
  call MakeGamma(Gamma10d)

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

  status=0
  inquire(file=intermediate_config,exist=status)
  if(.not.status) then
   print *, "ERROR: File ",trim(intermediate_config)," does not exist."
   return
  end if
  status=0
     
  !*************************************
  !*** Set the initial configuration ***
  !*************************************
  call initial_configuration(xmat,alpha,acceleration,itraj,init,&
       &iaccelerate,nfuzzy,input_config,acc_input,flux,mersenne_seed,ngauge)

  if(ngauge.EQ.1)then
     !ungauged
     alpha=0d0
     dtau_alpha=0d0
  end if
  !***********************************
  !******  Make the output file ******
  !***********************************
  call output_header(data_output,temperature,flux,&
       &ntau,nratio,dtau_xmat,dtaU_alpha,neig_max,neig_min,nbc,nbmn,&
       &init,input_config,output_config,iaccelerate,acc_input,acc_output,&
       &g_alpha,g_R,RCUT,upper_approx,max_err,max_iteration,CG_log,Pol_phase,Eigenval,&
       &isave,nsave,intermediate_config,imetropolis,ngauge,purebosonic)

  !*******************************************************
  !******  Make the intermediate configuration file ******
  !*******************************************************

  open(unit=unit_intermediate_config,&
          &file=intermediate_config,action='READ')
  ! everything not accessible from the updater is set to zero.
  nacceptance=0 !number of acceptance
  ntrial=0 !number of trial 
  ncv=0 !number of the constraint violation.
  n_bad_CG=0 !count how many times CG solver failed to converge.
  iteration=0
  ham_init=0
  ham_fin=0
  do while ((itraj.LE.ntraj).AND.(status.GE.0))
     call Read_Intermediate_Config(xmat,alpha,itraj,status)
     if((MOD(itraj,nskip).EQ.0).AND.(status.GE.0))then
        call hermitian_projection(xmat)
        call measurements(xmat,alpha,nbc,nbmn,temperature,flux,&
             &GAMMA10d,neig_max,neig_min,ham_init,ham_fin,itraj,ntrial,&
             iteration,max_err,max_iteration,ncv,n_bad_CG,nacceptance,ngauge,purebosonic)
     end if
     itraj=itraj+1
  end do
  close(unit_measurement)
  close(unit_intermediate_config)
  
  
end program BFSS
