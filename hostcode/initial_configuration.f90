  !*************************************
  !*************************************
  !*** Set the initial configuration ***
  !*************************************
  !*************************************

SUBROUTINE initial_configuration(xmat,alpha,acceleration,itraj,init,iaccelerate,nfuzzy,input_config,acc_input,flux,mersenne_seed)

  use mtmod !Mersenne twistor
  implicit none

  include '../staticparameters.f90'
  include '../unit_number.inc'
  !****** input ******
  integer init,iaccelerate,nfuzzy,mersenne_seed
  character(1000) input_config,acc_input
  double precision flux
  !****** output ******
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double precision alpha(1:nmat)
  double precision acceleration(1:nsite)
  integer itraj
  !********************

  double complex xmat_init(1:nmat,1:nmat,1:ndim,1:nsite)
  integer myrank,nprocs,ierr
  integer imat,jmat,ispin,jspin,ifuzzy,isite,idim
  double precision spin

  !init=0 -> continue from old config
  !init=1 -> new config, cold start  
  !init=2 -> fuzzy sphere; # of fuzzy spheres = nfuzzy
  if(init.EQ.0)then
     !continue from old config   
     open(unit=unit_input_config,status='OLD',file=input_config,action='READ')
     call mtgetu(unit_input_config)
     read(unit_input_config,*) itraj
     read(unit_input_config,*) xmat_init
     read(unit_input_config,*) alpha
     close(unit_input_config)
  else if(init.EQ.1)then
     !new config, cold start  
     itraj=1
     xmat_init=(0d0,0d0)
     !because alpha_1=...=alpha_N=0 is singular, 
     !we must slightly shift them. 
     do imat=1,nmat
        alpha(imat)=-1d0+2d0*dble(imat)/dble(nmat)
     end do
  else if(init.EQ.2)then
     !fuzzy sphere; # of fuzzy spheres = nfuzzy 
     !spin of the fuzzy sphere = (1/2)*(nmat/nfuzzy-1)
     itraj=1
     xmat_init=(0d0,0d0)
     spin=0.5d0*(dble(nmat)/dble(nfuzzy)-1d0)
     jspin=int(dble(nmat)/dble(nfuzzy)+0.1d0)!size of irreducible repr, 2*spin+1
     do isite=1,nsite
        do ifuzzy=1,nfuzzy
           do ispin=1,jspin
              imat=(ifuzzy-1)*jspin+ispin
              xmat_init(imat,imat,3,isite)=dcmplx(spin-dble(ispin-1))*dcmplx(flux)
           end do
           do ispin=1,jspin-1
              imat=(ifuzzy-1)*jspin+ispin
              xmat_init(imat,imat+1,1,isite)&
                   &=(0.5d0,0d0)*dcmplx(flux)&
                   &*dcmplx(dsqrt(spin*(spin+1d0)-(spin-dble(ispin))*(spin-dble(ispin-1))))
              xmat_init(imat+1,imat,1,isite)&
                   &=(0.5d0,0d0)*dcmplx(flux)&
                   &*dcmplx(dsqrt(spin*(spin+1d0)-(spin-dble(ispin))*(spin-dble(ispin-1))))
              xmat_init(imat,imat+1,2,isite)&
                   &=(0.0d0,-0.5d0)*dcmplx(flux)&
                   &*dcmplx(dsqrt(spin*(spin+1d0)-(spin-dble(ispin))*(spin-dble(ispin-1))))
              xmat_init(imat+1,imat,2,isite)&
                   &=(0.0d0,0.5d0)*dcmplx(flux)&
                   &*dcmplx(dsqrt(spin*(spin+1d0)-(spin-dble(ispin))*(spin-dble(ispin-1))))
           end do
        end do
     end do
     !because alpha_1=...=alpha_N=0 is singular, 
     !we must slightly shift them. 
     do imat=1,nmat
        alpha(imat)=-1d0+2d0*dble(imat)/dble(nmat)
     end do
  end if
  if(ngauge.EQ.0)then
     alpha=0d0
  end if
  !**************************
  !**************************
  !*** Set random numbers ***
  !**************************
  !**************************
  if(init.EQ.0)then
     !continue from old config. 
     !all nodes must access to the information about the random number. 
     open(unit=unit_input_config,status='OLD',file=input_config,action='READ')
     call mtgetu(unit_input_config)
     close(unit_input_config)
  else if(init.NE.0)then
     call sgrnd(mersenne_seed)!set a seed for the Mersenne twistor at each node.
  end if

  if(iaccelerate.EQ.0)then
     !read acceleration parameters from existing file.
     if(iaccelerate.EQ.0)then
        open(unit=unit_input_acc,status='OLD',file=acc_input,action='READ')
        read(unit_input_acc,*) acceleration
        close(unit_input_acc)     
     end if
  else if(iaccelerate.EQ.1)then
     !"naive choice" of the Fourier acceleration parameters. 
     call Fourier_acceleration_naive(acceleration)
  end if
  
  do isite=1,nsite  
     do idim=1,ndim
        do imat=1,nmat
           do jmat=1,nmat
              xmat(imat,jmat,idim,isite)= xmat_init(imat,jmat,idim,isite)
           end do
        end do
     end do
  end do
  !adjust margin.
  call Adjust_margin_xmat(xmat)

 
  return

END SUBROUTINE initial_configuration
