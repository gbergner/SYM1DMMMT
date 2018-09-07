  !*************************************
  !*************************************
  !*** Set the initial configuration ***
  !*************************************
  !*************************************

SUBROUTINE initial_configuration(xmat,alpha,acceleration,itraj,init,iaccelerate,nfuzzy,input_config,acc_input,flux,mersenne_seed)

  use mtmod !Mersenne twistor
  implicit none

  include 'size_parallel.h'
  include 'mpif.h'
  include 'unit_number.inc'
  !****** input ******
  integer init,iaccelerate,nfuzzy,mersenne_seed
  character(1000) input_config,acc_input
  double precision flux
  !****** output ******
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,-(nmargin-1):nsite_local+nmargin)
  double precision alpha(1:nmat_block*nblock)
  double precision acceleration(1:nsite_local)
  integer itraj
  !********************
  double complex, allocatable :: xmat_init(:,:,:,:)
  double precision, allocatable :: acceleration_init(:),acceleration_init_2(:)
  integer myrank,nprocs,ierr
  integer imat,ispin,jspin,ifuzzy,isite,isublat,i
  double precision spin

  call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS, IERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK, IERR)

  allocate(xmat_init(1:nmat_block*nblock,1:nmat_block*nblock,&
       &1:ndim,1:nsite_local*nsublat))
  if(MYRANK.EQ.0)then
     allocate(acceleration_init(1:nsite_local*nsublat*nblock*nblock))
     allocate(acceleration_init_2(1:nsite_local*nsublat))
     
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
        do imat=1,nmat_block*nblock
           alpha(imat)=-1d0+2d0*dble(imat)/dble(nmat_block*nblock)
        end do
     else if(init.EQ.2)then
        !fuzzy sphere; # of fuzzy spheres = nfuzzy 
        !spin of the fuzzy sphere = (1/2)*(nmat/nfuzzy-1)
        itraj=1
        xmat_init=(0d0,0d0)
        spin=0.5d0*(dble(nmat_block*nblock)/dble(nfuzzy)-1d0)
        jspin=int(dble(nmat_block*nblock)/dble(nfuzzy)+0.1d0)!size of irreducible repr, 2*spin+1
        do isite=1,nsite_local*nsublat
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
        do imat=1,nmat_block*nblock
            alpha(imat)=-1d0+2d0*dble(imat)/dble(nmat_block*nblock)
        end do

     end if
     !*******************************************************
     !*** Initial Configuration is now stored in myrank=0 ***
     !*******************************************************
  end if

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
     
     if(myrank.EQ.0)then
        open(unit=unit_input_acc,status='OLD',file=acc_input,action='READ')
        read(unit_input_acc,*) acceleration_init_2
        close(unit_input_acc)
        !We have to change the format
        do isite=1,nsite_local
           do isublat=1,nsublat
              do i=1,nblock*nblock
                 acceleration_init(isite+(i-1)*nsite_local&
                      &+(isublat-1)*nsite_local*nblock*nblock)=&
                      &acceleration_init_2(isite+(isublat-1)*nsite_local)
              end do
           end do
        end do
     end if
     !distribute the Fourier acceleration parameters to nodes. 
     call MPI_scatter(acceleration_init(1),nsite_local,&
          &MPI_DOUBLE_PRECISION,&
          &acceleration(1),nsite_local,MPI_DOUBLE_PRECISION,&
          &0,MPI_COMM_WORLD,IERR)
     
  else if(iaccelerate.EQ.1)then
     !"naive choice" of the Fourier acceleration parameters. 
     call Fourier_acceleration_naive(acceleration,myrank)
  end if
  
  !distribute the configuration to all ranks
  call matrix_format_change(xmat,xmat_init,myrank,2)

  deallocate(xmat_init)
  call MPI_Bcast(alpha(1),nmat_block*nblock,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_Bcast(itraj,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  
  !adjust the margin
  call Adjust_margin_xmat(xmat,myrank)

  if(MYRANK.EQ.0)then
     deallocate(acceleration_init)
     deallocate(acceleration_init_2)
  end if

  return

end SUBROUTINE initial_configuration
