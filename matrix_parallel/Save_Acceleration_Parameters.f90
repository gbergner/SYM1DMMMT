  !*************************************
  !*** Save the final configuration ****
  !*************************************

SUBROUTINE Save_Acceleration_Parameters(fluctuation,imeasure,acc_output)
  
  use mtmod !Mersenne twistor
  implicit none
  
  include 'size_parallel.h'
  include 'unit_number.inc'
  include 'mpif.h'
  !---------------------------------
  double complex xmat_mom(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local) 
  double precision fluctuation(1:nsite_local) 
  integer imeasure
  integer isite,isublat
  double precision, allocatable :: acceleration_fin(:)
  character(1000) acc_output

  integer myrank,nprocs,IERR,tag,ireq
  integer STATUS(MPI_STATUS_SIZE)

  call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS, IERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK, IERR)

  if(myrank.eq.0)then 
     allocate(acceleration_fin(1:nsite_local*nsublat))
  end if
    
  call Fourier_acceleration_optimize(xmat_mom,fluctuation,myrank,imeasure,2)
  

  if(myrank.EQ.0)then
     do isite=1,nsite_local
        acceleration_fin(isite)=fluctuation(isite)
     end do
  end if
  do isublat=2,nsublat
     tag=100
     if(myrank.EQ.(isublat-1)*nblock*nblock)then
        call MPI_Isend(fluctuation(1),nsite_local,&
             MPI_DOUBLE_PRECISION,&
             0,tag,MPI_COMM_WORLD,ireq,ierr)
        call MPI_Wait(ireq,status,ierr) 
     end if
     if(myrank.EQ.0)then
        call MPI_Recv(acceleration_fin((isublat-1)*nsite_local+1),&
             &nsite_local,&
             &MPI_DOUBLE_PRECISION,&
             &(isublat-1)*nblock*nblock,tag,MPI_COMM_WORLD,status,ierr)
     end if
     !call MPI_Wait(ireq,status,ierr) 
  end do
 

  if(myrank.eq.0)then
     
     open(UNIT=unit_output_acc, File = acc_output, STATUS = "REPLACE", ACTION = "WRITE")
     write(unit_output_acc,*) acceleration_fin
     close(unit_output_acc)
     deallocate(acceleration_fin)
  end if



  return

END SUBROUTINE Save_Acceleration_Parameters

