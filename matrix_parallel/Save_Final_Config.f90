  !*************************************
  !*** Save the final configuration ****
  !*************************************

SUBROUTINE Save_Final_Config(xmat,alpha,itraj,output_config)
  
  use mtmod !Mersenne twistor
  implicit none
  
  include 'size_parallel.h'
  include 'unit_number.inc'
  include 'mpif.h'
  !***** input *****
  double precision alpha(1:nmat_block*nblock)
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  integer itraj
  character(1000) output_config
  !*******************
  double complex, allocatable :: xmat_fin(:,:,:,:)
  integer myrank,IERR

  call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK, IERR)

  allocate(xmat_fin(1:nmat_block*nblock,1:nmat_block*nblock,&
       &1:ndim,1:nsite_local*nsublat))
    
  call matrix_format_change(xmat,xmat_fin,myrank,1)
  

  if(myrank.eq.0)then

     open(UNIT=unit_output_config, File = output_config,&
          &STATUS = "REPLACE", ACTION = "WRITE")
     call mtsaveu(unit_output_config)  
     write(unit_output_config,*) itraj
     write(unit_output_config,*) xmat_fin
     write(unit_output_config,*) alpha
     close(unit_output_config)
     
  end if

  deallocate(xmat_fin)


  return

END SUBROUTINE Save_Final_Config

