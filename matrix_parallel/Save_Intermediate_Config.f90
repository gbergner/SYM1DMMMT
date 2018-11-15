  !*************************************
  !*** Save the final configuration ****
  !*************************************

SUBROUTINE Save_Intermediate_Config(xmat,alpha,itraj)
  
  use mtmod !Mersenne twistor
  implicit none
  
  include 'size_parallel.h'
  include '../unit_number.inc'
  include 'mpif.h'
  !***** input *****
  double precision alpha(1:nmat_block*nblock)
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  integer itraj
  !*******************
  double complex, allocatable :: xmat_int(:,:,:,:)
  integer myrank,IERR

  call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK, IERR)

  allocate(xmat_int(1:nmat_block*nblock,1:nmat_block*nblock,&
       &1:ndim,1:nsite_local*nsublat))
    
  call matrix_format_change(xmat,xmat_int,myrank,1)
  
  if(myrank.eq.0)then
     call mtsaveu(unit_intermediate_config)  
     write(unit_intermediate_config,*) itraj
     write(unit_intermediate_config,*) xmat_int
     write(unit_intermediate_config,*) alpha
  end if

  deallocate(xmat_int)


  return

END SUBROUTINE Save_Intermediate_Config

