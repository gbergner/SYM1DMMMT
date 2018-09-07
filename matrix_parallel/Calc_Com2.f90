subroutine Calc_Com2(xmat,com2,myrank)

  implicit none

  include 'mpif.h'
  include 'size_parallel.h'
  !***** input *****
  integer myrank
  double complex xmat(1:nmat_block,1:nmat_block,&
       &1:ndim,-(nmargin-1):nsite_local+nmargin),&
       &com(1:nmat_block,1:nmat_block)
  !***** output *****
  double precision com2
  !*******************
  double precision com2_local
  integer isite,idim,jdim
  integer imat,jmat,kmat
  double complex xmat_row(1:nmat_block,1:nmat_block*nblock,1:ndim,1:nsite_local)
  double complex xmat_column(1:nmat_block*nblock,1:nmat_block,&
       &1:ndim,1:nsite_local)
  !***** For MPI *****
  integer IERR

  !move i-th row and j-th row of xmat to (i,j)-th node.
  call mpi_xmat_row(xmat,xmat_row,myrank)
  call mpi_xmat_column(xmat,xmat_column,myrank)

  com2_local=0d0
  do isite=1,nsite_local 
     do idim=1,ndim-1
        do jdim=idim+1,ndim
           com=(0d0,0d0)
!$omp parallel
!$omp do
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 do kmat=1,nmat_block*nblock
                    com(imat,jmat)=com(imat,jmat)&
                         &+xmat_row(imat,kmat,idim,isite)&
                         &*xmat_column(kmat,jmat,jdim,isite)&
                         &-xmat_row(imat,kmat,jdim,isite)&
                         &*xmat_column(kmat,jmat,idim,isite)
                 end do
              end do
           end do
!$omp end do
!$omp end parallel
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 com2_local=com2_local&
                      &+dble(com(imat,jmat)*dconjg(com(imat,jmat)))
              end do
           end do
        end do
     end do
  end do
  call MPI_Reduce(com2_local,com2,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,0,MPI_COMM_WORLD,IERR)
  if(myrank.EQ.0)then
     com2=com2/dble(nmat_block*nblock*nsite_local*nsublat)*2d0
  end if
  

  return

END subroutine Calc_Com2
