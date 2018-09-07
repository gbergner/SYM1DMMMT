! distribut X^{(i,j)} to X^{(i,k)}, k=1,...,nblcok
SUBROUTINE mpi_xmat_row(xmat,xmat_row,myrank)
 
  implicit none
  
  include 'mpif.h'
  include 'size_parallel.h'
  !***** input *****
  integer myrank
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  !***** output *****
  double complex xmat_row(1:nmat_block,1:nmat_block*nblock,1:ndim,1:nsite_local)
  !***** For MPI *****
  integer IERR,IREQ,send_rank,receive_rank,tag
  integer STATUS(MPI_STATUS_SIZE)
  double complex xmat_rcv(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)
  !*******************
  integer imat,jmat
  integer idim
  integer isite
  integer iblock,jblock,kblock_send,kblock_rcv
  integer isublat
  integer ishift
  

  call who_am_i(myrank,isublat,iblock,jblock)

!$omp parallel
!$omp do
  do imat=1,nmat_block
     do jmat=1,nmat_block
        do idim=1,ndim   
           do isite=1,nsite_local
              xmat_row(imat,jmat+(jblock-1)*nmat_block,idim,isite)&
                   =xmat(imat,jmat,idim,isite)
           end do
        end do
     end do
  end do
!$omp end do
!$omp end parallel  

  kblock_send=jblock
  kblock_rcv=jblock
  do ishift=1,nblock-1
     kblock_send=kblock_send+1
     kblock_rcv=kblock_rcv-1
     if(kblock_send.EQ.nblock+1)then
        kblock_send=1
     end if
     if(kblock_rcv.EQ.0)then
        kblock_rcv=nblock
     end if
     
     send_rank=(isublat-1)*nblock*nblock+(iblock-1)*nblock+kblock_send-1
     receive_rank=(isublat-1)*nblock*nblock+(iblock-1)*nblock+kblock_rcv-1

     tag=1
     call MPI_Isend(xmat(1,1,1,1),nmat_block*nmat_block*nsite_local*ndim,&
          &MPI_DOUBLE_COMPLEX,&
          &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
     call MPI_Recv(xmat_rcv(1,1,1,1),nmat_block*nmat_block*nsite_local*ndim,&
          &MPI_DOUBLE_COMPLEX,&
          &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
     call MPI_Wait(ireq,status,ierr)

!$omp parallel
!$omp do
     do imat=1,nmat_block
        do jmat=1,nmat_block
           do idim=1,ndim   
              do isite=1,nsite_local
                 xmat_row(imat,jmat+(kblock_rcv-1)*nmat_block,idim,isite)&
                      =xmat_rcv(imat,jmat,idim,isite)
              end do
           end do
        end do
     end do
!$omp end do
!$omp end parallel
     
  end do

  return
  
END SUBROUTINE mpi_xmat_row
