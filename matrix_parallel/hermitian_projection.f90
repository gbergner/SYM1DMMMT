SUBROUTINE hermitian_projection(xmat,myrank)
  
  implicit none
  
  include 'mpif.h'
  include 'size_parallel.h'
  !***** input *****
  integer myrank
  !***** input & output *****
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  !**************************
  double complex mat_rcv(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  integer idim,imat,jmat,isite
  integer isublat,iblock,jblock
  !***** for MPI *****
  integer send_rank,receive_rank,ireq,ierr,tag
  integer status(MPI_STATUS_SIZE)

  call who_am_i(myrank,isublat,iblock,jblock)

  if(iblock.EQ.jblock)then
     do idim=1,ndim
        do isite=-(nmargin-1),nsite_local+nmargin
!$omp parallel
!$omp do
           do imat=1,nmat_block-1
              do jmat=imat+1,nmat_block
                 xmat(jmat,imat,idim,isite)=&
                      &dconjg(xmat(imat,jmat,idim,isite))
                 
              end do
           end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do
           do imat=1,nmat_block
              xmat(imat,imat,idim,isite)=&
                   &dcmplx(dble(xmat(imat,imat,idim,isite)))
           end do
!$omp end do
!$omp end parallel       
        end do
     end do

  else 
     send_rank=(isublat-1)*nblock*nblock+(jblock-1)*nblock+iblock-1
     receive_rank=(isublat-1)*nblock*nblock+(jblock-1)*nblock+iblock-1

     tag=9
     call MPI_Isend(xmat(1,1,1,-(nmargin-1)),&
          nmat_block*nmat_block*(nsite_local+2*nmargin)*ndim,&
          &MPI_DOUBLE_COMPLEX,&
          &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
     call MPI_Recv(mat_rcv(1,1,1,-(nmargin-1)),&
          nmat_block*nmat_block*(nsite_local+2*nmargin)*ndim,&
          &MPI_DOUBLE_COMPLEX,&
          &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
     call MPI_Wait(ireq,status,ierr)
 
     do idim=1,ndim
        do isite=-(nmargin-1),nsite_local+nmargin
!$omp parallel
!$omp do
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 xmat(imat,jmat,idim,isite)=&
                      &dconjg(mat_rcv(jmat,imat,idim,isite))
              end do
           end do
!$omp end do
!$omp end parallel
        end do
     end do
     
     
  end if
  
  
  return
  
END SUBROUTINE 
