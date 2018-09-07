!iformat=1 -> usual matrix form (no-block) is srored at myrank=0
! i.e. block matrices to big matrix
!iformat=2 -> big matrix to block matrices
subroutine matrix_format_change(xmat,xmat_big,myrank,iformat)

  implicit none

  include 'mpif.h'
  include 'size_parallel.h'
  !***** input *****
  integer myrank,iformat
  !***** input when iformat=1, output when iformat=2*****
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
       -(nmargin-1):nsite_local+nmargin)
  !***** output when iformat=1, input when iformat=2*****
  double complex xmat_big(1:nmat_block*nblock,1:nmat_block*nblock,&
       &1:ndim,1:nsite_local*nsublat)
  !******************************************************
  double complex, allocatable :: xmat_intermediate(:,:,:,:)
  integer isite,idim
  integer imat,jmat
  integer isublat,ksublat,iblock,jblock,kblock,lblock
  !***** for MPI *****
  integer send_rank,receive_rank,ireq,ierr,tag
  integer status(MPI_STATUS_SIZE)
  double complex mat_rcv(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)

  call who_am_i(myrank,isublat,iblock,jblock)

  if(iformat.EQ.1)then
     if((iblock.EQ.1).AND.(jblock.EQ.1))then
        allocate(xmat_intermediate(1:nmat_block*nblock,1:nmat_block*nblock,&
             &1:ndim,1:nsite_local))
        do isite=1,nsite_local
           do idim=1,ndim
!$omp parallel
!$omp do
              do imat=1,nmat_block
                 do jmat=1,nmat_block     
                    xmat_intermediate(imat,jmat,idim,isite)&
                         &=xmat(imat,jmat,idim,isite)
                 end do
              end do
!$omp end do
!$omp end parallel
           end do
        end do
     end if
     do kblock=1,nblock
        do lblock=1,nblock
           if((kblock.NE.1).OR.(lblock.NE.1))then
              send_rank=(isublat-1)*nblock*nblock
              receive_rank=(isublat-1)*nblock*nblock+(kblock-1)*nblock+lblock-1
              tag=1
              if((kblock.EQ.iblock).AND.(lblock.EQ.jblock))then              
               !  call MPI_Isend(xmat(1,1,1,1),&
                 call MPI_Send(xmat(1,1,1,1),&
                      &nmat_block*nmat_block*nsite_local*ndim,&
                      &MPI_DOUBLE_COMPLEX,&
                      &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
              end if
              if((iblock.EQ.1).AND.(jblock.EQ.1))then
                 call MPI_Recv(mat_rcv(1,1,1,1),&
                      &nmat_block*nmat_block*nsite_local*ndim,&
                      &MPI_DOUBLE_COMPLEX,&
                      &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
                 do isite=1,nsite_local
                    do idim=1,ndim
!$omp parallel
!$omp do
                       do imat=1,nmat_block
                          do jmat=1,nmat_block     
                             xmat_intermediate(imat+(kblock-1)*nmat_block,&
                                  jmat+(lblock-1)*nmat_block,&
                                  idim,isite)&
                                  &=mat_rcv(imat,jmat,idim,isite)
  
                          end do
                       end do
!$omp end do
!$omp end parallel
                    end do
                 end do
              end if
           end if
           !  call MPI_Wait(ireq,status,ierr)
        end do
     end do
     
     
     
     if((iblock.EQ.1).AND.(jblock.EQ.1).AND.(isublat.EQ.1))then!i.e. myrank=0
        do isite=1,nsite_local
           do idim=1,ndim
!$omp parallel
!$omp do
              do imat=1,nmat_block*nblock
                 do jmat=1,nmat_block*nblock  
                    xmat_big(imat,jmat,idim,isite)=&
                         &xmat_intermediate(imat,jmat,idim,isite)
                 end do
              end do
!$omp end do
!$omp end parallel
           end do
        end do
     end if
     
     do ksublat=2,nsublat

        send_rank=0
        receive_rank=(ksublat-1)*nblock*nblock
        tag=ksublat
        if((iblock.EQ.1).AND.(jblock.EQ.1).AND.(isublat.EQ.ksublat))then
           !call MPI_Isend(xmat_intermediate(1,1,1,1),&
           call MPI_Send(xmat_intermediate(1,1,1,1),&
                &nmat_block*nblock*nmat_block*nblock*nsite_local*ndim,&
                &MPI_DOUBLE_COMPLEX,&
                &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
        end if
        if((iblock.EQ.1).AND.(jblock.EQ.1).AND.(isublat.EQ.1))then!i.e. myrank=0
           call MPI_Recv(xmat_big(1,1,1,(ksublat-1)*nsite_local+1),&
                &nmat_block*nblock*nmat_block*nblock*nsite_local*ndim,&
                &MPI_DOUBLE_COMPLEX,&
                &receive_rank,tag,MPI_COMM_WORLD,status,ierr)

        end if
        !call MPI_Wait(ireq,status,ierr)
     end do

     if((iblock.EQ.1).AND.(jblock.EQ.1))then
        deallocate(xmat_intermediate)
     end if

  else if(iformat.EQ.2)then

     call MPI_Bcast(xmat_big(1,1,1,1),&
          nmat_block*nblock*nmat_block*nblock*nsite_local*ndim*nsublat,&
          MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)

     do isite=1,nsite_local
        do idim=1,ndim
!$omp parallel
!$omp do
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 xmat(imat,jmat,idim,isite)=&
                      &xmat_big(imat+(iblock-1)*nmat_block,&
                      jmat+(jblock-1)*nmat_block,&
                      idim,isite+(isublat-1)*nsite_local)
              end do
           end do
!$omp end do
!$omp end parallel
        end do
     end do
     
     
  end if
  
  return

END subroutine matrix_format_change
