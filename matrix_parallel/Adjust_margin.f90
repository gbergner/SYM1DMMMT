!This subroutine adjust the margins of xmat (isite < 1 and isite > nsite_local),
!by communicating with the neigbouring processes.
SUBROUTINE Adjust_margin_xmat(xmat,myrank)

  implicit none
  include 'size_parallel.h'
  include 'mpif.h'
  !***** input *****
  integer myrank
  !***** input & output *****
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  integer iblock,jblock,isublat
  !***** for MPI *****
  integer send_rank,receive_rank,ireq,ierr,tag
  integer status(MPI_STATUS_SIZE)

  call who_am_i(myrank,isublat,iblock,jblock)
  if(isublat.ne.nsublat)then
     send_rank=myrank+nblock*nblock
  else
     send_rank=(iblock-1)*nblock+jblock-1
  end if
  if(isublat.ne.1)then
     receive_rank=myrank-nblock*nblock
  else
     receive_rank=(nsublat-1)*nblock*nblock+(iblock-1)*nblock+jblock-1
  end if
  tag=1
  call MPI_Isend(xmat(1,1,1,nsite_local-(nmargin-1)),&
       &nmat_block*nmat_block*ndim*nmargin,&
       &MPI_DOUBLE_COMPLEX,&
       &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
  call MPI_Recv(xmat(1,1,1,-(nmargin-1)),nmat_block*nmat_block*ndim*nmargin,&
       &MPI_DOUBLE_COMPLEX,&
       &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
  call MPI_Wait(ireq,status,ierr)
  !******
  if(isublat.ne.1)then
     send_rank=myrank-nblock*nblock
  else
     send_rank=(nsublat-1)*nblock*nblock+(iblock-1)*nblock+jblock-1
  end if
  if(isublat.ne.nsublat)then
     receive_rank=myrank+nblock*nblock
  else
     receive_rank=(iblock-1)*nblock+jblock-1
  end if
  tag=3
  call MPI_Isend(xmat(1,1,1,1),nmat_block*nmat_block*ndim*nmargin,&
       &MPI_DOUBLE_COMPLEX,&
       &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
  call MPI_Recv(xmat(1,1,1,nsite_local+1),nmat_block*nmat_block*ndim*nmargin,&
       &MPI_DOUBLE_COMPLEX,&
       &receive_rank,tag,MPI_COMM_WORLD,status,ierr) 
  call MPI_Wait(ireq,status,ierr)

  return

END SUBROUTINE Adjust_margin_xmat
!***************************************************************
!This subroutine adjust the margin of pseudo-fermion 
! (isite < 1 and isite > nsite_local), 
!by communicating with the neigbouring processes.
!It also properly sets the boundary condition (pbc or apbc).
SUBROUTINE Adjust_margin_and_bc_pf(pf,myrank,nbc)

  implicit none
  include 'size_parallel.h'
  include 'mpif.h'
  !***** input *****
  integer myrank,nbc
  !***** input & output *****
  double complex pf(1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin)
  !**************************
  integer isite,imat,jmat,ispin,isublat,iblock,jblock
  !***** for MPI *****
  integer send_rank,receive_rank,ireq,&
       ierr,tag
  integer status(MPI_STATUS_SIZE)

  call who_am_i(myrank,isublat,iblock,jblock)
  !*********************
  !*** adjust margin ***
  !*********************
  if(isublat.ne.nsublat)then
     send_rank=myrank+nblock*nblock
  else
     send_rank=(iblock-1)*nblock+jblock-1
  end if
  if(isublat.ne.1)then
     receive_rank=myrank-nblock*nblock
  else
     receive_rank=(nsublat-1)*nblock*nblock+(iblock-1)*nblock+jblock-1
  end if
  tag=1
  call MPI_Isend(pf(1,1,1,nsite_local-(nmargin-1)),&
       &nmat_block*nmat_block*nspin*nmargin,&
       &MPI_DOUBLE_COMPLEX,&
       &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
  call MPI_Recv(pf(1,1,1,-(nmargin-1)),nmat_block*nmat_block*nspin*nmargin,&
       &MPI_DOUBLE_COMPLEX,&
       &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
  call MPI_Wait(ireq,status,ierr)
  !*****
  if(isublat.ne.1)then
     send_rank=myrank-nblock*nblock
  else
     send_rank=(nsublat-1)*nblock*nblock+(iblock-1)*nblock+jblock-1
  end if
  if(isublat.ne.nsublat)then
     receive_rank=myrank+nblock*nblock
  else
     receive_rank=(iblock-1)*nblock+jblock-1
  end if
  tag=3
  call MPI_Isend(pf(1,1,1,1),nmat_block*nmat_block*nspin*nmargin,&
       &MPI_DOUBLE_COMPLEX,&
       &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
  call MPI_Recv(pf(1,1,1,nsite_local+1),nmat_block*nmat_block*nspin*nmargin,&
       &MPI_DOUBLE_COMPLEX,&
       &receive_rank,tag,MPI_COMM_WORLD,status,ierr) 
  call MPI_Wait(ireq,status,ierr)
  !***************************
  !*** boundary condition  ***
  !***************************
  if(nbc.EQ.0)then!pbc
     !don't do anything. 
  else if(nbc.EQ.1)then!apbc
     if(isublat.EQ.1)then
!$omp parallel
!$omp do
        do imat=1,nmat_block     
           do jmat=1,nmat_block
              do ispin=1,nspin
                 do isite=-(nmargin-1),0
                    pf(imat,jmat,ispin,isite)=&
                         &pf(imat,jmat,ispin,isite)*(-1d0,0d0)
                 end do
              end do
           end do
        end do
!$omp end do
!$omp end parallel
     else if(isublat.EQ.nsublat)then
!$omp parallel
!$omp do
        do imat=1,nmat_block     
           do jmat=1,nmat_block
              do ispin=1,nspin
                 do isite=nsite_local+1,nsite_local+nmargin
                    pf(imat,jmat,ispin,isite)=&
                         &pf(imat,jmat,ispin,isite)*(-1d0,0d0)
                 end do
              end do
           end do
        end do
!$omp end do
!$omp end parallel
     end if
  end if

  return

END SUBROUTINE Adjust_margin_and_bc_pf
 
!***************************************************************
!This subroutine adjust the margin of chi (isite < 1 and isite > nsite_local), 
!by communicating with the neigbouring processes.
!It also properly sets the boundary condition (pbc or apbc).
SUBROUTINE Adjust_margin_and_bc_Chi(Chi,myrank,nbc,nremez)

  implicit none
  include 'size_parallel.h'
  include 'mpif.h'
  !***** input *****
  integer myrank,nbc,nremez
  !***** input & output *****
  double complex Chi(1:nremez,1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin)
  !**************************
  integer isite,imat,jmat,ispin,iremez,isublat,iblock,jblock
  !***** for MPI *****
  integer send_rank,receive_rank,ireq,&
       ierr,tag
  integer status(MPI_STATUS_SIZE)

  call who_am_i(myrank,isublat,iblock,jblock)
  !*************************
  !*** adjust the margin ***
  !************************* 
   if(isublat.ne.nsublat)then
      send_rank=myrank+nblock*nblock
   else
      send_rank=(iblock-1)*nblock+jblock-1
   end if
   if(isublat.ne.1)then
      receive_rank=myrank-nblock*nblock
   else
      receive_rank=(nsublat-1)*nblock*nblock+(iblock-1)*nblock+jblock-1
   end if
  tag=1
  call MPI_Isend(Chi(1,1,1,1,nsite_local-(nmargin-1)),&
       &nremez*nmat_block*nmat_block*nspin*nmargin,&
       &MPI_DOUBLE_COMPLEX,&
       &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
  call MPI_Recv(Chi(1,1,1,1,-(nmargin-1)),&
       &nremez*nmat_block*nmat_block*nspin*nmargin,&
       &MPI_DOUBLE_COMPLEX,&
       &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
  call MPI_Wait(ireq,status,ierr)
  !******
  if(isublat.ne.1)then
     send_rank=myrank-nblock*nblock
  else
     send_rank=(nsublat-1)*nblock*nblock+(iblock-1)*nblock+jblock-1
  end if
  if(isublat.ne.nsublat)then
     receive_rank=myrank+nblock*nblock
  else
     receive_rank=(iblock-1)*nblock+jblock-1
  end if
  tag=3
  call MPI_Isend(Chi(1,1,1,1,1),nremez*nmat_block*nmat_block*nspin*nmargin,&
       &MPI_DOUBLE_COMPLEX,&
       &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
  call MPI_Recv(Chi(1,1,1,1,nsite_local+1),&
       &nremez*nmat_block*nmat_block*nspin*nmargin,&
       &MPI_DOUBLE_COMPLEX,&
       &receive_rank,tag,MPI_COMM_WORLD,status,ierr) 
  call MPI_Wait(ireq,status,ierr)
  !**************************
  !*** boundary condition *** 
  !**************************
  if(nbc.EQ.0)then!pbc
     !don't do anything. 
  else if(nbc.EQ.1)then!apbc
     if(isublat.EQ.1)then
        do ispin=1,nspin
!$omp parallel 
!$omp do
           do jmat=1,nmat_block
              do imat=1,nmat_block  
                 do iremez=1,nremez
                    
                    do isite=-(nmargin-1),0
                       Chi(iremez,imat,jmat,ispin,isite)=&
                            Chi(iremez,imat,jmat,ispin,isite)*(-1d0,0d0)
                    end do
                 end do
              end do
           end do
!$omp end do
!$omp end parallel
        end do
     else if(isublat.EQ.nsublat)then
        do ispin=1,nspin
!$omp parallel 
!$omp do
           do jmat=1,nmat_block
              do imat=1,nmat_block  
                 do iremez=1,nremez
                    do isite=nsite_local+1,nsite_local+nmargin
                       Chi(iremez,imat,jmat,ispin,isite)=&
                            Chi(iremez,imat,jmat,ispin,isite)*(-1d0,0d0)
                    end do
                 end do
              end do
           end do
!$omp end do
!$omp end parallel
        end do
     end if
  end if

  return

END SUBROUTINE Adjust_margin_and_bc_Chi
