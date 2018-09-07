!nxp=1 -> xmat to xmat_mom
!nxp=2 -> xmat_mom to xmat
subroutine Fourier_transform_xmat(xmat,xmat_mom,myrank,nxp)

  implicit none

  include 'mpif.h'
  include 'size_parallel.h'
  !***** input *****
  integer nxp
  integer myrank
  !***** input when nxp=1, output when nxp=2 *****
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,-(nmargin-1):nsite_local+nmargin)
  !***** output when nxp=1, input when nxp=2 *****
  double complex xmat_mom(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)
  !***********************************************
  double complex, allocatable :: xmat_all(:,:,:,:,:)
  double complex, allocatable :: xmat_mom_all(:,:,:,:,:)
  double complex phase
  double precision temp,pi
  integer isite,idim,imom
  integer imat,jmat
  integer isublat,ksublat
  !***** for MPI *****
  integer send_rank,receive_rank,ierr,tag,ireq
  integer status(MPI_STATUS_SIZE)
  integer iblock,jblock
  integer isublat_send,isublat_rcv,ishift

  pi=2d0*dasin(1d0)

  if (nxp.EQ.1)then

     allocate(xmat_all(1:nmat_block,1:nmat_block,1:ndim,&
          &1:nsite_local,1:nsublat))

     call who_am_i(myrank,isublat,iblock,jblock)
     
     do isite=1,nsite_local
        do idim=1,ndim
!$omp parallel
!$omp do
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 xmat_all(imat,jmat,idim,isite,isublat)&
                      &=xmat(imat,jmat,idim,isite)
              end do
           end do
!$omp end do
!$omp end parallel
        end do
     end do
     isublat_send=isublat
     isublat_rcv=isublat
     do ishift=1,nsublat-1        
        isublat_send=isublat_send+1
        isublat_rcv=isublat_rcv-1
        if(isublat_send.EQ.nsublat+1)then
           isublat_send=1
        end if
        if(isublat_rcv.EQ.0)then
           isublat_rcv=nsublat
        end if
        send_rank=(isublat_send-1)*nblock*nblock+(iblock-1)*nblock+jblock-1
        receive_rank=(isublat_rcv-1)*nblock*nblock+(iblock-1)*nblock+jblock-1
        
        tag=1
        call MPI_Isend(xmat(1,1,1,1),&
             &nmat_block*nmat_block*nsite_local*ndim,&
             &MPI_DOUBLE_COMPLEX,&
             &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
        call MPI_Recv(xmat_all(1,1,1,1,isublat_rcv),&
             &nmat_block*nmat_block*nsite_local*ndim,&
             &MPI_DOUBLE_COMPLEX,&
             &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
        call MPI_Wait(ireq,status,ierr)
        
     end do

     xmat_mom=(0d0,0d0)
     do imom=1,nsite_local
        do isite=1,nsite_local
           do ksublat=1,nsublat
              !calculate the phase.
              temp=2d0*pi*dble((imom+(isublat-1)*nsite_local)&
                   &*(isite+(ksublat-1)*nsite_local))&
                   &/dble(nsublat*nsite_local)
              phase=dcmplx(dcos(temp))-(0d0,1d0)*dcmplx(dsin(temp))
              phase=phase/dcmplx(dsqrt(dble(nsublat*nsite_local)))
              do idim=1,ndim
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       xmat_mom(imat,jmat,idim,imom)=&
                            &xmat_mom(imat,jmat,idim,imom)&
                            &+xmat_all(imat,jmat,idim,isite,ksublat)*phase
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
           end do
        end do
     end do
     
     deallocate(xmat_all)

  else if(nxp.EQ.2)then

     allocate(xmat_mom_all(1:nmat_block,1:nmat_block,1:ndim,&
          &1:nsite_local,1:nsublat))

     call who_am_i(myrank,isublat,iblock,jblock)

     do imom=1,nsite_local
        do idim=1,ndim
!$omp parallel 
!$omp do
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 xmat_mom_all(imat,jmat,idim,imom,isublat)&
                      &=xmat_mom(imat,jmat,idim,imom)
              end do
           end do
!$omp end do
!$omp end parallel
        end do
     end do
     isublat_send=isublat
     isublat_rcv=isublat
     do ishift=1,nsublat-1

        isublat_send=isublat_send+1
        isublat_rcv=isublat_rcv-1
        if(isublat_send.EQ.nsublat+1)then
           isublat_send=1
        end if
        if(isublat_rcv.EQ.0)then
           isublat_rcv=nsublat
        end if
        send_rank=(isublat_send-1)*nblock*nblock+(iblock-1)*nblock+jblock-1
        receive_rank=(isublat_rcv-1)*nblock*nblock+(iblock-1)*nblock+jblock-1

        tag=1
        call MPI_Isend(xmat_mom(1,1,1,1),&
             &nmat_block*nmat_block*nsite_local*ndim,&
             &MPI_DOUBLE_COMPLEX,&
             &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
        call MPI_Recv(xmat_mom_all(1,1,1,1,isublat_rcv),&
             &nmat_block*nmat_block*nsite_local*ndim,&
             &MPI_DOUBLE_COMPLEX,&
             &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
        call MPI_Wait(ireq,status,ierr)

     end do
     xmat=(0d0,0d0)
     do imom=1,nsite_local
        do isite=1,nsite_local
           do ksublat=1,nsublat
              !calculate the phase.
              temp=2d0*pi*dble((isite+(isublat-1)*nsite_local)&
                   &*(imom+(ksublat-1)*nsite_local))&
                   &/dble(nsublat*nsite_local)
              phase=dcmplx(dcos(temp))+(0d0,1d0)*dcmplx(dsin(temp))
              phase=phase/dcmplx(dsqrt(dble(nsublat*nsite_local)))
              do idim=1,ndim
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       xmat(imat,jmat,idim,isite)=&
                            &xmat(imat,jmat,idim,isite)&
                            &+xmat_mom_all(imat,jmat,idim,imom,ksublat)*phase
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
           end do
        end do
     end do

     deallocate(xmat_mom_all)

  end if

  return

END subroutine Fourier_transform_xmat
