!trace part of alpha is removed. 
!trace part of (¥int dt X) is removed. 
SUBROUTINE subtract_U1(xmat,alpha,myrank)
  
  implicit none
  include 'mpif.h'
  include 'size_parallel.h'
  integer IERR,myrank
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,-(nmargin-1):nsite_local+nmargin)
  double precision alpha(1:nmat_block*nblock)
  integer idim,imat,isite
  double complex trace(1:ndim),trace_local(1:ndim)
  double precision sum_alpha
  integer iblock,jblock
  integer isublat

  call who_am_i(myrank,isublat,iblock,jblock)
  trace_local=(0d0,0d0)
  trace=(0d0,0d0)
  if(iblock.EQ.jblock)then
     do idim=1,ndim
        do isite=1,nsite_local
           do imat=1,nmat_block
              trace_local(idim)=trace_local(idim)+xmat(imat,imat,idim,isite)
           end do
        end do
     end do
  end if
  !call MPI_Reduce(trace_local(1),trace(1),ndim,MPI_DOUBLE_COMPLEX,&
  !     MPI_SUM,0,MPI_COMM_WORLD,IERR)
  !if(myrank.EQ.0)then
  !   trace=trace/dcmplx(nsite_local*nsublat*nmat_block*nblock)
  !end if
  !call MPI_Bcast(trace(1),ndim,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
  call MPI_Allreduce(trace_local(1),trace(1),ndim,MPI_DOUBLE_COMPLEX,&
       MPI_SUM,MPI_COMM_WORLD,IERR)
  trace=trace/dcmplx(nsite_local*nsublat*nmat_block*nblock)

  if(iblock.EQ.jblock)then
     do idim=1,ndim
        !take care of the margin too. 
        do isite=-(nmargin-1),nsite_local+nmargin
!$omp parallel
!$omp do
           do imat=1,nmat_block
              xmat(imat,imat,idim,isite)=xmat(imat,imat,idim,isite)-trace(idim)
           end do
!$omp end do
!$omp end parallel
        end do
     end do
  end if
  sum_alpha=0d0
  do imat=1,nmat_block*nblock
     sum_alpha=sum_alpha+alpha(imat)
  end do
  sum_alpha=sum_alpha/dble(nmat_block*nblock)
!$omp parallel
!$omp do
  do imat=1,nmat_block*nblock
     alpha(imat)=alpha(imat)-sum_alpha
  end do
!$omp end do
!$omp end parallel
  !  end if
  !  call MPI_Bcast(alpha(1),nmat,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  
  return

END SUBROUTINE subtract_U1
