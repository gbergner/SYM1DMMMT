subroutine Calc_TrX2(xmat,trx2,myrank)

  implicit none

  include 'mpif.h'
  include 'size_parallel.h'
  integer myrank
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,-(nmargin-1):nsite_local+nmargin)
  double precision trx2,trx2_local

  integer isite,idim
  integer imat,jmat
  integer IERR

  trx2_local=0d0
  do isite=1,nsite_local
     !isite=0, nsite_local+1 are considered at neighboring nodes.
     do idim=1,ndim
        do imat=1,nmat_block
           do jmat=1,nmat_block
              trx2_local=trx2_local&
                   +dble(xmat(imat,jmat,idim,isite)&
                   *dconjg(xmat(imat,jmat,idim,isite)))
           end do
        end do
     end do
  end do

  call MPI_Reduce(trx2_local,trx2,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,0,MPI_COMM_WORLD,IERR)
  if(myrank.EQ.0)then
     trx2=trx2/dble(nmat_block*nblock*nsite_local*nsublat)
  end if


  return

END subroutine Calc_TrX2
!##################################################
subroutine Calc_TrX2_each(xmat,sum_trx2,myrank,trx2)

  implicit none

  include 'mpif.h'
  include 'size_parallel.h'
  integer myrank
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,-(nmargin-1):nsite_local+nmargin)
  double precision sum_trx2,trx2(1:ndim),trx2_local(1:ndim),trx2_global

  integer isite,idim
  integer imat,jmat
  integer IERR
  trx2_local=0d0
  do idim=1,ndim
     do isite=1,nsite_local
        !isite=0, nsite_local+1 are considered at neighboring nodes.
        do imat=1,nmat_block
           do jmat=1,nmat_block
              trx2_local(idim)=trx2_local(idim)&
                   +dble(xmat(imat,jmat,idim,isite)&
                   *dconjg(xmat(imat,jmat,idim,isite)))
           end do
        end do
     end do
     
  end do

  call MPI_Reduce(trx2_local,trx2,ndim,MPI_DOUBLE_PRECISION,&
          MPI_SUM,0,MPI_COMM_WORLD,IERR)
  if(myrank.EQ.0)then
     trx2=trx2/dble(nmat_block*nblock*nsite_local*nsublat)
  end if
  sum_trx2=0d0
  do idim=1,ndim
     sum_trx2=sum_trx2+trx2(idim)
  end do
  
  return

END subroutine Calc_TrX2_Each
