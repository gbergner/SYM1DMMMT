! Calculate the largest eigenvalue of (D^dag*D),  
! by multiplying (D^dag*D) many times to a random vector. 
SUBROUTINE Largest_eigenvalue(temperature,&
     xmat,alpha,GAMMA10d,neig,largest_eig,myrank,nbc,nbmn,flux)

  implicit none

  include 'mpif.h'
  include 'size_parallel.h'
  !***** input *****
  integer neig,nbc,nbmn,myrank
  double precision temperature,flux
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  double precision alpha(1:nmat_block*nblock)
  double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
  !***** output *****
  double precision largest_eig
  !******************
  double complex phi1(1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin)
  double complex phi2(1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin)
  !double complex trace(1:nspin*nsite_local*nsublat)
  !double complex trace_local(1:nspin*nsite_local*nsublat)
  double complex trace(1:nsublat,1:nspin,1:nsite_local)
  double complex trace_local(1:nsublat,1:nspin,1:nsite_local)
  double precision norm,norm_local
  integer imat,jmat
  integer ispin
  integer isite
  integer iblock,jblock
  integer isublat
  integer ieig
  integer i
  double complex xmat_row(1:nmat_block,1:nmat_block*nblock,1:ndim,1:nsite_local)
  double complex xmat_column(1:nmat_block*nblock,1:nmat_block,1:ndim,&
       &1:nsite_local)
  double precision r1,r2
  !***** For MPI *****  
  integer IERR

  call who_am_i(myrank,isublat,iblock,jblock)
  !move i-th row and j-th row of xmat to (i,j)-th node.
  call mpi_xmat_row(xmat,xmat_row,myrank)
  call mpi_xmat_column(xmat,xmat_column,myrank)
  !***********************************
  !**** generate a random vector. ****
  !***********************************

  if(myrank.GT.0)then
     !throw away some random numbers in order to avoid the communication
     do i=1,nmat_block*nmat_block*nspin*nsite_local*myrank
       call BoxMuller(r1,r2)
     end do
  end if

  norm_local=0d0
  do imat=1,nmat_block     
     do jmat=1,nmat_block
        do ispin=1,nspin
           do isite=1,nsite_local
              call BoxMuller(r1,r2)
              phi1(imat,jmat,ispin,isite)=&
                   &(dcmplx(r1)+dcmplx(r2)*(0D0,1D0))/dcmplx(dsqrt(2d0))
              norm_local=norm_local+(r1*r1+r2*r2)*0.5d0
           end do
        end do
     end do
  end do

  if(myrank.LT.nsublat*nblock*nblock-1)then
     !throw away some random numbers in order to avoid the communication
     do i=1,nmat_block*nmat_block*nspin*nsite_local*&
          (nsublat*nblock*nblock-1-myrank)
       call BoxMuller(r1,r2)
     end do
  end if
  !write(*,*)r1
  !****************************
  !*** traceless projection ***
  !****************************
  call who_am_i(myrank,isublat,iblock,jblock)
  trace_local=(0d0,0d0)
  if(iblock.EQ.jblock)then
     do ispin=1,nspin
        do isite=1,nsite_local
           do imat=1,nmat_block
              !trace_local((isublat-1)*nspin*nsite_local+&
              !     &(ispin-1)*nsite_local+isite)=&
             !      &trace_local((isublat-1)*nspin*nsite_local+&
             !      &(ispin-1)*nsite_local+isite)&
              !      &+phi1(imat,imat,ispin,isite)
              trace_local(isublat,ispin,isite)=trace_local(isublat,ispin,isite)&
                   &+phi1(imat,imat,ispin,isite)
           end do
        end do
     end do
  end if

    
  call MPI_Allreduce(trace_local,trace,nsublat*nspin*nsite_local,&
       MPI_DOUBLE_COMPLEX,&
       &MPI_SUM,MPI_COMM_WORLD,IERR)


     trace=trace/dcmplx(nmat_block*nblock)
 
  call who_am_i(myrank,isublat,iblock,jblock)
  if(iblock.EQ.jblock)then
     do ispin=1,nspin
        do isite=1,nsite_local
           do imat=1,nmat_block
              phi1(imat,imat,ispin,isite)=phi1(imat,imat,ispin,isite)&
                &-trace(isublat,ispin,isite)
           end do         
        end do
     end do
  end if
  !******************************
  !*** adjust margin and b.c. ***
  !******************************
  call Adjust_margin_and_bc_pf(phi1,myrank,nbc)
  !******************************************
  !*** random vector has been generated.  ***
  !******************************************

  do ieig=1,neig
     !********************************************
     !*** phi1 -> phi2=D*phi -> phi1=D^dag*phi2***
     !********************************************
     call Multiply_Dirac(temperature,xmat_row,xmat_column,&
          &alpha,phi1,phi2,GAMMA10d,nbmn,flux,myrank)
     !******************************
     !*** adjust margin and b.c. ***
     !******************************
     call Adjust_margin_and_bc_pf(phi2,myrank,nbc)
     
     call Multiply_Dirac_dagger(temperature,xmat_row,xmat_column,&
          &alpha,phi2,phi1,GAMMA10d,nbmn,flux,myrank)
     !******************************
     !*** adjust margin and b.c. ***
     !******************************
     call Adjust_margin_and_bc_pf(phi1,myrank,nbc)


     !************************************
     !*** calculate the norm of phi1.  ***
     !************************************
     norm_local=0d0
     do imat=1,nmat_block     
        do jmat=1,nmat_block
           do ispin=1,nspin
              do isite=1,nsite_local
                 norm_local=norm_local&
                      &+dble(phi1(imat,jmat,ispin,isite)&
                      &*dconjg(phi1(imat,jmat,ispin,isite)))
              end do
           end do
        end do
     end do
!     call MPI_Reduce(norm_local,norm,1,MPI_DOUBLE_PRECISION,&
!          &MPI_SUM,0,MPI_COMM_WORLD,IERR)
!     call MPI_Bcast(norm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
 call MPI_Allreduce(norm_local,norm,1,MPI_DOUBLE_PRECISION,&
      &MPI_SUM,MPI_COMM_WORLD,IERR)
     norm=dsqrt(norm)
     norm=1d0/norm
     phi1=phi1*norm
     norm=1d0/norm

     !write(*,*)norm
  end do
  
  largest_eig=norm

  return
  
END SUBROUTINE Largest_eigenvalue
