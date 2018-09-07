SUBROUTINE Multiply_Dirac_to_chi(temperature,xmat_row,xmat_column,alpha,&
chi,mchi,GAMMA10d,nbmn,flux,myrank)

  implicit none
  include 'size_parallel.h'
  !***** input *****
  integer nbmn
  integer myrank 
  double complex xmat_row(1:nmat_block,1:nmat_block*nblock,&
       &1:ndim,1:nsite_local)
  double complex xmat_column(1:nmat_block*nblock,1:nmat_block,&
       &1:ndim,1:nsite_local)
  double precision alpha(1:nmat_block*nblock)
  double precision temperature,flux
  double complex chi(1:nremez_md,1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin)
  double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
  !***** output *****
  double complex Mchi(1:nremez_md,1:nmat_block,1:nmat_block,&
       &1:nspin,1:nsite_local)
  !******************
  double complex pf1(1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin)
  double complex pf2(1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin)
  integer imat,jmat
  integer ispin
  integer isite
  integer iremez
 
  do iremez=1,nremez_md

!$omp parallel
!$omp do  
     do imat=1,nmat_block
        do jmat=1,nmat_block
           do ispin=1,nspin
              do isite=-nimprove,nsite_local+1+nimprove
                 pf1(imat,jmat,ispin,isite)=Chi(iremez,imat,jmat,ispin,isite)
              end do
           end do
        end do
     end do
!$omp end do
!$omp end parallel
     call Multiply_Dirac(temperature,xmat_row,xmat_column,&
          &alpha,pf1,pf2,GAMMA10d,nbmn,flux,myrank)
!$omp parallel
!$omp do
     do imat=1,nmat_block
        do jmat=1,nmat_block
           do ispin=1,nspin
              do isite=1,nsite_local
                MChi(iremez,imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)
             end do
          end do
        end do
     end do
!$omp end do
!$omp end parallel
     
  end do
  

  return

END SUBROUTINE Multiply_Dirac_to_chi
