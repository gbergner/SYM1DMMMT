!Derivative of (M*Chi)(t)_{ij} w.r.t. alpha_{k}  
SUBROUTINE Derivative_Dirac_alpha(alpha,chi,Deriv_Mchi_alpha,myrank)

  implicit none

  include 'size_parallel.h'
  double precision alpha(1:nmat_block*nblock)
  double complex chi(1:nremez_md,1:nmat_block,1:nmat_block,1:nspin,-(nmargin-1):nsite_local+nmargin)
  double complex phase,phase2
  double precision aij
  double complex Deriv_Mchi_alpha(1:nremez_md,1:nmat_block*nblock,&
       &1:nmat_block,1:nmat_block,1:nspin,1:nsite_local)

  integer imat,jmat,kmat
  integer ispin
  integer isite
  integer iremez

  integer iblock,jblock,isublat,myrank
  call who_am_i(myrank,isublat,iblock,jblock)
  !***************************
  !***  kinetic part only  ***
  !***************************
  Deriv_Mchi_alpha=(0d0,0d0)
  if(nimprove.EQ.0)then
     do isite=1,nsite_local
!$omp parallel
!$omp do
        do imat=1,nmat_block
           do jmat=1,nmat_block
              
              aij=alpha(imat+(iblock-1)*nmat_block)&
                   &-alpha(jmat+(jblock-1)*nmat_block)
              aij=aij/dble(nsite_local*nsublat)
              phase=dcmplx(dcos(aij))+(0d0,1d0)*dcmplx(dsin(aij))
              
              do iremez=1,nremez_md
                 do ispin=1,8!nspin
                    
                    kmat=jmat+(jblock-1)*nmat_block
                    Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite)=&
                         Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite)&
                         +chi(iremez,imat,jmat,ispin,isite+1)*phase&
                         /dcmplx(nsite_local*nsublat)*(0d0,-1d0)
                    
                    kmat=imat+(iblock-1)*nmat_block
                    Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite)=&
                         Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite)&
                         -chi(iremez,imat,jmat,ispin,isite+1)*phase&
                         /dcmplx(nsite_local*nsublat)*(0d0,-1d0)
                    
                    
                    kmat=jmat+(jblock-1)*nmat_block
                    Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite)=&
                         Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite)&
                         +chi(iremez,imat,jmat,ispin+8,isite-1)*dconjg(phase)&
                         /dcmplx(nsite_local*nsublat)*(0d0,-1d0)
                    
                    kmat=imat+(iblock-1)*nmat_block
                    Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite)=&
                         Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite)&
                         -chi(iremez,imat,jmat,ispin+8,isite-1)*dconjg(phase)&
                         /dcmplx(nsite_local*nsublat)*(0d0,-1d0)
                    
                 end do
              end do
           end do
        end do
!$omp end do
!$omp end parallel
     end do
     
  else if(nimprove.EQ.1)then
     do isite=1,nsite_local
!$omp parallel
!$omp do
        do imat=1,nmat_block
           do jmat=1,nmat_block
              
              aij=alpha(imat+(iblock-1)*nmat_block)&
                   &-alpha(jmat+(jblock-1)*nmat_block)
              aij=aij/dble(nsite_local*nsublat)
              phase=dcmplx(dcos(aij))+(0d0,1d0)*dcmplx(dsin(aij))
              phase2=phase*phase
              do iremez=1,nremez_md
                 do ispin=1,8!nspin
                    
                    kmat=jmat+(jblock-1)*nmat_block
                    Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite)=&
                         Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite)&
                         +chi(iremez,imat,jmat,ispin,isite+1)*phase&
                         /dcmplx(nsite_local*nsublat)*(0d0,-2d0)&
                         +chi(iremez,imat,jmat,ispin,isite+2)*phase2&
                         /dcmplx(nsite_local*nsublat)*(0d0,1d0)                    

                    kmat=imat+(iblock-1)*nmat_block
                    Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite)=&
                         Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite)&
                         -chi(iremez,imat,jmat,ispin,isite+1)*phase&
                         /dcmplx(nsite_local*nsublat)*(0d0,-2d0)&
                         -chi(iremez,imat,jmat,ispin,isite+2)*phase2&
                         /dcmplx(nsite_local*nsublat)*(0d0,1d0)
                    
                    
                    kmat=jmat+(jblock-1)*nmat_block
                    Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite)=&
                         Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite)&
                         +chi(iremez,imat,jmat,ispin+8,isite-1)*dconjg(phase)&
                         /dcmplx(nsite_local*nsublat)*(0d0,-2d0)&
                         +chi(iremez,imat,jmat,ispin+8,isite-2)*dconjg(phase2)&
                         /dcmplx(nsite_local*nsublat)*(0d0,1d0)                    

                    kmat=imat+(iblock-1)*nmat_block
                    Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite)=&
                         Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite)&
                         -chi(iremez,imat,jmat,ispin+8,isite-1)*dconjg(phase)&
                         /dcmplx(nsite_local*nsublat)*(0d0,-2d0)&
                         -chi(iremez,imat,jmat,ispin+8,isite-2)*dconjg(phase2)&
                         /dcmplx(nsite_local*nsublat)*(0d0,1d0)
                    
                 end do
              end do
           end do
        end do
!$omp end do
!$omp end parallel
     end do
  end if

  return

END SUBROUTINE Derivative_Dirac_alpha
