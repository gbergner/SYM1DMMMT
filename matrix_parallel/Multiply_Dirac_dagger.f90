! Pf2 = [(-i)*D/N]^dagger * Pf1
SUBROUTINE Multiply_Dirac_dagger(temperature,xmat_row_2,xmat_column_2,&
     &alpha,pf1,pf2,GAMMA10d,nbmn,flux,myrank)

  implicit none

  include 'mpif.h'
  include 'size_parallel.h'
  !***** input *****
  integer nbmn,myrank
  double precision temperature,flux
  double precision alpha(1:nmat_block*nblock)
  double complex pf1(1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin)
  double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
  double complex xmat_row_2(1:nmat_block,1:nmat_block*nblock,1:ndim,&
       &1:nsite_local)
  double complex xmat_column_2(1:nmat_block*nblock,1:nmat_block,1:ndim,&
       &1:nsite_local)
  !***** output *****
  !pf2=M*pf1, M: Dirac op
  double complex pf2(1:nmat_block,1:nmat_block,1:nspin,&
       -(nmargin-1):nsite_local+nmargin)
  !******************
  double complex xmat_row(1:nmat_block,1:nmat_block*nblock,1:ndim,1:nsite_local)
  double complex xmat_column(1:nmat_block*nblock,1:nmat_block,1:ndim,&
       &1:nsite_local)
  double complex Gam123(1:nspin,1:nspin),Gam12(1:nspin,1:nspin)
  double complex phase,phase2,GAM
  double precision aij
  double precision lattice_spacing
  integer imat,jmat,kmat
  integer idim
  integer ispin,jspin,kspin
  integer isite
  integer isublat
  integer iblock,jblock,kblock_send,kblock_rcv
  integer ishift
  !***** For MPI *****
  double complex mat_send(1:nmat_block,1:nmat_block,1:nspin,1:nsite_local)
  double complex mat_rcv(1:nmat_block,1:nmat_block,1:nspin,1:nsite_local)
  integer IERR,IREQ,send_rank,receive_rank,tag
  integer STATUS(MPI_STATUS_SIZE)

  do idim=1,ndim
!$omp parallel
!$omp do
     do imat=1,nmat_block
        do jmat=1,nmat_block*nblock
           do isite=1,nsite_local
              xmat_row(imat,jmat,idim,isite)=&
                   &dconjg(xmat_column_2(jmat,imat,idim,isite))
              xmat_column(jmat,imat,idim,isite)=&
                   &dconjg(xmat_row_2(imat,jmat,idim,isite))
           end do
        end do
     end do
!$omp end do
!$omp end parallel
  end do
  call who_am_i(myrank,isublat,iblock,jblock)
  lattice_spacing=1d0/temperature/dble(nsite_local*nsublat)
  !****************************************************
  !*****************************************************
  !***  kinetic part; no need for MPI communication  ***
  !*****************************************************
  !*****************************************************
  !pf2=(0d0,1d0)*pf1
  pf2=(0d0,0d0)
  if(nimprove.EQ.0)then
     !Naive action.
     do isite=1,nsite_local
!$omp parallel
!$omp do
        do jmat=1,nmat_block
           do imat=1,nmat_block
              aij=alpha(imat+(iblock-1)*nmat_block)&
                   &-alpha(jmat+(jblock-1)*nmat_block)
              aij=aij/dble(nsite_local*nsublat)
              phase=dcmplx(dcos(aij))+(0d0,1d0)*dcmplx(dsin(aij))
              do ispin=1,8!nspin
                 pf2(imat,jmat,ispin+8,isite)=&
                      pf2(imat,jmat,ispin+8,isite)&
                      -(1d0,0d0)*phase*pf1(imat,jmat,ispin,isite+1)&
                      +(1d0,0d0)*pf1(imat,jmat,ispin,isite)
                 
                 pf2(imat,jmat,ispin,isite)=&
                      pf2(imat,jmat,ispin,isite)&
                      +(1d0,0d0)*dconjg(phase)*pf1(imat,jmat,ispin+8,isite-1)&
                      -(1d0,0d0)*pf1(imat,jmat,ispin+8,isite)
              end do
           end do
        end do
!$omp end do
!$omp end parallel
     end do
  else if(nimprove.EQ.1)then
     !Improved action.
     do isite=1,nsite_local
!$omp parallel
!$omp do
        do jmat=1,nmat_block
           do imat=1,nmat_block
              aij=alpha(imat+(iblock-1)*nmat_block)&
                   &-alpha(jmat+(jblock-1)*nmat_block)
              aij=aij/dble(nsite_local*nsublat)
              phase=dcmplx(dcos(aij))+(0d0,1d0)*dcmplx(dsin(aij))
              phase2=phase*phase
              do ispin=1,8!nspin
                 pf2(imat,jmat,ispin+8,isite)=&
                      pf2(imat,jmat,ispin+8,isite)&
                      +(0.5d0,0d0)*phase2*pf1(imat,jmat,ispin,isite+2)&
                      -(2d0,0d0)*phase*pf1(imat,jmat,ispin,isite+1)&
                      +(1.5d0,0d0)*pf1(imat,jmat,ispin,isite)
                 
                 pf2(imat,jmat,ispin,isite)=&
                      pf2(imat,jmat,ispin,isite)&
                      -(0.5d0,0d0)*dconjg(phase2)*pf1(imat,jmat,ispin+8,isite-2)&
                      +(2d0,0d0)*dconjg(phase)*pf1(imat,jmat,ispin+8,isite-1)&
                      -(1.5d0,0d0)*pf1(imat,jmat,ispin+8,isite)
              end do
           end do
        end do
!$omp end do
!$omp end parallel
     end do
  end if
  !*****************************************************
  !*****************************************************
  !***  interaction part; MPI communication needed.  ***
  !*****************************************************
  !*****************************************************
  !###################################################################
  !### A part of xmat_column*pf does not involve MPI communication ###
  !###################################################################
!!$  do ispin=1,nspin
!!$     do jspin=1,nspin
!!$        do idim=1,ndim
!!$           if(abs(GAMMA10d(idim,ispin,jspin)).GT.1d-2)then
!!$              
!!$              do isite=1,nsite_local
!!$!$omp parallel
!!$!$omp do             
!!$                 do imat=1,nmat_block
!!$                    do jmat=1,nmat_block
!!$                       do kmat=1,nmat_block
!!$                          pf2(imat,jmat,ispin,isite)=&
!!$                               &pf2(imat,jmat,ispin,isite)-&
!!$                               &dcmplx(lattice_spacing)*&
!!$                               &dconjg(GAMMA10d(idim,jspin,ispin))&
!!$                               &*xmat_column((iblock-1)&
!!$                               &*nmat_block+imat,kmat,idim,isite)&
!!$                               &*pf1(kmat,jmat,jspin,isite)
!!$                       end do
!!$                    end do
!!$                 end do
!!$!$omp end do
!!$!$omp end parallel
!!$              end do
!!$           end if
!!$        end do
!!$     end do
!!$  end do
   include 'multiplication_dagger_direct_1.f'
  !############################################
  !### MPI-communication for xmat_column*pf ###
  !############################################
  kblock_send=iblock
  kblock_rcv=iblock
  do ishift=1,nblock-1
     kblock_send=kblock_send+1
     kblock_rcv=kblock_rcv-1
     if(kblock_send.EQ.nblock+1)then
        kblock_send=1
     end if
     if(kblock_rcv.EQ.0)then
        kblock_rcv=nblock
     end if  
     send_rank=(isublat-1)*nblock*nblock+(kblock_send-1)*nblock+jblock-1
     receive_rank=(isublat-1)*nblock*nblock+(kblock_rcv-1)*nblock+jblock-1     
     mat_send=(0d0,0d0)
!!$     do ispin=1,nspin
!!$        do jspin=1,nspin
!!$           do idim=1,ndim
!!$              if(abs(GAMMA10d(idim,ispin,jspin)).GT.1d-2)then
!!$                 do isite=1,nsite_local
!!$!$omp parallel
!!$!$omp do                 
!!$                    do imat=1,nmat_block
!!$                       do jmat=1,nmat_block
!!$                          do kmat=1,nmat_block
!!$                             mat_send(imat,jmat,ispin,isite)=&
!!$                                  &mat_send(imat,jmat,ispin,isite)-&
!!$                                  &dcmplx(lattice_spacing)*&
!!$                                  &dconjg(GAMMA10d(idim,jspin,ispin))&
!!$                                  &*xmat_column((kblock_send-1)&
!!$                                  &*nmat_block+imat,kmat,idim,isite)&
!!$                                  &*pf1(kmat,jmat,jspin,isite)        
!!$                          end do
!!$                       end do
!!$                    end do
!!$!$omp end do
!!$!$omp end parallel
!!$                 end do
!!$              end if
!!$           end do
!!$        end do
!!$     end do
  include 'multiplication_dagger_direct_1.f'
     tag=1
     call MPI_Isend(mat_send(1,1,1,1),&
          &nmat_block*nmat_block*nsite_local*nspin,&
          &MPI_DOUBLE_COMPLEX,&
          &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
     call MPI_Recv(mat_rcv(1,1,1,1),&
          &nmat_block*nmat_block*nsite_local*nspin,&
          &MPI_DOUBLE_COMPLEX,&
          &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
     call MPI_Wait(ireq,status,ierr)
     do ispin=1,nspin
        do isite=1,nsite_local
!$omp parallel
!$omp do
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 pf2(imat,jmat,ispin,isite)=&
                      &pf2(imat,jmat,ispin,isite)&
                      &+mat_rcv(imat,jmat,ispin,isite)
              end do
           end do
!$omp end do
!$omp end parallel
        end do
     end do
  end do
  !################################################################
  !### A part of pf*xmat_row does not involve MPI communication ###
  !################################################################
!!$  do ispin=1,nspin
!!$     do jspin=1,nspin
!!$        do idim=1,ndim
!!$           if(abs(GAMMA10d(idim,ispin,jspin)).GT.1d-2)then
!!$              do isite=1,nsite_local
!!$!$omp parallel 
!!$!$omp do             
!!$                 do imat=1,nmat_block
!!$                    do jmat=1,nmat_block
!!$                       do kmat=1,nmat_block
!!$                          pf2(imat,jmat,ispin,isite)=&
!!$                               &pf2(imat,jmat,ispin,isite)+&
!!$                                  &dcmplx(lattice_spacing)*&
!!$                                  &dconjg(GAMMA10d(idim,jspin,ispin))&
!!$                                  &*pf1(imat,kmat,jspin,isite)&
!!$                                  &*xmat_row(kmat,(jblock-1)&
!!$                                  &*nmat_block+jmat,idim,isite)
!!$                       end do
!!$                    end do
!!$                 end do
!!$!$omp end do
!!$!$omp end parallel
!!$              end do
!!$           end if
!!$        end do
!!$     end do
!!$  end do
  include 'multiplication_dagger_direct_1.f'
  !#########################################
  !### MPI-communication for pf*xmat_row ###
  !#########################################
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
     mat_send=(0d0,0d0)
!!$     do ispin=1,nspin
!!$        do jspin=1,nspin
!!$           do idim=1,ndim
!!$              if(abs(GAMMA10d(idim,ispin,jspin)).GT.1d-2)then
!!$                 do isite=1,nsite_local
!!$!$omp parallel
!!$!$omp do                 
!!$                    do imat=1,nmat_block
!!$                       do jmat=1,nmat_block
!!$                          do kmat=1,nmat_block
!!$                             mat_send(imat,jmat,ispin,isite)=&
!!$                                  &mat_send(imat,jmat,ispin,isite)+&
!!$                                  &dcmplx(lattice_spacing)*&
!!$                                  &dconjg(GAMMA10d(idim,jspin,ispin))&
!!$                                  &*pf1(imat,kmat,jspin,isite)&
!!$                                  &*xmat_row(kmat,(kblock_send-1)&
!!$                                  *nmat_block+jmat,idim,isite)
!!$                          end do
!!$                       end do
!!$                    end do
!!$!$omp end do
!!$!$omp end parallel
!!$                 end do
!!$              end if
!!$           end do
!!$        end do
!!$     end do
       include 'multiplication_dagger_direct_1.f'
     call MPI_Isend(mat_send(1,1,1,1),&
          &nmat_block*nmat_block*nsite_local*nspin,&
          &MPI_DOUBLE_COMPLEX,&
          &send_rank,tag,MPI_COMM_WORLD,ireq,ierr)
     call MPI_Recv(mat_rcv(1,1,1,1),&
          &nmat_block*nmat_block*nsite_local*nspin,&
          &MPI_DOUBLE_COMPLEX,&
          &receive_rank,tag,MPI_COMM_WORLD,status,ierr)
     call MPI_Wait(ireq,status,ierr)
     
     do ispin=1,nspin
        do isite=1,nsite_local
!$omp parallel
!$omp do
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 pf2(imat,jmat,ispin,isite)=&
                      &pf2(imat,jmat,ispin,isite)&
                      &+mat_rcv(imat,jmat,ispin,isite)
              end do
           end do
!$omp end do
!$omp end parallel
        end do
     end do
  end do
  !*************************************************************
  !*************************************************************
  !*** Plane wave deformation; no need for MPI communication ***
  !*************************************************************
  !*************************************************************
  if(nbmn.EQ.1)then
     Gam12=(0d0,0d0)
     Gam123=(0d0,0d0)
     do jspin=1,nspin
        do ispin=1,nspin
           do kspin=1,nspin
              Gam12(ispin,jspin)=Gam12(ispin,jspin)&
                   +Gamma10d(1,ispin,kspin)*Gamma10d(2,kspin,jspin)
           end do
        end do
     end do
     do jspin=1,nspin
        do ispin=1,nspin
           do kspin=1,nspin
              Gam123(ispin,jspin)=Gam123(ispin,jspin)&
                   +Gam12(ispin,kspin)*Gamma10d(3,kspin,jspin)
           end do
        end do
     end do
     Gam123=Gam123*dcmplx(flux)*(0d0,-0.75d0)*dcmplx(lattice_spacing)
     do isite=1,nsite_local
        do jspin=1,nspin
           do ispin=1,nspin
!$omp parallel
!$omp do
              do imat=1,nmat_block
                 do jmat=1,nmat_block
                    pf2(imat,jmat,ispin,isite)=&
                         pf2(imat,jmat,ispin,isite)&
                         +dconjg(Gam123(jspin,ispin))*pf1(imat,jmat,jspin,isite)
                 end do
              end do
!$omp end do
!$omp end parallel
           end do
        end do
     end do
  end if
  
  return

END SUBROUTINE Multiply_Dirac_dagger
