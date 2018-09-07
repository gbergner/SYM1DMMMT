!******************************************
!******************************************
!*** Calculate the action at each node. ***
!******************************************
!******************************************
subroutine Calc_action(temperature,xmat,xmat_row,xmat_column,alpha,action,myrank,ngauge)

  implicit none
  include 'size_parallel.h'
  integer myrank
  integer ngauge
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,-(nmargin-1):nsite_local+nmargin)
  double precision alpha(1:nmat_block*nblock)

  double precision action,kinetic,potential,gauge_fixing
  double precision temperature,lattice_spacing

  double complex commutator(1:nmat_block,1:nmat_block)
  double complex uxumx(1:nmat_block,1:nmat_block)
  integer isite,isite_p1
  integer idim,jdim
  integer imat,jmat,kmat
  double complex ei,ej

  double complex xmat_row(1:nmat_block,1:nmat_block*nblock,1:ndim,1:nsite_local)
  double complex xmat_column(1:nmat_block*nblock,1:nmat_block,1:ndim,1:nsite_local)

  integer iblock,jblock,isublat
  call who_am_i(myrank,isublat,iblock,jblock)

  !nprocs=nsublat*nmat_block*nmat_block
  lattice_spacing=1d0/temperature/dble(nsite_local*nsublat)
  !**********************
  !*** potential term ***
  !**********************
  potential=0d0
  do isite=1,nsite_local
     !isite=0, nsite_local+1 are considered at neighboring nodes.
     do idim=1,ndim-1
        do jdim=idim+1,ndim
           commutator=(0d0,0d0)
!$omp parallel 
!$omp do
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 do kmat=1,nmat_block*nblock
                    commutator(imat,jmat)=commutator(imat,jmat)&
                         +xmat_row(imat,kmat,idim,isite)*xmat_column(kmat,jmat,jdim,isite)&
                         -xmat_row(imat,kmat,jdim,isite)*xmat_column(kmat,jmat,idim,isite)
                 end do
              end do
           end do
!$omp end do
!$omp end parallel
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 potential=potential&
                      +dble(commutator(imat,jmat)*dconjg(commutator(imat,jmat)))
              end do
           end do           
        end do
     end do
  end do
  potential=potential*0.5d0*dble(nmat_block*nblock)*lattice_spacing
  !********************
  !*** kinetic term ***
  !********************
  kinetic=0d0
  if(nimprove.EQ.0)then
     do isite=1,nsite_local
        !isite=0, nsite_local+1 are considered at neighboring nodes.
        isite_p1=isite+1
        
        do idim=1,ndim
           !u(t)*x(t+a)*u^dagger(t) - x(t)
!$omp parallel
!$omp do           
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 
                 !exp(i*alpha_i)
                 ei=dcmplx(dcos(alpha(imat+(iblock-1)*nmat_block)&
                      &/dble(nsite_local*nsublat)))&
                      &+(0d0,1d0)*dcmplx(dsin(alpha(imat+(iblock-1)*nmat_block)&
                      &/dble(nsite_local*nsublat)))
                 !exp(-i*alpha_j)
                 ej=dcmplx(dcos(alpha(jmat+(jblock-1)*nmat_block)&
                      &/dble(nsite_local*nsublat)))&
                      -(0d0,1d0)*dcmplx(dsin(alpha(jmat+(jblock-1)*nmat_block)&
                      &/dble(nsite_local*nsublat)))
                 uxumx(imat,jmat)=&
                      ei*xmat(imat,jmat,idim,isite_p1)*ej&
                      -xmat(imat,jmat,idim,isite)
              end do
           end do
!$omp end do
!$omp end parallel        
           do imat=1,nmat_block
              do jmat=1,nmat_block 
                 !kinetic=kinetic+dble(uxumx(imat,jmat)*uxumx(jmat,imat))
                 kinetic=kinetic+dble(uxumx(imat,jmat)*dconjg(uxumx(imat,jmat)))
              end do
           end do
           
        end do
     end do

  else if(nimprove.EQ.1)then  

    do isite=1,nsite_local

        do idim=1,ndim
           !-0.5*u^2*x(t+2a)*(u^dagger)^2 + 2*u*x(t+a)*u^dagger - 1.5*x(t)
!$omp parallel
!$omp do            
           do imat=1,nmat_block
              do jmat=1,nmat_block
                 
                 !exp(i*alpha_i)
                 ei=dcmplx(dcos(alpha(imat+(iblock-1)*nmat_block)&
                      &/dble(nsite_local*nsublat)))&
                      &+(0d0,1d0)*dcmplx(dsin(alpha(imat+(iblock-1)*nmat_block)&
                      &/dble(nsite_local*nsublat)))
                 !exp(-i*alpha_j)
                 ej=dcmplx(dcos(alpha(jmat+(jblock-1)*nmat_block)&
                      &/dble(nsite_local*nsublat)))&
                      -(0d0,1d0)*dcmplx(dsin(alpha(jmat+(jblock-1)*nmat_block)&
                      &/dble(nsite_local*nsublat)))
                 uxumx(imat,jmat)=&
                      -(0.5d0,0d0)*ei*ei*xmat(imat,jmat,idim,isite+2)*ej*ej&
                      +(2d0,0d0)*ei*xmat(imat,jmat,idim,isite+1)*ej&
                      -(1.5d0,0d0)*xmat(imat,jmat,idim,isite)
              end do
           end do
!$omp end do
!$omp end parallel         
           do imat=1,nmat_block
              do jmat=1,nmat_block 
                 !kinetic=kinetic+dble(uxumx(imat,jmat)*uxumx(jmat,imat))
                 kinetic=kinetic+dble(uxumx(imat,jmat)*dconjg(uxumx(imat,jmat)))
              end do
           end do
   
     
        end do
     end do

  end if

  kinetic=kinetic*0.5d0*dble(nmat_block*nblock)/lattice_spacing
  gauge_fixing=0d0
  if(ngauge.eq.0)then
     if(myrank.EQ.0)then
        !*************************
        !*** gauge-fixing term ***
        !************************* 
        do imat=1,nmat_block*nblock-1
           do jmat=imat+1,nmat_block*nblock
              gauge_fixing=gauge_fixing&
                   -2d0*dlog(dabs(dsin(0.5d0*(alpha(imat)-alpha(jmat)))))
              !  write(*,*)gauge_fixing
           end do
        end do
     end if
  end if
  action=kinetic+potential+gauge_fixing

  return

END subroutine Calc_action
