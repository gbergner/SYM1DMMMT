!*******************************************************
!*** Calculate the action at each MPI process. *********
!*** Here, only kinetic+potential+FP is calculated. ****
!*******************************************************
subroutine Calc_action(temperature,xmat,alpha,action,ngauge)
  
  implicit none

  include '../staticparameters.f90'
  !***** input *****
  integer ngauge
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double precision temperature
  double precision alpha(1:nmat)
  !***** output *****
  double precision action
  !******************
  double precision kinetic,potential,gauge_fixing
  double precision lattice_spacing
  double complex commutator(1:nmat,1:nmat),uxumx(1:nmat,1:nmat)
  integer isite
  integer idim,jdim
  integer imat,jmat,kmat
  double complex ei,ej
  
  lattice_spacing=1d0/temperature/dble(nsite)
  !**********************
  !*** potential term ***
  !**********************
  potential=0d0
  do isite=1,nsite
     do idim=1,ndim-1
        do jdim=idim+1,ndim
           commutator=(0d0,0d0)
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    commutator(imat,jmat)=commutator(imat,jmat)&
                         &+xmat(imat,kmat,idim,isite)*xmat(kmat,jmat,jdim,isite)&
                         &-xmat(imat,kmat,jdim,isite)*xmat(kmat,jmat,idim,isite)
                 end do
              end do
           end do
           do imat=1,nmat
              do jmat=1,nmat
                 potential=potential&
                      &+dble(commutator(imat,jmat)*dconjg(commutator(imat,jmat)))
              end do
           end do
           
        end do
     end do
  end do
  potential=potential*0.5d0*dble(nmat)*lattice_spacing
  !********************
  !*** kinetic term ***
  !********************
  kinetic=0d0
  !********************
  !*** Naive Action ***
  !********************
  if(nimprove.EQ.0)then
     do isite=1,nsite
        do idim=1,ndim
           !u(t)*x(t+a)*u^dagger(t) - x(t)
           do imat=1,nmat
              do jmat=1,nmat
                 
                 !exp(i*alpha_i)
                 ei=dcmplx(dcos(alpha(imat)/dble(nsite)))&
                      &+(0d0,1d0)*dcmplx(dsin(alpha(imat)&
                      &/dble(nsite)))
                 !exp(-i*alpha_j)
                 ej=dcmplx(dcos(alpha(jmat)/dble(nsite)))&
                      &-(0d0,1d0)*dcmplx(dsin(alpha(jmat)&
                      &/dble(nsite)))
                 uxumx(imat,jmat)=&
                      &ei*xmat(imat,jmat,idim,isite+1)*ej&
                      &-xmat(imat,jmat,idim,isite)
              end do
           end do
           do imat=1,nmat
              do jmat=1,nmat
                 kinetic=kinetic+dble(uxumx(imat,jmat)*uxumx(jmat,imat))
              end do
           end do
           
        end do
     end do
     !***********************
     !*** Improved Action ***
     !***********************
  else if(nimprove.EQ.1)then

    do isite=1,nsite
        do idim=1,ndim
           !-0.5*u^2*x(t+a)*(u^dagger)^2 - 2*u*x(t+a)*u^dagger - 1.5*x(t)
           do imat=1,nmat
              do jmat=1,nmat
                 
                 !exp(i*alpha_i)
                 ei=dcmplx(dcos(alpha(imat)/dble(nsite)))&
                      &+(0d0,1d0)*dcmplx(dsin(alpha(imat)&
                      &/dble(nsite)))
                 !exp(-i*alpha_j)
                 ej=dcmplx(dcos(alpha(jmat)/dble(nsite)))&
                      &-(0d0,1d0)*dcmplx(dsin(alpha(jmat)&
                      &/dble(nsite)))
                 uxumx(imat,jmat)=&
                      &-(0.5d0,0d0)*ei*ei*xmat(imat,jmat,idim,isite+2)*ej*ej&
                      &+(2d0,0d0)*ei*xmat(imat,jmat,idim,isite+1)*ej&
                      &-(1.5d0,0d0)*xmat(imat,jmat,idim,isite)
              end do
           end do
           
           do imat=1,nmat
              do jmat=1,nmat
                 kinetic=kinetic+dble(uxumx(imat,jmat)*dconjg(uxumx(imat,jmat)))
              end do
           end do

           
        end do
     end do


  end if
  
  kinetic=kinetic*0.5d0*dble(nmat)/lattice_spacing
  !*************************
  !*** gauge-fixing term ***
  !*************************
  gauge_fixing=0d0
  if(ngauge.EQ.0)then
     do imat=1,nmat-1
        do jmat=imat+1,nmat
           gauge_fixing=gauge_fixing&
                -2d0*dlog(dabs(dsin(0.5d0*(alpha(imat)-alpha(jmat)))))
        end do
     end do
  end if
  action=kinetic+potential+gauge_fixing
  
  return

END subroutine Calc_action
