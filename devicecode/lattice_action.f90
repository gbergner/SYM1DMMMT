!*******************************************************
!*** Calculate the action at each MPI process. *********
!*** Here, only kinetic+potential+FP is calculated. ****
!*******************************************************
module lattice_action
contains
subroutine Calc_action_device(temperature,xmat,alpha,action,phase)

  use compiletimeconstants
  implicit none


  double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double precision, intent(in) :: temperature
  double complex, intent(in) :: phase(1:nmat,1:nmat,1:2)
  double precision, intent(in) :: alpha(1:nmat)
  !$acc declare present(xmat,temperature,phase,alpha)
  double precision, intent(out) :: action
  double precision kinetic,potential,gauge_fixing
  double precision lattice_spacing
  double complex uxumx,com
  integer isite
  integer idim,jdim
  integer imat,jmat,kmat
  double complex ei,ej
  
  lattice_spacing=1d0/temperature/dble(nsite)
  !**********************
  !*** potential term ***
  !**********************
  potential=0d0

  !$acc kernels
  do idim=1,ndim-1
   do jdim=idim+1,ndim
      do isite=1,nsite
           do imat=1,nmat
              do jmat=1,nmat
                 com=(0d0,0d0)
                 do kmat=1,nmat
                    com=com&
                         &+xmat(imat,kmat,idim,isite)*xmat(kmat,jmat,jdim,isite)&
                         &-xmat(imat,kmat,jdim,isite)*xmat(kmat,jmat,idim,isite)
                 end do
                 potential=potential&
                            &+dble(com*dconjg(com))
              end do
           end do
        end do
     end do
  end do
  !$acc end kernels

  potential=potential*0.5d0*dble(nmat)*lattice_spacing
  !********************
  !*** kinetic term ***
  !********************
  kinetic=0d0
  !********************
  !*** Naive Action ***
  !********************
  if(nimprove.EQ.0)then
  !$acc kernels
     do isite=1,nsite
        do idim=1,ndim
           !u(t)*x(t+a)*u^dagger(t) - x(t)
           do imat=1,nmat
              do jmat=1,nmat
                 !exp(i*alpha_i)
                 !exp(-i*alpha_j)
                 uxumx=&
                      &phase(imat,jmat,1)*xmat(imat,jmat,idim,isite+1)&
                      &-xmat(imat,jmat,idim,isite)
                 kinetic=kinetic+dble(uxumx*dconjg(uxumx))
                 !+dble(uxumx(imat,jmat)*uxumx(jmat,imat))
              end do
           end do
        end do
     end do
     !$acc end kernels
     !***********************
     !*** Improved Action ***
     !***********************
  else if(nimprove.EQ.1)then
    !$acc kernels
    do isite=1,nsite
        do idim=1,ndim
           !-0.5*u^2*x(t+a)*(u^dagger)^2 - 2*u*x(t+a)*u^dagger - 1.5*x(t)
           do imat=1,nmat
              do jmat=1,nmat
                 
                 !exp(i*alpha_i)
                 !exp(-i*alpha_j)
                 uxumx=&
                      &-(0.5d0,0d0)*phase(imat,jmat,2)*xmat(imat,jmat,idim,isite+2)&
                      &+(2d0,0d0)*phase(imat,jmat,1)*xmat(imat,jmat,idim,isite+1)&
                      &-(1.5d0,0d0)*xmat(imat,jmat,idim,isite)
                 kinetic=kinetic+dble(uxumx*dconjg(uxumx))
              end do
           end do
        end do
     end do
    !$acc end kernels

  end if
  
  kinetic=kinetic*0.5d0*dble(nmat)/lattice_spacing
  !*************************
  !*** gauge-fixing term ***
  !*************************
  gauge_fixing=0d0
  if(ngauge.EQ.1)then
    !$acc kernels
     do imat=1,nmat-1
        do jmat=imat+1,nmat
           gauge_fixing=gauge_fixing&
                -2d0*dlog(dabs(dsin(0.5d0*(alpha(imat)-alpha(jmat)))))
        end do
     end do
     !$acc end kernels
  end if
  action=kinetic+potential+gauge_fixing
  
  return

END subroutine Calc_action_device

! impose the constraint max(alpha_i-alpha_j) < 2*pi.
! info=0 -> OK, info=1 -> constraint is violated.
SUBROUTINE check_alpha_constraint_device(alpha,info)
  use compiletimeconstants
  implicit none
  !***** input *****
  double precision alpha(1:nmat)
  !$acc declare present(alpha)
  !***** output *****
  integer info
  !******************
  integer imat
  double precision mmax,mmin,pi

  !$acc kernels
  mmax=maxval(alpha)
  mmin=minval(alpha)
  !$acc end kernels
  pi=2d0*dasin(1d0)
  if(mmax-mmin.LT.2d0*pi)then
     info=0
  else
     info=1
  end if

  return

END SUBROUTINE check_alpha_constraint_device

end module lattice_action
