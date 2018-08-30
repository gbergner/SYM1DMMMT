!**********************************************************
!*** The internal energy is calculated, *******************
!*** by subtracting a constant from the bosonic action. ***
!*** note that "energy" is actualy E/N^2. *****************
!**********************************************************
subroutine Calc_energy(temperature,xmat,alpha,energy,nbmn,flux)

  implicit none

  include '../staticparameters.f90'
  !***** input *****
  integer nbmn
  double precision temperature,flux
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double precision alpha(1:nmat)
  !***** output *****
  double precision energy
  !******************
  double precision action,kinetic,potential,potential_BMN
  double precision lattice_spacing
  double complex commutator(1:nmat,1:nmat),uxumx(1:nmat,1:nmat)
  integer isite
  integer idim,jdim
  integer imat,jmat,kmat
  double complex ei,ej
  !**** For BMN deformation ****
  double complex x23(1:nmat,1:nmat),x32(1:nmat,1:nmat),&
       trx123,trx132
  double precision trx2_123,trx2_456789

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
                         +xmat(imat,kmat,idim,isite)*xmat(kmat,jmat,jdim,isite)&
                         -xmat(imat,kmat,jdim,isite)*xmat(kmat,jmat,idim,isite)
                 end do
              end do
           end do
           do imat=1,nmat
              do jmat=1,nmat
                 potential=potential&
                      +dble(commutator(imat,jmat)*dconjg(commutator(imat,jmat)))
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
                      +(0d0,1d0)*dcmplx(dsin(alpha(imat)/dble(nsite)))
                 !exp(-i*alpha_j)
                 ej=dcmplx(dcos(alpha(jmat)/dble(nsite)))&
                      -(0d0,1d0)*dcmplx(dsin(alpha(jmat)/dble(nsite)))
                 uxumx(imat,jmat)=&
                      ei*xmat(imat,jmat,idim,isite+1)*ej&
                      -xmat(imat,jmat,idim,isite)
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
                      -(0d0,1d0)*dcmplx(dsin(alpha(jmat)&
                      &/dble(nsite)))
                 uxumx(imat,jmat)=&
                      -(0.5d0,0d0)*ei*ei*xmat(imat,jmat,idim,isite+2)*ej*ej&
                      +(2d0,0d0)*ei*xmat(imat,jmat,idim,isite+1)*ej&
                      -(1.5d0,0d0)*xmat(imat,jmat,idim,isite)
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
  !******************************
  !*** plane wave deformation ***
  !******************************
  potential_BMN=0d0
  if(nbmn.EQ.1)then
     !*****************
     !*** mass term ***
     !*****************
     trx2_123=0d0
     trx2_456789=0d0
     do isite=1,nsite
        do imat=1,nmat
           do jmat=1,nmat
              do idim=1,3
                 trx2_123=trx2_123&
                      +dble(xmat(imat,jmat,idim,isite)&
                      *dconjg(xmat(imat,jmat,idim,isite)))
              end do
              do idim=4,9
                 trx2_456789=trx2_456789&
                      +dble(xmat(imat,jmat,idim,isite)&
                      *dconjg(xmat(imat,jmat,idim,isite)))
              end do
           end do
        end do
     end do
     potential_BMN=potential_BMN+flux*flux*(0.5d0*trx2_123+0.125d0*trx2_456789)&
          *lattice_spacing*dble(nmat)
     !******************
     !*** cubic term ***
     !******************
     trx123=(0d0,0d0)
     trx132=(0d0,0d0)
     do isite=1,nsite
        x23=(0d0,0d0)
        x32=(0d0,0d0)
        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 x23(imat,jmat)=x23(imat,jmat)&
                      +xmat(imat,kmat,2,isite)*xmat(kmat,jmat,3,isite)
                 x32(imat,jmat)=x32(imat,jmat)&
                      +xmat(imat,kmat,3,isite)*xmat(kmat,jmat,2,isite)
              end do
           end do
        end do
        do imat=1,nmat
           do jmat=1,nmat
              trx123=trx123+xmat(imat,jmat,1,isite)*x23(jmat,imat)
              trx132=trx132+xmat(imat,jmat,1,isite)*x32(jmat,imat)
           end do
        end do
     end do
     potential_BMN=potential_BMN+&
          dble((0d0,3d0)*(trx123-trx132))*flux*lattice_spacing*dble(nmat)
  end if

  !The gauge fixing term is not needed when we calculate the energy.
  action=kinetic+potential+potential_BMN
  
  energy=-3d0*temperature*action/dble(nmat*nmat)+1.5d0*temperature*dble(ndim*nsite)


  !after summing up contributions from all ranks, 
  !we must subtract the zero-mode of the U(1) part, 
  !1.5d0*temperature*dble(ndim)/dble(nmat*nmat). 
  energy=energy-1.5d0*temperature*dble(ndim)/dble(nmat*nmat)
  
  return

END subroutine Calc_energy
