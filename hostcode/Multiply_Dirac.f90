!pf2=(-i)*M*pf1, M: Dirac op
!Please be careful about the overall factor (-i). 
SUBROUTINE Multiply_Dirac(temperature,xmat,alpha,&
     pf1,pf2,GAMMA10d,nbmn,flux)
  
  implicit none
  include '../staticparameters.f90'
  !***** input *****
  integer nprocs,nbmn
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double precision alpha(1:nmat)
  double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
  double complex gam
  double complex pf1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
  !***** output *****
  double complex pf2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
  !******************
  double complex Gam123(1:nspin,1:nspin),Gam12(1:nspin,1:nspin)
  double complex phase,phase2
  double precision aij
  
  double precision temperature,lattice_spacing,flux
  integer imat,jmat,kmat
  integer idim
  integer ispin,jspin,kspin
  integer isite
  
  lattice_spacing=1d0/temperature/dble(nsite)
  !**********************
  !**********************
  !***  kinetic part  ***
  !**********************
  !**********************
  pf2=(0d0,0d0)

  !********************
  !*** Naive Action ***
  !********************
  if(nimprove.EQ.0)then
     do isite=1,nsite
        do jmat=1,nmat
           do imat=1,nmat
              aij=alpha(imat)-alpha(jmat)
              aij=aij/dble(nsite)
              phase=dcmplx(dcos(aij))+(0d0,1d0)*dcmplx(dsin(aij))
              do ispin=1,8
                 pf2(imat,jmat,ispin+8,isite)=&
                      pf2(imat,jmat,ispin+8,isite)&
                      +(1d0,0d0)*phase*pf1(imat,jmat,ispin,isite+1)&
                      -(1d0,0d0)*pf1(imat,jmat,ispin,isite)
                 
                 pf2(imat,jmat,ispin,isite)=&
                      pf2(imat,jmat,ispin,isite)&
                      -(1d0,0d0)*dconjg(phase)*pf1(imat,jmat,ispin+8,isite-1)&
                      +(1d0,0d0)*pf1(imat,jmat,ispin+8,isite)
                 
              end do
              
              
           end do
        end do
     end do
     !***********************
     !*** Improved Action ***
     !***********************
  else if(nimprove.EQ.1)then
     do isite=1,nsite
        do jmat=1,nmat
           do imat=1,nmat
              aij=alpha(imat)-alpha(jmat)
              aij=aij/dble(nsite)
              phase=dcmplx(dcos(aij))+(0d0,1d0)*dcmplx(dsin(aij))
              phase2=phase*phase
              do ispin=1,8
                 pf2(imat,jmat,ispin+8,isite)=&
                      pf2(imat,jmat,ispin+8,isite)&
                      -(0.5d0,0d0)*phase2*pf1(imat,jmat,ispin,isite+2)&
                      +(2d0,0d0)*phase*pf1(imat,jmat,ispin,isite+1)&
                      -(1.5d0,0d0)*pf1(imat,jmat,ispin,isite)
               
                 pf2(imat,jmat,ispin,isite)=&
                      pf2(imat,jmat,ispin,isite)&
                      +(0.5d0,0d0)*dconjg(phase2)*pf1(imat,jmat,ispin+8,isite-2)&
                      -(2d0,0d0)*dconjg(phase)*pf1(imat,jmat,ispin+8,isite-1)&
                      +(1.5d0,0d0)*pf1(imat,jmat,ispin+8,isite)
              end do
              
              
           end do
        end do
     end do
  end if
  !**************************
  !**************************
  !***  interaction part  ***
  !**************************
  !**************************
  include 'multiply_direct.f'
  !******************************
  !******************************
  !*** Plane wave deformation ***
  !******************************
  !******************************
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
     Gam123=lattice_spacing*Gam123*dcmplx(flux)*(0d0,-0.75d0)
     do isite=1,nsite
        do ispin=1,nspin
           do jspin=1,nspin
              do imat=1,nmat
                 do jmat=1,nmat
                    pf2(imat,jmat,ispin,isite)=&
                         pf2(imat,jmat,ispin,isite)&
                         +Gam123(ispin,jspin)*pf1(imat,jmat,jspin,isite)
                 end do
              end do
           end do
        end do
     end do
  end if
  
  return
  
END SUBROUTINE Multiply_Dirac
