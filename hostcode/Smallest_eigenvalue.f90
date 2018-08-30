! Calculate the smallest eigenvalue of (D^dag*D),  
! by multiplying (D^dag*D)^{-1} many times to a random vector. 
SUBROUTINE Smallest_eigenvalue(temperature,&
     &xmat,alpha,GAMMA10d,neig,smallest_eig,nbc,max_err,max_iteration,nbmn,&
     &flux,nsite,nmat,nremez_md,nremez_pf,nimprove,ngauge,nmargin)

  implicit none

  !include '../staticparameters.f90'

  integer nsite
  integer nmat
  integer nremez_md,nremez_pf
  integer nimprove
  integer ngauge
  integer npf
  parameter(npf=1)!only inside this routine
  integer ndim ! fix to 9
  parameter(ndim=9)
  integer nspin ! fix to 16
  parameter(nspin=16)
  integer nmargin!size of the margin
  integer neig,ieig,nbc,nremez,info,max_iteration,iteration,nbmn
  double precision max_err

  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double precision alpha(1:nmat)
  double precision bcoeff(1:1)!use multimass CG solver, with nremez=1
  double complex phi1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
  double complex phi2(1:1,1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)!use multimass CG solver, with nremez=1 and bcoeff=0
  double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
  double complex trace
  double precision smallest_eig,norm
  !pf2=M*pf1, M: Dirac op

  double precision temperature,flux
  integer imat,jmat
  integer ispin
  integer isite
  integer i

  integer nprocs,myrank,IERR

  double precision r1,r2

  !***********************************
  !**** generate a random vector. ****
  !***********************************
  norm=0d0
  do imat=1,nmat     
     do jmat=1,nmat
        do ispin=1,nspin
           do isite=1,nsite
              call BoxMuller(r1,r2)
              phi1(imat,jmat,ispin,isite,1)=&
                   &(dcmplx(r1)+dcmplx(r2)*(0D0,1D0))/dcmplx(dsqrt(2d0))
              norm=norm+(r1*r1+r2*r2)*0.5d0
           end do
        end do
     end do
  end do
  !****************************
  !*** traceless projection ***
  !****************************
  do ispin=1,nspin
     do isite=1,nsite
        trace=(0d0,0d0)
        do imat=1,nmat
           trace=trace+phi1(imat,imat,ispin,isite,1)
        end do
        trace=trace/dcmplx(nmat)
        do imat=1,nmat
           norm=norm-dble(phi1(imat,imat,ispin,isite,1)*dconjg(phi1(imat,imat,ispin,isite,1)))
           phi1(imat,imat,ispin,isite,1)=phi1(imat,imat,ispin,isite,1)-trace
           norm=norm+dble(phi1(imat,imat,ispin,isite,1)*dconjg(phi1(imat,imat,ispin,isite,1)))
        end do
        
     end do
  end do

  norm=dsqrt(norm)
  norm=1d0/norm
  phi1=phi1*norm
  !*********************************
  !*** adjust the margin and b.c.***
  !*********************************
  call Adjust_margin_and_bc_pf(phi1,nbc)
  !******************************************
  !*** random vector has been generated.  ***
  !******************************************
  nremez=1
  bcoeff(1)=0d0
  do ieig=1,neig
     !*****************************************************
     !*** phi1 -> phi2=(D^dag*D)^{-1}*phi1 -> phi1=phi2 ***
     !*****************************************************
     call solver_biCGm(nbc,nbmn,nremez,&
     xmat,alpha,phi1,phi2,GAMMA10d,&
     bcoeff,max_err,max_iteration,iteration,&
     temperature,flux,info)

!call solver_biCGm(nbc,nprocs,nremez,bcoeff,temperature,&
!     xmat,alpha,phi1,phi2,GAMMA10d,max_err,max_iteration,iteration,myrank,&
!     nbmn,flux,info)
     do imat=1,nmat     
        do jmat=1,nmat
           do ispin=1,nspin
              do isite=-(nmargin-1),nsite+nmargin
                 phi1(imat,jmat,ispin,isite,1)=phi2(1,imat,jmat,ispin,isite,1)
              end do
           end do
        end do
     end do
     !****************************************
     !*** traceless projection;            ***
     !*** to avoid the zero mode to appear ***
     !*** as numerical artifact            ***
     !****************************************
     ! (Practically, it does not seem to be necessary.) 
     !do ispin=1,nspin
     !   do isite=-nimprove,nsite_local+1+nimprove
     !      trace=(0d0,0d0)
     !      do imat=1,nmat
     !         trace=trace+phi1(imat,imat,ispin,isite)
     !      end do
     !      trace=trace/dcmplx(nmat)
     !      do imat=1,nmat
     !         phi1(imat,imat,ispin,isite)=phi1(imat,imat,ispin,isite)-trace
     !      end do
     !      
     !   end do
     !end do

     !************************************
     !*** calculate the norm of phi1.  ***
     !************************************
     norm=0d0
     do imat=1,nmat     
        do jmat=1,nmat
           do ispin=1,nspin
              do isite=1,nsite
                 norm=norm&
                      &+dble(dconjg(phi1(imat,jmat,ispin,isite,1)&
                      &*dconjg(phi1(imat,jmat,ispin,isite,1))))
              end do
           end do
        end do
     end do
     norm=dsqrt(norm)
     norm=1d0/norm
     phi1=phi1*norm
     norm=1d0/norm
    
     !write(*,*)1d0/norm
  end do
  
  smallest_eig=1d0/norm

  return
  
END SUBROUTINE Smallest_eigenvalue
