! Calculate the largest eigenvalue of (D^dag*D),  
! by multiplying (D^dag*D) many times to a random vector. 
SUBROUTINE Largest_eigenvalue(temperature,xmat,alpha,GAMMA10d,neig,largest_eig,nbc,nbmn,flux,&
     &nsite,nmat,nremez_md,nremez_pf,nimprove,ngauge,nmargin)

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



  !***** input *****
  integer neig,nbc,nbmn
  double precision temperature,flux
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double precision alpha(1:nmat)
  double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
  !***** output *****
  double precision largest_eig
  !******************
  double complex phi1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
  double complex phi2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
  double complex trace
  double precision norm
  integer imat,jmat
  integer ispin
  integer isite
  integer ieig,i
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
              phi1(imat,jmat,ispin,isite)=&
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
           trace=trace+phi1(imat,imat,ispin,isite)
        end do
        trace=trace/dcmplx(nmat)
        do imat=1,nmat
           norm=norm-dble(phi1(imat,imat,ispin,isite)*dconjg(phi1(imat,imat,ispin,isite)))
           phi1(imat,imat,ispin,isite)=phi1(imat,imat,ispin,isite)-trace
           norm=norm+dble(phi1(imat,imat,ispin,isite)*dconjg(phi1(imat,imat,ispin,isite)))
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

  do ieig=1,neig
     !********************************************
     !*** phi1 -> phi2=D*phi -> phi1=D^dag*phi2***
     !********************************************
     call Multiply_Dirac(temperature,&
          &xmat,alpha,phi1,phi2,GAMMA10d,nbmn,flux)
     !*********************************
     !*** adjust the margin and b.c.***
     !*********************************
     call Adjust_margin_and_bc_pf_2(phi2,nbc)
     
     call Multiply_Dirac_dagger(temperature,&
          &xmat,alpha,phi2,phi1,GAMMA10d,nbmn,flux)
     !*********************************
     !*** adjust the margin and b.c.***
     !*********************************
     call Adjust_margin_and_bc_pf_2(phi1,nbc)
     !************************************x
     !*** calculate the norm of phi1.  ***
     !************************************
     norm=0d0
     do imat=1,nmat     
        do jmat=1,nmat
           do ispin=1,nspin
              do isite=1,nsite
                 norm=norm&
                      &+dble(dconjg(phi1(imat,jmat,ispin,isite)&
                      &*dconjg(phi1(imat,jmat,ispin,isite))))
              end do
           end do
        end do
     end do
     norm=dsqrt(norm)
     norm=1d0/norm
     phi1=phi1*norm
     norm=1d0/norm
  end do
  
  largest_eig=norm

  return
  
END SUBROUTINE Largest_eigenvalue
