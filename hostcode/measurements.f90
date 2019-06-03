SUBROUTINE measurements(xmat,alpha,nbc,nbmn,temperature,flux,&
     &GAMMA10d,neig_max,neig_min,ham_init,ham_fin,itraj,ntrial,iteration,&
     &max_err,max_iteration,ncv,n_bad_CG,nacceptance,ngauge,purebosonic)

  implicit none
  include '../staticparameters.f90'
  include '../unit_number.inc'
  !input
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double precision alpha(1:nmat)
  integer nbc,nbmn,ntrial,iteration,itraj,nacceptance,ngauge,purebosonic
  integer ncv,n_bad_CG
  integer neig_max,neig_min
  double precision temperature,flux
  double complex Gamma10d(1:ndim,1:nspin,1:nspin)
  double complex eigval(1:ndim)
  double precision max_err
  integer max_iteration
  integer ivec
  double precision trx2,com2,Pol,largest_eig,smallest_eig,energy,ham_fin,ham_init


  !measurements. 
  call Calc_TrX2(xmat,trx2)
  call Calc_TrX2_eigenvalues(xmat,eigval)
  call Calc_Com2(xmat,com2)
  call Calc_energy(temperature,xmat,alpha,energy,nbmn,flux)
  call Calc_Polyakov(alpha,Pol)
  !largest eigenvalue of D^dag*D. 
  if(neig_max.GT.0)then
     call Largest_eigenvalue(temperature,&
          &xmat,alpha,GAMMA10d,neig_max,largest_eig,nbc,nbmn,flux,&
          &nsite,nmat,nremez_md,nremez_pf,nimprove,ngauge,nmargin)
  end if
  !smallest eigenvalue of D^dag*D. 
  if(neig_min.GT.0)then
     call Smallest_eigenvalue(temperature,&
          &xmat,alpha,GAMMA10d,neig_min,smallest_eig,nbc,&
          &max_err,max_iteration,nbmn,flux,nsite,nmat,nremez_md,&
          &nremez_pf,nimprove,ngauge,nmargin)
  end if

  !**************
  !*** output ***
  !**************
  if((neig_max.EQ.0).AND.(neig_min.EQ.0))then
     
     write(unit_measurement,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,5(1x,f15.7))')&
          &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
          &Pol,trx2,com2,dble(nacceptance)/dble(ntrial)
     
     write(*,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,5(1x,f15.7))')&
          &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
          &Pol,trx2,com2,dble(nacceptance)/dble(ntrial)
     
  else  if((neig_max.GT.0).AND.(neig_min.EQ.0))then
     
     write(unit_measurement,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,6(1x,f15.7))')&
          &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
          &Pol,trx2,com2,largest_eig,&
          &dble(nacceptance)/dble(ntrial)
     write(*,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,6(1x,f15.7))')&
             &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
             &Pol,trx2,com2,largest_eig,&
             &dble(nacceptance)/dble(ntrial)
     
  else  if((neig_max.EQ.0).AND.(neig_min.GT.0))then
     
     write(unit_measurement,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,6(1x,f15.7))')&
          &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
          &Pol,trx2,com2,smallest_eig,&
          &dble(nacceptance)/dble(ntrial)
     write(*,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,6(1x,f15.7))')&
          &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
          &Pol,trx2,com2,smallest_eig,&
          &dble(nacceptance)/dble(ntrial)
     
  else  if((neig_max.GT.0).AND.(neig_min.GT.0))then
     write(unit_measurement,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,7(1x,f15.7))')&
          &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
          &abs(Pol),trx2,com2,largest_eig,smallest_eig,&
          &dble(nacceptance)/dble(ntrial)
     write(*,'(I8,1x,f15.9,1x,I4,1x,I8,1x,I8,7(1x,f15.7))')&
          &itraj,ham_fin-ham_init,ncv,n_bad_CG,iteration,energy,&
          &abs(Pol),trx2,com2,largest_eig,smallest_eig,&
          &dble(nacceptance)/dble(ntrial)
  end if
  

  write(unit_Polyakov_phase,"(I)",advance="no") itraj
  do ivec=1,nmat
    write(unit_Polyakov_phase,"(G)",advance="no") alpha(ivec)
  end do
  write(unit_Polyakov_phase,*) ""

  write(unit_Eigenval,"(I)",advance="no") itraj
  do ivec=1,ndim
    write(unit_Eigenval,"(G)",advance="no") dble(eigval(ivec))
  end do
  write(unit_Eigenval,*) ""
  ! print out the imag part, intendet for checks
  !write(unit_Eigenval,*) eigval

  flush(unit_measurement)
  flush(unit_Polyakov_phase)
  flush(unit_Eigenval)

  return 

END SUBROUTINE measurements
