subroutine Calc_energy(temperature,xmat,alpha,energy,myrank,nbmn,flux,&
     &acoeff_md,bcoeff_md,acoeff_pf,bcoeff_pf,&
     &nbc,max_err,max_iteration)

  implicit none

  include 'mpif.h'
  include 'size_parallel.h'
  !***** input *****
  integer myrank,nbmn
  double precision temperature,flux
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
       -(nmargin-1):nsite_local+nmargin)
  double precision alpha(1:nmat_block*nblock)
  double precision acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)
  double precision acoeff_pf(0:nremez_pf),bcoeff_pf(1:nremez_pf)
  !***** output *****
  double precision energy
  !******************
  double precision action,kinetic,potential,potential_BMN,energy_local
  double precision lattice_spacing

  double complex commutator(1:nmat_block,1:nmat_block)
  double complex uxumx(1:nmat_block,1:nmat_block)
  integer isite!,isite_p1
  integer idim,jdim
  integer imat,jmat,kmat
  double complex ei,ej

  double complex xmat_row(1:nmat_block,1:nmat_block*nblock,1:ndim,1:nsite_local)
  double complex xmat_column(1:nmat_block*nblock,1:nmat_block,&
       &1:ndim,1:nsite_local)
  double complex x23(1:nmat_block,1:nmat_block),&
       x32(1:nmat_block,1:nmat_block),&
       trx123,trx132
  double precision trx2_123,trx2_456789

  integer iblock,jblock,isublat
  integer info,info_pf
  double complex pf(1:nmat_block,1:nmat_block,&
       1:nspin,-(nmargin-1):nsite_local+nmargin)
  double complex Chi(1:nremez_md,1:nmat_block,1:nmat_block,&
       1:nspin,-(nmargin-1):nsite_local+nmargin)
  double complex Mchi(1:nremez_md,1:nmat_block,1:nmat_block,&
       &1:nspin,1:nsite_local)
  
  double precision sum_pf
  double complex Gamma10d(1:ndim,1:nspin,1:nspin),gam123(1:nspin,1:nspin),gam12(1:nspin,1:nspin)
  integer ispin,jspin,kspin
  integer iremez
  integer nbc,max_iteration,iteration
  double precision max_err
  !***** For MPI *****
  integer IERR

  call who_am_i(myrank,isublat,iblock,jblock)
  !move i-th row and j-th row of xmat to (i,j)-th node.
  call mpi_xmat_row(xmat,xmat_row,myrank)
  call mpi_xmat_column(xmat,xmat_column,myrank)
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
                         &+xmat_row(imat,kmat,idim,isite)&
                         &*xmat_column(kmat,jmat,jdim,isite)&
                         &-xmat_row(imat,kmat,jdim,isite)&
                         &*xmat_column(kmat,jmat,idim,isite)
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
  potential=potential*0.75d0*dble(nmat_block*nblock)*lattice_spacing
  !******************************  
  !*** Plane wave deformation ***
  !******************************
  potential_BMN=0d0
  sum_pf=0d0
  if(nbmn.EQ.1)then
     !*****************
     !*** mass term ***
     !*****************
     trx2_123=0d0
     trx2_456789=0d0
     do isite=1,nsite_local
        do imat=1,nmat_block
           do jmat=1,nmat_block
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
     potential_BMN=potential_BMN+flux*flux*(trx2_123+0.25d0*trx2_456789)&
          *lattice_spacing*dble(nmat_block*nblock)
     !******************
     !*** cubic term ***
     !******************
     trx123=(0d0,0d0)
     trx132=(0d0,0d0)
     do isite=1,nsite_local
        x23=(0d0,0d0)
        x32=(0d0,0d0)
        do imat=1,nmat_block
           do jmat=1,nmat_block
              do kmat=1,nmat_block*nblock
                 x23(imat,jmat)=x23(imat,jmat)&
                      +xmat_row(imat,kmat,2,isite)&
                      *xmat_column(kmat,jmat,3,isite)
                 x32(imat,jmat)=x32(imat,jmat)&
                      +xmat_row(imat,kmat,3,isite)&
                      *xmat_column(kmat,jmat,2,isite)
              end do
           end do
        end do
        do imat=1,nmat_block
           do jmat=1,nmat_block
              trx123=trx123+x23(imat,jmat)*dconjg(xmat(imat,jmat,1,isite))
              trx132=trx132+x32(imat,jmat)*dconjg(xmat(imat,jmat,1,isite))
           end do
        end do
     end do
     potential_BMN=potential_BMN&
          &+dble((0d0,7.5d0)*(trx123-trx132))*flux*lattice_spacing&
          &*dble(nmat_block*nblock)

!!$
!!$     
!!$     call MakeGamma(Gamma10d)
!!$     call generate_pseudo_fermion_SUN(nbc,acoeff_pf,bcoeff_pf,temperature,&
!!$          &xmat,alpha,pf,GAMMA10d,max_err,max_iteration,iteration,myrank,&
!!$          &nbmn,flux,info_pf)
!!$     !info_pf=0 -> OK (CG solver converged)
!!$     !info_pf=1 -> error (CG solver did not converge)
!!$     !Calculate Chi_k = (D+bcoeff_md(k))^{-1}*pf by using multi-mass CG-solver
!!$     call solver_biCGm(nbc,nremez_md,bcoeff_md,temperature,&
!!$          &xmat,alpha,pf,chi,GAMMA10d,max_err,max_iteration,iteration,myrank,&
!!$          &nbmn,flux,info)
!!$
!!$     !Multiply (-i)*M
!!$     call Multiply_Dirac_to_chi(temperature,xmat_row,xmat_column,alpha,&
!!$          chi,mchi,GAMMA10d,nbmn,flux,myrank)
!!$
!!$     !Gamma10d is actually (-i)*gamma
!!$     !gam123 is actually (-i)^3*gamma_{123}
!!$     gam12=(0d0,0d0)
!!$     do ispin=1,nspin
!!$        do jspin=1,nspin
!!$           do kspin=1,nspin
!!$              gam12(ispin,jspin)=gam12(ispin,jspin)&
!!$                   &+Gamma10d(1,ispin,kspin)*Gamma10d(2,kspin,jspin)
!!$           end do
!!$        end do
!!$     end do
!!$     gam123=(0d0,0d0)
!!$     do ispin=1,nspin
!!$        do jspin=1,nspin
!!$           do kspin=1,nspin
!!$              gam123(ispin,jspin)=gam123(ispin,jspin)&
!!$                   &+gam12(ispin,kspin)*Gamma10d(3,kspin,jspin)
!!$           end do
!!$        end do
!!$     end do
!!$    
!!$     do iremez=1,nremez_md
!!$        do ispin=1,nspin        
!!$           do jspin=1,nspin
!!$              do imat=1,nmat_block
!!$                 do jmat=1,nmat_block
!!$                    do isite=1,nsite_local
!!$                       sum_pf=sum_pf+&
!!$                            &acoeff_md(iremez)&
!!$                            &*dble(dconjg(Chi(iremez,imat,jmat,ispin,isite))&
!!$                            &*gam123(ispin,jspin)&
!!$                            &*mchi(iremez,imat,jmat,jspin,isite)*(0d0,1d0))
!!$                    end do
!!$                 end do
!!$              end do
!!$           end do
!!$        end do
!!$     end do
!!$     !sum_pf=sum_pf*(-1.5d0*dble(nmat_block*nblock)*flux*temperature)
!!$     !Be careful about the normalization of pseudo fermion.
!!$     sum_pf=sum_pf*(-1.5d0*flux*temperature)*lattice_spacing
!!$     
     end if
!!$  !The gauge fixing term is not needed when we calculate the energy.
!!$  action=kinetic+potential+potential_BMN
!!$
!!$  energy_local=-3d0*temperature*action&
!!$       &/dble(nmat_block*nmat_block*nblock*nblock)&
!!$       &+1.5d0*temperature*dble(ndim*nsite_local)/dble(nblock*nblock)
!!$
!!$  !**************************************
!!$  !******** The sign was wrong!!  *******
!!$  sum_pf=-sum_pf
!!$  !**************************************
!!$  !**************************************
!!$  
  
  energy_local=energy_local+sum_pf/dble(nmat_block*nmat_block*nblock*nblock)

  call MPI_Reduce(energy_local,energy,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,0,MPI_COMM_WORLD,IERR)

  !after summing up contributions from all ranks, 
  !we must subtract the U(1) part, 1.5d0*temperature*dble(ndim)/dble(nmat*nmat)
  if(myrank.EQ.0)then

     energy=energy-1.5d0*temperature*dble(ndim)&
          /dble(nmat_block*nmat_block*nblock*nblock)
  end if

  return

END subroutine Calc_energy
