!##############################################################################
!######              BFSS matrix model on lattice                     #########
!######              Benchmark for CPU and GPU multiplications        #########
!######              written by Georg Bergner                         #########
!##############################################################################

! Setting some more or less random test vectors.
SUBROUTINE set_testvect(vect1)
    use compiletimeconstants
    implicit none
    double complex,intent(inout) :: vect1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
    integer isite,jmat,imat,ispin
    double precision tmp1,tmp2
    do isite=1,nsite
        do jmat=1,nmat
            do imat=1,nmat
                do ispin=1,nspin
                    tmp1=0.02*isite+0.1*(imat+jmat)
                    tmp2=ispin
                    vect1(imat,jmat,ispin,isite)=tmp1+tmp2
                end do
            end do
        end do
    end do
END SUBROUTINE set_testvect


program test_time_inverter
    use compiletimeconstants
    use dirac_operator
    use cublasinterface
    use cgm_solver
    implicit none
  
    double complex :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    double precision alpha(1:nmat)

    !Gamma matrices
    double complex Gamma10d(1:ndim,1:nspin,1:nspin)
    double precision temperature
    double precision flux

    integer :: nbc !boundary condition for fermions; 0 -> pbc, 1 -> apbc
    integer:: nbmn ! 0 -> BFSS, 1 -> BMN
  
    double complex :: phase(1:nmat,1:nmat,1:2)
    double complex :: Gam123(1:nspin,1:nspin)
    !$acc declare device_resident(phase,Gam123)
  
    double complex :: testvect0(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
    double complex :: testvect_d0(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
    double complex :: chi_d(1:nmat,1:nmat,&
        1:nspin,-(nmargin-1):nsite+nmargin,1:nremez_md,1:npf)
    !$acc declare device_resident(testvect_d0,chi_d)

    real :: start_time,stop_time,time1,time2,time1_cud,time2_cud
    integer :: isite,idim,imat,jmat,ipf
    double precision bcoeff_md(1:nremez_md),acoeff_md(0:nremez_md)
    double precision max_err
    integer max_iteration,info,iteration

    if(npf.EQ.1)then
     include 'remez_md_1.dat'
    else if(npf.EQ.2)then
     include 'remez_md_2.dat'
    end if

    nbc=1 !boundary condition for fermions; 0 -> pbc, 1 -> apbc
    nbmn=0 ! 0 -> BFSS, 1 -> BMN
    temperature=1d0*nbc
    flux=0.5d0*temperature
     
    call MakeGamma(Gamma10d)
 


    ! more or less random field xmat
    xmat=(0d0,0d0)
    do isite=-(nmargin-1),nsite+nmargin
        do idim=1,ndim
            do imat=1,nmat
                do jmat=1,nmat
                    xmat(imat,jmat,idim,isite)=CMPLX(imat+2.943*jmat,idim+0.02*isite)
                    xmat(jmat,imat,idim,isite)=conjg(xmat(imat,jmat,idim,isite))
                end do
            end do
        end do
    end do
    do imat=1,nmat
      alpha(imat)=imat*0.01
    end do

    do ipf=1,npf
      call set_testvect(testvect0(:,:,:,:,ipf))
    end do
    max_iteration=10000
    max_err=0.000000000000001d0
     print*, "vect in ", Sum(testvect0)
    !$acc data &
    !$acc copyin(testvect0,temperature,xmat,alpha,GAMMA10d,nbmn,flux,nbc,nbmn,bcoeff_md)


    call cudaSetupTimer()
    !$acc kernels
    testvect_d0=testvect0
    !$acc end kernels
    call  setup_data_device(alpha,flux,GAMMA10d,phase,Gam123,temperature)
    call cudaTimerStart()
    call cpu_time(start_time)
    call cgm_solver_device(nremez_md,bcoeff_md,nbmn,nbc,temperature,&
        max_err,max_iteration,xmat,phase,Gam123,testvect_d0,chi_d,info,iteration)
    call cpu_time(stop_time)
    call cudaTimerStop(time2_cud)
    time2=stop_time - start_time
    print *, "Device time:", &
        time2,"(",time2_cud, ")seconds ", iteration, " iterations"

    call cudaFinishTimer()
  !End test part
  !deallocate(xmat)
      !$acc end data
end program test_time_inverter
