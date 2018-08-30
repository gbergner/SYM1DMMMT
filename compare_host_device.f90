!##############################################################################
!######              Code for direct comparison between               #########
!######                        host and device code                   #########
!######                   comparison only of the part that has the    #########
!######                   same functionality on host and device       #########
!######                 written by Georg Bergner                      #########
!######                                                               #########
!##############################################################################

module compare_host_device
contains
    ! These helpler functions are not completely necessary, but they
    ! simplify the comparison of the results since the order of the indices
    ! etc. are changed in some cases.
    SUBROUTINE difference_vect(vect1,vect2,difference)
        use compiletimeconstants
        implicit none
        double complex vect1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex vect2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double precision difference
        integer isite,jmat,imat,ispin
        difference=0.0d0
        do isite=1,nsite
            do jmat=1,nmat
                do imat=1,nmat
                    do ispin=1,nspin
                        difference=difference+abs(vect1(imat,jmat,ispin,isite)-vect2(imat,jmat,ispin,isite))
                    end do
                end do
            end do
        end do
    END SUBROUTINE difference_vect

    ! Reorders the indices.
    SUBROUTINE set_chid(chid,chi)
        use compiletimeconstants
        implicit none
        double complex :: chid(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:nremez_md,1:npf)
        double complex chi(1:nremez_md,1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
        integer ipf,iremez,isite,jmat,imat,ispin
        !$acc declare device_resident(chid)
        !$acc declare present(chi)
        !$acc kernels
        do ipf=1,npf
            do iremez=1,nremez_md
                do isite=-(nmargin-1),nsite+nmargin
                    do jmat=1,nmat
                        do imat=1,nmat
                            do ispin=1,nspin
                                chid(imat,jmat,ispin,isite,iremez,ipf)=chi(iremez,imat,jmat,ispin,isite,ipf)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$acc end kernels
    END SUBROUTINE set_chid

    ! Difference, taking into account the different index ordering.
    SUBROUTINE difference_chid(chid,chi,difference)
        use compiletimeconstants
        implicit none
        double complex :: chid(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:nremez_md,1:npf)
        double complex chi(1:nremez_md,1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
        integer ipf,iremez,isite,jmat,imat,ispin
        double precision difference
         !$acc declare device_resident(chid)
        !$acc declare present(chi)
        difference=0.0d0
        !$acc kernels
        do ipf=1,npf
            do iremez=1,nremez_md
                do isite=-(nmargin-1),nsite+nmargin
                    !do isite=1,nsite
                    do jmat=1,nmat
                        do imat=1,nmat
                            do ispin=1,nspin
                                difference=difference+abs(chid(imat,jmat,ispin,isite,iremez,ipf)-chi(iremez,imat,jmat,ispin,isite,ipf))
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !$acc end kernels
    END SUBROUTINE difference_chid

    ! Norm without boundary parts (can be done much simpler).
    SUBROUTINE norm_vect(vect1,norm)
        use compiletimeconstants
        implicit none
        double complex vect1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double precision norm
        integer isite,jmat,imat,ispin
        norm=0.0d0
        do isite=1,nsite
            do jmat=1,nmat
                do imat=1,nmat
                    do ispin=1,nspin
                        norm=norm+real(vect1(imat,jmat,ispin,isite)*dconjg(vect1(imat,jmat,ispin,isite)))
                    end do
                end do
            end do
        end do
    END SUBROUTINE norm_vect

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

    SUBROUTINE check_err(errin)
        implicit none
        double precision errin
        double precision, save :: overallerr=0
        integer, save ::  largeerrcounter=0
        if(errin>1e-8) then
            largeerrcounter=largeerrcounter+1
            write (*,*) "!!! WARNING LARGE ERROR (see below)"
        end if
        overallerr=overallerr+errin
        write (*,*) "Error of test: ", errin, " sum of errors: ", overallerr, " counter: ",   largeerrcounter
    END SUBROUTINE check_err
 
    ! The following functions are just standardized output functions.
    SUBROUTINE print_difference(vect1,vect2,difference)
        use compiletimeconstants
        implicit none
        double complex vect1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex vect2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double precision difference, norm1, norm2
        call difference_vect(vect1,vect2,difference)
        call norm_vect(vect1,norm1)
        call norm_vect(vect2,norm2)
        write (*,*) "Difference of vectors=",difference, " norm1=",norm1, " norm2=", norm2
    END SUBROUTINE print_difference

    SUBROUTINE print_difference_chid(chid,chi)
        use compiletimeconstants
        implicit none
        double complex :: chid(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:nremez_md,1:npf)
        double complex chi(1:nremez_md,1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
        double precision difference
        call difference_chid(chid,chi,difference)
        write (*,*) "Difference of chi=",difference
        call check_err(difference)
    END SUBROUTINE print_difference_chid

    SUBROUTINE difference_vect_device(vect1,vect2)
        use compiletimeconstants
        implicit none
        double complex vect1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex vect2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double precision difference
        integer isite,jmat,imat,ispin
        difference=0.0d0
        !$acc kernels &
        !$acc present (vect1,vect2)
        do isite=1,nsite
            do jmat=1,nmat
                do imat=1,nmat
                    do ispin=1,nspin
                        difference=difference+abs(vect1(imat,jmat,ispin,isite)-vect2(imat,jmat,ispin,isite))
                    end do
                end do
            end do
        end do
        !$acc end kernels
        write (*,*) "Difference of device vectors=",difference
    END SUBROUTINE difference_vect_device

    SUBROUTINE printtime(time1,time2)
        implicit none
        real :: time1,time2
        write(*,*) "-------TIMING------"
        print *, "Device time:", &
            time2, "seconds"
        print *, "Host time:", &
            time1, "seconds"
        print *,"relative timing",time1/time2
    END SUBROUTINE printtime


    ! Main comparison is done here.
    SUBROUTINE compare_host_device(xmat,alpha,ncv,n_bad_CG,nacceptance,nbc,nbmn,&
        &temperature,flux,GAMMA10d,ntau,nratio,dtau_xmat,dtau_alpha,&
        &acceleration,g_alpha,g_R,RCUT,&
        &acoeff_md,bcoeff_md,acoeff_pf,bcoeff_pf,max_err,max_iteration,iteration,&
        &ham_init,ham_fin,ntrial,imetropolis)
        use compiletimeconstants
        use dirac_operator
        use fourier_transform
        use cgm_solver
        use mtmod !Mersenne twistor
        use outputstreamnumbers
        use RHMC_Updater
        use hmc_molecular_dynamics
        use hmc_hamiltonean
        use lattice_action
        use hmc_force
        use utils_measurements
        implicit none
        ! This include is required for the host serial code.
    include 'Fourier.inc'

        ! Lots of variables are initialized here since we need sometimes separate host and
        ! device memory of the same object.
        integer nbc,nbmn
        double precision temperature,flux
        double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
        double precision max_err
        integer max_iteration
        double precision acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)
        double precision acoeff_pf(0:nremez_pf),bcoeff_pf(1:nremez_pf)
        double precision g_alpha,g_R,RCUT
        integer ntau,nratio
        integer imetropolis
        double precision acceleration(1:nsite),dtau_alpha,dtau_xmat
        !input & output
        double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double precision alpha(1:nmat)
        double precision alpha_d(1:nmat)
        integer ncv,n_bad_CG,nacceptance,ntrial
        double precision ham_init,ham_fin
        !output
        integer iteration


        double precision metropolis
        integer info_pf,info_mol,info_CG_init,info_CG_fin,info_alpha,info_accept,info
        double complex backup_xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double precision backup_alpha(1:nmat)
        double complex P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision P_alpha(1:nmat)
        double complex pf(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
        double complex Chi(1:nremez_md,1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:npf)

       !double complex,dimension(:,:,:,:,:,:),allocatable :: Chi
        double complex P_xmat_d(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision P_alpha_d(1:nmat)
        !$acc declare device_resident(P_xmat_d,P_alpha_d)
       
        integer IERR,myrank,nprocs
        double complex testvect1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex testvect2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex testvect0(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
  
        double complex :: testvect_d0(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex :: testvect_d2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex :: phase(1:nmat,1:nmat,1:2)
        double complex :: Gam123(1:nspin,1:nspin)
        double complex :: pf_d(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
        double complex :: Chi_d(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:nremez_md,1:npf)
        !$acc declare device_resident(testvect_d0,testvect_d2,phase,Gam123,pf_d,Chi_d)
  
        double complex :: Chi_test(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:nremez_md,1:npf)
  
        real :: start_time,stop_time,time1,time2
        integer :: counter
        double precision tmp1,tmp2
        double complex  :: tmp1compl
  
        double complex  :: delh_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision :: delh_alpha(1:nmat)
  
        ! separate device and host part to test device data on host.
        double complex  :: delh_xmat_d(1:nmat,1:nmat,1:ndim,1:nsite)
        double complex  :: delh_xmat_h(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision :: delh_alpha_d(1:nmat)
        double precision :: delh_alpha_h(1:nmat)
        !$acc declare device_resident(delh_xmat_d,delh_alpha_d)
  
        double complex :: P_xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
        double complex :: P_xmat_mom2(1:nmat,1:nmat,1:ndim,1:nsite)
        double complex :: P_xmat_mom_d(1:nmat,1:nmat,1:ndim,1:nsite)
         !$acc declare device_resident(P_xmat_mom_d)
   
        double complex:: xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
        double complex:: xmat_mom2(1:nmat,1:nmat,1:ndim,1:nsite)
        double complex:: xmat_d(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double complex:: xmat_d2(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double complex:: xmat_mom_d(1:nmat,1:nmat,1:ndim,1:nsite)
        !$acc declare device_resident(xmat_mom_d,xmat_d)
        !!$acc declare present(xmat_d)
  
        integer imat,jmat,isite,idim
        integer info_alpha_d

        !allocate(Chi(1:nremez_md,1:nmat,1:nmat,&
         !   1:nspin,-(nmargin-1):nsite+nmargin,1:npf))
        ! The first part is just included from the standard HMC on the host
        ! in order to arrive at a representative random configuration for the test.
        info_pf=1
        info_mol=1
        info_CG_init=1
        info_CG_fin=1
        call generate_pseudo_fermion_SUN(pf,xmat,alpha,&
            &GAMMA10d,acoeff_pf,bcoeff_pf,max_err,max_iteration,iteration,&
            &nbc,nbmn,temperature,flux,info_pf)
        backup_xmat=xmat
        backup_alpha=alpha
        call Generate_P_xmat(P_xmat)
        call Generate_P_alpha(P_alpha)
        call Molecular_Dynamics(nbc,temperature,&
            &ntau,nratio,dtau_xmat,dtau_alpha,xmat,alpha,P_xmat,P_alpha,&
            &acoeff_md,bcoeff_md,pf,5,max_err,iteration,&
            &gamma10d,g_alpha,g_R,RCUT,acceleration,nbmn,flux,info_mol)
       
        ! Now the real test part starts with a printout of the paramters.
        max_iteration=10000
        max_err=0.000000000000001d0
        write (*,*) "nbmn=", nbmn, " sum(xmat)=", sum(xmat), " sum(P_xmat)=",sum(P_xmat)," sum(alpha)=",sum(alpha)," sum(P_alpha)=",sum(P_alpha)
        testvect0=1.0d0
        call set_testvect(testvect0)
        testvect1=0.0d0
        testvect2=0.0d0

        testvect_d2=testvect2
        !$acc data &
        !$acc copyin(testvect0,temperature,xmat,alpha,GAMMA10d,nbmn,flux,nbc,bcoeff_md,bcoeff_pf)&
        !$acc copyin(Chi,pf,g_R,RCUT,P_xmat,P_alpha,acoeff_md,g_alpha,acceleration,alpha_d)&
        !$acc copy (testvect2,info,Chi_test,phase,Gam123,delh_alpha,delh_xmat)
        !call norm_vect_device(testvect_d0,testnorm)
  
        !$acc kernels
        tmp1=sum(xmat)
        !$acc end kernels
        ! check basic communications.
        write (*,*)  " sum(xmat)_host=", sum(xmat)," device=",tmp1
        tmp1=abs(tmp1-sum(xmat))
        call check_err(tmp1)
        ! check multiplications
        call Multiply_Dirac(temperature,xmat,alpha,&
            testvect0,testvect1,GAMMA10d,nbmn,flux)
        write (*,*) "host multiply"
        !$acc kernels
        testvect_d0=testvect0
        !$acc end kernels
        call  setup_data_device(alpha,flux,GAMMA10d,phase,Gam123,temperature)
        call Multiply_Dirac_device(temperature,xmat,phase,Gam123,nbmn,&
            testvect_d0,testvect_d2)
        !$acc kernels
        testvect2=testvect_d2
        !$acc end kernels
        !$acc update host(testvect2)
        write (*,*) "device multiply"
        write (*,*) "-----> device multiply host data update"
        call print_difference(testvect1,testvect2,tmp1)
        call check_err(tmp1)
        call Multiply_Dirac_dagger(temperature,xmat,alpha,&
            testvect0,testvect1,GAMMA10d,nbmn,flux)
        call Multiply_Dirac_dagger_device(temperature,xmat,phase,Gam123,nbmn,&
            testvect_d0,testvect_d2)
        !$acc kernels
        testvect2=testvect_d2
        !$acc end kernels
        !$acc update host(testvect2)
        write (*,*) "------> device multiply host data update dagger"
        call print_difference(testvect1,testvect2,tmp1)
        call check_err(tmp1)
        write(*,*) "-----------------------------------------------------------"
        print *,"Testing CGM solver"
        !$acc kernels
        pf_d=pf
        !$acc end kernels
        call set_chid(Chi_d,chi)
        call cpu_time(start_time)
        call cgm_solver_device(nremez_md,bcoeff_md,nbmn,nbc,temperature,&
            max_err,max_iteration,xmat,phase,Gam123,pf_d,chi_d,info,iteration)
        call cpu_time(stop_time)
        time1=stop_time - start_time
        call cpu_time(start_time)
        call solver_biCGm(nbc,nbmn,nremez_md,&
            xmat,alpha,pf,chi,GAMMA10d,&
            bcoeff_md,max_err,max_iteration,iteration,&
            temperature,flux,info_CG_init)
        call cpu_time(stop_time)
        time2=stop_time - start_time
        !this time different order:
        call printtime(time2, time1)
        print*,'CGM info',info_CG_init, " ",info
        !$acc update device(chi)
        call print_difference_chid(chi_d,chi)
        !$acc kernels
        tmp1compl=Sum(chi_d)
        !$acc end kernels
        print*,"sum chi ",Sum(chi)," ",tmp1compl
        tmp1=abs(tmp1compl-Sum(chi))
        call check_err(tmp1)
        write(*,*) "-----------------------------------------------------------"
        !$acc kernels
        xmat_d=xmat
        alpha_d=alpha
        !$acc end kernels
        !This sets the host part to zero in order to ensure that the device part is used.
        xmat_d=(0d0,0d0)
        alpha_d=0d0
        call cpu_time(start_time)
        call Calc_action(temperature,xmat,alpha,tmp1)
        call cpu_time(stop_time)
        time1=stop_time-start_time
        call cpu_time(start_time)
        call Calc_action_device(temperature,xmat_d,alpha_d,tmp2,phase)
        call cpu_time(stop_time)
        time2=stop_time-start_time
        print *,"difference action calculation", abs(tmp1-tmp2), " ", tmp1, " ",tmp2
        call printtime(time1,time2)
        tmp1=abs(tmp1-tmp2)
        call check_err(tmp1)
        write(*,*) "-----------------------------------------------------------"
        call cpu_time(start_time)
        call Calc_Ham(temperature,xmat,alpha,P_xmat,P_alpha,tmp1,pf,chi,&
            &acoeff_md,g_R,RCUT,nbmn,flux)
        call cpu_time(stop_time)
        time1=stop_time-start_time
        !$acc kernels
        P_xmat_d=P_xmat
        !$acc end kernels
        !$acc kernels
        P_alpha_d=P_alpha
        !$acc end kernels
        call cpu_time(start_time)
        call Calc_Ham_device(temperature,xmat_d,alpha_d,P_xmat_d,P_alpha_d,tmp2,pf,Chi_d,&
            &acoeff_md,g_R,RCUT,nbmn,flux,phase)
        call cpu_time(stop_time)
        time2=stop_time-start_time
        print *,"difference ham calculation: ", abs(tmp1-tmp2), " ", tmp1, " ",tmp2
        call printtime(time1,time2)
        tmp1=abs(tmp1-tmp2)
        call check_err(tmp1)
        write(*,*) "-----------------------------------------------------------"
        call cpu_time(start_time)
        call Calc_Force_bosonic(delh_xmat,delh_alpha,xmat,alpha,chi,&
            GAMMA10d,g_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md)
        call cpu_time(stop_time)
        time1=stop_time-start_time
        call cpu_time(start_time)
        call Calc_Force_bosonic_device(delh_xmat_d,delh_alpha_d,xmat_d,alpha_d,Chi_d,&
            GAMMA10d,g_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md)
        call cpu_time(stop_time)
        time2=stop_time-start_time
        !$acc kernels
        delh_alpha_h=delh_alpha_d
        delh_xmat_h=delh_xmat_d
         !$acc end kernels
        print *,"difference of bosonic force alpha calculation ",Sum(abs(delh_alpha-delh_alpha_h))," ",Sum(abs(delh_alpha))," ",&
            & Sum(abs(delh_alpha_h))," ",Maxval(delh_alpha)-Maxval(delh_alpha_h)
        print *,"difference of bosonic force xmat calculation ",Sum(abs(delh_xmat-delh_xmat_h))," ",Sum(abs(delh_xmat))," ", Sum(abs(delh_xmat_h))
        call printtime(time1,time2)
        tmp1=Sum(abs(delh_alpha-delh_alpha_h))+Sum(abs(delh_xmat-delh_xmat_h))
        call check_err(tmp1)
        write(*,*) "-----------------------------------------------------------"
            !$acc kernels
        delh_alpha_d=0d0
        delh_xmat_d=0d0
        !$acc end kernels
        call cpu_time(start_time)
        call Calc_Force_fermionic(delh_xmat,delh_alpha,xmat,alpha,chi,&
            GAMMA10d,g_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md)
        call cpu_time(stop_time)
        time1=stop_time-start_time
        call cpu_time(start_time)
        call Add_Force_fermionic_device(1d0,delh_xmat_d,delh_alpha_d,xmat_d,chi_d,&
            g_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md,phase,Gam123,nbc)
        call cpu_time(stop_time)
        time2=stop_time-start_time
        !$acc kernels
        delh_alpha_h=delh_alpha_d
        delh_xmat_h=delh_xmat_d
         !$acc end kernels
        print *,"difference of fermionic force alpha calculation (add) ",Sum(abs(delh_alpha-delh_alpha_h))," ",Sum(abs(delh_alpha))," ", Sum(abs(delh_alpha_h))
        print *,"difference of fermonic force xmat calculation (add) ",Sum(abs(delh_xmat-delh_xmat_h))," ",Sum(delh_xmat)," ", Sum(delh_xmat_h)
        call printtime(time1,time2)
        tmp1=Sum(abs(delh_alpha-delh_alpha_h))+Sum(abs(delh_xmat-delh_xmat_h))
        call check_err(tmp1)
        do imat=1,nmat
          print*,"i ",imat," ",delh_alpha(imat)," ",delh_alpha_h(imat)
        end do
        write(*,*) "-----------------------------------------------------------"
        call cpu_time(start_time)
        call Fourier_transform_P_xmat(P_xmat,P_xmat_mom,1)
        call cpu_time(stop_time)
        time1=stop_time-start_time
        !$acc kernels
        P_xmat_d=P_xmat
        !$acc end kernels
        call cpu_time(start_time)
        call FT_P_xmat_device(P_xmat_d,P_xmat_mom_d)
        call cpu_time(stop_time)
        time2=stop_time-start_time
        !$acc kernels
        P_xmat_mom2=P_xmat_mom_d
        !$acc end kernels
        print *,"difference of FT P_xmat ",Sum(abs(P_xmat_mom-P_xmat_mom2))," ",Sum(abs(P_xmat_mom))," ", Sum(abs(P_xmat_mom2))
        call printtime(time1,time2)
        tmp1=Sum(abs(P_xmat_mom-P_xmat_mom2))
        call check_err(tmp1)
        write(*,*) "-----------------------------------------------------------"
        call Fourier_transform_P_xmat(P_xmat_mom2,P_xmat_mom,2)
        call cpu_time(stop_time)
        time1=stop_time-start_time
        !$acc kernels
        P_xmat_d=P_xmat_mom
        !$acc end kernels
        call cpu_time(start_time)
        call FTinv_P_xmat_device(P_xmat_mom_d,P_xmat_d)
        call cpu_time(stop_time)
        time2=stop_time-start_time
        !$acc kernels
        P_xmat_mom=P_xmat_mom_d
        !$acc end kernels
        print *,"difference of inv FT P_xmat ",Sum(abs(P_xmat_mom2-P_xmat_mom))," ",Sum(abs(P_xmat_mom2))," ", Sum(abs(P_xmat_mom))
        call printtime(time1,time2)
        tmp1=Sum(abs(P_xmat_mom2-P_xmat_mom))
        call check_err(tmp1)
        write(*,*) "-----------------------------------------------------------"
        call cpu_time(start_time)
        call Fourier_transform_xmat(xmat,xmat_mom,1)
        call cpu_time(stop_time)
        time1=stop_time-start_time
        !$acc kernels
        xmat_d=xmat
        !$acc end kernels
        call cpu_time(start_time)
        call FT_xmat_device(xmat_d,xmat_mom_d)
        call cpu_time(stop_time)
        time2=stop_time-start_time
        !$acc kernels
        xmat_mom2=xmat_mom_d
        !$acc end kernels
        print *,"difference of FT xmat ",Sum(abs(xmat_mom-xmat_mom2))," ",Sum(abs(xmat_mom))," ", Sum(abs(xmat_mom2))
        call printtime(time1,time2)
        tmp1=Sum(abs(xmat_mom-xmat_mom2))
        call check_err(tmp1)
        write(*,*) "-----------------------------------------------------------"
        call cpu_time(start_time)
        call Fourier_transform_xmat(xmat,xmat_mom,2)
        call cpu_time(stop_time)
        time1=stop_time-start_time
        !$acc kernels
        xmat_d=xmat
        !$acc end kernels
        call cpu_time(start_time)
        call FTinv_xmat_device(xmat_d,xmat_mom_d)
        call cpu_time(stop_time)
        time2=stop_time-start_time
        !$acc kernels
        xmat_mom2=xmat_mom_d
        !$acc end kernels
        print *,"difference of FT xmat ",Sum(abs(xmat_mom-xmat_mom2))," ",Sum(abs(xmat_mom))," ", Sum(abs(xmat_mom2))
        call printtime(time1,time2)
        tmp1=Sum(abs(xmat_mom-xmat_mom2))
        call check_err(tmp1)
        write(*,*) "-----------------------------------------------------------"
        !$acc kernels
        xmat_d=xmat
        !$acc end kernels
        print *,"before",Sum(abs(xmat))
        call cpu_time(start_time)
        call Adjust_margin_xmat(xmat)
        call cpu_time(stop_time)
        time1=stop_time-start_time
        print *,"after",Sum(abs(xmat))
        call cpu_time(start_time)
        call Adjust_margin_xmat_device(xmat_d)
        call cpu_time(stop_time)
        time2=stop_time-start_time
        !$acc kernels
        tmp1=Sum(abs(xmat_d))
        !$acc end kernels
        print *,"difference of AdjustMargins xmat ",abs(Sum(abs(xmat))-tmp1)," ",Sum(abs(xmat))," ", tmp1
        call printtime(time1,time2)
        tmp1=abs(Sum(abs(xmat))-tmp1)
        call check_err(tmp1)
        write(*,*) "-----------------------------------------------------------"
  
        xmat_d=xmat
        alpha_d=alpha
        !$acc update device(xmat_d,alpha_d)
        xmat_d=0d0
        alpha_d=0d0
        !$acc kernels
        P_xmat_d=P_xmat
        P_alpha_d=P_alpha
        !$acc end kernels
        !$acc kernels
        pf_d=pf
        !$acc end kernels
        print *, "dt ",dtau_xmat," ",dtau_alpha," ",nratio
        call cpu_time(start_time)
        call Molecular_Dynamics(nbc,temperature,&
            &1,nratio,dtau_xmat,dtau_alpha,xmat,alpha,P_xmat,P_alpha,&
            &acoeff_md,bcoeff_md,pf,max_iteration,max_err,iteration,&
            &gamma10d,g_alpha,g_R,RCUT,acceleration,nbmn,flux,info_mol)
        call cpu_time(stop_time)
        time1=stop_time-start_time
        call cpu_time(start_time)
        call Molecular_Dynamics_device(nbc,temperature,&
            &1,nratio,dtau_xmat,dtau_alpha,xmat_d,alpha_d,phase,P_xmat_d,P_alpha_d,&
            &acoeff_md,bcoeff_md,pf,max_iteration,max_err,iteration,&
            &gamma10d,gam123,g_alpha,g_R,RCUT,acceleration,nbmn,flux,info_mol)
        call cpu_time(stop_time)
        time2=stop_time-start_time
        !$acc update host(xmat_d,alpha_d)
        print *,"difference of xmat in molecular dynamics ",Sum(abs(xmat-xmat_d))," ",Sum(xmat)," ", Sum(xmat_d)," ",Sum(abs(xmat))
        print *,"difference of alpha in molecular dynamics ",Sum(alpha-alpha_d)," ",Sum(alpha)," ", Sum(alpha_d)
        call printtime(time1,time2)
        tmp1=Sum(abs(xmat-xmat_d))+Sum(alpha-alpha_d)
        call check_err(tmp1)
        write(*,*) "-----------------------------------------------------------"
        xmat_d=xmat
        alpha_d=alpha
        !$acc update device(xmat_d,alpha_d)
        call cpu_time(start_time)
        call subtract_U1(xmat,alpha)
        call cpu_time(stop_time)
        time1=stop_time-start_time
        call cpu_time(start_time)
        call subtract_U1_device(xmat_d,alpha_d)
        call cpu_time(stop_time)
        time2=stop_time-start_time
        !$acc update host(xmat_d,alpha_d)
        print *,"difference of xmat in subtr ",Sum(abs(xmat-xmat_d))," ",Sum(abs(xmat))," ", Sum(abs(xmat_d))
        print *,"difference of alpha in subtr ",Maxval(alpha-alpha_d)," ",Sum(alpha)," ", Sum(alpha_d)
        call printtime(time1,time2)
        tmp1=Sum(abs(xmat-xmat_d))+Maxval(abs(alpha-alpha_d))
        call check_err(tmp1)
        write(*,*) "-----------------------------------------------------------"
        xmat_d=xmat
        alpha_d=alpha
        !$acc update device(xmat_d,alpha_d)
        call cpu_time(start_time)
        call check_alpha_constraint(alpha,info_alpha)
        call cpu_time(stop_time)
        time1=stop_time-start_time
        call cpu_time(start_time)
        call check_alpha_constraint(alpha_d,info_alpha_d)
        call subtract_U1_device(xmat_d,alpha_d)
        call cpu_time(stop_time)
        time2=stop_time-start_time
        !$acc update host(xmat_d,alpha_d)
        print *,"difference of check alpha constraint ",info_alpha," ", info_alpha_d
        call printtime(time1,time2)
        tmp1=abs(info_alpha-info_alpha_d)
        call check_err(tmp1)
        write(*,*) "-----------------------------------------------------------"
        !$acc end data
 
        return

    END SUBROUTINE compare_host_device

END MODULE compare_host_device
