! alpha -> alpha + P_alpha*dtau_alpha
! xmat -> xmat + P_xmat*dtau_xmat
! P_alpha -> P_alpha - delh_alpha*dtau_alpha
! P_xmat -> P_xmat - delh_xmat*dtau_xmat
! delh_xmat(imat,jmat)=dH/dxmat(jmat,imat)
MODULE HMC_Molecular_Dynamics
contains

    SUBROUTINE Field_step(xmat_mom,P_xmat_mom,alpha,P_alpha,dtau_xmat,dtau_alpha,acceleration)

        use compiletimeconstants
        implicit none

        double complex,intent(inout) :: xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
        double complex,intent(in) :: P_xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision, intent(inout):: alpha(1:nmat)
        double precision, intent(in):: P_alpha(1:nmat)
        double precision, intent(in) :: acceleration(1:nsite)
        double precision, intent(in) :: dtau_xmat,dtau_alpha
        integer :: imom,idim,imat,jmat
         !$acc declare device_resident(xmat_mom,P_xmat_mom,P_alpha)
        !$acc declare present(acceleration,alpha)
        !$acc data copyin(dtau_xmat,dtau_alpha)
        !$acc kernels
        do imom=1,nsite
            do idim=1,ndim
                do imat=1,nmat
                    do jmat=1,nmat
                        xmat_mom(imat,jmat,idim,imom)=&
                            &xmat_mom(imat,jmat,idim,imom)&
                            &+P_xmat_mom(imat,jmat,idim,imom)&
                            &*dcmplx(dtau_xmat*acceleration(imom))
                    end do
                end do
            end do
        end do
        !$acc end kernels
        !$acc kernels
        alpha=alpha+P_alpha*dtau_alpha
      !$acc end kernels
      !$acc end data
    END SUBROUTINE Field_step


    SUBROUTINE Momentum_step(P_xmat_mom,xmat_mom,P_alpha,alpha,dtau_xmat,dtau_alpha,acceleration)

        use compiletimeconstants
        implicit none

        double complex,intent(in) :: xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
        double complex,intent(inout) :: P_xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision, intent(in):: alpha(1:nmat)
        double precision, intent(inout):: P_alpha(1:nmat)
        double precision, intent(in) :: acceleration(1:nsite)
        double precision, intent(in) :: dtau_xmat,dtau_alpha
        integer :: imom,idim,imat,jmat
        !$acc declare device_resident(xmat_mom,P_xmat_mom,P_alpha)
        !$acc declare present(alpha)
        !$acc declare present(acceleration)
        !$acc data copyin(dtau_xmat,dtau_alpha)
        !$acc kernels
        P_alpha=P_alpha-alpha*dtau_alpha
        !$acc end kernels
        !$acc kernels
        do imom=1,nsite
            do idim=1,ndim
                do jmat=1,nmat
                    do imat=1,nmat
                        P_xmat_mom(imat,jmat,idim,imom)=&
                            &P_xmat_mom(imat,jmat,idim,imom)&
                            &-xmat_mom(imat,jmat,idim,imom)&
                            &*dcmplx(dtau_xmat*acceleration(imom))
                    end do
                end do
            end do
        end do
         !$acc end kernels
         !$acc end data
    END SUBROUTINE Momentum_step


    subroutine Molecular_Dynamics_device(nbc,temperature,&
        &ntau,nratio,dtau_xmat,dtau_alpha,xmat,alpha,phase,P_xmat,P_alpha,&
        &acoeff_md,bcoeff_md,pf,max_iteration,max_err,&
        &iteration,gamma10d,gam123,g_alpha,g_R,RCUT,acceleration,&
        &nbmn,flux,info_CG,ngauge,purebosonic)
        use compiletimeconstants
        use dirac_operator
        use fourier_transform
        use cgm_solver
        use outputstreamnumbers
        use hmc_force
        use Adjust_margins
        use timer
        implicit none
        !***** input *****
        integer, intent(in) :: nbc,nbmn,ngauge,purebosonic
        integer, intent(in) :: max_iteration
        double precision, intent(in) :: max_err
        double precision, intent(in) :: g_alpha
        integer, intent(in) :: ntau,nratio
        double precision :: dtau_xmat,dtau_alpha
        double precision :: dtau_xmat_bos,dtau_alpha_bos
        !double precision dtau_xmat_pf,dtau_alpha_pf
        double precision, intent(in) :: g_R,RCUT
        double precision, intent(in) :: temperature,flux
        double precision, intent(in) :: acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)
        double complex, intent(in) :: pf(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
        double precision, intent(in) :: acceleration(1:nsite)
        double complex, intent(in) :: GAMMA10d(1:ndim,1:nspin,1:nspin)
        !***** input & output *****
        double precision alpha(1:nmat)
        double precision P_alpha(1:nmat)
        double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double complex xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
        double complex P_xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
        !***** output *****
        integer info_CG,iteration
        double complex P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
        !*************************
        double complex delh_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
        !double complex delh_xmat_pf(1:nmat,1:nmat,1:ndim,1:nsite)
        double complex delh_xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision delh_alpha(1:nmat)
        !double precision delh_alpha_pf(1:nmat)
        integer imat,jmat,imom,idim,step

        double complex :: phase(1:nmat,1:nmat,1:2)
        double complex :: Gam123(1:nspin,1:nspin),tmp
        !$acc declare device_resident(phase,Gam123)

        double complex Chi(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:nremez_md,1:npf)

        !$acc declare device_resident(xmat_mom,P_xmat,P_xmat_mom,P_alpha)
        !$acc declare present(xmat)

        !$acc declare present(alpha,g_R,RCUT,temperature,flux,acoeff_md,bcoeff_md,GAMMA10d,g_alpha)
        !$acc declare present(acceleration)

        !$acc declare device_resident(Chi)
        !$acc declare device_resident(delh_xmat,delh_xmat_mom,delh_alpha)
        !!$acc declare device_resident(delh_xmat,delh_xmat_mom,delh_alpha,delh_xmat_pf,delh_alpha_pf)

        call print_time_step("MD evolution start")
        !dtau_xmat_pf=dtau_xmat
        !dtau_alpha_pf=dtau_alpha
        dtau_xmat_bos=dtau_xmat/dble(2*nratio+1)
        dtau_alpha_bos=dtau_alpha/dble(2*nratio+1)
        !call  setup_data_device(alpha,flux,GAMMA10d,phase,Gam123)
        !*******************************
        !*** first step of leap frog ***
        !*******************************
        !Fourier transform from P_xmat to P_xmat_mom
        call FT_P_xmat_device(P_xmat,P_xmat_mom)
        !Fourier transform from xmat to xmat_mom
        call FT_xmat_device(xmat,xmat_mom)
        !move xmat_mom and alpha (1/2)-step forward.
        call Field_step(xmat_mom,P_xmat_mom,alpha,P_alpha,0.5d0*dtau_xmat_bos,0.5d0*dtau_alpha_bos,acceleration)
        !Fourier transform from xmat_mom to xmat
        call FTinv_xmat_device(xmat,xmat_mom)
  
        !************************************
        !*** second step,...,Ntau-th step ***
        !************************************
        info_CG=0
        !info_CG=0 -> CG solver correctly worked.
        !info_CG=1 -> CG solver returned error.
        step=1
        do while ((info_CG.EQ.0).AND.(step.LT.ntau*(2*nratio+1)))
     
            !adjust margin.
            call Adjust_margin_xmat_device(xmat)
            !Calculate Chi_k = (D+bcoeff_md(k))^{-1}*pf by using multi-mass CG-solver
     
            !calculate the force term in the coordinate space.
            !delh_xmat=dH/dX, delh_alpha=dH/(d alpha)
            call print_time_step("MD bos force calculation start")
            call Calc_Force_bosonic_device(delh_xmat,delh_alpha,xmat,alpha,&
                &g_alpha,g_R,RCUT,nbmn,flux,temperature,ngauge)
            call print_time_step("MD bos force calculation end")
            if(rhmc_force.EQ.1) then
                !$acc kernels
                tmp=Sum(abs(delh_xmat))
                !$acc end kernels
                print*,"bosonic force ",step," force ",real(tmp)
            end if
            if((mod(step,2*nratio+1).EQ.nratio+1).OR.(nratio.EQ.0))then
                call update_data_device(alpha,phase)
                call print_time_step("MD cgm solver start")
                if(purebosonic.eq.0) then
                    call cgm_solver_device(nremez_md,bcoeff_md,nbmn,nbc,temperature,&
                        max_err,max_iteration,xmat,phase,Gam123,pf,chi,info_CG,iteration)
                end if
                call print_time_step("MD cgm solver end")
                !Take CG_log
                write(unit_CG_log,*)"molecular evolution",iteration
                call print_time_step("MD ferm force calculation start")
                if(purebosonic.eq.0) then
                    call Add_Force_fermionic_device(dble(2*nratio+1),delh_xmat,delh_alpha,xmat,chi,&
                        g_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md,phase,Gam123,nbc,ngauge)
                end if
                if(rhmc_force.EQ.1) then
                    !$acc kernels
                    tmp=Sum(abs(delh_xmat))
                    !$acc end kernels
                    print*,"bosonic+fermionic force ",step," force ",real(tmp)
                end if
                call print_time_step("MD ferm force calculation end")
                !!$acc kernels
                !delh_xmat=delh_xmat+delh_xmat_pf*dcmplx(2*nratio+1)
                !delh_alpha=delh_alpha+delh_alpha_pf*dble(2*nratio+1)
                !!$acc end kernels
                !if(rhmc_verbose.EQ.1) then
                 !   !$acc kernels
                 !   tmp=Sum(abs(delh_xmat_pf))
                  !  !$acc end kernels
                  !  print*,"md d s ",step," fermion force ",real(tmp)
                !end if
            end if
            if(rhmc_verbose.EQ.1) then
                   !$acc kernels
                tmp=Sum(abs(delh_xmat))
                !$acc end kernels
                print*,"md d s ",step," force ",real(tmp)
            end if
            !Fourier transform from delh_xmat to delh_xmat_mom,
            !i.e. convert the force(delh_xmat) to the Fourier mode (delh_xmat_mom).
            call FT_P_xmat_device(delh_xmat,delh_xmat_mom)
            !move P_xmat_mom and P_alpha one step forward.
            call Momentum_step(P_xmat_mom,delh_xmat_mom,P_alpha,delh_alpha,dtau_xmat_bos,dtau_alpha_bos,acceleration)
            !move xmat_mom and alpha one step forward.
            call Field_step(xmat_mom,P_xmat_mom,alpha,P_alpha,dtau_xmat_bos,dtau_alpha_bos,acceleration)
            !Fourier transform from xmat_mom to xmat
            call FTinv_xmat_device(xmat,xmat_mom)

            if(rhmc_verbose.EQ.1) then
                !$acc kernels
                tmp=Sum(xmat)
                !$acc end kernels
                print*,"md d s ",step," xmat ",tmp
                !$acc kernels
                tmp=Sum(P_xmat)
                !$acc end kernels
                print*,"md d s ",step," mom ",tmp
            end if

            step=step+1
        end do
  
        !*****************
        !*** last step ***
        !*****************
        !info=0 -> so far so good, let's go to the last step.
        !info=1 -> error happened somewhere,
        !      just skip this last step and return info_CG=1.
        if(info_CG.EQ.0)then
            !adjust margin.
            call Adjust_margin_xmat_device(xmat)
            !Calculate Chi_k = (D+bcoeff_md(k))^{-1}*pf by using multi-mass CG-solver
     
            !calculate the force term in the coordinate space.
            !delh_xmat=dH/dX, delh_alpha=dH/(d alpha)
            call print_time_step("MD bos force calculation start")
            call Calc_Force_bosonic_device(delh_xmat,delh_alpha,xmat,alpha,&
                &g_alpha,g_R,RCUT,nbmn,flux,temperature,ngauge)
            call print_time_step("MD bos force calculation end")
            if(nratio.EQ.0)then
                call update_data_device(alpha,phase)
                call print_time_step("MD cgm solver start")
                if(purebosonic.eq.0) then
                    call cgm_solver_device(nremez_md,bcoeff_md,nbmn,nbc,temperature,&
                        max_err,max_iteration,xmat,phase,Gam123,pf,chi,info_CG,iteration)
                end if
                call print_time_step("MD cgm solver end")
                !Take CG_log
                write(unit_CG_log,*)"molecular evolution",iteration
                call print_time_step("MD ferm force calculation start")
                if(purebosonic.eq.0) then
                    call Add_Force_fermionic_device(dble(2*nratio+1),delh_xmat,delh_alpha,xmat,chi,&
                        g_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md,phase,Gam123,nbc,ngauge)
                end if
                call print_time_step("MD ferm force calculation end")
                !delh_xmat=delh_xmat+delh_xmat_pf*dcmplx(2*nratio+1)
                !delh_alpha=delh_alpha+delh_alpha_pf*dble(2*nratio+1)
                !if(rhmc_verbose.EQ.1) then
                !    !$acc kernels
                !    tmp=Sum(abs(delh_xmat_pf))
                !    !$acc end kernels
                !    print*,"md d s  ",step," fermion force ",real(tmp)
                !end if
            end if
            if(rhmc_verbose.EQ.1) then
                !$acc kernels
                tmp=Sum(abs(delh_xmat))
                !$acc end kernels
                print*,"md d s  ",step," complete force ",real(tmp)
            end if

            !Fourier transform from delh_xmat to delh_xmat_mom,
            !i.e. convert the force(delh_xmat) to the Fourier mode (delh_xmat_mom).
            call FT_P_xmat_device(delh_xmat,delh_xmat_mom)
            !move P_xmat_mom and P_alpha one step forward.
            call Momentum_step(P_xmat_mom,delh_xmat_mom,P_alpha,delh_alpha,dtau_xmat_bos,dtau_alpha_bos,acceleration)
            !move xmat_mom and alpha (1/2)-step forward.
            call Field_step(xmat_mom,P_xmat_mom,alpha,P_alpha,0.5d0*dtau_xmat_bos,0.5d0*dtau_alpha_bos,acceleration)
            !Fourier transform from xmat_mom to xmat
            call FTinv_xmat_device(xmat,xmat_mom)
            !adjust margin.
            call Adjust_margin_xmat_device(xmat)
            !Fourier transform from P_xmat_mom to P_xmat
            call FTinv_P_xmat_device(P_xmat,P_xmat_mom)

            if(rhmc_verbose.EQ.1) then
                !$acc kernels
                tmp=Sum(xmat)
                !$acc end kernels
                print*,"md d s ",step," xmat ",tmp
                !$acc kernels
                tmp=Sum(P_xmat)
                !$acc end kernels
                print*,"md d s ",step," mom ",tmp
            end if
        end if
     
        call print_time_step("MD evolution end")
        return

    END subroutine Molecular_Dynamics_device

    subroutine Molecular_Dynamics_device_SW(nbc,temperature,&
        &ntau,nratio,dtau_xmat,dtau_alpha,xmat,alpha,phase,P_xmat,P_alpha,&
        &acoeff_md,bcoeff_md,pf,max_iteration,max_err,&
        &iteration,gamma10d,gam123,g_alpha,g_R,RCUT,acceleration,&
        &nbmn,flux,info_CG,ngauge,purebosonic)
        use compiletimeconstants
        use dirac_operator
        use fourier_transform
        use cgm_solver
        use outputstreamnumbers
        use hmc_force
        use Adjust_margins
        use timer
        implicit none
        !***** input *****
        integer, intent(in) :: nbc,nbmn,ngauge,purebosonic
        integer, intent(in) :: max_iteration
        double precision, intent(in) :: max_err
        double precision, intent(in) :: g_alpha
        integer, intent(in) :: ntau,nratio
        double precision :: dtau_xmat,dtau_alpha
        double precision :: dtau_xmat_bos,dtau_alpha_bos
        !double precision dtau_xmat_pf,dtau_alpha_pf
        double precision, intent(in) :: g_R,RCUT
        double precision, intent(in) :: temperature,flux
        double precision, intent(in) :: acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)
        double complex, intent(in) :: pf(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
        double precision, intent(in) :: acceleration(1:nsite)
        double complex, intent(in) :: GAMMA10d(1:ndim,1:nspin,1:nspin)
        !***** input & output *****
        double precision alpha(1:nmat)
        double precision P_alpha(1:nmat)
        double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double complex xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
        double complex P_xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
        !***** output *****
        integer info_CG,iteration
        double complex P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
        !*************************
        double complex delh_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
        double complex delh_xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision delh_alpha(1:nmat)
        integer imat,jmat,imom,idim,step

        double complex :: phase(1:nmat,1:nmat,1:2)
        double complex :: Gam123(1:nspin,1:nspin),tmp
        !$acc declare device_resident(phase,Gam123)

        double complex Chi(1:nmat,1:nmat,&
            1:nspin,-(nmargin-1):nsite+nmargin,1:nremez_md,1:npf)

        !$acc declare device_resident(xmat_mom,P_xmat,P_xmat_mom,P_alpha)
        !$acc declare present(xmat)

        !$acc declare present(alpha,g_R,RCUT,temperature,flux,acoeff_md,bcoeff_md,GAMMA10d,g_alpha)
        !$acc declare present(acceleration)

        !$acc declare device_resident(Chi)
        !$acc declare device_resident(delh_xmat,delh_xmat_mom,delh_alpha)

        call print_time_step("MDSW outer evolution start")

        dtau_xmat_bos=dtau_xmat/dble(ntau)
        dtau_alpha_bos=dtau_alpha/dble(ntau)
        call FT_P_xmat_device(P_xmat,P_xmat_mom)

        step=0
        info_CG=0

        call Adjust_margin_xmat_device(xmat)
        call update_data_device(alpha,phase)
        call print_time_step("MD cgm solver start")
        if(purebosonic.eq.0) then
            call cgm_solver_device(nremez_md,bcoeff_md,nbmn,nbc,temperature,&
                max_err,max_iteration,xmat,phase,Gam123,pf,chi,info_CG,iteration)
        end if
        call print_time_step("MD cgm solver end")
        !Take CG_log
        write(unit_CG_log,*)"molecular evolution",iteration
        call print_time_step("MD ferm force calculation start")
        !$acc kernels
        delh_xmat=0d0
        delh_alpha=0d0
        !$acc end kernels
        if(purebosonic.eq.0) then
            call Add_Force_fermionic_device(1.0d0,delh_xmat,delh_alpha,xmat,chi,&
                g_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md,phase,Gam123,nbc,ngauge)
        end if
        if(rhmc_force.EQ.1) then
            !$acc kernels
            tmp=Sum(abs(delh_xmat))
            !$acc end kernels
            print*,"fermionic force ",step," force ",real(tmp)
        end if
        call FT_P_xmat_device(delh_xmat,delh_xmat_mom)
        call print_time_step("MD ferm force calculation end")

        call Momentum_step(P_xmat_mom,delh_xmat_mom,P_alpha,delh_alpha,dtau_xmat_bos/6.0d0,dtau_alpha_bos/6.0d0,acceleration)

        if(rhmc_force.EQ.1) then
            !$acc kernels
            tmp=Sum(abs(delh_xmat))
            !$acc end kernels
            print*,"fermionic force ",step," force ",real(tmp)
        end if

        call FT_xmat_device(xmat,xmat_mom)
        do while ((info_CG.EQ.0).AND.(step.LT.ntau))
            step=step+1
            call Molecular_Dynamics_device_SW_inner(nbc,temperature,&
                &nratio,0.5d0*dtau_xmat_bos,0.5d0*dtau_alpha_bos,xmat,xmat_mom,alpha,P_xmat_mom,P_alpha,&
                &g_alpha,g_R,RCUT,acceleration,&
                &nbmn,flux,ngauge,purebosonic)
            !call FTinv_xmat_device(xmat,xmat_mom)

            call Adjust_margin_xmat_device(xmat)
            call update_data_device(alpha,phase)
            call print_time_step("MD cgm solver start")
            if(purebosonic.eq.0) then
                call cgm_solver_device(nremez_md,bcoeff_md,nbmn,nbc,temperature,&
                    max_err,max_iteration,xmat,phase,Gam123,pf,chi,info_CG,iteration)
            end if
            call print_time_step("MD cgm solver end")
            !Take CG_log
            write(unit_CG_log,*)"molecular evolution",iteration
            call print_time_step("MD ferm force calculation start")
            !$acc kernels
            delh_xmat=0d0
            delh_alpha=0d0
            !$acc end kernels
            if(purebosonic.eq.0) then
                call Add_Force_fermionic_device(1.0d0,delh_xmat,delh_alpha,xmat,chi,&
                    g_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md,phase,Gam123,nbc,ngauge)
            end if
            if(rhmc_force.EQ.1) then
                !$acc kernels
                tmp=Sum(abs(xmat))
                !$acc end kernels
                print*,"fermionic force xmat ",step," force ",real(tmp)
                !$acc kernels
                tmp=Sum(abs(delh_xmat))
                !$acc end kernels
                print*,"fermionic force ",step," force ",real(tmp)
            end if
            call FT_P_xmat_device(delh_xmat,delh_xmat_mom)
            call print_time_step("MD ferm force calculation end")

            call Momentum_step(P_xmat_mom,delh_xmat_mom,P_alpha,delh_alpha,2.0d0*dtau_xmat_bos/3.0d0,2.0d0*dtau_alpha_bos/3.0d0,acceleration)
            call Molecular_Dynamics_device_SW_inner(nbc,temperature,&
                &nratio,0.5d0*dtau_xmat_bos,0.5d0*dtau_alpha_bos,xmat,xmat_mom,alpha,P_xmat_mom,P_alpha,&
                &g_alpha,g_R,RCUT,acceleration,&
                &nbmn,flux,ngauge,purebosonic)
            !call FTinv_xmat_device(xmat,xmat_mom)

            call Adjust_margin_xmat_device(xmat)
            call update_data_device(alpha,phase)
            call print_time_step("MD cgm solver start")
            if(purebosonic.eq.0) then
                call cgm_solver_device(nremez_md,bcoeff_md,nbmn,nbc,temperature,&
                    max_err,max_iteration,xmat,phase,Gam123,pf,chi,info_CG,iteration)
            end if
            call print_time_step("MD cgm solver end")
            !Take CG_log
            write(unit_CG_log,*)"molecular evolution",iteration
            call print_time_step("MD ferm force calculation start")
            !$acc kernels
            delh_xmat=0d0
            delh_alpha=0d0
            !$acc end kernels
            if(purebosonic.eq.0) then
                call Add_Force_fermionic_device(1.0d0,delh_xmat,delh_alpha,xmat,chi,&
                    g_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md,phase,Gam123,nbc,ngauge)
            end if
            if(rhmc_force.EQ.1) then
                !$acc kernels
                tmp=Sum(abs(delh_xmat))
                !$acc end kernels
                print*,"fermionic force ",step," force ",real(tmp)
            end if
            call FT_P_xmat_device(delh_xmat,delh_xmat_mom)
            call print_time_step("MD ferm force calculation end")

            if(step.EQ.ntau) then
                call Momentum_step(P_xmat_mom,delh_xmat_mom,P_alpha,delh_alpha,dtau_xmat_bos/6.0d0,dtau_alpha_bos/6.0d0,acceleration)
            else
                call Momentum_step(P_xmat_mom,delh_xmat_mom,P_alpha,delh_alpha,dtau_xmat_bos/3.0d0,dtau_alpha_bos/3.0d0,acceleration)
            end if
        end do
        call FTinv_xmat_device(xmat,xmat_mom)
        call Adjust_margin_xmat_device(xmat)
        call FTinv_P_xmat_device(P_xmat,P_xmat_mom)
        call print_time_step("MDSW Outer evolution end")
        return

    END subroutine Molecular_Dynamics_device_SW

    subroutine Molecular_Dynamics_device_SW_inner(nbc,temperature,&
        &ntau,dtau_xmat,dtau_alpha,xmat,xmat_mom,alpha,P_xmat_mom,P_alpha,&
        &g_alpha,g_R,RCUT,acceleration,&
        &nbmn,flux,ngauge,purebosonic)
        use compiletimeconstants
        use dirac_operator
        use fourier_transform
        use cgm_solver
        use outputstreamnumbers
        use hmc_force
        use Adjust_margins
        use timer
        implicit none
        !***** input *****
        integer, intent(in) :: nbc,nbmn,ngauge,purebosonic
        double precision, intent(in) :: g_alpha
        integer, intent(in) :: ntau
        double precision :: dtau_xmat,dtau_alpha
        double precision :: dtau_xmat_bos,dtau_alpha_bos
        double precision, intent(in) :: g_R,RCUT
        double precision, intent(in) :: temperature,flux
        double precision, intent(in) :: acceleration(1:nsite)
        !***** input & output *****
        double precision alpha(1:nmat)
        double precision P_alpha(1:nmat)
        double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double complex xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
        double complex P_xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
        !*************************
        double complex delh_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
        double complex delh_xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision delh_alpha(1:nmat)
        integer imat,jmat,imom,idim,step
        double complex :: tmp
        !$acc declare device_resident(xmat_mom,P_xmat_mom,P_alpha)
        !$acc declare present(xmat)

        !$acc declare present(alpha,g_R,RCUT,temperature,flux,g_alpha)
        !$acc declare present(acceleration)

        !$acc declare device_resident(delh_xmat,delh_xmat_mom,delh_alpha)

        call print_time_step("MDSW Inner evolution start")
        dtau_xmat_bos=dtau_xmat/dble(ntau)
        dtau_alpha_bos=dtau_alpha/dble(ntau)

        step=0

        call print_time_step("MD bos force calculation start")
        call Adjust_margin_xmat_device(xmat)
        call Calc_Force_bosonic_device(delh_xmat,delh_alpha,xmat,alpha,&
            &g_alpha,g_R,RCUT,nbmn,flux,temperature,ngauge)
        call FT_P_xmat_device(delh_xmat,delh_xmat_mom)
        if(rhmc_force.EQ.1) then
            !$acc kernels
            tmp=Sum(abs(delh_xmat))
            !$acc end kernels
            print*,"bosonic force ",step," force ",real(tmp)
        end if
        call print_time_step("MD bos force calculation end")

        call Momentum_step(P_xmat_mom,delh_xmat_mom,P_alpha,delh_alpha,dtau_xmat_bos/6.0d0,dtau_alpha_bos/6.0d0,acceleration)

        call FT_xmat_device(xmat,xmat_mom)
        do while (step.LT.ntau)
            step=step+1
            call Field_step(xmat_mom,P_xmat_mom,alpha,P_alpha,0.5d0*dtau_xmat_bos,0.5d0*dtau_alpha_bos,acceleration)
            call FTinv_xmat_device(xmat,xmat_mom)

            call print_time_step("MD bos force calculation start")
            call Adjust_margin_xmat_device(xmat)
            call Calc_Force_bosonic_device(delh_xmat,delh_alpha,xmat,alpha,&
                &g_alpha,g_R,RCUT,nbmn,flux,temperature,ngauge)
            call FT_P_xmat_device(delh_xmat,delh_xmat_mom)
            if(rhmc_force.EQ.1) then
                !$acc kernels
                tmp=Sum(abs(delh_xmat))
                !$acc end kernels
                print*,"bosonic force ",step," force ",real(tmp)
                !$acc kernels
                tmp=Sum(abs(P_xmat_mom))
                !$acc end kernels
                print*,"fermionic force xmat mom  ",step," force ",real(tmp)
                !$acc kernels
                tmp=Sum(abs(xmat))
                !$acc end kernels
                print*,"fermionic force xmat ",step," force ",real(tmp)
            end if
            call print_time_step("MD bos force calculation end")

            call Momentum_step(P_xmat_mom,delh_xmat_mom,P_alpha,delh_alpha,2.0d0*dtau_xmat_bos/3.0d0,2.0d0*dtau_alpha_bos/3.0d0,acceleration)
            call Field_step(xmat_mom,P_xmat_mom,alpha,P_alpha,0.5d0*dtau_xmat_bos,0.5d0*dtau_alpha_bos,acceleration)
            call FTinv_xmat_device(xmat,xmat_mom)

            call print_time_step("MD bos force calculation start")
            call Adjust_margin_xmat_device(xmat)
            call Calc_Force_bosonic_device(delh_xmat,delh_alpha,xmat,alpha,&
                &g_alpha,g_R,RCUT,nbmn,flux,temperature,ngauge)
            call FT_P_xmat_device(delh_xmat,delh_xmat_mom)
            if(rhmc_force.EQ.1) then
                !$acc kernels
                tmp=Sum(abs(delh_xmat))
                !$acc end kernels
                print*,"bosonic force ",step," force ",real(tmp)
            end if
            call print_time_step("MD bos force calculation end")

            if(step.EQ.ntau) then
                call Momentum_step(P_xmat_mom,delh_xmat_mom,P_alpha,delh_alpha,dtau_xmat_bos/6.0d0,dtau_alpha_bos/6.0d0,acceleration)
            else
                call Momentum_step(P_xmat_mom,delh_xmat_mom,P_alpha,delh_alpha,dtau_xmat_bos/3.0d0,dtau_alpha_bos/3.0d0,acceleration)
            end if
        end do
        call FTinv_xmat_device(xmat,xmat_mom)
        ! This might not be needed.
        call Adjust_margin_xmat_device(xmat)
        call print_time_step("MDSW Inner evolution end")
        return
    END subroutine Molecular_Dynamics_device_SW_inner

END MODULE HMC_Molecular_Dynamics

