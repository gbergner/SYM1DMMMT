! alpha -> alpha + P_alpha*dtau_alpha
! xmat -> xmat + P_xmat*dtau_xmat
! P_alpha -> P_alpha - delh_alpha*dtau_alpha
! P_xmat -> P_xmat - delh_xmat*dtau_xmat
! delh_xmat(imat,jmat)=dH/dxmat(jmat,imat)

subroutine Molecular_Dynamics(nbc,temperature,&
    &ntau,nratio,dtau_xmat,dtau_alpha,xmat,alpha,P_xmat,P_alpha,&
    &acoeff_md,bcoeff_md,pf,max_iteration,max_err,&
    &iteration,gamma10d,g_alpha,g_R,RCUT,acceleration,&
    &nbmn,flux,info_CG)

    implicit none

  include '../staticparameters.f90'
  include '../Fourier.inc'
  include '../unit_number.inc'
    !***** input *****
    integer nbc,nbmn
    integer max_iteration
    double precision max_err
    double precision g_alpha
    integer ntau,nratio
    double precision dtau_xmat,dtau_alpha
    double precision dtau_xmat_bos,dtau_alpha_bos
    !double precision dtau_xmat_pf,dtau_alpha_pf
    double precision g_R,RCUT
    doubleprecision temperature,flux
    double precision acoeff_md(0:nremez_md),bcoeff_md(1:nremez_md)
    double complex pf(1:nmat,1:nmat,&
        1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
    double precision acceleration(1:nsite)
    double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
    !***** input & output *****
    double precision alpha(1:nmat)
    double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    !***** output *****
    integer info_CG,iteration
    double complex P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
    double precision P_alpha(1:nmat)
    !*************************
    double complex xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
    double complex P_xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
    double complex delh_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
    double complex delh_xmat_pf(1:nmat,1:nmat,1:ndim,1:nsite)
    double complex delh_xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
    double precision delh_alpha(1:nmat)
    double precision delh_alpha_pf(1:nmat)
    integer imat,jmat
    integer imom
    integer idim
    integer step
    double complex Chi(1:nremez_md,1:nmat,1:nmat,&
        1:nspin,-(nmargin-1):nsite+nmargin,1:npf)

    !dtau_xmat_pf=dtau_xmat
    !dtau_alpha_pf=dtau_alpha
    dtau_xmat_bos=dtau_xmat/dble(2*nratio+1)
    dtau_alpha_bos=dtau_alpha/dble(2*nratio+1)
  
    !*******************************
    !*** first step of leap frog ***
    !*******************************
    !Fourier transform from P_xmat to P_xmat_mom
    call Fourier_transform_P_xmat(P_xmat,P_xmat_mom,x2p)
    !Fourier transform from xmat to xmat_mom
    call Fourier_transform_xmat(xmat,xmat_mom,x2p)
    !move xmat_mom and alpha (1/2)-step forward.
    do imom=1,nsite
        do idim=1,ndim
            do imat=1,nmat
                do jmat=1,nmat
                    xmat_mom(imat,jmat,idim,imom)=&
                        &xmat_mom(imat,jmat,idim,imom)&
                        &+P_xmat_mom(imat,jmat,idim,imom)&
                        &*dcmplx(0.5d0*dtau_xmat_bos*acceleration(imom))
                end do
            end do
        end do
    end do
    alpha=alpha+P_alpha*0.5d0*dtau_alpha_bos
    !Fourier transform from xmat_mom to xmat
    call Fourier_transform_xmat(xmat,xmat_mom,p2x)
    !************************************
    !*** second step,...,Ntau-th step ***
    !************************************
    info_CG=0
    !info_CG=0 -> CG solver correctly worked.
    !info_CG=1 -> CG solver returned error.
    step=1
    do while ((info_CG.EQ.0).AND.(step.LT.ntau*(2*nratio+1)))
     
        !adjust margin.
        call Adjust_margin_xmat(xmat)
        !Calculate Chi_k = (D+bcoeff_md(k))^{-1}*pf by using multi-mass CG-solver
     
        !calculate the force term in the coordinate space.
        !delh_xmat=dH/dX, delh_alpha=dH/(d alpha)
        call Calc_Force_bosonic(delh_xmat,delh_alpha,xmat,alpha,chi,&
            &GAMMA10d,g_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md)
        if((mod(step,2*nratio+1).EQ.nratio+1).OR.(nratio.EQ.0))then
            call solver_biCGm(nbc,nbmn,nremez_md,&
                &xmat,alpha,pf,chi,GAMMA10d,&
                &bcoeff_md,max_err,max_iteration,iteration,&
                &temperature,flux,info_CG)
            !Take CG_log
            write(unit_CG_log,*)"molecular evolution",iteration
            call Calc_Force_fermionic(delh_xmat_pf,delh_alpha_pf,xmat,alpha,chi,&
                &GAMMA10d,g_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md)
            !end if
            !if(mod(step,2*nratio+1).EQ.0)then
            delh_xmat=delh_xmat+delh_xmat_pf*dcmplx(2*nratio+1)
            delh_alpha=delh_alpha+delh_alpha_pf*dble(2*nratio+1)
            if(rhmc_verbose.EQ.1) then
                print*,"md h s ",step," fermion force ",Sum(abs(delh_xmat_pf))
            end if
        end if

        if(rhmc_verbose.EQ.1) then
            print*,"md h s ",step," force ",Sum(abs(delh_xmat))
        end if


        !Fourier transform from delh_xmat to delh_xmat_mom,
        !i.e. convert the force(delh_xmat) to the Fourier mode (delh_xmat_mom).
        call Fourier_transform_P_xmat(delh_xmat,delh_xmat_mom,1)
        !move P_xmat_mom and P_alpha one step forward.
        P_alpha=P_alpha-delh_alpha*dtau_alpha_bos
        do imom=1,nsite
            do idim=1,ndim
                do jmat=1,nmat
                    do imat=1,nmat
                        P_xmat_mom(imat,jmat,idim,imom)=&
                            &P_xmat_mom(imat,jmat,idim,imom)&
                            &-delh_xmat_mom(imat,jmat,idim,imom)&
                            &*dcmplx(dtau_xmat_bos*acceleration(imom))
                    end do
                end do
            end do
        end do
        !move xmat_mom and alpha one step forward.
        do imom=1,nsite
            do idim=1,ndim
                do jmat=1,nmat
                    do imat=1,nmat
                        xmat_mom(imat,jmat,idim,imom)=&
                            &xmat_mom(imat,jmat,idim,imom)&
                            &+P_xmat_mom(imat,jmat,idim,imom)&
                            &*dcmplx(dtau_xmat_bos*acceleration(imom))
                    end do
                end do
            end do
        end do
        alpha=alpha+P_alpha*dtau_alpha_bos
        !Fourier transform from xmat_mom to xmat
        call Fourier_transform_xmat(xmat,xmat_mom,p2x)

        if(rhmc_verbose.EQ.1) then
            print*,"md h s ",step," xmat ",Sum(xmat)
            print*,"md h s ",step," mom ",Sum(P_xmat)
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
        call Adjust_margin_xmat(xmat)
        !Calculate Chi_k = (D+bcoeff_md(k))^{-1}*pf by using multi-mass CG-solver
     
        !calculate the force term in the coordinate space.
        !delh_xmat=dH/dX, delh_alpha=dH/(d alpha)
        call Calc_Force_bosonic(delh_xmat,delh_alpha,xmat,alpha,chi,&
            &GAMMA10d,g_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md)
        if(nratio.EQ.0)then
            call solver_biCGm(nbc,nbmn,nremez_md,&
                &xmat,alpha,pf,chi,GAMMA10d,&
                &bcoeff_md,max_err,max_iteration,iteration,&
                &temperature,flux,info_CG)
            !Take CG_log
            write(unit_CG_log,*)"molecular evolution",iteration
            call Calc_Force_fermionic(delh_xmat_pf,delh_alpha_pf,xmat,alpha,chi,&
                &GAMMA10d,g_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md)
            !end if
            !if(mod(step,2*nratio+1).EQ.0)then
            delh_xmat=delh_xmat+delh_xmat_pf*dcmplx(2*nratio+1)
            delh_alpha=delh_alpha+delh_alpha_pf*dble(2*nratio+1)
             if(rhmc_verbose.EQ.1) then
                print*,"md h s ",step," fermion force ",Sum(abs(delh_xmat_pf))
            end if
        end if
        if(rhmc_verbose.EQ.1) then
            print*,"md h s ",step," force ",Sum(abs(delh_xmat))
        end if
        !Fourier transform from delh_xmat to delh_xmat_mom,
        !i.e. convert the force(delh_xmat) to the Fourier mode (delh_xmat_mom).
        call Fourier_transform_P_xmat(delh_xmat,delh_xmat_mom,x2p)
        !move P_xmat_mom and P_alpha one step forward.
        P_alpha=P_alpha-delh_alpha*dtau_alpha_bos
        do imom=1,nsite
            do idim=1,ndim
                do jmat=1,nmat
                    do imat=1,nmat
                        P_xmat_mom(imat,jmat,idim,imom)=&
                            &P_xmat_mom(imat,jmat,idim,imom)&
                            &-delh_xmat_mom(imat,jmat,idim,imom)&
                            &*dcmplx(dtau_xmat_bos*acceleration(imom))
                    end do
                end do
            end do
        end do
        !move xmat_mom and alpha (1/2)-step forward.
        do imom=1,nsite
            do idim=1,ndim
                do jmat=1,nmat
                    do imat=1,nmat
                        xmat_mom(imat,jmat,idim,imom)=&
                            &xmat_mom(imat,jmat,idim,imom)&
                            &+P_xmat_mom(imat,jmat,idim,imom)&
                            &*dcmplx(0.5d0*dtau_xmat_bos*acceleration(imom))
                    end do
                end do
            end do
        end do
        alpha=alpha+P_alpha*0.5d0*dtau_alpha_bos
        !Fourier transform from xmat_mom to xmat
        call Fourier_transform_xmat(xmat,xmat_mom,p2x)
        !adjust margin.
        call Adjust_margin_xmat(xmat)
        !Fourier transform from P_xmat_mom to P_xmat
        call Fourier_transform_P_xmat(P_xmat,P_xmat_mom,p2x)

        if(rhmc_verbose.EQ.1) then
            print*,"md h s ",step," xmat ",Sum(xmat)
            print*,"md h s ",step," mom ",Sum(P_xmat)
        end if

    end if
     
  
    return

END subroutine Molecular_Dynamics
