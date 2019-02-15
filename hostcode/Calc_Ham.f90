subroutine Calc_Ham(temperature,&
    xmat,alpha,P_xmat,P_alpha,ham,pf,chi,&
    acoeff_md,g_R,RCUT,nbmn,flux,ngauge,purebosonic)

    implicit none

  include '../staticparameters.f90'
    !***** input *****
    integer nbmn,ngauge,purebosonic
    double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
    double precision alpha(1:nmat)
    double complex P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
    double precision P_alpha(1:nmat)
    double complex pf(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
    double complex Chi(1:nremez_md,1:nmat,1:nmat,1:nspin,&
        &-(nmargin-1):nsite+nmargin,1:npf)
    double precision temperature,flux
    double precision g_R,RCUT
    double precision acoeff_md(0:nremez_md)
    !***** output *****
    double precision ham
    !*******************
    double precision alpha_max,alpha_min
    double precision pi
    double precision trx2
    integer isite
    integer ipf
    integer idim
    integer imat,jmat,kmat
    integer iremez,ispin
    double complex x23(1:nmat,1:nmat),x32(1:nmat,1:nmat),trx123,trx132
    double precision trx2_123,trx2_456789,lattice_spacing

    pi=2d0*dasin(1d0)
    lattice_spacing=1d0/temperature/dble(nsite)
    !**************************************
    !*** The bosonic part of the action ***
    !**************************************
    call Calc_action(temperature,xmat,alpha,ham,ngauge)
    if(rhmc_verbose.EQ.1) then
        print*, "hamilton action host ",ham
    end if
    !************************
    !*** 0.5*(momentum)^2 ***
    !************************
    do isite=1,nsite
        do idim=1,ndim
            do imat=1,nmat
                do jmat=1,nmat
                    ham=ham&
                        +dble(P_xmat(imat,jmat,idim,isite)&
                        *P_xmat(jmat,imat,idim,isite))*0.5d0
                end do
            end do
        end do
    end do
    if(ngauge.eq.0)then
        do imat=1,nmat
            ham=ham+dble(P_alpha(imat)*P_alpha(imat))*0.5d0
        end do
    end if
    if(rhmc_verbose.EQ.1) then
        print*, "hamilton momentum host ",ham
    end if
      !***************************
      !*** pseudo-fermion part ***
      !***************************
    !*** This part does not change during the Molecular evolution ***
    !!$  do ipf=1,npf
    !!$     do isite=1,nsite
    !!$        do ispin=1,nspin
    !!$           do imat=1,nmat
    !!$              do jmat=1,nmat
    !!$                 ham=ham+acoeff_md(0)*&
    !!$                      dble(pf(imat,jmat,ispin,isite,ipf)*&
    !!$                      dconjg(pf(imat,jmat,ispin,isite,ipf)))
    !!$              end do
    !!$           end do
    !!$        end do
    !!$     end do
    !!$  end do
    if(purebosonic.eq.0) then
        do ipf=1,npf
            do iremez=1,nremez_md
                do isite=1,nsite
                    do ispin=1,nspin
                        do imat=1,nmat
                            do jmat=1,nmat
                                ham=ham+acoeff_md(iremez)*&
                                    dble(Chi(iremez,imat,jmat,ispin,isite,ipf)*&
                                    dconjg(pf(imat,jmat,ispin,isite,ipf)))
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end if
    if(rhmc_verbose.EQ.1) then
        print*, "hamilton fermion host ",ham
    end if
    !****************************
    !*** constraint for alpha ***
    !****************************
    alpha_max=alpha(1)
    alpha_min=alpha(1)
    if(ngauge.eq.0)then
        do imat=2,nmat
            if(alpha(imat).GT.alpha_max)then
                alpha_max=alpha(imat)
            else if(alpha(imat).LT.alpha_min)then
                alpha_min=alpha(imat)
            end if
        end do
        if(alpha_max-alpha_min.LT.2d0*pi)then
            ham=ham-dlog(2d0*pi-(alpha_max-alpha_min))
        end if
    end if
    if(rhmc_verbose.EQ.1) then
        print*, "hamilton constraint alpha host ",ham
    end if
    !****************************
    !*** constraint for TrX^2 ***
    !****************************
    call Calc_TrX2(xmat,trx2)
    if(trx2.GE.RCUT)then
        ham=ham+g_R*(trx2-RCUT)*dble(NMAT)
    end if
    if(rhmc_verbose.EQ.1) then
        print*, "hamilton constraint R host ",ham
    end if
    !******************************
    !*** Plane wave deformation ***
    !******************************
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
        ham=ham+flux*flux*(0.5d0*trx2_123+0.125d0*trx2_456789)&
            *lattice_spacing*dble(nmat)
        if(rhmc_verbose.EQ.1) then
            print*, "hamilton mass term host ",ham
        end if
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
        ham=ham+dble((0d0,3d0)*(trx123-trx132))*flux*lattice_spacing*dble(nmat)
        if(rhmc_verbose.EQ.1) then
            print*, "hamilton cubic term host ",ham
        end if
    end if

    return

END subroutine Calc_Ham
