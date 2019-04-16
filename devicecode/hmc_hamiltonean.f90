
module hmc_hamiltonean
    implicit none
contains
    subroutine trace_Xmat2(xmat,trx2)
        use compiletimeconstants
        implicit none

        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double precision, intent(out) :: trx2
        !$acc declare present(xmat)
        integer isite,idim
        integer imat,jmat
        trx2=0d0
        !$acc parallel loop &
        !$acc collapse(4) &
        !$acc reduction(+:trx2)
        do isite=1,nsite
            do idim=1,ndim
                do jmat=1,nmat
                    do imat=1,nmat
                        trx2=trx2&
                            +dble(xmat(imat,jmat,idim,isite)&
                            *dconjg(xmat(imat,jmat,idim,isite)))
                    end do
                end do
            end do
        end do
        !$acc end parallel
        trx2=trx2/dble(nmat*nsite)
        return
    end subroutine trace_Xmat2

    subroutine Calc_Ham_device(temperature,&
        xmat,alpha,P_xmat,P_alpha,ham,pf,chi,&
        acoeff_md,g_R,RCUT,nbmn,flux,phase,ngauge,purebosonic)
        use compiletimeconstants
        use lattice_action
        implicit none

        !***** input *****
        integer, intent(in) :: nbmn,ngauge,purebosonic
        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double precision, intent(in) :: alpha(1:nmat)
        double complex, intent(in) :: P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision, intent(in) :: P_alpha(1:nmat)
        double complex, intent(in) :: pf(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
        double complex, intent(in) :: Chi(1:nmat,1:nmat,1:nspin,&
            &-(nmargin-1):nsite+nmargin,1:nremez_md,1:npf)
        double complex, intent(in) :: phase(1:nmat,1:nmat,1:2)
        double precision, intent(in) :: temperature,flux
        double precision, intent(in) :: g_R,RCUT
        double precision, intent(in) :: acoeff_md(0:nremez_md)
        !$acc declare present(nbmn,xmat,alpha,pf,temperature,flux,g_R,RCUT,acoeff_md,phase)
        !$acc declare device_resident(Chi,P_xmat,P_alpha)
        !***** output *****
        double precision, intent(out) :: ham
        !*******************
        double precision alpha_max,alpha_min
        double precision pi
        double precision trx2
        integer isite
        integer ipf
        integer idim
        integer imat,jmat,kmat
        integer iremez,ispin
        double complex x23,x32,trx123,trx132
        double precision trx2_123,trx2_456789,lattice_spacing
        pi=2d0*dasin(1d0)
        lattice_spacing=1d0/temperature/dble(nsite)
        !**************************************
        !*** The bosonic part of the action ***
        !**************************************
        ham=0.0d0
        call Calc_action_device(temperature,xmat,alpha,ham,phase,ngauge)
        !************************
        !*** 0.5*(momentum)^2 ***
        !************************
        if(rhmc_verbose.EQ.1) then
            print*, "hamilton action dev ",ham
        end if
        !$acc kernels
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
        !$acc end kernels
        if(ngauge.eq.0)then
            !$acc kernels
            do imat=1,nmat
                ham=ham+dble(P_alpha(imat)*P_alpha(imat))*0.5d0
            end do
        !$acc end kernels
        end if
        if(rhmc_verbose.EQ.1) then
            print*, "hamilton momentum dev ",ham
        end if
                !***************************
                !*** pseudo-fermion part ***
                !***************************
        if(purebosonic.eq.0) then
            do ipf=1,npf
                do iremez=1,nremez_md
                    !$acc kernels
                    do isite=1,nsite
                        do ispin=1,nspin
                            do imat=1,nmat
                                do jmat=1,nmat
                                    ham=ham+acoeff_md(iremez)*&
                                        dble(Chi(imat,jmat,ispin,isite,iremez,ipf)*&
                                        dconjg(pf(imat,jmat,ispin,isite,ipf)))
                                end do
                            end do
                        end do
                    end do
                     !$acc end kernels
                end do
            end do
        end if
        if(rhmc_verbose.EQ.1) then
            print*, "hamilton fermion dev ",ham
        end if
        !****************************
        !*** constraint for alpha ***
        !****************************
        if(ngauge.eq.0)then
            !$acc kernels
            alpha_max=Maxval(alpha)
            alpha_min=Minval(alpha)
            !$acc end kernels
            if(alpha_max-alpha_min.LT.2d0*pi)then
                ham=ham-dlog(2d0*pi-(alpha_max-alpha_min))
            end if
        end if
        if(rhmc_verbose.EQ.1) then
            print*, "hamilton constraint alpha dev ",ham
        end if

        !****************************
        !*** constraint for TrX^2 ***
        !****************************
        call trace_Xmat2(xmat,trx2)
        if(trx2.GE.RCUT)then
            ham=ham+g_R*(trx2-RCUT)*dble(NMAT)
        end if
        if(rhmc_verbose.EQ.1) then
            print*, "hamilton constraint R dev ",ham
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
            !$acc kernels
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
            !$acc end kernels
            ham=ham+flux*flux*(0.5d0*trx2_123+0.125d0*trx2_456789)&
                *lattice_spacing*dble(nmat)
            if(rhmc_verbose.EQ.1) then
                print*, "hamilton mass term dev ",ham
            end if
            !******************
            !*** cubic term ***
            !******************
            trx123=(0d0,0d0)
            trx132=(0d0,0d0)
            !$acc kernels
            do isite=1,nsite
                do imat=1,nmat
                    do jmat=1,nmat
                        !(jmat,imat)
                        x23=(0d0,0d0)
                        x32=(0d0,0d0)
                        do kmat=1,nmat
                            x23=x23&
                                +xmat(jmat,kmat,2,isite)*xmat(kmat,imat,3,isite)
                            x32=x32&
                                +xmat(jmat,kmat,3,isite)*xmat(kmat,imat,2,isite)
                        end do
                        trx123=trx123+xmat(imat,jmat,1,isite)*x23
                        trx132=trx132+xmat(imat,jmat,1,isite)*x32
                    end do
                end do
            end do
            !$acc end kernels
            ham=ham+dble((0d0,3d0)*(trx123-trx132))*flux*lattice_spacing*dble(nmat)
            if(rhmc_verbose.EQ.1) then
                print*, "hamilton cubic term dev ",ham
            end if
        end if

        return

    end subroutine Calc_Ham_device

end module hmc_hamiltonean
