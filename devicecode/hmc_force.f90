!****************************************************************
!***** Calculate the force term in the molecular evolution.            *****
!***** The fermionic force computation had to be completely rewritten. *****
!*****  Georg Bergner                                                  *****
!****************************************************************
! delh_alpha(imat) = dH/(d alpha_imat)
! delh_xmat(imat,jmat)=dH/(d xmat(jmat,imat)), 
!      be careful about the difference of ordering of indices.
! ******* Note that "Force" = (-1)*delh_xmat, (-1)*delh_alpha: ******
! alpha -> alpha + P_alpha*dtau_alpha
! xmat -> xmat + P_xmat*dtau_xmat
! P_alpha -> P_alpha - delh_alpha*dtau_alpha
! P_xmat -> P_xmat - delh_xmat*dtau_xmat

module hmc_force
    implicit none
contains
    SUBROUTINE Calc_Force_bosonic_device(delh_xmat,delh_alpha,xmat,alpha,gcoeff_alpha,g_R,RCUT,nbmn,flux,temperature)

        use compiletimeconstants
        implicit none


        !****** input *****
        integer, intent(in) :: nbmn
        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double precision, intent(in) :: alpha(1:nmat)
        double precision, intent(in) :: gcoeff_alpha
        double precision, intent(in) :: g_R,RCUT
        double precision, intent(in) :: temperature,flux
        !$acc declare present(g_R,RCUT,temperature,flux,gcoeff_alpha)
        !***** output *****
        double complex, intent(out) :: delh_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision, intent(out) :: delh_alpha(1:nmat)
        !$acc declare device_resident(delh_xmat,delh_alpha,xmat,alpha)
        !*******************
        double precision alpha_max,alpha_min,trx2
        double precision lattice_spacing,pi
        integer imat,jmat,kmat,lmat,imax,imin
        integer idim,jdim
        integer isite,isite_p1,isite_m1
        !double complex commutator(1:nmat,1:nmat,1:nsite,1:ndim,1:ndim),uxdu(1:nmat,1:nmat),&
         !    &udxu(1:nmat,1:nmat),xxx,ei,ej
        double complex commutator(1:nmat,1:nmat,1:ndim,1:ndim),uxdu(1:nmat,1:nmat),&
            &udxu(1:nmat,1:nmat),xxx,ei,ej
        !$acc declare device_resident(commutator)
        integer iremez
        integer ispin
        integer ipf
        double complex com12(1:nmat,1:nmat),com23(1:nmat,1:nmat),&
            &com31(1:nmat,1:nmat),temp,uxumx,tmp1,tmp2



        pi=2d0*dasin(1d0)
        lattice_spacing=1d0/temperature/dble(nsite)
        !*****************************
        !*** calculate delh_alpha ****
        !*****************************
        !$acc kernels
        delh_alpha=0d0
        !$acc end kernels
        if(ngauge.EQ.1)then
            ! delh_xmat=(0d0,0d0)
            if(nimprove.EQ.0)then
                !$acc kernels
                do isite=1,nsite
                    isite_p1=isite+1
           
                    do idim=1,ndim
                        do imat=1,nmat
                            do jmat=1,nmat
                                !exp(i*alpha_i)
                                ei=dcmplx(dcos(alpha(imat)/dble(nsite)))&
                                    +(0d0,1d0)*dcmplx(dsin(alpha(imat)/dble(nsite)))
                                !exp(-i*alpha_j)
                                ej=dcmplx(dcos(alpha(jmat)/dble(nsite)))&
                                    -(0d0,1d0)*dcmplx(dsin(alpha(jmat)/dble(nsite)))
                    
                                delh_alpha(imat)=delh_alpha(imat)&
                                    -dble(ei*xmat(imat,jmat,idim,isite_p1)&
                                    *ej*xmat(jmat,imat,idim,isite)*(0d0,1d0))
                                delh_alpha(jmat)=delh_alpha(jmat)&
                                    +dble(ei*xmat(imat,jmat,idim,isite_p1)&
                                    *ej*xmat(jmat,imat,idim,isite)*(0d0,1d0))
                    
                            end do
                        end do
                    end do
                end do
                !$acc end kernels
            else if(nimprove.EQ.1)then
                !$acc kernels
                do imat=1,nmat
                    do jmat=1,nmat
                        tmp1=(0d0,0d0)
                        tmp2=(0d0,0d0)
                        do isite=1,nsite
                            do idim=1,ndim
                                !-0.5*u^2*x(t+2a)*(u^dagger)^2 + 2*u*x(t+a)*u^dagger - 1.5*x(t)

                    
                                !exp(i*alpha_i)
                                ei=dcmplx(dcos(alpha(imat)&
                                    &/dble(nsite)))&
                                    &+(0d0,1d0)*dcmplx(dsin(alpha(imat)&
                                    &/dble(nsite)))
                                !exp(-i*alpha_j)
                                ej=dcmplx(dcos(alpha(jmat)&
                                    &/dble(nsite)))&
                                    -(0d0,1d0)*dcmplx(dsin(alpha(jmat)&
                                    &/dble(nsite)))
                                uxumx=&
                                    -(0.5d0,0d0)*ei*ei*xmat(imat,jmat,idim,isite+2)*ej*ej&
                                    +(2d0,0d0)*ei*xmat(imat,jmat,idim,isite+1)*ej&
                                    -(1.5d0,0d0)*xmat(imat,jmat,idim,isite)
                    
                                tmp1&
                                    &=tmp1&
                                    &+dble(((2d0,0d0)*ei*xmat(imat,jmat,idim,isite+1)*ej&
                                    &-ei*ei*xmat(imat,jmat,idim,isite+2)*ej*ej)&
                                    *dconjg(uxumx)*(0d0,1d0))
                    
                                tmp2&
                                    &=tmp2&
                                    &-dble(((2d0,0d0)*ei*xmat(imat,jmat,idim,isite+1)*ej&
                                    &-ei*ei*xmat(imat,jmat,idim,isite+2)*ej*ej)&
                                    *dconjg(uxumx)*(0d0,1d0))
                    
                            end do
                        end do
                        delh_alpha(imat)=delh_alpha(imat)+tmp1
                        delh_alpha(jmat)=delh_alpha(jmat)+tmp2
                    end do
                end do
               !$acc end kernels
            end if
     
            !$acc kernels
            delh_alpha=delh_alpha*dble(nmat)/lattice_spacing/dble(nsite)
            !$acc end kernels

            !gauge-fixing term
            !$acc kernels
            do imat=1,nmat
                do jmat=1,nmat
                    if(jmat.NE.imat)then
                        delh_alpha(imat)=delh_alpha(imat)&
                            -1d0/dtan(0.5d0*(alpha(imat)-alpha(jmat)))
                    end if
                end do
            end do
           !$acc end kernels
        end if

        !***************************
        !*** calculate delh_xmat ***
        !***************************
        !$acc kernels
        delh_xmat=(0d0,0d0)
        !$acc end kernels
        if(nimprove.EQ.0)then
            !$acc kernels
            do isite=1,nsite
                isite_p1=isite+1
                isite_m1=isite-1
        
                do idim=1,ndim
                    !u(t)*x(t+a)*u(t)^¥dagger
                    !u(t-a)^¥dagger*x(t-a)*u(t-a)
                    do imat=1,nmat
                        do jmat=1,nmat
                            !exp(i*alpha_i)
                            ei=dcmplx(dcos(alpha(imat)/dble(nsite)))&
                                +(0d0,1d0)*dcmplx(dsin(alpha(imat)/dble(nsite)))
                            !exp(-i*alpha_j)
                            ej=dcmplx(dcos(alpha(jmat)/dble(nsite)))&
                                -(0d0,1d0)*dcmplx(dsin(alpha(jmat)/dble(nsite)))
                            uxdu(imat,jmat)=&
                                ei*xmat(imat,jmat,idim,isite_p1)*ej
                            udxu(jmat,imat)=&
                                ej*xmat(jmat,imat,idim,isite_m1)*ei
                        end do
                    end do
           
                    do imat=1,nmat
                        do jmat=1,nmat
                            delh_xmat(imat,jmat,idim,isite)=&
                                dcmplx(2d0)*xmat(imat,jmat,idim,isite)
                 
                            delh_xmat(imat,jmat,idim,isite)=&
                                delh_xmat(imat,jmat,idim,isite)&
                                -uxdu(imat,jmat)&
                                -udxu(imat,jmat)
                        end do
                    end do
                end do
            end do
           !$acc end kernels
        else if(nimprove.EQ.1)then
            !$acc kernels
            do isite=1,nsite
                do idim=1,ndim
                    do imat=1,nmat
                        do jmat=1,nmat
                            !exp(i*alpha_i)
                            ei=dcmplx(dcos(alpha(imat)&
                                &/dble(nsite)))&
                                &+(0d0,1d0)*dcmplx(dsin(alpha(imat)&
                                &/dble(nsite)))
                            !exp(-i*alpha_j)
                            ej=dcmplx(dcos(alpha(jmat)&
                                &/dble(nsite)))&
                                &-(0d0,1d0)*dcmplx(dsin(alpha(jmat)&
                                &/dble(nsite)))
                 
                            delh_xmat(imat,jmat,idim,isite)=&
                                &delh_xmat(imat,jmat,idim,isite)&
                                &+(6.5d0,0d0)*xmat(imat,jmat,idim,isite)&
                                &-(4d0,0d0)*ei*xmat(imat,jmat,idim,isite+1)*ej&
                                &-(4d0,0d0)*dconjg(ei)*xmat(imat,jmat,idim,isite-1)&
                                &*dconjg(ej)&
                                &+(0.75d0,0d0)*ei*ei*xmat(imat,jmat,idim,isite+2)*ej*ej&
                                &+(0.75d0,0d0)*dconjg(ei*ei)&
                                &*xmat(imat,jmat,idim,isite-2)*dconjg(ej*ej)
                 
                        end do
                    end do
                end do
            end do
        !$acc end kernels


        end if

        !$acc kernels
        delh_xmat=delh_xmat*dcmplx(nmat)/dcmplx(lattice_spacing)
        !$acc end kernels

        !commutator term
          !commutator term
           !$acc kernels
        do isite=1,nsite
            commutator=(0d0,0d0)
            do idim=1,ndim-1
                do jdim=idim+1,ndim
                    do imat=1,nmat
                        do jmat=1,nmat
                            do kmat=1,nmat
                                commutator(imat,jmat,idim,jdim)=&
                                    commutator(imat,jmat,idim,jdim)&
                                    +xmat(imat,kmat,idim,isite)*xmat(kmat,jmat,jdim,isite)&
                                    -xmat(imat,kmat,jdim,isite)*xmat(kmat,jmat,idim,isite)
                            end do
                        end do
                    end do
                end do
            end do

            do idim=1,ndim
                do jdim=1,ndim
                    do imat=1,nmat
                        do jmat=1,nmat
                            xxx=(0d0,0d0)
                            if(jdim.GT.idim)then
                                do kmat=1,nmat
                                    xxx=xxx&
                                        +xmat(imat,kmat,jdim,isite)&
                                        *commutator(kmat,jmat,idim,jdim)&
                                        -commutator(imat,kmat,idim,jdim)&
                                        *xmat(kmat,jmat,jdim,isite)
                                end do
                            else if(jdim.LT.idim)then
                                do kmat=1,nmat
                                    xxx=xxx&
                                        -xmat(imat,kmat,jdim,isite)&
                                        *commutator(kmat,jmat,jdim,idim)&
                                        +commutator(imat,kmat,jdim,idim)&
                                        *xmat(kmat,jmat,jdim,isite)
                                end do
                            end if

                            delh_xmat(imat,jmat,idim,isite)=&
                                delh_xmat(imat,jmat,idim,isite)&
                                -dcmplx(nmat)*dcmplx(lattice_spacing)*xxx
                        end do
                    end do

                end do
            end do

        end do
          !$acc end kernels

        !****************************
        !*** constraint for alpha ***
        !****************************
        if(ngauge.EQ.1)then
            imax=1
            imin=1
            alpha_max=0d0
            alpha_min=0d0
            !$acc kernels
            do imat=1,nmat
                if(alpha(imat).GT.alpha_max)then
                    imax=imat
                    alpha_max=alpha(imat)
                else if(alpha(imat).LT.alpha_min)then
                    imin=imat
                    alpha_min=alpha(imat)
                end if
            end do
            !$acc end kernels
            if(alpha_max-alpha_min.LT.2d0*pi)then
                !$acc kernels
                delh_alpha(imax)=delh_alpha(imax)&
                    +1d0/(2d0*pi-(alpha_max-alpha_min)+1d0/gcoeff_alpha)
                delh_alpha(imin)=delh_alpha(imin)&
                    -1d0/(2d0*pi-(alpha_max-alpha_min)+1d0/gcoeff_alpha)
               !$acc end kernels
            else if(alpha_max-alpha_min.GE.2d0*pi)then
                  !$acc kernels
                delh_alpha(imax)=delh_alpha(imax)+gcoeff_alpha
                delh_alpha(imin)=delh_alpha(imin)-gcoeff_alpha
               !$acc end kernels
            end if
        end if
        !****************************
        !*** constraint for TrX^2 ***
        !****************************
        call Calc_TrX2(xmat,trx2)
        if(trx2.GE.RCUT)then
            !$acc kernels
            do isite=1,nsite
                do idim=1,ndim
                    do imat=1,nmat
                        do jmat=1,nmat
                            delh_xmat(imat,jmat,idim,isite)=&
                                delh_xmat(imat,jmat,idim,isite)&
                                +2d0*g_R*xmat(imat,jmat,idim,isite)&
                                /dble(nsite)
                        end do
                    end do
                end do
            end do
         !$acc end kernels
        end if
        !******************************
        !*** Plane wave deformation ***
        !******************************
        if(nbmn.EQ.1)then
            !***************************
            !*** deriv. of mass term ***
            !***************************
            !$acc kernels
            do isite=1,nsite
                do imat=1,nmat
                    do jmat=1,nmat
                        tmp1=dcmplx(flux*flux*lattice_spacing*dble(nmat))
                        do idim=1,3
                            delh_xmat(imat,jmat,idim,isite)=&
                                delh_xmat(imat,jmat,idim,isite)&
                                +xmat(imat,jmat,idim,isite)*tmp1
                        end do
                        tmp1=dcmplx(0.25d0*flux*flux*lattice_spacing*dble(nmat))
                        do idim=4,9
                            delh_xmat(imat,jmat,idim,isite)=&
                                delh_xmat(imat,jmat,idim,isite)&
                                +xmat(imat,jmat,idim,isite)*tmp1
                        end do
                    end do
                end do
            end do
            !$acc end kernels
            !******************
            !*** cubic term ***
            !******************
            !$acc kernels
            do isite=1,nsite
                com12=(0d0,0d0)
                com23=(0d0,0d0)
                com31=(0d0,0d0)
                do imat=1,nmat
                    do jmat=1,nmat
                        do kmat=1,nmat
                            com12(imat,jmat)=com12(imat,jmat)&
                                +xmat(imat,kmat,1,isite)*xmat(kmat,jmat,2,isite)&
                                -xmat(imat,kmat,2,isite)*xmat(kmat,jmat,1,isite)
                            com23(imat,jmat)=com23(imat,jmat)&
                                +xmat(imat,kmat,2,isite)*xmat(kmat,jmat,3,isite)&
                                -xmat(imat,kmat,3,isite)*xmat(kmat,jmat,2,isite)
                            com31(imat,jmat)=com31(imat,jmat)&
                                +xmat(imat,kmat,3,isite)*xmat(kmat,jmat,1,isite)&
                                -xmat(imat,kmat,1,isite)*xmat(kmat,jmat,3,isite)
                        end do
                    end do
                end do
                tmp1=(0d0,3d0)*dcmplx(flux*lattice_spacing*dble(nmat))
                do imat=1,nmat
                    do jmat=1,nmat
                        delh_xmat(imat,jmat,1,isite)=&
                            delh_xmat(imat,jmat,1,isite)+com23(imat,jmat)*tmp1
                        delh_xmat(imat,jmat,2,isite)=&
                            delh_xmat(imat,jmat,2,isite)+com31(imat,jmat)*tmp1
                        delh_xmat(imat,jmat,3,isite)=&
                            delh_xmat(imat,jmat,3,isite)+com12(imat,jmat)*tmp1
                    end do
                end do
            end do
             !$acc end kernels
        end if

        return
  
    END SUBROUTINE Calc_Force_Bosonic_device

        !**********************************************************
        ! This severly reduces the memory requirements!!
    SUBROUTINE Add_Force_fermionic_device(prefact,delh_xmat,delh_alpha,xmat,chi,&
        gcoeff_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md,phase,Gam123,nbc)

        use compiletimeconstants
        use dirac_operator
        use timer
        use gammamatrix
        implicit none


        !****** input *****
        integer nbmn,nbc
        double precision, intent(in) :: prefact
        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double complex, intent(in) :: Chi(1:nmat,1:nmat,1:nspin,&
            &-(nmargin-1):nsite+nmargin,1:nremez_md,1:npf)
        double precision, intent(in) :: gcoeff_alpha
        double precision, intent(in) :: g_R,RCUT
        double precision, intent(in) :: temperature,flux
        double precision, intent(in) :: acoeff_md(0:nremez_md)
        double complex, intent(in)  :: phase(1:nmat,1:nmat,1:2)
        double complex, intent(in)  :: Gam123(1:nspin,1:nspin)
        !$acc declare present(nbmn,nbc,xmat,gcoeff_alpha,g_R,RCUT,temperature,flux,acoeff_md)
        !$acc declare device_resident(Chi,phase,Gam123)
        !***** output *****
        double complex, intent(inout) :: delh_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
        double precision, intent(inout) :: delh_alpha(1:nmat)
        !$acc declare device_resident(delh_xmat,delh_alpha)
        !*******************
        double precision lattice_spacing,pi,tmpsum
        integer imat,jmat,kmat,lmat
        integer idim,jdim
        integer isite
        double complex pf1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        !$acc declare device_resident(pf1)
        integer iremez
        integer ispin
        integer ipf
        integer countr, jspin
        double complex :: tmpfact,gamtmp

        pi=2d0*dasin(1d0)
        lattice_spacing=1d0/temperature/dble(nsite)

        do ipf=1,npf
            do iremez=1,nremez_md
                call set_boundary_device(nbc,chi(:,:,:,:,iremez,ipf))
                call Multiply_Dirac_device(temperature,xmat,phase,Gam123,nbmn,&
                    chi(:,:,:,:,iremez,ipf),pf1)
                          !-acoeff*(M*chi)^dagger*(dM/dX)*chi
                do countr=1,144
                    ispin=gamispin(countr)
                    jspin=gamjspin(countr)
                    idim=gamidim(countr)
                    gamtmp=gamgam(countr)
                    tmpfact=-dcmplx(acoeff_md(iremez)*prefact*lattice_spacing)*gamtmp
                    !$acc kernels
                    do isite=1,nsite
                        do jmat=1,nmat
                            do kmat=1,nmat
                                do imat=1,nmat
                                    !pf2(imat,jmat,ispin,isite)=pf2(imat,jmat,ispin,isite)-dcmplx(lattice_spacing)*gamtmp&
                                    !*(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                                    !-xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
                                    ! note that the transpose has to be used.
                                    delh_xmat(kmat,jmat,idim,isite)=delh_xmat(kmat,jmat,idim,isite)-tmpfact&
                                        *chi(kmat,imat,jspin,isite,iremez,ipf)*dconjg(pf1(jmat,imat,ispin,isite))
                                    delh_xmat(kmat,jmat,idim,isite)=delh_xmat(kmat,jmat,idim,isite)+tmpfact&
                                        *chi(imat,jmat,jspin,isite,iremez,ipf)*dconjg(pf1(imat,kmat,ispin,isite))
                                    !take into account both
                                    !(M*chi)^dagger*(dM/dX)*chi and ((dM/dX)*chi)^dagger*(M*chi)
                                    delh_xmat(kmat,jmat,idim,isite)=delh_xmat(kmat,jmat,idim,isite)-dconjg(tmpfact&
                                        *chi(jmat,imat,jspin,isite,iremez,ipf)*dconjg(pf1(kmat,imat,ispin,isite)))
                                    delh_xmat(kmat,jmat,idim,isite)=delh_xmat(kmat,jmat,idim,isite)+dconjg(tmpfact&
                                        *chi(imat,kmat,jspin,isite,iremez,ipf)*dconjg(pf1(imat,jmat,ispin,isite)))
                                end do
                            end do
                        end do
                    end do
                   !$acc end kernels
                end do

                !********************
                !*** Naive Action ***
                !********************
                !take into account both
                !(M*chi)^dagger*(dM/d¥alpha)*chi and ((dM/d¥alpha)*chi)^dagger*(M*chi) --> factor 2
                if(ngauge.EQ.1)then
                    tmpfact=-dcmplx(prefact*acoeff_md(iremez))*(0d0,-2d0)/dcmplx(nsite)
                    if(nimprove.EQ.0)then
                           !aij=alpha(imat)-alpha(jmat)
                           !aij=aij/dble(nsite)
                           !phase=dcmplx(dcos(aij))+(0d0,1d0)*dcmplx(dsin(aij))=e^(i aij/nsite)
                           ! phase2=phase*phase=e^(i2 aij/nsite)
                        !$acc kernels
                        do isite=1,nsite
                            do ispin=1,8
                                do jmat=1,nmat
                                    do imat=1,nmat
                                        delh_alpha(jmat)=delh_alpha(jmat)&
                                            +dble(tmpfact*phase(imat,jmat,1)*chi(imat,jmat,ispin,isite+1,iremez,ipf)&
                                            *dconjg(pf1(imat,jmat,ispin+8,isite)))&
                                            +dble(tmpfact*dconjg(phase(imat,jmat,1))*chi(imat,jmat,ispin+8,isite-1,iremez,ipf)&
                                            *dconjg(pf1(imat,jmat,ispin,isite)))
                                        delh_alpha(imat)=delh_alpha(imat)&
                                            -dble(tmpfact*phase(imat,jmat,1)*chi(imat,jmat,ispin,isite+1,iremez,ipf)&
                                            *dconjg(pf1(imat,jmat,ispin+8,isite)))&
                                            -dble(tmpfact*dconjg(phase(imat,jmat,1))*chi(imat,jmat,ispin+8,isite-1,iremez,ipf)&
                                            *dconjg(pf1(imat,jmat,ispin,isite)))

                                        !pf2(imat,jmat,ispin+8,isite)=&
                                        !    pf2(imat,jmat,ispin+8,isite)&
                                        !    +(1d0,0d0)*phase(imat,jmat,1)*pf1(imat,jmat,ispin,isite+1)&
                                        !    -(1d0,0d0)*pf1(imat,jmat,ispin,isite)
                                        !pf2(imat,jmat,ispin,isite)=&
                                        !    pf2(imat,jmat,ispin,isite)&
                                        !   -(1d0,0d0)*dconjg(phase(imat,jmat,1))*pf1(imat,jmat,ispin+8,isite-1)&
                                        !    +(1d0,0d0)*pf1(imat,jmat,ispin+8,isite)
                                    end do
                                end do
                            end do
                        end do
                       !$acc end kernels
                       !***********************
                       !*** Improved Action ***
                       !***********************
                    else if(nimprove.EQ.1)then
                        !$acc kernels
                        do jmat=1,nmat
                            tmpsum=0d0
                            do isite=1,nsite
                                do ispin=1,8
                                    do imat=1,nmat
                                        tmpsum=tmpsum&
                                            +2d0*dble(tmpfact*phase(imat,jmat,1)*chi(imat,jmat,ispin,isite+1,iremez,ipf)&
                                            *dconjg(pf1(imat,jmat,ispin+8,isite)))&
                                            +2d0*dble(tmpfact*dconjg(phase(imat,jmat,1))*chi(imat,jmat,ispin+8,isite-1,iremez,ipf)&
                                            *dconjg(pf1(imat,jmat,ispin,isite)))&
                                            -dble(tmpfact*phase(imat,jmat,2)*chi(imat,jmat,ispin,isite+2,iremez,ipf)&
                                            *dconjg(pf1(imat,jmat,ispin+8,isite)))&
                                            -dble(tmpfact*dconjg(phase(imat,jmat,2))*chi(imat,jmat,ispin+8,isite-2,iremez,ipf)&
                                            *dconjg(pf1(imat,jmat,ispin,isite)))&
                                            -2d0*dble(tmpfact*phase(jmat,imat,1)*chi(jmat,imat,ispin,isite+1,iremez,ipf)&
                                            *dconjg(pf1(jmat,imat,ispin+8,isite)))&
                                            -2d0*dble(tmpfact*dconjg(phase(jmat,imat,1))*chi(jmat,imat,ispin+8,isite-1,iremez,ipf)&
                                            *dconjg(pf1(jmat,imat,ispin,isite)))&
                                            +dble(tmpfact*phase(jmat,imat,2)*chi(jmat,imat,ispin,isite+2,iremez,ipf)&
                                            *dconjg(pf1(jmat,imat,ispin+8,isite)))&
                                            +dble(tmpfact*dconjg(phase(jmat,imat,2))*chi(jmat,imat,ispin+8,isite-2,iremez,ipf)&
                                            *dconjg(pf1(jmat,imat,ispin,isite)))

                                    !                            pf2(imat,jmat,ispin+8,isite)=&
                                    !                                pf2(imat,jmat,ispin+8,isite)&
                                    !                               -(0.5d0,0d0)*phase(imat,jmat,2)*pf1(imat,jmat,ispin,isite+2)&
                                    !                              +(2d0,0d0)*phase(imat,jmat,1)*pf1(imat,jmat,ispin,isite+1)&
                                    !                             -(1.5d0,0d0)*pf1(imat,jmat,ispin,isite)

                                    !                        pf2(imat,jmat,ispin,isite)=&
                                    !                           pf2(imat,jmat,ispin,isite)&
                                    !                          +(0.5d0,0d0)*dconjg(phase(imat,jmat,2))*pf1(imat,jmat,ispin+8,isite-2)&
                                    !                         -(2d0,0d0)*dconjg(phase(imat,jmat,1))*pf1(imat,jmat,ispin+8,isite-1)&
                                    !                        +(1.5d0,0d0)*pf1(imat,jmat,ispin+8,isite)
                                    end do
                                end do
                            end do
                            delh_alpha(jmat)=delh_alpha(jmat)+tmpsum
                        end do
                       !$acc end kernels
                    end if
                end if !ngauge
            end do !iremez
        end do !ipf
        return

    END SUBROUTINE Add_Force_Fermionic_device

end module hmc_force
