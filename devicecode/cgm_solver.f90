! This is a multishift-CG, it is considerably differenent from the code in the
! serial host branch, but the result should be compatible.
! Despite the name in the serial host branch, the code is also a CG and not a BiCG.
! Georg Bergner

module cgm_solver
    implicit none
contains
    SUBROUTINE norm_vect_device(norm,pf)

        use compiletimeconstants
        implicit none

        double complex,intent(in) :: pf(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double precision, intent(out) :: norm
        integer :: isite,ispin,jspin,imat,jmat
        !$acc declare device_resident(pf)
        norm=0d0
        !$acc kernels
        do isite=1,nsite
            do ispin=1,nspin
                do jmat=1,nmat
                    do imat=1,nmat
                        norm=norm+real(pf(imat,jmat,ispin,isite)&
                            *dconjg((pf(imat,jmat,ispin,isite))))
                    end do
                end do
            end do
        end do
       !$acc end kernels
    END SUBROUTINE norm_vect_device


    SUBROUTINE check_cgm_solutions(nremez,bcoeff,nbmn,nbc,temperature,&
        xmat,phase,Gam123,pf_input,chi)

        use compiletimeconstants
        use dirac_operator
        implicit none

        integer, intent(in) :: nremez
        integer, intent(in) :: nbmn,nbc
        double precision,intent(in) :: bcoeff(1:nremez)
        double precision, intent(in) :: temperature
        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double complex, intent(in) :: phase(1:nmat,1:nmat,1:2)
        double complex, intent(in) :: Gam123(1:nspin,1:nspin)
        double complex, intent(in) :: pf_input(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
        double complex, intent(in) :: chi(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:nremez,1:npf)
        !$acc declare present(nbmn,nbc,bcoeff,temperature,xmat)
        !$acc declare device_resident(phase,Gam123,pf_input,chi)
        double complex ::  mp1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
            mp2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        !$acc declare device_resident(mp1,mp2)
        double precision :: tmpsum, bc, sumcompl
   
        integer imat,jmat
        integer ispin
        integer isite
        integer iremez
        integer :: ipf

        sumcompl=0d0

        do ipf=1,npf
            do iremez=1,nremez
                bc=bcoeff(iremez)
                  !$acc kernels
                do isite=1,nsite
                    do ispin=1,nspin
                        do jmat=1,nmat
                            do imat=1,nmat
                                mp1(imat,jmat,ispin,isite)=chi(imat,jmat,ispin,isite,iremez,ipf)
                            end do
                        end do
                    end do
                end do
                !$acc end kernels
                call set_boundary_device(nbc,mp1)
                call Multiply_Dirac_device(temperature,xmat,phase,Gam123,nbmn,mp1,mp2)
                call set_boundary_device(nbc,mp2)
                call Multiply_Dirac_dagger_device(temperature,xmat,phase,Gam123,nbmn,mp2,mp1)
     


                tmpsum=0.0d0
                !$acc kernels
                do isite=1,nsite
                    do ispin=1,nspin
                        do jmat=1,nmat
                            do imat=1,nmat
                                tmpsum=tmpsum+abs(mp1(imat,jmat,ispin,isite)+bc*chi(imat,jmat,ispin,isite,iremez,ipf)-pf_input(imat,jmat,ispin,isite,ipf))
                            end do
                        end do
                    end do
                end do
                !$acc end kernels
                write (*,*) "Error for ipf ",ipf," iremez ",iremez," is ",tmpsum
                sumcompl=sumcompl+tmpsum
            end do !loop nremez
        end do !loop npf
        write (*,*) "Error for complete inversion(sum) ",sumcompl
        if(sumcompl.GE.0.01d0) then
            write (*,*) "STOP: ERROR: large miss in CGM inversion!"
            call exit(-1)
        end if
    END SUBROUTINE check_cgm_solutions

    SUBROUTINE cgm_solver_device(nremez,bcoeff,nbmn,nbc,temperature,&
        max_err,max_iteration,xmat,phase,Gam123,pf_input,chi,info,iteration)
        use compiletimeconstants
        use dirac_operator
        use cublasinterface
        implicit none

        integer, intent(in) :: nremez
        integer, intent(in) :: nbmn,nbc
        double precision,intent(in) :: bcoeff(1:nremez)
        double precision, intent(in) :: temperature
        integer, intent(in) :: max_iteration
        double precision, intent(in) :: max_err
        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double complex, intent(in) :: phase(1:nmat,1:nmat,1:2)
        double complex, intent(in) :: Gam123(1:nspin,1:nspin)
        integer, intent(out) :: info
    
        double complex, intent(in) :: pf_input(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
        double complex, intent(out) :: chi(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:nremez,1:npf)
        !$acc declare present(nbmn,nbc,bcoeff,temperature,xmat)
        !$acc declare device_resident(phase,Gam123,pf_input,chi)

        double complex ::  ps(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:nremez),&
            r(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
            p(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
            mp1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
            mp2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        !$acc declare device_resident(ps,r,p,mp1,mp2)
       
        integer imat,jmat
        integer ispin
        integer isite
        integer iremez
        integer,intent(out):: iteration
        integer :: ipf
        double precision :: normPF2,NormR
        double precision :: numerator,denominator
        double precision :: a,aold
        double precision :: bc,bsc,asc
        double precision  :: Zsc(0:2,1:nremez)
        !$acc declare device_resident(Zsc)
        double precision :: error(1:nremez)
        integer :: hantei(1:nremez)

        type(c_devptr), device :: xptr_d(nsite,ndim)
        type(c_devptr), device :: pptr_d(nsite,nspin)
        type(c_devptr), device :: mp1ptr_d(nsite,nspin)
        type(c_devptr), device :: mp2ptr_d(nsite,nspin)
  

        if(cublascgm==1) then
            call setup_cublas()
            call setup_cublas_pointers_xmat(xmat,xptr_d)
            call setup_cublas_pointers_pf(p,pptr_d)
            call setup_cublas_pointers_pf(mp1,mp1ptr_d)
            call setup_cublas_pointers_pf(mp2,mp2ptr_d)
        end if

        do ipf=1,npf

            call norm_vect_device(normPF2,pf_input(:,:,:,:,ipf))
            normPF2=dsqrt(normPF2)
            if(solver_verbose.EQ.3) then
                write(*,*) "norm device input vector (",ipf,")",normPF2
            end if
            !initial condition
            !$acc kernels
            do iremez=1,nremez
                do isite=1,nsite
                    do ispin=1,nspin
                        do jmat=1,nmat
                            do imat=1,nmat
                                chi(imat,jmat,ispin,isite,iremez,ipf)=(0d0,0d0)
                                ps(imat,jmat,ispin,isite,iremez)=pf_input(imat,jmat,ispin,isite,ipf)
                            end do
                        end do
                    end do
                end do
            end do
            !$acc end kernels
            !$acc kernels
            do isite=1,nsite
                do ispin=1,nspin
                    do jmat=1,nmat
                        do imat=1,nmat
                            r(imat,jmat,ispin,isite)=pf_input(imat,jmat,ispin,isite,ipf)
                            p(imat,jmat,ispin,isite)=pf_input(imat,jmat,ispin,isite,ipf)
                        end do
                    end do
                end do
            end do
             !$acc end kernels
             !$acc kernels
            do iremez=1,nremez
                Zsc(0,iremez)=(1d0,0d0)
                Zsc(1,iremez)=(1d0,0d0)
            end do
            !$acc end kernels
            aold=(1d0,0d0)
            bc=(1d0,0d0)
            !loop
            iteration=1
            info=1
            hantei=1
            call norm_vect_device(numerator,r)
     
            !$acc data copyin(hantei,a,aold,bc,error)
            do while((iteration.LT.max_iteration).AND.(info.EQ.1))
                call set_boundary_device(nbc,p)
                if(cublascgm==1) then
                    call Multiply_Dirac_device_cuda(temperature,xmat,phase,Gam123,nbmn,p,mp1,xptr_d,pptr_d,mp1ptr_d)
                else
                    call Multiply_Dirac_device(temperature,xmat,phase,Gam123,nbmn,p,mp1)
                end if
                call norm_vect_device(denominator,mp1)
        
                a=numerator/denominator
                !$acc update device(a)
                !ps update
                !$acc wait(2)
                !hantei update device
                !$acc wait(1)
                !update Xs
                !$acc kernels async(1)
                !$acc loop independent private(asc)
                do iremez=1,nremez
                    if(hantei(iremez).EQ.1)then
                        Zsc(2,iremez)=&
                            Zsc(1,iremez)*Zsc(0,iremez)*aold&
                            /(aold*Zsc(0,iremez)*(1d0&
                            +a*bcoeff(iremez))&
                            +a*bc*(Zsc(0,iremez)&
                            -Zsc(1,iremez)))
                        asc=a*Zsc(2,iremez)&
                            /Zsc(1,iremez)
                        do isite=1,nsite
                            do ispin=1,nspin
                                do jmat=1,nmat
                                    do imat=1,nmat
                                        chi(imat,jmat,ispin,isite,iremez,ipf)=&
                                            chi(imat,jmat,ispin,isite,iremez,ipf)&
                                            +asc&
                                            *ps(imat,jmat,ispin,isite,iremez)
                                    end do
                                end do
                            end do
                        end do
                    end if
                end do
                !$acc end kernels

                call set_boundary_device(nbc,mp1)
                if(cublascgm==1) then
                    call Multiply_Dirac_dagger_device_cuda(temperature,xmat,phase,Gam123,nbmn,mp1,mp2,xptr_d,mp1ptr_d,mp2ptr_d)
                else
                    call Multiply_Dirac_dagger_device(temperature,xmat,phase,Gam123,nbmn,mp1,mp2)
                end if

                denominator=numerator
                numerator=(0d0,0d0)
                !$acc kernels
                do isite=1,nsite
                    do ispin=1,nspin
                        do jmat=1,nmat
                            do imat=1,nmat
                                r(imat,jmat,ispin,isite)=r(imat,jmat,ispin,isite)-a*mp2(imat,jmat,ispin,isite)
                                numerator=numerator&
                                    +real(r(imat,jmat,ispin,isite)&
                                    *dconjg(r(imat,jmat,ispin,isite)))
                            end do
                        end do
                    end do
                end do
                !$acc end kernels
 
                bc=numerator/denominator
                 !$acc update device(bc)

                NormR=dsqrt(numerator)
                !xs and zsc update
                !$acc wait(1)
                !$acc kernels
                do iremez=1,nremez
                    if(hantei(iremez).EQ.1)then
                        ERROR(iremez)=Zsc(2,iremez)&
                            *NormR/normPF2
                    end if
                end do
                !$acc end kernels
                !$acc update host(error) async(1)
         
                !$acc kernels
                do isite=1,nsite
                    do ispin=1,nspin
                        do jmat=1,nmat
                            do imat=1,nmat
                                p(imat,jmat,ispin,isite)=r(imat,jmat,ispin,isite)+bc*p(imat,jmat,ispin,isite)
                            end do
                        end do
                    end do
                end do
                !$acc end kernels
        
                !complete xs update and Zsc calculation


                !$acc kernels async(2)
                !$acc loop independent private(bsc)
                do iremez=1,nremez
                    if(hantei(iremez).EQ.1)then
                        bsc=bc&
                            *(Zsc(2,iremez)&
                            /Zsc(1,iremez))*(Zsc(2,iremez)&
                            /Zsc(1,iremez))
                        do isite=1,nsite
                            do ispin=1,nspin
                                do jmat=1,nmat
                                    do imat=1,nmat
                                        ps(imat,jmat,ispin,isite,iremez)&
                                            =Zsc(2,iremez)*r(imat,jmat,ispin,isite)&
                                            +bsc*ps(imat,jmat,ispin,isite,iremez)
                                    end do
                                end do

                            end do
                        end do
                    end if
                end do
                !$acc end kernels
                !
                !         !Evaluation of err.
                !         !**evaluate ERROR for each iremez.
                !         !**if ERROR is smaller than Max_Err, then stop updating Zs.
                !         !**This procedure is crucial!!

                !complete error update
                !$acc wait(1)
                do iremez=1,nremez
                    if((ERROR(iremez).LT.max_err).AND.(hantei(iremez).EQ.1))then
                        hantei(iremez)=0
                        if(solver_verbose.EQ.2) then
                            write(*,*)  "finished value ",ipf," ",iremez," ", ERROR(iremez), " ",Zsc(2,iremez), " ",NormR," ",normPF2
                        end if
                    end if
                end do
                !$acc update device(hantei) async(1)
                if(solver_verbose.EQ.3) then
                    write(*,*)  "iteration",iteration," nromR ",NormR
                    write(*,*)  "iteration",iteration," errorfull",Error
                    write(*,*) "hantei",hantei
                    NormR=0.0
                    !$acc kernels
                    do isite=1,nsite
                        do ispin=1,nspin
                            do jmat=1,nmat
                                do imat=1,nmat
                                    NormR=NormR+chi(imat,jmat,ispin,isite,1,ipf)
                                end do
                            end do
                        end do
                    end do
                    !$acc end kernels

                    write(*,*) "xsnorm device",NormR
                end if
                iteration=iteration+1
                info=0
                do iremez=1,nremez
                    if(hantei(iremez).EQ.1)then
                        info=1
                    end if
                end do
                aold=a
                !$acc update  device(aold)
                !$acc kernels
                Zsc(0,:)=Zsc(1,:)
                Zsc(1,:)=Zsc(2,:)
                !$acc end kernels
            end do !loop iteration

               !$acc wait
               !$acc end data

              !boundary condition
            do iremez=1,nremez
                call set_boundary_device(nbc,chi(:,:,:,:,iremez,ipf))
                if(solver_verbose.EQ.2) then
                    write(*,*)  "error after inversion=",ERROR(iremez)," ",iremez," ",ipf
                end if
            end do

        end do !loop pf fields
        if(solver_verbose.EQ.1) then
            write(*,*)  "check after iteration=",iteration
            call check_cgm_solutions(nremez,bcoeff,nbmn,nbc,temperature,&
                xmat,phase,Gam123,pf_input,chi)
        end if
        if(cublascgm==1) then
            call finish_cublas()
        end if
        return

    END SUBROUTINE cgm_solver_device

    SUBROUTINE cg_solver_device(nbmn,nbc,temperature,&
        max_err,max_iteration,xmat,phase,Gam123,pf_input,pf_sol,info,iteration)
        use compiletimeconstants
        use dirac_operator
        use cublasinterface
        implicit none

        integer, intent(in) :: nbmn,nbc
        double precision, intent(in) :: temperature
        integer, intent(in) :: max_iteration
        double precision, intent(in) :: max_err
        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double complex, intent(in) :: phase(1:nmat,1:nmat,1:2)
        double complex, intent(in) :: Gam123(1:nspin,1:nspin)
        integer, intent(out) :: info

        double complex, intent(in) :: pf_input(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex, intent(out) :: pf_sol(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        !$acc declare present(nbmn,nbc,temperature,xmat)
        !$acc declare device_resident(phase,Gam123,pf_input,pf_sol)

        double complex :: r(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
            p(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
            mp1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
            mp2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        !$acc declare device_resident(r,p,mp1,mp2)

        integer imat,jmat
        integer ispin
        integer isite
        integer,intent(out):: iteration
        double precision :: normPF2,NormR
        double precision :: numerator,denominator
        double precision :: a,aold
        double precision :: bc
        double precision :: error

        type(c_devptr), device :: xptr_d(nsite,ndim)
        type(c_devptr), device :: pptr_d(nsite,nspin)
        type(c_devptr), device :: mp1ptr_d(nsite,nspin)
        type(c_devptr), device :: mp2ptr_d(nsite,nspin)


        if(cublascgm==1) then
            call setup_cublas()
            call setup_cublas_pointers_xmat(xmat,xptr_d)
            call setup_cublas_pointers_pf(p,pptr_d)
            call setup_cublas_pointers_pf(mp1,mp1ptr_d)
            call setup_cublas_pointers_pf(mp2,mp2ptr_d)
        end if

        call norm_vect_device(normPF2,pf_input)
        normPF2=dsqrt(normPF2)
        if(solver_verbose.EQ.3) then
            write(*,*) "norm device input vector ",normPF2
        end if
        !initial condition
        !$acc kernels
        pf_sol=(0d0,0d0)
        p=pf_input
        r=pf_input
        !$acc end kernels
        aold=(1d0,0d0)
        bc=(1d0,0d0)
        !loop
        iteration=1
        info=1
        call norm_vect_device(numerator,r)

        !$acc data copyin(a,aold,bc,error)
        do while((iteration.LT.max_iteration).AND.(info.EQ.1))
            call set_boundary_device(nbc,p)
            if(cublascgm==1) then
                call Multiply_Dirac_device_cuda(temperature,xmat,phase,Gam123,nbmn,p,mp1,xptr_d,pptr_d,mp1ptr_d)
            else
                call Multiply_Dirac_device(temperature,xmat,phase,Gam123,nbmn,p,mp1)
            end if
            call norm_vect_device(denominator,mp1)

            a=numerator/denominator
            !$acc update device(a)

            call set_boundary_device(nbc,mp1)
            if(cublascgm==1) then
                call Multiply_Dirac_dagger_device_cuda(temperature,xmat,phase,Gam123,nbmn,mp1,mp2,xptr_d,mp1ptr_d,mp2ptr_d)
            else
                call Multiply_Dirac_dagger_device(temperature,xmat,phase,Gam123,nbmn,mp1,mp2)
            end if

            denominator=numerator
            numerator=(0d0,0d0)
            !$acc kernels
            do isite=1,nsite
                do ispin=1,nspin
                    do jmat=1,nmat
                        do imat=1,nmat
                            pf_sol(imat,jmat,ispin,isite)=pf_sol(imat,jmat,ispin,isite)+a*p(imat,jmat,ispin,isite)
                            r(imat,jmat,ispin,isite)=r(imat,jmat,ispin,isite)-a*mp2(imat,jmat,ispin,isite)
                            numerator=numerator&
                                +real(r(imat,jmat,ispin,isite)&
                                *dconjg(r(imat,jmat,ispin,isite)))
                        end do
                    end do
                end do
            end do
            !$acc end kernels

            bc=numerator/denominator
            !$acc update device(bc)

            NormR=dsqrt(numerator)

            error=NormR/normPF2
            if((error.LT.max_err))then
                info=0
            end if


            !$acc kernels
            do isite=1,nsite
                do ispin=1,nspin
                    do jmat=1,nmat
                        do imat=1,nmat
                            p(imat,jmat,ispin,isite)=r(imat,jmat,ispin,isite)+bc*p(imat,jmat,ispin,isite)
                        end do
                    end do
                end do
            end do
            !$acc end kernels

            if(solver_verbose.EQ.3) then
                write(*,*)  "iteration",iteration," nromR ",NormR
                write(*,*)  "iteration",iteration," errorfull",Error
                NormR=0.0
                !$acc kernels
                do isite=1,nsite
                    do ispin=1,nspin
                        do jmat=1,nmat
                            do imat=1,nmat
                                NormR=NormR+pf_sol(imat,jmat,ispin,isite)
                            end do
                        end do
                    end do
                end do
                !$acc end kernels

                write(*,*) "xsnorm device",NormR
            end if
            iteration=iteration+1
            aold=a
            !$acc update  device(aold)
        end do !loop iteration

            !$acc end data
        call set_boundary_device(nbc,pf_sol)
          !boundary condition

        if(solver_verbose.EQ.2) then
            write(*,*)  "error after inversion=",error
        end if

        if(solver_verbose.EQ.1) then
            write(*,*)  "check after iteration=",iteration
            call Multiply_Dirac_device(temperature,xmat,phase,Gam123,nbmn,pf_sol,mp1)
            call set_boundary_device(nbc,mp1)
            call Multiply_Dirac_dagger_device(temperature,xmat,phase,Gam123,nbmn,mp1,mp2)
            error=0.0d0
            !$acc kernels
            do isite=1,nsite
                do ispin=1,nspin
                    do jmat=1,nmat
                        do imat=1,nmat
                            error=error+abs(mp2(imat,jmat,ispin,isite)-pf_input(imat,jmat,ispin,isite))
                        end do
                    end do
                end do
            end do
            !$acc end kernels
            write (*,*) "Error for is ",error
        end if
        if(cublascgm==1) then
            call finish_cublas()
        end if
        return

    END SUBROUTINE cg_solver_device

end module cgm_solver
