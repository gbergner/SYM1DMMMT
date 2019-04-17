!************************************************
!*** Calculate trx2 = (1/N)*¥int dt Tr(X_I^2) ***
!************************************************
module utils_measurements
    implicit none
contains
    subroutine Calc_TrX2_device(xmat,trx2)
  
        use compiletimeconstants
        implicit none

        double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        !$acc declare present(xmat)
        double precision trx2
  
        integer isite,idim
        integer imat,jmat

        trx2=0d0
        !$acc kernels
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
        !$acc end kernels
        trx2=trx2/dble(nmat*nsite)
  
        return
  
    END subroutine Calc_TrX2_device

    !trace part of (¥int dt X) and alpha are removed.
    SUBROUTINE subtract_U1_device(xmat,alpha)

        use compiletimeconstants
        implicit none

        !***** input & output *****
        double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double precision alpha(1:nmat)
        !$acc declare present(xmat,alpha)
        !****************************
        integer idim,imat,isite
        double complex trace(1:ndim),tmp
        double precision sum_alpha

        !***********************************************
        !**** trace part of (¥int dt X) is removed. ****
        !***********************************************
        do idim=1,ndim
            tmp=(0d0,0d0)
            !$acc kernels
            do isite=1,nsite
                do imat=1,nmat
                    tmp=tmp+xmat(imat,imat,idim,isite)
                end do
            end do
            !$acc end kernels
            tmp=tmp/dcmplx(nsite*nmat)
            !take care of the margin too.
            !$acc kernels
            do isite=-(nmargin-1),nsite+nmargin
                do imat=1,nmat
                    xmat(imat,imat,idim,isite)=xmat(imat,imat,idim,isite)-tmp
                end do
            end do
            !$acc end kernels
            trace(idim)=tmp/dcmplx(nsite*nmat)
        end do
        !*********************************************************
        !*** Trace part of alpha is removed. *********************
        !*** Do the same at all nodes, to avoid communication. ***
        !*********************************************************
        !$acc kernels
        sum_alpha=Sum(alpha)
        !$acc end kernels
        sum_alpha=sum_alpha/dble(nmat)
        !$acc kernels
        do imat=1,nmat
            alpha(imat)=alpha(imat)-sum_alpha
        end do
        !$acc end kernels
        return

    END SUBROUTINE subtract_U1_device
    SUBROUTINE hermitian_projection_device(xmat)

        use compiletimeconstants
        implicit none

        double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        !$acc declare present(xmat)
        integer idim,imat,jmat,isite

        !$acc kernels
        do idim=1,ndim
            do isite=-(nmargin-1),nsite+nmargin
                do imat=1,nmat-1
                    do jmat=imat+1,nmat
                        xmat(jmat,imat,idim,isite)=dconjg(xmat(imat,jmat,idim,isite))
                    end do
                end do
                do imat=1,nmat
                    xmat(imat,imat,idim,isite)=dcmplx(dble(xmat(imat,imat,idim,isite)))
                end do
            end do
        end do
        !$acc end kernels
        return

    END SUBROUTINE hermitian_projection_device

    SUBROUTINE Largest_eigenvalue_device(temperature,xmat,phase,Gam123,nbmn,neig,largest_eig,nbc)

        use compiletimeconstants
        use dirac_operator
        implicit none

        !***** input *****
        integer, intent(in) :: neig,nbc,nbmn
        double precision, intent(in) :: temperature
        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double complex, intent(in) :: phase(1:nmat,1:nmat,1:2)
        double complex, intent(in) :: Gam123(1:nspin,1:nspin)
        !***** output *****
        double precision largest_eig
        !******************
        double complex phi1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex phi2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        !$acc declare device_resident(phi1,phi2)
        double complex trace
        double precision norm
        integer imat,jmat
        integer ispin
        integer isite
        integer ieig,i
        double precision r1,r2
        !$acc declare present(temperature,xmat,nbmn,phase,Gam123)

        !***********************************
        !**** generate a random vector. ****
        !***********************************
        norm=0d0
        do imat=1,nmat
            do jmat=1,nmat
                do ispin=1,nspin
                    do isite=1,nsite
                        call BoxMuller(r1,r2)
                        phi1(imat,jmat,ispin,isite)=&
                            &(dcmplx(r1)+dcmplx(r2)*(0D0,1D0))/dcmplx(dsqrt(2d0))
                        norm=norm+(r1*r1+r2*r2)*0.5d0
                    end do
                end do
            end do
        end do
        !****************************
        !*** traceless projection ***
        !****************************
        do ispin=1,nspin
            do isite=1,nsite
                trace=(0d0,0d0)
                do imat=1,nmat
                    trace=trace+phi1(imat,imat,ispin,isite)
                end do
                trace=trace/dcmplx(nmat)
                do imat=1,nmat
                    norm=norm-dble(phi1(imat,imat,ispin,isite)*dconjg(phi1(imat,imat,ispin,isite)))
                    phi1(imat,imat,ispin,isite)=phi1(imat,imat,ispin,isite)-trace
                    norm=norm+dble(phi1(imat,imat,ispin,isite)*dconjg(phi1(imat,imat,ispin,isite)))
                end do

            end do
        end do

        norm=dsqrt(norm)
        norm=1d0/norm
        phi1=phi1*norm
        !*********************************
        !*** adjust the margin and b.c.***
        !*********************************
        call Adjust_margin_and_bc_pf(phi1,nbc)
        !******************************************
        !*** random vector has been generated.  ***
        !******************************************

        do ieig=1,neig
            !********************************************
            !*** phi1 -> phi2=D*phi -> phi1=D^dag*phi2***
            !********************************************
            call Multiply_Dirac(temperature,xmat,phase,Gam123,nbmn,&
                phi1,phi2)
            !*********************************
            !*** adjust the margin and b.c.***
            !*********************************
            call Adjust_margin_and_bc_pf_2(phi2,nbc)

            call Multiply_Dirac_dagger(temperature,xmat,phase,Gam123,nbmn,&
                phi2,phi1)
            !*********************************
            !*** adjust the margin and b.c.***
            !*********************************
            call Adjust_margin_and_bc_pf_2(phi1,nbc)
            !************************************x
            !*** calculate the norm of phi1.  ***
            !************************************
            norm=0d0
            do imat=1,nmat
                do jmat=1,nmat
                    do ispin=1,nspin
                        do isite=1,nsite
                            norm=norm&
                                &+dble(dconjg(phi1(imat,jmat,ispin,isite)&
                                &*dconjg(phi1(imat,jmat,ispin,isite))))
                        end do
                    end do
                end do
            end do
            norm=dsqrt(norm)
            norm=1d0/norm
            phi1=phi1*norm
            norm=1d0/norm
        end do

        largest_eig=norm

        return

    END SUBROUTINE Largest_eigenvalue_device

    ! Calculate the smallest eigenvalue of (D^dag*D),
    ! by multiplying (D^dag*D)^{-1} many times to a random vector.
    SUBROUTINE Smallest_eigenvalue_device(temperature,xmat,phase,Gam123,nbmn,neig,smallest_eig,nbc,max_err,max_iteration)

        use compiletimeconstants
        use cgm_solver
        implicit none

        !***** input *****
        integer, intent(in) :: neig,nbc,nbmn
        double precision, intent(in) :: temperature
        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double complex, intent(in) :: phase(1:nmat,1:nmat,1:2)
        double complex, intent(in) :: Gam123(1:nspin,1:nspin)
        !******************
        double complex phi1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex phi2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        !$acc declare device_resident(phi1,phi2)
        double complex trace
        double precision norm
        integer imat,jmat
        integer ispin
        integer isite
        integer ieig,i
        double precision r1,r2
        !$acc declare present(temperature,xmat,nbmn,phase,Gam123)
        integer info,max_iteration,iteration
        double precision max_err
        double precision smallest_eig
        !pf2=M*pf1, M: Dirac op

        !***********************************
        !**** generate a random vector. ****
        !***********************************
        norm=0d0
        do imat=1,nmat
            do jmat=1,nmat
                do ispin=1,nspin
                    do isite=1,nsite
                        call BoxMuller(r1,r2)
                        phi1(imat,jmat,ispin,isite)=&
                            &(dcmplx(r1)+dcmplx(r2)*(0D0,1D0))/dcmplx(dsqrt(2d0))
                        norm=norm+(r1*r1+r2*r2)*0.5d0
                    end do
                end do
            end do
        end do
        !****************************
        !*** traceless projection ***
        !****************************
        do ispin=1,nspin
            do isite=1,nsite
                trace=(0d0,0d0)
                do imat=1,nmat
                    trace=trace+phi1(imat,imat,ispin,isite)
                end do
                trace=trace/dcmplx(nmat)
                do imat=1,nmat
                    norm=norm-dble(phi1(imat,imat,ispin,isite)*dconjg(phi1(imat,imat,ispin,isite)))
                    phi1(imat,imat,ispin,isite)=phi1(imat,imat,ispin,isite)-trace
                    norm=norm+dble(phi1(imat,imat,ispin,isite)*dconjg(phi1(imat,imat,ispin,isite)))
                end do

            end do
        end do

        norm=dsqrt(norm)
        norm=1d0/norm
        phi1=phi1*norm
        !*********************************
        !*** adjust the margin and b.c.***
        !*********************************
        call Adjust_margin_and_bc_pf(phi1,nbc)
        !******************************************
        !*** random vector has been generated.  ***
        !******************************************
        do ieig=1,neig
            !*****************************************************
            !*** phi1 -> phi2=(D^dag*D)^{-1}*phi1 -> phi1=phi2 ***
            !*****************************************************
            call cg_solver_device(nbmn,nbc,temperature,&
                max_err,max_iteration,xmat,phase,Gam123,phi1,phi2,info,iteration)

            do imat=1,nmat
                do jmat=1,nmat
                    do ispin=1,nspin
                        do isite=-(nmargin-1),nsite+nmargin
                            phi1(imat,jmat,ispin,isite)=phi2(imat,jmat,ispin,isite)
                        end do
                    end do
                end do
            end do
            !****************************************
            !*** traceless projection;            ***
            !*** to avoid the zero mode to appear ***
            !*** as numerical artifact            ***
            !****************************************
            ! (Practically, it does not seem to be necessary.)
            !do ispin=1,nspin
            !   do isite=-nimprove,nsite_local+1+nimprove
            !      trace=(0d0,0d0)
            !      do imat=1,nmat
            !         trace=trace+phi1(imat,imat,ispin,isite)
            !      end do
            !      trace=trace/dcmplx(nmat)
            !      do imat=1,nmat
            !         phi1(imat,imat,ispin,isite)=phi1(imat,imat,ispin,isite)-trace
            !      end do
            !
            !   end do
            !end do

            !************************************
            !*** calculate the norm of phi1.  ***
            !************************************
            norm=0d0
            do imat=1,nmat
                do jmat=1,nmat
                    do ispin=1,nspin
                        do isite=1,nsite
                            norm=norm&
                                &+dble(dconjg(phi1(imat,jmat,ispin,isite)&
                                &*dconjg(phi1(imat,jmat,ispin,isite))))
                        end do
                    end do
                end do
            end do
            norm=dsqrt(norm)
            norm=1d0/norm
            phi1=phi1*norm
            norm=1d0/norm

           !write(*,*)1d0/norm
        end do

        smallest_eig=1d0/norm

        return

    END SUBROUTINE Smallest_eigenvalue_device

    SUBROUTINE save_checkpoint(xmat,alpha,itraj)
        use compiletimeconstants
        use mtmod !Mersenne twistor
        implicit none

        double precision alpha(1:nmat)
        double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        integer itraj,isite,idim,imat,jmat
        open(UNIT=42, File = "chp.dat", STATUS = "REPLACE", ACTION = "WRITE")
        call mtsaveu(42)
        write(42,*) itraj
        flush(42)
        !        close(42)
        !        open(UNIT=42, File = "chp.dat", STATUS = "REPLACE", ACTION = "WRITE",ACCESS='direct',recl=)
        do isite=1,nsite
            do idim=1,ndim
                do jmat=1,nmat
                    do imat=1,nmat
                        write(42,'(E22.16)') xmat(imat,jmat,idim,isite)%re
                        write(42,'(E22.16)') xmat(imat,jmat,idim,isite)%im
                    end do
                end do
            end do
        end do
        do imat=1,nmat
            write(42,'(E22.16)') alpha(imat)
        end do
        flush(42)
        close(42)
        return
    END SUBROUTINE save_checkpoint

    SUBROUTINE read_checkpoint(xmat,alpha,itraj,init,mersenne_seed)
        use compiletimeconstants
        use mtmod !Mersenne twistor
        use Adjust_margins
        implicit none
        !****** output ******
        double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        double precision alpha(1:nmat)
        double precision tmp1,tmp2
        integer itraj,isite,idim,imat,jmat
        integer stat,init
        integer mersenne_seed
        logical exist

        inquire(file="chp.dat",exist=exist)
        if(exist) then
            open(UNIT=42, File = "chp.dat", STATUS = "OLD", ACTION = "READ")
            call mtgetu(42)
            read(42,*) itraj
            do isite=1,nsite
                do idim=1,ndim
                    do jmat=1,nmat
                        do imat=1,nmat
                            read(42,'(E22.16)') tmp1
                            read(42,'(E22.16)') tmp2
                            xmat(imat,jmat,idim,isite)=CMPLX(tmp1,tmp2)
                        end do
                    end do
                end do
            end do
            do imat=1,nmat
                read(42,'(E22.16)') alpha(imat)
            end do
            close(42)

            if(init.NE.4)then
                call sgrnd(mersenne_seed)
                print*,"Reset random seed."
            end if
        else
            print*, "WARNING: FILE chp.dat DOES NOT EXIST, NO CHECKPOINT"
            print*, "RESTART WITH MERSENNE SEED FROM CONFIG FILE"
            call sgrnd(mersenne_seed)
        end if
        !call Adjust_margin_xmat_device(xmat) This must be called after update of device
    END SUBROUTINE read_checkpoint

end module utils_measurements
