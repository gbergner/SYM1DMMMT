! This module offers the possiblility to call the native CUBLAS functions
! in the multiplication of the Dirac operator. This might lead to a
! better preformance in rare cases. Georg Bergner

module cublasinterface
    use cudafor
    use openacc_cublas
    implicit none

    save
     type(cublasHandle), value :: cublas_handle
     integer :: cublas_status
     type(cudaEvent) :: startEvent, stopEvent
     integer :: timer_status

contains
    subroutine cudaSetupTimer()
      use cudafor
      implicit none
      timer_status = cudaEventCreate(startEvent)
      timer_status = cudaEventCreate(stopEvent)
    end subroutine cudaSetupTimer

    subroutine cudaFinishTimer()
      use cudafor
      implicit none
      timer_status = cudaEventDestroy(startEvent)
      timer_status = cudaEventDestroy(stopEvent)
    end subroutine cudaFinishTimer

    subroutine cudaTimerStart()
      use cudafor
      implicit none
      timer_status = cudaEventRecord(startEvent, 0)
    end subroutine cudaTimerStart

    subroutine cudaTimerStop(time)
          use cudafor
      implicit none
      real,intent(out) :: time
      timer_status = cudaEventRecord(stopEvent, 0)
      timer_status = cudaEventSynchronize(stopEvent)
      timer_status = cudaEventElapsedTime(time, startEvent, stopEvent)
    end subroutine cudaTimerStop

    subroutine  callZgemmBatched(pf2,xmat,pf1,alpha,beta,ispin,jspin,idim)
        use cudafor
        use cublas
        use iso_c_binding
        use compiletimeconstants
        implicit none
        integer, intent(in) :: ispin,jspin,idim
        complex(c_double_complex), value,intent(in) :: alpha, beta
        double complex, intent(in) :: pf1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex, intent(out) :: pf2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        !$acc declare device_resident(pf1,pf2)
        !$acc declare present(xmat)

        type(c_devptr) :: Aptr(nsite)
        type(c_devptr), device :: Aptr_d(nsite)
        type(c_devptr) :: Bptr(nsite)
        type(c_devptr), device :: Bptr_d(nsite)
        type(c_devptr) :: Cptr(nsite)
        type(c_devptr), device :: Cptr_d(nsite)
        type(cublasHandle) :: h1
        integer :: isite, istat
        ! create array of device pointers on host and copy it to device
        do isite = 1,nsite
            Aptr(isite) = c_devloc(xmat(1,1,idim,isite))
            Bptr(isite) = c_devloc(pf1(1,1,jspin,isite))
            Cptr(isite) = c_devloc(pf2(1,1,ispin,isite))
        end do
        Aptr_d=Aptr
        Bptr_d=Bptr
        Cptr_d=Cptr
        istat = cublasCreate(h1)
        if (istat /= CUBLAS_STATUS_SUCCESS) write(*,*) 'cublasCreate failed'
            !C [ i ] = α op ( A [ i ] ) op ( B [ i ] ) + β C [ i ] ,  for i  ∈ [ 0 , b a t c h C o u n t − 1 ]
        istat= cublasZgemmBatched(h1, 'n', 'n', &
            nmat, nmat, nmat, -alpha, Aptr_d, nmat, Bptr_d, nmat, &
            beta, Cptr_d, nmat, nsite)
        istat= cublasZgemmBatched(h1, 'n', 'n', &
            nmat, nmat, nmat, alpha, Bptr_d, nmat, Aptr_d, nmat, &
            beta, Cptr_d, nmat, nsite)
        if (istat /= CUBLAS_STATUS_SUCCESS) write(*,*) 'cublas failed: ', istat
        istat = cublasDestroy(h1)
        if (istat /= CUBLAS_STATUS_SUCCESS) write(*,*) 'cublasDestroy failed'
    end subroutine

    subroutine  setup_cublas_pointers_xmat(xmat,xptr_d)
        use cudafor
        use cublas
        use iso_c_binding
        use compiletimeconstants
        implicit none
        double complex, intent(in) :: xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        !$acc declare present(xmat)

        type(c_devptr) :: Aptr(nsite)
        type(c_devptr), device, intent(out) :: xptr_d(nsite,ndim)
        integer :: isite, idim
        ! create array of device pointers on host and copy it to device
        do idim = 1,ndim
            do isite = 1,nsite
                Aptr(isite) = c_devloc(xmat(1,1,idim,isite))
            end do
            xptr_d(:,idim)=Aptr
        end do
    end subroutine

    subroutine  setup_cublas_pointers_pf(pf,pfptr_d)
        use cudafor
        use cublas
        use iso_c_binding
        use compiletimeconstants
        implicit none
        double complex, intent(in) :: pf(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
        !$acc declare device_resident(pf)
        type(c_devptr) :: Aptr(nsite)
        type(c_devptr), device, intent(out) :: pfptr_d(nsite,nspin)
        integer :: isite, ispin
        ! create array of device pointers on host and copy it to device
        do ispin = 1,nspin
            do isite = 1,nsite
                Aptr(isite) = c_devloc(pf(1,1,ispin,isite))
            end do
            pfptr_d(:,ispin)=Aptr
        end do
    end subroutine

    subroutine setup_cublas()
      use cublas
      implicit none
      cublas_status = cublasCreate(cublas_handle)
        if (cublas_status /= CUBLAS_STATUS_SUCCESS) write(*,*) 'cublasCreate failed'
    end subroutine setup_cublas

    subroutine finish_cublas()
      use cublas
      implicit none
      cublas_status = cublasDestroy(cublas_handle)
        if (cublas_status /= CUBLAS_STATUS_SUCCESS) write(*,*) 'cublasDestroy failed'
    end subroutine finish_cublas

    subroutine check_cublas()
      implicit none
      if (cublas_status /= CUBLAS_STATUS_SUCCESS) write(*,*) 'cublas error'
    end subroutine check_cublas

    subroutine  multiply_cublas_pointer(Aptr_d,Bptr_d,Cptr_d,alpha,beta,ispin,jspin,idim)
        use cudafor
        use cublas
        use iso_c_binding
        use compiletimeconstants
        implicit none
        integer, intent(in) :: ispin,jspin,idim
        complex(c_double_complex), value,intent(in) :: alpha, beta
        type(c_devptr), device :: Aptr_d(nsite,ndim)
        type(c_devptr), device :: Bptr_d(nsite,nspin)
        type(c_devptr), device :: Cptr_d(nsite,nspin)
            !C [ i ] = α op ( A [ i ] ) op ( B [ i ] ) + β C [ i ] ,  for i  ∈ [ 0 , b a t c h C o u n t − 1 ]
        cublas_status= cublasZgemmBatched(cublas_handle, 'n', 'n', &
            nmat, nmat, nmat, -alpha, Aptr_d(:,idim), nmat, Bptr_d(:,jspin), nmat, &
            beta, Cptr_d(:,ispin), nmat, nsite)
        cublas_status=  cublasZgemmBatched(cublas_handle, 'n', 'n', &
            nmat, nmat, nmat, alpha, Bptr_d(:,jspin), nmat, Aptr_d(:,idim), nmat, &
            beta, Cptr_d(:,ispin), nmat, nsite)
        !cublas_status=cudaDeviceSynchronize()
    end subroutine

    subroutine  multiply_cublas_pointer_streams(Aptr_d,Bptr_d,Cptr_d,prefact,dag)
        use cudafor
        use cublas
        use iso_c_binding
        use compiletimeconstants
        use gammamatrix
        implicit none
        complex(c_double_complex), value :: alpha, beta
        type(c_devptr), device :: Aptr_d(nsite,ndim)
        type(c_devptr), device :: Bptr_d(nsite,nspin)
        type(c_devptr), device :: Cptr_d(nsite,nspin)
        integer(kind=cuda_stream_kind) :: streams(nspin)
        integer(kind=cuda_stream_kind) :: oldstream
        double complex :: gamtmp,prefact
        integer :: imat,jmat,kmat
        integer :: idim
        integer :: ispin,jspin,kspin
        integer :: isite
        integer :: countr
        integer :: status
        logical,intent(in) :: dag
        status=cublasGetStream(cublas_handle,oldstream)
        do ispin=1,nspin
          status = cudaStreamCreate(streams(ispin))
        end do
        do countr=1,144
            ispin=gamispin(countr)
            jspin=gamjspin(countr)
            idim=gamidim(countr)
            if(dag) then
             gamtmp=dconjg(gamgam(countr))
            else
             gamtmp=gamgam(countr)
            end if
            status=cublasSetStream(cublas_handle,streams(ispin))
            alpha=gamtmp*prefact
            beta=dcmplx(1.d0)
           cublas_status= cublasZgemmBatched(cublas_handle, 'n', 'n', &
            nmat, nmat, nmat, -alpha, Aptr_d(:,idim), nmat, Bptr_d(:,jspin), nmat, &
            beta, Cptr_d(:,ispin), nmat, nsite)
           cublas_status=  cublasZgemmBatched(cublas_handle, 'n', 'n', &
            nmat, nmat, nmat, alpha, Bptr_d(:,jspin), nmat, Aptr_d(:,idim), nmat, &
            beta, Cptr_d(:,ispin), nmat, nsite)
        end do
        do ispin=1,nspin
           status = cudaStreamSynchronize(streams(ispin))
        end do
        do ispin=1,nspin
          status = cudaStreamDestroy(streams(ispin))
        end do
        status = cublasSetStream(cublas_handle,oldstream)
    end subroutine
end module cublasinterface
