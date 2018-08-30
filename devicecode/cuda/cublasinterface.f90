! This module offers the possiblility to call the native CUBLAS functions
! in the multiplication of the Dirac operator. This might lead to a
! better preformance in rare cases. Georg Bergner

module cublasinterface
    use cudafor
    use openacc_cublas
    implicit none
    interface
        subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc ) &
            bind(c,name='cublasZgemm')
            use iso_c_binding
            integer(c_int), value :: m, n, k, lda, ldb, ldc
            complex(c_double_complex), device, dimension(m,n) :: a, b, c
            complex(c_double_complex), value :: alpha, beta
            character(kind=c_char), value :: transa, transb
        end subroutine zgemm
    end interface

    interface cublasSgemmBatched
        integer(c_int) function &
            cublasSgemmBatched_hpm(h, transa, transb, &
            m, n, k, alpha, Aarray, lda, Barray, ldb, &
            beta, Carray, ldc, batchCount) &
            bind(c,name='cublasSgemmBatched')
            use iso_c_binding
            use cublas
            type(cublasHandle), value :: h
            integer, value :: transa
            integer, value :: transb
            integer(c_int), value :: m, n, k
            real :: alpha
            type(c_devptr), device :: Aarray(*)
            integer(c_int), value :: lda
            type(c_devptr), device :: Barray(*)
            integer(c_int), value :: ldb
            real :: beta
            type(c_devptr), device :: Carray(*)
            integer(c_int), value :: ldc
            integer(c_int), value :: batchCount
        end function cublasSgemmBatched_hpm

        integer(c_int) function &
            cublasSgemmBatched_dpm(h, transa, transb, &
            m, n, k, alpha, Aarray, lda, Barray, ldb, &
            beta, Carray, ldc, batchCount) &
            bind(c,name='cublasSgemmBatched')
            use iso_c_binding
            use cublas
            type(cublasHandle), value :: h
            integer, value :: transa
            integer, value :: transb
            integer(c_int), value :: m, n, k
            real, device :: alpha
            type(c_devptr), device :: Aarray(*)
            integer(c_int), value :: lda
            type(c_devptr), device :: Barray(*)
            integer(c_int), value :: ldb
            real, device :: beta
            type(c_devptr), device :: Carray(*)
            integer(c_int), value :: ldc
            integer(c_int), value :: batchCount
        end function cublasSgemmBatched_dpm
    end interface cublasSgemmBatched

    interface
        integer(c_int) function &
            cublasZgemmBatched(h, transa, transb, &
            m, n, k, alpha, Aarray, lda, Barray, ldb, &
            beta, Carray, ldc, batchCount) &
            bind(c,name='cublasSgemmBatched')
            use iso_c_binding
            use cublas
            type(cublasHandle), value :: h
            integer, value :: transa
            integer, value :: transb
            integer(c_int), value :: m, n, k
            complex(c_double_complex) :: alpha
            type(c_devptr), device :: Aarray(*)
            integer(c_int), value :: lda
            type(c_devptr), device :: Barray(*)
            integer(c_int), value :: ldb
            complex(c_double_complex) :: beta
            type(c_devptr), device :: Carray(*)
            integer(c_int), value :: ldc
            integer(c_int), value :: batchCount
        end function cublasZgemmBatched
    end interface

contains

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
end module cublasinterface
