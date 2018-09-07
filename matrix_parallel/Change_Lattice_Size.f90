!##############################################################################
!######              BFSS matrix model on lattice                     #########
!######                                                               #########
!######                 written by Masanori Hanada                    #########
!######                                                               #########
!######            ver.1 Change the lattice size.                     ######### 
!##############################################################################
!Mersenne twister
include 'mt19937.f90'

program ChangeLatticeSize

  use mtmod
  implicit none
  !---------------------------------
  integer nsite_new,nsite_old,nmat,ndim
  parameter(ndim=9)!do not change
  parameter(nmat=3)
  parameter(nsite_new=8)
  parameter(nsite_old=4)

  integer itraj

  double complex xmat_old(1:nmat,1:nmat,1:ndim,1:nsite_old)
  double complex xmat_new(1:nmat,1:nmat,1:ndim,1:nsite_new)
  double precision alpha(1:nmat)
  character(150)input_config,output_config
  !*************************
  !**** read parameters ****
  !*************************
  open(unit=10,status='OLD',file='input_parallel_v1.dat',action='READ')
  read(10,*) input_config
  read(10,*) output_config
  close(10)
  !*************************************
  !*************************************
  !*** Set the initial configuration ***
  !*************************************
  !*************************************

 
  open(unit=9,status='OLD',file=input_config,action='READ')
  call mtgetu(9)
  read(9,*) itraj
  read(9,*) xmat_old
  read(9,*) alpha
  close(9)
  
  call Fourier_transform(nmat,ndim,nsite_old,nsite_new,xmat_old,xmat_new)
  
  open(UNIT = 22, File = output_config, STATUS = "REPLACE", ACTION = "WRITE")
  call mtsaveu(22)  
  write(22,*) itraj
  write(22,*) xmat_new
  write(22,*) alpha
  close(22)
  
  !*******************************************
  !*** It's the end, have a good weekend! ****
  !*******************************************


end program ChangeLatticeSize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Fourier_transform(nmat,ndim,nsite_old,nsite_new,xmat_old,xmat_new)
  
  implicit none
  

  integer nmat,ndim,nsite_new,nsite_old
  double complex xmat_old(1:nmat,1:nmat,1:ndim,1:nsite_old)
  double complex xmat_mom_old(1:nmat,1:nmat,1:ndim,1:nsite_old)
  double complex xmat_new(1:nmat,1:nmat,1:ndim,1:nsite_new)
  double complex xmat_mom_new(1:nmat,1:nmat,1:ndim,1:nsite_new)
  
  double complex phase
  double precision temp,pi
  
  integer isite,idim,imom
  integer imat,jmat

  pi=2d0*dasin(1d0)
  
  
  xmat_mom_old=(0d0,0d0)
  do idim=1,ndim
!$omp parallel
!$omp do
     do imat=1,nmat
        do jmat=1,nmat
           do imom=1,nsite_old
              do isite=1,nsite_old
                 temp=2d0*pi*dble(imom*isite)/dble(nsite_old)
                 phase=dcmplx(dcos(temp))-(0d0,1d0)*dcmplx(dsin(temp))
                 phase=phase/dcmplx(dsqrt(dble(nsite_old)))
                 xmat_mom_old(imat,jmat,idim,imom)=&
                      xmat_mom_old(imat,jmat,idim,imom)&
                      +xmat_old(imat,jmat,idim,isite)*phase
              end do
           end do
        end do
     end do
!$omp end do
!$omp end parallel
  end do

  xmat_mom_new=(0d0,0d0)
  do idim=1,ndim
!$omp parallel
!$omp do
     do imat=1,nmat
        do jmat=1,nmat
           do imom=1,int(dble(nsite_old)/2d0+0.1d0)
              xmat_mom_new(imat,jmat,idim,imom)=&
                   xmat_mom_old(imat,jmat,idim,imom)
              xmat_mom_new(imat,jmat,idim,nsite_new-imom)=&
                   dconjg(xmat_mom_old(jmat,imat,idim,imom))
                   
           end do

           xmat_mom_new(imat,jmat,idim,nsite_new)=&
                xmat_mom_old(imat,jmat,idim,nsite_old)
        end do
     end do
!$omp end do
!$omp end parallel
  end do

  xmat_new=(0d0,0d0)
  do idim=1,ndim
!$omp parallel
!$omp do
     do imat=1,nmat
        do jmat=1,nmat
           do imom=1,nsite_new
              do isite=1,nsite_new
                 temp=2d0*pi*dble(imom*isite)/dble(nsite_new)
                 phase=dcmplx(dcos(temp))+(0d0,1d0)*dcmplx(dsin(temp))
                 !be careful about the normalization!!
                 !phase=phase/dcmplx(dsqrt(dble(nsite_new)))
                 phase=phase/dcmplx(dsqrt(dble(nsite_old)))
                 xmat_new(imat,jmat,idim,isite)=&
                      xmat_new(imat,jmat,idim,isite)&
                      +xmat_mom_new(imat,jmat,idim,imom)*phase
              end do
           end do
        end do
     end do
!$omp end do
!$omp end parallel
  end do

  return

END subroutine Fourier_transform
