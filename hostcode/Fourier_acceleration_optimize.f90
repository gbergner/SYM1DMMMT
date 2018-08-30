! info=1 -> measure fluctuation
! info=2 -> set the acceleration parameter
subroutine Fourier_acceleration_optimize(xmat_mom,fluctuation,imeasure,info)

  implicit none
  
  include '../staticparameters.f90'
  !***** input *****
  integer info 
  !***** input, when info=1 *****
  double complex xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite)
  !***** input & output *****
  double precision fluctuation(1:nsite)
  !***** input & output when info=1, input when info=2 *****
  integer imeasure !number of measurement performed.
  !************************
  double precision temp
  integer imom,imat,jmat,idim
  
  
  if(info.EQ.1)then
     !******************************************************************
     !*** measure fluctuation and add to the sum of previous values. ***
     !******************************************************************
     !fluctuation(p)=|Xmat_mom(p)|
     do imom=1,nsite
        do imat=1,nmat
           do jmat=1,nmat
              do idim=1,ndim
                 fluctuation(imom)=fluctuation(imom)&
                      +dble(xmat_mom(imat,jmat,idim,imom)&
                      *dconjg(xmat_mom(imat,jmat,idim,imom)))
              end do
           end do
        end do
     end do
     imeasure=imeasure+1

  else if(info.EQ.2)then

     fluctuation=fluctuation/dble(nmat*imeasure)
     !*************************************************************
     !*** fluctuation(p) and fluctuation(-p) should be averaged ***
     !*************************************************************
     do imom=1,int(dble(nsite)*0.5d0+0.01d0)-1
        temp=fluctuation(imom)+fluctuation(nsite-imom)
        temp=temp*0.5d0
        fluctuation(imom)=temp
        fluctuation(nsite-imom)=temp
     end do
     !**********************************************************
     !*** Fourier acceleration parameter = sqrt(fluctuation) ***
     !**********************************************************
     fluctuation=dsqrt(fluctuation)
     
        
  end if

  return

END subroutine Fourier_acceleration_optimize
