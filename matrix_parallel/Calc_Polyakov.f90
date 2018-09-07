!**************************************
!*** Polaakov loop (absolute value) ***
!************************************** 
SUBROUTINE Calc_Polyakov(nmat,alpha,Pol)
  
  implicit none
  
  integer nmat
  double precision alpha(1:nmat)
  double precision Pol,re_pol,im_pol
  integer imat
  
  re_pol=0d0
  im_pol=0d0
 
  do imat=1,nmat
     re_pol=re_pol+dcos(alpha(imat))
     im_pol=im_pol+dsin(alpha(imat))
  end do

  re_pol=re_pol/dble(nmat)
  im_pol=im_pol/dble(nmat)
 
  Pol=dsqrt(re_pol**2d0+im_pol**2d0)
  
  return
  
END SUBROUTINE Calc_Polyakov
