! Generate P_alpha with Gaussian weight.
! We do not have to care about the traceless condition. 
SUBROUTINE Generate_P_alpha(nmat,P_alpha)
  
  !use mtmod !Mersenne twistor
  implicit none
  integer nmat
  integer imat
  double precision r1,r2

  double precision P_alpha(1:nmat)
  
  do imat=1,nmat
     call BoxMuller(r1,r2)
     P_alpha(imat)=r1
  end do
  
  !write(*,*)r1

  return
  
END SUBROUTINE Generate_P_alpha
