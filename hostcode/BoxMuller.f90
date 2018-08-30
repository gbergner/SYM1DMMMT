!****************************************************************
!****************************************************************
!***  Box-Muller method for generating Gaussian random number ***
!****************************************************************
!****************************************************************
SUBROUTINE BoxMuller(p,q)  

  use mtmod !Mersenne twistor
  implicit none 
  !***** output *****
  doubleprecision p,q
  !******************
  double precision r,s,Pi

  Pi=2d0*DASIN(1d0)
  !uniform random numbers between 0 and 1
  r=grnd()
  s=grnd()
  !Gaussian random numbers, 
  !with weights proportional to e^{-p^2/2} and e^{-q^2/2}
  p=dsqrt(-2d0*dlog(r))*DSIN(2d0*Pi*s)
  q=dsqrt(-2d0*dlog(r))*DCOS(2d0*Pi*s)

  return

END SUBROUTINE BoxMuller
