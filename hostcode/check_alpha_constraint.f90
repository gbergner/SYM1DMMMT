! impose the constraint max(alpha_i-alpha_j) < 2*pi. 
! info=0 -> OK, info=1 -> constraint is violated. 
SUBROUTINE check_alpha_constraint(alpha,info)
  
  implicit none
  include '../staticparameters.f90'
  !***** input *****
  double precision alpha(1:nmat)
  !***** output *****
  integer info
  !******************
  integer imat
  double precision max,min,pi
 
  
  max=alpha(1)
  min=alpha(1)
  do imat=2,nmat
     if(max.LT.alpha(imat))then
        max=alpha(imat)
     else if(min.GT.alpha(imat))then
        min=alpha(imat)
     end if
  end do
  pi=2d0*dasin(1d0)
  if(max-min.LT.2d0*pi)then
     info=0
  else
     info=1
  end if
  
  return
  
END SUBROUTINE check_alpha_constraint
