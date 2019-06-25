  !*************************************
  !*** Save the final configuration ****
  !*************************************


SUBROUTINE Read_Intermediate_Config(xmat,alpha,itraj,status)
  
  use mtmod !Mersenne twistor
  implicit none
  
  include '../staticparameters.f90'
  include '../unit_number.inc'
  !---------------------------------
  double precision alpha(1:nmat)
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  integer itraj,status
  if(status.EQ.0) then
   call mtgetus(unit_intermediate_config,status)
  end if
  if(status.EQ.0) then
   read(unit_intermediate_config,*,IOSTAT=status) itraj
  end if
  if(status.EQ.0) then
   read(unit_intermediate_config,*,IOSTAT=status) xmat
  end if
  if(status.EQ.0) then
   read(unit_intermediate_config,*,IOSTAT=status) alpha
  end if
  return
END SUBROUTINE Read_Intermediate_Config
