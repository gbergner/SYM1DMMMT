  !*************************************
  !*** Save the final configuration ****
  !*************************************


SUBROUTINE Save_Intermediate_Config(xmat,alpha,itraj)
  
  use mtmod !Mersenne twistor
  implicit none
  
  include '../staticparameters.f90'
  include '../unit_number.inc'
  !---------------------------------
  double precision alpha(1:nmat)
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  integer itraj


  call mtsaveu(unit_intermediate_config)  
  write(unit_intermediate_config,*) itraj
  write(unit_intermediate_config,*) xmat
  write(unit_intermediate_config,*) alpha
  
  return

END SUBROUTINE Save_Intermediate_Config
