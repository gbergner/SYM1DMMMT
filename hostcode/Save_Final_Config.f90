  !*************************************
  !*** Save the final configuration ****
  !*************************************


SUBROUTINE Save_Final_Config(xmat,alpha,itraj,output_config)
  
  use mtmod !Mersenne twistor
  implicit none
  
  include '../staticparameters.f90'
  include '../unit_number.inc'
  !---------------------------------
  double precision alpha(1:nmat)
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  integer itraj
  character(1000) output_config


  open(UNIT=unit_output_config, File = output_config, STATUS = "REPLACE", ACTION = "WRITE")
  call mtsaveu(unit_output_config)  
  write(unit_output_config,*) itraj
  write(unit_output_config,*) xmat
  write(unit_output_config,*) alpha
  close(unit_output_config)

  return

END SUBROUTINE Save_Final_Config
