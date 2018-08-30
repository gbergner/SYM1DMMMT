  !*************************************
  !*** Save the final configuration ****
  !*************************************


SUBROUTINE Save_Acceleration_Parameters(fluctuation,imeasure,acc_output)
  
  use mtmod !Mersenne twistor
  implicit none
  
  include '../staticparameters.f90'
  include '../unit_number.inc'
  !---------------------------------
  double complex xmat_mom(1:nmat,1:nmat,1:ndim,1:nsite) 
  double precision fluctuation(1:nsite) 
  integer imeasure
  character(1000) acc_output


  call Fourier_acceleration_optimize(xmat_mom,fluctuation,&
       &imeasure,2)



  open(UNIT=unit_output_acc, File = acc_output, STATUS = "REPLACE", ACTION = "WRITE")
  write(unit_output_acc,*) fluctuation
  close(unit_output_acc)  

  return

END SUBROUTINE Save_Acceleration_Parameters
