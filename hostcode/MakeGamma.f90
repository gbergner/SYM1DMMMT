!****************************************************
!*** 10d Gamma matrices, written by J. Nishimura. ***
!****************************************************
!sigma_{1,2,3}: Pauli matrices
!sigma_4 = identity
!
!gamma_1=sigma_3 X sigma_4 X sigma_4 X sigma_4
!gamma_2=sigma_2 X sigma_2 X sigma_2 X sigma_2
!gamma_3=sigma_2 X sigma_2 X sigma_4 X sigma_1
!gamma_4=sigma_2 X sigma_2 X sigma_4 X sigma_3
!gamma_5=sigma_2 X sigma_1 X sigma_2 X sigma_4
!gamma_6=sigma_2 X sigma_3 X sigma_2 X sigma_4
!gamma_7=sigma_2 X sigma_4 X sigma_1 X sigma_2
!gamma_8=sigma_2 X sigma_4 X sigma_3 X sigma_2
!gamma_9=sigma_2 X sigma_4 X sigma_4 X sigma_4*(0d0,-1d0)
!
!Here "X" is a tensor product.
SUBROUTINE MakeGamma(Gamma10d)

  implicit none
  include '../staticparameters.f90'

  double complex Gamma10d(1:ndim,1:nspin,1:nspin)
  double complex sigma(1:2,1:2,1:4)
  INTEGER I1,J1,I2,J2,I3,J3,I4,J4
  INTEGER I,J

  !Pauli matrices

  SIGMA(1,1,1)=(0d0,0d0)
  SIGMA(2,1,1)=(1d0,0d0)
  SIGMA(1,2,1)=(1d0,0d0)
  SIGMA(2,2,1)=(0d0,0d0)
  
  SIGMA(1,1,2)=(0d0,0d0)
  SIGMA(2,1,2)=(0d0,-1d0)
  SIGMA(1,2,2)=(0d0,1d0)
  SIGMA(2,2,2)=(0d0,0d0)
  
  SIGMA(1,1,3)=(1d0,0d0)
  SIGMA(2,1,3)=(0d0,0d0)
  SIGMA(1,2,3)=(0d0,0d0)
  SIGMA(2,2,3)=(-1d0,0d0)
  
  SIGMA(1,1,4)=(1d0,0d0)
  SIGMA(2,1,4)=(0d0,0d0)
  SIGMA(1,2,4)=(0d0,0d0)
  SIGMA(2,2,4)=(1d0,0d0)
  
  do i1 = 1,2
     do i2 = 1,2
        do i3 = 1,2
           do i4 = 1,2
              I=8*(i1-1)+4*(i2-1)+2*(i3-1)+i4
              do j1 = 1,2
                 do j2 = 1,2
                    do j3 = 1,2
                       do j4 = 1,2
                          J=8*(j1-1)+4*(j2-1)+2*(j3-1)+j4
                          
                          Gamma10d(1,I,J)=SIGMA(I1,J1,3)*SIGMA(I2,J2,4)&
                               *SIGMA(I3,J3,4)*SIGMA(I4,J4,4)
                          Gamma10d(2,I,J)=SIGMA(I1,J1,2)*SIGMA(I2,J2,2)&
                               *SIGMA(I3,J3,2)*SIGMA(I4,J4,2)
                          Gamma10d(3,I,J)=SIGMA(I1,J1,2)*SIGMA(I2,J2,2)&
                               *SIGMA(I3,J3,4)*SIGMA(I4,J4,1)
                          Gamma10d(4,I,J)=SIGMA(I1,J1,2)*SIGMA(I2,J2,2)&
                               *SIGMA(I3,J3,4)*SIGMA(I4,J4,3)
                          Gamma10d(5,I,J)=SIGMA(I1,J1,2)*SIGMA(I2,J2,1)&
                               *SIGMA(I3,J3,2)*SIGMA(I4,J4,4)
                          Gamma10d(6,I,J)=SIGMA(I1,J1,2)*SIGMA(I2,J2,3)&
                               *SIGMA(I3,J3,2)*SIGMA(I4,J4,4)
                          Gamma10d(7,I,J)=SIGMA(I1,J1,2)*SIGMA(I2,J2,4)&
                               *SIGMA(I3,J3,1)*SIGMA(I4,J4,2)
                          Gamma10d(8,I,J)=SIGMA(I1,J1,2)*SIGMA(I2,J2,4)&
                               *SIGMA(I3,J3,3)*SIGMA(I4,J4,2)
                          !Gamma10d(9,I,J)=SIGMA(I1,J1,1)*SIGMA(I2,J2,4)&
                          !     *SIGMA(I3,J3,4)*SIGMA(I4,J4,4)
                          Gamma10d(9,I,J)=(0d0,-1d0)*SIGMA(I1,J1,4)*SIGMA(I2,J2,4)&
                               *SIGMA(I3,J3,4)*SIGMA(I4,J4,4)
                          
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do  

  !******* A factor (-i) is needed for a historical reason.... *****
  !******* It is a very bad notation. Sorry. ***********************
  Gamma10d=Gamma10d*(0d0,-1d0)

  RETURN

END SUBROUTINE MakeGamma
