!##############################################################
!######    Estimation of auto-correlation                 #####
!######    written by Masanori Hanada                     #####     
!##############################################################
program JackKnife

  implicit none
  !----------------------------------------
  !     Iteration
  integer iteration
  parameter(iteration=5000)
  !---------------------------------------- 
  integer Binsize,BinStart,BinEnd,BinStep
  parameter(BinStart=5)
  parameter(BinEnd=101)
  parameter(BinStep=5)
  !----------------------------------------
  character(150) input,output,aho
  integer i,j

  doubleprecision readdata(1:7),average(1:7),err(1:7),BinAverage(1:7)
 ! doubleprecision readdata(1:5),average(1:5),err(1:5),BinAverage(1:5)
  integer read1
  !*************************
  !*** evaluation of err *** 
  !*************************
  write(input,"('N4S16T20_apbc.txt')")
  BinSize=BinStart
  do while(Binsize.LT.BinEnd)
     
     open(UNIT = 12, File = input, STATUS = "OLD", ACTION = "READ")
     do i=1,9
        read(12,*)aho
     end do

     average=0d0
     err=0d0
     do i=1,INT(dble(iteration)/dble(BinSize))
        BinAverage=0d0
        do j=1,Binsize
           read(12,*)read1,aho,readdata
           BinAverage=BinAverage+readdata
        end do
        BinAverage=BinAverage/dble(Binsize)
        average=average+BinAverage
        err=err+BinAverage**2d0

     end do
     close(12)
     average=average/dble(INT(dble(iteration)/dble(BinSize)))
     err=err/dble(INT(dble(iteration)/dble(BinSize)))
     err=err-average**2d0
     err=sqrt(err/dble(INT(dble(iteration)/dble(BinSize))-1))
     write(*,*)dble(binsize),&
          average(3),err(3),average(4),err(4),average(5),err(5),average(6),err(6),average(7),err(7)
     
     Binsize=Binsize+BinStep
  end do
  




 


end program JackKnife
