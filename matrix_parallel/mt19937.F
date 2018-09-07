C     // initialization of random number generator =============================
      subroutine init_rand(iseed)

      implicit none

      integer, parameter :: N = 624

C     // rundom number seed
      integer :: iseed 

C     // random number sequence
      integer ir_here, irand(0:N-1)
      common /rng_mt19937/ ir_here, irand

C     // generate initial numbers by linear congruential method 
      irand(0)= iand(iseed,-1)
      do ir_here=1,N-1
         irand(ir_here) = iand(69069 * irand(ir_here-1), -1)
      enddo

      return
      end


C     // Mersenne Twister 19937 random number generator : double in [0,1) ======
C     // Ref: ACM Transactions on Modeling and Computer Simulation
C     //      Vol. 8, Issue 1, 1998, 3-30 

      real(8) function rand_1double()

      implicit none

C     // parameters and variables for MT19937
      integer, parameter :: N = 624, M = 397
      integer, parameter :: mask_up = -2147483648, mask_low =  2147483647
      integer, parameter :: Amat(0:1) = (/0, -1727483681/)
      integer :: ishift1,ishift2,ishift3,ishift4
      integer, parameter :: mask_b = -1658038656, mask_c = -272236544

C     // random number and its indices
      integer :: ir_here, irand(0:N-1), ir, irand_tmp
      common /rng_mt19937/ ir_here, irand


C     // generate N numbers at once, only when ir_here exceeds N 
      if(ir_here .ge. N) then

C        // update earlier N-M numbers
         do ir=0,N-M-1
C           // mix (w-r) and r bits
            irand_tmp = ior( iand(irand(ir),  mask_up),
     .                       iand(irand(ir+1),mask_low) )
C           // multiply A and ieor 
            irand(ir) = ieor( ieor( irand(ir+M),ishft(irand_tmp,-1) ),
     .                        Amat(iand(irand_tmp,1)) )
         enddo

C        // update next M-1 numbers
         do ir=N-M,N-2
            irand_tmp = ior( iand(irand(ir),  mask_up),
     .                       iand(irand(ir+1),mask_low) )
            irand(ir) = ieor( ieor( irand(ir+(M-N)),ishft(irand_tmp,-1) ),
     .                        Amat(iand(irand_tmp,1)))
         enddo

C        // update latest one
         irand_tmp  = ior( iand(irand(N-1),mask_up),
     .                     iand(irand(0),  mask_low))
         irand(N-1) = ieor( ieor( irand(M-1),ishft(irand_tmp,-1) ),
     .                      Amat(iand(irand_tmp,1)))
         ir_here = 0

      endif

C     // tempering
      irand_tmp = irand(ir_here)
      irand_tmp = ieor( irand_tmp, ishft(irand_tmp,-11) )
      irand_tmp = ieor( irand_tmp, iand( ishft(irand_tmp,  7), mask_b ) )
      irand_tmp = ieor( irand_tmp, iand( ishft(irand_tmp, 15), mask_c ) )
      irand_tmp = ieor( irand_tmp, ishft(irand_tmp,-18) )

      ir_here = ir_here+1

c       rand_1double = (dble(irand_tmp)+dble(2.0d0**31))/dble(2.0d0**32)

      if( irand_tmp .lt. 0 ) then
         rand_1double = dble(irand_tmp)/(2.0d0**32) + 1.e0;
      else
         rand_1double = dble(irand_tmp)/(2.0d0**32)
      endif

      return
      end


C     // Mersenne Twister 19937 random number generator : 32 bit integer =======
C     // Ref: ACM Transactions on Modeling and Computer Simulation
C     //      Vol. 8, Issue 1, 1998, 3-30 

      integer function rand_1int()

      implicit none

C     // parameters and variables for MT19937
      integer, parameter :: N = 624, M = 397
      integer, parameter :: mask_up = -2147483648, mask_low =  2147483647
      integer, parameter :: Amat(0:1) = (/0, -1727483681/)
      integer :: ishift1,ishift2,ishift3,ishift4
      integer, parameter :: mask_b = -1658038656, mask_c = -272236544

C     // random number and its indices
      integer :: ir_here, irand(0:N-1), ir, irand_tmp
      common /rng_mt19937/ ir_here, irand


C     // generate N numbers at once, only when ir_here exceeds N 
      if(ir_here .ge. N) then

C        // update earlier N-M numbers
         do ir=0,N-M-1
C           // mix (w-r) and r bits
            irand_tmp = ior( iand(irand(ir),  mask_up),
     .                       iand(irand(ir+1),mask_low) )
C           // multiply A and ieor 
            irand(ir) = ieor( ieor( irand(ir+M),ishft(irand_tmp,-1) ),
     .                        Amat(iand(irand_tmp,1)) )
         enddo

C        // update next M-1 numbers
         do ir=N-M,N-2
            irand_tmp = ior( iand(irand(ir),  mask_up),
     .                       iand(irand(ir+1),mask_low) )
            irand(ir) = ieor( ieor( irand(ir+(M-N)),ishft(irand_tmp,-1) ),
     .                        Amat(iand(irand_tmp,1)))
         enddo

C        // update latest one
         irand_tmp  = ior( iand(irand(N-1),mask_up),
     .                     iand(irand(0),  mask_low))
         irand(N-1) = ieor( ieor( irand(M-1),ishft(irand_tmp,-1) ),
     .                      Amat(iand(irand_tmp,1)))
         ir_here = 0

      endif

C     // tempering
      irand_tmp = irand(ir_here)
      irand_tmp = ieor( irand_tmp, ishft(irand_tmp,-11) )
      irand_tmp = ieor( irand_tmp, iand( ishft(irand_tmp,  7), mask_b ) )
      irand_tmp = ieor( irand_tmp, iand( ishft(irand_tmp, 15), mask_c ) )
      irand_tmp = ieor( irand_tmp, ishft(irand_tmp,-18) )

      ir_here = ir_here+1

      rand_1int = irand_tmp

      return
      end


C     // save random numbers ===================================================
      subroutine save_rand(ifile)

      implicit none

      integer, parameter :: N = 624

C     // random number sequence
      integer ir_here, irand(0:N-1)
      common /rng_mt19937/ ir_here, irand

      integer :: ifile

      write(ifile,*) ir_here
      write(ifile,*) irand

      return
      end


C     // load random numbers ===================================================
      subroutine load_rand(ifile)

      implicit none

      integer, parameter :: N = 624

C     // random number sequence
      integer ir_here, irand(0:N-1)
      common /rng_mt19937/ ir_here, irand

      integer :: ifile

      read(ifile,*) ir_here
      read(ifile,*) irand

      return
      end


C     // generate and store random numbers into an array =======================

      subroutine rand_int_array(irand, nr)

      implicit none

      integer :: irand, ir, nr, rand_1int
      dimension irand(nr)

C     // generate random numbers 
      do ir=1,nr
         irand(ir)=rand_1int()
      enddo

      return
      end


C     // generate and store random numbers into an array =======================

      subroutine rand_double_array(rand, nr)

      implicit none

      integer :: irand, ir, nr
      real(8) :: rand, rand_1double
      dimension rand(nr)

C     // generate random numbers 
      do ir=1,nr
         rand(ir)=rand_1double()
      enddo

      return
      end


