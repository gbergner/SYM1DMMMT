! Generate P_xmat with Gaussian weight.
! We do not have to care about the traceless condition. 
SUBROUTINE Generate_P_xmat(P_xmat,myrank)
  
  implicit none
  include 'size_parallel.h'
  !***** input *****
  integer myrank
  !***** output *****
  double complex P_xmat(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)
  !******************
  integer imat,jmat,idim,isite
  double precision r1,r2
  integer isublat,iblock,jblock,i,j,s
  
  call who_am_i(myrank,isublat,iblock,jblock)
  !**************************
  !*** off-diagonal block ***
  !**************************  
  do s=1,nsublat
     do i=1,nblock-1
        do j=i+1,nblock
           
           if ((s.EQ.isublat).AND.(i.EQ.iblock).AND.(j.EQ.jblock))then
              
              do isite=1,nsite_local
                 do idim=1,ndim
                    do imat=1,nmat_block
                       do jmat=1,nmat_block
                          call BoxMuller(r1,r2)
                          P_xmat(imat,jmat,idim,isite)=&
                               dcmplx(r1/dsqrt(2d0))+&
                               dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
                       end do
                    end do
                 end do
              end do
              
           else if((s.EQ.isublat).AND.(j.EQ.iblock).AND.(i.EQ.jblock))then
              do isite=1,nsite_local
                 do idim=1,ndim
                    do imat=1,nmat_block
                       do jmat=1,nmat_block
                          call BoxMuller(r1,r2)
                          P_xmat(jmat,imat,idim,isite)=&
                               dcmplx(r1/dsqrt(2d0))-&
                               dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
                       end do
                    end do
                 end do
              end do
              
           else
              !throw away some random numbers 
              !in order to avoid the communication
              do isite=1,nsite_local
                 do idim=1,ndim
                    do imat=1,nmat_block
                       do jmat=1,nmat_block
                          call BoxMuller(r1,r2)
                       end do
                    end do
                    
                 end do
              end do
           end if
           
        end do
     end do
  end do
  !**********************
  !*** diagonal block ***
  !**********************  
  do s=1,nsublat
     do i=1,nblock
        if((s.EQ.isublat).AND.(i.EQ.iblock).AND.(i.EQ.jblock))then
           do isite=1,nsite_local
              do idim=1,ndim
                 do imat=1,nmat_block-1
                    do jmat=imat+1,nmat_block
                       call BoxMuller(r1,r2)
                       P_xmat(imat,jmat,idim,isite)=&
                            dcmplx(r1/dsqrt(2d0))&
                            &+dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
                       P_xmat(jmat,imat,idim,isite)=&
                            dcmplx(r1/dsqrt(2d0))&
                            &-dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
                    end do
                 end do
                 do imat=1,nmat_block
                    call BoxMuller(r1,r2)
                    P_xmat(imat,imat,idim,isite)=dcmplx(r1)
                 end do
              end do
           end do
        else
           !throw away some random numbers 
           !in order to avoid the communication
           do isite=1,nsite_local
              do idim=1,ndim
                 do imat=1,nmat_block-1
                    do jmat=imat+1,nmat_block
                       call BoxMuller(r1,r2)
                    end do
                 end do
                 do imat=1,nmat_block
                    call BoxMuller(r1,r2)
                 end do
              end do
           end do
        end if
     end do
  end do

  return
  
END SUBROUTINE Generate_P_xmat
