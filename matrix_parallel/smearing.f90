!******  1-step smearing  ******
SUBROUTINE smearing_one_step(xmat,xmat_smeared,s,myrank)

  implicit none
  include 'size_parallel.h'
  !***** input *****
  double precision s
  integer myrank
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  !***** output *****
  double complex xmat_smeared(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  !******************
  integer imat,jmat,idim,isite
  double complex p,q
  
  p=dcmplx(1d0-2d0*s)
  q=dcmplx(s)
  
  do isite=1,nsite_local
     do idim=1,ndim
        do jmat=1,nmat_block
           do imat=1,nmat_block
              xmat_smeared(imat,jmat,idim,isite)=p*xmat(imat,jmat,idim,isite)&
                   &+q*(xmat(imat,jmat,idim,isite+1)+xmat(imat,jmat,idim,isite-1))
           end do
        end do
     end do
  end do
  
  call Adjust_margin_xmat(xmat_smeared,myrank)
  
  return
  
END SUBROUTINE smearing_one_step
!********************************
!******  smearing of xmat  ******
!********************************
SUBROUTINE smearing_xmat(xmat,xmat_smeared,s,myrank,nsmear)

  implicit none
  include 'size_parallel.h'
  !***** input *****
  double precision s
  integer myrank
  integer nsmear
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  !***** output *****
  double complex xmat_smeared(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  !**************************
  integer ismear
  double complex temp(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)

  if(nsmear.GT.0)then
     temp=xmat
     
     do ismear=1,nsmear
        call smearing_one_step(temp,xmat_smeared,s,myrank)
        temp=xmat_smeared
     end do
  else if(nsmear.EQ.0)then
     xmat_smeared=xmat
  end if
  
  return
  
END SUBROUTINE smearing_xmat
!************************************
!******  smearing of dS_PF/dX  ******
!************************************
SUBROUTINE smearing_delPF(delPF_xmat_smeared,delPF_xmat,s,myrank,nsmear)

  implicit none
  include 'size_parallel.h'
  !***** input *****
  double precision s
  integer myrank
  integer nsmear
  !***** input *****
  double complex delPF_xmat_smeared(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)

  !***** output *****
  double complex delPF_xmat(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)

  !**************************
  double complex temp1(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  double complex temp2(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  integer ismear
  integer imat,jmat,idim,isite

  if(nsmear.GT.0)then
     temp1=(0d0,0d0)
     do isite=1,nsite_local
        do idim=1,ndim
           do jmat=1,nmat_block
              do imat=1,nmat_block
                 temp1(imat,jmat,idim,isite)=delPF_xmat_smeared(imat,jmat,idim,isite)
              end do
           end do
        end do
     end do
     call Adjust_margin_xmat(temp1,myrank)
     
     do ismear=1,nsmear
        call smearing_one_step(temp1,temp2,s,myrank)
        temp1=temp2
     end do
     do isite=1,nsite_local
        do idim=1,ndim
           do jmat=1,nmat_block
              do imat=1,nmat_block
                 delPF_xmat(imat,jmat,idim,isite)=temp2(imat,jmat,idim,isite)
              end do
           end do
        end do
     end do
     
  else if(nsmear.EQ.0)then
     delPF_xmat=delPF_xmat_smeared
  end if

     
  return
  
END SUBROUTINE smearing_delPF



