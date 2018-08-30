SUBROUTINE Adjust_margin_xmat(xmat)
  
  implicit none
  include '../staticparameters.f90'
  !***** input & output *****
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  !**************************
  integer imat,jmat
  integer idim
  integer isite
  
  do isite=-(nmargin-1),0
     do idim=1,ndim
        do imat=1,nmat
           do jmat=1,nmat
              xmat(imat,jmat,idim,isite)=xmat(imat,jmat,idim,nsite+isite)
           end do
        end do
     end do
  end do
  do isite=nsite+1,nsite+nmargin
     do idim=1,ndim
        do imat=1,nmat
           do jmat=1,nmat
              xmat(imat,jmat,idim,isite)=xmat(imat,jmat,idim,isite-nsite)
           end do
        end do
     end do
  end do
  
  return
  
END SUBROUTINE Adjust_margin_xmat
!***************************************************************
SUBROUTINE Adjust_margin_and_bc_pf(pf,nbc)
  
  implicit none
  include '../staticparameters.f90'
  !***** input *****
  integer nbc
  !***** input & output *****
  double complex pf(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
  !*******************
  integer isite,imat,jmat,ispin,ipf
  
  do ipf=1,npf
     if(nbc.EQ.0)then
        !pbc
        do imat=1,nmat     
           do jmat=1,nmat
              do ispin=1,nspin
                 do isite=-(nmargin-1),0
                    pf(imat,jmat,ispin,isite,ipf)=&
                         &pf(imat,jmat,ispin,nsite+isite,ipf)
                 end do
              end do
           enddo
        end do
        do imat=1,nmat     
           do jmat=1,nmat
              do ispin=1,nspin
                 do isite=nsite+1,nsite+nmargin
                    pf(imat,jmat,ispin,isite,ipf)=&
                         &pf(imat,jmat,ispin,isite-nsite,ipf)
                 end do
              end do
           end do
        end do
     else if(nbc.EQ.1)then!apbc
        do imat=1,nmat     
           do jmat=1,nmat
              do ispin=1,nspin
                 do isite=-(nmargin-1),0
                    pf(imat,jmat,ispin,isite,ipf)=&
                         &pf(imat,jmat,ispin,nsite+isite,ipf)*(-1d0,0d0)
                 end do
              end do
           enddo
        end do
        do imat=1,nmat     
           do jmat=1,nmat
              do ispin=1,nspin
                 do isite=nsite+1,nsite+nmargin
                    pf(imat,jmat,ispin,isite,ipf)=&
                         &pf(imat,jmat,ispin,isite-nsite,ipf)*(-1d0,0d0)
                 end do
              end do
           end do
        end do
     end if

  end do
     
  return

END SUBROUTINE Adjust_margin_and_bc_pf
 
!***************************************************************
SUBROUTINE Adjust_margin_and_bc_Chi(Chi,nbc,nremez)

  implicit none
  include '../staticparameters.f90'
  !***** input *****
  integer nbc,nremez
  !***** input & output *****
  double complex Chi(1:nremez,1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
  !**************************
  integer isite,imat,jmat,ispin,iremez,ipf

  do ipf=1,npf
     if(nbc.EQ.0)then!pbc
        do ispin=1,nspin
           do jmat=1,nmat
              do imat=1,nmat  
                 do iremez=1,nremez          
                    do isite=-(nmargin-1),0
                       Chi(iremez,imat,jmat,ispin,isite,ipf)=&
                            Chi(iremez,imat,jmat,ispin,isite+nsite,ipf)
                    end do
                 end do
              end do
           end do
        end do
        do ispin=1,nspin
           do jmat=1,nmat
              do imat=1,nmat  
                 do iremez=1,nremez
                    do isite=nsite+1,nsite+nmargin
                       Chi(iremez,imat,jmat,ispin,isite,ipf)=&
                            Chi(iremez,imat,jmat,ispin,isite-nsite,ipf)
                    end do
                 end do
              end do
           end do
        end do
     else if(nbc.EQ.1)then!apbc
        do ispin=1,nspin
           do jmat=1,nmat
              do imat=1,nmat  
                 do iremez=1,nremez          
                    do isite=-(nmargin-1),0
                       Chi(iremez,imat,jmat,ispin,isite,ipf)=&
                            Chi(iremez,imat,jmat,ispin,isite+nsite,ipf)*(-1d0,0d0)
                    end do
                 end do
              end do
           end do
        end do
        do ispin=1,nspin
           do jmat=1,nmat
              do imat=1,nmat  
                 do iremez=1,nremez
                    do isite=nsite+1,nsite+nmargin
                       Chi(iremez,imat,jmat,ispin,isite,ipf)=&
                            Chi(iremez,imat,jmat,ispin,isite-nsite,ipf)*(-1d0,0d0)
                    end do
                 end do
              end do
           end do
        end do
     end if
  end do
  
  return
  
END SUBROUTINE Adjust_margin_and_bc_Chi

!***************************************************************
SUBROUTINE Adjust_margin_and_bc_pf_2(pf,nbc)
  
  implicit none
  include '../staticparameters.f90'
  !***** input *****
  integer nbc
  !***** input & output *****
  double complex pf(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
  !*******************
  integer isite,imat,jmat,ispin
  
  if(nbc.EQ.0)then
     !pbc
     do imat=1,nmat     
        do jmat=1,nmat
           do ispin=1,nspin
              do isite=-(nmargin-1),0
                 pf(imat,jmat,ispin,isite)=&
                      &pf(imat,jmat,ispin,nsite+isite)
              end do
           end do
        enddo
     end do
     do imat=1,nmat     
        do jmat=1,nmat
           do ispin=1,nspin
              do isite=nsite+1,nsite+nmargin
                 pf(imat,jmat,ispin,isite)=&
                      &pf(imat,jmat,ispin,isite-nsite)
              end do
           end do
        end do
     end do
  else if(nbc.EQ.1)then!apbc
     do imat=1,nmat     
        do jmat=1,nmat
           do ispin=1,nspin
              do isite=-(nmargin-1),0
                 pf(imat,jmat,ispin,isite)=&
                      &pf(imat,jmat,ispin,nsite+isite)*(-1d0,0d0)
              end do
           end do
        enddo
     end do
     do imat=1,nmat     
        do jmat=1,nmat
           do ispin=1,nspin
              do isite=nsite+1,nsite+nmargin
                 pf(imat,jmat,ispin,isite)=&
                      &pf(imat,jmat,ispin,isite-nsite)*(-1d0,0d0)
              end do
           end do
        end do
     end do
  end if
      
  return

END SUBROUTINE Adjust_margin_and_bc_pf_2
