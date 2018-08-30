!Derivative of (M*pf) w.r.t. X_{ji}(t)  (not by X_{ij}(t))
SUBROUTINE Derivative_Dirac_X(iremez,ipf,temperature,&
     chi,Deriv_Mchi_X,GAMMA10d)
  
  implicit none
  include '../staticparameters.f90'
  !***** input *****
  integer iremez
  double precision temperature
  double complex chi(1:nremez_md,1:nmat,1:nmat,1:nspin,&
       &-(nmargin-1):nsite+nmargin,1:npf)
  double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
  !***** output *****
  !number of indices must be less than 8...
  double complex Deriv_MChi_X(1:nmat,1:nmat,1:ndim,&
       1:nmat,1:nmat,1:nspin,1:nsite)
  !******************
  double precision lattice_spacing
  integer imat,jmat,kmat,lmat
  integer idim
  integer ispin,jspin
  integer isite
  integer ipf

  lattice_spacing=1d0/temperature/dble(nsite)
  
  Deriv_MChi_X=(0d0,0d0)

  do isite=1,nsite
     do imat=1,nmat
        do jmat=1,nmat
           do ispin=1,nspin
              do idim=1,ndim
                 do jspin=1,nspin
                    lmat=jmat
                    do kmat=1,nmat
                       
                       Deriv_Mchi_X(lmat,kmat,idim,imat,jmat,ispin,isite)=&
                            Deriv_Mchi_X(lmat,kmat,idim,imat,jmat,ispin,isite)&
                            +GAMMA10d(idim,ispin,jspin)&
                            *chi(iremez,imat,kmat,jspin,isite,ipf)&
                            *dcmplx(lattice_spacing)
                       
                    end do
                    
                    kmat=imat
                    do lmat=1,nmat
                       
                       Deriv_Mchi_X(lmat,kmat,idim,imat,jmat,ispin,isite)=&
                            Deriv_Mchi_X(lmat,kmat,idim,imat,jmat,ispin,isite)&
                            -GAMMA10d(idim,ispin,jspin)&
                            *chi(iremez,lmat,jmat,jspin,isite,ipf)&
                            *dcmplx(lattice_spacing)
                       
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do
  
  return
  
END SUBROUTINE Derivative_Dirac_X
