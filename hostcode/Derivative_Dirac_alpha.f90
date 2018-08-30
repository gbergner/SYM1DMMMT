!Derivative of (M*Chi)(t)_{ij} w.r.t. alpha_{k}  
SUBROUTINE Derivative_Dirac_alpha(alpha,&
     chi,Deriv_Mchi_alpha)
  
  implicit none
  include '../staticparameters.f90'
  
  double precision alpha(1:nmat)
  double complex chi(1:nremez_md,1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)
  double complex phase,phase2
  double precision aij
  double complex Deriv_Mchi_alpha(1:nremez_md,1:nmat,1:nmat,1:nmat,1:nspin,1:nsite,1:npf)

  integer imat,jmat,kmat
  integer ispin
  integer isite
  integer iremez
  integer ipf
  !***************************
  !***************************
  !***  kinetic part only  ***
  !***************************
  !***************************
  Deriv_Mchi_alpha=(0d0,0d0)
  do ipf=1,npf
     !********************
     !*** Naive Action ***
     !********************
     if(nimprove.EQ.0)then
        do isite=1,nsite
           do imat=1,nmat
              do jmat=1,nmat
                 
                 aij=alpha(imat)-alpha(jmat)
                 aij=aij/dble(nsite)
                 phase=dcmplx(dcos(aij))+(0d0,1d0)*dcmplx(dsin(aij))
                 
                 do iremez=1,nremez_md
                    do ispin=1,8
                       
                       kmat=jmat
                       Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite,ipf)=&
                            Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite,ipf)&
                         +chi(iremez,imat,jmat,ispin,isite+1,ipf)*phase&
                         /dcmplx(nsite)*(0d0,-1d0)
                       
                       kmat=imat
                       Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite,ipf)=&
                            Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite,ipf)&
                            -chi(iremez,imat,jmat,ispin,isite+1,ipf)*phase&
                         /dcmplx(nsite)*(0d0,-1d0)
                    
                    
                       kmat=jmat
                       Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite,ipf)=&
                            Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite,ipf)&
                            +chi(iremez,imat,jmat,ispin+8,isite-1,ipf)*dconjg(phase)&
                            /dcmplx(nsite)*(0d0,-1d0)
                       
                      kmat=imat
                       Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite,ipf)=&
                            Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite,ipf)&
                            -chi(iremez,imat,jmat,ispin+8,isite-1,ipf)*dconjg(phase)&
                            /dcmplx(nsite)*(0d0,-1d0)
                       
                    end do
                 end do
              end do
           end do
        end do
        !***********************
        !*** Improved Action ***
        !***********************
     else if(nimprove.EQ.1)then
        do isite=1,nsite
           do imat=1,nmat
              do jmat=1,nmat
                 
                 aij=alpha(imat)-alpha(jmat)
                 aij=aij/dble(nsite)
                 phase=dcmplx(dcos(aij))+(0d0,1d0)*dcmplx(dsin(aij))
                 phase2=phase*phase
                 do iremez=1,nremez_md
                    do ispin=1,8
                       
                       kmat=jmat
                       Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite,ipf)=&
                            Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite,ipf)&
                            +chi(iremez,imat,jmat,ispin,isite+1,ipf)*phase&
                            /dcmplx(nsite)*(0d0,-2d0)&
                            +chi(iremez,imat,jmat,ispin,isite+2,ipf)*phase2&
                            /dcmplx(nsite)*(0d0,1d0)  
                       
                       kmat=imat
                       Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite,ipf)=&
                            Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin+8,isite,ipf)&
                            -chi(iremez,imat,jmat,ispin,isite+1,ipf)*phase&
                            /dcmplx(nsite)*(0d0,-2d0)&
                            -chi(iremez,imat,jmat,ispin,isite+2,ipf)*phase2&
                            /dcmplx(nsite)*(0d0,1d0)
                       
                    
                       kmat=jmat
                       Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite,ipf)=&
                            Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite,ipf)&
                            +chi(iremez,imat,jmat,ispin+8,isite-1,ipf)*dconjg(phase)&
                            /dcmplx(nsite)*(0d0,-2d0)&
                            +chi(iremez,imat,jmat,ispin+8,isite-2,ipf)*dconjg(phase2)&
                            /dcmplx(nsite)*(0d0,1d0)
                       
                       kmat=imat
                       Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite,ipf)=&
                            Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite,ipf)&
                            -chi(iremez,imat,jmat,ispin+8,isite-1,ipf)*dconjg(phase)&
                            /dcmplx(nsite)*(0d0,-2d0)&
                            -chi(iremez,imat,jmat,ispin+8,isite-2,ipf)*dconjg(phase2)&
                            /dcmplx(nsite)*(0d0,1d0)
                       
                    end do
                 end do
              end do
           end do
        end do
     end if
  end do
  
  return
  
END SUBROUTINE Derivative_Dirac_alpha
