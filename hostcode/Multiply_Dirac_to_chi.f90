!Mchi=M*chi, M: Dirac op
SUBROUTINE Multiply_Dirac_to_chi(temperature,xmat,alpha,&
chi,mchi,GAMMA10d,nbmn,flux)

  implicit none
  include '../staticparameters.f90'
  !***** input *****
  integer nprocs,nbmn
  double precision temperature,flux
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double precision alpha(1:nmat)
  double complex chi(1:nremez_md,1:nmat,1:nmat,1:nspin,&
       &-(nmargin-1):nsite+nmargin,1:npf)
  double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
  !***** output *****
  double complex Mchi(1:nremez_md,1:nmat,1:nmat,1:nspin,&
       &-(nmargin-1):nsite+nmargin,1:npf)
  !******************
  double complex pf1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
  double complex pf2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
  integer imat,jmat
  integer ispin
  integer isite
  integer iremez
  integer ipf

  do ipf=1,npf
     do iremez=1,nremez_md
        
        do imat=1,nmat
           do jmat=1,nmat
              do ispin=1,nspin
                 do isite=-(nmargin-1),nsite+nmargin
                    pf1(imat,jmat,ispin,isite)=Chi(iremez,imat,jmat,ispin,isite,ipf)
                 end do
              end do
           end do
        end do
        call Multiply_Dirac(temperature,xmat,alpha,&
             pf1,pf2,GAMMA10d,nbmn,flux)
        do imat=1,nmat
           do jmat=1,nmat
              do ispin=1,nspin
                 do isite=-(nmargin-1),nsite+nmargin
                    MChi(iremez,imat,jmat,ispin,isite,ipf)=pf2(imat,jmat,ispin,isite)
                 end do
              end do
           end do
        end do
        
     end do
  end do
  
  return

END SUBROUTINE Multiply_Dirac_to_chi
