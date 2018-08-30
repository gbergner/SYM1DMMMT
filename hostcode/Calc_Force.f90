!****************************************************************
!***** Calculate the force term in the molecular evolution. *****
!****************************************************************
! delh_alpha(imat) = dH/(d alpha_imat)
! delh_xmat(imat,jmat)=dH/(d xmat(jmat,imat)), 
!      be careful about the difference of ordering of indices.
! ******* Note that "Force" = (-1)*delh_xmat, (-1)*delh_alpha: ******
! alpha -> alpha + P_alpha*dtau_alpha
! xmat -> xmat + P_xmat*dtau_xmat
! P_alpha -> P_alpha - delh_alpha*dtau_alpha
! P_xmat -> P_xmat - delh_xmat*dtau_xmat

SUBROUTINE Calc_Force_bosonic(delh_xmat,delh_alpha,xmat,alpha,chi,&
     GAMMA10d,gcoeff_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md)

  implicit none

  include '../staticparameters.f90'
  !****** input *****
  integer nbmn
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double precision alpha(1:nmat)
  double complex Chi(1:nremez_md,1:nmat,1:nmat,1:nspin,&
       &-(nmargin-1):nsite+nmargin,1:npf)
  double precision gcoeff_alpha
  double precision g_R,trx2,RCUT
  double precision temperature,flux
  double precision acoeff_md(0:nremez_md)
  double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
  !***** output *****
  double complex delh_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double precision delh_alpha(1:nmat)
  !*******************
  double precision alpha_max,alpha_min
  double precision lattice_spacing,pi
  integer imat,jmat,kmat,lmat,imax,imin
  integer idim,jdim
  integer isite,isite_p1,isite_m1
  double complex commutator(1:nmat,1:nmat,1:ndim,1:ndim),uxdu(1:nmat,1:nmat),&
       &udxu(1:nmat,1:nmat),xxx,ei,ej
  integer iremez
  integer ispin
  integer ipf
  double complex com12(1:nmat,1:nmat),com23(1:nmat,1:nmat),&
       &com31(1:nmat,1:nmat),temp,uxumx(1:nmat,1:nmat)



  pi=2d0*dasin(1d0)
  lattice_spacing=1d0/temperature/dble(nsite)
  !*****************************
  !*** calculate delh_alpha ****
  !*****************************
  delh_alpha=0d0
  if(ngauge.EQ.1)then
     ! delh_xmat=(0d0,0d0)
     if(nimprove.EQ.0)then
        do isite=1,nsite
           isite_p1=isite+1
           
           do idim=1,ndim
              do imat=1,nmat
                 do jmat=1,nmat
                    !exp(i*alpha_i)
                    ei=dcmplx(dcos(alpha(imat)/dble(nsite)))&
                         +(0d0,1d0)*dcmplx(dsin(alpha(imat)/dble(nsite)))
                    !exp(-i*alpha_j)
                    ej=dcmplx(dcos(alpha(jmat)/dble(nsite)))&
                         -(0d0,1d0)*dcmplx(dsin(alpha(jmat)/dble(nsite)))
                    
                    delh_alpha(imat)=delh_alpha(imat)&
                         -dble(ei*xmat(imat,jmat,idim,isite_p1)&
                         *ej*xmat(jmat,imat,idim,isite)*(0d0,1d0))
                    delh_alpha(jmat)=delh_alpha(jmat)&
                         +dble(ei*xmat(imat,jmat,idim,isite_p1)&                      
                         *ej*xmat(jmat,imat,idim,isite)*(0d0,1d0))
                    
                 end do
              end do
           end do
        end do
        
     else if(nimprove.EQ.1)then
        do isite=1,nsite
           do idim=1,ndim
              !-0.5*u^2*x(t+2a)*(u^dagger)^2 + 2*u*x(t+a)*u^dagger - 1.5*x(t)
              do imat=1,nmat
                 do jmat=1,nmat
                    
                    !exp(i*alpha_i)
                    ei=dcmplx(dcos(alpha(imat)&
                         &/dble(nsite)))&
                         &+(0d0,1d0)*dcmplx(dsin(alpha(imat)&
                         &/dble(nsite)))
                    !exp(-i*alpha_j)
                    ej=dcmplx(dcos(alpha(jmat)&
                         &/dble(nsite)))&
                         -(0d0,1d0)*dcmplx(dsin(alpha(jmat)&
                         &/dble(nsite)))
                    uxumx(imat,jmat)=&
                         -(0.5d0,0d0)*ei*ei*xmat(imat,jmat,idim,isite+2)*ej*ej&
                         +(2d0,0d0)*ei*xmat(imat,jmat,idim,isite+1)*ej&
                         -(1.5d0,0d0)*xmat(imat,jmat,idim,isite)
                    
                    delh_alpha(imat)&
                         &=delh_alpha(imat)&
                         &+dble(((2d0,0d0)*ei*xmat(imat,jmat,idim,isite+1)*ej&
                         &-ei*ei*xmat(imat,jmat,idim,isite+2)*ej*ej)&
                         *dconjg(uxumx(imat,jmat))*(0d0,1d0))
                    
                    delh_alpha(jmat)&
                         &=delh_alpha(jmat)&
                         &-dble(((2d0,0d0)*ei*xmat(imat,jmat,idim,isite+1)*ej&
                         &-ei*ei*xmat(imat,jmat,idim,isite+2)*ej*ej)&
                         *dconjg(uxumx(imat,jmat))*(0d0,1d0))
                    
                 end do
              end do
           end do
        end do
        
     end if
     
     delh_alpha=delh_alpha*dble(nmat)/lattice_spacing/dble(nsite)

     !gauge-fixing term
     do imat=1,nmat
        do jmat=1,nmat
           if(jmat.NE.imat)then
              delh_alpha(imat)=delh_alpha(imat)& 
                   -1d0/dtan(0.5d0*(alpha(imat)-alpha(jmat)))
           end if
        end do
     end do

  end if

  !***************************
  !*** calculate delh_xmat ***
  !***************************
  delh_xmat=(0d0,0d0)
  if(nimprove.EQ.0)then
     do isite=1,nsite
        isite_p1=isite+1
        isite_m1=isite-1
        
        do idim=1,ndim
           !u(t)*x(t+a)*u(t)^¥dagger
           !u(t-a)^¥dagger*x(t-a)*u(t-a)
           do imat=1,nmat
              do jmat=1,nmat
                 !exp(i*alpha_i)
                 ei=dcmplx(dcos(alpha(imat)/dble(nsite)))&
                      +(0d0,1d0)*dcmplx(dsin(alpha(imat)/dble(nsite)))
                 !exp(-i*alpha_j)
                 ej=dcmplx(dcos(alpha(jmat)/dble(nsite)))&
                      -(0d0,1d0)*dcmplx(dsin(alpha(jmat)/dble(nsite))) 
                 uxdu(imat,jmat)=&
                      ei*xmat(imat,jmat,idim,isite_p1)*ej
                 udxu(jmat,imat)=&
                      ej*xmat(jmat,imat,idim,isite_m1)*ei
              end do
           end do
           
           do imat=1,nmat
              do jmat=1,nmat
                 delh_xmat(imat,jmat,idim,isite)=&
                      dcmplx(2d0)*xmat(imat,jmat,idim,isite)
                 
                 delh_xmat(imat,jmat,idim,isite)=&
                      delh_xmat(imat,jmat,idim,isite)&
                      -uxdu(imat,jmat)&
                      -udxu(imat,jmat)
              end do
           end do
        end do
     end do
  else if(nimprove.EQ.1)then

     do isite=1,nsite
        do idim=1,ndim           
           do imat=1,nmat
              do jmat=1,nmat
                 !exp(i*alpha_i)
                 ei=dcmplx(dcos(alpha(imat)&
                      &/dble(nsite)))&
                      &+(0d0,1d0)*dcmplx(dsin(alpha(imat)&
                      &/dble(nsite)))
                 !exp(-i*alpha_j)
                 ej=dcmplx(dcos(alpha(jmat)&
                      &/dble(nsite)))&
                      &-(0d0,1d0)*dcmplx(dsin(alpha(jmat)&
                      &/dble(nsite)))
                 
                 delh_xmat(imat,jmat,idim,isite)=&
                      &delh_xmat(imat,jmat,idim,isite)&
                      &+(6.5d0,0d0)*xmat(imat,jmat,idim,isite)&
                      &-(4d0,0d0)*ei*xmat(imat,jmat,idim,isite+1)*ej&
                      &-(4d0,0d0)*dconjg(ei)*xmat(imat,jmat,idim,isite-1)&
                      &*dconjg(ej)&
                      &+(0.75d0,0d0)*ei*ei*xmat(imat,jmat,idim,isite+2)*ej*ej&
                      &+(0.75d0,0d0)*dconjg(ei*ei)&
                      &*xmat(imat,jmat,idim,isite-2)*dconjg(ej*ej)
                 
              end do
           end do
        end do
     end do
 


  end if

  delh_xmat=delh_xmat*dcmplx(nmat)/dcmplx(lattice_spacing)


  !commutator term
  do isite=1,nsite
     commutator=(0d0,0d0)
     do idim=1,ndim-1
        do jdim=idim+1,ndim
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    commutator(imat,jmat,idim,jdim)=&
                         commutator(imat,jmat,idim,jdim)&
                         +xmat(imat,kmat,idim,isite)*xmat(kmat,jmat,jdim,isite)&
                         -xmat(imat,kmat,jdim,isite)*xmat(kmat,jmat,idim,isite)
                 end do
              end do
           end do
        end do
     end do
     
     do idim=1,ndim
        do jdim=1,ndim
           do imat=1,nmat
              do jmat=1,nmat
                 xxx=(0d0,0d0)
                 if(jdim.GT.idim)then
                    do kmat=1,nmat
                       xxx=xxx&
                            +xmat(imat,kmat,jdim,isite)&
                            *commutator(kmat,jmat,idim,jdim)&
                            -commutator(imat,kmat,idim,jdim)&
                            *xmat(kmat,jmat,jdim,isite)
                    end do
                 else if(jdim.LT.idim)then
                    do kmat=1,nmat
                       xxx=xxx&
                            -xmat(imat,kmat,jdim,isite)&
                            *commutator(kmat,jmat,jdim,idim)&
                            +commutator(imat,kmat,jdim,idim)&
                            *xmat(kmat,jmat,jdim,isite)
                    end do
                 end if
                 
                 delh_xmat(imat,jmat,idim,isite)=&
                      delh_xmat(imat,jmat,idim,isite)&
                      -dcmplx(nmat)*dcmplx(lattice_spacing)*xxx
              end do
           end do
           
        end do
     end do
     
  end do 

  !****************************  
  !*** constraint for alpha ***
  !****************************
  if(ngauge.EQ.1)then
     imax=1
     imin=1
     alpha_max=alpha(1)
     alpha_min=alpha(1)
     do imat=2,nmat
        if(alpha(imat).GT.alpha_max)then
           imax=imat
           alpha_max=alpha(imat)
        else if(alpha(imat).LT.alpha_min)then
           imin=imat
           alpha_min=alpha(imat)
        end if
     end do
     if(alpha_max-alpha_min.LT.2d0*pi)then
        delh_alpha(imax)=delh_alpha(imax)&
             +1d0/(2d0*pi-(alpha_max-alpha_min)+1d0/gcoeff_alpha)
        delh_alpha(imin)=delh_alpha(imin)&
             -1d0/(2d0*pi-(alpha_max-alpha_min)+1d0/gcoeff_alpha)
     else if(alpha_max-alpha_min.GE.2d0*pi)then
        delh_alpha(imax)=delh_alpha(imax)+gcoeff_alpha
        delh_alpha(imin)=delh_alpha(imin)-gcoeff_alpha
     end if
     
   end if
  !****************************  
  !*** constraint for TrX^2 ***
  !****************************
  call Calc_TrX2(xmat,trx2)
  if(trx2.GE.RCUT)then
     do isite=1,nsite
        do idim=1,ndim
           do imat=1,nmat
              do jmat=1,nmat
                 delh_xmat(imat,jmat,idim,isite)=&
                      delh_xmat(imat,jmat,idim,isite)&
                      +2d0*g_R*xmat(imat,jmat,idim,isite)&
                      /dble(nsite)
              end do
           end do
        end do
     end do
  end if
  !******************************  
  !*** Plane wave deformation ***
  !******************************
  if(nbmn.EQ.1)then
     !***************************
     !*** deriv. of mass term ***
     !***************************
     do isite=1,nsite
        do imat=1,nmat
           do jmat=1,nmat
              temp=dcmplx(flux*flux*lattice_spacing*dble(nmat))
              do idim=1,3
                 delh_xmat(imat,jmat,idim,isite)=&
                      delh_xmat(imat,jmat,idim,isite)&
                      +xmat(imat,jmat,idim,isite)*temp
              end do
              temp=dcmplx(0.25d0*flux*flux*lattice_spacing*dble(nmat))
              do idim=4,9
                 delh_xmat(imat,jmat,idim,isite)=&
                      delh_xmat(imat,jmat,idim,isite)&
                      +xmat(imat,jmat,idim,isite)*temp
              end do
           end do
        end do
     end do
     !******************
     !*** cubic term ***
     !******************
     do isite=1,nsite
        com12=(0d0,0d0)
        com23=(0d0,0d0)
        com31=(0d0,0d0)
        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 com12(imat,jmat)=com12(imat,jmat)&
                      +xmat(imat,kmat,1,isite)*xmat(kmat,jmat,2,isite)&
                      -xmat(imat,kmat,2,isite)*xmat(kmat,jmat,1,isite)
                 com23(imat,jmat)=com23(imat,jmat)&
                      +xmat(imat,kmat,2,isite)*xmat(kmat,jmat,3,isite)&
                      -xmat(imat,kmat,3,isite)*xmat(kmat,jmat,2,isite)
                 com31(imat,jmat)=com31(imat,jmat)&
                      +xmat(imat,kmat,3,isite)*xmat(kmat,jmat,1,isite)&
                      -xmat(imat,kmat,1,isite)*xmat(kmat,jmat,3,isite)
              end do
           end do
        end do
        temp=(0d0,3d0)*dcmplx(flux*lattice_spacing*dble(nmat))
        do imat=1,nmat
           do jmat=1,nmat  
              delh_xmat(imat,jmat,1,isite)=&
                   delh_xmat(imat,jmat,1,isite)+com23(imat,jmat)*temp
              delh_xmat(imat,jmat,2,isite)=&
                   delh_xmat(imat,jmat,2,isite)+com31(imat,jmat)*temp
              delh_xmat(imat,jmat,3,isite)=&
                   delh_xmat(imat,jmat,3,isite)+com12(imat,jmat)*temp
           end do
        end do
     end do
  
  end if

  return
  
END SUBROUTINE Calc_Force_Bosonic
!**********************************************************
SUBROUTINE Calc_Force_fermionic(delh_xmat,delh_alpha,xmat,alpha,chi,&
     GAMMA10d,gcoeff_alpha,g_R,RCUT,nbmn,flux,temperature,acoeff_md)

  implicit none

  include '../staticparameters.f90'
  !****** input *****
  integer nbmn
  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double precision alpha(1:nmat)
  double complex Chi(1:nremez_md,1:nmat,1:nmat,1:nspin,&
       &-(nmargin-1):nsite+nmargin,1:npf)
  double precision gcoeff_alpha
  double precision g_R,trx2,RCUT
  double precision temperature,flux
  double precision acoeff_md(0:nremez_md)
  double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
  !***** output *****
  double complex delh_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double precision delh_alpha(1:nmat)
  !*******************
  double precision alpha_max,alpha_min
  double precision lattice_spacing,pi
  integer imat,jmat,kmat,lmat,imax,imin
  integer idim,jdim
  integer isite,isite_p1,isite_m1
  double complex MChi(1:nremez_md,1:nmat,1:nmat,1:nspin,&
       &-(nmargin-1):nsite+nmargin,1:npf)
  double complex delPF_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double precision delPF_alpha(1:nmat)
  integer iremez
  integer ispin
  integer ipf

  double complex Deriv_Mchi_alpha(1:nremez_md,1:nmat,1:nmat,1:nmat,&
       &1:nspin,1:nsite,1:npf)
  !number of indices must be less than 8...
  double complex Deriv_Mchi_X(1:nmat,1:nmat,1:ndim,&
       &1:nmat,1:nmat,1:nspin,1:nsite)

  pi=2d0*dasin(1d0)
  lattice_spacing=1d0/temperature/dble(nsite)

  delh_xmat=(0d0,0d0)
  delh_alpha=0d0
  !*******************************
  !***** pseudo-fermion part *****
  !*******************************
  if(ngauge.EQ.1)then
     !(dM/d¥alpha)*chi
     call Derivative_Dirac_alpha(alpha,&
          chi,Deriv_Mchi_alpha)
  end if
  !M*chi
  call Multiply_Dirac_to_chi(temperature,xmat,alpha,&
       chi,mchi,GAMMA10d,nbmn,flux)

  delPF_xmat=(0d0,0d0)
  delPF_alpha=0d0
  

  do ipf=1,npf
     do iremez=1,nremez_md 
        !(dM/dX)*chi  
        call Derivative_Dirac_X(iremez,ipf,temperature,&
             chi,Deriv_Mchi_X,GAMMA10d)  
        !-acoeff*(M*chi)^dagger*(dM/dX)*chi
        do kmat=1,nmat
           do lmat=1,nmat
              do idim=1,ndim
                 do isite=1,nsite
                    
                    do imat=1,nmat
                       do jmat=1,nmat
                          do ispin=1,nspin
                             delPF_xmat(kmat,lmat,idim,isite)=&
                                  delPF_xmat(kmat,lmat,idim,isite)&
                                  -dcmplx(acoeff_md(iremez))&
                                  *dconjg(Mchi(iremez,imat,jmat,ispin,isite,ipf))&
                                  *Deriv_Mchi_X(kmat,lmat,idim,imat,jmat,ispin,isite)

                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do
  !take into account both 
  !(M*chi)^dagger*(dM/dX)*chi and ((dM/dX)*chi)^dagger*(M*chi)
  do kmat=1,nmat
     do lmat=1,nmat
        do idim=1,ndim
           do isite=1,nsite
              delh_xmat(kmat,lmat,idim,isite)=&
                   delh_xmat(kmat,lmat,idim,isite)&
                   +delPF_xmat(kmat,lmat,idim,isite)&
                   +dconjg(delPF_xmat(lmat,kmat,idim,isite))
           end do
        end do
     end do
  end do
  if(ngauge.EQ.1)then
     !-acoeff*(M*chi)^dagger*(dM/alpha)*chi
     do iremez=1,nremez_md
        do ipf=1,npf
           do kmat=1,nmat
              do isite=1,nsite
                 do imat=1,nmat
                    do jmat=1,nmat
                       do ispin=1,nspin
                          delPF_alpha(kmat)=delPF_alpha(kmat)&
                               -acoeff_md(iremez)&
                               *dble(dconjg(Mchi(iremez,imat,jmat,ispin,isite,ipf))&
                               *Deriv_Mchi_alpha(iremez,kmat,imat,jmat,ispin,isite,ipf))
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
     !take into account both 
     !(M*chi)^dagger*(dM/d¥alpha)*chi and ((dM/d¥alpha)*chi)^dagger*(M*chi)
     do kmat=1,nmat
        delh_alpha(kmat)=delh_alpha(kmat)+delPF_alpha(kmat)*2d0
     end do
  end if


  return
  
END SUBROUTINE Calc_Force_Fermionic
