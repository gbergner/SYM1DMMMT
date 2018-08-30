! gives error when pf=0.
!**************************************************
!**************************************************
!**************   Multi-mass BiCG   ***************
!**************************************************
!**************************************************
SUBROUTINE solver_biCGm(nbc,nbmn,nremez,&
     xmat,alpha,pf_input,chi,GAMMA10d,&
     bcoeff,max_err,max_iteration,iteration,&
     temperature,flux,info)

  implicit none

  include '../staticparameters.f90'

  integer nbmn,nremez,nbc

  integer info_temp

  double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
  double precision alpha(1:nmat)
  double complex pf_input(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)!input
  double complex pf(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)
  double complex chi(1:nremez,1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin,1:npf)!output
  double complex GAMMA10d(1:ndim,1:nspin,1:nspin)

  double precision temperature,flux
  integer imat,jmat
  integer ispin
  integer isite
  integer iremez

  integer max_iteration,iteration,info
  double precision bcoeff(1:nremez)!bcoeff(1) is the smallest. 
  double precision max_err
  double precision ERROR(1:nremez),NormR 

  integer hantei(1:nremez)

  double complex Xs0(1:nremez,1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       R0(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       P0(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       RB0(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       PB0(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       Ps0(1:nremez,1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       PBs0(1:nremez,1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       Xs1(1:nremez,1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       R1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       P1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       RB1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       PB1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       Ps1(1:nremez,1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       PBs1(1:nremez,1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       Xs2(1:nremez,1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       R2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       P2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       RB2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       PB2(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       Ps2(1:nremez,1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       PBs2(1:nremez,1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin),&
       A(0:max_iteration),As(0:max_iteration,1:nremez),&
       B(0:max_iteration),Bs(0:max_iteration,1:nremez),&
       Zs(0:max_iteration,1:nremez),&
       numerator,denominator,bcoeffcmp(1:nremez),&
       mp1(1:nmat,1:nmat,1:nspin,-(nmargin-1):nsite+nmargin)

  double precision NormPF2
  integer ipf
  bcoeffcmp=dcmplx(bcoeff)

  do ipf=1,npf
     do isite=-(nmargin-1),nsite+nmargin
        do ispin=1,nspin
           do jmat=1,nmat
              do imat=1,nmat 
                 pf(imat,jmat,ispin,isite)=pf_input(imat,jmat,ispin,isite,ipf)
              end do
           end do
        end do
     end do
  
     hantei=1
     normPF2=0d0
     do isite=1,nsite
        do ispin=1,nspin
           do jmat=1,nmat
              do imat=1,nmat
                 normPF2=normPF2+dble(pf(imat,jmat,ispin,isite)&
                      *dconjg((pf(imat,jmat,ispin,isite))))
              end do
           end do
        end do
     end do
     normPF2=dsqrt(normPF2)
     if (solver_verbose.EQ.3) then
       write(*,*) "norm host(",ipf,")",normPF2
     end if
     !initial condition
     do isite=-(nmargin-1),nsite+nmargin
        do ispin=1,nspin
           !$omp parallel do
           do jmat=1,nmat       
              do imat=1,nmat
                 
                 do iremez=1,nremez             
                    Xs1(iremez,imat,jmat,ispin,isite)=(0d0,0d0)
                    Ps0(iremez,imat,jmat,ispin,isite)=(0d0,0d0)
                    PBs0(iremez,imat,jmat,ispin,isite)=(0d0,0d0)
                    Ps1(iremez,imat,jmat,ispin,isite)=PF(imat,jmat,ispin,isite)
                    PBs1(iremez,imat,jmat,ispin,isite)=dconjg(PF(imat,jmat,ispin,isite))
                 end do
                 
                 P0(imat,jmat,ispin,isite)=(0d0,0d0)
                 PB0(imat,jmat,ispin,isite)=(0d0,0d0)
                 R0(imat,jmat,ispin,isite)=pf(imat,jmat,ispin,isite)
                 RB0(imat,jmat,ispin,isite)=dconjg(pf(imat,jmat,ispin,isite))
                 R1(imat,jmat,ispin,isite)=pf(imat,jmat,ispin,isite)
                 RB1(imat,jmat,ispin,isite)=dconjg(pf(imat,jmat,ispin,isite))
                 P1(imat,jmat,ispin,isite)=pf(imat,jmat,ispin,isite)
                 PB1(imat,jmat,ispin,isite)=dconjg(pf(imat,jmat,ispin,isite))
                 
              end do
           end do
        end do
     end do
     
     
     do iremez=1,nremez 
        Zs(0,iremez)=(1d0,0d0)
        Zs(1,iremez)=(1d0,0d0)
        As(0,iremez)=(1d0,0d0)
        Bs(0,iremez)=(1d0,0d0)
     end do
     A(0)=(1d0,0d0)
     B(0)=(1d0,0d0)
     !loop
     iteration=1
     info=1
     do while((iteration.LT.max_iteration).AND.(info.EQ.1))
        !update A
        denominator=(0d0,0d0)
        numerator=(0d0,0d0)
        do isite=1,nsite
           do ispin=1,nspin
              do jmat=1,nmat
                 do imat=1,nmat
                    numerator=numerator+R1(imat,jmat,ispin,isite)&
                         *dconjg(R1(imat,jmat,ispin,isite))
                 end do
              end do
           end do
        end do
        
        !boundary condition  
        
        if(nbc.EQ.0)then!pbc
           do ispin=1,nspin
              do jmat=1,nmat     
                 do imat=1,nmat
                    do isite=-(nmargin-1),0
                       P1(imat,jmat,ispin,isite)=&
                            P1(imat,jmat,ispin,nsite+isite)
                    end do
                 end do
              end do
           end do
           do ispin=1,nspin
              do jmat=1,nmat     
                 do imat=1,nmat
                    do isite=nsite+1,nsite+nmargin
                       P1(imat,jmat,ispin,isite)=&
                            P1(imat,jmat,ispin,isite-nsite)
                    end do
                 end do
              end do
           end do
           
        else if(nbc.EQ.1)then!apbc
           do ispin=1,nspin
              do jmat=1,nmat     
                 do imat=1,nmat
                    do isite=-(nmargin-1),0
                       P1(imat,jmat,ispin,isite)=&
                            P1(imat,jmat,ispin,nsite+isite)*(-1d0,0d0)
                    end do
                 end do
              end do
           end do
           do ispin=1,nspin
              do jmat=1,nmat     
                 do imat=1,nmat
                    do isite=nsite+1,nsite+nmargin
                       P1(imat,jmat,ispin,isite)=&
                            P1(imat,jmat,ispin,isite-nsite)*(-1d0,0d0)
                    end do
                 end do
              end do
           end do
        end if
        
        call Multiply_Dirac(temperature,xmat,alpha,P1,mp1,GAMMA10d,nbmn,flux)
        
        
        !boundary condition  
        
        if(nbc.EQ.0)then!pbc
           do ispin=1,nspin
              do jmat=1,nmat     
                 do imat=1,nmat
                    do isite=-(nmargin-1),0
                       mp1(imat,jmat,ispin,isite)=&
                            mp1(imat,jmat,ispin,isite+nsite)
                    end do
                 end do
              end do
           end do
           do ispin=1,nspin
              do jmat=1,nmat     
                 do imat=1,nmat
                    do isite=nsite+1, nsite+nmargin
                       mp1(imat,jmat,ispin,isite)=&
                            mp1(imat,jmat,ispin,isite-nsite)
                    end do
                 end do
              end do
           end do
        else if(nbc.EQ.1)then!apbc
           do ispin=1,nspin
              do jmat=1,nmat     
                 do imat=1,nmat
                    do isite=-(nmargin-1),0
                       mp1(imat,jmat,ispin,isite)=&
                            mp1(imat,jmat,ispin,isite+nsite)*(-1d0,0d0)
                    end do
                 end do
              end do
           end do
           do ispin=1,nspin
              do jmat=1,nmat     
                 do imat=1,nmat
                    do isite=nsite+1, nsite+nmargin
                       mp1(imat,jmat,ispin,isite)=&
                            mp1(imat,jmat,ispin,isite-nsite)*(-1d0,0d0)
                    end do
                 end do
              end do
           end do
        end if
        do isite=1,nsite
           do ispin=1,nspin
              do jmat=1,nmat     
                 do imat=1,nmat
                    denominator=denominator&
                         +mp1(imat,jmat,ispin,isite)*&
                         dconjg(mp1(imat,jmat,ispin,isite))
                 end do
              end do
           end do
        end do
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        A(iteration)=numerator/denominator
        !write(*,*)denominator
        !write(*,*)A(iteration)
        !update Zs and As

        do iremez=1,nremez
           if(hantei(iremez).EQ.0)then!not update. already finished. 
              Zs(iteration+1,iremez)=Zs(iteration,iremez)
              As(iteration,iremez)=A(iteration)
           else
              Zs(iteration+1,iremez)=&
                   Zs(iteration,iremez)*Zs(iteration-1,iremez)*A(iteration-1)&
                   /(A(iteration-1)*Zs(iteration-1,iremez)*((1d0,0d0)&
                   +A(iteration)*bcoeffcmp(iremez))&
                   +A(iteration)*B(iteration-1)*(Zs(iteration-1,iremez)&
                   -Zs(iteration,iremez)))
              
              As(iteration,iremez)=A(iteration)*Zs(iteration+1,iremez)&
                   /Zs(iteration,iremez)
           end if
        end do
        

        
        
        !update Xs
        do iremez=1,nremez
           if(hantei(iremez).EQ.1)then 
              do isite=1,nsite
                 do ispin=1,nspin
                    do jmat=1,nmat
                       do imat=1,nmat
                          Xs2(iremez,imat,jmat,ispin,isite)=&
                               Xs1(iremez,imat,jmat,ispin,isite)&
                               +As(iteration,iremez)&
                               *Ps1(iremez,imat,jmat,ispin,isite)   
                       end do
                    end do
                 end do
              end do
           end if
        end do
        
!!!!!!!!!!!!!!!!!!!!!!!
        !update R and RB
        
        call Multiply_Dirac_dagger(temperature,xmat,alpha,mp1,R2,GAMMA10d,nbmn,flux)

        
        R2=R2*A(iteration)*(-1d0,0d0)+R1
        RB2=dconjg(R2)
        
        !update B
        denominator=(0d0,0d0)
        numerator=(0d0,0d0)
        do isite=1,nsite
           do ispin=1,nspin
              do jmat=1,nmat
                 do imat=1,nmat  
                    denominator=denominator&
                         +R1(imat,jmat,ispin,isite)&
                         *dconjg(R1(imat,jmat,ispin,isite))
                    numerator=numerator&
                         +R2(imat,jmat,ispin,isite)&
                         *dconjg(R2(imat,jmat,ispin,isite))
                 end do
              end do
           end do
        end do
        
        B(iteration)=numerator/denominator
                                        
        !update Bs
        do iremez=1,nremez
           if(hantei(iremez).EQ.1)then
              ! DO NOT write as "(Zs(ite+1,i)/Zs(ite,i))**2d0". It causes err.
              Bs(iteration,iremez)=B(iteration)&
                   *(Zs(iteration+1,iremez)&
                   /Zs(iteration,iremez))*(Zs(iteration+1,iremez)&
                   /Zs(iteration,iremez))
           end if
        end do
        !update P, PB, Ps and PBs. 
        
        P2=R2+B(iteration)*P1
        PB2=dconjg(P2)
        
        do iremez=1,nremez
           if(hantei(iremez).EQ.1)then
              do isite=1,nsite
                 do ispin=1,nspin
                    !$omp parallel do
                    do jmat=1,nmat
                       do imat=1,nmat
                          
                          !   Ps(ite+1,j,i)=Zs(ite+1,j)*R(ite+1,i)+Bs(ite,j)*Ps(ite,j,i)
                          Ps2(iremez,imat,jmat,ispin,isite)&
                               =Zs(iteration+1,iremez)*R2(imat,jmat,ispin,isite)&
                               +Bs(iteration,iremez)*Ps1(iremez,imat,jmat,ispin,isite)    
                          
                          PBs2(iremez,imat,jmat,ispin,isite)&
                               =dconjg(Ps2(iremez,imat,jmat,ispin,isite))
                       end do
                    end do
                    
                 end do
              end do
           end if
        end do
        
        !Evaluation of err.
        !**evaluate ERROR for each iremez.
        !**if ERROR is smaller than Max_Err, then stop updating Zs.  
        !**This procedure is crucial!!
        NormR=0d0
        do isite=1,nsite
           do ispin=1,nspin
              do jmat=1,nmat
                 do imat=1,nmat
                    NormR=NormR+dble(R2(imat,jmat,ispin,isite)&
                         *dconjg(R2(imat,jmat,ispin,isite)))
                 end do
              end do
           end do
        end do
        NormR=dsqrt(NormR)

        !write(*,*)hantei
        do iremez=1,nremez
           if(hantei(iremez).EQ.1)then
              ERROR(iremez)=dsqrt(dble(Zs(iteration+1,iremez)&
                   *dconjg(Zs(iteration+1,iremez))))&
                   *NormR/normPF2
              !write(*,*)normPF2
           end if
        end do
        do iremez=1,nremez
           if((ERROR(iremez).LT.max_err).AND.(hantei(iremez).EQ.1))then
              do isite=1,nsite
                 do ispin=1,nspin
                    do jmat=1,nmat
                       do imat=1,nmat
                          !   Chi(i,j)=Xs(ite+1,i,j)!final result.
                          Chi(iremez,imat,jmat,ispin,isite,ipf)&
                               =Xs2(iremez,imat,jmat,ispin,isite)!final result.
                       end do
                    end do
                 end do
              end do
              hantei(iremez)=0
           end if
        end do
        if (solver_verbose.EQ.3) then
        write(*,*) "iteration ",iteration," normr host ",normr
                  write(*,*)  "iteration",iteration," errorfull",Error
          write(*,*) "hantei",hantei

          NormR=0.0
              do isite=1,nsite
                 do ispin=1,nspin
                    do jmat=1,nmat
                       do imat=1,nmat
                               NormR=NormR+xs2(1,imat,jmat,ispin,isite)
                       end do
                    end do
                 end do
              end do
          
                   write(*,*) "xsnorm device",NormR
        end if
        iteration=iteration+1
        
        !write(*,*)iteration,Error(1)
        
        Xs0=Xs1
        R0=R1
        P0=P1
        RB0=RB1
        PB0=PB1
        Ps0=Ps1
        PBs0=PBs1  
        
        Xs1=Xs2
        R1=R2
        P1=P2
        RB1=RB2
        PB1=PB2
        Ps1=Ps2
        PBs1=PBs2
        
        info_temp=0
        do iremez=1,nremez
           if(hantei(iremez).EQ.1)then
              info_temp=1
           end if
        end do
        info=info_temp
        
     end do
     
     
     !boundary condition  
     
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

END SUBROUTINE Solver_biCGm
