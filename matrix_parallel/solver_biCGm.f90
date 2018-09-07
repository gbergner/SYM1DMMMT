! gives error when pf=0.
!**************************************************
!**************************************************
!**************   Multi-mass BiCG   ***************
!**************************************************
!**************************************************
SUBROUTINE solver_biCGm(nbc,nremez,bcoeff,temperature,&
     &xmat,alpha,pf,chi,GAMMA10d,max_err,max_iteration,iteration,myrank,&
     &nbmn,flux,info)

  implicit none

  include 'mpif.h'
  include 'size_parallel.h'
  !***** input *****
  integer nbmn,nbc,nremez
  integer myrank
  double complex xmat(1:nmat_block,1:nmat_block,1:ndim,&
       &-(nmargin-1):nsite_local+nmargin)
  double precision alpha(1:nmat_block*nblock)
  double complex GAMMA10d(1:ndim,1:nspin,1:nspin)
  double precision temperature,flux
  double complex pf(1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin)
  integer max_iteration
  double precision bcoeff(1:nremez)!bcoeff(1) is the smallest. 
  double precision max_err
  !***** output *****
  double complex chi(1:nremez,1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin)
  integer iteration,info 
  !*******************
  integer info_temp
  integer imat,jmat
  integer ispin
  integer isite
  integer iremez
  double precision ERROR(1:nremez),NormR_local,NormR 
  integer hantei(1:nremez)

  double complex Xs0(1:nremez,1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin),&
       R0(1:nmat_block,1:nmat_block,1:nspin,-(nmargin-1):nsite_local+nmargin),&
       P0(1:nmat_block,1:nmat_block,1:nspin,-(nmargin-1):nsite_local+nmargin),&
       RB0(1:nmat_block,1:nmat_block,1:nspin,-(nmargin-1):nsite_local+nmargin),&
       PB0(1:nmat_block,1:nmat_block,1:nspin,-(nmargin-1):nsite_local+nmargin),&
       Ps0(1:nremez,1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin),&
       PBs0(1:nremez,1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin),&
       Xs1(1:nremez,1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin),&
       R1(1:nmat_block,1:nmat_block,1:nspin,-(nmargin-1):nsite_local+nmargin),&
       P1(1:nmat_block,1:nmat_block,1:nspin,-(nmargin-1):nsite_local+nmargin),&
       RB1(1:nmat_block,1:nmat_block,1:nspin,-(nmargin-1):nsite_local+nmargin),&
       PB1(1:nmat_block,1:nmat_block,1:nspin,-(nmargin-1):nsite_local+nmargin),&
       Ps1(1:nremez,1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin),&
       PBs1(1:nremez,1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin),&
       Xs2(1:nremez,1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin),&
       R2(1:nmat_block,1:nmat_block,1:nspin,-(nmargin-1):nsite_local+nmargin),&
       P2(1:nmat_block,1:nmat_block,1:nspin,-(nmargin-1):nsite_local+nmargin),&
       RB2(1:nmat_block,1:nmat_block,1:nspin,-(nmargin-1):nsite_local+nmargin),&
       PB2(1:nmat_block,1:nmat_block,1:nspin,-(nmargin-1):nsite_local+nmargin),&
       Ps2(1:nremez,1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin),&
       PBs2(1:nremez,1:nmat_block,1:nmat_block,1:nspin,&
       &-(nmargin-1):nsite_local+nmargin),&
       A(0:max_iteration),As(0:max_iteration,1:nremez),&
       B(0:max_iteration),Bs(0:max_iteration,1:nremez),&
       Zs(0:max_iteration,1:nremez),&
       numerator,denominator,bcoeffcmp(1:nremez),&
       numerator_local,denominator_local,&
       mp1(1:nmat_block,1:nmat_block,1:nspin,-(nmargin-1):nsite_local+nmargin)

  double precision NormPF2,normPF2_local
  double complex xmat_row(1:nmat_block,1:nmat_block*nblock,1:ndim,1:nsite_local)
  double complex xmat_column(1:nmat_block*nblock,1:nmat_block,1:ndim,&
       &1:nsite_local)
  !***** For MPI *****
  integer IERR

  bcoeffcmp=dcmplx(bcoeff)
  hantei=1
  normPF2_local=0d0
  !move i-th row and j-th row of xmat to (i,j)-th node.
  call mpi_xmat_row(xmat,xmat_row,myrank)
  call mpi_xmat_column(xmat,xmat_column,myrank)
  

  do isite=1,nsite_local
     do ispin=1,nspin
        do jmat=1,nmat_block
           do imat=1,nmat_block
              normPF2_local=normPF2_local+dble(pf(imat,jmat,ispin,isite)&
                   *dconjg((pf(imat,jmat,ispin,isite))))
           end do
        end do
     end do
  end do
  !collect normPF2_local to myrank=0 and calculate total value normPF2
  call MPI_Allreduce(normPF2_local,normPF2,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,MPI_COMM_WORLD,IERR)
  normPF2=dsqrt(normPF2)
  !write(*,*)normPF2_local
  !write(*,*)normPF2

  !initial condition
  do isite=0,nsite_local+1
     do ispin=1,nspin
!$omp parallel 
!$omp do
        do jmat=1,nmat_block       
           do imat=1,nmat_block
              
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
!$omp end do
!$omp end parallel
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
     denominator_local=(0d0,0d0)
     numerator_local=(0d0,0d0)
     do isite=1,nsite_local
        do ispin=1,nspin
           do jmat=1,nmat_block
              do imat=1,nmat_block
                 numerator_local=numerator_local+R1(imat,jmat,ispin,isite)&
                      *dconjg(R1(imat,jmat,ispin,isite))
              end do
           end do
        end do
     end do
     !collect numerator_local to myrank=0 and calculate total value numerator
     call MPI_Allreduce(numerator_local,numerator,1,MPI_DOUBLE_COMPLEX,&
          MPI_SUM,MPI_COMM_WORLD,IERR)
     !write(*,*)numerator
     !******************************
     !*** adjust margin and b.c. ***
     !******************************
     call Adjust_margin_and_bc_pf(P1,myrank,nbc)

     call Multiply_Dirac(temperature,xmat_row,xmat_column,&
          &alpha,P1,mp1,GAMMA10d,nbmn,flux,myrank)

     !******************************
     !*** adjust margin and b.c. ***
     !******************************
     call Adjust_margin_and_bc_pf(mp1,myrank,nbc)

     do isite=1,nsite_local
        do ispin=1,nspin
           do jmat=1,nmat_block      
              do imat=1,nmat_block 
                 denominator_local=denominator_local&
                      +mp1(imat,jmat,ispin,isite)*&
                      dconjg(mp1(imat,jmat,ispin,isite))
              end do
           end do
        end do
     end do 
     !collect denominator_local to myrank=0 and calculate total value denominator
     call MPI_Allreduce(denominator_local,denominator,1,MPI_DOUBLE_COMPLEX,&
          MPI_SUM,MPI_COMM_WORLD,IERR)

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
           do isite=1,nsite_local
              do ispin=1,nspin
!$omp parallel 
!$omp do
                 do jmat=1,nmat_block 
                    do imat=1,nmat_block 
                       Xs2(iremez,imat,jmat,ispin,isite)=&
                            Xs1(iremez,imat,jmat,ispin,isite)&
                            +As(iteration,iremez)&
                            *Ps1(iremez,imat,jmat,ispin,isite)   
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
           end do
        end if
     end do
!!!!!!!!!!!!!!!!!!!!!!!
     !update R and RB
     call Multiply_Dirac_dagger(temperature,xmat_row,xmat_column,&
          &alpha,mp1,R2,GAMMA10d,nbmn,flux,myrank)

     R2=R2*A(iteration)*(-1d0,0d0)+R1
     RB2=dconjg(R2)

     !update B
     denominator_local=(0d0,0d0)
     numerator_local=(0d0,0d0)
     do isite=1,nsite_local
        do ispin=1,nspin
           do jmat=1,nmat_block
              do imat=1,nmat_block  
                 denominator_local=denominator_local&
                      +R1(imat,jmat,ispin,isite)&
                      *dconjg(R1(imat,jmat,ispin,isite))
                 numerator_local=numerator_local&
                      +R2(imat,jmat,ispin,isite)&
                      *dconjg(R2(imat,jmat,ispin,isite))
              end do
           end do
        end do
     end do


     !collect denominator_local to myrank=0 and calculate total value denominator

     call MPI_Allreduce(denominator_local,denominator,1,MPI_DOUBLE_COMPLEX,&
          MPI_SUM,MPI_COMM_WORLD,IERR)
     !collect numerator_local to myrank=0 and calculate total value numerator
     call MPI_Allreduce(numerator_local,numerator,1,MPI_DOUBLE_COMPLEX,&
          MPI_SUM,MPI_COMM_WORLD,IERR)

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
           do isite=1,nsite_local
              do ispin=1,nspin
!$omp parallel 
!$omp do
                 do jmat=1,nmat_block
                    do imat=1,nmat_block
                 
                       !   Ps(ite+1,j,i)=Zs(ite+1,j)*R(ite+1,i)+Bs(ite,j)*Ps(ite,j,i)
                       Ps2(iremez,imat,jmat,ispin,isite)&
                            =Zs(iteration+1,iremez)*R2(imat,jmat,ispin,isite)&
                            +Bs(iteration,iremez)*Ps1(iremez,imat,jmat,ispin,isite)    

                       PBs2(iremez,imat,jmat,ispin,isite)&
                            =dconjg(Ps2(iremez,imat,jmat,ispin,isite))
                    end do
                 end do
!$omp end do
!$omp end parallel

              end do
           end do
        end if
     end do

     !Evaluation of err.
     !**evaluate ERROR for each iremez.
     !**if ERROR is smaller than Max_Err, then stop updating Zs.  
     !**This procedure is crucial!!
     NormR_local=0d0
     do isite=1,nsite_local
        do ispin=1,nspin
           do jmat=1,nmat_block
              do imat=1,nmat_block
                 NormR_local=NormR_local+dble(R2(imat,jmat,ispin,isite)&
                      *dconjg(R2(imat,jmat,ispin,isite)))
              end do
           end do
        end do
     end do

     call MPI_Reduce(NormR_local,NormR,1,MPI_DOUBLE_PRECISION,&
          MPI_SUM,0,MPI_COMM_WORLD,IERR)
     NormR=dsqrt(NormR) 

     call MPI_Bcast(NormR,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

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
           do isite=1,nsite_local
              do ispin=1,nspin
!$omp parallel 
!$omp do
                 do jmat=1,nmat_block
                    do imat=1,nmat_block
                       !   Chi(i,j)=Xs(ite+1,i,j)!final result.
                       Chi(iremez,imat,jmat,ispin,isite)&
                            =Xs2(iremez,imat,jmat,ispin,isite)!final result.
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
           end do
           hantei(iremez)=0
        end if
     end do
     iteration=iteration+1

     !if(myrank.EQ.0)then
     !   write(*,*)iteration,Error(1)
     !end if

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

   !******************************
   !*** adjust margin and b.c. ***
   !******************************
   call Adjust_margin_and_bc_Chi(Chi,myrank,nbc,nremez)



  return

END SUBROUTINE Solver_biCGm
