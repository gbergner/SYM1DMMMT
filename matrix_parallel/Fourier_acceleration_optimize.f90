! info=1 -> measure fluctuation
! info=2 -> set the acceleration parameter
subroutine Fourier_acceleration_optimize(xmat_mom,fluctuation,myrank,imeasure,info)

  implicit none

  include 'mpif.h'
  include 'size_parallel.h'
  integer IERR

  integer myrank,info
  double precision fluctuation(1:nsite_local)
  double precision, allocatable :: fluctuation_all_1(:),fluctuation_all_2(:)
  double precision temp
  integer imom,imat,jmat,idim
  integer isublat,iblock,jblock,isite
  double complex xmat_mom(1:nmat_block,1:nmat_block,1:ndim,1:nsite_local)

  integer imeasure !number of measurement performed. 

  if(info.EQ.1)then
     !fluctuation(p)=|Xmat_mom(p)|
     do imom=1,nsite_local
        do imat=1,nmat_block
           do jmat=1,nmat_block
              do idim=1,ndim
                 fluctuation(imom)=fluctuation(imom)&
                      +dble(xmat_mom(imat,jmat,idim,imom)&
                      *dconjg(xmat_mom(imat,jmat,idim,imom)))
              end do
           end do
        end do
     end do
     imeasure=imeasure+1

  else if(info.EQ.2)then
     fluctuation=fluctuation/dble(nmat_block*nblock*imeasure)
     if(myrank.EQ.0)then
        allocate(fluctuation_all_1(1:nsite_local*nsublat*nblock*nblock))
        allocate(fluctuation_all_2(1:nsite_local*nsublat))
     end if

     call MPI_Gather(fluctuation(1),nsite_local,&
          MPI_DOUBLE_PRECISION,&
          fluctuation_all_1(1),nsite_local,MPI_DOUBLE_PRECISION,&
          0,MPI_COMM_WORLD,IERR)
     if(myrank.EQ.0)then
        fluctuation_all_2=0d0
        do isite=1,nsite_local
           do isublat=1,nsublat
              do iblock=1,nblock
                 do jblock=1,nblock
                    fluctuation_all_2(isite+(isublat-1)*nsite_local)=&
                         &fluctuation_all_2(isite+(isublat-1)*nsite_local)&
                         &+fluctuation_all_1(isite&
                         &+((iblock-1)*nblock+jblock-1)*nsite_local&
                         &+(isublat-1)*nsite_local*nblock*nblock)
                 end do
              end do
           end do
        end do
     end if

     if(myrank.EQ.0)then
        do imom=1,int(dble(nblock*nsite_local)*0.5d0+0.01d0)-1
           temp=fluctuation_all_2(imom)+&
                &fluctuation_all_2(nsublat*nsite_local-imom)
           temp=temp*0.5d0
           fluctuation_all_2(imom)=temp
           fluctuation_all_2(nsublat*nsite_local-imom)=temp
        end do
        fluctuation_all_2=dsqrt(fluctuation_all_2)


        do isite=1,nsite_local
           do isublat=1,nsublat
              do iblock=1,nblock
                 do jblock=1,nblock
                    fluctuation_all_1(isite&
                         &+((iblock-1)*nblock+jblock-1)*nsite_local&
                         &+(isublat-1)*nsite_local*nblock*nblock)=&
                         &fluctuation_all_2(isite+(isublat-1)*nsite_local)
                 end do
           end do
        end do
     end do
     

     end if
         call MPI_scatter(fluctuation_all_1(1),nsite_local,&
              MPI_DOUBLE_PRECISION,&
              fluctuation(1),nsite_local,MPI_DOUBLE_PRECISION,&
              0,MPI_COMM_WORLD,IERR)
     if(myrank.EQ.0)then
        deallocate(fluctuation_all_1)
        deallocate(fluctuation_all_2)
     end if

  end if

  return

END subroutine Fourier_acceleration_optimize
