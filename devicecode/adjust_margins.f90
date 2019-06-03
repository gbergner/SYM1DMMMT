! alpha -> alpha + P_alpha*dtau_alpha
! xmat -> xmat + P_xmat*dtau_xmat
! P_alpha -> P_alpha - delh_alpha*dtau_alpha
! P_xmat -> P_xmat - delh_xmat*dtau_xmat
! delh_xmat(imat,jmat)=dH/dxmat(jmat,imat)
MODULE Adjust_margins
contains
    SUBROUTINE Adjust_margin_xmat_device(xmat)

        use compiletimeconstants
        implicit none

        !***** input & output *****
        double complex xmat(1:nmat,1:nmat,1:ndim,-(nmargin-1):nsite+nmargin)
        !$acc declare device_resident(xmat)
        !**************************
        integer imat,jmat
        integer idim
        integer isite
  
        !$acc kernels
        do isite=-(nmargin-1),0
            do idim=1,ndim
                do imat=1,nmat
                    do jmat=1,nmat
                        xmat(imat,jmat,idim,isite)=xmat(imat,jmat,idim,nsite+isite)
                    end do
                end do
            end do
        end do
        !$acc end kernels
        !$acc kernels
        do isite=nsite+1,nsite+nmargin
            do idim=1,ndim
                do imat=1,nmat
                    do jmat=1,nmat
                        xmat(imat,jmat,idim,isite)=xmat(imat,jmat,idim,isite-nsite)
                    end do
                end do
            end do
        end do
        !$acc end kernels
        return
  
    END SUBROUTINE Adjust_margin_xmat_device
END MODULE Adjust_margins

