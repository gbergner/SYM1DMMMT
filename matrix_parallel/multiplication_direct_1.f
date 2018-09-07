!           1 -th nonzero term
 idim=           1
 ispin=           1
 jspin=           1
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
   do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !           2 -th nonzero term
 idim=           1
 ispin=           2
 jspin=           2
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !           3 -th nonzero term
 idim=           1
 ispin=           3
 jspin=           3
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !           4 -th nonzero term
 idim=           1
 ispin=           4
 jspin=           4
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !           5 -th nonzero term
 idim=           1
 ispin=           5
 jspin=           5
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !           6 -th nonzero term
 idim=           1
 ispin=           6
 jspin=           6
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !           7 -th nonzero term
 idim=           1
 ispin=           7
 jspin=           7
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !           8 -th nonzero term
 idim=           1
 ispin=           8
 jspin=           8
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !           9 -th nonzero term
 idim=           1
 ispin=           9
 jspin=           9
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          10 -th nonzero term
 idim=           1
 ispin=          10
 jspin=          10
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          11 -th nonzero term
 idim=           1
 ispin=          11
 jspin=          11
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          12 -th nonzero term
 idim=           1
 ispin=          12
 jspin=          12
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          13 -th nonzero term
 idim=           1
 ispin=          13
 jspin=          13
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          14 -th nonzero term
 idim=           1
 ispin=          14
 jspin=          14
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          15 -th nonzero term
 idim=           1
 ispin=          15
 jspin=          15
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          16 -th nonzero term
 idim=           1
 ispin=          16
 jspin=          16
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          17 -th nonzero term
 idim=           2
 ispin=           1
 jspin=          16
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          18 -th nonzero term
 idim=           2
 ispin=           2
 jspin=          15
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          19 -th nonzero term
 idim=           2
 ispin=           3
 jspin=          14
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          20 -th nonzero term
 idim=           2
 ispin=           4
 jspin=          13
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          21 -th nonzero term
 idim=           2
 ispin=           5
 jspin=          12
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          22 -th nonzero term
 idim=           2
 ispin=           6
 jspin=          11
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          23 -th nonzero term
 idim=           2
 ispin=           7
 jspin=          10
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          24 -th nonzero term
 idim=           2
 ispin=           8
 jspin=           9
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          25 -th nonzero term
 idim=           2
 ispin=           9
 jspin=           8
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          26 -th nonzero term
 idim=           2
 ispin=          10
 jspin=           7
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          27 -th nonzero term
 idim=           2
 ispin=          11
 jspin=           6
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          28 -th nonzero term
 idim=           2
 ispin=          12
 jspin=           5
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          29 -th nonzero term
 idim=           2
 ispin=          13
 jspin=           4
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          30 -th nonzero term
 idim=           2
 ispin=          14
 jspin=           3
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          31 -th nonzero term
 idim=           2
 ispin=          15
 jspin=           2
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          32 -th nonzero term
 idim=           2
 ispin=          16
 jspin=           1
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          33 -th nonzero term
 idim=           3
 ispin=           1
 jspin=          14
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          34 -th nonzero term
 idim=           3
 ispin=           2
 jspin=          13
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          35 -th nonzero term
 idim=           3
 ispin=           3
 jspin=          16
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          36 -th nonzero term
 idim=           3
 ispin=           4
 jspin=          15
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          37 -th nonzero term
 idim=           3
 ispin=           5
 jspin=          10
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          38 -th nonzero term
 idim=           3
 ispin=           6
 jspin=           9
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          39 -th nonzero term
 idim=           3
 ispin=           7
 jspin=          12
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          40 -th nonzero term
 idim=           3
 ispin=           8
 jspin=          11
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          41 -th nonzero term
 idim=           3
 ispin=           9
 jspin=           6
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          42 -th nonzero term
 idim=           3
 ispin=          10
 jspin=           5
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          43 -th nonzero term
 idim=           3
 ispin=          11
 jspin=           8
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          44 -th nonzero term
 idim=           3
 ispin=          12
 jspin=           7
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          45 -th nonzero term
 idim=           3
 ispin=          13
 jspin=           2
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          46 -th nonzero term
 idim=           3
 ispin=          14
 jspin=           1
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          47 -th nonzero term
 idim=           3
 ispin=          15
 jspin=           4
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          48 -th nonzero term
 idim=           3
 ispin=          16
 jspin=           3
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          49 -th nonzero term
 idim=           4
 ispin=           1
 jspin=          13
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          50 -th nonzero term
 idim=           4
 ispin=           2
 jspin=          14
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          51 -th nonzero term
 idim=           4
 ispin=           3
 jspin=          15
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          52 -th nonzero term
 idim=           4
 ispin=           4
 jspin=          16
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          53 -th nonzero term
 idim=           4
 ispin=           5
 jspin=           9
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          54 -th nonzero term
 idim=           4
 ispin=           6
 jspin=          10
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          55 -th nonzero term
 idim=           4
 ispin=           7
 jspin=          11
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          56 -th nonzero term
 idim=           4
 ispin=           8
 jspin=          12
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          57 -th nonzero term
 idim=           4
 ispin=           9
 jspin=           5
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          58 -th nonzero term
 idim=           4
 ispin=          10
 jspin=           6
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          59 -th nonzero term
 idim=           4
 ispin=          11
 jspin=           7
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          60 -th nonzero term
 idim=           4
 ispin=          12
 jspin=           8
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          61 -th nonzero term
 idim=           4
 ispin=          13
 jspin=           1
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          62 -th nonzero term
 idim=           4
 ispin=          14
 jspin=           2
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          63 -th nonzero term
 idim=           4
 ispin=          15
 jspin=           3
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          64 -th nonzero term
 idim=           4
 ispin=          16
 jspin=           4
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          65 -th nonzero term
 idim=           5
 ispin=           1
 jspin=          15
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          66 -th nonzero term
 idim=           5
 ispin=           2
 jspin=          16
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          67 -th nonzero term
 idim=           5
 ispin=           3
 jspin=          13
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          68 -th nonzero term
 idim=           5
 ispin=           4
 jspin=          14
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          69 -th nonzero term
 idim=           5
 ispin=           5
 jspin=          11
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          70 -th nonzero term
 idim=           5
 ispin=           6
 jspin=          12
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          71 -th nonzero term
 idim=           5
 ispin=           7
 jspin=           9
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          72 -th nonzero term
 idim=           5
 ispin=           8
 jspin=          10
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          73 -th nonzero term
 idim=           5
 ispin=           9
 jspin=           7
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          74 -th nonzero term
 idim=           5
 ispin=          10
 jspin=           8
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          75 -th nonzero term
 idim=           5
 ispin=          11
 jspin=           5
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          76 -th nonzero term
 idim=           5
 ispin=          12
 jspin=           6
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          77 -th nonzero term
 idim=           5
 ispin=          13
 jspin=           3
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          78 -th nonzero term
 idim=           5
 ispin=          14
 jspin=           4
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          79 -th nonzero term
 idim=           5
 ispin=          15
 jspin=           1
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          80 -th nonzero term
 idim=           5
 ispin=          16
 jspin=           2
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          81 -th nonzero term
 idim=           6
 ispin=           1
 jspin=          11
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          82 -th nonzero term
 idim=           6
 ispin=           2
 jspin=          12
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          83 -th nonzero term
 idim=           6
 ispin=           3
 jspin=           9
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          84 -th nonzero term
 idim=           6
 ispin=           4
 jspin=          10
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          85 -th nonzero term
 idim=           6
 ispin=           5
 jspin=          15
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          86 -th nonzero term
 idim=           6
 ispin=           6
 jspin=          16
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          87 -th nonzero term
 idim=           6
 ispin=           7
 jspin=          13
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          88 -th nonzero term
 idim=           6
 ispin=           8
 jspin=          14
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          89 -th nonzero term
 idim=           6
 ispin=           9
 jspin=           3
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          90 -th nonzero term
 idim=           6
 ispin=          10
 jspin=           4
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          91 -th nonzero term
 idim=           6
 ispin=          11
 jspin=           1
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          92 -th nonzero term
 idim=           6
 ispin=          12
 jspin=           2
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          93 -th nonzero term
 idim=           6
 ispin=          13
 jspin=           7
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          94 -th nonzero term
 idim=           6
 ispin=          14
 jspin=           8
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          95 -th nonzero term
 idim=           6
 ispin=          15
 jspin=           5
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          96 -th nonzero term
 idim=           6
 ispin=          16
 jspin=           6
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          97 -th nonzero term
 idim=           7
 ispin=           1
 jspin=          12
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          98 -th nonzero term
 idim=           7
 ispin=           2
 jspin=          11
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !          99 -th nonzero term
 idim=           7
 ispin=           3
 jspin=          10
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         100 -th nonzero term
 idim=           7
 ispin=           4
 jspin=           9
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         101 -th nonzero term
 idim=           7
 ispin=           5
 jspin=          16
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         102 -th nonzero term
 idim=           7
 ispin=           6
 jspin=          15
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         103 -th nonzero term
 idim=           7
 ispin=           7
 jspin=          14
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         104 -th nonzero term
 idim=           7
 ispin=           8
 jspin=          13
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         105 -th nonzero term
 idim=           7
 ispin=           9
 jspin=           4
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         106 -th nonzero term
 idim=           7
 ispin=          10
 jspin=           3
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         107 -th nonzero term
 idim=           7
 ispin=          11
 jspin=           2
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         108 -th nonzero term
 idim=           7
 ispin=          12
 jspin=           1
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         109 -th nonzero term
 idim=           7
 ispin=          13
 jspin=           8
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         110 -th nonzero term
 idim=           7
 ispin=          14
 jspin=           7
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         111 -th nonzero term
 idim=           7
 ispin=          15
 jspin=           6
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         112 -th nonzero term
 idim=           7
 ispin=          16
 jspin=           5
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         113 -th nonzero term
 idim=           8
 ispin=           1
 jspin=          10
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         114 -th nonzero term
 idim=           8
 ispin=           2
 jspin=           9
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         115 -th nonzero term
 idim=           8
 ispin=           3
 jspin=          12
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         116 -th nonzero term
 idim=           8
 ispin=           4
 jspin=          11
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         117 -th nonzero term
 idim=           8
 ispin=           5
 jspin=          14
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         118 -th nonzero term
 idim=           8
 ispin=           6
 jspin=          13
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         119 -th nonzero term
 idim=           8
 ispin=           7
 jspin=          16
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         120 -th nonzero term
 idim=           8
 ispin=           8
 jspin=          15
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         121 -th nonzero term
 idim=           8
 ispin=           9
 jspin=           2
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         122 -th nonzero term
 idim=           8
 ispin=          10
 jspin=           1
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         123 -th nonzero term
 idim=           8
 ispin=          11
 jspin=           4
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         124 -th nonzero term
 idim=           8
 ispin=          12
 jspin=           3
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         125 -th nonzero term
 idim=           8
 ispin=          13
 jspin=           6
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         126 -th nonzero term
 idim=           8
 ispin=          14
 jspin=           5
 Gam= ( -0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         127 -th nonzero term
 idim=           8
 ispin=          15
 jspin=           8
 Gam= (  0.0000000000000000     ,  1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         128 -th nonzero term
 idim=           8
 ispin=          16
 jspin=           7
 Gam= (  0.0000000000000000     , -1.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         129 -th nonzero term
 idim=           9
 ispin=           1
 jspin=           1
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         130 -th nonzero term
 idim=           9
 ispin=           2
 jspin=           2
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         131 -th nonzero term
 idim=           9
 ispin=           3
 jspin=           3
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         132 -th nonzero term
 idim=           9
 ispin=           4
 jspin=           4
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         133 -th nonzero term
 idim=           9
 ispin=           5
 jspin=           5
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         134 -th nonzero term
 idim=           9
 ispin=           6
 jspin=           6
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         135 -th nonzero term
 idim=           9
 ispin=           7
 jspin=           7
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         136 -th nonzero term
 idim=           9
 ispin=           8
 jspin=           8
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         137 -th nonzero term
 idim=           9
 ispin=           9
 jspin=           9
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         138 -th nonzero term
 idim=           9
 ispin=          10
 jspin=          10
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         139 -th nonzero term
 idim=           9
 ispin=          11
 jspin=          11
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         140 -th nonzero term
 idim=           9
 ispin=          12
 jspin=          12
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         141 -th nonzero term
 idim=           9
 ispin=          13
 jspin=          13
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         142 -th nonzero term
 idim=           9
 ispin=          14
 jspin=          14
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         143 -th nonzero term
 idim=           9
 ispin=          15
 jspin=          15
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
 !         144 -th nonzero term
 idim=           9
 ispin=          16
 jspin=          16
 Gam= ( -1.0000000000000000     , -0.0000000000000000     )
    do isite=1,nsite_local
!$omp parallel
!$omp do
                 do imat=1,nmat_block
                    do jmat=1,nmat_block
                       do kmat=1,nmat_block
                          pf2(imat,jmat,ispin,isite)=&
                               &pf2(imat,jmat,ispin,isite)-&
                               &dcmplx(lattice_spacing)*&
                               &GAM&
                               &*xmat_column((iblock-1)&
                               &*nmat_block+imat,kmat,idim,isite)&
                               &*pf1(kmat,jmat,jspin,isite)

                       end do
                    end do
                 end do
!$omp end do
!$omp end parallel
              end do
