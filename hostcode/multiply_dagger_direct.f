 !           1 -th
 ispin=           1
 jspin=           1
 idim=           1
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

  do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !           2 -th
 ispin=           1
 jspin=           1
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !           3 -th
 ispin=           1
 jspin=          10
 idim=           8
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !           4 -th
 ispin=           1
 jspin=          11
 idim=           6
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !           5 -th
 ispin=           1
 jspin=          12
 idim=           7
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !           6 -th
 ispin=           1
 jspin=          13
 idim=           4
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !           7 -th
 ispin=           1
 jspin=          14
 idim=           3
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !           8 -th
 ispin=           1
 jspin=          15
 idim=           5
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !           9 -th
 ispin=           1
 jspin=          16
 idim=           2
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          10 -th
 ispin=           2
 jspin=           2
 idim=           1
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          11 -th
 ispin=           2
 jspin=           2
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          12 -th
 ispin=           2
 jspin=           9
 idim=           8
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          13 -th
 ispin=           2
 jspin=          11
 idim=           7
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          14 -th
 ispin=           2
 jspin=          12
 idim=           6
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          15 -th
 ispin=           2
 jspin=          13
 idim=           3
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          16 -th
 ispin=           2
 jspin=          14
 idim=           4
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          17 -th
 ispin=           2
 jspin=          15
 idim=           2
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          18 -th
 ispin=           2
 jspin=          16
 idim=           5
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          19 -th
 ispin=           3
 jspin=           3
 idim=           1
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          20 -th
 ispin=           3
 jspin=           3
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          21 -th
 ispin=           3
 jspin=           9
 idim=           6
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          22 -th
 ispin=           3
 jspin=          10
 idim=           7
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          23 -th
 ispin=           3
 jspin=          12
 idim=           8
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          24 -th
 ispin=           3
 jspin=          13
 idim=           5
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          25 -th
 ispin=           3
 jspin=          14
 idim=           2
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          26 -th
 ispin=           3
 jspin=          15
 idim=           4
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          27 -th
 ispin=           3
 jspin=          16
 idim=           3
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          28 -th
 ispin=           4
 jspin=           4
 idim=           1
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          29 -th
 ispin=           4
 jspin=           4
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          30 -th
 ispin=           4
 jspin=           9
 idim=           7
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          31 -th
 ispin=           4
 jspin=          10
 idim=           6
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          32 -th
 ispin=           4
 jspin=          11
 idim=           8
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          33 -th
 ispin=           4
 jspin=          13
 idim=           2
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          34 -th
 ispin=           4
 jspin=          14
 idim=           5
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          35 -th
 ispin=           4
 jspin=          15
 idim=           3
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          36 -th
 ispin=           4
 jspin=          16
 idim=           4
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          37 -th
 ispin=           5
 jspin=           5
 idim=           1
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          38 -th
 ispin=           5
 jspin=           5
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          39 -th
 ispin=           5
 jspin=           9
 idim=           4
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          40 -th
 ispin=           5
 jspin=          10
 idim=           3
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          41 -th
 ispin=           5
 jspin=          11
 idim=           5
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          42 -th
 ispin=           5
 jspin=          12
 idim=           2
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          43 -th
 ispin=           5
 jspin=          14
 idim=           8
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          44 -th
 ispin=           5
 jspin=          15
 idim=           6
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          45 -th
 ispin=           5
 jspin=          16
 idim=           7
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          46 -th
 ispin=           6
 jspin=           6
 idim=           1
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          47 -th
 ispin=           6
 jspin=           6
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          48 -th
 ispin=           6
 jspin=           9
 idim=           3
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          49 -th
 ispin=           6
 jspin=          10
 idim=           4
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          50 -th
 ispin=           6
 jspin=          11
 idim=           2
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          51 -th
 ispin=           6
 jspin=          12
 idim=           5
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          52 -th
 ispin=           6
 jspin=          13
 idim=           8
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          53 -th
 ispin=           6
 jspin=          15
 idim=           7
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          54 -th
 ispin=           6
 jspin=          16
 idim=           6
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          55 -th
 ispin=           7
 jspin=           7
 idim=           1
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          56 -th
 ispin=           7
 jspin=           7
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          57 -th
 ispin=           7
 jspin=           9
 idim=           5
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          58 -th
 ispin=           7
 jspin=          10
 idim=           2
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          59 -th
 ispin=           7
 jspin=          11
 idim=           4
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          60 -th
 ispin=           7
 jspin=          12
 idim=           3
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          61 -th
 ispin=           7
 jspin=          13
 idim=           6
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          62 -th
 ispin=           7
 jspin=          14
 idim=           7
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          63 -th
 ispin=           7
 jspin=          16
 idim=           8
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          64 -th
 ispin=           8
 jspin=           8
 idim=           1
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          65 -th
 ispin=           8
 jspin=           8
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          66 -th
 ispin=           8
 jspin=           9
 idim=           2
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          67 -th
 ispin=           8
 jspin=          10
 idim=           5
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          68 -th
 ispin=           8
 jspin=          11
 idim=           3
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          69 -th
 ispin=           8
 jspin=          12
 idim=           4
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          70 -th
 ispin=           8
 jspin=          13
 idim=           7
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          71 -th
 ispin=           8
 jspin=          14
 idim=           6
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          72 -th
 ispin=           8
 jspin=          15
 idim=           8
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          73 -th
 ispin=           9
 jspin=           2
 idim=           8
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          74 -th
 ispin=           9
 jspin=           3
 idim=           6
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          75 -th
 ispin=           9
 jspin=           4
 idim=           7
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          76 -th
 ispin=           9
 jspin=           5
 idim=           4
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          77 -th
 ispin=           9
 jspin=           6
 idim=           3
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          78 -th
 ispin=           9
 jspin=           7
 idim=           5
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          79 -th
 ispin=           9
 jspin=           8
 idim=           2
 gamconjg= ( -0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          80 -th
 ispin=           9
 jspin=           9
 idim=           1
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          81 -th
 ispin=           9
 jspin=           9
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          82 -th
 ispin=          10
 jspin=           1
 idim=           8
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          83 -th
 ispin=          10
 jspin=           3
 idim=           7
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          84 -th
 ispin=          10
 jspin=           4
 idim=           6
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          85 -th
 ispin=          10
 jspin=           5
 idim=           3
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          86 -th
 ispin=          10
 jspin=           6
 idim=           4
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          87 -th
 ispin=          10
 jspin=           7
 idim=           2
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          88 -th
 ispin=          10
 jspin=           8
 idim=           5
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          89 -th
 ispin=          10
 jspin=          10
 idim=           1
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          90 -th
 ispin=          10
 jspin=          10
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          91 -th
 ispin=          11
 jspin=           1
 idim=           6
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          92 -th
 ispin=          11
 jspin=           2
 idim=           7
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          93 -th
 ispin=          11
 jspin=           4
 idim=           8
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          94 -th
 ispin=          11
 jspin=           5
 idim=           5
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          95 -th
 ispin=          11
 jspin=           6
 idim=           2
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          96 -th
 ispin=          11
 jspin=           7
 idim=           4
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          97 -th
 ispin=          11
 jspin=           8
 idim=           3
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          98 -th
 ispin=          11
 jspin=          11
 idim=           1
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !          99 -th
 ispin=          11
 jspin=          11
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         100 -th
 ispin=          12
 jspin=           1
 idim=           7
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         101 -th
 ispin=          12
 jspin=           2
 idim=           6
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         102 -th
 ispin=          12
 jspin=           3
 idim=           8
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         103 -th
 ispin=          12
 jspin=           5
 idim=           2
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         104 -th
 ispin=          12
 jspin=           6
 idim=           5
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         105 -th
 ispin=          12
 jspin=           7
 idim=           3
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         106 -th
 ispin=          12
 jspin=           8
 idim=           4
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         107 -th
 ispin=          12
 jspin=          12
 idim=           1
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         108 -th
 ispin=          12
 jspin=          12
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         109 -th
 ispin=          13
 jspin=           1
 idim=           4
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         110 -th
 ispin=          13
 jspin=           2
 idim=           3
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         111 -th
 ispin=          13
 jspin=           3
 idim=           5
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         112 -th
 ispin=          13
 jspin=           4
 idim=           2
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         113 -th
 ispin=          13
 jspin=           6
 idim=           8
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         114 -th
 ispin=          13
 jspin=           7
 idim=           6
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         115 -th
 ispin=          13
 jspin=           8
 idim=           7
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         116 -th
 ispin=          13
 jspin=          13
 idim=           1
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         117 -th
 ispin=          13
 jspin=          13
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         118 -th
 ispin=          14
 jspin=           1
 idim=           3
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         119 -th
 ispin=          14
 jspin=           2
 idim=           4
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         120 -th
 ispin=          14
 jspin=           3
 idim=           2
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         121 -th
 ispin=          14
 jspin=           4
 idim=           5
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         122 -th
 ispin=          14
 jspin=           5
 idim=           8
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         123 -th
 ispin=          14
 jspin=           7
 idim=           7
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         124 -th
 ispin=          14
 jspin=           8
 idim=           6
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         125 -th
 ispin=          14
 jspin=          14
 idim=           1
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         126 -th
 ispin=          14
 jspin=          14
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         127 -th
 ispin=          15
 jspin=           1
 idim=           5
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         128 -th
 ispin=          15
 jspin=           2
 idim=           2
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         129 -th
 ispin=          15
 jspin=           3
 idim=           4
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         130 -th
 ispin=          15
 jspin=           4
 idim=           3
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         131 -th
 ispin=          15
 jspin=           5
 idim=           6
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         132 -th
 ispin=          15
 jspin=           6
 idim=           7
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         133 -th
 ispin=          15
 jspin=           8
 idim=           8
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         134 -th
 ispin=          15
 jspin=          15
 idim=           1
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         135 -th
 ispin=          15
 jspin=          15
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         136 -th
 ispin=          16
 jspin=           1
 idim=           2
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         137 -th
 ispin=          16
 jspin=           2
 idim=           5
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         138 -th
 ispin=          16
 jspin=           3
 idim=           3
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         139 -th
 ispin=          16
 jspin=           4
 idim=           4
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         140 -th
 ispin=          16
 jspin=           5
 idim=           7
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         141 -th
 ispin=          16
 jspin=           6
 idim=           6
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         142 -th
 ispin=          16
 jspin=           7
 idim=           8
 gamconjg= (  0.0000000000000000     ,  1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         143 -th
 ispin=          16
 jspin=          16
 idim=           1
 gamconjg= (  0.0000000000000000     , -1.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
  
 !         144 -th
 ispin=          16
 jspin=          16
 idim=           9
 gamconjg= ( -1.0000000000000000     ,  0.0000000000000000     )

   do isite=1,nsite
     do jmat=1,nmat
        do imat=1,nmat
           do kmat=1,nmat
              pf2(imat,jmat,ispin,isite)=&
                   pf2(imat,jmat,ispin,isite)&
                   -dcmplx(lattice_spacing)&
                   *gamconjg&
                   *(xmat(imat,kmat,idim,isite)*pf1(kmat,jmat,jspin,isite)&
                   -xmat(kmat,jmat,idim,isite)*pf1(imat,kmat,jspin,isite))
              
           end do
        end do
     end do
  end do
