SUBROUTINE who_am_i(myrank,isublat,iblock,jblock)

  implicit none
  include 'size_parallel.h'
  integer myrank,isublat,iblock,jblock,i

  i=myrank+1
  jblock=mod(i,nblock)
  if(jblock.EQ.0)then
     jblock=nblock
  end if
  i=int(dble(i-jblock)/dble(nblock)+0.01d0)+1
  iblock=mod(i,nblock)
  if(iblock.EQ.0)then
     iblock=nblock
  end if
  isublat=int(dble(i-iblock)/dble(nblock)+0.01d0)+1

  return

END SUBROUTINE who_am_i
