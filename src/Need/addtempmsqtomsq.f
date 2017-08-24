      subroutine addtemptomsq(msq,temp,ime,jme,it,jt,j,thescale)
      implicit none
      include 'constants.f'
      integer ime,jme,it,jt,j
      double precision msq(fn:nf,fn:nf),temp(fn:nf,fn:nf),thescale
      msq(ime,jme)=msq(ime,jme)+temp(it,jt)*thescale
      if(temp(it,jt) .ne. 0d0)
     & write(6,*) "Adding temp(",it,jt,") = ",temp(it,jt),
     & " to msq(",ime,jme,") at j=",j
      end
      subroutine addtempwtomsq(msq,tempw,ime,jme,it,jt,j,thescale)
      implicit none
      include 'constants.f'
      integer ime,jme,it,jt,j
      double precision msq(fn:nf,fn:nf),tempw(fn:nf,fn:nf),thescale
      msq(ime,jme)=msq(ime,jme)+tempw(it,jt)*thescale
      if(tempw(it,jt) .ne. 0d0)
     & write(6,*) "Adding tempw(",it,jt,") = ",tempw(it,jt),
     & " to msq(",ime,jme,") at j=",j
      end
