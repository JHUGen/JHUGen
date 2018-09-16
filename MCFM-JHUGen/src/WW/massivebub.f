      subroutine massivebub(k1,k4,k2,k3,k5,k6,za,zb,bub)
      implicit none
      include 'types.f'
      
c----Bubble coefficients extracted from BDK 11.5, 11.8
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'scale.f'
      include 'blabels.f'
      include 'docheck.f'
      include 'first.f'
      character*9 st,st1
      integer:: j,k1,k2,k3,k4,k5,k6,h1,h2,e
      complex(dp):: b(2,2,7),bcoeff(8),bub(2,2,-2:0),Bint(7,-2:0),
     & qlI2,tmp
      real(dp):: mtsq,s12,s34,s56,s134,s156
      mtsq=mt**2

c--- QCDLoop already initialized from call to massivebox
      if (first) then
      first=.false. 
      call qlinit
c      write(6,*) 'mtsq',mtsq
c      write(6,*) 'musq',musq
c      pause
      endif

c--- note: these look unnatural due to the
c---       (k1,k4,k2,k3,k5,k6) ordering in the call
      s134=s(k1,k2)+s(k1,k3)+s(k2,k3)
      s156=s(k1,k5)+s(k1,k6)+s(k5,k6)
      s12=s(k1,k4)
      s34=s(k2,k3)
      s56=s(k5,k6)

      do h1=1,2
      do h2=1,2

      if ((h1==1) .and. (h2==1)) st='q+qb-g-g-' 
      if ((h1==1) .and. (h2==2)) st='q+qb-g-g+' 
      if ((h1==2) .and. (h2==1)) st='q+qb-g+g-' 
      if ((h1==2) .and. (h2==2)) st='q+qb-g+g+' 
      if (st == 'q+qb-g-g-') then
          st1='q+qb-g+g+'
          call mbc(st1,k4,k1,k3,k2,k6,k5,zb,za,bcoeff)
c--- this call interchanges roles of 1,2 so change coeffs. accordingly
          tmp=bcoeff(b134)
        bcoeff(b134)=bcoeff(b234)
        bcoeff(b234)=tmp
c--- (-,-) amplitudes
c      write(6,*) '--'
c      write(6,*) 'bcoeff(b12)',
c     & bcoeff(b12),abs(bcoeff(b12))
c      write(6,*) 'bcoeff(b12zm)',
c     & bcoeff(b12zm),abs(bcoeff(b12zm))
c      write(6,*) 'bcoeff(b234)',
c     & bcoeff(b234),abs(bcoeff(b234))
c      write(6,*) 'bcoeff(b56)',
c     & bcoeff(b56),abs(bcoeff(b56))
c      write(6,*) 'bcoeff(b134)',
c     & bcoeff(b134),abs(bcoeff(b134))
c      write(6,*) 'bcoeff(b34)',
c     & bcoeff(b34),abs(bcoeff(b34))
c      write(6,*) 'bcoeff(b12m)',
c     & bcoeff(b12m),abs(bcoeff(b12m))
c      write(6,*) 'sum',
c     & +bcoeff(b12)+bcoeff(b34)+bcoeff(b56)+bcoeff(b134)
c     & +bcoeff(b234)+bcoeff(b12m)
c      write(6,*) 'bcoeff(rat)',
c     & bcoeff(rat),abs(bcoeff(rat))
c      write(6,*)
        
      elseif (st=='q+qb-g-g+') then
          st1='q+qb-g+g-'
          call mbc(st1,k4,k1,k2,k3,k5,k6,za,zb,bcoeff)
c--- this call interchanges roles of 1,2 so change coeffs. accordingly
          tmp=bcoeff(b134)
        bcoeff(b134)=bcoeff(b234)
        bcoeff(b234)=tmp
c--- (-,+) amplitudes
c      write(6,*) '-+'
c      write(6,*) 'bcoeff(b12)',
c     & bcoeff(b12),abs(bcoeff(b12))
c      write(6,*) 'bcoeff(b12zm)',
c     & bcoeff(b12zm),abs(bcoeff(b12zm))
c      write(6,*) 'bcoeff(b234)',
c     & bcoeff(b234),abs(bcoeff(b234))
c      write(6,*) 'bcoeff(b56)',
c     & bcoeff(b56),abs(bcoeff(b56))
c      write(6,*) 'bcoeff(b134)',
c     & bcoeff(b134),abs(bcoeff(b134))
c      write(6,*) 'bcoeff(b34)',
c     & bcoeff(b34),abs(bcoeff(b34))
c      write(6,*) 'bcoeff(b12m)',
c     & bcoeff(b12m),abs(bcoeff(b12m))
c      write(6,*) 'sum',
c     & +bcoeff(b12)+bcoeff(b34)+bcoeff(b56)+bcoeff(b134)
c     & +bcoeff(b234)+bcoeff(b12m)
c      write(6,*) 'bcoeff(rat)',
c     & bcoeff(rat),abs(bcoeff(rat))
c      write(6,*)
        
      elseif (st=='q+qb-g+g+') then
          call mbc(st,k1,k4,k2,k3,k5,k6,za,zb,bcoeff)
c--- (+,+) amplitudes
c      write(6,*) '++'
c      write(6,*) 'bcoeff(b12)',
c     & bcoeff(b12),abs(bcoeff(b12))
c      write(6,*) 'bcoeff(b12zm)',
c     & bcoeff(b12zm),abs(bcoeff(b12zm))
c      write(6,*) 'bcoeff(b234)',
c     & bcoeff(b234),abs(bcoeff(b234))
c      write(6,*) 'bcoeff(b56)',
c     & bcoeff(b56),abs(bcoeff(b56))
c      write(6,*) 'bcoeff(b134)',
c     & bcoeff(b134),abs(bcoeff(b134))
c      write(6,*) 'bcoeff(b34)',
c     & bcoeff(b34),abs(bcoeff(b34))
c      write(6,*) 'bcoeff(b12m)',
c     & bcoeff(b12m),abs(bcoeff(b12m))
c      write(6,*) 'sum',
c     & +bcoeff(b12)+bcoeff(b34)+bcoeff(b56)+bcoeff(b134)
c     & +bcoeff(b234)+bcoeff(b12m)
c      write(6,*) 'bcoeff(rat)',
c     & bcoeff(rat),abs(bcoeff(rat))
c      write(6,*)
        
      elseif (st=='q+qb-g+g-') then
          call mbc(st,k1,k4,k2,k3,k5,k6,za,zb,bcoeff)
c--- (+,-) amplitudes
c      write(6,*) '+-'
c      write(6,*) 'bcoeff(b12)',
c     & bcoeff(b12),abs(bcoeff(b12))
c      write(6,*) 'bcoeff(b12zm)',
c     & bcoeff(b12zm),abs(bcoeff(b12zm))
c      write(6,*) 'bcoeff(b234)',
c     & bcoeff(b234),abs(bcoeff(b234))
c      write(6,*) 'bcoeff(b56)',
c     & bcoeff(b56),abs(bcoeff(b56))
c      write(6,*) 'bcoeff(b34)',
c     & bcoeff(b34),abs(bcoeff(b34))
c      write(6,*) 'bcoeff(b134)',
c     & bcoeff(b134),abs(bcoeff(b134))
c      write(6,*) 'bcoeff(b12m)',
c     & bcoeff(b12m),abs(bcoeff(b12m))
c      write(6,*) 'sum',
c     & +bcoeff(b12)+bcoeff(b34)+bcoeff(b56)+bcoeff(b134)
c     & +bcoeff(b234)+bcoeff(b12m)
c      write(6,*) 'bcoeff(rat)',
c     & bcoeff(rat),abs(bcoeff(rat))
c      write(6,*)
        
      endif

      do j=1,7
      b(h1,h2,j)=bcoeff(j)          
      enddo
      
      enddo
      enddo

c--- compare with numerical code
      if (docheck) then
        call comparenum(2,7,b)
      endif

      do h1=1,2
      do h2=1,2
      bub(h1,h2,-2)=czip
      enddo
      enddo
      

      do e=-1,0
      Bint(1,e)=qlI2(s12,0._dp,0._dp,musq,e)
      enddo
      do e=-1,0
      Bint(2,e)=qlI2(s34,0._dp,mtsq,musq,e)
      enddo
      do e=-1,0
      Bint(3,e)=qlI2(s56,0._dp,mtsq,musq,e)
      enddo
      do e=-1,0
      Bint(4,e)=qlI2(s134,0._dp,mtsq,musq,e)
      enddo
      do e=-1,0
      Bint(5,e)=qlI2(s156,0._dp,mtsq,musq,e)
      enddo
      do e=-1,0
      Bint(6,e)=qlI2(s12,mtsq,mtsq,musq,e)
      enddo
C-----rational part
      Bint(7,-1)=czip
      Bint(7, 0)=cone

c--- Note: Bint(8) is the massless result, should not be included in the sum
      do h1=1,2
      do h2=1,2
      do e=-2,0
      bub(h1,h2,e)=czip
      do j=1,7
      bub(h1,h2,e)=bub(h1,h2,e)+b(h1,h2,j)*Bint(j,e)
      enddo
      enddo
      enddo
      enddo

      if (docheck) then
        do j=1,7
        write(6,*) 'j,Bint(j)',j,Bint(j,0)*b(2,1,j)
        enddo
      endif
      
      return
      end

