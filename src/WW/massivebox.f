      subroutine massivebox(k1,k2,k3,k4,k5,k6,za,zb,box)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'scale.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'docheck.f'
      include 'first.f'
      double complex d(2,2,6),box(2,2,-2:0),qlI4,Dint(6,-2:0)
      double precision s12,s34,s56,s134,s156,mtsq
      integer j,k1,k2,k3,k4,k5,k6,h1,h2,e
      common/transferbox/d
!$omp threadprivate(/transferbox/)

      mtsq=mt**2
      if (first) then
      first=.false. 
c      write(6,*) 'mtsq',mtsq
c      write(6,*) 'musq',musq
      call qlinit
c      pause
      endif

      s134=s(k1,k3)+s(k1,k4)+s(k3,k4)
      s156=s(k1,k5)+s(k1,k6)+s(k5,k6)
      s12=s(k1,k2)
      s34=s(k3,k4)
      s56=s(k5,k6)

c--- Initialize all box coefficients to zero
      do h1=1,2
      do h2=1,2
      do j=1,6
      d(h1,h2,j)=czip
      enddo
      enddo
      enddo

c--- Boxes 1 and 2      
      call box1(k1,k2,k3,k4,k5,k6,za,zb,d(2,2,1),d(2,1,1),
     &                                  d(1,2,1),d(1,1,1))
      call box1(k1,k2,k6,k5,k4,k3,zb,za,d(1,1,2),d(1,2,2),
     &                                  d(2,1,2),d(2,2,2))

c--- Box 3 and 4                                                                                                                                                            
      call box3(k1,k2,k3,k4,k5,k6,za,zb,d(2,2,3),d(2,1,3),
     &                                  d(1,2,3),d(1,1,3))
      call box3(k2,k1,k6,k5,k4,k3,zb,za,d(1,1,4),d(2,1,4),
     &                                  d(1,2,4),d(2,2,4))

c--- Box 4
c      call box4(k1,k2,k3,k4,k5,k6,za,zb,d(2,2,4),d(2,1,4),
c     &                                  d(1,2,4),d(1,1,4))
c      write(6,*) 'in massivebox:or:d(1,1,4)',d(1,1,4)
c      write(6,*) 'in massivebox:or:d(1,2,4)',d(1,2,4)
c      write(6,*) 'in massivebox:or:d(2,1,4)',d(2,1,4)
c      write(6,*) 'in massivebox:or:d(2,2,4)',d(2,2,4)
c      pause

c--- Boxes 5 and 6      
      call box5(k1,k2,k3,k4,k5,k6,za,zb,d(2,2,5),d(2,1,5),
     &                                  d(1,2,5),d(1,1,5))
      call box5(k1,k2,k6,k5,k4,k3,zb,za,d(1,1,6),d(1,2,6),
     &                                  d(2,1,6),d(2,2,6))

c--- compare with numerical code
      if (docheck) then
        call comparenum(4,6,d)
      endif 

      do e=-2,0
      Dint(1,e)=qlI4(s34,0d0,0d0,s56,s134,s12,mtsq,0d0,0d0,0d0,musq,e)
      enddo
      do e=-2,0
      Dint(2,e)=qlI4(s56,0d0,0d0,s34,s156,s12,mtsq,0d0,0d0,0d0,musq,e)
      enddo
      do e=-2,0
      Dint(3,e)=qlI4(s56,0d0,s34,0d0,s156,s134,0d0,mtsq,mtsq,0d0,musq,e)
      enddo
      do e=-2,0
      Dint(4,e)=qlI4(s56,0d0,s34,0d0,s156,s134,mtsq,0d0,0d0,mtsq,musq,e)
      enddo
      do e=-2,0
      Dint(5,e)=qlI4(s34,0d0,0d0,s56,s134,s12,0d0,mtsq,mtsq,mtsq,musq,e)
      enddo
      do e=-2,0
      Dint(6,e)=qlI4(s56,0d0,0d0,s34,s156,s12,0d0,mtsq,mtsq,mtsq,musq,e)
      enddo

      do h1=1,2
      do h2=1,2
      do e=-2,0
      box(h1,h2,e)=czip
      do j=1,6
      box(h1,h2,e)=box(h1,h2,e)+d(h1,h2,j)*Dint(j,e)
      enddo
      enddo
      enddo
      enddo

      if (docheck) then
        do j=1,6
        write(6,*) 'j,Dint(j)',j,Dint(j,0)*d(2,1,j)
        enddo
      endif
      
      return


      end
