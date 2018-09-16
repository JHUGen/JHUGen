      subroutine massivetri6(k1,k2,k3,k4,k5,k6,za,zb,triang)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'scale.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'docheck.f'
      complex(dp):: c(2,2,12),d(2,2,6),Cint(12,-2:0),
     & triang(2,2,-2:0),qlI3,tmp
      real(dp):: s12,s34,s56,s134,s156,mtsq,Delta,shift,
     & cred13,cred23
      integer:: j,k1,k2,k3,k4,k5,k6,e,h1,h2
      common/transferbox/d
!$omp threadprivate(/transferbox/)

      mtsq=mt**2

c--- QCDLoop already initialized from call to massivebox
c      if (first) then
c      first=.false. 
c      write(6,*) 'mtsq',mtsq
c      write(6,*) 'musq',musq
c      call qlinit
cc      pause
c      endif

      s134=s(k1,k3)+s(k1,k4)+s(k3,k4)
      s156=s(k2,k3)+s(k2,k4)+s(k3,k4)
      s12=s(k1,k2)
      s34=s(k3,k4)
      s56=s(k5,k6)

c--- Initialize all triangle coefficients to zero
      do h1=1,2
      do h2=1,2
      do j=1,12
      c(h1,h2,j)=czip
      enddo     
      enddo     
      enddo     

c--- Triangle 2  ! NO LONGER NEEDED
c      call triangle2(1,2,3,4,5,6,za,zb,c(2,2,2),c(2,1,2),
c     &                                 c(1,2,2),c(1,1,2))
     
c--- Triangle 6
      call triangle6(1,2,3,4,5,6,za,zb,c(2,2,6),c(2,1,6))
      call triangle6(2,1,3,4,5,6,za,zb,c(1,1,6),c(1,2,6))

c--- Triangles 7 and 9
c      call triangle7(1,2,3,4,5,6,za,zb,c(2,2,7),c(2,1,7))
c      call triangle7(1,2,6,5,4,3,zb,za,c(1,1,9),c(1,2,9))
      call triangle7new(1,2,3,4,5,6,za,zb,c(2,2,7),c(2,1,7))
      call triangle7new(1,2,6,5,4,3,zb,za,c(1,1,9),c(1,2,9))
      
c      call triangle9(1,2,3,4,5,6,za,zb,c(2,2,9),c(2,1,9))
c      call triangle9(1,2,6,5,4,3,zb,za,c(1,1,7),c(1,2,7))
      call triangle9new(1,2,3,4,5,6,za,zb,c(2,2,9),c(2,1,9))
      call triangle9new(1,2,6,5,4,3,zb,za,c(1,1,7),c(1,2,7))

c--- Triangles 8 and 10
c      call triangle9(2,1,3,4,5,6,za,zb,c(2,2,8),c(1,2,8))
c      call triangle9(2,1,6,5,4,3,zb,za,c(1,1,10),c(2,1,10))
      call triangle9new(2,1,3,4,5,6,za,zb,c(2,2,8),c(1,2,8))
      call triangle9new(2,1,6,5,4,3,zb,za,c(1,1,10),c(2,1,10))
      
c      call triangle7(2,1,3,4,5,6,za,zb,c(2,2,10),c(1,2,10))
c      call triangle7(2,1,6,5,4,3,zb,za,c(1,1,8),c(2,1,8))
      call triangle7new(2,1,3,4,5,6,za,zb,c(2,2,10),c(1,2,10))
      call triangle7new(2,1,6,5,4,3,zb,za,c(1,1,8),c(2,1,8))

c--- Triangle 11
c      call triangle11(1,2,3,4,5,6,za,zb,c(2,2,11),c(2,1,11))
c--- This works for mp
c      call triangle11(2,1,3,4,5,6,za,zb,tmp,c(1,2,11))
c--- This works for mm
c      call triangle11(2,1,4,3,6,5,zb,za,c(1,1,11),tmp)
      call triangle11new(1,2,3,4,5,6,za,zb,c(2,2,11),c(2,1,11))
c--- This works for mp
      call triangle11new(2,1,3,4,5,6,za,zb,tmp,c(1,2,11))
c--- This works for mm
      call triangle11new(2,1,4,3,6,5,zb,za,c(1,1,11),tmp)

c--- Triangle 12
      call triangle12(1,2,3,4,5,6,za,zb,c(2,2,12),c(2,1,12))
      call triangle12(1,2,6,5,4,3,zb,za,c(1,1,12),c(1,2,12))
c      call triangle12new(1,2,3,4,5,6,za,zb,c(2,2,12),c(2,1,12))
c      call triangle12new(1,2,6,5,4,3,zb,za,c(1,1,12),c(1,2,12))


      

c--- Perform shifting of triangle coefficients 
      Delta=(mtsq-s134)*(mtsq-s156)-(mtsq-s34)*(mtsq-s56)
      shift=(2._dp*s12*mtsq-Delta)/Delta**2
c      cred31=2._dp*(s134-s34)*shift
c      cred34=2._dp*(s156-s56)*shift
c      cred42=2._dp*(s156-s34)*shift
c      cred43=2._dp*(s134-s56)*shift
       cred13=2._dp*((s134+mtsq)*(s56+s34-s12)-2._dp*s134*mtsq-2._dp*s34*s56)
     & /(s12*(s134-mtsq)**2)
       cred23=2._dp*((s156+mtsq)*(s56+s34-s12)-2._dp*s156*mtsq-2._dp*s34*s56)
     & /(s12*(s156-mtsq)**2)


      do h1=1,2
      do h2=1,2
      c(h1,h2,6)=c(h1,h2,6)
     & -0.5_dp*cred13*d(h1,h2,1)-0.5_dp*cred23*d(h1,h2,2)
      c(h1,h2,7) = c(h1,h2,7)-(s134-s34)*shift*d(h1,h2,3)
      c(h1,h2,8) = c(h1,h2,8)-(s134-s56)*shift*d(h1,h2,4)
      c(h1,h2,9) = c(h1,h2,9)-(s156-s56)*shift*d(h1,h2,3)
      c(h1,h2,10)=c(h1,h2,10)-(s156-s34)*shift*d(h1,h2,4)
      enddo
      enddo
      
c--- compare with numerical code
      if (docheck) then
        call comparenum(3,12,c)    
      endif

c      do e=-2,0
c      Cint(1,e)=qlI3(0._dp,0._dp,s12,0._dp,0._dp,0._dp,musq,e)
c      enddo
c      do e=-2,0
c      Cint(2,e)=qlI3(0._dp,s134,s34,0._dp,0._dp,mtsq,musq,e)
c      enddo
c      do e=-2,0
c      Cint(3,e)=qlI3(0._dp,s56,s134,0._dp,0._dp,mtsq,musq,e)
c      enddo
c      do e=-2,0
c      Cint(4,e)=qlI3(0._dp,s156,s56,0._dp,0._dp,mtsq,musq,e)
c      enddo
c      do e=-2,0
c      Cint(5,e)=qlI3(0._dp,s34,s156,0._dp,0._dp,mtsq,musq,e)
c      enddo
      do e=-2,0
      Cint(6,e)=qlI3(s12,s56,s34,0._dp,0._dp,mtsq,musq,e)
      enddo
      do e=-2,0
      Cint(7,e)=qlI3(s134,0._dp,s34,0._dp,mtsq,mtsq,musq,e)
      enddo
      do e=-2,0
      Cint(8,e)=qlI3(s56,0._dp,s134,0._dp,mtsq,mtsq,musq,e)
      enddo
      do e=-2,0
      Cint(9,e)=qlI3(s156,0._dp,s56,0._dp,mtsq,mtsq,musq,e)
      enddo
      do e=-2,0
      Cint(10,e)=qlI3(s34,0._dp,s156,0._dp,mtsq,mtsq,musq,e)
      enddo
      do e=-2,0
         Cint(11,e)=qlI3(s56,s12,s34,0._dp,mtsq,mtsq,musq,e)
      enddo
      do e=-2,0
      Cint(12,e)=qlI3(s12,0._dp,0._dp,mtsq,mtsq,mtsq,musq,e)
      enddo

      do h1=1,2
      do h2=1,2
      do e=-2,0
      triang(h1,h2,e)=czip
      do j=6,12
      triang(h1,h2,e)=triang(h1,h2,e)+c(h1,h2,j)*Cint(j,e)
      enddo
      enddo
      enddo
      enddo

      if (docheck) then
        do j=1,12
        write(6,*) 'j,Cint(j)',j,Cint(j,0)*c(2,1,j)
        enddo
      endif
      
      return
      end

