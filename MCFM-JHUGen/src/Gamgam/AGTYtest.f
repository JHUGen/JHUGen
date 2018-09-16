
      subroutine AGTYtest(ss,tt,uu)
!====== C.Williams July 2015
!==== routine for testing the various functions needed to
!==== build the two-loop GamGam and DirGam results
!===== Check is against independent Mathematica routines

!==== note that the various pieces must be tested for the
!==== physical region, s > 0, t< 0, u < 0
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      include 'scale.f'
      real(dp):: ss,tt,uu
      real(dp):: ddilog,Li3,Li4
      real(dp):: BigX,BigY,BigS,BigU,BigT,x,y,z
!==== routines
      real(dp):: Asx,Bsx,Csx,D1sx,D2sx,E1sx,E2sx,E3sx
      real(dp):: AGTYAs,AGTYBs,AGTYCs,AGTYD1s,AGTYD2s,AGTYE1s,AGTYE2s,
     & AGTYE3s
      real(dp) :: F1sx,F2sx
      real(dp) :: AGTYF1s,AGTYF2s
      real(dp):: test
      real(dp) Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,Li2z,Li3z,Li4z,Li4zinv

!==== initializtaion
!===== set scale
      musq=1_dp

!----- define various pieces
      BigX=log(-tt/ss)
      BigY=log(-uu/ss)
      BigS=log(ss/musq)
      BigU=log(-uu/musq)
      BigT=log(-tt/musq)

      x=-tt/ss
      y=-uu/ss
      z=-uu/tt

      Li2x=ddilog(x)
      Li3x=Li3(x)
      Li4x=Li4(x)

      Li2y=ddilog(y)
      Li3y=Li3(y)
      Li4y=Li4(y)

      Li2z=ddilog(z)
      Li3z=Li3(z)
      Li4z=Li4(z)

      Li4zinv=Li4(1_dp/z)
      write(6,*) '********** AGTY TEST **************'
      write(6,*) '* points :                         *'
      write(6,*) 'ss = ',ss
      write(6,*) 'tt = ',tt
      write(6,*) 'uu = ',uu
      write(6,*) '* Big Letters :                         *'
      write(6,*) 'Big X= ',BigX
      write(6,*) 'Big Y = ',BigY
      write(6,*) 'Big S = ',BigS
      write(6,*) 'Big U = ',BigU
      write(6,*) 'Big T = ',BigT
      write(6,*) '* lil Letters :                         *'
      write(6,*) 'x  =',x
      write(6,*) 'y = ',y
      write(6,*) 'z = ',z
      write(6,*) '* PolyLogs :                         *'
      write(6,*) 'Li4(x)  =',Li4x
      write(6,*) 'Li3(x) = ',Li3x
      write(6,*) 'Li2(x) = ',Li2x
      write(6,*) 'Li4(y)  =',Li4y
      write(6,*) 'Li3(y) = ',Li3y
      write(6,*) 'Li2(y) = ',Li2y
      write(6,*) 'Li4(z)  =',Li4z
      write(6,*) 'Li3(z) = ',Li3z
      write(6,*) 'Li2(z) = ',Li2z
      write(6,*) '************************************'


 !      test=Asx(s,t,u,BigX,Ly,Li2x,Li3x,Li3y,Li4x,Li4y,Li4z)
      write(6,*) 'A routines '
      write(6,*) 'Asx ',
     & Asx(ss,tt,uu,BigX,BigY,Li2x,Li3x,Li3y,Li4x,Li4y,Li4z)
      write(6,*) 'As ',
     & AGTYAs(ss,tt,uu,BigX,BigY,Li2x,Li2y,Li3x,Li3y,Li4x,Li4y,
     & Li4z,Li4zinv)
      write(6,*) 'B routines '
      write(6,*) 'Bsx ',
     & Bsx(ss,tt,uu,BigX,BigY,BigS,Li2x,Li3x,Li2y,Li3y,Li4y,Li4z)
      write(6,*) 'Bs ',
     & AGTYBs(ss,tt,uu,Bigx,Bigy,Bigs,Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,
     & Li4z,Li4zinv)
      write(6,*) 'C routines '
      write(6,*) 'Csx ',
     & Csx(ss,tt,uu,Bigx,Bigy,Bigs,Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,Li4z)
      write(6,*) 'Cs ',
     & AGTYCs(ss,tt,uu,Bigx,Bigy,Bigs,Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,
     & Li4z,Li4zinv)
      write(6,*) 'D 1 routines'
      write(6,*) 'D1s ',AGTYD1s(ss,tt,uu,BigX,Bigy,Bigs,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z,Li4zinv)
      write(6,*) 'D1sx ',D1sx(ss,tt,uu,Bigx,Bigy,Bigs,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z)
      write(6,*) 'D 2 routines'
      write(6,*) 'D2s ',AGTYD2s(ss,tt,uu,BigX,Bigy,Bigs,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z,Li4zinv)
      write(6,*) 'D2sx ',D2sx(ss,tt,uu,Bigx,Bigy,Bigs,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z)
      write(6,*) 'E 1 routines'
      write(6,*) 'E1s ',AGTYE1s(tt,uu,Bigx,Bigy,Bigs)
      write(6,*) 'E1sx ',E1sx(tt,uu,Bigx,Bigy,Bigs)
      write(6,*) 'E 2 routines'
      write(6,*) 'E2s ',AGTYE2s(ss,tt,uu,Bigx,BigY,Bigs,Li2x,Li2y
     &     ,Li3x,Li3y)
      write(6,*) 'E2sx ',E2sx(ss,tt,uu,Bigx,Bigy,Bigs,Li2x,Li3x,Li3y)
      write(6,*) 'E 3 routines'
      write(6,*) 'E3s ',AGTYE3s(tt,uu,BigX,Bigy,Bigs)
      write(6,*) 'E3sx ',E3sx(tt,uu,Bigx,Bigy,Bigs)

      write(6,*) 'F 1 routines'
      write(6,*) 'F1s ',AGTYF1s(tt,uu,BigX,BigY,BigS)
      write(6,*) 'F1sx ',F1sx(tt,uu,Bigx,Bigy,Bigs)

      return
       end


!---- Mathematica Checked routines (name, formula)
!     As (A.1)
!     Bs
!     Cs
!     E1s
!     E2s

!!===== E4s is missing!
!==== F2s and u are missing!



!----- Ds needed fix to get (t<->u) in Fortran
