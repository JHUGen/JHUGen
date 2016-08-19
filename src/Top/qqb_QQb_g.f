      subroutine qqb_QQb_g(p,msq) 
      implicit none

************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2002.                                                  *
*     calculate the element squared
*     for the process                                                  *
c---   Notation                                                       *
C      Quark antiquark annihilation in order alfa_s^3
C      Q(P1) + Qbar(P2) --> q(-P3) + qbar(-P4) + g(-P5)
************************************************************************
      include 'constants.f'
      include 'sprods_com.f'
      
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision ttbqqbg,ttbggg
      double precision wtqqb,wtqbq,wtqg,wtgq,wtqbg,wtgqb,wtgg

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      call dotem(5,p,s)

      wtqqb=ttbqqbg(4,3,2,1,5)
      wtqbq=ttbqqbg(4,3,1,2,5)

      wtgq=-ttbqqbg(4,3,5,2,1)
      wtqg=-ttbqqbg(4,3,5,1,2)
      wtgqb=-ttbqqbg(4,3,2,5,1)
      wtqbg=-ttbqqbg(4,3,1,5,2)
      wtgg=ttbggg(4,3,1,2,5)

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      do k=-nf,nf
      if ((j .eq. 0) .and. (k.eq.0)) then
          msq(j,k)=avegg*wtgg

      elseif ((j .gt. 0) .and. (j .eq. -k)) then
          msq(j,k)=aveqq*wtqqb
      elseif ((j .lt. 0) .and. (j .eq. -k)) then
          msq(j,k)=aveqq*wtqbq

      elseif ((j .lt. 0) .and. (k.eq.0)) then
          msq(j,k)=aveqg*wtqbg
      elseif ((j .gt. 0) .and. (k.eq.0)) then
          msq(j,k)=aveqg*wtqg
      elseif ((j .eq. 0) .and. (k.lt.0)) then
          msq(j,k)=aveqg*wtgqb
      elseif ((j .eq. 0) .and. (k.gt.0)) then
          msq(j,k)=aveqg*wtgq
      endif
      enddo
      enddo

      return
      end

      double precision function ttbqqbg(i1,i2,i3,i4,i5)
      implicit none
C      Q(P1) + Qbar(P2) --> q(-P3) + qb(-P4) + g(-P5)
C      Taken from:-
C      %\cite{Ellis:1986ef}
C      \bibitem{Ellis:1986ef}
C      R.~K.~Ellis and J.~C.~Sexton,
C      %``Explicit Formulae For Heavy Flavor Production,''
C      Nucl.\ Phys.\ B {\bf 282}, 642 (1987).
C      %%CITATION = NUPHA,B282,642;%%
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'breit.f'
      integer i1,i2,i3,i4,i5
      double precision xm2,DL1,DL2,DL3,DL4,res1,res2
      double precision P13,P14,P15,P23,P24,P25,P12,P35,P34,P45,S12

      xm2=mass2**2
      S12=2d0*xm2+S(i1,i2)
      P12=+0.5d0*s(i1,i2)

      P13=+0.5d0*s(i1,i3)
      P14=+0.5d0*s(i1,i4)
      P15=+0.5d0*s(i1,i5)

      P23=+0.5d0*s(i2,i3)
      P24=+0.5d0*s(i2,i4)
      P25=+0.5d0*s(i2,i5)

      P35=+0.5d0*s(i3,i5)
      P34=+0.5d0*s(i3,i4)
      P45=+0.5d0*s(i4,i5)

      DL1 = P13/P25-2d0*P35/s12
      DL2 = P14/P25-2d0*P45/S12
      DL3 = P23/P15-2d0*P35/S12
      DL4 = P24/P15-2d0*P45/S12

      res1 =(P13**2+P23**2+P14**2+P24**2
     . +XM2*(P12+P34+XM2))/2D0/S12/P34
     . *(4D0*V**2/xn*(P13/P15/P35+P24/P25/P45)
     . +4D0*V/xn*(2D0*P14/P15/P45+2D0*P23/P25/P35
     . -P13/P15/P35-P24/P25/P45-P12/P15/P25-P34/P35/P45) )
      res2 = 
     . -V*(xn**2-4D0)/xn
     . *2D0*XM2/S12/P34*((P13-P14)/P25-(P23-P24)/P15)
     . +4D0*V**2/xn*XM2*( (P35**2+P45**2)/P35/P45/S12**2 
     . -0.5d0*(1D0/P15+1D0/P25+1D0/P35+1D0/P45)/S12
     . -0.25d0*(1D0/P15+1D0/P25+XM2/P15**2+XM2/P25**2+4D0/S12)/P34
     . -(DL1**2+DL2**2+DL3**2+DL4**2)/4D0/P34**2)
     . -2D0*V/xn
     . *XM2/S12/P34*(1D0+2D0*P34/S12+XM2/P15+XM2/P25
     . +(P35**2+P45**2)/P15/P25     
     . +(P13-P14)*(DL1-DL2)/P34+(P23-P24)*(DL3-DL4)/P34 )     


      ttbqqbg=gsq**3*(res1+res2)

      return
      end

      double precision function ttbggg(i1,i2,i3,i4,i5)
      double precision Bsexton
      integer i1,i2,i3,i4,i5
      include 'qcdcouple.f'

       ttbggg=gsq**3*(
     . +Bsexton(i1,i2,i3,i4,i5)+Bsexton(i2,i1,i3,i4,i5)
     . +Bsexton(i1,i2,i4,i5,i3)+Bsexton(i2,i1,i4,i5,i3)
     . +Bsexton(i1,i2,i5,i3,i4)+Bsexton(i2,i1,i5,i3,i4)

     . +Bsexton(i1,i2,i5,i4,i3)+Bsexton(i2,i1,i5,i4,i3)
     . +Bsexton(i1,i2,i4,i3,i5)+Bsexton(i2,i1,i4,i3,i5)
     . +Bsexton(i1,i2,i3,i5,i4)+Bsexton(i2,i1,i3,i5,i4))

      return
      end

      double precision function Bsexton(i1,i2,i3,i4,i5)
C      This is the four dimensional result for
C      Quark antiquark annihilation into three gluons in order alfa_s^3
C      Q(P1) + Qbar(P2) --> g(-P3) + g(-P4) + g(-P5)
C      Result requires symmetrisation.   
C      No average over initial spins and colours
C      Taken from:-
C      %\cite{Ellis:1986ef}
C      \bibitem{Ellis:1986ef}
C      R.~K.~Ellis and J.~C.~Sexton,
C      %``Explicit Formulae For Heavy Flavor Production,''
C      Nucl.\ Phys.\ B {\bf 282}, 642 (1987).
C      %%CITATION = NUPHA,B282,642;%%
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'breit.f'
      include 'first.f'
      integer i1,i2,i3,i4,i5
      double precision xm2,xm4,res
      double precision P13,P14,P15,P23,P24,P25,P12,P35,P34,P45,S12

      if (first) then
      first=.false.
      write(6,*) 'Heavy Quark mass:',mass2
      endif 


      xm2=mass2**2
      XM4=XM2**2
      S12=2d0*xm2+S(i1,i2)
      P12=+0.5d0*s(i1,i2)

      P13=+0.5d0*s(i1,i3)
      P14=+0.5d0*s(i1,i4)
      P15=+0.5d0*s(i1,i5)

      P23=+0.5d0*s(i2,i3)
      P24=+0.5d0*s(i2,i4)
      P25=+0.5d0*s(i2,i5)

      P35=+0.5d0*s(i3,i5)
      P34=+0.5d0*s(i3,i4)
      P45=+0.5d0*s(i4,i5)
      RES=
     .   V*(xn**2+1D0)/4D0/xn**2
     .   *(P12*(P15**2+P25**2+2D0*XM2*(P35+P45))/P13/P14/P23/P24 )
     . -V*(xn**2+1D0)*XM2/xn**2
     . *(S12**2+2D0*P13*P14+2D0*P23*P24+4D0*XM2*P12)/P13/P14/P23/P24/8D0 
     . -V*(
     . + P13*(P13**2+P23**2+2D0*XM2*(P34+P35))/P34/P14/P15/P25
     . + (P15**2+P25**2+2D0*XM2*(P35+P45))/2D0/P34/P14/P23 ) 
     . +V*XM2
     . *( (P12-3D0*XM2)/P34/P15/P25
     .   +(3D0*P13*P23+5D0*P13*P25-2D0*P13**2)/P13/P14/P23/P24 )
     . +V*XM4
     . *( +(P13**2+P23**2+P15**2+P15*P25-S12*P13)/P34/P13/P24/P15/P25
     .    +(P15+P24)/P13/P23/P14/P25 )
      RES = RES
     . +2D0*V*xn**2 
     . *P23*(P13**2+P23**2+2D0*XM2*(P34+P35))/S12/P34/P45/P25
     . +V*xn**2
     . *(P14*P24*(P14**2+P24**2+2D0*XM2*(P34+P45))/S12/P34/P45/P13/P25)
      RES = RES
     . +2D0*V*xn**2*XM2
     . *(4D0*P15*P25+P13*P23+P14*P24
     . -0.125d0*S12*(S12+2D0*P34))/S12/P34/P15/P25
     . +V*xn**2*XM2 
     . *((P34-2D0*XM2)/S12/P13/P24-S12/4D0/P34/P15/P25)
     . -2D0*V*xn**2*XM4*(0.25d0*S12**2+P34**2+P45**2)
     . /S12/P34/P45/P13/P25
      RES = RES
     . +4D0*V*xn**2*XM2
     . *((P34**2+P35**2+P45**2)/S12**2/P34**2+2D0*P23**2/S12/P34**2/P15)
     . +V**2*XM2*
     . ( (P13/P25-P23/P15)**2/P34**2+2D0*XM2/P34/P15**2
     . +XM2*(P34+P45)/P34/P45/P13/P25 )
      RES = RES
     . +V**2/xn**2*XM2*0.25d0
     . *( (P13**2+P23**2-(P13+P23)*P45)/P13/P23/P14/P25
     .  +2D0*XM2*S12/P13/P23/P14/P25-4D0*XM4/P13/P14/P25**2)
      RES = RES
     . +V**3/xn**2*XM2*0.5D0
     . *((XM4+4D0*XM2*P24-2D0*P24*P25)/P13**2/P24**2
     . +(S12*P34+2D0*P15*P25-P14*P24)/P13/P23/P14/P24)
     . +V/xn**2*XM2 
     . *(XM4/P13/P23/P14/P25
     . -(5D0*S12*P13+8D0*P13*P23+6D0*P13*P25-4D0*P13**2)
     .  /4D0/P13/P14/P23/P24)
      BSEXTON=RES
      RETURN
      end 
