      subroutine qqb_QQb_g(p,msq)
      implicit none
      include 'types.f'


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
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'

      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: ttbqqbg,ttbggg
      real(dp):: wtqqb,wtqbq,wtqg,wtgq,wtqbg,wtgqb,wtgg

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
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
      if ((j == 0) .and. (k==0)) then
          msq(j,k)=avegg*wtgg

      elseif ((j > 0) .and. (j == -k)) then
          msq(j,k)=aveqq*wtqqb
      elseif ((j < 0) .and. (j == -k)) then
          msq(j,k)=aveqq*wtqbq

      elseif ((j < 0) .and. (k==0)) then
          msq(j,k)=aveqg*wtqbg
      elseif ((j > 0) .and. (k==0)) then
          msq(j,k)=aveqg*wtqg
      elseif ((j == 0) .and. (k<0)) then
          msq(j,k)=aveqg*wtgqb
      elseif ((j == 0) .and. (k>0)) then
          msq(j,k)=aveqg*wtgq
      endif
      enddo
      enddo

      return
      end

      function ttbqqbg(i1,i2,i3,i4,i5)
      implicit none
      include 'types.f'
      real(dp):: ttbqqbg

C      Q(P1) + Qbar(P2) --> q(-P3) + qb(-P4) + g(-P5)
C      Taken from:-
C      %\cite{Ellis:1986ef}
C      \bibitem{Ellis:1986ef}
C      R.~K.~Ellis and J.~C.~Sexton,
C      %``Explicit Formulae For Heavy Flavor Production,''
C      Nucl.\ Phys.\ B {\bf 282}, 642 (1987).
C      %%CITATION = NUPHA,B282,642;%%
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'breit.f'
      integer:: i1,i2,i3,i4,i5
      real(dp):: xm2,DL1,DL2,DL3,DL4,res1,res2
      real(dp):: P13,P14,P15,P23,P24,P25,P12,P35,P34,P45,S12

      xm2=mass2**2
      S12=2._dp*xm2+S(i1,i2)
      P12=+0.5_dp*s(i1,i2)

      P13=+0.5_dp*s(i1,i3)
      P14=+0.5_dp*s(i1,i4)
      P15=+0.5_dp*s(i1,i5)

      P23=+0.5_dp*s(i2,i3)
      P24=+0.5_dp*s(i2,i4)
      P25=+0.5_dp*s(i2,i5)

      P35=+0.5_dp*s(i3,i5)
      P34=+0.5_dp*s(i3,i4)
      P45=+0.5_dp*s(i4,i5)

      DL1 = P13/P25-2._dp*P35/s12
      DL2 = P14/P25-2._dp*P45/S12
      DL3 = P23/P15-2._dp*P35/S12
      DL4 = P24/P15-2._dp*P45/S12

      res1 =(P13**2+P23**2+P14**2+P24**2
     & +XM2*(P12+P34+XM2))/2._dp/S12/P34
     & *(4._dp*V**2/xn*(P13/P15/P35+P24/P25/P45)
     & +4._dp*V/xn*(2._dp*P14/P15/P45+2._dp*P23/P25/P35
     & -P13/P15/P35-P24/P25/P45-P12/P15/P25-P34/P35/P45) )
      res2 =
     & -V*(xn**2-4._dp)/xn
     & *2._dp*XM2/S12/P34*((P13-P14)/P25-(P23-P24)/P15)
     & +4._dp*V**2/xn*XM2*( (P35**2+P45**2)/P35/P45/S12**2
     & -0.5_dp*(1._dp/P15+1._dp/P25+1._dp/P35+1._dp/P45)/S12
     & -0.25_dp*(1._dp/P15+1._dp/P25+XM2/P15**2+XM2/P25**2+4._dp/S12)/P34
     & -(DL1**2+DL2**2+DL3**2+DL4**2)/4._dp/P34**2)
     & -2._dp*V/xn
     & *XM2/S12/P34*(1._dp+2._dp*P34/S12+XM2/P15+XM2/P25
     & +(P35**2+P45**2)/P15/P25
     & +(P13-P14)*(DL1-DL2)/P34+(P23-P24)*(DL3-DL4)/P34 )


      ttbqqbg=gsq**3*(res1+res2)

      return
      end

      function ttbggg(i1,i2,i3,i4,i5)
      implicit none
      include 'types.f'
      real(dp):: ttbggg
      real(dp):: Bsexton
      integer:: i1,i2,i3,i4,i5
      include 'qcdcouple.f'

       ttbggg=gsq**3*(
     & +Bsexton(i1,i2,i3,i4,i5)+Bsexton(i2,i1,i3,i4,i5)
     & +Bsexton(i1,i2,i4,i5,i3)+Bsexton(i2,i1,i4,i5,i3)
     & +Bsexton(i1,i2,i5,i3,i4)+Bsexton(i2,i1,i5,i3,i4)

     & +Bsexton(i1,i2,i5,i4,i3)+Bsexton(i2,i1,i5,i4,i3)
     & +Bsexton(i1,i2,i4,i3,i5)+Bsexton(i2,i1,i4,i3,i5)
     & +Bsexton(i1,i2,i3,i5,i4)+Bsexton(i2,i1,i3,i5,i4))

      return
      end

      function Bsexton(i1,i2,i3,i4,i5)
      implicit none
      include 'types.f'
      real(dp):: Bsexton
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

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'breit.f'
      include 'first.f'
      integer:: i1,i2,i3,i4,i5
      real(dp):: xm2,xm4,res
      real(dp):: P13,P14,P15,P23,P24,P25,P12,P35,P34,P45,S12

      if (first) then
      first=.false.
      write(6,*) 'Heavy Quark mass:',mass2
      endif


      xm2=mass2**2
      XM4=XM2**2
      S12=2._dp*xm2+S(i1,i2)
      P12=+0.5_dp*s(i1,i2)

      P13=+0.5_dp*s(i1,i3)
      P14=+0.5_dp*s(i1,i4)
      P15=+0.5_dp*s(i1,i5)

      P23=+0.5_dp*s(i2,i3)
      P24=+0.5_dp*s(i2,i4)
      P25=+0.5_dp*s(i2,i5)

      P35=+0.5_dp*s(i3,i5)
      P34=+0.5_dp*s(i3,i4)
      P45=+0.5_dp*s(i4,i5)
      RES=
     &   V*(xn**2+1._dp)/4._dp/xn**2
     &   *(P12*(P15**2+P25**2+2._dp*XM2*(P35+P45))/P13/P14/P23/P24 )
     & -V*(xn**2+1._dp)*XM2/xn**2
     & *(S12**2+2._dp*P13*P14+2._dp*P23*P24+4._dp*XM2*P12)/P13/P14/P23/P24/8._dp
     & -V*(
     & + P13*(P13**2+P23**2+2._dp*XM2*(P34+P35))/P34/P14/P15/P25
     & + (P15**2+P25**2+2._dp*XM2*(P35+P45))/2._dp/P34/P14/P23 )
     & +V*XM2
     & *( (P12-3._dp*XM2)/P34/P15/P25
     &   +(3._dp*P13*P23+5._dp*P13*P25-2._dp*P13**2)/P13/P14/P23/P24 )
     & +V*XM4
     & *( +(P13**2+P23**2+P15**2+P15*P25-S12*P13)/P34/P13/P24/P15/P25
     &    +(P15+P24)/P13/P23/P14/P25 )
      RES = RES
     & +2._dp*V*xn**2
     & *P23*(P13**2+P23**2+2._dp*XM2*(P34+P35))/S12/P34/P45/P25
     & +V*xn**2
     & *(P14*P24*(P14**2+P24**2+2._dp*XM2*(P34+P45))/S12/P34/P45/P13/P25)
      RES = RES
     & +2._dp*V*xn**2*XM2
     & *(4._dp*P15*P25+P13*P23+P14*P24
     & -0.125_dp*S12*(S12+2._dp*P34))/S12/P34/P15/P25
     & +V*xn**2*XM2
     & *((P34-2._dp*XM2)/S12/P13/P24-S12/4._dp/P34/P15/P25)
     & -2._dp*V*xn**2*XM4*(0.25_dp*S12**2+P34**2+P45**2)
     & /S12/P34/P45/P13/P25
      RES = RES
     & +4._dp*V*xn**2*XM2
     & *((P34**2+P35**2+P45**2)/S12**2/P34**2+2._dp*P23**2/S12/P34**2/P15)
     & +V**2*XM2*
     & ( (P13/P25-P23/P15)**2/P34**2+2._dp*XM2/P34/P15**2
     & +XM2*(P34+P45)/P34/P45/P13/P25 )
      RES = RES
     & +V**2/xn**2*XM2*0.25_dp
     & *( (P13**2+P23**2-(P13+P23)*P45)/P13/P23/P14/P25
     &  +2._dp*XM2*S12/P13/P23/P14/P25-4._dp*XM4/P13/P14/P25**2)
      RES = RES
     & +V**3/xn**2*XM2*0.5_dp
     & *((XM4+4._dp*XM2*P24-2._dp*P24*P25)/P13**2/P24**2
     & +(S12*P34+2._dp*P15*P25-P14*P24)/P13/P23/P14/P24)
     & +V/xn**2*XM2
     & *(XM4/P13/P23/P14/P25
     & -(5._dp*S12*P13+8._dp*P13*P23+6._dp*P13*P25-4._dp*P13**2)
     &  /4._dp/P13/P14/P23/P24)
      BSEXTON=RES
      RETURN
      end
