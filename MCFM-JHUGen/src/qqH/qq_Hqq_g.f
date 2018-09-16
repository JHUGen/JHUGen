      subroutine qq_Hqq_g(p,msq)
      implicit none
      include 'types.f'

c---Matrix element squared averaged over initial colors and spins
c
c     q(-p1)+q(-p2) -->  H(p3,p4)+q(p5)+q(p6)+g(p7)
c                           |
c                           |
c                           |
c                           ---> b(p3)+bbar(p4)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'hdecaymode.f'
      integer:: j,k,m,n
      real(dp):: p(mxpart,4),facqq,facqg,s34
      real(dp):: msq(-nf:nf,-nf:nf),hdecay,
     & msqll,msqlr,msqzzin,msqwzin,msqwl,
     & msxll,msxlr,msxzzin,msxwzin,msxwl,msqgamgam
      real(dp):: msqx(fn:nf,fn:nf,fn:nf,fn:nf)
      common/msq_all/msqx
      logical:: includeall

      integer,parameter::pn(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

c--- This flag decides whether or not to include all types of diagrams:
c---  FALSE --> only diagrams that look like WBF
c---  TRUE --> extra "associated" diagrams e.g. qqb -> W --> W(->qqb)H
c--- The comparison with Madgraph verified by JMC on 2/18/03, as follows:
c---    FALSE --> exact agreement
c---    TRUE ---> untested, no point in including
      includeall=.false.

      call dotem(7,p,s)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      do m=-nf,nf
      do n=-nf,nf
      msqx(j,k,m,n)=0._dp
      enddo
      enddo
      enddo
      enddo

      s34=(p(3,4)+p(4,4))**2
     & -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(hmass)
      else
      write(6,*) 'Unimplemented process in gg_hgg_v'
      stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      facqq=0.25_dp*gsq*Cf*gwsq**3*hdecay
      facqg=-facqq*(aveqg/aveqq)
C Extra factor of gsq*Cf compared to lowest order

C q-q and qbar-qbar
      call msq_gpieces(1,2,5,6,7,msqll,msqlr,msqzzin,msqwzin,msqwl)
      call msq_gpieces(1,2,6,5,7,msxll,msxlr,msxzzin,msxwzin,msxwl)

      do j=1,nf
      do k=1,nf
      do m=1,nf
      do n=1,nf

      if ((j==m) .and. (k==n)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqq*(
     & +msqll*((L(j)*L(k))**2+(R(j)*R(k))**2)
     & +msqlr*((L(j)*R(k))**2+(R(j)*L(k))**2))
      endif
      if ((j==n) .and. (k==m)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqq*(
     & +msxll*((L(j)*L(k))**2+(R(j)*R(k))**2)
     & +msxlr*((L(j)*R(k))**2+(R(j)*L(k))**2))
      endif
      if ((j==k) .and. (j==m) .and. (j==n)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqq*(
     & +msqzzin*((L(j)*L(k))**2+(R(j)*R(k))**2))
      endif
      if (  (pn(j)+pn(m) == +3) .and. (pn(k)+pn(n) == +3)
     ..and. (pn(j)+pn(k) == +3) ) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqq*Vsq(j,-m)*Vsq(k,-n)*(
     & +msxwzin*L(j)*L(k)+msxwl)
      endif
      if (  (pn(j)+pn(n) == +3) .and. (pn(k)+pn(m) == +3)
     ..and. (pn(j)+pn(k) == +3) ) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqq*Vsq(j,-n)*Vsq(k,-m)*(
     & +msqwzin*L(j)*L(k)+msqwl)
      endif
      if (j == k) then
        msqx(j,k,m,n)=msqx(j,k,m,n)/2._dp
      endif

      msqx(-j,-k,-m,-n)=msqx(j,k,m,n)

      enddo
      enddo
      enddo
      enddo

C q-qb
      call msq_gpieces(1,6,5,2,7,msqll,msqlr,msqzzin,msqwzin,msqwl)
      call msq_gpieces(1,6,2,5,7,msxll,msxlr,msxzzin,msxwzin,msxwl)

      do j=1,nf
      do k=-nf,-1
      do m=1,nf
      do n=-nf,-1

      if ((j==m) .and. (k==n)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqq*(
     & +msqll*((L(j)*L(-k))**2+(R(j)*R(-k))**2)
     & +msqlr*((L(j)*R(-k))**2+(R(j)*L(-k))**2))
       if (includeall) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqq*Vsq(j,k)*Vsq(m,n)*(
     & +msqwzin*L(j)*L(-k)+msqwl)
         if (j == -k) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqq*(
     & +msxll*((L(j)*L(j))**2+(R(j)*R(j))**2)
     & +msxlr*((L(j)*R(j))**2+(R(j)*L(j))**2)
     & +msqzzin*((L(j)*L(j))**2+(R(j)*R(j))**2))
         endif
       endif
      endif

      if (   (pn(j)+pn(m) == +3) .and. (pn(k)+pn(n) == -3)
     & .and. (pn(j)+pn(k) == pn(m)+pn(n)) ) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqq*Vsq(j,-m)*Vsq(k,-n)*msxwl
       if (includeall) then
         if ((j == -k) .and. (m == -n)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqq*(
     & +msxll*((L(j)*L(m))**2+(R(j)*R(m))**2)
     & +msxlr*((L(j)*R(m))**2+(R(j)*L(m))**2)
     & +msxwzin*L(j)*L(m)*Vsq(j,-m)*Vsq(k,-n)
     & )
         endif
       endif
      endif

      msqx(-j,-k,-m,-n)=msqx(j,k,m,n)

      enddo
      enddo
      enddo
      enddo

C qb-q
      call msq_gpieces(6,2,1,5,7,msqll,msqlr,msqzzin,msqwzin,msqwl)
      call msq_gpieces(6,2,5,1,7,msxll,msxlr,msxzzin,msxwzin,msxwl)

      do j=-nf,1
      do k=1,nf
      do m=1,nf
      do n=-nf,-1

      if ((j==n) .and. (k==m)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqq*(
     & +msqll*((L(-j)*L(k))**2+(R(-j)*R(k))**2)
     & +msqlr*((L(-j)*R(k))**2+(R(-j)*L(k))**2))
       if (includeall) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqq*Vsq(j,k)*Vsq(m,n)*(
     & +msqwzin*L(-j)*L(k)+msqwl)
         if (j == -k) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqq*(
     & +msxll*((L(-j)*L(-j))**2+(R(-j)*R(-j))**2)
     & +msxlr*((L(-j)*R(-j))**2+(R(-j)*L(-j))**2)
     & +msqzzin*((L(-j)*L(-j))**2+(R(-j)*R(-j))**2))
         endif
       endif
      endif

      if (   (pn(j)+pn(n) == -3) .and. (pn(k)+pn(m) == +3)
     & .and. (pn(j)+pn(k) == pn(m)+pn(n)) ) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqq*Vsq(j,-n)*Vsq(k,-m)*msxwl
       if (includeall) then
         if ((j == -k) .and. (m == -n)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqq*(
     & +msxll*((L(-j)*L(m))**2+(R(-j)*R(m))**2)
     & +msxlr*((L(-j)*R(m))**2+(R(-j)*L(m))**2)
     & +msxwzin*L(-j)*L(m)*Vsq(j,-n)*Vsq(k,-m)
     & )
         endif
       endif
      endif

      msqx(-j,-k,-m,-n)=msqx(j,k,m,n)

      enddo
      enddo
      enddo
      enddo

C q-g and qbar-g
      call msq_gpieces(1,7,5,6,2,msqll,msqlr,msqzzin,msqwzin,msqwl)
      call msq_gpieces(1,7,6,5,2,msxll,msxlr,msxzzin,msxwzin,msxwl)

      k=0
      do j=1,nf
      do m=1,nf
      do n=1,nf

      if ((j==m) .and. (n>0)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqg*(
     & +msqll*((L(j)*L(n))**2+(R(j)*R(n))**2)
     & +msqlr*((L(j)*R(n))**2+(R(j)*L(n))**2))
      endif
      if ((j==n) .and. (m>0)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqg*(
     & +msxll*((L(j)*L(m))**2+(R(j)*R(m))**2)
     & +msxlr*((L(j)*R(m))**2+(R(j)*L(m))**2))
      endif
      if ((j==m) .and. (j==n) .and. (m>0)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqg*(
     & +msqzzin*((L(j)*L(m))**2+(R(j)*R(m))**2))
      endif
      if (  (pn(j)+pn(m) == +3) .and. (pn(m)+pn(n) == +3)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqg*Vsq(j,-m)*Vsum(n)*(
     & +msxwzin*L(j)*L(m)+msxwl)
      endif
      if (  (pn(j)+pn(n) == +3) .and. (pn(m)+pn(n) == +3)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqg*Vsq(j,-n)*Vsum(m)*(
     & +msqwzin*L(j)*L(n)+msqwl)
      endif
      if (m == n) then
        msqx(j,k,m,n)=msqx(j,k,m,n)/2._dp
      endif

      msqx(-j,-k,-m,-n)=msqx(j,k,m,n)

      enddo
      enddo
      enddo

C g-q and g-qbar
      call msq_gpieces(7,2,5,6,1,msqll,msqlr,msqzzin,msqwzin,msqwl)
      call msq_gpieces(7,2,6,5,1,msxll,msxlr,msxzzin,msxwzin,msxwl)

      j=0
      do k=1,nf
      do m=1,nf
      do n=1,nf

      if ((k==n) .and. (m>0)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqg*(
     & +msqll*((L(k)*L(m))**2+(R(k)*R(m))**2)
     & +msqlr*((L(k)*R(m))**2+(R(k)*L(m))**2))
      endif
      if ((k==m) .and. (n>0)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqg*(
     & +msxll*((L(k)*L(n))**2+(R(k)*R(n))**2)
     & +msxlr*((L(k)*R(n))**2+(R(k)*L(n))**2))
      endif
      if ((k==m) .and. (k==n) .and. (m>0)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqg*(
     & +msqzzin*((L(k)*L(m))**2+(R(k)*R(m))**2))
      endif
      if (  (pn(k)+pn(n) == +3) .and. (pn(m)+pn(n) == +3)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqg*Vsq(k,-n)*Vsum(m)*(
     & +msxwzin*L(k)*L(n)+msxwl)
      endif
      if (  (pn(k)+pn(m) == +3) .and. (pn(m)+pn(n) == +3)) then
      msqx(j,k,m,n)=msqx(j,k,m,n)+facqg*Vsq(k,-m)*Vsum(n)*(
     & +msqwzin*L(k)*L(m)+msqwl)
      endif
      if (m == n) then
        msqx(j,k,m,n)=msqx(j,k,m,n)/2._dp
      endif

      msqx(-j,-k,-m,-n)=msqx(j,k,m,n)

      enddo
      enddo
      enddo

      do j=-nf,nf
      do k=-nf,nf
      do m=-nf,nf
      do n=-nf,nf
c--- add in to total for msq(j,k)
c        if (m>=n) then
c          msq(j,k) = msq(j,k)+msqx(j,k,m,n)
c        endif
c        if     ((j > 0) .and. (k < 0)) then
c          if (m>=n) msq(j,k) = msq(j,k)+msqx(j,k,m,n)
c        elseif ((j < 0) .and. (k > 0)) then
c          if (m<=n) msq(j,k) = msq(j,k)+msqx(j,k,m,n)
c        elseif ((j > 0) .and. (k > 0)) then
c          if ((pn(j)==pn(n)) .and. (pn(k)==pn(m)))
c     &                msq(j,k) = msq(j,k)+msqx(j,k,m,n)
c        elseif ((j < 0) .and. (k < 0)) then
c          if ((pn(j)==pn(n)) .and. (pn(k)==pn(m)))
c     &                msq(j,k) = msq(j,k)+msqx(j,k,m,n)
c        else
c          if (m>=n) msq(j,k) = msq(j,k)+msqx(j,k,m,n)
c        endif
        if ((j==m) .and. (k == n)) then
          msq(j,k) = msq(j,k)+msqx(j,k,m,n)
        endif
        if ((j==0) .and. (m*k > 0) .and. (n*k < 0)) then
          msq(j,k) = msq(j,k)+msqx(j,k,m,n)
        endif
        if ((k==0) .and. (m*j > 0) .and. (n*j > 0)) then
          msq(j,k) = msq(j,k)+msqx(j,k,m,n)
        endif
      enddo
      enddo
      enddo
      enddo

      return
      end


      subroutine msq_gpieces(i1,i2,i5,i6,i7,zll,zLR,zzLL,wzLL,wll)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      real(dp):: zll,zLR,zzLL,wzll,wll,htheta
      real(dp):: msqll_a,msqlr_a,msqinterf_a
      real(dp):: msqll_b,msqlr_b,msqinterf_b,msqll_c,msqll_d
      real(dp):: msqinterfx_a,msqinterfx_b
      real(dp):: propw,propz,x
      real(dp):: msq_gsamehel,msq_gopphel,msq_ginterf,msq_ginterfx
      real(dp):: s15,s16,s25,s26,s157,s167,s257,s267
      integer:: i1,i2,i5,i6,i7
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      htheta(x)=half+sign(half,x)
      propw(s15)=sign(one,(s15-wmass**2))*sqrt(
     .((s15-wmass**2)**2+htheta(s15)*(wmass*wwidth)**2)/wmass)
      propz(s15)=sign(one,(s15-zmass**2))
     & *sqrt(sqrt(1._dp-xw)/xw/2._dp/zmass
     & *((s15-zmass**2)**2+htheta(s15)*(zmass*zwidth)**2))

c--- Calculate some invariants
      s15=s(i1,i5)
      s16=s(i1,i6)
      s25=s(i2,i5)
      s26=s(i2,i6)
      s157=s15+s(i1,i7)+s(i5,i7)
      s167=s16+s(i1,i7)+s(i6,i7)
      s257=s25+s(i2,i7)+s(i5,i7)
      s267=s26+s(i2,i7)+s(i6,i7)

c--- Note that we have to add up both the radiation from the (15) line
c--- and that from the (26) line. Colour prohibits any interference -
c--- except in the identical case (msq_interfx below).
c--- Obtain (26) radiation from (15) by 1<->2, 5<->6.
      msqll_a=msq_gsamehel(i1,i2,i5,i6,i7)
      msqll_b=msq_gsamehel(i2,i1,i6,i5,i7)
      msqll_c=msq_gsamehel(i1,i2,i6,i5,i7)
      msqll_d=msq_gsamehel(i2,i1,i5,i6,i7)
      msqlr_a=msq_gopphel(i1,i2,i5,i6,i7)
      msqlr_b=msq_gopphel(i2,i1,i6,i5,i7)
      msqinterf_a=msq_ginterf(i1,i2,i5,i6,i7)
      msqinterf_b=msq_ginterf(i2,i1,i6,i5,i7)
      msqinterfx_a=msq_ginterfx(i1,i2,i5,i6,i7)
      msqinterfx_b=msq_ginterfx(i2,i1,i6,i5,i7)

c--- catch the unwanted diagrams for the gluon-quark processes
      if (i7 == 1) then
        msqll_b=0._dp
        msqlr_b=0._dp
        msqinterf_b=0._dp
        msqinterfx_b=0._dp
        msqinterfx_a=0._dp
        msqll_d=0._dp
      endif
      if (i7 == 2) then
        msqll_a=0._dp
        msqlr_a=0._dp
        msqinterf_a=0._dp
        msqinterfx_a=0._dp
        msqinterfx_b=0._dp
        msqll_c=0._dp
      endif

      zll=msqll_a/(propz(s157)*propz(s26))**2
     &   +msqll_b/(propz(s267)*propz(s15))**2

      zlr=msqlr_a/(propz(s157)*propz(s26))**2
     &   +msqlr_b/(propz(s267)*propz(s15))**2

      zzll=-2._dp/xn*msqinterf_a/(propz(s157)*propz(s26))
     &                        /(propz(s167)*propz(s25))
     &     -2._dp/xn*msqinterf_b/(propz(s267)*propz(s15))
     &                        /(propz(s257)*propz(s16))
     &     -2._dp/xn*msqinterfx_a/(propz(s167)*propz(s25))
     &                         /(propz(s267)*propz(s15))
     &     -2._dp/xn*msqinterfx_b/(propz(s157)*propz(s26))
     &                         /(propz(s257)*propz(s16))

      wzll=-2._dp/xn*msqinterf_a/(propz(s157)*propz(s26))
     &                        /(propw(s167)*propw(s25))
     &     -2._dp/xn*msqinterf_b/(propz(s267)*propz(s15))
     &                        /(propw(s257)*propw(s16))
     &     -2._dp/xn*msqinterfx_a/(propw(s167)*propw(s25))
     &                         /(propz(s267)*propz(s15))
     &     -2._dp/xn*msqinterfx_b/(propz(s157)*propz(s26))
     &                         /(propw(s257)*propw(s16))

      wll=msqll_c/(propw(s167)*propw(s25))**2
     &   +msqll_d/(propw(s257)*propw(s16))**2

      return
      end

      function msq_gsamehel(j1,j2,j5,j6,j7)
      implicit none
      include 'types.f'
      real(dp):: msq_gsamehel

      integer:: j1,j2,j5,j6,j7
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'

      msq_gsamehel=
     &    + 4._dp*s(j1,j2)*s(j1,j5)/s(j1,j7)*s(j5,j6)/s(j5,j7)
     &    + 2._dp*s(j1,j2)*s(j1,j5)/s(j1,j7)/s(j5,j7)*s(j6,j7)
     &    - 2._dp*s(j1,j2)*s(j1,j6)/s(j1,j7)
     &    + 2._dp*s(j1,j2)/s(j1,j7)*s(j5,j6)
     &    + 2._dp*s(j1,j2)*s(j5,j6)/s(j5,j7)
     &    + 2._dp*s(j1,j2)/s(j5,j7)*s(j6,j7)
     &    + 2._dp*s(j1,j5)/s(j1,j7)*s(j2,j7)*s(j5,j6)/s(j5,j7)
     &    + 2._dp/s(j1,j7)*s(j2,j7)*s(j5,j6)
     &    - 2._dp*s(j2,j5)*s(j5,j6)/s(j5,j7)

      return
      end

      function msq_gopphel(j1,j2,j5,j6,j7)
      implicit none
      include 'types.f'
      real(dp):: msq_gopphel

      integer:: j1,j2,j5,j6,j7
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'

      msq_gopphel=
     &    - 2._dp*s(j1,j2)*s(j1,j6)/s(j1,j7)
     &    + 4._dp*s(j1,j5)*s(j1,j6)/s(j1,j7)*s(j2,j5)/s(j5,j7)
     &    + 2._dp*s(j1,j5)*s(j1,j6)/s(j1,j7)*s(j2,j7)/s(j5,j7)
     &    + 2._dp*s(j1,j5)/s(j1,j7)*s(j2,j5)/s(j5,j7)*s(j6,j7)
     &    + 2._dp*s(j1,j6)/s(j1,j7)*s(j2,j5)
     &    + 2._dp*s(j1,j6)*s(j2,j5)/s(j5,j7)
     &    + 2._dp*s(j1,j6)*s(j2,j7)/s(j5,j7)
     &    + 2._dp/s(j1,j7)*s(j2,j5)*s(j6,j7)
     &    - 2._dp*s(j2,j5)*s(j5,j6)/s(j5,j7)

      return
      end

      function msq_ginterf(j1,j2,j5,j6,j7)
      implicit none
      include 'types.f'
      real(dp):: msq_ginterf

      integer:: j1,j2,j5,j6,j7
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'

      msq_ginterf=
     &    - 2._dp*s(j1,j2)*s(j1,j5)/s(j1,j7)*s(j5,j6)/s(j5,j7)
     &    - s(j1,j2)*s(j1,j5)/s(j1,j7)/s(j5,j7)*s(j6,j7)
     &    + s(j1,j2)*s(j1,j5)/s(j1,j7)
     &    - 2._dp*s(j1,j2)*s(j1,j6)/s(j1,j7)*s(j5,j6)/s(j6,j7)
     &    - s(j1,j2)*s(j1,j6)/s(j1,j7)*s(j5,j7)/s(j6,j7)
     &    + s(j1,j2)*s(j1,j6)/s(j1,j7)
     &    - 2._dp*s(j1,j2)/s(j1,j7)*s(j5,j6)
     &    + s(j1,j2)*s(j5,j6)/s(j5,j7)
     &    + s(j1,j2)*s(j5,j6)/s(j6,j7)
     &    + 2._dp*s(j1,j2)*s(j5,j6)**2/s(j5,j7)/s(j6,j7)
     &    + 2._dp*s(j1,j2)
     &    - s(j1,j5)/s(j1,j7)*s(j2,j7)*s(j5,j6)/s(j5,j7)
     &    - s(j1,j6)/s(j1,j7)*s(j2,j7)*s(j5,j6)/s(j6,j7)
     &    - 2._dp/s(j1,j7)*s(j2,j7)*s(j5,j6)
     &    + s(j2,j5)*s(j5,j6)/s(j5,j7)
     &    + s(j2,j6)*s(j5,j6)/s(j6,j7)

      return
      end

      function msq_ginterfx(j1,j2,j5,j6,j7)
      implicit none
      include 'types.f'
      real(dp):: msq_ginterfx

      integer:: j1,j2,j5,j6,j7
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'

      msq_ginterfx=
     &    + s(j1,j2)*s(j1,j5)/s(j1,j7)
     &    - 2._dp*s(j1,j2)*s(j1,j6)/s(j1,j7)*s(j5,j6)/s(j6,j7)
     &    - s(j1,j2)*s(j1,j6)/s(j1,j7)*s(j5,j7)/s(j6,j7)
     &    + s(j1,j2)/s(j1,j7)*s(j5,j6)
     &    + s(j1,j2)*s(j2,j5)/s(j2,j7)
     &    - 2._dp*s(j1,j2)*s(j2,j6)/s(j2,j7)*s(j5,j6)/s(j6,j7)
     &    - s(j1,j2)*s(j2,j6)/s(j2,j7)*s(j5,j7)/s(j6,j7)
     &    + s(j1,j2)/s(j2,j7)*s(j5,j6)
     &    - 2._dp*s(j1,j2)*s(j5,j6)/s(j6,j7)
     &    - 2._dp*s(j1,j2)*s(j5,j7)/s(j6,j7)
     &    + 2._dp*s(j1,j2)**2/s(j1,j7)/s(j2,j7)*s(j5,j6)
     &    - s(j1,j6)/s(j1,j7)*s(j2,j7)*s(j5,j6)/s(j6,j7)
     &    + s(j1,j6)*s(j5,j6)/s(j6,j7)
     &    - s(j1,j7)*s(j2,j6)/s(j2,j7)*s(j5,j6)/s(j6,j7)
     &    + s(j2,j6)*s(j5,j6)/s(j6,j7)
     &    + 2._dp*s(j5,j6)

      return
      end

