      subroutine checkscheme(j1,j2,j3,j4)
      implicit none
      include 'types.f'
c--- this routine checks the (universal) transition rules between
c--- amplitudes in the dred and tH-V schemes

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'scheme.f'
      include 'deltar.f'
      include 'zprods_com.f'
      integer:: j1,j2,j3,j4,h1,h2,h3,h4
      complex(dp):: Alo(3,2,2,2,2),
     & AvirtDR(3,2,2,2,2),AvirtHV(3,2,2,2,2),
     & A0ab(2,2,2),A0ba(2,2,2),
     & A41abDR(2,2,2),A41baDR(2,2,2),A43abDR(2,2,2),A43baDR(2,2,2),
     & A41abHV(2,2,2),A41baHV(2,2,2),A43abHV(2,2,2),A43baHV(2,2,2),
     & Amploqarb(2,2),AmpvqarbDR(2,2),AmpvqarbHV(2,2),
     & AmpvqaqaDR(2,2),AmpvqaqaHV(2,2),
     & A0Hqarbmppm,A0Hqarbmpmp,A41Hqarbmppm,A41Hqarbmpmp,
     & A42Hqarbmppm,A42Hqarbmpmp
      real(dp):: renDR,renHV,H4prenorm,diff

c--- 4-GLUON CHECK

      call   Amplo(j1,j2,j3,j4,Alo)

c--- dimensional reduction scheme
      scheme='dred'
      deltar=0._dp
      call Ampvirt_gggg(j1,j2,j3,j4,AvirtDR)
      renDR=H4prenorm()

c--- 't Hooft-Veltman scheme
      scheme='tH-V'
      deltar=1._dp
      call Ampvirt_gggg(j1,j2,j3,j4,AvirtHV)
      renHV=H4prenorm()

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      diff=real((AvirtDR(1,h1,h2,h3,h4)+renDR*Alo(1,h1,h2,h3,h4)
     &         -(AvirtHV(1,h1,h2,h3,h4)+renHV*Alo(1,h1,h2,h3,h4)))
     &          /Alo(1,h1,h2,h3,h4))
     &          *xn
     &    -(4._dp*xn/6._dp) ! expected answer, units of gsq/(4*pi)^2
      write(6,'(a6,4i3,f14.9)') 'GGGG',h1,h2,h3,h4,diff
      enddo
      enddo
      enddo
      enddo


c--- 2-QUARK, 2-GLUON CHECK

      call   Amplo_AQgg(j1,j2,j3,j4,A0ab,A0ba)

c--- dimensional reduction scheme
      scheme='dred'
      deltar=0._dp
      call Ampvirt_AQgg(j1,j2,j3,j4,A41abDR,A41baDR,A43abDR,A43baDR)
      renDR=H4prenorm()

c--- 't Hooft-Veltman scheme
      scheme='tH-V'
      deltar=1._dp
      call Ampvirt_AQgg(j1,j2,j3,j4,A41abHV,A41baHV,A43abHV,A43baHV)
      renHV=H4prenorm()

      do h1=1,2
      do h2=1,2
      do h3=1,2
      diff=real((A41abDR(h1,h2,h3)+renDR*A0ab(h1,h2,h3)
     &         -(A41abHV(h1,h2,h3)+renHV*A0ab(h1,h2,h3)))
     &          /A0ab(h1,h2,h3))
     &          *xn
     &    -(Cf+xn/3._dp) ! expected answer, units of gsq/(4*pi)^2
      write(6,'(a6,4i3,f14.9)') 'AQGG',h1,h2,h3,h4,diff
      diff=real((A43abDR(h1,h2,h3)
     &         -(A43abHV(h1,h2,h3)))
     &          /A0ab(h1,h2,h3))
     &          *xn
     &    -zip         ! expected answer is zero
      write(6,'(a6,4i3,f14.9)') 'AQgg',h1,h2,h3,h4,diff
      enddo
      enddo
      enddo



c--- 4-QUARK (NON-IDENTICAL) CHECK

      amploqarb(1,2)=A0Hqarbmppm(j1,j2,j3,j4,za,zb)
      amploqarb(1,1)=A0Hqarbmpmp(j1,j2,j3,j4,za,zb)

c--- dimensional reduction scheme
      scheme='dred'
      deltar=0._dp
      ampvqarbDR(1,2) =A41Hqarbmppm(j1,j2,j3,j4,za,zb)
      ampvqarbDR(1,1) =A41Hqarbmpmp(j1,j2,j3,j4,za,zb)
      renDR=H4prenorm()

c--- 't Hooft-Veltman scheme
      scheme='tH-V'
      deltar=1._dp
      ampvqarbHV(1,2) =A41Hqarbmppm(j1,j2,j3,j4,za,zb)
      ampvqarbHV(1,1) =A41Hqarbmpmp(j1,j2,j3,j4,za,zb)
      renHV=H4prenorm()

      do h2=1,2
      diff=real((ampvqarbDR(1,h2)+renDR*amploqarb(1,h2)
     &         -(ampvqarbHV(1,h2)+renHV*amploqarb(1,h2)))
     &          /amploqarb(1,h2))
     &          *xn
     &    -(4._dp*Cf/2._dp) ! expected answer, units of gsq/(4*pi)^2
      write(6,'(a6,2i3,f14.9)') 'QARB',1,h2,diff
      enddo


c--- 4-QUARK (IDENTICAL) CHECK

      amploqarb(1,2)=A0Hqarbmppm(j1,j2,j3,j4,za,zb)
      amploqarb(1,1)=A0Hqarbmpmp(j1,j2,j3,j4,za,zb)

c--- dimensional reduction scheme
      scheme='dred'
      deltar=0._dp
      ampvqaqaDR(1,2) =A42Hqarbmppm(j1,j2,j3,j4,za,zb)
      ampvqaqaDR(1,1) =A42Hqarbmpmp(j1,j2,j3,j4,za,zb)
      renDR=H4prenorm()

c--- 't Hooft-Veltman scheme
      scheme='tH-V'
      deltar=1._dp
      ampvqaqaHV(1,2) =A42Hqarbmppm(j1,j2,j3,j4,za,zb)
      ampvqaqaHV(1,1) =A42Hqarbmpmp(j1,j2,j3,j4,za,zb)
      renHV=H4prenorm()

c--- now convert to A4;2 as defined in EGZ
c---     D & S : (N*d14*d23*A41+d12*d34*A42)
c---     EGZ   : (N*d14*d23-d12*d34)*A41+d12*d34*A42)
c--- so that A42(EGZ)=A42(DS)+A41
      ampvqaqaDR(1,2)=ampvqaqaDR(1,2)+ampvqarbDR(1,2)
      ampvqaqaDR(1,1)=ampvqaqaDR(1,1)+ampvqarbDR(1,1)
      ampvqaqaHV(1,2)=ampvqaqaHV(1,2)+ampvqarbHV(1,2)
      ampvqaqaHV(1,1)=ampvqaqaHV(1,1)+ampvqarbHV(1,1)

      do h2=1,2
      diff=real((ampvqaqaDR(1,h2)
     &         -(ampvqaqaHV(1,h2)))
     &          /amploqarb(1,h2)
     &          )
     &          *xn
     &    -zip          ! expected answer is zero
      write(6,'(a6,2i3,f14.9)') 'QAQA',1,h2,diff
      enddo


c      pause

      return
      end

