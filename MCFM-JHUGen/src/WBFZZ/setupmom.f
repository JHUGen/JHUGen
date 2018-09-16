      subroutine setupmom(p)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      real(dp):: al,be,ro,si,th,ph,p(mxpart,4),E1,E2,E7,y

C----setup p
      al=1.04225677d0
      be=0.5564225677d0
      ro=0.86454322d0
      si=1.04225677d0
      th=0.63425677d0
      ph=0.145563425677d0
      zmass=91.188d0
      p(:,:)=0d0

      E1=343.43678d0
      E2=436.68843678d0
      y=0.5d0*(1d0+cos(th))
      E7=-(zmass**2-zmass*E1-zmass*E2+E1*E2)/(zmass-E1*y-E2*(1d0-y))

      p(1,4)=-E1
      p(1,1)=0d0
      p(1,2)=0d0
      p(1,3)=+E1

      p(2,4)=-E2
      p(2,1)=0d0
      p(2,2)=0d0
      p(2,3)=-E2

      p(3,4)=+0.5d0*zmass
      p(3,1)=+0.5d0*zmass*sin(al)*cos(ro)
      p(3,2)=+0.5d0*zmass*sin(al)*sin(ro)
      p(3,3)=+0.5d0*zmass*cos(al)

      p(4,4)=+0.5d0*zmass
      p(4,1)=-0.5d0*zmass*sin(al)*cos(ro)
      p(4,2)=-0.5d0*zmass*sin(al)*sin(ro)
      p(4,3)=-0.5d0*zmass*cos(al)

      p(5,4)=+0.5d0*zmass
      p(5,1)=+0.5d0*zmass*sin(be)*cos(si)
      p(5,2)=+0.5d0*zmass*sin(be)*sin(si)
      p(5,3)=+0.5d0*zmass*cos(be)

      p(6,4)=+0.5d0*zmass
      p(6,1)=-0.5d0*zmass*sin(be)*cos(si)
      p(6,2)=-0.5d0*zmass*sin(be)*sin(si)
      p(6,3)=-0.5d0*zmass*cos(be)

      p(7,4)=E7
      p(7,1)=E7*sin(th)*sin(ph)
      p(7,2)=E7*sin(th)*cos(ph)
      p(7,3)=E7*cos(th)

      p(8,:)=-p(1,:)-p(2,:)-p(3,:)-p(4,:)-p(5,:)-p(6,:)-p(7,:)

      return
      end
