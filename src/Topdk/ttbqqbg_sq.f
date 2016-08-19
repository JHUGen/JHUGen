      double precision function ttbqqbg_sq(i1,i2,i9,i5,i3,i6,i8)
c--- returns the summed squared helicity amplitudes for the
c--- ttbqqbg amplitudes
      implicit none
      include 'constants.f'
      integer i1,i2,i9,i5,i3,i6,i8
      double complex sq,tq,qq,rq,
     . ttbqqbsqpp,ttbqqbsqpm,ttbqqbsqmp,ttbqqbsqmm,
     . ttbqqbtqpp,ttbqqbtqpm,ttbqqbtqmp,ttbqqbtqmm,
     . ttbqqbqqpp,ttbqqbqqpm,ttbqqbqqmp,ttbqqbqqmm,
     . ttbqqbrqpp,ttbqqbrqpm,ttbqqbrqmp,ttbqqbrqmm
      double precision appsq,apmsq,ampsq,ammsq
      
c-- for q-qb , there are four colour amplitudes:
c---    sq        proportional to Ta(it,i1)*delta(i2,ib)
c---    tq        proportional to Ta(i2,ib)*delta(it,i1)
c---    qq        proportional to Ta(i2,i1)*delta(it,ib)
c---    rq        proportional to Ta(it,ib)*delta(i2,i1)
      sq=ttbqqbsqpp(i1,i2,i9,i5,i3,i6,i8)
      tq=ttbqqbtqpp(i1,i2,i9,i5,i3,i6,i8)
      qq=ttbqqbqqpp(i1,i2,i9,i5,i3,i6,i8)
      rq=ttbqqbrqpp(i1,i2,i9,i5,i3,i6,i8)
      appsq=xn*cf*(
     . +xn*cdabs(sq)**2+xn*cdabs(tq)**2
     . +xn*cdabs(qq)**2+xn*cdabs(rq)**2
     . +2d0*Dble((sq+tq)*dconjg(rq+qq)))
      
      sq=ttbqqbsqpm(i1,i2,i9,i5,i3,i6,i8)
      tq=ttbqqbtqpm(i1,i2,i9,i5,i3,i6,i8)
      qq=ttbqqbqqpm(i1,i2,i9,i5,i3,i6,i8)
      rq=ttbqqbrqpm(i1,i2,i9,i5,i3,i6,i8)
      apmsq=xn*cf*(
     . +xn*cdabs(sq)**2+xn*cdabs(tq)**2
     . +xn*cdabs(qq)**2+xn*cdabs(rq)**2
     . +2d0*Dble((sq+tq)*dconjg(rq+qq)))
      
      sq=ttbqqbsqmp(i1,i2,i9,i5,i3,i6,i8)
      tq=ttbqqbtqmp(i1,i2,i9,i5,i3,i6,i8)
      qq=ttbqqbqqmp(i1,i2,i9,i5,i3,i6,i8)
      rq=ttbqqbrqmp(i1,i2,i9,i5,i3,i6,i8)
      ampsq=xn*cf*(
     . +xn*cdabs(sq)**2+xn*cdabs(tq)**2
     . +xn*cdabs(qq)**2+xn*cdabs(rq)**2
     . +2d0*Dble((sq+tq)*dconjg(rq+qq)))
      
      sq=ttbqqbsqmm(i1,i2,i9,i5,i3,i6,i8)
      tq=ttbqqbtqmm(i1,i2,i9,i5,i3,i6,i8)
      qq=ttbqqbqqmm(i1,i2,i9,i5,i3,i6,i8)
      rq=ttbqqbrqmm(i1,i2,i9,i5,i3,i6,i8)
      ammsq=xn*cf*(
     . +xn*cdabs(sq)**2+xn*cdabs(tq)**2
     . +xn*cdabs(qq)**2+xn*cdabs(rq)**2
     . +2d0*Dble((sq+tq)*dconjg(rq+qq)))
      
      ttbqqbg_sq=appsq+apmsq+ampsq+ammsq

      return
      end
      
