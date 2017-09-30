      subroutine qqb_higgs_odd(p,msq)
      implicit none
c---Matrix element squared averaged over initial colors and spins
c     f(-p1) + f(-p2) --> A + f(p5)
c                         |
c                         --> b(p3)+bbar(p4)
c   
c   Here A is pseudoscalar (CP-odd) Higgs  
c   Note that the approximate forms (for large mbsq) are exactly
c   the same as in the scalar case, apart from an overall factor                   
c--all momenta incoming

      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'hdecaymode.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),gg,qg,gq,qq,hdecay
      double precision ehsvm3_odd,ehsvm4_odd,origmbsq,s34,msqgamgam

      s34=(p(3,4)+p(4,4))**2
     &   -(p(3,1)+p(4,1))**2
     &   -(p(3,2)+p(4,2))**2
     &   -(p(3,3)+p(4,3))**2

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

      origmbsq=mbsq
      mbsq=mt**2
      gg=+avegg*ehsvm3_odd(s(1,2),s(1,5),s(2,5))*hdecay
      qq=+aveqq*ehsvm4_odd(s(1,2),s(1,5),s(2,5))*hdecay
      qg=-aveqg*ehsvm4_odd(s(1,5),s(1,2),s(2,5))*hdecay
      gq=-aveqg*ehsvm4_odd(s(2,5),s(1,5),s(1,2))*hdecay
      mbsq=origmbsq

      do j=-nf,nf    
      do k=-nf,nf
      msq(j,k)=0d0

      if ((j.eq. 0) .or. (k.eq.0)) then
           if ((j.eq. 0) .and. (k.eq.0)) then
                msq(j,k)=gg
           elseif ((j.eq.0).and.(k.ne.0)) then
                msq(j,k)=gq
           elseif ((j.ne.0).and.(k.eq.0)) then
                msq(j,k)=qg
           endif
      elseif ((j.eq.-k).and. (j.ne.0)) then
           msq(j,k)=qq
      endif

      enddo
      enddo

      return
      end
      
      
      double precision function ehsvm3_odd(s,t,u)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
c---Matrix element squared Eqn 2.2 of EHSV
      double complex ehsva2_odd,ehsva4_odd
      double precision s,t,u
      logical approx
      parameter(approx=.false.)
c--- approx TRUE uses the heavy fermion approximation to Msq

      if (approx) then
      ehsvm3_odd=gwsq/pi*as**3*xn*V/9d0*(
     .        hmass**8+s**4+t**4+u**4)/s/t/u/wmass**2
      else
      ehsvm3_odd=
     . abs(ehsva2_odd(s,t,u))**2+abs(ehsva2_odd(u,s,t))**2
     .+abs(ehsva2_odd(t,u,s))**2+abs(ehsva4_odd(s,t,u))**2 
      ehsvm3_odd=gwsq/pi*as**3*xn*V*hmass**8/(s*t*u*wmass**2)*ehsvm3_odd
      endif
c--- effective Lagrangian increased by factor of (3/2) wrt. scalar case,
c--- hence that factor squared in the rate
      ehsvm3_odd=ehsvm3_odd*(3d0/2d0)**2
      
      return
      end
      
      
      double precision function ehsvm4_odd(s,t,u)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
c---Matrix element squared Eqn 2.6 of EHSV
      double complex ehsva5_odd
      double precision s,t,u
      logical approx
      parameter(approx=.false.)
c--- approx TRUE uses the heavy fermion approximation to Msq

      if (approx) then
      ehsvm4_odd=gwsq/pi*as**3*V/18d0*(u**2+t**2)/s/wmass**2
      else
      ehsvm4_odd=abs(ehsva5_odd(s,t,u))**2
      ehsvm4_odd=gwsq/(4d0*pi)*as**3*V/2d0*(u**2+t**2)/(s*wmass**2)
     . *hmass**4/(u+t)**2*ehsvm4_odd
      endif
c--- effective Lagrangian increased by factor of (3/2) wrt. scalar case,
c--- hence that factor squared in the rate
      ehsvm4_odd=ehsvm4_odd*(3d0/2d0)**2
      
      return
      end
      
      
 
