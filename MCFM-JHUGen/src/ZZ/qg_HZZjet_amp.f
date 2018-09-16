      subroutine qg_HZZjet_amp(i1,i2,i7,za,zb,amp_t,amp_b) 
      implicit none
      include 'types.f'
       
!===== C. Williams October 2013
!=== this is the amplitude for q(qb) g => H  + q(qb) needed for the 
!=== ZZ interference study 
!=== particle ordering 
!====== 0 --> q(i1)^- qb(i2)^+ g(i7) H( Z(34) Z(56) ) 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'masses.f' 
      include 'zprods_decl.f' 
      include 'zcouple.f'
      include 'sprods_com.f' 
      include 'ewcouple.f' 
      include 'scale.f' 
      include 'qcdcouple.f'
      include 'anom_higgs.f' 
      integer:: i2,i1,i7,hq,hg,h34,h56
      complex(dp):: amp_t(2,2,2,2),amp_b(2,2,2,2),
     & amp_th(2,2),amp_bh(2,2),
     & C0mt,C0mb,B12mt,B127mt,B12mb,B127mb,funcmt,funcmb,
     & qlI2,qlI3,fac_pro,fac_dec,mybmb,mybmt,
     & mtsfac,mbsfac,prop127,prop34,prop56,H4l(2,2),higgsprop
      real(dp):: s12,s127,mt2,mb2,sinthw,rescale

      s12=s(i1,i2)
      s127=s(i1,i2)+s(i1,i7)+s(i2,i7)
      mt2=mt**2
      mb2=mb**2
      sinthw=sqrt(xw)
      mtsfac=cplx1(mt2/(s127-s12))
      mbsfac=cplx1(mb2/(s127-s12))

!==== for width studies rescale by appropriate factor 
      if((keep_smhiggs_norm).and.(anom_higgs)) then 
         rescale=chi_higgs**2 
      else
         rescale=one
      endif


C----order for amp_t is hq,hg,h34,h56
      amp_t(:,:,:,:)=czip
      amp_b(:,:,:,:)=czip

!---- factor of (ta)^ij extracted
      fac_pro=im/16._dp/pi**2/rt2*4._dp*sqrt(gsq)**3*sqrt(esq)

!===== propagators 
      prop127=higgsprop(s127)
      prop34=cone/cplx2(s(3,4)-zmass**2,zmass*zwidth)
      prop56=cone/cplx2(s(5,6)-zmass**2,zmass*zwidth)
c--- Amplitudes for production 
      C0mt=qlI3(zip,s127,s12,mt2,mt2,mt2,musq,0)
      C0mb=qlI3(zip,s127,s12,mb2,mb2,mb2,musq,0)

      B12mt=qlI2(s12,mt2,mt2,musq,0)
      B12mb=qlI2(s12,mb2,mb2,musq,0)
      
      B127mt=qlI2(s127,mt2,mt2,musq,0)
      B127mb=qlI2(s127,mb2,mb2,musq,0)

      mybmt=B12mt-B127mt
      mybmb=B12mb-B127mb
       
c--- Amplitudes for decay
      H4l(1,1)=za(3,5)*zb(4,6)*l1*l2
     &        *wmass/(sinthw*(one-xw))*prop34*prop56
      H4l(2,1)=za(4,5)*zb(3,6)*r1*l2
     &        *wmass/(sinthw*(one-xw))*prop34*prop56
      H4l(1,2)=za(3,6)*zb(4,5)*l1*r2
     &        *wmass/(sinthw*(one-xw))*prop34*prop56
      H4l(2,2)=za(4,6)*zb(3,5)*r1*r2
     &        *wmass/(sinthw*(one-xw))*prop34*prop56

      fac_dec=-two*im*sqrt(esq)**3

c---- Indices of amp_th are (helicity of quark line, helicity of gluon)   

!===== Top loop amplitudes for production  
      funcmt=(C0mt*mt2/s12*(half-two*mtsfac)
     & +mtsfac/(s127-s12)*mybmt-mtsfac/s12)/(sinthw*wmass)

!===== -ve quark, +ve gluon 
      amp_th(1,2)=za(i2,i1)*zb(i2,i7)**2*funcmt
!===== +ve quark, +ve gluon by exchanging 1 and 2 
      amp_th(2,2)=za(i1,i2)*zb(i1,i7)**2*funcmt
!===== -ve quark, -ve gluon 
      amp_th(1,1)=za(i1,i7)**2*zb(i2,i1)*funcmt
!===== +ve quark, -ve gluon 
      amp_th(2,1)=za(i2,i7)**2*zb(i1,i2)*funcmt
 
 
!===== Bottom loop amplitudes for production  
      funcmb=(C0mb*mb2/s12*(half-two*mbsfac)
     & +mbsfac/(s127-s12)*mybmb-mbsfac/s12)/(sinthw*wmass)
!===== -ve quark, +ve gluon 
      amp_bh(1,2)=za(i2,i1)*zb(i2,i7)**2*funcmb
!===== +ve quark, +ve gluon by exchanging 1 and 2 
      amp_bh(2,2)=za(i1,i2)*zb(i1,i7)**2*funcmb
!===== -ve quark, -ve gluon 
      amp_bh(1,1)=za(i1,i7)**2*zb(i2,i1)*funcmb
!===== +ve quark, -ve gluon 
      amp_bh(2,1)=za(i2,i7)**2*zb(i1,i2)*funcmb
 

C---include im for Higgs propagator
      do hq=1,2
      do hg=1,2
      do h34=1,2
      do h56=1,2
      amp_t(hq,hg,h34,h56)=
     & amp_th(hq,hg)*H4l(h34,h56)*im*prop127*fac_pro*fac_dec*rescale
      amp_b(hq,hg,h34,h56)=
     & amp_bh(hq,hg)*H4l(h34,h56)*im*prop127*fac_pro*fac_dec*rescale
      enddo
      enddo
      enddo
      enddo

      return 
      end


