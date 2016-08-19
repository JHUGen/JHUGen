*      $Id: qpff.h,v 1.1 1995/12/12 10:03:48 gj Exp $
*      -------------------------------------------------------------
*      INCLUDE FILE FOR THE QPFF ROUTINES.  
*                              Geert Jan van Oldenborgh.
*      -------------------------------------------------------------
*       please do not change, and recompile _everything_ when you do.
*      -------------------------------------------------------------
*
*      this parameter determines how far the scalar npoint functions
*      will look back to find the same parameters (when lmem is true)
*
      integer memory
      parameter(memory=12)
*
*              if .TRUE. then                        default (ffinit)
*      lwrite: give debug output                     .FALSE.
*      ltest:  test consistency internally (slow)    .TRUE.
*      l4also: in C0 (and higher), also consider the algorithm with 16
*              dilogs                                .TRUE.
*      ldc3c4: in D0 (and higher), also consider possible cancellations
*              between the C0s                      .TRUE.
*      lmem:   before computing the C0 and higher, first check whether
*              it has already been done recently     .FALSE.
*      lwarn:  give warning messages (concerning numerical stability)
*                                                    .TRUE.
*      ldot:   leave the dotproducts and some determinants in common
*                                                    .FALSE.
*      onshel: (in ffz?0 only): use onshell momenta  .TRUE.
*      lsmug:  internal use
*      lnasty: internal use
*
      logical lwrite,ltest,l4also,ldc3c4,lmem,lwarn,ldot,onshel,lsmug,
     +            lnasty
*
*      nwidth: number of widths within which the complex mass is used
*      nschem: scheme to handle the complex mass (see ffinit.f)
*      idot:      internal flags to signal that some of the dotproducts 
*            are input: 0: none; 1: external pi.pj, 2: external + 
*            kinematical determinant, 3: all dotproducts + kindet.
*
      integer nwidth,nschem,idot
*
*      xloss:  factor that the final result of a subtraction can be
*              smaller than the terms without warning (default 1/8)
*      precx:  precision of real numbers, determined at runtime by
*              ffinit (IEEE: 4.e-16)
*      precc:  same for complex numbers
*      xalogm: smallest real number of which a log can be taken,
*              determined at runtime by ffinit (IEEE: 2.e-308)
*      xclogm: same for complex.
*      xalog2: xalogm**2
*      xclog2: xclogm**2
*      reqprc: not used
*      x[0124]:0,1,2,4
*      x05:      1/2
*      pi:      pi
*      pi6:      pi**2/6
*      pi12:      pi**2/12
*      xlg2:      log(2)
*      bf:      factors in the expansion of dilog (~Bernouilli numbers)
*      xninv:      1/n
*      xn2inv:      1/n**2
*      xinfac:      1/n!
*      fpij2:      vi.vj for 2point function 1-2: si, 3-3:  pi
*      fpij3:      vi.vj for 3point function 1-3: si, 4-6:  pi
*      fpij4:      vi.vj for 4point function 1-4: si, 5-10: pi
*      fpij5:      vi.vj for 5point function 1-5: si, 6-15: pi
*      fpij6:      vi.vj for 6point function 1-6: si, 7-21: pi
*      fdel2:      del2 = delta_(p1,p2)^(p1,p2) = p1^2.p2^2 - p1.p2^2 in C0
*      fdel3:      del3 = delta_(p1,p2,p3)^(p1,p2,p3) in D0
*      fdel4s:      del4s = delta_(s1,s2,s3,s4)^(s1,s2,s3,s4) in D0
*      fdel4:      del4 = delta_(p1,p2,p3,p4)^(p1,p2,p3,p4) in E0
*      fdl3i:      del3i = delta_(pj,pk,pl)^(pj,pk,pl) in E0, D0 without si
*      fdl4si:      dl4si = del4s in E0, D0 without si
*      fdl3ij:      same in F0 without si and sj.
*      fd4sij:      dl4si = del4s in E0, D0 without si
*      fdl4i:      delta4 in F0 without si.
*      fodel2:      same offshell (in case of complex or z-functions)
*      fodel3:      -"-
*      cfdl4s:      -"-
*      fodel4:      -"-
*      fodl3i:      -"-
*      fod3ij:      -"-
*      fodl4i:      -"-
*      fidel3:      ier of del3 (is not included in D0)
*      fidel4:      ier of del4 (is not included in E0)
*      fidl3i:      ier of dl3i (is not included in E0)
*      fid3ij:      ier of dl3ij (is not included in F0)
*      fidl4i:      ier of dl4i (is not included in F0)
*
      real(kind=qp) xloss,precx,precc,xalogm,xclogm,xalog2,xclog2,
     +            reqprc,x0,x05,x1,x2,x4,pi,pi6,pi12,xlg2,bf(20),
     +            xninv(30),xn2inv(30),xinfac(30),
     +            fpij2(3,3),fpij3(6,6),fpij4(10,10),fpij5(15,15),
     +            fpij6(21,21),fdel2,fdel3,fdel4s,fdel4,fdl3i(5),
     +            fdl4si(5),fdl3ij(6,6),fd4sij(6,6),fdl4i(6),fodel2,
     +            fodel3,fodel4,fodl3i(5),fod3ij(6,6),fodl4i(6)
      real(kind=qp) zip,zero,one,half,two,ten
      integer fidel3,fidel4,fidl3i(5),fid3ij(6,6),fidl4i(6)
*
*      c[0124]:0,1,2,4 complex
*      c05:      1/2 complex
*      c2ipi:      2*i*pi
*      cipi2:      i*pi**2
*      cfp..:      complex version of fp..., only defined in ff[cz]*
*      cmipj:      (internal only) mi^2 - pj^2 in C0
*      c2sisj:      (internal only) 2*si.sj in D0
*      cfdl4s:      del4s in complex case (D0)
*      ca1:      (internal only) complex A1
*      csdl2p: (internal only) complex transformed sqrt(del2)
*
      complex(kind=qp) c0,c05,c1,c2,c4,c2ipi,cipi2,
     +            cfpij2(3,3),cfpij3(6,6),cfpij4(10,10),cfpij5(15,15),
     +            cfpij6(21,21),cmipj(3,3),c2sisj(4,4),cfdl4s,ca1
*
*      nevent:      number in integration loop (to be updated by user)
*      ner:      can be used to signal numerical problems (see ffrcvr)
*      id:      identifier of scalar function (to be set by user)
*      idsub:      internal identifier to pinpoint errors
*      inx:      in D0: p(inx(i,j)) = isgn(i,j)*(s(i)-s(j))
*      inx5:      in E0: p(inx5(i,j)) = isgn5(i,j)*(s(i)-s(j))
*      inx6:      in F0: p(inx6(i,j)) = isgn6(i,j)*(s(i)-s(j))
*      isgn:      see inx
*      isgn5:      see inx5
*      isgn6:      see inx6
*      iold:      rotation matrix for 4point function
*      isgrot:      signs to iold
*      isgn34:      +1 or -1: which root to choose in the transformation (D0)
*      isgnal:      +1 or -1: which root to choose in the alpha-trick (C0)
*      irota3:      save the number of positions the C0 configuration has been 
*            rotated over
*      irota4:      same for the D0
*      irota5:      same for the E0
*      irota6:      same for the F0
*
      integer nevent,ner,id,idsub,inx(4,4),isgn(4,4),inx5(5,5),
     +            isgn5(5,5),inx6(6,6),isgn6(6,6),isgn34,isgnal,iold(13,
     +            12),isgrot(10,12),irota3,irota4,irota5,irota6
      integer idum93(2)
*
*      parameters

      parameter(zip=0._qp,zero=0._qp,one=1._qp,half=0.5_qp,two=2._qp,
     & ten=10._qp)
      parameter(x0=0._qp,x1=1._qp,x05=.5_qp,x2=2._qp,x4=4._qp,
     +           c0=(0._qp,0._qp),c05=(.5_qp,0._qp),c1=(1._qp,0._qp),
     +           c2=(2._qp,0._qp),c4=(4._qp,0._qp))
      parameter(
     + c2ipi=cmplx(0._qp,6.28318530717958647692528676655896_qp,qp),
     + cipi2=cmplx(0._qp,9.869604401089358618834490999876_qp,qp),
     + pi =3.14159265358979323846264338327948_qp,
     + pi6=1.644934066848226436472415166646_qp,
     + pi12=.822467033424113218236207583323_qp,
     + xlg2=.6931471805599453094172321214581_qp)
*
*      common
*
      common /qpffsign/isgn34,isgnal
      common /qpffprec/ xloss,precx,precc,xalogm,xclogm,xalog2,xclog2,
     +            reqprc
      common /qpffflag/ lwrite,ltest,l4also,ldc3c4,lmem,lwarn,ldot,
     +            nevent,ner,id,idsub,nwidth,nschem,onshel,idot
      common /qpffcnst/ bf,xninv,xn2inv,xinfac,inx,isgn,iold,isgrot,
     +            inx5,isgn5,inx6,isgn6
      common /qpffrota/ irota3,irota4,irota5,irota6
      common /qpffdot/ fpij2,fpij3,fpij4,fpij5,fpij6
      common /qpffdel/ fdel2,fdel3,fdel4s,fdel4,fdl3i,fdl4si,fdl3ij,
     +            fd4sij,fdl4i
      common /qpffidel/ fidel3,fidel4,fidl3i,fid3ij,fidl4i
      common /qpffcdot/ cfpij2,cfpij3,cfpij4,cfpij5,cfpij6
      common /qpffcdel/ fodel2,fodel3,cfdl4s,fodel4,fodl3i,fod3ij,fodl4i
      common /qpffsmug/ lsmug,lnasty,idum93,cmipj,c2sisj,ca1
!$omp threadprivate(/qpffsign/,/qpffprec/,/qpffflag/,/qpffcnst/,/qpffrota/,
!$omp& /qpffdot/,/qpffdel/,/qpffidel/,/qpffcdot/,/qpffsmug/)
