      subroutine writetable(
     & A6treemp,A6treepm,A6treemm,A6treepp,
     & ALCmp,ALCpm,ALCmm,ALCpp,
     & ASLmp,ASLpm,ASLmm,ASLpp,
     & Afmp,Afpm,Afmm,Afpp)
      implicit none
      include 'constants.f'
      integer j
      double complex 
     & A6treemp,A6treepm,A6treemm,A6treepp,
     & ALCmp(-2:0),ALCpm(-2:0),ALCmm(-2:0),ALCpp(-2:0),
     & ASLmp(-2:0),ASLpm(-2:0),ASLmm(-2:0),ASLpp(-2:0),
     & Afmp(-2:0),Afpm(-2:0),Afmm(-2:0),Afpp(-2:0)
      double precision p(mxpart,4),wmass,mb
      
      do j=-2,0
      ALCmp(j)=ALCmp(j)/A6treemp
      ALCpm(j)=ALCpm(j)/A6treepm
      ALCmm(j)=ALCmm(j)/A6treemm
      ALCpp(j)=ALCpp(j)/A6treepp
      if (abs(dimag(ALCmp(j))).lt.1d-15) 
     & ALCmp(j)=Dcmplx(dreal(ALCmp(j)),0d0)
      if (abs(dimag(ALCpm(j))).lt.1d-15) 
     & ALCpm(j)=Dcmplx(dreal(ALCpm(j)),0d0)
      if (abs(dimag(ALCmm(j))).lt.1d-15) 
     & ALCmm(j)=Dcmplx(dreal(ALCmm(j)),0d0)
      if (abs(dimag(ALCpp(j))).lt.1d-15) 
     & ALCpp(j)=Dcmplx(dreal(ALCpp(j)),0d0)

      ASLmp(j)=ASLmp(j)/A6treemp
      ASLpm(j)=ASLpm(j)/A6treepm
      ASLmm(j)=ASLmm(j)/A6treemm
      ASLpp(j)=ASLpp(j)/A6treepp
      if (abs(dimag(ASLmp(j))).lt.1d-15) 
     & ASLmp(j)=Dcmplx(dreal(ASLmp(j)),0d0)
      if (abs(dimag(ASLpm(j))).lt.1d-15) 
     & ASLpm(j)=Dcmplx(dreal(ASLpm(j)),0d0)
      if (abs(dimag(ASLmm(j))).lt.1d-15) 
     & ASLmm(j)=Dcmplx(dreal(ASLmm(j)),0d0)
      if (abs(dimag(ASLpp(j))).lt.1d-15) 
     & ASLpp(j)=Dcmplx(dreal(ASLpp(j)),0d0)

      Afmp(j)=Afmp(j)/A6treemp
      Afpm(j)=Afpm(j)/A6treepm
      Afmm(j)=Afmm(j)/A6treemm
      Afpp(j)=Afpp(j)/A6treepp
      if (abs(dimag(Afmp(j))).lt.1d-15) 
     & Afmp(j)=Dcmplx(dreal(Afmp(j)),0d0)
      if (abs(dimag(Afpm(j))).lt.1d-15) 
     & Afpm(j)=Dcmplx(dreal(Afpm(j)),0d0)
      if (abs(dimag(Afmm(j))).lt.1d-15) 
     & Afmm(j)=Dcmplx(dreal(Afmm(j)),0d0)
      if (abs(dimag(Afpp(j))).lt.1d-15) 
     & Afpp(j)=Dcmplx(dreal(Afpp(j)),0d0)
      enddo
      open(unit=66,file='table1.tex',status='unknown')
      open(unit=67,file='table2.tex',status='unknown')

      write(66,*) '\\begin{table}'
      write(66,*) '\\begin{tabular}{|c||c|c|c||c|c|c|}'
      write(66,*) '\\hline'
      write(66,*) '&\\multicolumn{3}{c||}{',
     & '$A_6(1_q^-,2_\\Qb^+,3_Q^-,4_',
     & '\\qb^+,5_{\\bar{\\ell}}^+,6_{\\ell}^-)$}'
      write(66,*) '&\\multicolumn{3}{c||}{',
     & '$A_6(1_q^-,2_\Qb^-,3_Q^+,4_\qb^+,',
     & '5_{\\bar{\\ell}}^+,6_{\\ell}^-)$}' 
      write(66,*) '\\\\'
      write(66,*) '\\hline'
      write(66,*) '\\hline'
      write(66,*) '& $1/\\e^2$ & $1/\\e$ & $\\e^0$'
      write(66,*) '&$1/\\e^2$ & $1/\e$ & $\\e^0$ \\\\'
      write(66,*) '\\hline'
      write(66,*) '\\hline $A_6^{\\rm tree}$ '
      write(66,88) 0d0,0d0,dreal(A6treemp),0d0,0d0,dreal(A6treepm),
     &             0d0,0d0,dimag(A6treemp),0d0,0d0,dimag(A6treepm)
      write(66,*) '\\hline $A/A_6^{\\rm tree}$'
      write(66,88) 
     & dreal(ALCmp(-2)),
     & dreal(ALCmp(-1)),
     & dreal(ALCmp( 0)),
     & dreal(ALCpm(-2)),
     & dreal(ALCpm(-1)),
     & dreal(ALCpm( 0)),
     & dimag(ALCmp(-2)),
     & dimag(ALCmp(-1)),
     & dimag(ALCmp( 0)),
     & dimag(ALCpm(-2)),
     & dimag(ALCpm(-1)),
     & dimag(ALCpm( 0))
      write(66,*) '\\hline $A^{sl}/A_6^{\\rm tree}$'
      write(66,88) 
     & dreal(ASLmp(-2)),
     & dreal(ASLmp(-1)),
     & dreal(ASLmp( 0)),
     & dreal(ASLpm(-2)),
     & dreal(ASLpm(-1)),
     & dreal(ASLpm( 0)),
     & dimag(ASLmp(-2)),
     & dimag(ASLmp(-1)),
     & dimag(ASLmp( 0)),
     & dimag(ASLpm(-2)),
     & dimag(ASLpm(-1)),
     & dimag(ASLpm( 0))

      write(66,*) '\\hline $A^f/A_6^{\\rm tree}$ '
      write(66,88) 
     & dreal(Afmp(-2)),
     & dreal(Afmp(-1)),
     & dreal(Afmp( 0)),
     & dreal(Afpm(-2)),
     & dreal(Afpm(-1)),
     & dreal(Afpm( 0)),
     & dimag(Afmp(-2)),
     & dimag(Afmp(-1)),
     & dimag(Afmp( 0)),
     & dimag(Afpm(-2)),
     & dimag(Afpm(-1)),
     & dimag(Afpm( 0))

      write(66,*) '\\hline'

   88 format(6('& $',f19.11,'   $ '),'\\\\ \n',
     .                   6('& $',f19.11,'\\,i$ '),'\\\\ \n')   
      write(66,*) '\\hline'
      write(66,*) '\\end{tabular}'
      write(66,*) '\\caption{Numerical values',
     & ' of primitive amplitudes'
      write(66,*) 'at the kinematic point defined in'
      write(66,*) 'Eq.~(\\ref{kinpoint}).'
      write(66,*) '\\label{numres}}'
      write(66,*) '\\end{table}'


      write(67,*) '\\begin{table}'
      write(67,*) '\\begin{tabular}{|c||c|c|c||c|c|c|}'
      write(67,*) '\\hline'
      write(67,*) '&\\multicolumn{3}{c||}{',
     & '$A_6(1_q^-,2_\\Qb^-,3_Q^-,4_',
     & '\\qb^+,5_{\\bar{\\ell}}^+,6_{\\ell}^-)$}'
      write(67,*) '&\\multicolumn{3}{c||}{',
     & '$A_6(1_q^-,2_\Qb^+,3_Q^+,4_\qb^+,',
     & '5_{\\bar{\\ell}}^+,6_{\\ell}^-)$}' 
      write(67,*) '\\\\'
      write(67,*) '\\hline'
      write(67,*) '\\hline'
      write(67,*) '& $1/\\e^2$ & $1/\\e$ & $\\e^0$'
      write(67,*) '&$1/\\e^2$ & $1/\e$ & $\\e^0$ \\\\'
      write(67,*) '\\hline'
      write(67,*) '\\hline $A_6^{\\rm tree}$ '
      write(67,88) 0d0,0d0,dreal(A6treemm),0d0,0d0,dreal(A6treepp),
     &             0d0,0d0,dimag(A6treemm),0d0,0d0,dimag(A6treepp)
      write(67,*) '\\hline $A/A_6^{\\rm tree}$'
      write(67,88) 
     & dreal(ALCmm(-2)),
     & dreal(ALCmm(-1)),
     & dreal(ALCmm( 0)),
     & dreal(ALCpp(-2)),
     & dreal(ALCpp(-1)),
     & dreal(ALCpp( 0)),
     & dimag(ALCmm(-2)),
     & dimag(ALCmm(-1)),
     & dimag(ALCmm( 0)),
     & dimag(ALCpp(-2)),
     & dimag(ALCpp(-1)),
     & dimag(ALCpp( 0))
      write(67,*) '\\hline $A^{sl}/A_6^{\\rm tree}$'
      write(67,88) 
     & dreal(ASLmm(-2)),
     & dreal(ASLmm(-1)),
     & dreal(ASLmm( 0)),
     & dreal(ASLpp(-2)),
     & dreal(ASLpp(-1)),
     & dreal(ASLpp( 0)),
     & dimag(ASLmm(-2)),
     & dimag(ASLmm(-1)),
     & dimag(ASLmm( 0)),
     & dimag(ASLpp(-2)),
     & dimag(ASLpp(-1)),
     & dimag(ASLpp( 0))

      write(67,*) '\\hline $A^f/A_6^{\\rm tree}$ '
      write(67,88) 
     & dreal(Afmm(-2)),
     & dreal(Afmm(-1)),
     & dreal(Afmm( 0)),
     & dreal(Afpp(-2)),
     & dreal(Afpp(-1)),
     & dreal(Afpp( 0)),
     & dimag(Afmm(-2)),
     & dimag(Afmm(-1)),
     & dimag(Afmm( 0)),
     & dimag(Afpp(-2)),
     & dimag(Afpp(-1)),
     & dimag(Afpp( 0))

      write(67,*) '\\hline'

      write(67,*) '\\hline'
      write(67,*) '\\end{tabular}'
      write(67,*) '\\caption{Numerical values',
     & ' of primitive amplitudes'
      write(67,*) 'at the kinematic point defined in'
      write(67,*) 'Eq.~(\\ref{kinpoint}).'
      write(67,*) '\\label{numres}}'
      write(67,*) '\\end{table}'

      include 'MCFMpoint.f'
      
      write(6,*) 'mb   =',mb
      write(6,*) 'wmass=',wmass
      do j=1,6
      write(6,77) j,p(j,4),p(j,1),p(j,2),p(j,3)
      enddo
      
      stop
      return
      
   77 format('p_',i1,' &=& (',3(f16.11,','),f19.11,'), \\nn \\\\')   
      
      end
