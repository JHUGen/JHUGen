      subroutine computepdfuncertainty(PDFarray,PDFcentral,
     & PDFperror,PDFnerror,PDFerror)
      implicit none
      include 'types.f'
c--- Routine to compute PDF uncertainty given an array of results
c--- for different PDF uncertainty sets (PDFarray).
c--- Returns central value (PDFcentral), positive and negative
c--- excursions (PDFperror, PDFnerror) and symmetric uncertainty (PDFerror).
c--- Note that the appropriate method for computing uncertainties
c--- is chosen according to the value of "PDFname" used in LHAPDF
c--- (as a result the logic may need to be updated in the future,
c---  it is current as of 9/2013)
c---
c--- This routine is similar to the one provided natively in LHAPDFv5.8.5
c--- onwards ("uncertainties.f", written by Graeme Watt)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'lhapdf.f'
      include 'PDFerrors.f'
      integer:: j
      real(dp):: PDFarray(0:1000),PDFcentral,PDFperror,PDFnerror,
     & PDFerror,sum1,sum2,PDFMCav,PDFMCer
      logical:: first
      data first/.true./
      save first

c--- NNPDF: just compute average and standard deviation    
      if ((index(PDFname,'NNPDF') > 0) .or.
     &    (index(PDFname,'nnpdf') > 0)) then
        if (first) then
        write(6,*)'****************************************************'
        write(6,*)'*    Using MC prescription for PDF uncertainties,  *'
        write(6,*)'*    appropriate for NNPDF sets                    *'
        write(6,*)'*                                                  *'
        write(6,*)'*   (for details and references, see Eqn. (158)    *'
        write(6,*)'*    in Appendix B of arXiv:0808.1231 [hep-ph])    *'
        write(6,*)'****************************************************'
        first=.false.
        endif
        sum1=zip
        sum2=zip
        do j=1,maxPDFsets
           if (PDFarray(j) .ne. 0.) then
             sum1=sum1+PDFarray(j)
             sum2=sum2+PDFarray(j)**2._dp
           endif
        enddo
        PDFMCav = sum1/maxPDFsets
        PDFMCer = sqrt(sum2/maxPDFsets -  PDFMCav**2._dp)

        PDFperror=PDFMCer
        PDFnerror=PDFMCer
        PDFerror=PDFMCer
        PDFcentral=PDFMCav
        return
      endif      

c--- Alekhin et al: Hessian approach (symmetric)     
      if ((index(PDFname,'ABM') > 0) .or.
     &    (index(PDFname,'ABKM') > 0) .or.
     &    (index(PDFname,'A02M') > 0).or.
     &    (index(PDFname,'abm') > 0) .or.
     &    (index(PDFname,'abkm') > 0) .or.
     &    (index(PDFname,'a02m') > 0)) then
        if (first) then
        write(6,*)'****************************************************'
        write(6,*)'*   Using symmetric Hessian prescription for PDF   *'
        write(6,*)'*   uncertainties, appropriate for Alekhin et al.  *'
        write(6,*)'****************************************************'
        first=.false.
        endif
        PDFerror=zip
        do j=1,maxPDFsets
           if (PDFarray(j) .ne. 0.) then
             PDFerror=PDFerror+(PDFarray(j)-PDFarray(0))**2._dp
           endif
        enddo
        PDFerror=sqrt(PDFerror)
        PDFperror=PDFerror
        PDFnerror=PDFerror
        PDFcentral=PDFarray(0)
        return
      endif      

c--- everyone else (CTEQ and MSTW): Hessian approach (asymmetric)     
        if (first) then
        write(6,*)'****************************************************'
        write(6,*)'*   Using asymmetric Hessian prescription for PDF  *'
        write(6,*)'*   uncertainties, appropriate for CTEQ, MSTW      *'
        write(6,*)'*                                                  *'
        write(6,*)'*        (see, for example Eqn. (43) of            *'
        write(6,*)'*         J.Campbell, J.Huston, W.J.Stirling,      *'
        write(6,*)'*         Rep. Prog. Phys. 70 (2007) 89)           *'
        write(6,*)'****************************************************'
        first=.false.
        endif
        PDFperror=zip
        PDFnerror=zip
        PDFerror=zip
        do j=1,maxPDFsets-1,2
           if (PDFarray(j) .ne. 0.) then
              PDFperror=PDFperror+max(zip,
     &        PDFarray(j)-PDFarray(0),PDFarray(j+1)-PDFarray(0))**2
              PDFnerror=PDFnerror+max(zip,
     &        PDFarray(0)-PDFarray(j),PDFarray(0)-PDFarray(j+1))**2
              PDFerror=PDFerror+(PDFarray(j)-PDFarray(j+1))**2
           endif
        enddo
        PDFerror=half*sqrt(PDFerror)
        PDFperror=sqrt(PDFperror)
        PDFnerror=sqrt(PDFnerror)
        PDFcentral=PDFarray(0)
        return

      return
      end
      
      
      
      
