//#########################################################################
//#                                                                       #
//#                               FROOT.C                                 #
//#                         Version June 2, 2008                          #
//#                                                                       #
//#########################################################################
//#                                                                       #
//# A simple interface to write numerical output from a standalone        #
//# Fortran or C program into a ROOT ntuple and perform various           #
//# manipulations on the output                                           #
//#                                                                       #
//#########################################################################
//# Send comments to Pavel Nadolsky (nadolsky@pa.msu.edu)


using namespace std;

// These are the C functions accessible from Fortran. They must be 
// described as 'extern "C"'

extern "C" {
  void initrootnt_(const char *title, const char *access, int ltitle, int laccess);
  void reinitrootnt_(const char *access, int laccess);
  void addntbranch_(float *element, const char *chtag, int ltag);
  void fillntbranch_(const char *chtag, int ltag);
  int getnumbranches_();
  void rootntoutp_();
  void printnt_();
  void teststr_(const char *str, int lstr);

}//extern "C"

//////////////////////////////////////////////////////////////////////////
// To link to the ROOT libraries, compile with a preprocessor option MYROOT
// (-DMYROOT) specified in Makefile. If MYROOT is not defined, the code 
// will be linked to dummy functions defined below

#ifdef MYROOT

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"


//Pointers to the ROOT application for writing the ntuple, 
//file with the ntuple, and the ntuple itself
TApplication* app;
TFile *outfile;   
TTree *nt;

//=========================================================== String routines
char* strtoupper(char *str)
  /* Converts a string str of the length strlen to uppercase letters */
{ int m=0;
  do
    str[m] = toupper(str[m]);
  while (str[m++]);
  return &(str[0]);
}



//--------------------------------------------------------------------------
char *truestr(const char *str, int len)
     // Converts a Fortran character string str of length len into a properly
     // formed C character string with '\0' at the end; trims blank spaces
     // at the end of the string
{
 int tlen=0;
 char tem;
/*Counts non-blank charachters in a string str until a first blank character
  or the end of the string is met*/

 while (str[tlen] !=' ' && (tlen < len ) && (tem=str[tlen++],tem)) ;


 char *tstr = new char[tlen+1];
 strncpy(tstr,str,tlen);
 tstr[tlen]='\0';

 return tstr;
}//truestr->


void teststr_(const char *str,int lstr)
{
  cout << "In C:"<<truestr(str,lstr) <<"==========="<<endl;
  return;
}



//------------------------------------------------------------------------

void initrootnt_(const char *title, const char *access, int ltitle, int laccess)
//Initialize the ROOT file for a given access mode 
//Inputs: title -- the name of the file (e.g., 'resbos.root')
//       access -- a string specifying the type of the access
//       (equivalent to the string 'option' in the TFile constructor)
//If access = 'NEW' or 'CREATE'   create a new file and open it 
//                           for writing. If the file already exists 
//                           the file is not opened.
              //           = 'RECREATE'      create a new file. If the file already
//                             exists, it will be overwritten.
//           = 'UPDATE'        open an existing file for writing.
//                             If no file exists, it is created.
//           = 'READ'          open an existing file for reading (default).
{
  char *appName="FROOT.C";
  int apprgc=1; 
  char *appargv[]={appName};
  app=new TApplication("App", &apprgc, appargv);

  outfile = new TFile(truestr(title,ltitle),strtoupper(truestr(access,laccess)));

  nt = (TTree*)outfile->Get("h10");
  if (nt == 0x0 ) 
    nt=new TTree("h10","h10");
    
  return;
} //initrootnt_ ->

//------------------------------------------------------------------------

void reinitrootnt_(const char *access, int laccess)
//Re-initialize the ROOT ntuple using a different access mode
//This instruction is used to change the access mode to the ntuple before
// it was closed by calling rootntoutp_(). See the example of the call of
// reinitrootnt_ in taste_froot.f. The access modes are described in
// the header of the subroutine initrootnt_(). 
{

  outfile->ReOpen(strtoupper(truestr(access,laccess)));
    
  nt = (TNtuple*)outfile->Get("h10");  

  if (nt == 0x0) 
    {
      cout << "STOP in reinitrootnt_: error in opening the ntuple" << endl;
      exit(10);
    }

  return;
} //reinitrootnt_ ->


//------------------------------------------------------------------------
  void addntbranch_(float *element, const char *chtag, int ltag)
    //Add an ntuple branch pointing at a Fortran variable element
{
  nt->Branch(truestr(chtag,ltag),element, truestr(chtag,ltag));

  return;
  }//addntbranch_ ->


//------------------------------------------------------------------------
void fillntbranch_(const char *chtag, int ltag)
 //Fill the variable into the branch with the name chtag, or all branches
 //if chtag='all'
{
  if (strcmp(truestr(chtag,ltag),"all") == 0)
      nt->Fill();
  else
     (nt->GetBranch(truestr(chtag,ltag)))->Fill();

  return;
}//fillntbranch_ ->

//------------------------------------------------------------------------
void rootntoutp_()
  //Close the ntuple file
{
  nt->Write("",TObject::kOverwrite);
  outfile->Close();
  return;
}//rootntoutp_ ->

//------------------------------------------------------------------------
int getnumbranches_()
//return the number of the branches in the ROOT ntuple
{
  return nt->GetNbranches();
}//getnumbranches_

//------------------------------------------------------------------------
void printnt_()
//Print the list of branches in the ntuple
{
  nt->Print();
}//print_



//------------------------------------------------------------------------

#else //MYROOT is not defined; link to dummy functions instead of the actual ROOT functions
//////////////////////////////////////////////////////////////////////////
 void initrootnt_(const char *title, const char *access, int ltitle, int laccess){return;}
 void reinitrootnt_(const char *access, int laccess){return;}
 void addntbranch_(float *element, const char *chtag, int ltag){return;}
 void fillntbranch_(const char *chtag, int ltag){return;}
 int getnumbranches_(){return;}
 void printnt_(){return;}
 void rootntoutp_(){return;}

#endif //MYROOT

