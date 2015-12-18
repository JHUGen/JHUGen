#define WorkPath "/home/schulze/1git_JHUGen/VH_NLO/"
#define setP1P1 "0"
#define setP2P2 "0"
#define setP3P3 "MZ^2"
#define setP4P4 "MH^2"
#define setDST   "4"
#define setDSTm4 "0"
#define InterfereDiags1 "1,1"
#define InterfereDiags2 "1,1"
#include /home/schulze/lib/FeynArtsToForm/header2.frm


** original FeynArt output
*#include `WorkPath'Tree_qqbVH_input.frm

** my modified anomalous version

#include `WorkPath'Tree_qqbVH_BSM_input.frm



* up quark masses
id MU=0;
argument;
   id MU=0;
endargument;

id IndexDelta(Col1?,Col2?) = 1;
id DID(Col1?)=1;


#include `WorkPath'Tree_qqbVH_calc.frm





Format Mathematica;
Print;
Bracket a1HVV,a2HVV,a3HVV,vev;
.sort;


#write <`WorkPath'Tree_qqbVH_output.dat> "NumAmps = `NumAmps';\n";
#write <`WorkPath'Tree_qqbVH_output.dat> "AmpList = \"  `AmpList'\";\n";

#do TheAmp = `AmpList'
  #write <`WorkPath'Tree_qqbVH_output.dat> "M`TheAmp' = (%E);\n ", `TheAmp';
#enddo




.end
