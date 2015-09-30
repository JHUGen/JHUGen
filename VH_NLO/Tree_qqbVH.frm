#define WorkPath "/afs/cern.ch/user/m/maschulz/1git_JHUGen/VH_NLO/"
#define setP1P1 "0"
#define setP2P2 "0"
#define setP3P3 "MZ^2"
#define setP4P4 "MH^2"
#define setDST   "4"
#define setDSTm4 "0"
#define InterfereDiags1 "1,1"
#define InterfereDiags2 "1,1"


#include /afs/cern.ch/user/m/maschulz/lib/FeynArtsToForm/header2.frm

#include `WorkPath'Tree_qqbVH_input.frm

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
Bracket IZ,SME;
.sort;





.end
