InitPath =  $HomeDirectory<>"/lib/MathLib/";
FeynArtsPath =  $HomeDirectory<>"/lib/FeynArts-3.9/";
FeynArtsToFormPath =  $HomeDirectory<>"/lib/FeynArtsToForm/";
ProjectPath =  $HomeDirectory<>"/1git_JHUGen/VH_NLO/";
Unprotect[Diagonal];
Clear[Diagonal];

Get[ FeynArtsPath<>"FeynArts.m" ];
Get[ FeynArtsToFormPath<>"ExtractAmp.m" ];


TheProcess = {V[5],V[5]} -> {V[2],S[1]};

TheLevel = {Classes};
TheGenericModel = "Lorentz";
TheModel = SMQCD;
PaintIt  = True;
WriteAmp = True;

SetOptions[CreateTopologies, ExcludeTopologies -> {Tadpoles,WFCorrections}];


SetOptions[InsertFields,
    	InsertionLevel -> TheLevel,
    	GenericModel -> TheGenericModel ,
    	Model -> TheModel,
        ExcludeParticles -> { SV[2],SV[3],S[2],U[1],U[2],U[3] },
        LastSelections-> {},
        Restrictions -> { NoLightFHCoupling }
];

SetOptions[Paint, Numbering->Simple, ColumnsXRows->{5,7} ];
CKM = IndexDelta;



(* generate Tree level diagrams *)
TheTreeTop  = CreateTopologies[ 1, 2 -> 2 ];



TheTreeDiag = InsertFields[ TheTreeTop, TheProcess ];

If[PaintIt,
   Paint[TheTreeDiag, DisplayFunction->(Display[ProjectPath<>"Tree_ggVH.ps",#]&) ];
];

TheTreeAmp   = CreateFeynAmp[ TheTreeDiag,   PreFactor-> 1 ];

Do[  TheTreeAmp[[i,3]]  *=DID[i] ,{i,1,Length[ TheTreeAmp]} ];

(* removing diagrams with closed dn-quark loops *)
(* TheTreeAmp2 = Delete[TheTreeAmp,{{5},{6},{14}}]; *)
TheTreeAmp2 = TheTreeAmp;


MyTreeAmp   = ExtractAmp[ TheTreeAmp2 ];
If[WriteAmp,
   AmpToFORM[ Join[MyTreeAmp] ,ProjectPath<>"Tree_ggVH_input.frm"];
];

