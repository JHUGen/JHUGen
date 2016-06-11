{
  gSystem->Load("libgfortran.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->AddIncludePath("-I$ROOFITSYS/include/");
  gSystem->Load("libmcfm_7p0.so");
  gSystem->Load("libME.so");
}
