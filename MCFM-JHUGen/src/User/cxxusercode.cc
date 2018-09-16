#include <iostream>
#include <cmath>
#include <string>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
using namespace fastjet;
using namespace std;

const int mxpart = 12;
const int maxd   = 40;

extern "C" {

// integer nextnplot
// common/plotindex/nextnplot
extern struct {
  int nextnplot;
} plotindex_;


// integer maxd,ndmax
// parameter (maxd=40)
// double precision ptilde(0:maxd,mxpart,4)
// double precision ptildejet(0:maxd,mxpart,4)
// common/ptildes/ptilde,ptildejet,ndmax
extern struct {
  double ptilde[4][mxpart][maxd+1];
  double ptildejet[4][mxpart][maxd+1];
  int ndmax;
} ptildes_;


// integer npart
// common/npart/npart
extern struct {
  int npart;
} npart_;


// character*2 plabel(mxpart)
// common/plabel/plabel
extern struct {
  char plabel[mxpart][2];
} plabel_;


void bookplot_(const int & iplot,  //< the index of the plot 
               const char * tag,   //< either "book" or "plot"
               const char * title, //< the plot title
               const double & var,
               const double & wt,
               const double & wt2,
               const double & xmax,
               const double & xmin,
               const double & dx,
               const char * llplot, // "lin" or "log"
               const int tag_len,   // pass string lengths by value!
               const int title_len  // pass string lengths by value!
               );

} // extern "C"


//----------------------------------------------------------------------
/// Interface to carry out histogramming with the MCFM mbook package.
///
/// When calling from C++, use this interface and not the f77
/// bookplot_ one defined above, because this one correctly handles the
/// the string aspects
void bookplot( const int & iplot,         //< the index of the plot 
               const std::string & tag,   //< either "book" or "plot"
               const std::string & title, //< the plot title
               const double & var,        //< the value of the variable
               const double & wt,         //< event weight
               const double & wt2,        //< squared event weight
               const double & xmax,       //< lower limit of histogram
               const double & xmin,       //< upper limit of histogram
               const double & dx,         //< bin spacing
               const std::string & llplot // "lin" or "log"
               ) {
  bookplot_(iplot, tag.c_str(), title.c_str(), 
            var, wt, wt2, xmax, xmin, dx, llplot.c_str(), 
            tag.size(),
            title.size());
}


//----------------------------------------------------------------------
/// return the scalar sum of the pt's of all particles
double ht_particles(int nd) {
  double ht = 0.0;
  for (int i = 2; i < 2+npart_.npart; i++) {
    ht += sqrt(pow(ptildes_.ptilde[0][i][nd],2) +
               pow(ptildes_.ptilde[1][i][nd],2));
  }
  return ht;
}

//----------------------------------------------------------------------
/// returns a vector of the particles in the event; if the offset
/// is zero, then it corresponds to all outgoing particles; 
///
/// if the offset is higher, then the first "offset" particles are
/// discarded;
///
/// if the offset is -2, then the incoming particles are the first two
///
vector<PseudoJet> mcfm_particles(int nd, int offset = 0) {
  vector<PseudoJet> result(npart_.npart-offset);
  for (int i = 2 + offset; i < 2+npart_.npart; i++) {
    result[i-(2+offset)].reset(ptildes_.ptilde[0][i][nd],
                               ptildes_.ptilde[1][i][nd],
                               ptildes_.ptilde[2][i][nd],
                               ptildes_.ptilde[3][i][nd]);
  }
  return result;
}



//----------------------------------------------------------------------
// the two user routines called by MCFM
extern "C" {

bool userincludedipole_(const int & nd, 
                        const double * ppart , 
                        const bool & mcfm_result) {
  bool result = mcfm_result;

  if (ht_particles(nd) < 1000.0) return false;
  
  return result;
}


//----------------------------------------------------------------------
void userplotter_( const double * pjet, 
                   const double & wt,
                   const double & wt2,
                   const int    & nd) {

  static bool first = true;
  string tag;
  if (first) {
    first = false;
    tag = "book";
  } else {
    tag = "plot";
  }

  double ht = ht_particles(nd);

  vector<PseudoJet> all_particles = mcfm_particles(nd,0);
  vector<PseudoJet> hadrons       = mcfm_particles(nd,2);
  
  JetDefinition jet_def(antikt_algorithm, 0.7);
  ClusterSequence cs_all(all_particles, jet_def);
  ClusterSequence cs_hadrons(hadrons, jet_def);
  vector<PseudoJet> jets_all       = sorted_by_pt(cs_all.    inclusive_jets());
  vector<PseudoJet> jets_hadrons   = sorted_by_pt(cs_hadrons.inclusive_jets());

  int iplot = plotindex_.nextnplot;
  
  double ptmax = 5000;
  double ptbin = 100;

  bookplot(iplot, tag, "UserHT", ht, wt, wt2, 0.0, 2*ptmax, ptbin, "lin");
  iplot++;

  // get the inclusive jet distribution, made up of all particles
  if (tag == "plot") {
    for (unsigned i = 0; i < jets_all.size(); i++) {
      bookplot(iplot, tag, "incl_all", jets_all[i].perp(), wt, wt2, 
               0.0, ptmax, ptbin, "lin");
    }
  } else {
    bookplot(iplot, tag, "incl_all", 0, wt, wt2, 0.0, ptmax, ptbin, "lin");
  }
  iplot++;

  // get the inclusive jet distribution, made up of hadrons
  if (tag == "plot") {
    for (unsigned i = 0; i < jets_hadrons.size(); i++) {
      bookplot(iplot, tag, "incl_hadrons", jets_hadrons[i].perp(), wt, wt2, 
               0.0, ptmax, ptbin, "lin");
    }
  } else {
    bookplot(iplot, tag, "incl_hadrons", 0, wt, wt2, 0.0, ptmax, ptbin, "lin");
  }
  iplot++;

}


} // extern "C"
