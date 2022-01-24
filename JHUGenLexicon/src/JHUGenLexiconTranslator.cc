#include <cassert>
#include <exception>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "JHUGenLexiconCouplings.h"
#include "JHUGenLexiconTranslator.h"
#include "JHUGenLexiconHelperFunctions.h"


using namespace std;
using namespace JHUGenLexiconHelperFunctions;
using namespace JHUGenLexiconIOHelpers;
using namespace JHUGenLexiconCouplings;


JHUGenLexiconTranslator::JHUGenLexiconTranslator(JHUGenLexiconOptionParser const& opts_) :
  opts(opts_)
{
  translate();
}
void JHUGenLexiconTranslator::translate(){
  auto const& input_flags = opts.getInputFlags();
  auto const& input_parameters = opts.getInputParameters();
  auto const& input_couplings = opts.getInputCouplings();
  auto const& basis_input = opts.getInputBasis();
  auto const& basis_output = opts.getOutputBasis();
  bool include_triple_quartic_gauge; getValueWithDefault<std::string, bool>(input_flags, "include_triple_quartic_gauge", include_triple_quartic_gauge, false); 
  bool HW_couplings_only; getValueWithDefault<std::string, bool>(input_flags, "HW_couplings_only", HW_couplings_only, false);
  bool HZ_couplings_only; getValueWithDefault<std::string, bool>(input_flags, "HZ_couplings_only", HZ_couplings_only, false);
  bool custodial_symmetry; getValueWithDefault<std::string, bool>(input_flags, "custodial_symmetry", custodial_symmetry, false);
  double vev_lam; getValueWithDefault<std::string, double>(input_parameters, "vev_lam", vev_lam, DEFVAL_VEV_LAM);
  double delta_v; getValueWithDefault<std::string, double>(input_parameters, "delta_v" , delta_v ,DEFVAL_DELTA_V);
  double delta_m; getValueWithDefault<std::string, double>(input_parameters, "delta_m" , delta_m ,DEFVAL_DELTA_M);
  double sw; getValueWithDefault<std::string, double>(input_parameters, "sin2ThetaW", sw, DEFVAL_SW);
  double cw = 1.0 - sw;
  double alpha; getValueWithDefault<std::string, double>(input_parameters, "alpha", alpha, DEFVAL_ALPHA);
  const double pi = 3.14159265358979323846;
  double e = sqrt(4*pi*alpha);
  double MZ; getValueWithDefault<std::string, double>(input_parameters, "MZ", MZ, DEFVAL_MZ);
  double MW; getValueWithDefault<std::string, double>(input_parameters, "MW", MW, DEFVAL_MW);
  // Get the translation matrix
  std::vector<std::vector<double>> tmatrix = getTranslationMatrix(basis_input, basis_output, input_flags, input_parameters);
  // Get the input vector
  std::vector< std::pair<double, double> > vinput = getOrderedInputCouplings(basis_input, input_flags, input_parameters, input_couplings);
  // Assign the output vector
  std::vector< std::pair<double, double> > voutput(tmatrix.size(), std::pair<double, double>(0, 0));

  // Do the matrix multiplication
  for (size_t i=0; i<tmatrix.size(); i++){
    auto const& trow = tmatrix.at(i);
    for (size_t j=0; j<trow.size(); j++){
      voutput.at(i).first += trow.at(j) * vinput.at(j).first;
      voutput.at(i).second += trow.at(j) * vinput.at(j).second;
    }
  }
  // Fix output by subracting by constant vectors

  // Offsets for input Amplitude JHUGen
  if (basis_input == bAmplitude_JHUGen){
    if(basis_output == bHiggsBasis){
      voutput.at(coupl_hbasis_dCz).first -= 1.0;
      voutput.at(coupl_hbasis_dCw).first -= 1.0;
    }
  }
  // Offsets for EFT_JHUGen inputs
  else if (basis_input == bEFT_JHUGen){
    if (basis_output == bAmplitude_JHUGen){
      if (include_triple_quartic_gauge){
        voutput.at(coupl_ampjhutrip_dV_Z).first += 1.0;
        voutput.at(coupl_ampjhutrip_dV_A).first += 1.0;
        voutput.at(coupl_ampjhutrip_dP_Z).first += 1.0;
        voutput.at(coupl_ampjhutrip_dM_Z).first += 1.0;
        voutput.at(coupl_ampjhutrip_dZZWpWm).first += cw/sw;
        voutput.at(coupl_ampjhutrip_dZAWpWm).first += sqrt(cw)/sqrt(sw);
        // forcing the d2 gamma d3.0 gamma to be 1.0
        voutput.at(coupl_ampjhutrip_dP_A).first = 1.0;
        voutput.at(coupl_ampjhutrip_dM_A).first = 1.0;
        voutput.at(coupl_ampjhutrip_dAAWpWm).first = 1.0;
      }
    }
    if (basis_output == bEFT_JHUGen){
      if(include_triple_quartic_gauge){
        voutput.at(coupl_eftjhutrip_dV_Z).first += 1.0;
        voutput.at(coupl_eftjhutrip_dV_A).first += 1.0;
        voutput.at(coupl_eftjhutrip_dP_Z).first += 1.0;
        voutput.at(coupl_eftjhutrip_dM_Z).first += 1.0;
        voutput.at(coupl_eftjhutrip_dZZWpWm).first += cw/sw;
        voutput.at(coupl_eftjhutrip_dZAWpWm).first += sqrt(cw)/sqrt(sw);
        // forcing the d2 gamma d3.0 gamma to be 1.0 as well as dAAWpWm
        voutput.at(coupl_eftjhutrip_dP_A).first = 1.0;
        voutput.at(coupl_eftjhutrip_dM_A).first = 1.0;
        voutput.at(coupl_eftjhutrip_dAAWpWm).first = 1.0;         
      }
    }
    if (basis_output == bHiggsBasis){
      if (custodial_symmetry){
        if(include_triple_quartic_gauge){
          voutput.at(coupl_hbasistrip_dCz).first -= 1.0;
          voutput.at(coupl_hbasistrip_dCw).first -= 1.0;
        }
        else{
          voutput.at(coupl_hbasis_dCz).first -= 1.0;
          voutput.at(coupl_hbasis_dCw).first -= 1.0;
        }
      }
      else{
        if(include_triple_quartic_gauge){
          voutput.at(coupl_hbasistrip_dCz).first -= 1.0;
          voutput.at(coupl_hbasistrip_dCw).first -= 1.0 -4.0 * delta_m;
        }
        else{
          voutput.at(coupl_hbasis_dCz).first -= 1.0;
          voutput.at(coupl_hbasis_dCw).first -= 1.0 -4.0 * delta_m;
        }
      }
    }
    if (basis_output == bEFT_HiggsBasis){
      if(include_triple_quartic_gauge){
        voutput.at(coupl_efthbasistrip_dCz).first -= 1.0;
      }
      else{
        voutput.at(coupl_efthbasis_dCz).first -= 1.0;
      }
    }
    if (basis_output == bWarsawBasis){
      if (custodial_symmetry){
        voutput.at(coupl_warsaw_cHbx).first -= 1.0/(vev_lam);   
      }
      else{
        voutput.at(coupl_warsaw_cHbx).first += 1.0/(vev_lam)* ((2.0 * delta_m) - 1.0);
        voutput.at(coupl_warsaw_cHD).first -= 1.0/(vev_lam)* (4.0 * delta_m);
      }
    }
  }
  // Offsets for Higgs Basis inputs
  else if (basis_input == bHiggsBasis && basis_output == bAmplitude_JHUGen){
    voutput.at(coupl_ampjhu_ghz1).first += 2;
    voutput.at(coupl_ampjhu_ghw1).first += 2;
  }
  // Offsets for EFT Higgs Basis inputs
  else if (basis_input ==bEFT_HiggsBasis){
    if (basis_output == bAmplitude_JHUGen){
      if (custodial_symmetry){
        if (include_triple_quartic_gauge){
          voutput.at(coupl_ampjhutrip_ghz1).first += 2.0;
          voutput.at(coupl_ampjhutrip_ghw1).first += 2.0;
          voutput.at(coupl_ampjhutrip_dV_Z).first += 1.0;
          voutput.at(coupl_ampjhutrip_dV_A).first += 1.0;
          voutput.at(coupl_ampjhutrip_dP_Z).first += 1.0;
          voutput.at(coupl_ampjhutrip_dM_Z).first += 1.0;
          voutput.at(coupl_ampjhutrip_dZZWpWm).first += cw/sw;
          voutput.at(coupl_ampjhutrip_dZAWpWm).first += sqrt(cw)/sqrt(sw);
          voutput.at(coupl_ampjhutrip_dP_A).first = 1.0;
          voutput.at(coupl_ampjhutrip_dM_A).first = 1.0;
          voutput.at(coupl_ampjhutrip_dAAWpWm).first = 1.0;
        }
        else{
          voutput.at(coupl_ampjhu_ghz1).first += 2;
          voutput.at(coupl_ampjhu_ghw1).first += 2;
        }
      }
      else{
        if (include_triple_quartic_gauge){
          voutput.at(coupl_ampjhutrip_ghz1).first += 2.0;
          voutput.at(coupl_ampjhutrip_ghw1).first += 2.0 + 8.0*delta_m;
          voutput.at(coupl_ampjhutrip_dV_Z).first += 1.0;
          voutput.at(coupl_ampjhutrip_dV_A).first += 1.0;
          voutput.at(coupl_ampjhutrip_dP_Z).first += 1.0;
          voutput.at(coupl_ampjhutrip_dM_Z).first += 1.0;
          voutput.at(coupl_ampjhutrip_dZZWpWm).first += cw/sw;
          voutput.at(coupl_ampjhutrip_dZAWpWm).first += sqrt(cw)/sqrt(sw);
          voutput.at(coupl_ampjhutrip_dP_A).first = 1.0;
          voutput.at(coupl_ampjhutrip_dM_A).first = 1.0;
          voutput.at(coupl_ampjhutrip_dAAWpWm).first = 1.0;
        }
        else{
          voutput.at(coupl_ampjhu_ghz1).first += 2;
          voutput.at(coupl_ampjhu_ghw1).first += 2 + 8.0*delta_m;
        } 
      }   
    }
    if (basis_output == bEFT_JHUGen){
      if (include_triple_quartic_gauge){
        voutput.at(coupl_eftjhutrip_ghz1).first += 2.0;
        voutput.at(coupl_eftjhutrip_dV_Z).first += 1.0;
        voutput.at(coupl_eftjhutrip_dV_A).first += 1.0;
        voutput.at(coupl_eftjhutrip_dP_Z).first += 1.0;
        voutput.at(coupl_eftjhutrip_dM_Z).first += 1.0;
        voutput.at(coupl_eftjhutrip_dZZWpWm).first += cw/sw;
        voutput.at(coupl_eftjhutrip_dZAWpWm).first += sqrt(cw)/sqrt(sw);
        voutput.at(coupl_eftjhutrip_dP_A).first = 1.0;
        voutput.at(coupl_eftjhutrip_dM_A).first = 1.0;
        voutput.at(coupl_eftjhutrip_dAAWpWm).first = 1.0;
      }
      else{
	voutput.at(coupl_eftjhu_ghz1).first += 2;	
      }
    }
    if (basis_output == bWarsawBasis){
      if (!custodial_symmetry){
 	voutput.at(coupl_warsaw_cHbx).first += (1.0/vev_lam)* 2.0 * delta_m;
        voutput.at(coupl_warsaw_cHD).first += -(1.0/vev_lam)* 4.0 * delta_m;
      }
    }
  }
  // Offsets for Warsaw Basis input
  else if (basis_input == bWarsawBasis){
    if (custodial_symmetry){
      if (basis_output == bAmplitude_JHUGen){
        if (include_triple_quartic_gauge){
            voutput.at(coupl_ampjhutrip_ghz1).first += 2.0;
            voutput.at(coupl_ampjhutrip_ghw1).first += 2.0;
            voutput.at(coupl_ampjhutrip_dV_Z).first += 1.0;
            voutput.at(coupl_ampjhutrip_dV_A).first += 1.0;
            voutput.at(coupl_ampjhutrip_dP_Z).first += 1.0;
            voutput.at(coupl_ampjhutrip_dM_Z).first += 1.0;
            voutput.at(coupl_ampjhutrip_dZZWpWm).first += cw/sw;
            voutput.at(coupl_ampjhutrip_dZAWpWm).first += sqrt(cw)/sqrt(sw);
            voutput.at(coupl_ampjhutrip_dP_A).first = 1.0;
            voutput.at(coupl_ampjhutrip_dM_A).first = 1.0;
            voutput.at(coupl_ampjhutrip_dAAWpWm).first = 1.0;
        }
        else {
          voutput.at(coupl_ampjhu_ghz1).first += 2.0;
          voutput.at(coupl_ampjhu_ghw1).first += 2.0;
        }
      }
      if (basis_output==bEFT_JHUGen){
        if (include_triple_quartic_gauge){
          voutput.at(coupl_eftjhutrip_ghz1).first += 2.0;
          voutput.at(coupl_eftjhutrip_dV_Z).first += 1.0;
          voutput.at(coupl_eftjhutrip_dV_A).first += 1.0;
          voutput.at(coupl_eftjhutrip_dP_Z).first += 1.0;
          voutput.at(coupl_eftjhutrip_dM_Z).first += 1.0;
          voutput.at(coupl_eftjhutrip_dZZWpWm).first += cw/sw;
          voutput.at(coupl_eftjhutrip_dZAWpWm).first += sqrt(cw)/sqrt(sw);
          voutput.at(coupl_eftjhutrip_dP_A).first = 1.0;
          voutput.at(coupl_eftjhutrip_dM_A).first = 1.0;
          voutput.at(coupl_eftjhutrip_dAAWpWm).first = 1.0;
        }
        else{
          voutput.at(coupl_eftjhu_ghz1).first +=2.0;
        }
      }
    }
    else{
      if (basis_output == bHiggsBasis){
	if(include_triple_quartic_gauge){
	  voutput.at(coupl_hbasistrip_dCz).first += (vev_lam)*-3.0*delta_v;
	  voutput.at(coupl_hbasistrip_dCw).first += (vev_lam)*-(sw)/(cw-sw) *delta_v;
	  voutput.at(coupl_hbasistrip_Czbx).first += (vev_lam)*2.0*pow(e,2)/(sw)*delta_v; /// Edited!!!
	  voutput.at(coupl_hbasistrip_Cwbx).first += (vev_lam)*(pow(MW,2)*sw)/(2.0*pow(e,2)*pow(MZ,2)*(1-2.0*sw))*delta_v;
	  voutput.at(coupl_hbasistrip_Cabx).first += (vev_lam)*-(4.0*cw*sw)/(pow(e,2)*(cw-sw)) *delta_v;
	  voutput.at(coupl_hbasistrip_dKz).first += (vev_lam)*4.0/(16.0*sw-8.0) * delta_v;
	  voutput.at(coupl_hbasistrip_dg1z).first += (vev_lam)*-4.0/(4.0-8.0*sw) * delta_v;
	}
	else{
	  voutput.at(coupl_hbasis_dCz).first +=(vev_lam)* -3.0*delta_v;
          voutput.at(coupl_hbasis_dCw).first +=(vev_lam)* -(sw)/(cw-sw) *delta_v;
          voutput.at(coupl_hbasis_Czbx).first +=(vev_lam)* 4.0*(sw)/(2.0*pow(e,2))*delta_v;
  	  voutput.at(coupl_hbasis_Cwbx).first +=(vev_lam)* (pow(MW,2)*sw)/(2.0*pow(e,2)*pow(MZ,2)*(1-2.0*sw))*delta_v;
          voutput.at(coupl_hbasis_Cabx).first +=(vev_lam)* -(4.0*cw*sw)/(pow(e,2)*(cw-sw)) *delta_v;
	}
      }
      if (basis_output ==bEFT_HiggsBasis){
	if(include_triple_quartic_gauge){
          voutput.at(coupl_efthbasistrip_dCz).first +=(vev_lam)* -3.0*delta_v;
          voutput.at(coupl_efthbasistrip_Czbx).first +=(vev_lam)* 4.0*(sw)/(2.0*pow(e,2))*delta_v;
          voutput.at(coupl_efthbasistrip_dKz).first +=(vev_lam)* 4.0/(16.0*sw-8.0) * delta_v;
          voutput.at(coupl_efthbasistrip_dg1z).first +=(vev_lam)* -4.0/(4.0-8.0*sw) * delta_v;
        }
        else{
          voutput.at(coupl_efthbasis_dCz).first +=(vev_lam)* -3.0*delta_v;
          voutput.at(coupl_efthbasis_Czbx).first +=(vev_lam)* 2.0*(sw)/(pow(e,2))*delta_v;
        }
      }
      if (basis_output == bAmplitude_JHUGen){
        if (include_triple_quartic_gauge){
          voutput.at(coupl_ampjhutrip_ghz1).first += (vev_lam)*( -6.0 * delta_v + 2.0);
          voutput.at(coupl_ampjhutrip_ghw1).first += (2.0 -(vev_lam)*(sw)/(cw-sw) *delta_v) ;
          voutput.at(coupl_ampjhutrip_ghz1_prime2).first += (vev_lam)*2.0 *delta_v ;
          voutput.at(coupl_ampjhutrip_ghw1_prime2).first += (vev_lam)* pow(MW,2)/(2.0*pow(MZ,2)*(1.0-2.0*sw)) * delta_v;
          voutput.at(coupl_ampjhutrip_ghzgs1_prime2).first +=(vev_lam)* - (4.0*sqrt(sw)*sqrt(cw))/(cw-sw) * delta_v;
          voutput.at(coupl_ampjhutrip_dV_Z).first += 1.0 + (vev_lam)* 4.0/(16.0*sw-8.0) * delta_v;
          voutput.at(coupl_ampjhutrip_dV_A).first += 1.0;
          voutput.at(coupl_ampjhutrip_dP_Z).first += 1.0 - (vev_lam)*4.0/(4.0-8.0*sw) * delta_v;
          voutput.at(coupl_ampjhutrip_dM_Z).first += 1.0;
          voutput.at(coupl_ampjhutrip_dZZWpWm).first += cw/sw *(1.0 -(vev_lam)*4.0/(4.0-8.0*sw) * delta_v);
          voutput.at(coupl_ampjhutrip_dZAWpWm).first += sqrt(cw)/sqrt(sw) * (1.0 - (vev_lam)*4.0/(4.0-8.0*sw) * delta_v);
          voutput.at(coupl_ampjhutrip_dP_A).first = 1.0;
          voutput.at(coupl_ampjhutrip_dM_A).first = 1.0;
          voutput.at(coupl_ampjhutrip_dAAWpWm).first = 1.0;
        }
        else {
          voutput.at(coupl_ampjhu_ghz1).first += -(vev_lam)*6.0 * delta_v + 2.0;
          voutput.at(coupl_ampjhu_ghw1).first += 2.0 -(vev_lam)*(sw)/(cw-sw) *delta_v ;
          voutput.at(coupl_ampjhu_ghz1_prime2).first += (vev_lam)*2.0 *delta_v ;
          voutput.at(coupl_ampjhu_ghw1_prime2).first += (vev_lam)*pow(MW,2)/(2.0*pow(MZ,2)*(1.0-2.0*sw)) * delta_v;
          voutput.at(coupl_ampjhu_ghzgs1_prime2).first += - (vev_lam)*(4.0*sqrt(sw)*sqrt(cw))/(cw-sw) * delta_v;
    	}
      }
      if (basis_output==bEFT_JHUGen){
        if (include_triple_quartic_gauge){
          voutput.at(coupl_eftjhutrip_ghz1).first += -6.0 *(vev_lam)* delta_v + 2.0;
          voutput.at(coupl_eftjhutrip_ghz1_prime2).first += 2.0*(vev_lam) * delta_v ;
          voutput.at(coupl_eftjhutrip_dV_Z).first += (-4.0+8*sw)/(-4.0+8*sw) + (vev_lam)*4.0/(-4.0+8.0*sw) * delta_v;
          voutput.at(coupl_eftjhutrip_dV_A).first += 1.0;
          voutput.at(coupl_eftjhutrip_dP_Z).first += (-4.0+8*sw)/(-4.0+8.0*sw) + (vev_lam)*4.0/(-4.0+8.0*sw) * delta_v;
          voutput.at(coupl_eftjhutrip_dM_Z).first += (-4.0+8*sw)/(-4.0+8.0*sw) + (vev_lam)*4.0/(-4.0+8.0*sw) * delta_v;
          voutput.at(coupl_eftjhutrip_dZZWpWm).first += (cw)/(sw*(cw-sw)) - (2.0*cw)/(cw-sw) - vev_lam * (2.0*cw)/(sw*(cw-sw)) * delta_v ;
          voutput.at(coupl_eftjhutrip_dZAWpWm).first += (2*sqrt(cw)*sqrt(sw))/(sw-cw) - (4*sqrt(cw))/(4*sqrt(sw)*(sw-cw)) + vev_lam * sqrt(cw)/(sqrt(sw)*(sw-cw)) * delta_v;
          voutput.at(coupl_eftjhutrip_dP_A).first = 1.0;
          voutput.at(coupl_eftjhutrip_dM_A).first = 1.0;
          voutput.at(coupl_eftjhutrip_dAAWpWm).first = 1.0;
	}
        else{
          voutput.at(coupl_eftjhu_ghz1).first += -6.0 * (vev_lam)*delta_v + 2.0;
          voutput.at(coupl_eftjhu_ghz1_prime2).first += 2.0*(vev_lam)*delta_v ;
        }
      }
    }
  }
  // Format for either CC or NC output
  
  // Formatting for Charged Current
  if (HW_couplings_only){
    if (basis_output==bAmplitude_JHUGen){
      if(include_triple_quartic_gauge){
        std::vector< std::pair<double, double> > tempoutput(17, std::pair<double, double>(0, 0));
        tempoutput.at(coupl_ampjhutripcc_ghw1) = voutput.at(coupl_ampjhutrip_ghw1);
        tempoutput.at(coupl_ampjhutripcc_ghw1_prime2) = voutput.at(coupl_ampjhutrip_ghw1_prime2);
        tempoutput.at(coupl_ampjhutripcc_ghw2) = voutput.at(coupl_ampjhutrip_ghw2);
        tempoutput.at(coupl_ampjhutripcc_ghw4) = voutput.at(coupl_ampjhutrip_ghw4);
        tempoutput.at(coupl_ampjhutripcc_ghg2) = voutput.at(coupl_ampjhutrip_ghg2);
        tempoutput.at(coupl_ampjhutripcc_ghg4) = voutput.at(coupl_ampjhutrip_ghg4);
        tempoutput.at(coupl_ampjhutripcc_dV_Z) = voutput.at(coupl_ampjhutrip_dV_Z);
        tempoutput.at(coupl_ampjhutripcc_dV_A) = voutput.at(coupl_ampjhutrip_dV_A);
        tempoutput.at(coupl_ampjhutripcc_dP_Z) = voutput.at(coupl_ampjhutrip_dP_Z);
        tempoutput.at(coupl_ampjhutripcc_dP_A) = voutput.at(coupl_ampjhutrip_dP_A);
        tempoutput.at(coupl_ampjhutripcc_dM_Z) = voutput.at(coupl_ampjhutrip_dM_Z);
        tempoutput.at(coupl_ampjhutripcc_dM_A) = voutput.at(coupl_ampjhutrip_dM_A);
        tempoutput.at(coupl_ampjhutripcc_dFour_Z) = voutput.at(coupl_ampjhutrip_dFour_Z);
        tempoutput.at(coupl_ampjhutripcc_dFour_A) = voutput.at(coupl_ampjhutrip_dFour_A);
        tempoutput.at(coupl_ampjhutripcc_dZZWpWm) = voutput.at(coupl_ampjhutrip_dZZWpWm);
        tempoutput.at(coupl_ampjhutripcc_dZAWpWm) = voutput.at(coupl_ampjhutrip_dZAWpWm);          
        tempoutput.at(coupl_ampjhutripcc_dAAWpWm) = voutput.at(coupl_ampjhutrip_dAAWpWm);
        voutput = tempoutput;
      }
      else {
        std::vector< std::pair<double, double> > tempoutput(6, std::pair<double, double>(0, 0));
        tempoutput.at(coupl_ampjhucc_ghw1) = voutput.at(coupl_ampjhu_ghw1);
        tempoutput.at(coupl_ampjhucc_ghw1_prime2) = voutput.at(coupl_ampjhu_ghw1_prime2);
        tempoutput.at(coupl_ampjhucc_ghw2) = voutput.at(coupl_ampjhu_ghw2);
        tempoutput.at(coupl_ampjhucc_ghw4) = voutput.at(coupl_ampjhu_ghw4);
        tempoutput.at(coupl_ampjhucc_ghg2) = voutput.at(coupl_ampjhutrip_ghg2);
        tempoutput.at(coupl_ampjhucc_ghg4) = voutput.at(coupl_ampjhutrip_ghg4);
        voutput = tempoutput;
      }
    }
  }
  else if (HZ_couplings_only){
    if (basis_output==bAmplitude_JHUGen){
      if(include_triple_quartic_gauge){
        std::vector< std::pair<double, double> > tempoutput(23, std::pair<double, double>(0, 0));
        tempoutput.at(coupl_ampjhutripnc_ghz1) = voutput.at(coupl_ampjhutrip_ghz1);
        tempoutput.at(coupl_ampjhutripnc_ghz1_prime2) = voutput.at(coupl_ampjhutrip_ghz1_prime2);
        tempoutput.at(coupl_ampjhutripnc_ghz2) = voutput.at(coupl_ampjhutrip_ghz2);
        tempoutput.at(coupl_ampjhutripnc_ghz4) = voutput.at(coupl_ampjhutrip_ghz4);
        tempoutput.at(coupl_ampjhutripnc_ghz4) = voutput.at(coupl_ampjhutrip_ghz4);
        tempoutput.at(coupl_ampjhutripnc_ghzgs1_prime2) = voutput.at(coupl_ampjhutrip_ghzgs1_prime2);
        tempoutput.at(coupl_ampjhutripnc_ghzgs2) = voutput.at(coupl_ampjhutrip_ghzgs2);
        tempoutput.at(coupl_ampjhutripnc_ghzgs4) = voutput.at(coupl_ampjhutrip_ghzgs4);
        tempoutput.at(coupl_ampjhutripnc_ghgsgs2) = voutput.at(coupl_ampjhutrip_ghgsgs2);
        tempoutput.at(coupl_ampjhutripnc_ghgsgs4) = voutput.at(coupl_ampjhutrip_ghgsgs4);
        tempoutput.at(coupl_ampjhutripnc_ghg2) = voutput.at(coupl_ampjhutrip_ghg2);
        tempoutput.at(coupl_ampjhutripnc_ghg4) = voutput.at(coupl_ampjhutrip_ghg4);
        tempoutput.at(coupl_ampjhutripnc_dV_Z) = voutput.at(coupl_ampjhutrip_dV_Z);
        tempoutput.at(coupl_ampjhutripnc_dV_A) = voutput.at(coupl_ampjhutrip_dV_A);
        tempoutput.at(coupl_ampjhutripnc_dP_Z) = voutput.at(coupl_ampjhutrip_dP_Z);
        tempoutput.at(coupl_ampjhutripnc_dP_A) = voutput.at(coupl_ampjhutrip_dP_A);
        tempoutput.at(coupl_ampjhutripnc_dM_Z) = voutput.at(coupl_ampjhutrip_dM_Z);
        tempoutput.at(coupl_ampjhutripnc_dM_A) = voutput.at(coupl_ampjhutrip_dM_A);
        tempoutput.at(coupl_ampjhutripnc_dFour_Z) = voutput.at(coupl_ampjhutrip_dFour_Z);
        tempoutput.at(coupl_ampjhutripnc_dFour_A) = voutput.at(coupl_ampjhutrip_dFour_A);
        tempoutput.at(coupl_ampjhutripnc_dZZWpWm) = voutput.at(coupl_ampjhutrip_dZZWpWm);
        tempoutput.at(coupl_ampjhutripnc_dZAWpWm) = voutput.at(coupl_ampjhutrip_dZAWpWm);
        tempoutput.at(coupl_ampjhutripnc_dAAWpWm) = voutput.at(coupl_ampjhutrip_dAAWpWm);
        voutput = tempoutput;
      }
      else {
        std::vector< std::pair<double, double> > tempoutput(12, std::pair<double, double>(0, 0));
        tempoutput.at(coupl_ampjhunc_ghz1) = voutput.at(coupl_ampjhutrip_ghz1);
        tempoutput.at(coupl_ampjhunc_ghz1_prime2) = voutput.at(coupl_ampjhutrip_ghz1_prime2);
        tempoutput.at(coupl_ampjhunc_ghz2) = voutput.at(coupl_ampjhutrip_ghz2);
        tempoutput.at(coupl_ampjhunc_ghz4) = voutput.at(coupl_ampjhutrip_ghz4);
        tempoutput.at(coupl_ampjhunc_ghz4) = voutput.at(coupl_ampjhutrip_ghz4);
        tempoutput.at(coupl_ampjhunc_ghzgs1_prime2) = voutput.at(coupl_ampjhutrip_ghzgs1_prime2);
        tempoutput.at(coupl_ampjhunc_ghzgs2) = voutput.at(coupl_ampjhutrip_ghzgs2);
        tempoutput.at(coupl_ampjhunc_ghzgs4) = voutput.at(coupl_ampjhutrip_ghzgs4);
        tempoutput.at(coupl_ampjhunc_ghgsgs2) = voutput.at(coupl_ampjhutrip_ghgsgs2);
        tempoutput.at(coupl_ampjhunc_ghgsgs4) = voutput.at(coupl_ampjhutrip_ghgsgs4);
        tempoutput.at(coupl_ampjhunc_ghg2) = voutput.at(coupl_ampjhutrip_ghg2);
        tempoutput.at(coupl_ampjhunc_ghg4) = voutput.at(coupl_ampjhutrip_ghg4);
        voutput = tempoutput;
      }
    }
  } 
  // Set the results
  interpretOutputCouplings(basis_output, input_flags, input_parameters, voutput);
}


std::vector<std::vector<double>> JHUGenLexiconTranslator::getTranslationMatrix(
  JHUGenLexiconIOHelpers::IOBasisType const& basis_input, JHUGenLexiconIOHelpers::IOBasisType const& basis_output,
  std::unordered_map<std::string, bool> const& input_flags,
  std::unordered_map<std::string, double> const& input_parameters
) const{
  bool include_triple_quartic_gauge; getValueWithDefault<std::string, bool>(input_flags, "include_triple_quartic_gauge", include_triple_quartic_gauge, false);
  bool custodial_symmetry; getValueWithDefault<std::string, bool>(input_flags, "custodial_symmetry", custodial_symmetry, false);
  double alpha; getValueWithDefault<std::string, double>(input_parameters, "alpha", alpha, DEFVAL_ALPHA);
  double sw; getValueWithDefault<std::string, double>(input_parameters, "sin2ThetaW", sw, DEFVAL_SW);
  double vev_lam; getValueWithDefault<std::string, double>(input_parameters, "vev_lam", vev_lam, DEFVAL_VEV_LAM);
  const double pi = 3.14159265358979323846;
  double e = sqrt(4*pi*alpha);
  double cw = 1.0-sw;
  double MZ; getValueWithDefault<std::string, double>(input_parameters, "MZ", MZ, DEFVAL_MZ);
  double MW; getValueWithDefault<std::string, double>(input_parameters, "MW", MW, DEFVAL_MW);
  double Lambda_z1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_z1", Lambda_z1, DEFVAL_LAMBDA_VI);
  double Lambda_w1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_w1", Lambda_w1, DEFVAL_LAMBDA_VI);
  double Lambda_zgs1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_zgs1", Lambda_zgs1, DEFVAL_LAMBDA_VI);
  std::vector<std::vector<double>> res;

  // This is where the translation matrixes should be coded carefully.
  // See also the getOrderedInputCouplings and interpretOutputCouplings functions for what the input/output vectors expect for JHUGen conventions (e.g. ghz1_prime2 scaled already by MZ^2/L1.0ZZ^2, so no need to put that for example).
  if (basis_input == bAmplitude_JHUGen){
    if (basis_output == bAmplitude_JHUGen){
      res.assign(nAmplitude_JHUGen_CouplingTypes, std::vector<double>(nAmplitude_JHUGen_CouplingTypes, 0));//Assign Correct Dimensions of translation matrix
      for (size_t i=0; i<(size_t) nAmplitude_JHUGen_CouplingTypes; i++) res.at(i).at(i)=1.0;
    }
    else if (basis_output == bEFT_JHUGen){
      res.assign(nEFT_JHUGen_CouplingTypes, std::vector<double>(nAmplitude_JHUGen_CouplingTypes, 0));
    }
    else if (basis_output == bHiggsBasis){
      res.assign(nHiggsBasis_CouplingTypes, std::vector<double>(nAmplitude_JHUGen_CouplingTypes, 0));
      res[coupl_hbasis_dCz][coupl_ampjhu_ghz1]=.5;
      res[coupl_hbasis_Czz][coupl_ampjhu_ghz2]=-(2.0*cw*sw)/pow(e,2);
      res[coupl_hbasis_Czbx][coupl_ampjhu_ghz1_prime2]= (sw)/(pow(e,2));
      res[coupl_hbasis_tCzz][coupl_ampjhu_ghz4]= (-2*cw*sw)/pow(e,2);
      res[coupl_hbasis_dCw][coupl_ampjhu_ghw1]= .5;
      res[coupl_hbasis_Cww][coupl_ampjhu_ghw2]= (-2.0*sw)/pow(e,2);
      res[coupl_hbasis_Cwbx][coupl_ampjhu_ghw1_prime2]= (sw)/pow(e,2);
      res[coupl_hbasis_tCww][coupl_ampjhu_ghw4]= (-2.0*sw)/pow(e,2);
      res[coupl_hbasis_Cza][coupl_ampjhu_ghzgs2]= (-2.0*sqrt(sw)*sqrt(cw))/pow(e,2);
      res[coupl_hbasis_tCza][coupl_ampjhu_ghzgs4]= (-2.0*sqrt(sw)*sqrt(cw))/pow(e,2);
      res[coupl_hbasis_Cabx][coupl_ampjhu_ghzgs1_prime2]= (sqrt(cw)*sqrt(sw))/pow(e,2);
      res[coupl_hbasis_Caa][coupl_ampjhu_ghgsgs2]= -2.0/pow(e,2);
      res[coupl_hbasis_tCaa][coupl_ampjhu_ghgsgs4]= -2.0/pow(e,2);
      res[coupl_hbasis_Cgg][coupl_ampjhu_ghg2]= -2.0/1.0;
      res[coupl_hbasis_tCgg][coupl_ampjhu_ghg4]= -2.0/1.0;
    }
    else if (basis_output == bEFT_HiggsBasis){
      res.assign(nEFT_HiggsBasis_CouplingTypes, std::vector<double>(nAmplitude_JHUGen_CouplingTypes, 0));
    }
  }
  else if (basis_input == bEFT_JHUGen){
    if (basis_output == bAmplitude_JHUGen){
      if (include_triple_quartic_gauge){
        res.assign(nAmplitude_JHUGen_Include_Triple_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
        res[coupl_ampjhutrip_ghz1][coupl_eftjhu_ghz1]= 1.0;
        res[coupl_ampjhutrip_ghz1_prime2][coupl_eftjhu_ghz1_prime2]= 1.0;
        res[coupl_ampjhutrip_ghz2][coupl_eftjhu_ghz2]= 1.0;
        res[coupl_ampjhutrip_ghz4][coupl_eftjhu_ghz4]= 1.0;
        res[coupl_ampjhutrip_ghw1][coupl_eftjhu_ghz1]= 1.0;
        res[coupl_ampjhutrip_ghw1_prime2][coupl_eftjhu_ghz1_prime2]= 1.0/pow(MZ,2) * pow(MW,2)/(cw - sw);
        res[coupl_ampjhutrip_ghw1_prime2][coupl_eftjhu_ghz2]= (-2*sw)/pow(MZ,2) * pow(MW,2)/(cw - sw);
        res[coupl_ampjhutrip_ghw1_prime2][coupl_eftjhu_ghzgs2]= (2*sqrt(sw))/sqrt(cw) *(cw - sw)/pow(MZ,2)*pow(MW,2)/(cw - sw);
        res[coupl_ampjhutrip_ghw1_prime2][coupl_eftjhu_ghgsgs2]= (2*sw)/pow(MZ,2)*pow(MW,2)/(cw - sw);
        res[coupl_ampjhutrip_ghw2][coupl_eftjhu_ghz2]= cw;
        res[coupl_ampjhutrip_ghw2][coupl_eftjhu_ghzgs2]= 2*sqrt(sw)*sqrt(cw);
        res[coupl_ampjhutrip_ghw2][coupl_eftjhu_ghgsgs2]= sw;
        res[coupl_ampjhutrip_ghw4][coupl_eftjhu_ghz4] = cw;
        res[coupl_ampjhutrip_ghw4][coupl_eftjhu_ghzgs4] = 2*sqrt(sw)*sqrt(cw);
        res[coupl_ampjhutrip_ghw4][coupl_eftjhu_ghgsgs4] = sw;
        res[coupl_ampjhutrip_ghzgs1_prime2][coupl_eftjhu_ghz1_prime2]= (2*sqrt(sw)*sqrt(cw))/(cw - sw);
        res[coupl_ampjhutrip_ghzgs1_prime2][coupl_eftjhu_ghz2]= (-2*sqrt(sw)*sqrt(cw))/(cw - sw);
        res[coupl_ampjhutrip_ghzgs1_prime2][coupl_eftjhu_ghzgs2]= 2;
        res[coupl_ampjhutrip_ghzgs1_prime2][coupl_eftjhu_ghgsgs2]= (2*sqrt(sw)*sqrt(cw))/(cw - sw);
        res[coupl_ampjhutrip_ghzgs2][coupl_eftjhu_ghzgs2]= 1.0;
        res[coupl_ampjhutrip_ghzgs4][coupl_eftjhu_ghzgs4]= 1.0;
        res[coupl_ampjhutrip_ghgsgs2][coupl_eftjhu_ghgsgs2]= 1.0;
        res[coupl_ampjhutrip_ghgsgs4][coupl_eftjhu_ghgsgs4]= 1.0;
        res[coupl_ampjhutrip_ghg2][coupl_eftjhu_ghg2]= 1.0;
        res[coupl_ampjhutrip_ghg4][coupl_eftjhu_ghg4]= 1.0;
        res[coupl_ampjhutrip_dV_Z][coupl_eftjhu_ghz1_prime2]= -1.0/(2.0*(cw - sw));
        res[coupl_ampjhutrip_dV_Z][coupl_eftjhu_ghz2]= (2.0*sw*cw)/(cw - sw);
        res[coupl_ampjhutrip_dV_Z][coupl_eftjhu_ghzgs2]= -2.0*sqrt(sw)*sqrt(cw);
        res[coupl_ampjhutrip_dV_Z][coupl_eftjhu_ghgsgs2]=  -(2.0*sw*cw)/(cw - sw); //check!
        res[coupl_ampjhutrip_dV_A][coupl_eftjhu_ghz2]= -cw;
        res[coupl_ampjhutrip_dV_A][coupl_eftjhu_ghzgs2]= (sqrt(cw)/sqrt(sw) - 2*sqrt(sw)*sqrt(cw));
        res[coupl_ampjhutrip_dV_A][coupl_eftjhu_ghgsgs2]= cw;
        res[coupl_ampjhutrip_dP_Z][coupl_eftjhu_ghz1_prime2]= -1.0/(2.0*(cw - sw));
        res[coupl_ampjhutrip_dP_Z][coupl_eftjhu_ghz2]= sw/(cw - sw);
        res[coupl_ampjhutrip_dP_Z][coupl_eftjhu_ghzgs2]= -sqrt(sw)/sqrt(cw);
        res[coupl_ampjhutrip_dP_Z][coupl_eftjhu_ghgsgs2]= -sw/(cw-sw);
        //res[coupl_ampjhutrip_dP_A][coupl_eftjhu_ghz1]= 0;  // dP_A is fixed by any relation and should be kept as 1.0 to preserve symmetry
        res[coupl_ampjhutrip_dM_Z][coupl_eftjhu_ghz1_prime2]= -1.0/(2.0*(cw - sw)) ;
        res[coupl_ampjhutrip_dM_Z][coupl_eftjhu_ghz2]= sw/(cw-sw);
        res[coupl_ampjhutrip_dM_Z][coupl_eftjhu_ghzgs2]= -sqrt(sw)/sqrt(cw);
        res[coupl_ampjhutrip_dM_Z][coupl_eftjhu_ghgsgs2]= -sw/(cw-sw);
        //res[coupl_ampjhutrip_dM_A][coupl_eftjhu_ghz1]= 0; // dM_A is fixed by any relation and should be kept as 1.0 to preserve symmetry
        res[coupl_ampjhutrip_dFour_Z][coupl_eftjhu_ghz4]= sw;
        res[coupl_ampjhutrip_dFour_Z][coupl_eftjhu_ghzgs4] = (-sqrt(sw)/sqrt(cw) + 2*pow(sw,3.0/2)/sqrt(cw));
        res[coupl_ampjhutrip_dFour_Z][coupl_eftjhu_ghgsgs4]= -sw;
        res[coupl_ampjhutrip_dFour_A][coupl_eftjhu_ghz4]= -cw;
        res[coupl_ampjhutrip_dFour_A][coupl_eftjhu_ghzgs4]= (sqrt(cw)/sqrt(sw) - 2*sqrt(sw)*sqrt(cw));
        res[coupl_ampjhutrip_dFour_A][coupl_eftjhu_ghgsgs4] = cw;
        res[coupl_ampjhutrip_dZZWpWm][coupl_eftjhu_ghz1_prime2]=  -(cw)/(sw*(cw - sw));
        res[coupl_ampjhutrip_dZZWpWm][coupl_eftjhu_ghz2]= (2.0*cw)/(cw-sw);
        res[coupl_ampjhutrip_dZZWpWm][coupl_eftjhu_ghzgs2]= (-2.0*sqrt(cw))/sqrt(sw);
        res[coupl_ampjhutrip_dZZWpWm][coupl_eftjhu_ghgsgs2]= -(2.0*cw)/(cw-sw);
        res[coupl_ampjhutrip_dZAWpWm][coupl_eftjhu_ghz1_prime2]= -(sqrt(cw))/(2*sqrt(sw)*(cw - sw));
        res[coupl_ampjhutrip_dZAWpWm][coupl_eftjhu_ghz2]= (sqrt(sw)*sqrt(cw))/(cw-sw);
        res[coupl_ampjhutrip_dZAWpWm][coupl_eftjhu_ghzgs2]= -1.0;
        res[coupl_ampjhutrip_dZAWpWm][coupl_eftjhu_ghgsgs2]=-(sqrt(sw)*sqrt(cw))/(cw-sw);
        //res[coupl_ampjhutrip_dAAWpWm][coupl_eftjhu_ghz1]= 0; ?
      }
      else{
        res.assign(nAmplitude_JHUGen_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
        res[coupl_ampjhu_ghz1][coupl_eftjhu_ghz1]= 1.0;
        res[coupl_ampjhu_ghz1_prime2][coupl_eftjhu_ghz1_prime2]= 1.0;
        res[coupl_ampjhu_ghz2][coupl_eftjhu_ghz2]= 1.0;
        res[coupl_ampjhu_ghz4][coupl_eftjhu_ghz4]= 1.0;
        res[coupl_ampjhu_ghw1][coupl_eftjhu_ghz1]= 1.0;
        res[coupl_ampjhu_ghw1_prime2][coupl_eftjhu_ghz1_prime2]= 1.0/pow(MZ,2) * pow(MW,2)/(cw - sw);
        res[coupl_ampjhu_ghw1_prime2][coupl_eftjhu_ghz2]= (-2*sw)/pow(MZ,2) * pow(MW,2)/(cw - sw);
        res[coupl_ampjhu_ghw1_prime2][coupl_eftjhu_ghzgs2]= (2*sqrt(sw))/sqrt(cw) *(cw - sw)/pow(MZ,2)*pow(MW,2)/(cw - sw);
        res[coupl_ampjhu_ghw1_prime2][coupl_eftjhu_ghgsgs2]= (2*sw)/pow(MZ,2)*pow(MW,2)/(cw - sw);
        res[coupl_ampjhu_ghw2][coupl_eftjhu_ghz2]= cw;
        res[coupl_ampjhu_ghw2][coupl_eftjhu_ghzgs2]= 2*sqrt(sw)*sqrt(cw);
        res[coupl_ampjhu_ghw2][coupl_eftjhu_ghgsgs2]= sw;
        res[coupl_ampjhu_ghw4][coupl_eftjhu_ghz4] = cw;
        res[coupl_ampjhu_ghw4][coupl_eftjhu_ghzgs4] = 2*sqrt(sw)*sqrt(cw);
        res[coupl_ampjhu_ghw4][coupl_eftjhu_ghgsgs4] = sw;
        res[coupl_ampjhu_ghzgs1_prime2][coupl_eftjhu_ghz1_prime2]= (2*sqrt(sw)*sqrt(cw))/(cw - sw);
        res[coupl_ampjhu_ghzgs1_prime2][coupl_eftjhu_ghz2]= (-2*sqrt(sw)*sqrt(cw))/(cw - sw);
        res[coupl_ampjhu_ghzgs1_prime2][coupl_eftjhu_ghzgs2]= 2;
        res[coupl_ampjhu_ghzgs1_prime2][coupl_eftjhu_ghgsgs2]= (2*sqrt(sw)*sqrt(cw))/(cw - sw);
        res[coupl_ampjhu_ghzgs2][coupl_eftjhu_ghzgs2]= 1.0;
        res[coupl_ampjhu_ghzgs4][coupl_eftjhu_ghzgs4]= 1.0;
        res[coupl_ampjhu_ghgsgs2][coupl_eftjhu_ghgsgs2]= 1.0;
        res[coupl_ampjhu_ghgsgs4][coupl_eftjhu_ghgsgs4]= 1.0;
        res[coupl_ampjhu_ghg2][coupl_eftjhu_ghg2]= 1.0;
        res[coupl_ampjhu_ghg4][coupl_eftjhu_ghg4]= 1.0;
      }
    }
    else if (basis_output == bEFT_JHUGen){
      if (include_triple_quartic_gauge){
        res.assign(nEFT_JHUGen_Include_Triple_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
        res[coupl_eftjhutrip_ghz1][coupl_eftjhu_ghz1]=1.0;
        res[coupl_eftjhutrip_ghz1_prime2][coupl_eftjhu_ghz1_prime2]=1.0;
        res[coupl_eftjhutrip_ghz2][coupl_eftjhu_ghz2]=1.0;
        res[coupl_eftjhutrip_ghz4][coupl_eftjhu_ghz4]=1.0;
        res[coupl_eftjhutrip_ghzgs2][coupl_eftjhu_ghzgs2]=1.0;
        res[coupl_eftjhutrip_ghzgs4][coupl_eftjhu_ghzgs4]=1.0;
        res[coupl_eftjhutrip_ghgsgs2][coupl_eftjhu_ghgsgs2]=1.0;
        res[coupl_eftjhutrip_ghgsgs4][coupl_eftjhu_ghgsgs4]=1.0;
        res[coupl_eftjhutrip_ghg2][coupl_eftjhu_ghg2]=1.0;
        res[coupl_eftjhutrip_ghg4][coupl_eftjhu_ghg4]=1.0;
        res[coupl_eftjhutrip_dV_Z][coupl_eftjhu_ghz1_prime2]= -1.0/(2.0*(cw - sw));
        res[coupl_eftjhutrip_dV_Z][coupl_eftjhu_ghz2]= (2.0*sw*cw)/(cw - sw);
        res[coupl_eftjhutrip_dV_Z][coupl_eftjhu_ghzgs2]= -2.0*sqrt(sw)*sqrt(cw);
        res[coupl_eftjhutrip_dV_Z][coupl_eftjhu_ghgsgs2]=  -(2.0*sw*cw)/(cw - sw); //check!
        res[coupl_eftjhutrip_dV_A][coupl_eftjhu_ghz2]= -cw;
        res[coupl_eftjhutrip_dV_A][coupl_eftjhu_ghzgs2]= (sqrt(cw)/sqrt(sw) - 2*sqrt(sw)*sqrt(cw));
        res[coupl_eftjhutrip_dV_A][coupl_eftjhu_ghgsgs2]= cw;
        res[coupl_eftjhutrip_dP_Z][coupl_eftjhu_ghz1_prime2]= -1.0/(2.0*(cw - sw));
        res[coupl_eftjhutrip_dP_Z][coupl_eftjhu_ghz2]= sw/(cw - sw);
        res[coupl_eftjhutrip_dP_Z][coupl_eftjhu_ghzgs2]= -sqrt(sw)/sqrt(cw);
        res[coupl_eftjhutrip_dP_Z][coupl_eftjhu_ghgsgs2]= -sw/(cw-sw);
        //res[coupl_ampjhutrip_dP_A][coupl_eftjhu_ghz1]= 0;  // dP_A is fixed by any relation and should be kept as 1.0 to preserve symmetry
        res[coupl_eftjhutrip_dM_Z][coupl_eftjhu_ghz1_prime2]= -1.0/(2.0*(cw - sw)) ;
        res[coupl_eftjhutrip_dM_Z][coupl_eftjhu_ghz2]= sw/(cw-sw);
        res[coupl_eftjhutrip_dM_Z][coupl_eftjhu_ghzgs2]= -sqrt(sw)/sqrt(cw);
        res[coupl_eftjhutrip_dM_Z][coupl_eftjhu_ghgsgs2]= -sw/(cw-sw);
        //res[coupl_ampjhutrip_dM_A][coupl_eftjhu_ghz1]= 0; // dM_A is fixed by any relation and should be kept as 1.0 to preserve symmetry
        res[coupl_eftjhutrip_dFour_Z][coupl_eftjhu_ghz4]= sw;
        res[coupl_eftjhutrip_dFour_Z][coupl_eftjhu_ghzgs4] = (-sqrt(sw)/sqrt(cw) + 2*pow(sw,3.0/2)/sqrt(cw));
        res[coupl_eftjhutrip_dFour_Z][coupl_eftjhu_ghgsgs4]= -sw;
        res[coupl_eftjhutrip_dFour_A][coupl_eftjhu_ghz4]= -cw;
        res[coupl_eftjhutrip_dFour_A][coupl_eftjhu_ghzgs4]= (sqrt(cw)/sqrt(sw) - 2*sqrt(sw)*sqrt(cw));
        res[coupl_eftjhutrip_dFour_A][coupl_eftjhu_ghgsgs4] = cw;
        res[coupl_eftjhutrip_dZZWpWm][coupl_eftjhu_ghz1_prime2]=  -(cw)/(sw*(cw - sw));
        res[coupl_eftjhutrip_dZZWpWm][coupl_eftjhu_ghz2]=  (2.0*cw)/(cw-sw);
        res[coupl_eftjhutrip_dZZWpWm][coupl_eftjhu_ghzgs2]= (-2.0*sqrt(cw))/sqrt(sw);
        res[coupl_eftjhutrip_dZZWpWm][coupl_eftjhu_ghgsgs2]= -(2.0*cw)/(cw-sw);
        res[coupl_eftjhutrip_dZAWpWm][coupl_eftjhu_ghz1_prime2]= -(sqrt(cw))/(2*sqrt(sw)*(cw - sw));
        res[coupl_eftjhutrip_dZAWpWm][coupl_eftjhu_ghz2]= (sqrt(sw)*sqrt(cw))/(cw-sw);
        res[coupl_eftjhutrip_dZAWpWm][coupl_eftjhu_ghzgs2]= -1.0;
        res[coupl_eftjhutrip_dZAWpWm][coupl_eftjhu_ghgsgs2]=-(sqrt(sw)*sqrt(cw))/(cw-sw);
        //res[coupl_ampjhutrip_dAAWpWm][coupl_eftjhu_ghz1]= 0;
        cout << "here" ;
      }
      else{
        res.assign(nEFT_JHUGen_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
        for (size_t i=0; i<(size_t) nEFT_JHUGen_CouplingTypes; i++) res.at(i).at(i)=1.0;
      }
    }
    else if (basis_output == bHiggsBasis){
      if (include_triple_quartic_gauge){
        res.assign(nHiggsBasis_Include_Triple_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
        res[coupl_hbasistrip_dCz][coupl_eftjhu_ghz1] = .5;
        res[coupl_hbasistrip_Czz][coupl_eftjhu_ghz2] = -((2*cw*sw)/pow(e,2));
        res[coupl_hbasistrip_Czbx][coupl_eftjhu_ghz1_prime2] = sw/pow(e,2);
        res[coupl_hbasistrip_tCzz][coupl_eftjhu_ghz4] = -((2*cw*sw)/pow(e,2));
        res[coupl_hbasistrip_dCw][coupl_eftjhu_ghz1] = .5;
        res[coupl_hbasistrip_Cww][coupl_eftjhu_ghz2] =  -(2*cw*sw)/pow(e,2);
        res[coupl_hbasistrip_Cww][coupl_eftjhu_ghzgs2] = -((4*sqrt(cw)*pow(sw,3.0/2))/pow(e,2));
        res[coupl_hbasistrip_Cww][coupl_eftjhu_ghgsgs2] =  -((2*pow(sw,2))/pow(e,2));
        res[coupl_hbasistrip_Cwbx][coupl_eftjhu_ghz1_prime2] = (pow(MW,2)*sw)/(pow(e,2)*pow(MZ,2)*(cw - sw));
        res[coupl_hbasistrip_Cwbx][coupl_eftjhu_ghz2] = -((2*pow(MW,2)*pow(sw,2))/(pow(e,2)*pow(MZ,2)*(cw - sw)));
        res[coupl_hbasistrip_Cwbx][coupl_eftjhu_ghzgs2] = ((2*pow(MW,2)*pow(sw,3.0/2.0))/(sqrt(cw)*pow(e,2)*pow(MZ,2)));
        res[coupl_hbasistrip_Cwbx][coupl_eftjhu_ghgsgs2] = (2*pow(MW,2)*pow(sw,2))/(pow(e,2)*pow(MZ,2)*(cw - sw));
        res[coupl_hbasistrip_tCww][coupl_eftjhu_ghz4] = -((2*cw*sw)/pow(e,2));
        res[coupl_hbasistrip_tCww][coupl_eftjhu_ghzgs4] = -((4*sqrt(cw)*pow(sw,3.0/2))/pow(e,2));
        res[coupl_hbasistrip_tCww][coupl_eftjhu_ghgsgs4] =  -((2*pow(sw,2))/pow(e,2));
        res[coupl_hbasistrip_Cza][coupl_eftjhu_ghzgs2] = -((2*sqrt(cw)*sqrt(sw))/pow(e,2));
        res[coupl_hbasistrip_tCza][coupl_eftjhu_ghzgs4] = -((2*sqrt(cw)*sqrt(sw))/pow(e,2));
        res[coupl_hbasistrip_Cabx][coupl_eftjhu_ghz1_prime2] = (2*cw*sw)/(pow(e,2)*(cw - sw));
        res[coupl_hbasistrip_Cabx][coupl_eftjhu_ghz2] = -((2*cw*sw)/(pow(e,2)*(cw - sw)));
        res[coupl_hbasistrip_Cabx][coupl_eftjhu_ghzgs2] = (2*sqrt(sw)*sqrt(cw))/pow(e,2);
        res[coupl_hbasistrip_Cabx][coupl_eftjhu_ghgsgs2] = (2*cw*sw)/(pow(e,2)*(cw - sw));
        res[coupl_hbasistrip_Caa][coupl_eftjhu_ghgsgs2] = -(2/pow(e,2));
        res[coupl_hbasistrip_tCaa][coupl_eftjhu_ghgsgs4] = -(2/pow(e,2));
        res[coupl_hbasistrip_Cgg][coupl_eftjhu_ghg2] = -2.0/1.0;
        res[coupl_hbasistrip_tCgg][coupl_eftjhu_ghg4] = -2.0/1.0;
        res[coupl_hbasistrip_dKa][coupl_eftjhu_ghz2] = -cw/2;
        res[coupl_hbasistrip_dKa][coupl_eftjhu_ghzgs2] = (1.0/2.0)*(sqrt(cw)*(1.0-2.0*sw))/(sqrt(sw));
        res[coupl_hbasistrip_dKa][coupl_eftjhu_ghgsgs2] = cw/2;
        res[coupl_hbasistrip_tKa][coupl_eftjhu_ghz4] = -cw/2;
        res[coupl_hbasistrip_tKa][coupl_eftjhu_ghzgs4] = 1.0/2*(sqrt(cw)/sqrt(sw) - 2*sqrt(cw)*sqrt(sw));
        res[coupl_hbasistrip_tKa][coupl_eftjhu_ghgsgs4] = cw/2;
        res[coupl_hbasistrip_dKz][coupl_eftjhu_ghz1_prime2] = -(1.0/(4*(cw - sw)));
        res[coupl_hbasistrip_dKz][coupl_eftjhu_ghz2] = (cw*sw)/(cw - sw);
        res[coupl_hbasistrip_dKz][coupl_eftjhu_ghzgs2] = -sqrt(cw)*sqrt(sw);
        res[coupl_hbasistrip_dKz][coupl_eftjhu_ghgsgs2] = -((cw*sw)/(cw - sw));
        res[coupl_hbasistrip_tKz][coupl_eftjhu_ghz4] = sw/2;
        res[coupl_hbasistrip_tKz][coupl_eftjhu_ghzgs4] = 1.0/2 *(-(sqrt(sw)/sqrt(cw)) + (2*pow(sw,3.0/2))/sqrt(cw));
        res[coupl_hbasistrip_tKz][coupl_eftjhu_ghgsgs4] = -sw/2;
        res[coupl_hbasistrip_dg1z][coupl_eftjhu_ghz1_prime2] = (1.0/(2.0*(2.0*sw-1.0)));
        res[coupl_hbasistrip_dg1z][coupl_eftjhu_ghz2] = -sw/(2.0*sw-1.0);
        res[coupl_hbasistrip_dg1z][coupl_eftjhu_ghzgs2] = -(sqrt(sw)/sqrt(cw));
        res[coupl_hbasistrip_dg1z][coupl_eftjhu_ghgsgs2] = sw/(2.0*sw-1.0);
      }
      else{
        res.assign(nHiggsBasis_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
        res[coupl_hbasis_dCz][coupl_eftjhu_ghz1] = .5;
        res[coupl_hbasis_Czz][coupl_eftjhu_ghz2] = -((2*cw*sw)/pow(e,2));
        res[coupl_hbasis_Czbx][coupl_eftjhu_ghz1_prime2] = sw/pow(e,2);
        res[coupl_hbasis_tCzz][coupl_eftjhu_ghz4] = -((2*cw*sw)/pow(e,2));
        res[coupl_hbasis_dCw][coupl_eftjhu_ghz1] = .5;
        res[coupl_hbasis_Cww][coupl_eftjhu_ghz2] =  -(2*cw*sw)/pow(e,2);
        res[coupl_hbasis_Cww][coupl_eftjhu_ghzgs2] = -((4*sqrt(cw)*pow(sw,3.0/2))/pow(e,2));
        res[coupl_hbasis_Cww][coupl_eftjhu_ghgsgs2] =  -((2*pow(sw,2))/pow(e,2));
        res[coupl_hbasis_Cwbx][coupl_eftjhu_ghz1_prime2] = (pow(MW,2)*sw)/(pow(e,2)*pow(MZ,2)*(cw - sw));
        res[coupl_hbasis_Cwbx][coupl_eftjhu_ghz2] = -((2*pow(MW,2)*pow(sw,2))/(pow(e,2)*pow(MZ,2)*(cw - sw)));
        res[coupl_hbasis_Cwbx][coupl_eftjhu_ghzgs2] = ((2*pow(MW,2)*pow(sw,3.0/2.0))/(sqrt(cw)*pow(e,2)*pow(MZ,2)));
        res[coupl_hbasis_Cwbx][coupl_eftjhu_ghgsgs2] = (2*pow(MW,2)*pow(sw,2))/(pow(e,2)*pow(MZ,2)*(cw - sw));
        res[coupl_hbasis_tCww][coupl_eftjhu_ghz4] = -((2*cw*sw)/pow(e,2));
        res[coupl_hbasis_tCww][coupl_eftjhu_ghzgs4] = -((4*sqrt(cw)*pow(sw,3.0/2))/pow(e,2));
        res[coupl_hbasis_tCww][coupl_eftjhu_ghgsgs4] =  -((2*pow(sw,2))/pow(e,2));
        res[coupl_hbasis_Cza][coupl_eftjhu_ghzgs2] = -((2*sqrt(cw)*sqrt(sw))/pow(e,2));
        res[coupl_hbasis_tCza][coupl_eftjhu_ghzgs4] = -((2*sqrt(cw)*sqrt(sw))/pow(e,2));
        res[coupl_hbasis_Cabx][coupl_eftjhu_ghz1_prime2] = (2*cw*sw)/(pow(e,2)*(cw - sw));
        res[coupl_hbasis_Cabx][coupl_eftjhu_ghz2] = -((2*cw*sw)/(pow(e,2)*(cw - sw)));
        res[coupl_hbasis_Cabx][coupl_eftjhu_ghzgs2] = (2*sqrt(cw)*sqrt(sw))/pow(e,2);
        res[coupl_hbasis_Cabx][coupl_eftjhu_ghgsgs2] = (2*cw*sw)/(pow(e,2)*(cw - sw));
        res[coupl_hbasis_Caa][coupl_eftjhu_ghgsgs2] = -(2/pow(e,2));
        res[coupl_hbasis_tCaa][coupl_eftjhu_ghgsgs4] = -(2/pow(e,2));
        res[coupl_hbasis_Cgg][coupl_eftjhu_ghg2] = -2.0/1.0;
        res[coupl_hbasis_tCgg][coupl_eftjhu_ghg4] = -2.0/1.0;
      }
    }
    else if (basis_output == bEFT_HiggsBasis){
      if (include_triple_quartic_gauge){
        res.assign(nEFT_HiggsBasis_Include_Triple_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
        res[coupl_efthbasistrip_dCz][coupl_eftjhu_ghz1] = .5;
        res[coupl_efthbasistrip_Czz][coupl_eftjhu_ghz2] = -((2*cw*sw)/pow(e,2));
        res[coupl_efthbasistrip_Czbx][coupl_eftjhu_ghz1_prime2] = sw/pow(e,2);
        res[coupl_efthbasistrip_tCzz][coupl_eftjhu_ghz4] = -((2*cw*sw)/pow(e,2));
        res[coupl_efthbasistrip_Cza][coupl_eftjhu_ghzgs2] = -((2*sqrt(cw)*sqrt(sw))/pow(e,2));
        res[coupl_efthbasistrip_tCza][coupl_eftjhu_ghzgs4] = -((2*sqrt(cw)*sqrt(sw))/pow(e,2));
        res[coupl_efthbasistrip_Caa][coupl_eftjhu_ghgsgs2] = -(2/pow(e,2));
        res[coupl_efthbasistrip_tCaa][coupl_eftjhu_ghgsgs4] = -(2/pow(e,2));
        res[coupl_efthbasistrip_Cgg][coupl_eftjhu_ghg2] = -2.0/1.0;
        res[coupl_efthbasistrip_tCgg][coupl_eftjhu_ghg4] = -2.0/1.0;
        res[coupl_efthbasistrip_dKa][coupl_eftjhu_ghz2] = -cw/2.0;
        res[coupl_efthbasistrip_dKa][coupl_eftjhu_ghzgs2] = (1.0/2.0)*(sqrt(cw)*(1.0-2.0*sw))/(sqrt(sw));
        res[coupl_efthbasistrip_dKa][coupl_eftjhu_ghgsgs2] = cw/2;
        res[coupl_efthbasistrip_tKa][coupl_eftjhu_ghz4] = -cw/2;
        res[coupl_efthbasistrip_tKa][coupl_eftjhu_ghzgs4] = 1.0/2*(sqrt(cw)/sqrt(sw) - 2*sqrt(cw)*sqrt(sw));
        res[coupl_efthbasistrip_tKa][coupl_eftjhu_ghgsgs4] = cw/2;
        res[coupl_efthbasistrip_dKz][coupl_eftjhu_ghz1_prime2] = -(1.0/(2.0*(2.0-4.0*sw)));
        res[coupl_efthbasistrip_dKz][coupl_eftjhu_ghz2] = (sw*(sw-1.0))/(2.0*sw-1.0);
        res[coupl_efthbasistrip_dKz][coupl_eftjhu_ghzgs2] = -sqrt(cw)*sqrt(sw);
        res[coupl_efthbasistrip_dKz][coupl_eftjhu_ghgsgs2] = -(sw*(sw-1.0))/(2.0*sw-1.0);
        res[coupl_efthbasistrip_tKz][coupl_eftjhu_ghz4] = sw/2;
        res[coupl_efthbasistrip_tKz][coupl_eftjhu_ghzgs4] = 1.0/2 *(-(sqrt(sw)/sqrt(cw)) + (2*pow(sw,3.0/2))/sqrt(cw));
        res[coupl_efthbasistrip_tKz][coupl_eftjhu_ghgsgs4] = -sw/2;
        res[coupl_efthbasistrip_dg1z][coupl_eftjhu_ghz1_prime2] = (1.0/(2.0*(2.0*sw-1.0)));
        res[coupl_efthbasistrip_dg1z][coupl_eftjhu_ghz2] = -sw/(2.0*sw-1.0);
        res[coupl_efthbasistrip_dg1z][coupl_eftjhu_ghzgs2] = -(sqrt(sw)/sqrt(cw));
        res[coupl_efthbasistrip_dg1z][coupl_eftjhu_ghgsgs2] = sw/(2.0*sw-1.0);
      }
      else{
        res.assign(nEFT_HiggsBasis_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
        res[coupl_efthbasis_dCz][coupl_eftjhu_ghz1] = .5;
        res[coupl_efthbasis_Czbx][coupl_eftjhu_ghz1_prime2] = sw/pow(e,2);
        res[coupl_efthbasis_Czz][coupl_eftjhu_ghz2] = -((2*cw*sw)/pow(e,2));
        res[coupl_efthbasis_tCzz][coupl_eftjhu_ghz4] = -((2*cw*sw)/pow(e,2));
        res[coupl_efthbasis_Cza][coupl_eftjhu_ghzgs2] = -((2*sqrt(cw)*sqrt(sw))/pow(e,2));
        res[coupl_efthbasis_tCza][coupl_eftjhu_ghzgs4] = -((2*sqrt(cw)*sqrt(sw))/pow(e,2));
        res[coupl_efthbasis_Caa][coupl_eftjhu_ghgsgs2] = -(2/pow(e,2));
        res[coupl_efthbasis_tCaa][coupl_eftjhu_ghgsgs4] = -(2/pow(e,2));
        res[coupl_efthbasis_Cgg][coupl_eftjhu_ghg2] = -2.0/1.0;
        res[coupl_efthbasis_tCgg][coupl_eftjhu_ghg4] = -2/1.0;
      }
    }
    else if (basis_output == bWarsawBasis){
      if(custodial_symmetry){
        res.assign(nWarsawBasis_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
        res[coupl_warsaw_cHbx][coupl_eftjhu_ghz1]=(1.0/vev_lam)*1.0/2.0;
        res[coupl_warsaw_cHbx][coupl_eftjhu_ghz2]=(1.0/vev_lam)*-(2.0*cw*sw)/(cw-sw);
        res[coupl_warsaw_cHbx][coupl_eftjhu_ghz1_prime2]=(1.0/vev_lam)*(3.0*(cw-sw))/(2.0*(cw-sw));
        res[coupl_warsaw_cHbx][coupl_eftjhu_ghzgs2]=(1.0/vev_lam)*2.0*sqrt(cw)*sqrt(sw);
        res[coupl_warsaw_cHbx][coupl_eftjhu_ghgsgs2]=(1.0/vev_lam)*(2.0*cw*sw)/(cw-sw);
        res[coupl_warsaw_cHG][coupl_eftjhu_ghg2]=(1.0/vev_lam)*-1.0/(2.0);
        res[coupl_warsaw_cHW][coupl_eftjhu_ghz2]=(1.0/vev_lam)*-cw/(2.0);
        res[coupl_warsaw_cHW][coupl_eftjhu_ghzgs2]=(1.0/vev_lam)*-sqrt(cw)*sqrt(sw);
        res[coupl_warsaw_cHW][coupl_eftjhu_ghgsgs2]=(1.0/vev_lam)*-(sw)/2.0;
        res[coupl_warsaw_cHWB][coupl_eftjhu_ghz2]=(1.0/vev_lam)*-sqrt(cw)*sqrt(sw) ;
        res[coupl_warsaw_cHWB][coupl_eftjhu_ghzgs2]=(1.0/vev_lam)*(cw-sw);
        res[coupl_warsaw_cHWB][coupl_eftjhu_ghgsgs2]=(1.0/vev_lam)*sqrt(cw)*sqrt(sw) ;
        res[coupl_warsaw_cHB][coupl_eftjhu_ghz2]=(1.0/vev_lam)*-sw/(2.0);
        res[coupl_warsaw_cHB][coupl_eftjhu_ghzgs2]=(1.0/vev_lam)*sqrt(cw)*sqrt(sw);
        res[coupl_warsaw_cHB][coupl_eftjhu_ghgsgs2]=(1.0/vev_lam)*-(cw)/2.0;
        res[coupl_warsaw_cHD][coupl_eftjhu_ghz2]=(1.0/vev_lam)*(4.0*cw*sw)/(cw-sw);
        res[coupl_warsaw_cHD][coupl_eftjhu_ghz1_prime2]=(1.0/vev_lam)*(-2.0*sw)/((cw-sw));
        res[coupl_warsaw_cHD][coupl_eftjhu_ghzgs2]=(1.0/vev_lam)*-4.0*sqrt(sw)*sqrt(cw);
        res[coupl_warsaw_cHD][coupl_eftjhu_ghgsgs2]=(1.0/vev_lam)*(-4.0*cw*sw)/(cw-sw);
        res[coupl_warsaw_tcHB][coupl_eftjhu_ghz4]=(1.0/vev_lam)*-sw/(2.0);
        res[coupl_warsaw_tcHB][coupl_eftjhu_ghzgs4]=(1.0/vev_lam)*sqrt(cw)*sqrt(sw);
        res[coupl_warsaw_tcHB][coupl_eftjhu_ghgsgs4]=(1.0/vev_lam)*-cw/2.0;
        res[coupl_warsaw_tcHW][coupl_eftjhu_ghz4]=(1.0/vev_lam)*-cw/2.0;
        res[coupl_warsaw_tcHW][coupl_eftjhu_ghzgs4]=(1.0/vev_lam)*-sqrt(cw)*sqrt(sw);
        res[coupl_warsaw_tcHW][coupl_eftjhu_ghgsgs4]=(1.0/vev_lam)*-sw/2.0;
        res[coupl_warsaw_tcHWB][coupl_eftjhu_ghz4]=(1.0/vev_lam)*-sqrt(cw)*sqrt(sw);
        res[coupl_warsaw_tcHWB][coupl_eftjhu_ghzgs4]=(1.0/vev_lam)*(cw-sw);
        res[coupl_warsaw_tcHWB][coupl_eftjhu_ghgsgs4]=(1.0/vev_lam)*sqrt(cw)*sqrt(sw);
        res[coupl_warsaw_tcHG][coupl_eftjhu_ghg4]=(1.0/vev_lam)*-1.0/2.0;
      }
      else{
	res.assign(nWarsawBasis_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
        res[coupl_warsaw_cHbx][coupl_eftjhu_ghz1]=(1.0/vev_lam)*1.0/2.0;
	res[coupl_warsaw_cHbx][coupl_eftjhu_ghz1_prime2]=(1.0/vev_lam)*(-3.0*cw+sw)/(2.0*(sw-cw));
        res[coupl_warsaw_cHbx][coupl_eftjhu_ghz2]=(1.0/vev_lam)*(-2.0*cw*sw)/(cw-sw);
        res[coupl_warsaw_cHbx][coupl_eftjhu_ghzgs2]=(1.0/vev_lam)*(2.0*sqrt(sw)*sqrt(cw));
        res[coupl_warsaw_cHbx][coupl_eftjhu_ghgsgs2]=(1.0/vev_lam)*(2.0*cw*sw)/(cw-sw);
        res[coupl_warsaw_cHG][coupl_eftjhu_ghg2]=(1.0/vev_lam)*-1.0/2.0;
        res[coupl_warsaw_cHW][coupl_eftjhu_ghz2]=(1.0/vev_lam)*-cw/(2.0);
        res[coupl_warsaw_cHW][coupl_eftjhu_ghzgs2]=(1.0/vev_lam)*-sqrt(cw)*sqrt(sw);
        res[coupl_warsaw_cHW][coupl_eftjhu_ghgsgs2]=(1.0/vev_lam)*-(sw)/2.0;
        res[coupl_warsaw_cHWB][coupl_eftjhu_ghz2]=(1.0/vev_lam)*-sqrt(cw)*sqrt(sw);
        res[coupl_warsaw_cHWB][coupl_eftjhu_ghzgs2]=(1.0/vev_lam)*(cw-sw);
        res[coupl_warsaw_cHWB][coupl_eftjhu_ghgsgs2]=(1.0/vev_lam)*sqrt(cw)*sqrt(sw);
        res[coupl_warsaw_cHB][coupl_eftjhu_ghz2]=(1.0/vev_lam)*-sw/(2.0);
        res[coupl_warsaw_cHB][coupl_eftjhu_ghzgs2]=(1.0/vev_lam)*sqrt(cw)*sqrt(sw);
        res[coupl_warsaw_cHB][coupl_eftjhu_ghgsgs2]=(1.0/vev_lam)*-(cw)/2.0;
        res[coupl_warsaw_cHD][coupl_eftjhu_ghz1_prime2]=(1.0/vev_lam)*-(2.0*sw)/(cw-sw);
        res[coupl_warsaw_cHD][coupl_eftjhu_ghz2]=(1.0/vev_lam)*(4.0*sw*cw)/(cw-sw);
        res[coupl_warsaw_cHD][coupl_eftjhu_ghzgs2]=(1.0/vev_lam)*(-4.0*sqrt(sw)*sqrt(cw));
        res[coupl_warsaw_cHD][coupl_eftjhu_ghgsgs2]=(1.0/vev_lam)*(-4.0*cw*sw)/(cw-sw);
        res[coupl_warsaw_tcHB][coupl_eftjhu_ghz4]=(1.0/vev_lam)*-sw/(2.0);
        res[coupl_warsaw_tcHB][coupl_eftjhu_ghzgs4]=(1.0/vev_lam)*sqrt(cw)*sqrt(sw);
        res[coupl_warsaw_tcHB][coupl_eftjhu_ghgsgs4]=(1.0/vev_lam)*-cw/2.0;
        res[coupl_warsaw_tcHW][coupl_eftjhu_ghz4]=(1.0/vev_lam)*-cw/2.0;
        res[coupl_warsaw_tcHW][coupl_eftjhu_ghzgs4]=(1.0/vev_lam)*-sqrt(cw)*sqrt(sw);
        res[coupl_warsaw_tcHW][coupl_eftjhu_ghgsgs4]=(1.0/vev_lam)*-sw/2.0;
        res[coupl_warsaw_tcHWB][coupl_eftjhu_ghz4]=(1.0/vev_lam)*-sqrt(cw)*sqrt(sw);
        res[coupl_warsaw_tcHWB][coupl_eftjhu_ghzgs4]=(1.0/vev_lam)*(cw-sw);
        res[coupl_warsaw_tcHWB][coupl_eftjhu_ghgsgs4]=(1.0/vev_lam)*sqrt(cw)*sqrt(sw);
        res[coupl_warsaw_tcHG][coupl_eftjhu_ghg4]=(1.0/vev_lam)*-1.0/2.0;
      }
    }
  }
  else if (basis_input == bHiggsBasis){
    if (basis_output == bAmplitude_JHUGen){
      res.assign(nAmplitude_JHUGen_CouplingTypes, std::vector<double>(nHiggsBasis_CouplingTypes, 0));
      res[coupl_ampjhu_ghz1][coupl_hbasis_dCz] = 2.0;
      res[coupl_ampjhu_ghz1_prime2][coupl_hbasis_Czbx] = (pow(e,2))/(sw);
      res[coupl_ampjhu_ghz2][coupl_hbasis_Czz] = -(pow(e,2)/(2.0*cw*sw));
      res[coupl_ampjhu_ghz4][coupl_hbasis_tCzz] = -(pow(e,2)/(2.0*cw*sw));
      res[coupl_ampjhu_ghw1][coupl_hbasis_dCw] = 2.0;
      res[coupl_ampjhu_ghw1_prime2][coupl_hbasis_Cwbx] = (pow(e,2))/(sw);
      res[coupl_ampjhu_ghw2][coupl_hbasis_Cww] = -(pow(e,2)/(2.0*sw));
      res[coupl_ampjhu_ghw4][coupl_hbasis_tCww] = -(pow(e,2)/(2.0*sw));
      res[coupl_ampjhu_ghzgs1_prime2][coupl_hbasis_Cabx] = (pow(e,2))/(sqrt(cw)*sqrt(sw));
      res[coupl_ampjhu_ghzgs2][coupl_hbasis_Cza] = -(pow(e,2)/(2.0*sqrt(cw)*sqrt(sw)));
      res[coupl_ampjhu_ghzgs4][coupl_hbasis_tCza] = -(pow(e,2)/(2.0*sqrt(cw)*sqrt(sw)));
      res[coupl_ampjhu_ghgsgs2][coupl_hbasis_Caa] =  -(pow(e,2)/2.0);
      res[coupl_ampjhu_ghgsgs4][coupl_hbasis_tCaa] =  -(pow(e,2)/2.0);
      res[coupl_ampjhu_ghg2][coupl_hbasis_Cgg] = -(1.0/2.0);
      res[coupl_ampjhu_ghg4][coupl_hbasis_tCgg] = -(1.0/2.0);
    }
    else if (basis_output == bEFT_JHUGen){
      res.assign(nEFT_JHUGen_CouplingTypes, std::vector<double>(nHiggsBasis_CouplingTypes, 0));
    }
    else if (basis_output == bEFT_HiggsBasis){
      res.assign(nEFT_HiggsBasis_CouplingTypes, std::vector<double>(nHiggsBasis_CouplingTypes, 0));
    }
    else if (basis_output == bHiggsBasis){
      res.assign(nHiggsBasis_CouplingTypes, std::vector<double>(nHiggsBasis_CouplingTypes, 0));
      for (size_t i=0; i<(size_t) nHiggsBasis_CouplingTypes; i++) res.at(i).at(i)=1.0;
    }
  }
  else if (basis_input == bEFT_HiggsBasis){
    if (basis_output == bAmplitude_JHUGen){
      if (include_triple_quartic_gauge){
        res.assign(nAmplitude_JHUGen_Include_Triple_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
        res[coupl_ampjhutrip_ghz1][coupl_efthbasis_dCz] = 2;
        res[coupl_ampjhutrip_ghz1_prime2][coupl_efthbasis_Czbx] = (pow(e,2))/(sw);
        res[coupl_ampjhutrip_ghz2][coupl_efthbasis_Czz] = -(pow(e,2)/(2*cw*sw));
        res[coupl_ampjhutrip_ghz4][coupl_efthbasis_tCzz] = -(pow(e,2))/(2*cw*sw);
        res[coupl_ampjhutrip_ghw1][coupl_efthbasis_dCz] = 2;
        res[coupl_ampjhutrip_ghw1_prime2][coupl_efthbasis_Czbx] = (pow(e,2)*pow(MW,2))/(pow(MZ,2)*sw*(cw - sw));
        res[coupl_ampjhutrip_ghw1_prime2][coupl_efthbasis_Czz] = (pow(e,2)*pow(MW,2))/(cw*pow(MZ,2)*(cw - sw));
        res[coupl_ampjhutrip_ghw1_prime2][coupl_efthbasis_Cza] = -((pow(e,2)*pow(MW,2))/(cw*pow(MZ,2)));
        res[coupl_ampjhutrip_ghw1_prime2][coupl_efthbasis_Caa] = -((pow(e,2)*pow(MW,2)*sw)/(pow(MZ,2)*(cw - sw)));
        res[coupl_ampjhutrip_ghw2][coupl_efthbasis_Czz] = -(pow(e,2)/(2*sw));
        res[coupl_ampjhutrip_ghw2][coupl_efthbasis_Cza] = -pow(e,2);
        res[coupl_ampjhutrip_ghw2][coupl_efthbasis_Caa] = -(1.0/2)*pow(e,2)*sw;
        res[coupl_ampjhutrip_ghw4][coupl_efthbasis_tCzz] = -(pow(e,2))/(2*sw);
        res[coupl_ampjhutrip_ghw4][coupl_efthbasis_tCza] = -pow(e,2);
        res[coupl_ampjhutrip_ghw4][coupl_efthbasis_tCaa] = -(1.0/2)*pow(e,2)*sw;
        res[coupl_ampjhutrip_ghzgs1_prime2][coupl_efthbasis_Czbx] = (2*sqrt(cw)*pow(e,2))/(sqrt(sw)*(cw - sw));
        res[coupl_ampjhutrip_ghzgs1_prime2][coupl_efthbasis_Czz] = (pow(e,2))/(sqrt(cw)*sqrt(sw)*(cw - sw));
        res[coupl_ampjhutrip_ghzgs1_prime2][coupl_efthbasis_Cza] = -((pow(e,2))/(sqrt(cw)*sqrt(sw)));
        res[coupl_ampjhutrip_ghzgs1_prime2][coupl_efthbasis_Caa] = -((sqrt(cw)*pow(e,2)*sqrt(sw))/(cw - sw));
        res[coupl_ampjhutrip_ghzgs2][coupl_efthbasis_Cza] = -(pow(e,2)/(2*sqrt(cw)*sqrt(sw)));
        res[coupl_ampjhutrip_ghzgs4][coupl_efthbasis_tCza] = -(pow(e,2)/(2*sqrt(cw)*sqrt(sw)));
        res[coupl_ampjhutrip_ghgsgs2][coupl_efthbasis_Caa] = -(pow(e,2)/2);
        res[coupl_ampjhutrip_ghgsgs4][coupl_efthbasis_tCaa] = -(pow(e,2)/2);
        res[coupl_ampjhutrip_ghg2][coupl_efthbasis_Cgg] = -(1.0/2.0);
        res[coupl_ampjhutrip_ghg4][coupl_efthbasis_tCgg] = -(1.0/2.0);
        res[coupl_ampjhutrip_dV_Z][coupl_efthbasis_Czbx] = -(pow(e,2))/(2*sw*(cw - sw));
        res[coupl_ampjhutrip_dV_Z][coupl_efthbasis_Czz] = (pow(e,2)/(-cw + sw));
        res[coupl_ampjhutrip_dV_Z][coupl_efthbasis_Cza] = pow(e,2);
        res[coupl_ampjhutrip_dV_Z][coupl_efthbasis_Caa] = -(cw*pow(e,2)*sw)/(-cw + sw);
        res[coupl_ampjhutrip_dV_A][coupl_efthbasis_Czz] = (pow(e,2)/(2*sw));
        res[coupl_ampjhutrip_dV_A][coupl_efthbasis_Cza] = ((sw-cw)*pow(e,2))/(2.0*sw);
        res[coupl_ampjhutrip_dV_A][coupl_efthbasis_Caa] = -(1.0/2.0)*cw*pow(e,2);
        res[coupl_ampjhutrip_dP_Z][coupl_efthbasis_Czbx] = -((pow(e,2))/(2*sw*(cw - sw)));
        res[coupl_ampjhutrip_dP_Z][coupl_efthbasis_Czz] = -(pow(e,2)/(2*cw*(cw - sw)));
        res[coupl_ampjhutrip_dP_Z][coupl_efthbasis_Cza] = pow(e,2)/(2.0*cw);
        res[coupl_ampjhutrip_dP_Z][coupl_efthbasis_Caa] = ((pow(e,2)*sw)/(2*(cw - sw)));
        //res[coupl_ampjhutrip_dP_A][coupl_efthbasis_dCz] = 1.0;
        res[coupl_ampjhutrip_dM_Z][coupl_efthbasis_Czbx] = -((pow(e,2))/(2*sw*(cw - sw)));
        res[coupl_ampjhutrip_dM_Z][coupl_efthbasis_Czz] = -(pow(e,2)/(2*cw*(cw - sw)));
        res[coupl_ampjhutrip_dM_Z][coupl_efthbasis_Cza] = pow(e,2)/(2.0*cw);
        res[coupl_ampjhutrip_dM_Z][coupl_efthbasis_Caa] = ((pow(e,2)*sw)/(2*(cw - sw)));
        //res[coupl_ampjhutrip_dM_A][coupl_efthbasis_dCz] = 1.0;
        res[coupl_ampjhutrip_dFour_Z][coupl_efthbasis_tCzz] = -(pow(e,2)/(2*cw));
        res[coupl_ampjhutrip_dFour_Z][coupl_efthbasis_tCza] = -((pow(e,2)*(-(sqrt(sw)/sqrt(cw)) + (2*pow(sw,3.0/2))/sqrt(cw)))/(2*sqrt(cw)*sqrt(sw)));
        res[coupl_ampjhutrip_dFour_Z][coupl_efthbasis_tCaa] = (pow(e,2)*sw)/2;
        res[coupl_ampjhutrip_dFour_A][coupl_efthbasis_tCzz] = pow(e,2)/(2*sw);
        res[coupl_ampjhutrip_dFour_A][coupl_efthbasis_tCza] = -((pow(e,2)*(sqrt(cw)/sqrt(sw) - 2*sqrt(cw)*sqrt(sw)))/(2*sqrt(cw)*sqrt(sw)));
        res[coupl_ampjhutrip_dFour_A][coupl_efthbasis_tCaa] = -(1.0/2)*cw*pow(e,2);
        res[coupl_ampjhutrip_dZZWpWm][coupl_efthbasis_Czbx] = -((pow(e,2)*cw)/(pow(sw,2)*(cw - sw)));
        res[coupl_ampjhutrip_dZZWpWm][coupl_efthbasis_Czz] = -(pow(e,2)/(sw*(cw - sw)));
        res[coupl_ampjhutrip_dZZWpWm][coupl_efthbasis_Cza] = (pow(e,2))/(sw);
        res[coupl_ampjhutrip_dZZWpWm][coupl_efthbasis_Caa] = ((pow(e,2)*cw)/(cw - sw));
        res[coupl_ampjhutrip_dZAWpWm][coupl_efthbasis_Czbx] = -((pow(e,2)*sqrt(cw))/(2*pow(sw,3.0/2.0)*(cw - sw)));
        res[coupl_ampjhutrip_dZAWpWm][coupl_efthbasis_Czz] = -(pow(e,2)/(2*sqrt(cw)*sqrt(sw)*(cw - sw)));
        res[coupl_ampjhutrip_dZAWpWm][coupl_efthbasis_Cza] = pow(e,2)/(2*sqrt(sw)*sqrt(cw));
        res[coupl_ampjhutrip_dZAWpWm][coupl_efthbasis_Caa] = ((pow(e,2)*sqrt(sw)*sqrt(cw))/(2*(cw - sw)));
        // res[coupl_ampjhutrip_dZZWpWm][coupl_efthbasis_dCz] = 1.0;
      }//
      else{
        res.assign(nAmplitude_JHUGen_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
        res[coupl_ampjhu_ghz1][coupl_efthbasis_dCz] = 2;
        res[coupl_ampjhu_ghz1_prime2][coupl_efthbasis_Czbx] = (pow(e,2))/(sw);
        res[coupl_ampjhu_ghz2][coupl_efthbasis_Czz] = -(pow(e,2)/(2*cw*sw));
        res[coupl_ampjhu_ghz4][coupl_efthbasis_tCzz] = -(pow(e,2))/(2*cw*sw);
        res[coupl_ampjhu_ghw1][coupl_efthbasis_dCz] = 2;
        res[coupl_ampjhu_ghw1_prime2][coupl_efthbasis_Czbx] = (pow(e,2)*pow(MW,2))/(pow(MZ,2)*sw*(cw - sw));
        res[coupl_ampjhu_ghw1_prime2][coupl_efthbasis_Czz] = (pow(e,2)*pow(MW,2))/(cw*pow(MZ,2)*(cw - sw));
        res[coupl_ampjhu_ghw1_prime2][coupl_efthbasis_Cza] = -((pow(e,2)*pow(MW,2))/(cw*pow(MZ,2)));
        res[coupl_ampjhu_ghw1_prime2][coupl_efthbasis_Caa] = -((pow(e,2)*pow(MW,2)*sw)/(pow(MZ,2)*(cw - sw)));
        res[coupl_ampjhu_ghw2][coupl_efthbasis_Czz] = -(pow(e,2)/(2*sw));
        res[coupl_ampjhu_ghw2][coupl_efthbasis_Cza] = -pow(e,2);
        res[coupl_ampjhu_ghw2][coupl_efthbasis_Caa] = -(1.0/2)*pow(e,2)*sw;
        res[coupl_ampjhu_ghw4][coupl_efthbasis_tCzz] = -(pow(e,2))/(2*sw);
        res[coupl_ampjhu_ghw4][coupl_efthbasis_tCza] = -pow(e,2);
        res[coupl_ampjhu_ghw4][coupl_efthbasis_tCaa] = -(1.0/2)*pow(e,2)*sw;
        res[coupl_ampjhu_ghzgs1_prime2][coupl_efthbasis_Czbx] = (2*sqrt(cw)*pow(e,2))/(sqrt(sw)*(cw - sw));
        res[coupl_ampjhu_ghzgs1_prime2][coupl_efthbasis_Czz] = (pow(e,2))/(sqrt(cw)*sqrt(sw)*(cw - sw));
        res[coupl_ampjhu_ghzgs1_prime2][coupl_efthbasis_Cza] = -((pow(e,2))/(sqrt(cw)*sqrt(sw)));
        res[coupl_ampjhu_ghzgs1_prime2][coupl_efthbasis_Caa] = -((sqrt(cw)*pow(e,2)*sqrt(sw))/(cw - sw));
        res[coupl_ampjhu_ghzgs2][coupl_efthbasis_Cza] = -(pow(e,2)/(2*sqrt(cw)*sqrt(sw)));
        res[coupl_ampjhu_ghzgs4][coupl_efthbasis_tCza] = -(pow(e,2)/(2*sqrt(cw)*sqrt(sw)));
        res[coupl_ampjhu_ghgsgs2][coupl_efthbasis_Caa] = -(pow(e,2)/2);
        res[coupl_ampjhu_ghgsgs4][coupl_efthbasis_tCaa] = -(pow(e,2)/2);
        res[coupl_ampjhu_ghg2][coupl_efthbasis_Cgg] = -(1.0/2);
        res[coupl_ampjhu_ghg4][coupl_efthbasis_tCgg] = -(1.0/2);
      }//
    }//
    else if (basis_output == bEFT_JHUGen){
      if(include_triple_quartic_gauge){
        res.assign(nEFT_JHUGen_Include_Triple_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
        res[coupl_eftjhutrip_ghz1][coupl_efthbasis_dCz]= 2;
        res[coupl_eftjhutrip_ghz1_prime2][coupl_efthbasis_Czbx]= (pow(e,2))/(sw);
        res[coupl_eftjhutrip_ghz2][coupl_efthbasis_Czz]= -(pow(e,2))/(2*cw*sw);
        res[coupl_eftjhutrip_ghz4][coupl_efthbasis_tCzz]= -pow(e,2)/(2*cw*sw);
        res[coupl_eftjhutrip_ghzgs2][coupl_efthbasis_Cza]= -(pow(e,2)/(2*sqrt(cw)*sqrt(sw)));
        res[coupl_eftjhutrip_ghzgs4][coupl_efthbasis_tCza]= -(pow(e,2)/(2*sqrt(cw)*sqrt(sw)));
        res[coupl_eftjhutrip_ghgsgs2][coupl_efthbasis_Caa]=  -(pow(e,2)/2);
        res[coupl_eftjhutrip_ghgsgs4][coupl_efthbasis_tCaa]= -(pow(e,2)/2);
        res[coupl_eftjhutrip_ghg2][coupl_efthbasis_Cgg]= -(1.0/2);
        res[coupl_eftjhutrip_ghg4][coupl_efthbasis_tCgg]= -(1.0/2);
        res[coupl_eftjhutrip_dV_Z][coupl_efthbasis_Czbx] = -(pow(e,2))/(2*sw*(cw - sw));
        res[coupl_eftjhutrip_dV_Z][coupl_efthbasis_Czz] = (pow(e,2)/(-cw + sw));
        res[coupl_eftjhutrip_dV_Z][coupl_efthbasis_Cza] = pow(e,2);
        res[coupl_eftjhutrip_dV_Z][coupl_efthbasis_Caa] = -(cw*pow(e,2)*sw)/(-cw + sw);
        res[coupl_eftjhutrip_dV_A][coupl_efthbasis_Czz] = (pow(e,2)/(2*sw));
        res[coupl_eftjhutrip_dV_A][coupl_efthbasis_Cza] = ((sw-cw)*pow(e,2))/(2.0*sw);
        res[coupl_eftjhutrip_dV_A][coupl_efthbasis_Caa] = -(1.0/2.0)*cw*pow(e,2);
        res[coupl_eftjhutrip_dP_Z][coupl_efthbasis_Czbx] = -((pow(e,2))/(2*sw*(cw - sw)));
        res[coupl_eftjhutrip_dP_Z][coupl_efthbasis_Czz] = -(pow(e,2)/(2*cw*(cw - sw)));
        res[coupl_eftjhutrip_dP_Z][coupl_efthbasis_Cza] = pow(e,2)/(2.0*cw);
        res[coupl_eftjhutrip_dP_Z][coupl_efthbasis_Caa] = ((pow(e,2)*sw)/(2*(cw - sw)));
        //res[coupl_ampjhutrip_dP_A][coupl_efthbasis_dCz] = 1.0;
        res[coupl_eftjhutrip_dM_Z][coupl_efthbasis_Czbx] = -((pow(e,2))/(2*sw*(cw - sw)));
        res[coupl_eftjhutrip_dM_Z][coupl_efthbasis_Czz] = -(pow(e,2)/(2*cw*(cw - sw)));
        res[coupl_eftjhutrip_dM_Z][coupl_efthbasis_Cza] = pow(e,2)/(2.0*cw);
        res[coupl_eftjhutrip_dM_Z][coupl_efthbasis_Caa] = ((pow(e,2)*sw)/(2*(cw - sw)));
        //res[coupl_ampjhutrip_dM_A][coupl_efthbasis_dCz] = 1.0;
        res[coupl_eftjhutrip_dFour_Z][coupl_efthbasis_tCzz] = -(pow(e,2)/(2*cw));
        res[coupl_eftjhutrip_dFour_Z][coupl_efthbasis_tCza] = -((pow(e,2)*(-(sqrt(sw)/sqrt(cw)) + (2*pow(sw,3.0/2))/sqrt(cw)))/(2*sqrt(cw)*sqrt(sw)));
        res[coupl_eftjhutrip_dFour_Z][coupl_efthbasis_tCaa] = (pow(e,2)*sw)/2;
        res[coupl_eftjhutrip_dFour_A][coupl_efthbasis_tCzz] = pow(e,2)/(2*sw);
        res[coupl_eftjhutrip_dFour_A][coupl_efthbasis_tCza] = -((pow(e,2)*(sqrt(cw)/sqrt(sw) - 2*sqrt(cw)*sqrt(sw)))/(2*sqrt(cw)*sqrt(sw)));
        res[coupl_eftjhutrip_dFour_A][coupl_efthbasis_tCaa] = -(1.0/2)*cw*pow(e,2);
        res[coupl_eftjhutrip_dZZWpWm][coupl_efthbasis_Czbx] = -((pow(e,2)*cw)/(pow(sw,2)*(cw - sw)));
        res[coupl_eftjhutrip_dZZWpWm][coupl_efthbasis_Czz] = -(pow(e,2)/(sw*(cw - sw)));
        res[coupl_eftjhutrip_dZZWpWm][coupl_efthbasis_Cza] = (pow(e,2))/(sw);
        res[coupl_eftjhutrip_dZZWpWm][coupl_efthbasis_Caa] = ((pow(e,2)*cw)/(cw - sw));
        res[coupl_eftjhutrip_dZAWpWm][coupl_efthbasis_Czbx] = -((pow(e,2)*sqrt(cw))/(2*pow(sw,3.0/2.0)*(cw - sw)));
        res[coupl_eftjhutrip_dZAWpWm][coupl_efthbasis_Czz] = -(pow(e,2)/(2*sqrt(cw)*sqrt(sw)*(cw - sw)));
        res[coupl_eftjhutrip_dZAWpWm][coupl_efthbasis_Cza] = pow(e,2)/(2*sqrt(sw)*sqrt(cw));
        res[coupl_eftjhutrip_dZAWpWm][coupl_efthbasis_Caa] = ((pow(e,2)*sqrt(sw)*sqrt(cw))/(2*(cw - sw)));
        //res[coupl_ampjhutrip_dZZWpWm][coupl_efthbasis_dCz] = 1.0;
      }
      else{
        res.assign(nEFT_JHUGen_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
        res[coupl_eftjhu_ghz1][coupl_efthbasis_dCz]= 2;
        res[coupl_eftjhu_ghz1_prime2][coupl_efthbasis_Czbx]= (pow(e,2))/(sw);
        res[coupl_eftjhu_ghz2][coupl_efthbasis_Czz]= -(pow(e,2))/(2*cw*sw);
        res[coupl_eftjhu_ghz4][coupl_efthbasis_tCzz]= -pow(e,2)/(2*cw*sw);
        res[coupl_eftjhu_ghzgs2][coupl_efthbasis_Cza]= -(pow(e,2)/(2*sqrt(cw)*sqrt(sw)));
        res[coupl_eftjhu_ghzgs4][coupl_efthbasis_tCza]= -(pow(e,2)/(2*sqrt(cw)*sqrt(sw)));
        res[coupl_eftjhu_ghgsgs2][coupl_efthbasis_Caa]=  -(pow(e,2)/2);
        res[coupl_eftjhu_ghgsgs4][coupl_efthbasis_tCaa]= -(pow(e,2)/2);
        res[coupl_eftjhu_ghg2][coupl_efthbasis_Cgg]= -1.0/2;
        res[coupl_eftjhu_ghg4][coupl_efthbasis_tCgg]= -1.0/2;
      }
    }
    else if (basis_output == bHiggsBasis){
      if (include_triple_quartic_gauge){
        res.assign(nHiggsBasis_Include_Triple_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
        res[coupl_hbasistrip_dCz][coupl_efthbasis_dCz]= 1.0;
        res[coupl_hbasistrip_Czz][coupl_efthbasis_Czz]= 1.0;
        res[coupl_hbasistrip_Czbx][coupl_efthbasis_Czbx]= 1.0;
        res[coupl_hbasistrip_tCzz][coupl_efthbasis_tCzz]= 1.0;
        res[coupl_hbasistrip_dCw][coupl_efthbasis_dCz]= 1.0;
        res[coupl_hbasistrip_Cww][coupl_efthbasis_Czz]= 1.0;
        res[coupl_hbasistrip_Cww][coupl_efthbasis_Cza]= 2*sw;
        res[coupl_hbasistrip_Cww][coupl_efthbasis_Caa]= pow(sw,2);
        res[coupl_hbasistrip_Cwbx][coupl_efthbasis_Czbx]= pow(MW,2)/(pow(MZ,2)*(cw - sw));
        res[coupl_hbasistrip_Cwbx][coupl_efthbasis_Czz]= (pow(MW,2)*sw)/(cw*pow(MZ,2)*(cw - sw));
        res[coupl_hbasistrip_Cwbx][coupl_efthbasis_Cza]= -((pow(MW,2)*sw)/(cw*pow(MZ,2)));
        res[coupl_hbasistrip_Cwbx][coupl_efthbasis_Caa]= -((pow(MW,2)*pow(sw,2))/(pow(MZ,2)*(cw - sw)));
        res[coupl_hbasistrip_tCww][coupl_efthbasis_tCzz]= 1.0;
        res[coupl_hbasistrip_tCww][coupl_efthbasis_tCza]= 2*sw;
        res[coupl_hbasistrip_tCww][coupl_efthbasis_tCaa]= pow(sw,2);
        res[coupl_hbasistrip_Cza][coupl_efthbasis_Cza]= 1.0;
        res[coupl_hbasistrip_tCza][coupl_efthbasis_tCza]= 1.0;
        res[coupl_hbasistrip_Cabx][coupl_efthbasis_Czbx]= (2*cw)/(cw - sw);
        res[coupl_hbasistrip_Cabx][coupl_efthbasis_Czz]= 1.0/(cw - sw);
        res[coupl_hbasistrip_Cabx][coupl_efthbasis_Cza]= -1.0;
        res[coupl_hbasistrip_Cabx][coupl_efthbasis_Caa]= -((cw*sw)/(cw - sw));
        res[coupl_hbasistrip_Caa][coupl_efthbasis_Caa]= 1.0;
        res[coupl_hbasistrip_tCaa][coupl_efthbasis_tCaa]= 1.0;
        res[coupl_hbasistrip_Cgg][coupl_efthbasis_Cgg]= 1.0;
        res[coupl_hbasistrip_tCgg][coupl_efthbasis_tCgg]= 1.0;
        res[coupl_hbasistrip_dKa][coupl_efthbasis_Czz]= (pow(e,2)/(4*sw));
        res[coupl_hbasistrip_dKa][coupl_efthbasis_Cza]= (pow(e,2)*(2.0*sw-1.0))/(4.0*sw);
        res[coupl_hbasistrip_dKa][coupl_efthbasis_Caa]= -(1.0/4)*cw*pow(e,2);
        res[coupl_hbasistrip_tKa][coupl_efthbasis_tCzz]= pow(e,2)/(4*sw);
        res[coupl_hbasistrip_tKa][coupl_efthbasis_tCza]= -((pow(e,2)*(sqrt(cw)/sqrt(sw) - 2*sqrt(cw)*sqrt(sw)))/(4*sqrt(cw)*sqrt(sw)));
        res[coupl_hbasistrip_tKa][coupl_efthbasis_tCaa]= -(1.0/4)*cw*pow(e,2);
        res[coupl_hbasistrip_dKz][coupl_efthbasis_Czbx]= -(pow(e,2)/(4*sw*(cw - sw)));
        res[coupl_hbasistrip_dKz][coupl_efthbasis_Czz]= (pow(e,2)/(2*(-cw + sw)));
        res[coupl_hbasistrip_dKz][coupl_efthbasis_Cza]= pow(e,2)/2;
        res[coupl_hbasistrip_dKz][coupl_efthbasis_Caa]= -(cw*pow(e,2)*sw)/(2*(-cw + sw));
        res[coupl_hbasistrip_tKz][coupl_efthbasis_tCzz]= -(pow(e,2)/(4*cw));
        res[coupl_hbasistrip_tKz][coupl_efthbasis_tCza]= -((pow(e,2)*(-(sqrt(sw)/sqrt(cw)) + (2*pow(sw,3.0/2))/sqrt(cw)))/(4*sqrt(cw)*sqrt(sw)));
        res[coupl_hbasistrip_tKz][coupl_efthbasis_tCaa]= (pow(e,2)*sw)/4;
        res[coupl_hbasistrip_dg1z][coupl_efthbasis_Czbx]= -(pow(e,2)/(2*sw*(cw - sw)));
        res[coupl_hbasistrip_dg1z][coupl_efthbasis_Czz]= -(pow(e,2)/(2*cw*(cw - sw)));
        res[coupl_hbasistrip_dg1z][coupl_efthbasis_Cza]= pow(e,2)/(2*cw);
        res[coupl_hbasistrip_dg1z][coupl_efthbasis_Caa] = ((pow(e,2)*sw)/(2*(cw - sw)));
      }
      else{
        res.assign(nHiggsBasis_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
        res[coupl_hbasis_dCz][coupl_efthbasis_dCz]= 1.0;
        res[coupl_hbasis_Czz][coupl_efthbasis_Czz]= 1.0;
        res[coupl_hbasis_Czbx][coupl_efthbasis_Czbx]= 1.0;
        res[coupl_hbasis_tCzz][coupl_efthbasis_tCzz]= 1.0;
        res[coupl_hbasis_dCw][coupl_efthbasis_dCz]= 1.0;
        res[coupl_hbasis_Cww][coupl_efthbasis_Czz]= 1.0;
        res[coupl_hbasis_Cww][coupl_efthbasis_Cza]= 2*sw;
        res[coupl_hbasis_Cww][coupl_efthbasis_Caa]= pow(sw,2);
        res[coupl_hbasis_Cwbx][coupl_efthbasis_Czbx]= pow(MW,2)/(pow(MZ,2)*(cw - sw));
        res[coupl_hbasis_Cwbx][coupl_efthbasis_Czz]= (sw*pow(MW,2))/(cw*pow(MZ,2)*(cw - sw));
        res[coupl_hbasis_Cwbx][coupl_efthbasis_Cza]= -((pow(MW,2)*sw)/(cw*pow(MZ,2)));
        res[coupl_hbasis_Cwbx][coupl_efthbasis_Caa]= -((pow(MW,2)*pow(sw,2))/(pow(MZ,2)*(cw - sw)));
        res[coupl_hbasis_tCww][coupl_efthbasis_tCzz]= 1.0;
        res[coupl_hbasis_tCww][coupl_efthbasis_tCza]= 2*sw;
        res[coupl_hbasis_tCww][coupl_efthbasis_tCaa]= pow(sw,2);
        res[coupl_hbasis_Cza][coupl_efthbasis_Cza]= 1.0;
        res[coupl_hbasis_tCza][coupl_efthbasis_tCza]= 1.0;
        res[coupl_hbasis_Cabx][coupl_efthbasis_Czbx]= (2*cw)/(cw - sw);
        res[coupl_hbasis_Cabx][coupl_efthbasis_Czz]= 1.0/(cw - sw);
        res[coupl_hbasis_Cabx][coupl_efthbasis_Cza]= -1.0;
        res[coupl_hbasis_Cabx][coupl_efthbasis_Caa]= -((cw*sw)/(cw - sw));
        res[coupl_hbasis_Caa][coupl_efthbasis_Caa]= 1.0;
        res[coupl_hbasis_tCaa][coupl_efthbasis_tCaa]= 1.0;
        res[coupl_hbasis_Cgg][coupl_efthbasis_Cgg]= 1.0;
        res[coupl_hbasis_tCgg][coupl_efthbasis_tCgg]= 1.0;
      }
    }
    else if (basis_output == bEFT_HiggsBasis){
      if (include_triple_quartic_gauge){
        res.assign(nEFT_HiggsBasis_Include_Triple_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
        res[coupl_efthbasistrip_dCz][coupl_efthbasis_dCz]= 1.0;
        res[coupl_efthbasistrip_Czz][coupl_efthbasis_Czz]= 1.0;
        res[coupl_efthbasistrip_Czbx][coupl_efthbasis_Czbx]= 1.0;
        res[coupl_efthbasistrip_tCzz][coupl_efthbasis_tCzz]= 1.0;
        res[coupl_efthbasistrip_Cza][coupl_efthbasis_Cza]= 1.0;
        res[coupl_efthbasistrip_tCza][coupl_efthbasis_tCza]= 1.0;
        res[coupl_efthbasistrip_Caa][coupl_efthbasis_Caa]= 1.0;
        res[coupl_efthbasistrip_tCaa][coupl_efthbasis_tCaa]= 1.0;
        res[coupl_efthbasistrip_Cgg][coupl_efthbasis_Cgg]= 1.0;
        res[coupl_efthbasistrip_tCgg][coupl_efthbasis_tCgg]= 1.0;
        res[coupl_efthbasistrip_dKa][coupl_efthbasis_Czz]= (pow(e,2)/(4*sw));
        res[coupl_efthbasistrip_dKa][coupl_efthbasis_Cza]= (pow(e,2)*(2.0*sw-1.0))/(4.0*sw);
        res[coupl_efthbasistrip_dKa][coupl_efthbasis_Caa]= -(1.0/4)*cw*pow(e,2);
        res[coupl_efthbasistrip_tKa][coupl_efthbasis_tCzz]= pow(e,2)/(4*sw);
        res[coupl_efthbasistrip_tKa][coupl_efthbasis_tCza]= -((pow(e,2)*(sqrt(cw)/sqrt(sw) - 2*sqrt(cw)*sqrt(sw)))/(4*sqrt(cw)*sqrt(sw)));
        res[coupl_efthbasistrip_tKa][coupl_efthbasis_tCaa]= -(1.0/4)*cw*pow(e,2);
        res[coupl_efthbasistrip_dKz][coupl_efthbasis_Czbx]= -(pow(e,2)/(4*sw*(cw - sw)));
        res[coupl_efthbasistrip_dKz][coupl_efthbasis_Czz]= (pow(e,2)/(2*(-cw + sw)));
        res[coupl_efthbasistrip_dKz][coupl_efthbasis_Cza]= pow(e,2)/2;
        res[coupl_efthbasistrip_dKz][coupl_efthbasis_Caa]= -(cw*pow(e,2)*sw)/(2*(-cw + sw));
        res[coupl_efthbasistrip_tKz][coupl_efthbasis_tCzz]= -(pow(e,2)/(4*cw));
        res[coupl_efthbasistrip_tKz][coupl_efthbasis_tCza]= -((pow(e,2)*(-(sqrt(sw)/sqrt(cw)) + (2*pow(sw,3.0/2))/sqrt(cw)))/(4*sqrt(cw)*sqrt(sw)));
        res[coupl_efthbasistrip_tKz][coupl_efthbasis_tCaa]= (pow(e,2)*sw)/4;
        res[coupl_efthbasistrip_dg1z][coupl_efthbasis_Czbx]= -(pow(e,2)/(2*sw*(cw - sw)));
        res[coupl_efthbasistrip_dg1z][coupl_efthbasis_Czz]= -(pow(e,2)/(2*cw*(cw - sw)));
        res[coupl_efthbasistrip_dg1z][coupl_efthbasis_Cza]= pow(e,2)/(2*cw);
        res[coupl_efthbasistrip_dg1z][coupl_efthbasis_Caa] = ((pow(e,2)*sw)/(2*(cw - sw)));
      }
      else{
        res.assign(nEFT_HiggsBasis_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
        for (size_t i=0; i<(size_t) nEFT_HiggsBasis_CouplingTypes; i++) res.at(i).at(i)=1.0;
      }
    }
    else if (basis_output == bWarsawBasis){
      if (custodial_symmetry){
        res.assign(nWarsawBasis_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
	res[coupl_warsaw_cHbx][coupl_efthbasis_dCz]=(1.0/vev_lam)*1.0;
        res[coupl_warsaw_cHbx][coupl_efthbasis_Czbx]=(1.0/vev_lam)*(pow(e,2)*(-3.0*cw+sw))/(2.0*sw*(sw-cw));
        res[coupl_warsaw_cHbx][coupl_efthbasis_Czz]=(1.0/vev_lam)*-pow(e,2)/(sw-cw);
        res[coupl_warsaw_cHbx][coupl_efthbasis_Cza]=(1.0/vev_lam)*(-pow(e,2));
        res[coupl_warsaw_cHbx][coupl_efthbasis_Caa]=(1.0/vev_lam)*(cw*sw*pow(e,2))/(sw-cw);
        res[coupl_warsaw_cHG][coupl_efthbasis_Cgg]=(1.0/vev_lam)/(4.0);
        res[coupl_warsaw_cHW][coupl_efthbasis_Czz]=(1.0/vev_lam)*pow(e,2)/(4.0*sw);
        res[coupl_warsaw_cHW][coupl_efthbasis_Cza]=(1.0/vev_lam)*pow(e,2.0)/2.0;
        res[coupl_warsaw_cHW][coupl_efthbasis_Caa]=(1.0/vev_lam)*(sw*pow(e,2.0))/4.0;
        res[coupl_warsaw_cHWB][coupl_efthbasis_Czz]=(1.0/vev_lam)*pow(e,2.0)/(2.0*sqrt(cw)*sqrt(sw)) ;
        res[coupl_warsaw_cHWB][coupl_efthbasis_Cza]=(1.0/vev_lam)*(pow(e,2.0)*(sw-cw))/(2.0*sqrt(cw)*sqrt(sw)) ;
        res[coupl_warsaw_cHWB][coupl_efthbasis_Caa]=(1.0/vev_lam)*-(sqrt(cw)*sqrt(sw)*pow(e,2.0))/2.0 ;
        res[coupl_warsaw_cHB][coupl_efthbasis_Czz]=(1.0/vev_lam)*pow(e,2)/(4.0*cw);
        res[coupl_warsaw_cHB][coupl_efthbasis_Cza]=(1.0/vev_lam)*-pow(e,2)/2.0;
        res[coupl_warsaw_cHB][coupl_efthbasis_Caa]=(1.0/vev_lam)*(cw*pow(e,2))/(4.0);
	res[coupl_warsaw_cHD][coupl_efthbasis_Czbx]=(1.0/vev_lam)*(-2.0*pow(e,2))/(cw-sw);
        res[coupl_warsaw_cHD][coupl_efthbasis_Czz]=(1.0/vev_lam)*(-2.0*pow(e,2))/(cw-sw);
        res[coupl_warsaw_cHD][coupl_efthbasis_Cza]=(1.0/vev_lam)*(2.0*pow(e,2));
        res[coupl_warsaw_cHD][coupl_efthbasis_Caa]=(1.0/vev_lam)*(2.0*cw*sw*pow(e,2))/(cw-sw);
        res[coupl_warsaw_tcHB][coupl_efthbasis_tCzz]=(1.0/vev_lam)*pow(e,2)/(4.0*cw);
        res[coupl_warsaw_tcHB][coupl_efthbasis_tCza]=(1.0/vev_lam)*-pow(e,2)/2.0;
        res[coupl_warsaw_tcHB][coupl_efthbasis_tCaa]=(1.0/vev_lam)*(cw*pow(e,2))/4.0;
        res[coupl_warsaw_tcHW][coupl_efthbasis_tCzz]=(1.0/vev_lam)*pow(e,2)/(4.0*sw);
        res[coupl_warsaw_tcHW][coupl_efthbasis_tCza]=(1.0/vev_lam)*pow(e,2)/(2.0);
        res[coupl_warsaw_tcHW][coupl_efthbasis_tCaa]=(1.0/vev_lam)*(pow(e,2)*sw)/4.0;
        res[coupl_warsaw_tcHWB][coupl_efthbasis_tCzz]=(1.0/vev_lam)*pow(e,2)/(2.0*sqrt(sw)*sqrt(cw));
        res[coupl_warsaw_tcHWB][coupl_efthbasis_tCza]=(1.0/vev_lam)*(-pow(e,2)*(cw-sw))/(2.0*sqrt(sw)*sqrt(cw));
        res[coupl_warsaw_tcHWB][coupl_efthbasis_tCaa]=(1.0/vev_lam)*-(pow(e,2)*sqrt(sw)*sqrt(cw))/2.0;
        res[coupl_warsaw_tcHG][coupl_efthbasis_tCgg]=(1.0/vev_lam)/4.0;
      }
      else{
        res.assign(nWarsawBasis_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
	res[coupl_warsaw_cHbx][coupl_efthbasis_dCz]=(1.0/vev_lam)*1.0;
        res[coupl_warsaw_cHbx][coupl_efthbasis_Czbx]=(1.0/vev_lam)*(pow(e,2)*(-3.0*cw+sw))/(2.0*sw*(sw-cw));
        res[coupl_warsaw_cHbx][coupl_efthbasis_Czz]=(1.0/vev_lam)*-pow(e,2)/(sw-cw);
	res[coupl_warsaw_cHbx][coupl_efthbasis_Cza]=(1.0/vev_lam)*(-pow(e,2));
        res[coupl_warsaw_cHbx][coupl_efthbasis_Caa]=(1.0/vev_lam)*(cw*sw*pow(e,2))/(sw-cw);
        res[coupl_warsaw_cHG][coupl_efthbasis_Cgg]=(1.0/vev_lam)/(4.0);
        res[coupl_warsaw_cHW][coupl_efthbasis_Czz]=(1.0/vev_lam)*pow(e,2)/(4.0*sw);
        res[coupl_warsaw_cHW][coupl_efthbasis_Cza]=(1.0/vev_lam)*pow(e,2.0)/2.0;
        res[coupl_warsaw_cHW][coupl_efthbasis_Caa]=(1.0/vev_lam)*(sw*pow(e,2.0))/4.0;
        res[coupl_warsaw_cHWB][coupl_efthbasis_Czz]=(1.0/vev_lam)*pow(e,2.0)/(2.0*sqrt(cw)*sqrt(sw)) ;
        res[coupl_warsaw_cHWB][coupl_efthbasis_Cza]=(1.0/vev_lam)*(pow(e,2.0)*(sw-cw))/(2.0*sqrt(cw)*sqrt(sw)) ;
        res[coupl_warsaw_cHWB][coupl_efthbasis_Caa]=(1.0/vev_lam)*-(sqrt(cw)*sqrt(sw)*pow(e,2.0))/2.0 ;
        res[coupl_warsaw_cHB][coupl_efthbasis_Czz]=(1.0/vev_lam)*pow(e,2)/(4.0*cw);
        res[coupl_warsaw_cHB][coupl_efthbasis_Cza]=(1.0/vev_lam)*-pow(e,2)/2.0;
        res[coupl_warsaw_cHB][coupl_efthbasis_Caa]=(1.0/vev_lam)*(cw*pow(e,2))/(4.0);
        res[coupl_warsaw_cHD][coupl_efthbasis_Czbx]=(1.0/vev_lam)*(-2.0*pow(e,2))/(cw-sw);
	res[coupl_warsaw_cHD][coupl_efthbasis_Czz]=(1.0/vev_lam)*(-2.0*pow(e,2))/(cw-sw);
        res[coupl_warsaw_cHD][coupl_efthbasis_Cza]=(1.0/vev_lam)*(2.0*pow(e,2));
        res[coupl_warsaw_cHD][coupl_efthbasis_Caa]=(1.0/vev_lam)*(2.0*cw*sw*pow(e,2))/(cw-sw);
        res[coupl_warsaw_tcHB][coupl_efthbasis_tCzz]=(1.0/vev_lam)*pow(e,2)/(4.0*cw);
        res[coupl_warsaw_tcHB][coupl_efthbasis_tCza]=(1.0/vev_lam)*-pow(e,2)/2.0;
        res[coupl_warsaw_tcHB][coupl_efthbasis_tCaa]=(1.0/vev_lam)*(cw*pow(e,2))/4.0;
        res[coupl_warsaw_tcHW][coupl_efthbasis_tCzz]=(1.0/vev_lam)*pow(e,2)/(4.0*sw);
        res[coupl_warsaw_tcHW][coupl_efthbasis_tCza]=(1.0/vev_lam)*pow(e,2)/(2.0);
        res[coupl_warsaw_tcHW][coupl_efthbasis_tCaa]=(1.0/vev_lam)*(pow(e,2)*sw)/4.0;
        res[coupl_warsaw_tcHWB][coupl_efthbasis_tCzz]=(1.0/vev_lam)*pow(e,2)/(2.0*sqrt(sw)*sqrt(cw));
        res[coupl_warsaw_tcHWB][coupl_efthbasis_tCza]=(1.0/vev_lam)*(-pow(e,2)*(cw-sw))/(2.0*sqrt(sw)*sqrt(cw));
        res[coupl_warsaw_tcHWB][coupl_efthbasis_tCaa]=(1.0/vev_lam)*-(pow(e,2)*sqrt(sw)*sqrt(cw))/2.0;
        res[coupl_warsaw_tcHG][coupl_efthbasis_tCgg]=(1.0/vev_lam)/4.0;
      }
    }
  }
  else if (basis_input == bWarsawBasis){
    if (basis_output == bHiggsBasis){
      if (custodial_symmetry){
        if (include_triple_quartic_gauge){
          res.assign(nHiggsBasis_Include_Triple_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_hbasistrip_dCz][coupl_warsaw_cHbx]=(vev_lam)*1.0;
          res[coupl_hbasistrip_dCz][coupl_warsaw_cHWB]=(vev_lam)*3.0*pow(e,2)/sw;
          res[coupl_hbasistrip_dCz][coupl_warsaw_cHD]=(vev_lam)*(-(1.0/4.0)+(3.0*cw)/(4.0*sw));
          res[coupl_hbasistrip_Czz][coupl_warsaw_cHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_hbasistrip_Czz][coupl_warsaw_cHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_hbasistrip_Czz][coupl_warsaw_cHB]=(vev_lam)*(4.0*pow(sw,2)*cw)/pow(e,2);
          res[coupl_hbasistrip_Czbx][coupl_warsaw_cHWB]=(vev_lam)*-2.0;
          res[coupl_hbasistrip_Czbx][coupl_warsaw_cHD]=(vev_lam)*(sw-cw)/(2.0*pow(e,2));
          res[coupl_hbasistrip_tCzz][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw*pow(sw,2))/pow(e,2);
          res[coupl_hbasistrip_tCzz][coupl_warsaw_tcHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_hbasistrip_tCzz][coupl_warsaw_tcHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_hbasistrip_dCw][coupl_warsaw_cHbx]=(vev_lam)*1.0;
          res[coupl_hbasistrip_dCw][coupl_warsaw_cHWB]=(vev_lam)*3.0*pow(e,2)/sw;
          res[coupl_hbasistrip_dCw][coupl_warsaw_cHD]=(vev_lam)*(-(1.0/4.0)+(3.0*cw)/(4.0*sw));
          res[coupl_hbasistrip_Cww][coupl_warsaw_cHW]=(vev_lam)*(4.0*sw)/pow(e,2);
          res[coupl_hbasistrip_Cwbx][coupl_warsaw_cHD]=(vev_lam)*-(pow(MW,2.0)/(2.0*pow(e,2)*pow(MZ,2)*(cw-sw)))*(1-2*sw);
          res[coupl_hbasistrip_Cwbx][coupl_warsaw_cHWB]=(vev_lam)*-(pow(MW,2.0)/(2.0*pow(e,2)*pow(MZ,2)*(cw-sw))) * (4.0*(pow(e,2)-(pow(sw,3.0/2.0)/sqrt(cw))));
          res[coupl_hbasistrip_tCww][coupl_warsaw_tcHW]=(vev_lam)*4.0*sw/pow(e,2);
          res[coupl_hbasistrip_Cza][coupl_warsaw_cHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_hbasistrip_Cza][coupl_warsaw_cHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_hbasistrip_Cza][coupl_warsaw_cHB]=(vev_lam)*-(4.0*cw*sw)/pow(e,2);
          res[coupl_hbasistrip_tCza][coupl_warsaw_tcHB]=(vev_lam)*-(4.0*sw*cw)/pow(e,2);
          res[coupl_hbasistrip_tCza][coupl_warsaw_tcHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_hbasistrip_tCza][coupl_warsaw_tcHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_hbasistrip_Caa][coupl_warsaw_cHW]=(vev_lam)*(4.0*sw)/(pow(e,2));
          res[coupl_hbasistrip_Caa][coupl_warsaw_cHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/(pow(e,2));
          res[coupl_hbasistrip_Caa][coupl_warsaw_cHB]=(vev_lam)*(4.0*cw)/(pow(e,2));
          res[coupl_hbasistrip_Cabx][coupl_warsaw_cHWB]=(vev_lam)*(4.0*cw*pow(e,2)-2.0*sqrt(sw)*sqrt(cw))/(pow(e,2)*(sw-cw));
          res[coupl_hbasistrip_Cabx][coupl_warsaw_cHD]=(vev_lam)* (1.0-3.0*sw+2.0*pow(sw,2.0))/(pow(e,2)*(sw-cw));
          res[coupl_hbasistrip_tCaa][coupl_warsaw_tcHW]=(vev_lam)*(4.0*sw)/pow(e,2);
          res[coupl_hbasistrip_tCaa][coupl_warsaw_tcHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/pow(e,2);
          res[coupl_hbasistrip_tCaa][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw)/pow(e,2);
          res[coupl_hbasistrip_Cgg][coupl_warsaw_cHG]=(vev_lam)*4.0;
          res[coupl_hbasistrip_tCgg][coupl_warsaw_tcHG]=(vev_lam)*4.0;
          res[coupl_hbasistrip_dKa][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(cw))/(2.0*sqrt(sw));
          res[coupl_hbasistrip_tKa][coupl_warsaw_tcHWB]=(vev_lam)*(sqrt(cw))/(2.0*sqrt(sw));
          res[coupl_hbasistrip_dKz][coupl_warsaw_cHWB]=(vev_lam)*-(sqrt(cw)*pow(e,2)-2*pow(sw,3.0/2.0)+2*pow(sw,5.0/2.0))/(2.0*sw*sqrt(cw)*(sw-cw));
          res[coupl_hbasistrip_dKz][coupl_warsaw_cHD]=(vev_lam)*1.0/(8.0*sw);
          res[coupl_hbasistrip_tKz][coupl_warsaw_tcHWB]=(vev_lam)*(-sqrt(sw))/(2.0*sqrt(cw));
          res[coupl_hbasistrip_dg1z][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(cw)*pow(e,2)-pow(sw,3.0/2.0))/(sqrt(cw)*sw*(cw-sw));
          res[coupl_hbasistrip_dg1z][coupl_warsaw_cHD]=(vev_lam)*1.0/(4.0*sw);
        }
        else{
          res.assign(nHiggsBasis_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_hbasis_dCz][coupl_warsaw_cHbx]=(vev_lam)*1.0;
          res[coupl_hbasis_dCz][coupl_warsaw_cHWB]=(vev_lam)*3.0*pow(e,2)/sw;
          res[coupl_hbasis_dCz][coupl_warsaw_cHD]=(vev_lam)*(-(1.0/4.0)+(3.0*cw)/(4.0*sw));
          res[coupl_hbasis_Czz][coupl_warsaw_cHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_hbasis_Czz][coupl_warsaw_cHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_hbasis_Czz][coupl_warsaw_cHB]=(vev_lam)*(4.0*pow(sw,2)*cw)/pow(e,2);
          res[coupl_hbasis_Czbx][coupl_warsaw_cHWB]=(vev_lam)*-2.0;
          res[coupl_hbasis_Czbx][coupl_warsaw_cHD]=(vev_lam)*(sw-cw)/(2.0*pow(e,2));
          res[coupl_hbasis_tCzz][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw*pow(sw,2))/pow(e,2);
          res[coupl_hbasis_tCzz][coupl_warsaw_tcHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_hbasis_tCzz][coupl_warsaw_tcHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_hbasis_dCw][coupl_warsaw_cHbx]=(vev_lam)*1.0;
          res[coupl_hbasis_dCw][coupl_warsaw_cHWB]=(vev_lam)*(3.0*pow(e,2))/sw;
          res[coupl_hbasis_dCw][coupl_warsaw_cHD]=(vev_lam)*(-(1.0/4.0)+(3.0*cw)/(4.0*sw));
          res[coupl_hbasis_Cww][coupl_warsaw_cHW]=(vev_lam)*1.0;
          res[coupl_hbasis_Cwbx][coupl_warsaw_cHB]=(vev_lam)*-2.0;
          res[coupl_hbasis_Cwbx][coupl_warsaw_cHD]=(vev_lam)*-cw/(2*pow(e,2));
          res[coupl_hbasis_tCww][coupl_warsaw_tcHW]=(vev_lam)*4.0*sw/pow(e,2);
          res[coupl_hbasis_Cza][coupl_warsaw_cHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_hbasis_Cza][coupl_warsaw_cHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_hbasis_Cza][coupl_warsaw_cHB]=(vev_lam)*-(4.0*cw*sw)/pow(e,2);
          res[coupl_hbasis_tCza][coupl_warsaw_tcHB]=(vev_lam)*-(4.0*sw*cw)/pow(e,2);
          res[coupl_hbasis_tCza][coupl_warsaw_tcHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_hbasis_tCza][coupl_warsaw_tcHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_hbasis_Caa][coupl_warsaw_cHW]=(vev_lam)*(4.0*sw)/(pow(e,2));
          res[coupl_hbasis_Caa][coupl_warsaw_cHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/(pow(e,2));
          res[coupl_hbasis_Caa][coupl_warsaw_cHB]=(vev_lam)*(4.0*cw)/(pow(e,2));
          res[coupl_hbasis_Cabx][coupl_warsaw_cHWB]=(vev_lam)*(4.0*cw*pow(e,2)-2.0*sqrt(sw)*sqrt(cw))/(pow(e,2)*(cw-sw));
          res[coupl_hbasis_Cabx][coupl_warsaw_cHD]=(vev_lam)* (1.0-3.0*sw+2.0*pow(sw,2.0))/(pow(e,2)*(cw-sw));
          res[coupl_hbasis_tCaa][coupl_warsaw_tcHW]=(vev_lam)*(4.0*sw)/pow(e,2);
          res[coupl_hbasis_tCaa][coupl_warsaw_tcHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/pow(e,2);
          res[coupl_hbasis_tCaa][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw)/pow(e,2);
          res[coupl_hbasis_Cgg][coupl_warsaw_cHG]=(vev_lam)*4.0;
          res[coupl_hbasis_tCgg][coupl_warsaw_tcHG]=(vev_lam)*4.0;
        }
      }
      else {
        if (include_triple_quartic_gauge){
          res.assign(nHiggsBasis_Include_Triple_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_hbasistrip_dCz][coupl_warsaw_cHbx]=(vev_lam)*1.0;
          res[coupl_hbasistrip_dCz][coupl_warsaw_cHD]=(vev_lam)*-(1.0/4.0);
          res[coupl_hbasistrip_Czz][coupl_warsaw_cHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_hbasistrip_Czz][coupl_warsaw_cHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_hbasistrip_Czz][coupl_warsaw_cHB]=(vev_lam)*(4.0*pow(sw,2)*cw)/pow(e,2);
          res[coupl_hbasistrip_Czbx][coupl_warsaw_cHD]=(vev_lam)*(sw)/(2.0*pow(e,2));
          res[coupl_hbasistrip_tCzz][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw*pow(sw,2))/pow(e,2);
          res[coupl_hbasistrip_tCzz][coupl_warsaw_tcHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_hbasistrip_tCzz][coupl_warsaw_tcHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_hbasistrip_dCw][coupl_warsaw_cHbx]=(vev_lam)*1.0;
          res[coupl_hbasistrip_dCw][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(sw)*sqrt(cw))/(cw-sw);
          res[coupl_hbasistrip_dCw][coupl_warsaw_cHD]=(vev_lam)*(sw/(4.0*(cw-sw))-2.0/(4.0*(cw-sw)));
          res[coupl_hbasistrip_Cww][coupl_warsaw_cHW]=(vev_lam)*(4.0*sw)/pow(e,2);
          res[coupl_hbasistrip_Cwbx][coupl_warsaw_cHD]=(vev_lam)*(pow(MW,2)*sw)/(2.0*pow(e,2)*pow(MZ,2)*(cw-sw));          
	  res[coupl_hbasistrip_Cwbx][coupl_warsaw_cHWB]=(vev_lam)*(2.0*pow(sw,3.0/2.0)*pow(MW,2))/(sqrt(cw)*(cw-sw)*pow(e,2)*pow(MZ,2));
          res[coupl_hbasistrip_tCww][coupl_warsaw_tcHW]=(vev_lam)*4.0*sw/pow(e,2);
          res[coupl_hbasistrip_Cza][coupl_warsaw_cHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_hbasistrip_Cza][coupl_warsaw_cHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_hbasistrip_Cza][coupl_warsaw_cHB]=(vev_lam)*-(4.0*cw*sw)/pow(e,2);
          res[coupl_hbasistrip_tCza][coupl_warsaw_tcHB]=(vev_lam)*-(4.0*sw*cw)/pow(e,2);
          res[coupl_hbasistrip_tCza][coupl_warsaw_tcHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_hbasistrip_tCza][coupl_warsaw_tcHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_hbasistrip_Caa][coupl_warsaw_cHW]=(vev_lam)*(4.0*sw)/(pow(e,2));
          res[coupl_hbasistrip_Caa][coupl_warsaw_cHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/(pow(e,2));
          res[coupl_hbasistrip_Caa][coupl_warsaw_cHB]=(vev_lam)*(4.0*cw)/(pow(e,2));
          res[coupl_hbasistrip_Cabx][coupl_warsaw_cHWB]=(vev_lam)*(-2.0*sqrt(sw)*sqrt(cw))/(pow(e,2)*(cw-sw));
          res[coupl_hbasistrip_Cabx][coupl_warsaw_cHD]=(vev_lam)* (-cw*sw)/(pow(e,2)*(cw-sw));
          res[coupl_hbasistrip_tCaa][coupl_warsaw_tcHW]=(vev_lam)*(4.0*sw)/pow(e,2);
          res[coupl_hbasistrip_tCaa][coupl_warsaw_tcHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/pow(e,2);
          res[coupl_hbasistrip_tCaa][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw)/pow(e,2);
          res[coupl_hbasistrip_Cgg][coupl_warsaw_cHG]=(vev_lam)*4.0;
          res[coupl_hbasistrip_tCgg][coupl_warsaw_tcHG]=(vev_lam)*4.0;
          res[coupl_hbasistrip_dKa][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(cw))/(2.0*sqrt(sw));
          res[coupl_hbasistrip_tKa][coupl_warsaw_tcHWB]=(vev_lam)*(sqrt(cw))/(2.0*sqrt(sw));
          res[coupl_hbasistrip_dKz][coupl_warsaw_cHWB]=(vev_lam)*1.0/2.0 * (8.0*sqrt(sw)*sqrt(cw))/(8.0*sw-4);
          res[coupl_hbasistrip_dKz][coupl_warsaw_cHD]=(vev_lam)*1.0/2.0 * 1.0/(8.0*sw - 4.0);
          res[coupl_hbasistrip_tKz][coupl_warsaw_tcHWB]=(vev_lam)*(-sqrt(sw))/(2.0*sqrt(cw));
          res[coupl_hbasistrip_dg1z][coupl_warsaw_cHWB]=(vev_lam)*(4.0*sqrt(sw))/(sqrt(cw)*(8.0*sw-4.0));
          res[coupl_hbasistrip_dg1z][coupl_warsaw_cHD]=(vev_lam)*1.0/(8.0*sw-4.0);
        }
        else{
          res.assign(nHiggsBasis_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_hbasis_dCz][coupl_warsaw_cHbx]=(vev_lam)*1.0;
          res[coupl_hbasis_dCz][coupl_warsaw_cHD]=(vev_lam)*-(1.0/4.0);
          res[coupl_hbasis_Czz][coupl_warsaw_cHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_hbasis_Czz][coupl_warsaw_cHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_hbasis_Czz][coupl_warsaw_cHB]=(vev_lam)*(4.0*pow(sw,2)*cw)/pow(e,2);
          res[coupl_hbasis_Czbx][coupl_warsaw_cHD]=(vev_lam)*(sw/(2.0*pow(e,2)));
          res[coupl_hbasis_tCzz][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw*pow(sw,2))/pow(e,2);
          res[coupl_hbasis_tCzz][coupl_warsaw_tcHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_hbasis_tCzz][coupl_warsaw_tcHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_hbasis_dCw][coupl_warsaw_cHbx]=(vev_lam)*1.0;
          res[coupl_hbasis_dCw][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(sw)*sqrt(cw))/(cw-sw);
          res[coupl_hbasis_dCw][coupl_warsaw_cHD]=(vev_lam)*(sw/(4.0*(cw-sw))-2.0/(4.0*(cw-sw)));
          res[coupl_hbasis_Cww][coupl_warsaw_cHW]=(vev_lam)*(4.0*sw)/pow(e,2);
          res[coupl_hbasis_Cwbx][coupl_warsaw_cHD]=(vev_lam)*(pow(MW,2)*sw)/(2.0*pow(e,2)*pow(MZ,2)*(cw-sw));
          res[coupl_hbasis_Cwbx][coupl_warsaw_cHWB]=(vev_lam)*(2.0*pow(sw,3.0/2.0)*pow(MW,2))/(sqrt(cw)*(cw-sw)*pow(e,2)*pow(MZ,2));
          res[coupl_hbasis_tCww][coupl_warsaw_tcHW]=(vev_lam)*4.0*sw/pow(e,2);
          res[coupl_hbasis_Cza][coupl_warsaw_cHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_hbasis_Cza][coupl_warsaw_cHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_hbasis_Cza][coupl_warsaw_cHB]=(vev_lam)*-(4.0*cw*sw)/pow(e,2);
          res[coupl_hbasis_tCza][coupl_warsaw_tcHB]=(vev_lam)*-(4.0*sw*cw)/pow(e,2);
          res[coupl_hbasis_tCza][coupl_warsaw_tcHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_hbasis_tCza][coupl_warsaw_tcHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_hbasis_Caa][coupl_warsaw_cHW]=(vev_lam)*(4.0*sw)/(pow(e,2));
          res[coupl_hbasis_Caa][coupl_warsaw_cHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/(pow(e,2));
          res[coupl_hbasis_Caa][coupl_warsaw_cHB]=(vev_lam)*(4.0*cw)/(pow(e,2));
          res[coupl_hbasis_Cabx][coupl_warsaw_cHWB]=(vev_lam)*(-2.0*sqrt(sw)*sqrt(cw))/(pow(e,2)*(cw-sw));
          res[coupl_hbasis_Cabx][coupl_warsaw_cHD]=(vev_lam)* (-cw*sw)/(pow(e,2)*(cw-sw));
          res[coupl_hbasis_tCaa][coupl_warsaw_tcHW]=(vev_lam)*(4.0*sw)/pow(e,2);
          res[coupl_hbasis_tCaa][coupl_warsaw_tcHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/pow(e,2);
          res[coupl_hbasis_tCaa][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw)/pow(e,2);
          res[coupl_hbasis_Cgg][coupl_warsaw_cHG]=(vev_lam)*4.0;
          res[coupl_hbasis_tCgg][coupl_warsaw_tcHG]=(vev_lam)*4.0;
        }
      }
    } 
    if (basis_output == bEFT_HiggsBasis){
      if (custodial_symmetry){
        if (include_triple_quartic_gauge){
          res.assign(nEFT_HiggsBasis_Include_Triple_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_efthbasistrip_dCz][coupl_warsaw_cHbx]=(vev_lam)*1.0;
          res[coupl_efthbasistrip_dCz][coupl_warsaw_cHWB]=(vev_lam)*3.0*pow(e,2)/sw;
          res[coupl_efthbasistrip_dCz][coupl_warsaw_cHD]=(vev_lam)*(-(1.0/4.0)+(3.0*cw)/(4.0*sw));
          res[coupl_efthbasistrip_Czz][coupl_warsaw_cHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_efthbasistrip_Czz][coupl_warsaw_cHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_efthbasistrip_Czz][coupl_warsaw_cHB]=(vev_lam)*(4.0*pow(sw,2)*cw)/pow(e,2);
          res[coupl_efthbasistrip_Czbx][coupl_warsaw_cHWB]=(vev_lam)*-2.0;
          res[coupl_efthbasistrip_Czbx][coupl_warsaw_cHD]=(vev_lam)*(sw-cw)/(2.0*pow(e,2));
          res[coupl_efthbasistrip_tCzz][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw*pow(sw,2))/pow(e,2);
          res[coupl_efthbasistrip_tCzz][coupl_warsaw_tcHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_efthbasistrip_tCzz][coupl_warsaw_tcHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_efthbasistrip_Cza][coupl_warsaw_cHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_efthbasistrip_Cza][coupl_warsaw_cHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_efthbasistrip_Cza][coupl_warsaw_cHB]=(vev_lam)*-(4.0*cw*sw)/pow(e,2);
          res[coupl_efthbasistrip_tCza][coupl_warsaw_tcHB]=(vev_lam)*-(4.0*sw*cw)/pow(e,2);
          res[coupl_efthbasistrip_tCza][coupl_warsaw_tcHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_efthbasistrip_tCza][coupl_warsaw_tcHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_efthbasistrip_Caa][coupl_warsaw_cHW]=(vev_lam)*(4.0*sw)/(pow(e,2));
          res[coupl_efthbasistrip_Caa][coupl_warsaw_cHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/(pow(e,2));
          res[coupl_efthbasistrip_Caa][coupl_warsaw_cHB]=(vev_lam)*(4.0*cw)/(pow(e,2));
          res[coupl_efthbasistrip_tCaa][coupl_warsaw_tcHW]=(vev_lam)*(4.0*sw)/pow(e,2);
          res[coupl_efthbasistrip_tCaa][coupl_warsaw_tcHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/pow(e,2);
          res[coupl_efthbasistrip_tCaa][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw)/pow(e,2);
          res[coupl_efthbasistrip_Cgg][coupl_warsaw_cHG]=(vev_lam)*4.0;
          res[coupl_efthbasistrip_tCgg][coupl_warsaw_tcHG]=(vev_lam)*4.0;
          res[coupl_efthbasistrip_dKa][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(cw))/(2.0*sqrt(sw));
          res[coupl_efthbasistrip_tKa][coupl_warsaw_tcHWB]=(vev_lam)*(sqrt(cw))/(2.0*sqrt(sw));
          res[coupl_efthbasistrip_dKz][coupl_warsaw_cHWB]=(vev_lam)*-(sqrt(cw)*pow(e,2)-2*pow(sw,3.0/2.0)+2*pow(sw,5.0/2.0))/(2.0*sw*sqrt(cw)*(sw-cw));
          res[coupl_efthbasistrip_dKz][coupl_warsaw_cHD]=(vev_lam)*1.0/(8.0*sw);
          res[coupl_efthbasistrip_tKz][coupl_warsaw_tcHWB]=(vev_lam)*(-sqrt(sw))/(2.0*sqrt(cw));
	  res[coupl_efthbasistrip_dg1z][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(cw)*pow(e,2)-pow(sw,3.0/2.0))/(sqrt(cw)*sw*(cw-sw));
          res[coupl_efthbasistrip_dg1z][coupl_warsaw_cHD]=(vev_lam)*1.0/(4.0*sw);
        }
        else{
          res.assign(nEFT_HiggsBasis_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_efthbasis_dCz][coupl_warsaw_cHbx]=(vev_lam)*1.0;
          res[coupl_efthbasis_dCz][coupl_warsaw_cHWB]=(vev_lam)*3.0*sqrt(cw)/sqrt(sw);
          res[coupl_efthbasis_dCz][coupl_warsaw_cHD]=(vev_lam)*(-(1.0/4.0)+(3.0*cw)/(4.0*sw));
          res[coupl_efthbasis_Czz][coupl_warsaw_cHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_efthbasis_Czz][coupl_warsaw_cHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_efthbasis_Czz][coupl_warsaw_cHB]=(vev_lam)*(4.0*pow(sw,2)*cw)/pow(e,2);
          res[coupl_efthbasis_Czbx][coupl_warsaw_cHWB]=(vev_lam)*-2.0*(sqrt(sw)/(pow(e,2)));
          res[coupl_efthbasis_Czbx][coupl_warsaw_cHD]=(vev_lam)*(sw-cw)/(2.0*sqrt(cw)*pow(e,2));
          res[coupl_efthbasis_tCzz][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw*pow(sw,2))/pow(e,2);
          res[coupl_efthbasis_tCzz][coupl_warsaw_tcHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_efthbasis_tCzz][coupl_warsaw_tcHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_efthbasis_Cza][coupl_warsaw_cHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_efthbasis_Cza][coupl_warsaw_cHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_efthbasis_Cza][coupl_warsaw_cHB]=(vev_lam)*-(4.0*cw*sw)/pow(e,2);
          res[coupl_efthbasis_tCza][coupl_warsaw_tcHB]=(vev_lam)*-(4.0*sw*cw)/pow(e,2);
          res[coupl_efthbasis_tCza][coupl_warsaw_tcHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_efthbasis_tCza][coupl_warsaw_tcHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_efthbasis_Caa][coupl_warsaw_cHW]=(vev_lam)*(4.0*sw)/(pow(e,2));
          res[coupl_efthbasis_Caa][coupl_warsaw_cHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/(pow(e,2));
          res[coupl_efthbasis_Caa][coupl_warsaw_cHB]=(vev_lam)*(4.0*cw)/(pow(e,2));
          res[coupl_efthbasis_tCaa][coupl_warsaw_tcHW]=(vev_lam)*(4.0*sw)/pow(e,2);
          res[coupl_efthbasis_tCaa][coupl_warsaw_tcHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/pow(e,2);
          res[coupl_efthbasis_tCaa][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw)/pow(e,2);
          res[coupl_efthbasis_Cgg][coupl_warsaw_cHG]=(vev_lam)*4.0;
          res[coupl_efthbasis_tCgg][coupl_warsaw_tcHG]=(vev_lam)*4.0;
        }
      }
      else {    
        if (include_triple_quartic_gauge){
          res.assign(nEFT_HiggsBasis_Include_Triple_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_efthbasistrip_dCz][coupl_warsaw_cHbx]=(vev_lam)*1.0;
          res[coupl_efthbasistrip_dCz][coupl_warsaw_cHD]=(vev_lam)*-(1.0/4.0);
          res[coupl_efthbasistrip_Czz][coupl_warsaw_cHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_efthbasistrip_Czz][coupl_warsaw_cHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_efthbasistrip_Czz][coupl_warsaw_cHB]=(vev_lam)*(4.0*pow(sw,2)*cw)/pow(e,2);
          res[coupl_efthbasistrip_Czbx][coupl_warsaw_cHD]=(vev_lam)*(sw/(2.0*pow(e,2)));
          res[coupl_efthbasistrip_tCzz][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw*pow(sw,2))/pow(e,2);
          res[coupl_efthbasistrip_tCzz][coupl_warsaw_tcHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_efthbasistrip_tCzz][coupl_warsaw_tcHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_efthbasistrip_Cza][coupl_warsaw_cHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_efthbasistrip_Cza][coupl_warsaw_cHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_efthbasistrip_Cza][coupl_warsaw_cHB]=(vev_lam)*-(4.0*cw*sw)/pow(e,2);
          res[coupl_efthbasistrip_tCza][coupl_warsaw_tcHB]=(vev_lam)*-(4.0*sw*cw)/pow(e,2);
          res[coupl_efthbasistrip_tCza][coupl_warsaw_tcHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_efthbasistrip_tCza][coupl_warsaw_tcHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_efthbasistrip_Caa][coupl_warsaw_cHW]=(vev_lam)*(4.0*sw)/(pow(e,2));
          res[coupl_efthbasistrip_Caa][coupl_warsaw_cHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/(pow(e,2));
          res[coupl_efthbasistrip_Caa][coupl_warsaw_cHB]=(vev_lam)*(4.0*cw)/(pow(e,2));
          res[coupl_efthbasistrip_tCaa][coupl_warsaw_tcHW]=(vev_lam)*(4.0*sw)/pow(e,2);
          res[coupl_efthbasistrip_tCaa][coupl_warsaw_tcHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/pow(e,2);
          res[coupl_efthbasistrip_tCaa][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw)/pow(e,2);
          res[coupl_efthbasistrip_Cgg][coupl_warsaw_cHG]=(vev_lam)*4.0;
          res[coupl_efthbasistrip_tCgg][coupl_warsaw_tcHG]=(vev_lam)*4.0;
          res[coupl_efthbasistrip_dKa][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(cw))/(2.0*sqrt(sw));
          res[coupl_efthbasistrip_tKa][coupl_warsaw_tcHWB]=(vev_lam)*(sqrt(cw))/(2.0*sqrt(sw));
          res[coupl_efthbasistrip_dKz][coupl_warsaw_cHWB]=(vev_lam)*1.0/2.0 * (8.0*sqrt(sw)*sqrt(cw))/(8.0*sw-4);
          res[coupl_efthbasistrip_dKz][coupl_warsaw_cHD]=(vev_lam)*1.0/2.0 * 1.0/(8.0*sw - 4.0);
          res[coupl_efthbasistrip_tKz][coupl_warsaw_tcHWB]=(vev_lam)*(-sqrt(sw))/(2.0*sqrt(cw));
          res[coupl_efthbasistrip_dg1z][coupl_warsaw_cHWB]=(vev_lam)*(4.0*sqrt(sw))/(sqrt(cw)*(8.0*sw-4.0));
          res[coupl_efthbasistrip_dg1z][coupl_warsaw_cHD]=(vev_lam)*1.0/(8.0*sw-4.0);
	}
        else{
          res.assign(nEFT_HiggsBasis_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_efthbasis_dCz][coupl_warsaw_cHbx]=(vev_lam)*1.0;
          res[coupl_efthbasis_dCz][coupl_warsaw_cHD]=(vev_lam)*-(1.0/4.0);
          res[coupl_efthbasis_Czz][coupl_warsaw_cHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_efthbasis_Czz][coupl_warsaw_cHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_efthbasis_Czz][coupl_warsaw_cHB]=(vev_lam)*(4.0*pow(sw,2)*cw)/pow(e,2);
          res[coupl_efthbasis_Czbx][coupl_warsaw_cHD]=(vev_lam)*(sw)/(2.0*pow(e,2));
          res[coupl_efthbasis_tCzz][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw*pow(sw,2))/pow(e,2);
          res[coupl_efthbasis_tCzz][coupl_warsaw_tcHW]=(vev_lam)*(4.0*pow(cw,2)*sw)/pow(e,2);
          res[coupl_efthbasis_tCzz][coupl_warsaw_tcHWB]=(vev_lam)*(4.0*pow(cw,3.0/2.0)*pow(sw,3.0/2.0))/pow(e,2);
          res[coupl_efthbasis_Cza][coupl_warsaw_cHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_efthbasis_Cza][coupl_warsaw_cHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_efthbasis_Cza][coupl_warsaw_cHB]=(vev_lam)*-(4.0*cw*sw)/pow(e,2);
          res[coupl_efthbasis_tCza][coupl_warsaw_tcHB]=(vev_lam)*-(4.0*sw*cw)/pow(e,2);
          res[coupl_efthbasis_tCza][coupl_warsaw_tcHW]=(vev_lam)*(4.0*cw*sw)/pow(e,2);
          res[coupl_efthbasis_tCza][coupl_warsaw_tcHWB]=(vev_lam)*-(2.0*sqrt(cw)*sqrt(sw)*(cw-sw))/(pow(e,2));
          res[coupl_efthbasis_Caa][coupl_warsaw_cHW]=(vev_lam)*(4.0*sw)/(pow(e,2));
          res[coupl_efthbasis_Caa][coupl_warsaw_cHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/(pow(e,2));
          res[coupl_efthbasis_Caa][coupl_warsaw_cHB]=(vev_lam)*(4.0*cw)/(pow(e,2));
          res[coupl_efthbasis_tCaa][coupl_warsaw_tcHW]=(vev_lam)*(4.0*sw)/pow(e,2);
          res[coupl_efthbasis_tCaa][coupl_warsaw_tcHWB]=(vev_lam)*-(4.0*sqrt(cw)*sqrt(sw))/pow(e,2);
          res[coupl_efthbasis_tCaa][coupl_warsaw_tcHB]=(vev_lam)*(4.0*cw)/pow(e,2);
          res[coupl_efthbasis_Cgg][coupl_warsaw_cHG]=(vev_lam)*4.0;
          res[coupl_efthbasis_tCgg][coupl_warsaw_tcHG]=(vev_lam)*4.0;
        }
      }
    } 
    else if (basis_output == bAmplitude_JHUGen){
      if(custodial_symmetry){
        if (include_triple_quartic_gauge){
          res.assign(nAmplitude_JHUGen_Include_Triple_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_ampjhutrip_ghz1][coupl_warsaw_cHbx]=(vev_lam)*2.0;
          res[coupl_ampjhutrip_ghz1][coupl_warsaw_cHWB]=(vev_lam)*6.0*pow(e,2)/sw;
          res[coupl_ampjhutrip_ghz1][coupl_warsaw_cHD]=(vev_lam)*(-(2.0/4.0)+(6.0*cw)/(4.0*sw));
          res[coupl_ampjhutrip_ghz2][coupl_warsaw_cHB]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhutrip_ghz2][coupl_warsaw_cHW]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhutrip_ghz2][coupl_warsaw_cHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
          res[coupl_ampjhutrip_ghz1_prime2][coupl_warsaw_cHWB]=(vev_lam)*-(2.0*pow(e,2))/(sw);
          res[coupl_ampjhutrip_ghz1_prime2][coupl_warsaw_cHD]=(vev_lam)*(1.0-1.0/(2.0*sw));
          res[coupl_ampjhutrip_ghz4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhutrip_ghz4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhutrip_ghz4][coupl_warsaw_tcHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
          res[coupl_ampjhutrip_ghw1][coupl_warsaw_cHbx]=(vev_lam)*2.0;
          res[coupl_ampjhutrip_ghw1][coupl_warsaw_cHWB]=(vev_lam)*6.0*pow(e,2)/sw;
          res[coupl_ampjhutrip_ghw1][coupl_warsaw_cHD]=(vev_lam)*(-(2.0/4.0)+(6.0*cw)/(4.0*sw));
          res[coupl_ampjhutrip_ghw1_prime2][coupl_warsaw_cHWB]=(vev_lam)*-(pow(MW,2)/(2.0*pow(MZ,2)*sw*(1-2.0*sw))) * 4.0*(pow(e,2)-(pow(sw,3.0/2.0)/sqrt(cw)));
          res[coupl_ampjhutrip_ghw1_prime2][coupl_warsaw_cHD]=(vev_lam)*-(pow(MW,2)/(2.0*pow(MZ,2)*sw*(1-2.0*sw)))*(1-2.0*sw);
          res[coupl_ampjhutrip_ghw2][coupl_warsaw_cHW]=(vev_lam)*-2.0;
          res[coupl_ampjhutrip_ghw4][coupl_warsaw_tcHW]=(vev_lam)*-2.0;
          res[coupl_ampjhutrip_ghzgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_ampjhutrip_ghzgs2][coupl_warsaw_cHWB]=(vev_lam)*(cw-sw);
          res[coupl_ampjhutrip_ghzgs2][coupl_warsaw_cHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_ampjhutrip_ghzgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_ampjhutrip_ghzgs4][coupl_warsaw_tcHWB]=(vev_lam)*(cw-sw);
          res[coupl_ampjhutrip_ghzgs4][coupl_warsaw_tcHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_ampjhutrip_ghzgs1_prime2][coupl_warsaw_cHWB]=(vev_lam)*(2*sqrt(cw)*sqrt(sw)-4.0*cw*pow(e,2))/(sqrt(sw)*sqrt(cw)*(cw-sw));
          res[coupl_ampjhutrip_ghzgs1_prime2][coupl_warsaw_cHD]=(vev_lam)*(-2.0*pow(sw,2)+3.0*sw-1.0)/(sqrt(sw)*sqrt(cw)*(cw-sw));
          res[coupl_ampjhutrip_ghgsgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhutrip_ghgsgs2][coupl_warsaw_cHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_ampjhutrip_ghgsgs2][coupl_warsaw_cHB]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhutrip_ghgsgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhutrip_ghgsgs4][coupl_warsaw_tcHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_ampjhutrip_ghgsgs4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhutrip_ghg2][coupl_warsaw_cHG]=(vev_lam)*-2.0;
          res[coupl_ampjhutrip_ghg4][coupl_warsaw_tcHG]=(vev_lam)*-2.0;
          res[coupl_ampjhutrip_dV_Z][coupl_warsaw_cHD]=(vev_lam)*1.0/(4*sw);
          res[coupl_ampjhutrip_dV_Z][coupl_warsaw_cHWB]=(vev_lam)*-(sqrt(cw)*pow(e,2)-2.0*pow(sw,3.0/2.0)+2.0*pow(sw,5.0/2.0))/(sqrt(cw)*sw*(sw-cw));
          res[coupl_ampjhutrip_dV_A][coupl_warsaw_cHWB]=(vev_lam)*sqrt(cw)/sqrt(sw);
          res[coupl_ampjhutrip_dP_Z][coupl_warsaw_cHD]=(vev_lam)*1.0/(4.0*sw);
          res[coupl_ampjhutrip_dP_Z][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(cw)*pow(e,2)-pow(sw,3.0/2.0))/(sqrt(cw)*sw*(cw-sw));
          //res[coupl_ampjhutrip_dP_A][coupl_warsaw_tcHG]=(vev_lam)*0;
          res[coupl_ampjhutrip_dM_Z][coupl_warsaw_cHD]=(vev_lam)*1.0/(4.0*sw);
          res[coupl_ampjhutrip_dM_Z][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(cw)*pow(e,2)-pow(sw,3.0/2.0))/(sqrt(cw)*sw*(cw-sw));
          //res[coupl_ampjhutrip_dM_A][coupl_warsaw_tcHG]=(vev_lam)*0;
          res[coupl_ampjhutrip_dFour_Z][coupl_warsaw_tcHWB]=(vev_lam)*-(sqrt(sw)/sqrt(cw));
          res[coupl_ampjhutrip_dFour_A][coupl_warsaw_tcHWB]=(vev_lam)*sqrt(cw)/sqrt(sw);
          res[coupl_ampjhutrip_dZZWpWm][coupl_warsaw_cHD]=(vev_lam)*(-2.0*pow(sw,2)+3.0*sw-1.0)/(2.0*pow(sw,2)*(sw-cw));
          res[coupl_ampjhutrip_dZZWpWm][coupl_warsaw_cHWB]=(vev_lam)*(2.0*sqrt(cw)*(-sqrt(cw)*pow(e,2)+pow(sw,3.0/2.0)))/(pow(sw,2)*(sw-cw));
          res[coupl_ampjhutrip_dZAWpWm][coupl_warsaw_cHD]=(vev_lam)*(sqrt(cw))/(4.0*pow(sw,3.0/2.0));
          res[coupl_ampjhutrip_dZAWpWm][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(cw)*pow(e,2)-pow(sw,3.0/2.0))/(pow(sw,3.0/2.0)*(cw-sw));
          //res[coupl_ampjhutrip_dAAwpWm][coupl_warsaw_tcHG]=(vev_lam)*0;
        }
        else{
          res.assign(nAmplitude_JHUGen_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_ampjhu_ghz1][coupl_warsaw_cHbx]=(vev_lam)*2.0;
          res[coupl_ampjhu_ghz1][coupl_warsaw_cHWB]=(vev_lam)*6.0*pow(e,2)/sw;
          res[coupl_ampjhu_ghz1][coupl_warsaw_cHD]=(vev_lam)*(-(2.0/4.0)+(6.0*cw)/(4.0*sw));
          res[coupl_ampjhu_ghz2][coupl_warsaw_cHB]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhu_ghz2][coupl_warsaw_cHW]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhu_ghz2][coupl_warsaw_cHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
          res[coupl_ampjhu_ghz1_prime2][coupl_warsaw_cHWB]=(vev_lam)*-(2.0*pow(e,2))/(sw);
          res[coupl_ampjhu_ghz1_prime2][coupl_warsaw_cHD]=(vev_lam)*(1.0-1.0/(2.0*sw));
          res[coupl_ampjhu_ghz4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhu_ghz4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhu_ghz4][coupl_warsaw_tcHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
          res[coupl_ampjhu_ghw1][coupl_warsaw_cHbx]=(vev_lam)*2.0;
          res[coupl_ampjhu_ghw1][coupl_warsaw_cHWB]=(vev_lam)*6.0*pow(e,2)/sw;
          res[coupl_ampjhu_ghw1][coupl_warsaw_cHD]=(vev_lam)*(-(2.0/4.0)+(6.0*cw)/(4.0*sw));
          res[coupl_ampjhu_ghw1_prime2][coupl_warsaw_cHWB]=(vev_lam)*-(pow(MW,2)/(2.0*pow(MZ,2)*sw*(1-2.0*sw))) * 4.0*(pow(e,2)-(pow(sw,3.0/2.0)/sqrt(cw)));
          res[coupl_ampjhu_ghw1_prime2][coupl_warsaw_cHD]=(vev_lam)*-(pow(MW,2)/(2.0*pow(MZ,2)*sw*(1-2.0*sw)))*(1-2.0*sw);
          res[coupl_ampjhu_ghw2][coupl_warsaw_cHW]=(vev_lam)*-2.0;
          res[coupl_ampjhu_ghw4][coupl_warsaw_tcHW]=(vev_lam)*-2.0;
          res[coupl_ampjhu_ghzgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_ampjhu_ghzgs2][coupl_warsaw_cHWB]=(vev_lam)*(cw-sw);
          res[coupl_ampjhu_ghzgs2][coupl_warsaw_cHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_ampjhu_ghzgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_ampjhu_ghzgs4][coupl_warsaw_tcHWB]=(vev_lam)*(cw-sw);
          res[coupl_ampjhu_ghzgs4][coupl_warsaw_tcHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_ampjhu_ghzgs1_prime2][coupl_warsaw_cHWB]=(vev_lam)*(2*sqrt(cw)*sqrt(sw)-4.0*cw*pow(e,2))/(sqrt(sw)*sqrt(cw)*(cw-sw));
          res[coupl_ampjhu_ghzgs1_prime2][coupl_warsaw_cHD]=(vev_lam)*(-2.0*pow(sw,2)+3.0*sw-1.0)/(sqrt(sw)*sqrt(cw)*(cw-sw));
          res[coupl_ampjhu_ghgsgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhu_ghgsgs2][coupl_warsaw_cHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_ampjhu_ghgsgs2][coupl_warsaw_cHB]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhu_ghgsgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhu_ghgsgs4][coupl_warsaw_tcHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_ampjhu_ghgsgs4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhu_ghg2][coupl_warsaw_cHG]=(vev_lam)*-2.0;
          res[coupl_ampjhu_ghg4][coupl_warsaw_tcHG]=(vev_lam)*-2.0;
        } 
      }
      else{
        if (include_triple_quartic_gauge){
          res.assign(nAmplitude_JHUGen_Include_Triple_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_ampjhutrip_ghz1][coupl_warsaw_cHbx]=(vev_lam)*2.0;
          res[coupl_ampjhutrip_ghz1][coupl_warsaw_cHD]=(vev_lam)*-(2.0/4.0);
          res[coupl_ampjhutrip_ghz2][coupl_warsaw_cHB]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhutrip_ghz2][coupl_warsaw_cHW]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhutrip_ghz2][coupl_warsaw_cHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
          res[coupl_ampjhutrip_ghz1_prime2][coupl_warsaw_cHD]=(vev_lam)*(1.0/2.0);
          res[coupl_ampjhutrip_ghz4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhutrip_ghz4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhutrip_ghz4][coupl_warsaw_tcHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
          res[coupl_ampjhutrip_ghw1][coupl_warsaw_cHbx]=(vev_lam)*(4.0*(cw-sw))/(2.0*cw-sw);
          res[coupl_ampjhutrip_ghw1][coupl_warsaw_cHWB]=(vev_lam)*(4.0*sqrt(cw)*sqrt(cw))/(2.0*(cw-sw));
          res[coupl_ampjhutrip_ghw1][coupl_warsaw_cHD]=(vev_lam)*(sw-2.0*cw)/(2.0*(cw-sw));
          res[coupl_ampjhutrip_ghw1_prime2][coupl_warsaw_cHWB]=(vev_lam)*(2.0*sqrt(sw)*pow(MW,2))/(sqrt(cw)*(cw-sw)*pow(MZ,2));
          res[coupl_ampjhutrip_ghw1_prime2][coupl_warsaw_cHD]=(vev_lam)*(pow(MW,2))/(2.0*pow(MZ,2)*(cw-sw));
          res[coupl_ampjhutrip_ghw2][coupl_warsaw_cHW]=(vev_lam)*-2.0;
          res[coupl_ampjhutrip_ghw4][coupl_warsaw_tcHW]=(vev_lam)*-2.0;
          res[coupl_ampjhutrip_ghzgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_ampjhutrip_ghzgs2][coupl_warsaw_cHWB]=(vev_lam)*(cw-sw);
          res[coupl_ampjhutrip_ghzgs2][coupl_warsaw_cHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_ampjhutrip_ghzgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_ampjhutrip_ghzgs4][coupl_warsaw_tcHWB]=(vev_lam)*(cw-sw);
          res[coupl_ampjhutrip_ghzgs4][coupl_warsaw_tcHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_ampjhutrip_ghzgs1_prime2][coupl_warsaw_cHWB]=(vev_lam)*(-2.0)/(cw-sw);
          res[coupl_ampjhutrip_ghzgs1_prime2][coupl_warsaw_cHD]=(vev_lam)*(-sqrt(cw)*sqrt(sw))/(cw-sw);
          res[coupl_ampjhutrip_ghgsgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhutrip_ghgsgs2][coupl_warsaw_cHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_ampjhutrip_ghgsgs2][coupl_warsaw_cHB]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhutrip_ghgsgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhutrip_ghgsgs4][coupl_warsaw_tcHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_ampjhutrip_ghgsgs4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhutrip_ghg2][coupl_warsaw_cHG]=(vev_lam)*-2.0;
          res[coupl_ampjhutrip_ghg4][coupl_warsaw_tcHG]=(vev_lam)*-2.0;
          res[coupl_ampjhutrip_dV_Z][coupl_warsaw_cHD]=(vev_lam)*(2.0)/(16*sw-8.0);
          res[coupl_ampjhutrip_dV_Z][coupl_warsaw_cHWB]=(vev_lam)*(8.0*sqrt(cw)*sqrt(sw))/(16*sw-8.0);
          res[coupl_ampjhutrip_dV_A][coupl_warsaw_cHWB]=(vev_lam)*sqrt(cw)/sqrt(sw);
          res[coupl_ampjhutrip_dP_Z][coupl_warsaw_cHD]=(vev_lam)*(4.0*sqrt(sw))/(cw*(8.0*sw-4.0));
          res[coupl_ampjhutrip_dP_Z][coupl_warsaw_cHWB]=(vev_lam)*(1.0)/(8.0*sw-4.0);
          //res[coupl_ampjhutrip_dP_A][coupl_warsaw_tcHG]=(vev_lam)*0;
          res[coupl_ampjhutrip_dM_Z][coupl_warsaw_cHD]=(vev_lam)*1.0/(4.0*sw);
          res[coupl_ampjhutrip_dM_Z][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(cw)*pow(e,2)-pow(sw,3.0/2.0))/(sqrt(cw)*sw*(cw-sw));
          //res[coupl_ampjhutrip_dM_A][coupl_warsaw_tcHG]=(vev_lam)*0;
          res[coupl_ampjhutrip_dFour_Z][coupl_warsaw_tcHWB]=(vev_lam)*-(sqrt(sw)/sqrt(cw));
          res[coupl_ampjhutrip_dFour_A][coupl_warsaw_tcHWB]=(vev_lam)*sqrt(cw)/sqrt(sw);
          res[coupl_ampjhutrip_dZZWpWm][coupl_warsaw_cHD]=(vev_lam)*(2.0*cw)/(sw*(8.0-4.0));
          res[coupl_ampjhutrip_dZZWpWm][coupl_warsaw_cHWB]=(vev_lam)*(8.0*cw*sqrt(sw))/(sw*sqrt(cw)*(8.0*sw-4.0));
          res[coupl_ampjhutrip_dZAWpWm][coupl_warsaw_cHD]=(vev_lam)*(sqrt(sw)/sqrt(cw))*(1.0/(8.0*sw-4.0));
          res[coupl_ampjhutrip_dZAWpWm][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(sw)/sqrt(cw)) * (4.0*sqrt(sw))/(sqrt(cw)*(8.0*sw-4.0));
          //res[coupl_ampjhutrip_dAAwpWm][coupl_warsaw_tcHG]=(vev_lam)*0;
        } 
        else{
          res.assign(nAmplitude_JHUGen_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_ampjhu_ghz1][coupl_warsaw_cHbx]=(vev_lam)*2.0;
          res[coupl_ampjhu_ghz1][coupl_warsaw_cHD]=(vev_lam)*(-(2.0/4.0));
          res[coupl_ampjhu_ghz2][coupl_warsaw_cHB]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhu_ghz2][coupl_warsaw_cHW]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhu_ghz2][coupl_warsaw_cHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
          res[coupl_ampjhu_ghz1_prime2][coupl_warsaw_cHD]=(vev_lam)*(1.0/2.0);
          res[coupl_ampjhu_ghz4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhu_ghz4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhu_ghz4][coupl_warsaw_tcHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
	  res[coupl_ampjhu_ghw1][coupl_warsaw_cHbx]=(vev_lam)*(4.0*(cw-sw))/(2.0*cw-sw);
          res[coupl_ampjhu_ghw1][coupl_warsaw_cHWB]=(vev_lam)*(4.0*sqrt(cw)*sqrt(cw))/(2.0*(cw-sw));
          res[coupl_ampjhu_ghw1][coupl_warsaw_cHD]=(vev_lam)*(sw-2.0*cw)/(2.0*(cw-sw));
          res[coupl_ampjhu_ghw1][coupl_warsaw_cHD]=(vev_lam)*(-(2.0/4.0)+(6.0*cw)/(4.0*sw));
          res[coupl_ampjhu_ghw1_prime2][coupl_warsaw_cHWB]=(vev_lam)*-(pow(MW,2)/(2.0*pow(MZ,2)*sw*(1-2.0*sw))) * 4.0*(pow(e,2)-(pow(sw,3.0/2.0)/sqrt(cw)));
          res[coupl_ampjhu_ghw1_prime2][coupl_warsaw_cHD]=(vev_lam)*-(pow(MW,2)/(2.0*pow(MZ,2)*sw*(1-2.0*sw)))*(1-2.0*sw);
          res[coupl_ampjhu_ghw2][coupl_warsaw_cHW]=(vev_lam)*-2.0;
          res[coupl_ampjhu_ghw4][coupl_warsaw_tcHW]=(vev_lam)*-2.0;
          res[coupl_ampjhu_ghzgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_ampjhu_ghzgs2][coupl_warsaw_cHWB]=(vev_lam)*(cw-sw);
          res[coupl_ampjhu_ghzgs2][coupl_warsaw_cHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_ampjhu_ghzgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_ampjhu_ghzgs4][coupl_warsaw_tcHWB]=(vev_lam)*(cw-sw);
          res[coupl_ampjhu_ghzgs4][coupl_warsaw_tcHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_ampjhu_ghzgs1_prime2][coupl_warsaw_cHWB]=(vev_lam)*(-2.0)/(cw-sw);
          res[coupl_ampjhu_ghzgs1_prime2][coupl_warsaw_cHD]=(vev_lam)*(-sqrt(cw)*sqrt(sw))/(cw-sw);
          res[coupl_ampjhu_ghgsgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhu_ghgsgs2][coupl_warsaw_cHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_ampjhu_ghgsgs2][coupl_warsaw_cHB]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhu_ghgsgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sw;
          res[coupl_ampjhu_ghgsgs4][coupl_warsaw_tcHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_ampjhu_ghgsgs4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*cw;
          res[coupl_ampjhu_ghg2][coupl_warsaw_cHG]=(vev_lam)*-2.0;
          res[coupl_ampjhu_ghg4][coupl_warsaw_tcHG]=(vev_lam)*-2.0;
        }
      }
    }
    else if (basis_output == bEFT_JHUGen){
      if(custodial_symmetry){
        if (include_triple_quartic_gauge){
          res.assign(nEFT_JHUGen_Include_Triple_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_eftjhutrip_ghz1][coupl_warsaw_cHbx]=(vev_lam)*2.0;
          res[coupl_eftjhutrip_ghz1][coupl_warsaw_cHWB]=(vev_lam)*6.0*pow(e,2)/sw;
          res[coupl_eftjhutrip_ghz1][coupl_warsaw_cHD]=(vev_lam)*(-(2.0/4.0)+(6.0*cw)/(4.0*sw));
          res[coupl_eftjhutrip_ghz2][coupl_warsaw_cHB]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhutrip_ghz2][coupl_warsaw_cHW]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhutrip_ghz2][coupl_warsaw_cHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
          res[coupl_eftjhutrip_ghz1_prime2][coupl_warsaw_cHWB]=(vev_lam)*-(2.0*pow(e,2))/(sw);
          res[coupl_eftjhutrip_ghz1_prime2][coupl_warsaw_cHD]=(vev_lam)*(1.0-1.0/(2.0*sw));
          res[coupl_eftjhutrip_ghz4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhutrip_ghz4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhutrip_ghz4][coupl_warsaw_tcHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
          res[coupl_eftjhutrip_ghzgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_eftjhutrip_ghzgs2][coupl_warsaw_cHWB]=(vev_lam)*(cw-sw);
          res[coupl_eftjhutrip_ghzgs2][coupl_warsaw_cHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_eftjhutrip_ghzgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_eftjhutrip_ghzgs4][coupl_warsaw_tcHWB]=(vev_lam)*(cw-sw);
          res[coupl_eftjhutrip_ghzgs4][coupl_warsaw_tcHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_eftjhutrip_ghgsgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhutrip_ghgsgs2][coupl_warsaw_cHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_eftjhutrip_ghgsgs2][coupl_warsaw_cHB]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhutrip_ghgsgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhutrip_ghgsgs4][coupl_warsaw_tcHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_eftjhutrip_ghgsgs4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhutrip_ghg2][coupl_warsaw_cHG]=(vev_lam)*-2.0;
          res[coupl_eftjhutrip_ghg4][coupl_warsaw_tcHG]=(vev_lam)*-2.0;
          res[coupl_eftjhutrip_dV_Z][coupl_warsaw_cHD]=(vev_lam)*1.0/(4*sw);
          res[coupl_eftjhutrip_dV_Z][coupl_warsaw_cHWB]=(vev_lam)*-(sqrt(cw)*pow(e,2)-2.0*pow(sw,3.0/2.0)+2.0*pow(sw,5.0/2.0))/(sqrt(cw)*sw*(sw-cw));
          res[coupl_eftjhutrip_dV_A][coupl_warsaw_cHWB]=(vev_lam)*sqrt(cw)/sqrt(sw);
          res[coupl_eftjhutrip_dP_Z][coupl_warsaw_cHD]=(vev_lam)*1.0/(4.0*sw);
          res[coupl_eftjhutrip_dP_Z][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(cw)*pow(e,2)-pow(sw,3.0/2.0))/(sqrt(cw)*sw*(cw-sw));
	  //res[coupl_ampjhutrip_dP_A][coupl_warsaw_tcHG]=(vev_lam)*0;
          res[coupl_eftjhutrip_dM_Z][coupl_warsaw_cHD]=(vev_lam)*1.0/(4.0*sw);
          res[coupl_eftjhutrip_dM_Z][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(cw)*pow(e,2)-pow(sw,3.0/2.0))/(sqrt(cw)*sw*(cw-sw));
          //res[coupl_ampjhutrip_dM_A][coupl_warsaw_tcHG]=(vev_lam)*0;
	  res[coupl_eftjhutrip_dFour_Z][coupl_warsaw_tcHWB]=(vev_lam)*-(sqrt(sw)/sqrt(cw));
          res[coupl_eftjhutrip_dFour_A][coupl_warsaw_tcHWB]=(vev_lam)*sqrt(cw)/sqrt(sw);
          res[coupl_eftjhutrip_dZZWpWm][coupl_warsaw_cHD]=(vev_lam)*(2.0*cw)/(sw*(8.0-4.0));
          res[coupl_eftjhutrip_dZZWpWm][coupl_warsaw_cHWB]=(vev_lam)*(8.0*cw*sqrt(sw))/(sw*sqrt(cw)*(8.0*sw-4.0));
          res[coupl_eftjhutrip_dZAWpWm][coupl_warsaw_cHD]=(vev_lam)*(sqrt(sw)/sqrt(cw))*(1.0/(8.0*sw-4.0));
          res[coupl_eftjhutrip_dZAWpWm][coupl_warsaw_cHWB]=(vev_lam)*(sqrt(sw)/sqrt(cw)) * (4.0*sqrt(sw))/(sqrt(cw)*(8.0*sw-4.0));
          //res[coupl_ampjhutrip_dAAwpWm][coupl_warsaw_tcHG]=(vev_lam)*0; 
	}
        else{
          res.assign(nEFT_JHUGen_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_eftjhutrip_ghz1][coupl_warsaw_cHbx]=(vev_lam)*2.0;
          res[coupl_eftjhutrip_ghz1][coupl_warsaw_cHWB]=(vev_lam)*6.0*pow(e,2)/sw;
          res[coupl_eftjhutrip_ghz1][coupl_warsaw_cHD]=(vev_lam)*(-(2.0/4.0)+(6.0*cw)/(4.0*sw));
          res[coupl_eftjhutrip_ghz2][coupl_warsaw_cHB]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhutrip_ghz2][coupl_warsaw_cHW]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhutrip_ghz2][coupl_warsaw_cHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
          res[coupl_eftjhutrip_ghz1_prime2][coupl_warsaw_cHWB]=(vev_lam)*-(2.0*pow(e,2))/(sw);
          res[coupl_eftjhutrip_ghz1_prime2][coupl_warsaw_cHD]=(vev_lam)*(1.0-1.0/(2.0*sw));
          res[coupl_eftjhutrip_ghz4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhutrip_ghz4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhutrip_ghz4][coupl_warsaw_tcHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
          res[coupl_eftjhutrip_ghzgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_eftjhutrip_ghzgs2][coupl_warsaw_cHWB]=(vev_lam)*(cw-sw);
          res[coupl_eftjhutrip_ghzgs2][coupl_warsaw_cHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_eftjhutrip_ghzgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_eftjhutrip_ghzgs4][coupl_warsaw_tcHWB]=(vev_lam)*(cw-sw);
          res[coupl_eftjhutrip_ghzgs4][coupl_warsaw_tcHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_eftjhutrip_ghgsgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhutrip_ghgsgs2][coupl_warsaw_cHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_eftjhutrip_ghgsgs2][coupl_warsaw_cHB]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhutrip_ghgsgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhutrip_ghgsgs4][coupl_warsaw_tcHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_eftjhutrip_ghgsgs4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhutrip_ghg2][coupl_warsaw_cHG]=(vev_lam)*-2.0;
          res[coupl_eftjhutrip_ghg4][coupl_warsaw_tcHG]=(vev_lam)*-2.0;	
	}
      }
      else{
        if (include_triple_quartic_gauge){
          res.assign(nEFT_JHUGen_Include_Triple_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_eftjhutrip_ghz1][coupl_warsaw_cHbx]=(vev_lam)*2.0;
          res[coupl_eftjhutrip_ghz1][coupl_warsaw_cHD]=(vev_lam)*-(2.0/4.0);
          res[coupl_eftjhutrip_ghz2][coupl_warsaw_cHB]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhutrip_ghz2][coupl_warsaw_cHW]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhutrip_ghz2][coupl_warsaw_cHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
          res[coupl_eftjhutrip_ghz1_prime2][coupl_warsaw_cHD]=(vev_lam)*(1.0/2.0);
          res[coupl_eftjhutrip_ghz4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhutrip_ghz4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhutrip_ghz4][coupl_warsaw_tcHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
          res[coupl_eftjhutrip_ghzgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_eftjhutrip_ghzgs2][coupl_warsaw_cHWB]=(vev_lam)*(cw-sw);
          res[coupl_eftjhutrip_ghzgs2][coupl_warsaw_cHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_eftjhutrip_ghzgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_eftjhutrip_ghzgs4][coupl_warsaw_tcHWB]=(vev_lam)*(cw-sw);
          res[coupl_eftjhutrip_ghzgs4][coupl_warsaw_tcHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_eftjhutrip_ghgsgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhutrip_ghgsgs2][coupl_warsaw_cHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_eftjhutrip_ghgsgs2][coupl_warsaw_cHB]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhutrip_ghgsgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhutrip_ghgsgs4][coupl_warsaw_tcHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_eftjhutrip_ghgsgs4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhutrip_ghg2][coupl_warsaw_cHG]=(vev_lam)*-2.0;
          res[coupl_eftjhutrip_ghg4][coupl_warsaw_tcHG]=(vev_lam)*-2.0;
          res[coupl_eftjhutrip_dV_Z][coupl_warsaw_cHD]=(vev_lam)*1.0/(-4.0+8.0*sw);
          res[coupl_eftjhutrip_dV_Z][coupl_warsaw_cHWB]=(vev_lam)*(8.0*sqrt(cw)*sqrt(sw))/(-4.0+8.0*sw);
          res[coupl_eftjhutrip_dV_A][coupl_warsaw_cHWB]=(vev_lam)*sqrt(cw)/sqrt(sw);
          res[coupl_eftjhutrip_dP_Z][coupl_warsaw_cHD]=(vev_lam)*1.0/(4.0*(-1.0+2*sw));
          res[coupl_eftjhutrip_dP_Z][coupl_warsaw_cHWB]=(vev_lam)*4.0*sqrt(sw)/(sqrt(cw)*4*(-1.0+2*sw));
          //res[coupl_ampjhutrip_dP_A][coupl_warsaw_tcHG]=(vev_lam)*0;
          res[coupl_eftjhutrip_dM_Z][coupl_warsaw_cHD]=(vev_lam)*1.0/(4.0*(-1.0+2*sw));
          res[coupl_eftjhutrip_dM_Z][coupl_warsaw_cHWB]=(vev_lam)*4.0*sqrt(sw)/(sqrt(cw)*4*(-1.0+2*sw));
          //res[coupl_ampjhutrip_dM_A][coupl_warsaw_tcHG]=(vev_lam)*0;
          res[coupl_eftjhutrip_dFour_Z][coupl_warsaw_tcHWB]=(vev_lam)*-(sqrt(sw)/sqrt(cw));
          res[coupl_eftjhutrip_dFour_A][coupl_warsaw_tcHWB]=(vev_lam)*sqrt(cw)/sqrt(sw);
          res[coupl_eftjhutrip_dZZWpWm][coupl_warsaw_cHWB]=(vev_lam)*-(2*sqrt(cw))/(sqrt(sw)*(cw-sw));
          res[coupl_eftjhutrip_dZZWpWm][coupl_warsaw_cHD]=(vev_lam)*-(cw)/(2*sw*(cw-sw));
          res[coupl_eftjhutrip_dZAWpWm][coupl_warsaw_cHWB]=(vev_lam)*1/(sw-cw);
          res[coupl_eftjhutrip_dZAWpWm][coupl_warsaw_cHD]=(vev_lam)*sqrt(cw)/(4*sqrt(sw)*(sw-cw));
        }
        else{
	  res.assign(nEFT_JHUGen_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
          res[coupl_eftjhu_ghz1][coupl_warsaw_cHbx]=(vev_lam)*2.0;
          res[coupl_eftjhu_ghz1][coupl_warsaw_cHD]=(vev_lam)*(-(2.0/4.0));
          res[coupl_eftjhu_ghz2][coupl_warsaw_cHB]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhu_ghz2][coupl_warsaw_cHW]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhu_ghz2][coupl_warsaw_cHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
          res[coupl_eftjhu_ghz1_prime2][coupl_warsaw_cHD]=(vev_lam)*(1.0/2.0);
          res[coupl_eftjhu_ghz4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhu_ghz4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhu_ghz4][coupl_warsaw_tcHWB]=(vev_lam)*-2.0*sqrt(sw*cw);
          res[coupl_eftjhu_ghzgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_eftjhu_ghzgs2][coupl_warsaw_cHWB]=(vev_lam)*(cw-sw);
          res[coupl_eftjhu_ghzgs2][coupl_warsaw_cHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_eftjhu_ghzgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sqrt(cw*sw);
          res[coupl_eftjhu_ghzgs4][coupl_warsaw_tcHWB]=(vev_lam)*(cw-sw);
          res[coupl_eftjhu_ghzgs4][coupl_warsaw_tcHB]=(vev_lam)*2.0*sqrt(cw)*sqrt(sw);
          res[coupl_eftjhu_ghgsgs2][coupl_warsaw_cHW]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhu_ghgsgs2][coupl_warsaw_cHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_eftjhu_ghgsgs2][coupl_warsaw_cHB]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhu_ghgsgs4][coupl_warsaw_tcHW]=(vev_lam)*-2.0*sw;
          res[coupl_eftjhu_ghgsgs4][coupl_warsaw_tcHWB]=(vev_lam)*2.0*sqrt(cw*sw);
          res[coupl_eftjhu_ghgsgs4][coupl_warsaw_tcHB]=(vev_lam)*-2.0*cw;
          res[coupl_eftjhu_ghg2][coupl_warsaw_cHG]=(vev_lam)*-2.0;
          res[coupl_eftjhu_ghg4][coupl_warsaw_tcHG]=(vev_lam)*-2.0;	
        }    
      } 
    }
    else if (basis_output == bWarsawBasis){
      res.assign(nWarsawBasis_CouplingTypes, std::vector<double>(nWarsawBasis_CouplingTypes, 0));
      for (size_t i=0; i<(size_t) nWarsawBasis_CouplingTypes; i++) res.at(i).at(i)=1.0;
    }
  }
  if (res.empty()){
    cerr << "JHUGenLexiconTranslator::getTranslationMatrix: Translation from input basis " << basis_input << " to output basis " << basis_output << " is not implemented." << endl;
    assert(0);
  }
  return res;
}

std::vector< std::pair<double, double> > JHUGenLexiconTranslator::getOrderedInputCouplings(
  JHUGenLexiconIOHelpers::IOBasisType const& basis_input,
  std::unordered_map<std::string, bool> const& input_flags,
  std::unordered_map<std::string, double> const& input_parameters,
  std::unordered_map<std::string, std::pair<double, double> > const& input_couplings
) const{
  bool useMCFMAtInput; getValueWithDefault<std::string, bool>(input_flags, "useMCFMAtInput", useMCFMAtInput, false);
  bool includeTripleGauge; getValueWithDefault<std::string, bool>(input_flags, "includeTripleGauge", includeTripleGauge, false);
  bool distinguish_HWWcouplings; getValueWithDefault<std::string, bool>(input_flags, "distinguish_HWWcouplings", distinguish_HWWcouplings, false);
  bool switch_convention; getValueWithDefault<std::string, bool>(input_flags, "switch_convention", switch_convention, false);
  double MZ; getValueWithDefault<std::string, double>(input_parameters, "MZ", MZ, DEFVAL_MZ);
  double MW; getValueWithDefault<std::string, double>(input_parameters, "MW", MW, DEFVAL_MW);
  double sw; getValueWithDefault<std::string, double>(input_parameters, "sw", sw, DEFVAL_SW);
  double cw=1.0-sw;
  double Lambda_z1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_z1", Lambda_z1, DEFVAL_LAMBDA_VI);
  double Lambda_w1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_w1", Lambda_w1, DEFVAL_LAMBDA_VI);
  double Lambda_zgs1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_zgs1", Lambda_zgs1, DEFVAL_LAMBDA_VI);

  std::vector< std::pair<double, double> > res;
// The last line dFour_Z checks and scales the output to what is expected by jhugen
#define COUPLING_COMMAND(NAME, PREFIX, DEFVAL) \
  getValueWithDefault<std::string, std::pair<double, double>>(input_couplings, #NAME, res.at(coupl_##PREFIX##_##NAME), std::pair<double, double>(DEFVAL, 0)); \
  if (useMCFMAtInput && (std::string(#NAME).find("ghz")!=std::string::npos || std::string(#NAME).find("ghw")!=std::string::npos)){ res.at(coupl_##PREFIX##_##NAME).first *= 2.; res.at(coupl_##PREFIX##_##NAME).second *= 2.; } \
  if (std::string(#NAME).find("ghzgs")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ res.at(coupl_##PREFIX##_##NAME).first *= std::pow(MZ/Lambda_zgs1, 2); res.at(coupl_##PREFIX##_##NAME).second *= std::pow(MZ/Lambda_zgs1, 2); } \
  else if (std::string(#NAME).find("ghz")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ res.at(coupl_##PREFIX##_##NAME).first *= std::pow(MZ/Lambda_z1, 2); res.at(coupl_##PREFIX##_##NAME).second *= std::pow(MZ/Lambda_z1, 2); } \
  else if (std::string(#NAME).find("ghw")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ res.at(coupl_##PREFIX##_##NAME).first *= std::pow(MW/Lambda_w1, 2); res.at(coupl_##PREFIX##_##NAME).second *= std::pow(MW/Lambda_w1, 2); } \
  if (switch_convention && ((std::string(#NAME).find("ghzgs2")!=std::string::npos) || (std::string(#NAME).find("ghzgs4")!=std::string::npos) || (std::string(#NAME).find("ghzgs1_prime2")!=std::string::npos) || (std::string(#NAME).find("Cza")!=std::string::npos) || (std::string(#NAME).find("tCza")!=std::string::npos) || (std::string(#NAME).find("Cabx")!=std::string::npos))){ res.at(coupl_##PREFIX##_##NAME).first *= -1.0; res.at(coupl_##PREFIX##_##NAME).second *= -1.0; }
  switch (basis_input){
  case bAmplitude_JHUGen:
  {
    res.assign(nAmplitude_JHUGen_CouplingTypes, std::pair<double, double>(0, 0));
    AMPLITUDE_JHUGEN_COUPLING_COMMANDS;
    break;
  }
  case bEFT_JHUGen:
  {
    res.assign(nEFT_JHUGen_CouplingTypes, std::pair<double, double>(0, 0));
    EFT_JHUGEN_COUPLING_COMMANDS;
    break;
  }
  case bHiggsBasis:
  {
    res.assign(nHiggsBasis_CouplingTypes, std::pair<double, double>(0, 0));
    HIGGSBASIS_COUPLING_COMMANDS;
    break;
  }
  case bEFT_HiggsBasis:
  {
    res.assign(nEFT_HiggsBasis_CouplingTypes, std::pair<double, double>(0, 0));
    EFT_HIGGSBASIS_COUPLING_COMMANDS;
    break;
  }
  case bWarsawBasis:
  {
    res.assign(nWarsawBasis_CouplingTypes, std::pair<double, double>(0, 0));
    WARSAWBASIS_COUPLING_COMMANDS;
    break;
  }
  default:
    cerr << "JHUGenLexiconTranslator::getOrderedInputCouplings: Input basis " << basis_input << " is not implemented." << endl;
    assert(0);
  }

#undef COUPLING_COMMAND

  return res;
}

void JHUGenLexiconTranslator::interpretOutputCouplings(
  JHUGenLexiconIOHelpers::IOBasisType const& basis_output,
  std::unordered_map<std::string, bool> const& input_flags,
  std::unordered_map<std::string, double> const& input_parameters,
  std::vector< std::pair<double, double> >& output_vector
){
  auto const& basis_input = opts.getInputBasis();
  bool useMCFMAtOutput; getValueWithDefault<std::string, bool>(input_flags, "useMCFMAtOutput", useMCFMAtOutput, false);
  bool include_triple_quartic_gauge; getValueWithDefault<std::string, bool>(input_flags, "include_triple_quartic_gauge", include_triple_quartic_gauge, false);
  bool HW_couplings_only; getValueWithDefault<std::string, bool>(input_flags, "HW_couplings_only", HW_couplings_only, false);
  bool HZ_couplings_only; getValueWithDefault<std::string, bool>(input_flags, "HZ_couplings_only", HZ_couplings_only, false);
  bool switch_convention; getValueWithDefault<std::string, bool>(input_flags, "switch_convention", switch_convention, false);
  double MZ; getValueWithDefault<std::string, double>(input_parameters, "MZ", MZ, DEFVAL_MZ);
  double MW; getValueWithDefault<std::string, double>(input_parameters, "MW", MW, DEFVAL_MW);
  double Lambda_z1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_z1", Lambda_z1, DEFVAL_LAMBDA_VI);
  double Lambda_w1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_w1", Lambda_w1, DEFVAL_LAMBDA_VI);
  double Lambda_zgs1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_zgs1", Lambda_zgs1, DEFVAL_LAMBDA_VI);
#define COUPLING_COMMAND(NAME, PREFIX, DEFVAL) \
  if (useMCFMAtOutput && (std::string(#NAME).find("ghz")!=std::string::npos || std::string(#NAME).find("ghw")!=std::string::npos)){ output_vector.at(coupl_##PREFIX##_##NAME).first /= 2.; output_vector.at(coupl_##PREFIX##_##NAME).second /= 2.; } \
  if (std::string(#NAME).find("ghzgs")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ output_vector.at(coupl_##PREFIX##_##NAME).first /= std::pow(MZ/Lambda_zgs1, 2); output_vector.at(coupl_##PREFIX##_##NAME).second /= std::pow(MZ/Lambda_zgs1, 2); } \
  else if (std::string(#NAME).find("ghz")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ output_vector.at(coupl_##PREFIX##_##NAME).first /= std::pow(MZ/Lambda_z1, 2); output_vector.at(coupl_##PREFIX##_##NAME).second /= std::pow(MZ/Lambda_z1, 2); } \
  else if (std::string(#NAME).find("ghw")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ output_vector.at(coupl_##PREFIX##_##NAME).first /= std::pow(MW/Lambda_w1, 2); output_vector.at(coupl_##PREFIX##_##NAME).second /= std::pow(MW/Lambda_w1, 2); } \
  else if (basis_input == bWarsawBasis && switch_convention && ((std::string(#NAME).find("ghzgs2")!=std::string::npos) || (std::string(#NAME).find("ghzgs4")!=std::string::npos) || (std::string(#NAME).find("ghzgs1_prime2")!=std::string::npos) || (std::string(#NAME).find("Cza")!=std::string::npos) || (std::string(#NAME).find("tCza")!=std::string::npos) || (std::string(#NAME).find("Cabx")!=std::string::npos))){ output_vector.at(coupl_##PREFIX##_##NAME).first *= -1.0; output_vector.at(coupl_##PREFIX##_##NAME).second *= -1.0; } \
  result_couplings[#NAME] = output_vector.at(coupl_##PREFIX##_##NAME);  
  switch (basis_output){
    case bAmplitude_JHUGen:
    {
      if (include_triple_quartic_gauge){
	if (HW_couplings_only){
          AMPLITUDE_JHUGEN_INCLUDE_TRIPLE_CHARGED_CURRENT_COUPLING_COMMANDS
        }
	else if (HZ_couplings_only){
          AMPLITUDE_JHUGEN_INCLUDE_TRIPLE_NEUTRAL_CURRENT_COUPLING_COMMANDS
        }
        else{
          AMPLITUDE_JHUGEN_INCLUDE_TRIPLE_COUPLING_COMMANDS;
        }
      }
      else{
        if (HW_couplings_only){
          AMPLITUDE_JHUGEN_CHARGED_CURRENT_COUPLING_COMMANDS
        }
        else if (HZ_couplings_only){
          AMPLITUDE_JHUGEN_NEUTRAL_CURRENT_COUPLING_COMMANDS
        }
        else{
          AMPLITUDE_JHUGEN_COUPLING_COMMANDS;
        }
      }
      break;
    }
    case bEFT_JHUGen:
    {
      if (include_triple_quartic_gauge){
        EFT_JHUGEN_INCLUDE_TRIPLE_COUPLING_COMMANDS;
      }
      else{
        EFT_JHUGEN_COUPLING_COMMANDS;
      }
      break;
    }
    case bHiggsBasis:
    {
      if (include_triple_quartic_gauge){
        HIGGSBASIS_INCLUDE_TRIPLE_COUPLING_COMMANDS;
      }
      else{
        HIGGSBASIS_COUPLING_COMMANDS;
      }
      break;
    }
    case bEFT_HiggsBasis:
    {
      if  (include_triple_quartic_gauge){
        EFT_HIGGSBASIS_INCLUDE_TRIPLE_COUPLING_COMMANDS;
      }
      else{
        EFT_HIGGSBASIS_COUPLING_COMMANDS;
      }
      break;
    }
    case bWarsawBasis:
    {
      WARSAWBASIS_COUPLING_COMMANDS;
      break;
    }
    default:
      cerr << "JHUGenLexiconTranslator::interpretOutputCouplings: Output basis " << basis_output << " is not implemented." << endl;
      assert(0);
  }

#undef COUPLING_COMMAND
}
