#ifndef __pctScatteringWEPLFunction_h
#define __pctScatteringWEPLFunction_h

#include "pctBetheBlochFunctor.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "TMath.h"

namespace pct
{

namespace Functor
{

namespace ScatteringWEPL
{
std::vector<double> sigma2_lut;
std::vector<int> thickness_lut;
unsigned int range = 0;
const double proton_mass_c2 = CLHEP::proton_mass_c2;

// Equation 8 of Schulte Med Phys 2008
class Highland
{
public:
  static double GetSigma2(const double thickness, const double invbeta2p2)
    {
    const double E0 = 13.6*CLHEP::MeV;
    const double X0 = 38.08*CLHEP::cm;
    return E0*E0/X0 * (1+0.038*vcl_log(thickness/X0))* (1+0.038*vcl_log(thickness/X0)) * invbeta2p2;
    }
};

// Equations 39 and 40 of Gottschalk Med Phys 2010
class DifferentialMoliere
{
public:
  static double GetSigma2(const double eIn, const double invbeta2p2)
    {
    const double p1v1 = 1./((eIn*CLHEP::MeV + proton_mass_c2) / (eIn*CLHEP::MeV + 2*proton_mass_c2) / eIn*CLHEP::MeV);
    const double Xs = 46.88*CLHEP::cm; 
    const double pv = TMath::Sqrt(1./invbeta2p2);
    double fdm = (0.5244 + 0.1975 * TMath::Log10(1-(pv/p1v1)*(pv/p1v1)) + 0.2320 * TMath::Log10(pv) - 0.0098 * TMath::Log10(pv) * TMath::Log10(1-(pv/p1v1)*(pv/p1v1)));
    return fdm * 15.* 15. / pv / pv / Xs;
    }
};

class ScatteringLUT
{
public:
  // LUT from analytical scattering models   
  static double SetLUT(const double eIn)
    {
    double m_IonizationPotential = 68.9984*CLHEP::eV; 
    pct::Functor::IntegratedBetheBlochProtonStoppingPowerInverse<float, double> bethe( m_IonizationPotential, 500 * CLHEP::MeV );
    unsigned const int max_waterthickness = 500*CLHEP::mm;
    double invbeta2p2  = 0.;
    double energycurrent  = 0.;
    double sigma2 = 0.;

    for(unsigned int waterthickness = 1; waterthickness<max_waterthickness; waterthickness++)
      {    
      energycurrent = bethe.GetEnergy(waterthickness*CLHEP::mm, eIn*CLHEP::MeV);
      if(energycurrent<0.001) break;        
      invbeta2p2 = (energycurrent+proton_mass_c2)*(energycurrent+proton_mass_c2) / 
                  ((energycurrent+2*proton_mass_c2)*(energycurrent+2*proton_mass_c2) * energycurrent*energycurrent); 
      //sigma2 += Highland::GetSigma2(waterthickness,invbeta2p2); 
      sigma2 += DifferentialMoliere::GetSigma2(eIn,invbeta2p2); // Closer to the G4 values than Highland's 
      //std::cout<< sigma2 <<std::endl;   
      sigma2_lut.push_back( sigma2 );
      thickness_lut.push_back(waterthickness);
      range = waterthickness; 
      }
    //return sigma2_lut, thickness_lut ; 
    }

  static double SetG4LUT(const double eIn)
    {
    // Reads LUT from a txt file
    char inputfilename[200];
    sprintf(inputfilename, "/home/ctquinones/src/pct/database/ScatteringWaterLut_%iMeV.txt", (int)eIn);
    std::ifstream inputfile;
    inputfile.open(inputfilename);

    double waterthickness = 0.;
    double sigma2 = 0.;    

    if(!inputfile) 
      {
      std::cout<<"File not found"<<std::endl;
      }

    while (!inputfile.eof()) 
      {
      inputfile >> waterthickness >> sigma2;
      // Stores the values 
      sigma2_lut.push_back( sigma2 );
      thickness_lut.push_back(waterthickness);
      }
      range = waterthickness; 
    }// end SetG4LUT

}; // end class Scattering LUT

class ConvertToScatteringWEPL
{
public:
  static double GetValue(const double sigma2)
    {
    double value = 0., x0 = 0., y0 = 0., x1 = 0., y1 = 0.;
    for(unsigned int i = 0; i < range; i++)
      {
      if(sigma2 < sigma2_lut[i])
      continue;
        {
        x0 = sigma2_lut[i];
        x1 = sigma2_lut[i+1];
        y0 = thickness_lut[i];
        y1 = thickness_lut[i+1];  
        }
       } 

    // calculate wepl using linear interpolation      
    value = y0 + (y1-y0) * (sigma2-x0) / (x1 - x0); 
    if(value!=value)
      value = 0.;
  
    return value;
    }

}; // end class ConvertToScatteringWEPL

}; // end namespace ScatteringWEPL

}; // end namespace Functor

} // end namespace pct

#endif
