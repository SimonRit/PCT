#include "pctschulte_ggo.h"

#include <rtkMacro.h>
#include <rtkGgoFunctions.h>

#include "pctSchulteMLPFunction.h"
#include "pctBetheBlochFunctor.h"
#include "pctEnergyStragglingFunctor.h"

#include <fstream>

int main(int argc, char * argv[])
{
  GGO(pctschulte, args_info);

  // For position and angle straggling, compute sigma1 in Schulte 2008
  double u = args_info.length_arg;

  const double intForSigmaSqTheta1  = pct::Functor::SchulteMLP::IntegralForSigmaSqTheta ::GetValue(u);
  const double intForSigmaSqTTheta1 = pct::Functor::SchulteMLP::IntegralForSigmaSqTTheta::GetValue(u);
  const double intForSigmaSqT1      = pct::Functor::SchulteMLP::IntegralForSigmaSqT     ::GetValue(u);

  itk::Matrix<double, 2, 2> Sigma1;
  Sigma1(1,1) = intForSigmaSqTheta1;
  Sigma1(0,1) = u * Sigma1(1,1) - intForSigmaSqTTheta1;
  Sigma1(1,0) = Sigma1(0,1);
  Sigma1(0,0) = u * ( 2*Sigma1(0,1) - u*Sigma1(1,1) ) + intForSigmaSqT1;
  Sigma1 *= pct::Functor::SchulteMLP::ConstantPartOfIntegrals::GetValue(0.,u);

  // For mean energy
  pct::Functor::IntegratedBetheBlochProtonStoppingPowerInverse<float, double> bethe( args_info.ionpot_arg * CLHEP::eV, 500 * CLHEP::MeV );

  if(args_info.dEdx_given)
    {
    std::ofstream os(args_info.dEdx_arg);
    pct::Functor::BetheBlochProtonStoppingPower<double, double> bb;
    os << "E (MeV), dEdX (MeV.cm²/g)" << std::endl;
    for(double E=250; E>=5.; E-=5.)
      {
      os << E/CLHEP::MeV;
      os << ',';
      os << bb.GetValue(E*CLHEP::MeV, args_info.ionpot_arg*CLHEP::eV)/(1000.*CLHEP::kg/CLHEP::m3) * CLHEP::g/(CLHEP::MeV*CLHEP::cm2);
      os << std::endl;
      }
    os.close();
    }

  pct::SchulteMLPFunction::Pointer mlp;
  mlp = pct::SchulteMLPFunction::New();
  itk::Vector<double,3> p1,p2,d;
  p1.Fill(0.);
  p2 = p1;
  d = p1;
  p1[2] = -110.;
  p2[2] = 110.;
  d[2] = 1.;
  mlp->Init(p1, p2, d, d);
  itk::Matrix<double, 2, 2> error;

  switch(args_info.parameter_arg)
    {
    case(parameter_arg_energyMean):
    if(args_info.energy_given)
      std::cout << bethe.GetEnergy(args_info.length_arg*CLHEP::mm, args_info.energy_arg*CLHEP::MeV)/CLHEP::MeV << std::endl;
    else
      std::cout << bethe.GetEnergy(args_info.length_arg*CLHEP::mm)/CLHEP::MeV << std::endl;
    return EXIT_SUCCESS;

    case(parameter_arg_energySD):
    std::cout << pct::Functor::EnergyStragglingFunctor<float, double>::GetValue(args_info.length_arg*CLHEP::mm)/CLHEP::MeV << std::endl;
    return EXIT_SUCCESS;

    case(parameter_arg_positionSD):
    std::cout << sqrt(Sigma1(0,0)) << std::endl;
    return EXIT_SUCCESS;

    case(parameter_arg_mlpSD):
    mlp->EvaluateError(args_info.length_arg, error);
    std::cout << sqrt(error[0][0]) << std::endl;
    return EXIT_SUCCESS;

    default:
    return EXIT_FAILURE;
    }

//  std::cout << u/CLHEP::mm << '\t' << sqrt(Sigma1(0,0))/CLHEP::mm << std::endl;

  std::cerr << "Unknow parameter required" << std::endl;
  return EXIT_FAILURE;
}
