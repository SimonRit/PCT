#include "pctschulte_ggo.h"

#include <rtkMacro.h>
#include <rtkGgoFunctions.h>

#include "itkSchulteMLPFunction.h"
#include "pctBetheBlochFunctor.h"

int main(int argc, char * argv[])
{
  GGO(pctschulte, args_info);

  // For position and angle straggling, compute sigma1 in Schulte 2008
  double u = args_info.length_arg;

  const double intForSigmaSqTheta1  = itk::Functor::SchulteMLP::IntegralForSigmaSqTheta ::GetValue(u);
  const double intForSigmaSqTTheta1 = itk::Functor::SchulteMLP::IntegralForSigmaSqTTheta::GetValue(u);
  const double intForSigmaSqT1      = itk::Functor::SchulteMLP::IntegralForSigmaSqT     ::GetValue(u);

  itk::Matrix<double, 2, 2> Sigma1;
  Sigma1(1,1) = intForSigmaSqTheta1;
  Sigma1(0,1) = u * Sigma1(1,1) - intForSigmaSqTTheta1;
  Sigma1(1,0) = Sigma1(0,1);
  Sigma1(0,0) = u * ( 2*Sigma1(0,1) - u*Sigma1(1,1) ) + intForSigmaSqT1;
  Sigma1 *= itk::Functor::SchulteMLP::ConstantPartOfIntegrals::GetValue(0.,u);

  // For mean energy
  pct::Functor::IntegratedBetheBlochProtonStoppingPowerInverse<float, double> bethe;

  switch(args_info.parameter_arg)
    {
    case(parameter_arg_energyMean):
    if(args_info.energy_given)
      std::cout << bethe.GetEnergy(args_info.length_arg*CLHEP::mm, args_info.energy_arg*CLHEP::MeV)/CLHEP::MeV << std::endl;
    else
      std::cout << bethe.GetEnergy(args_info.length_arg*CLHEP::mm)/CLHEP::MeV << std::endl;
    return EXIT_SUCCESS;

    case(parameter_arg_energySD):
    //TODO

    case(parameter_arg_positionSD):
    std::cout << sqrt(Sigma1(0,0)) << std::endl;
    return EXIT_SUCCESS;
    }

//  std::cout << u/CLHEP::mm << '\t' << sqrt(Sigma1(0,0))/CLHEP::mm << std::endl;

  std::cerr << "Unknow parameter required" << std::endl;
  return EXIT_FAILURE;
}
