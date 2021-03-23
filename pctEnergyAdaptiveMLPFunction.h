#ifndef __pctEnergyAdaptiveMLPFunction_h
#define __pctEnergyAdaptiveMLPFunction_h

#include "CLHEP/Units/PhysicalConstants.h"

#include "pctMostLikelyPathFunction.h"

// #include <cmath>

namespace pct
{

namespace Functor
{

namespace FlexibleMLP
{

void
GetLinearCoefficients( itk::Vector<double, 2>& ab, const double E_in, const double E_out, const double deltaU )
{
  double inverseScatteringPower_in, inverseScatteringPower_out;

  /** Calculate 1/T = beta^2 p^2 / Omega0^2 * X0
  * 938.3 is the proton rest mass; It would be better to avoid hard coding
  * the contants 13.6 and X0=36.1 actually cancel out when evaluating the MLP,
  * but it is better to keep them for physics debugging.
  * inverseScatteringPower is unitless.
  */
  inverseScatteringPower_in = (E_in+2*938.3) * E_in / (E_in+938.3) / 13.6;
  inverseScatteringPower_in *= inverseScatteringPower_in * 36.1 * CLHEP::cm;
  inverseScatteringPower_out = (E_out+2*938.3) * E_out / (E_out+938.3) / 13.6;
  inverseScatteringPower_out *= inverseScatteringPower_out * 36.1*CLHEP::cm;

  ab[1] = inverseScatteringPower_in;
  ab[0] = (inverseScatteringPower_out - inverseScatteringPower_in) / deltaU;
}

// class which calculates factor A, B, C, D as in eq. 19 in [Krah 2019, PMB]
class FactorsABCD
{
public:
  static double GetA(const double uOut, itk::Vector<double, 2> ab)
    {
    double A = std::log(uOut*ab[0]/ab[1] + 1.) / ab[0];
    return A;
    }

  static double GetB(const double uOut, itk::Vector<double, 2> ab)
    {
    double uOut_tilde = uOut*ab[0]/ab[1];
    double B = (uOut_tilde - std::log(uOut_tilde + 1)) * ab[1] / std::pow(ab[0], 2);
    return B;
    }

  static double GetC(const double uOut, itk::Vector<double, 2> ab)
    {
    double uOut_tilde = uOut*ab[0]/ab[1];
    double C = ((uOut_tilde + 1)*std::log(uOut_tilde + 1) - uOut_tilde) * ab[1] / std::pow(ab[0], 2);
    return C;
    }

  static double GetD(const double uOut, itk::Vector<double, 2> ab)
    {
    double uOut_tilde = uOut*ab[0]/ab[1];
    double D = (0.5*uOut_tilde*uOut_tilde + uOut_tilde - (uOut_tilde + 1)*std::log(uOut_tilde + 1)) * std::pow(ab[1], 2) / std::pow(ab[0], 3);
    return D;
    }
};

class CoefficientsC
{
public:
  static void GetValue(itk::Vector<double, 2>& CoefficientsC, const double uOut, const itk::Vector<double, 2> vIn, const itk::Vector<double, 2> vOut, const double A, const double B, const double C, const double D)
  {
    CoefficientsC[0] = (-B * (vOut[0] - vIn[0] - vIn[1] * uOut) + D * (vOut[1] - vIn[1])) / (A*D - B*C);
    CoefficientsC[1] = (A * (vOut[0] - vIn[0] - vIn[1] * uOut) - C * (vOut[1] - vIn[1])) / (A*D - B*C);
  }
};


} // end namespace EnergyAdaptiveMLPFunction

} // end namespace Functor

/** \class EnergyAdaptiveMLPFunction
 *
 * \ingroup Functions
 */
class ITK_EXPORT EnergyAdaptiveMLPFunction:
    public MostLikelyPathFunction<double>
{
public:
  /** Standard class typedefs. */
  typedef EnergyAdaptiveMLPFunction                       Self;
  typedef MostLikelyPathFunction<double>           Superclass;
  typedef itk::SmartPointer<Self>                  Pointer;
  typedef itk::SmartPointer<const Self>            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Useful defines. */
  typedef Superclass::VectorType VectorType;

  /* Initialize the MLP estimator */
  virtual void Init(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut, double eIn, double eOut) override;

  /* Evaluate MLP in multiple depths u */ 
  virtual void Evaluate( std::vector<double> u, std::vector<double> &x, std::vector<double> &y ) override;

  /** Evaluate the error (x,y) (equation 27) at depth z. */
  void EvaluateError( const double u1, itk::Matrix<double, 2, 2> &error);

#ifdef MLP_TIMING
  /** Print timing information */
  virtual void PrintTiming(std::ostream& os) override;
#endif

protected:

  /// Constructor
  EnergyAdaptiveMLPFunction();
  EnergyAdaptiveMLPFunction(const int polydeg);

  /// Destructor
  ~EnergyAdaptiveMLPFunction(){}

private:
  EnergyAdaptiveMLPFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  itk::Vector<double, 1> m_ScalarTest;

  // vectors holding the linear coefficients a and b of the inverse kinematic term
  itk::Vector<double, 2> m_ab;

  // vectors holding the factors c0 and c1; constant per MLP
  itk::Vector<double, 2> m_c_x;
  itk::Vector<double, 2> m_c_y;

  // regrouped coefficients for MLP calculation; constant per MLP
  itk::Vector<double, 4> m_dm_x;
  itk::Vector<double, 4> m_dm_y;

  // Depth position at entrance and exit, only u1 is variable
  double m_uOrigin;
  double m_u0;
  double m_u2;

  // Entrance and exit parameters (equation 1)
  itk::Vector<double, 2> m_x0;
  itk::Vector<double, 2> m_x2;
  itk::Vector<double, 2> m_y0;
  itk::Vector<double, 2> m_y2;


#ifdef MLP_TIMING
  itk::TimeProbe m_EvaluateProbe1;
  itk::TimeProbe m_EvaluateProbe2;
#endif

};

} // end namespace pct

#include "pctEnergyAdaptiveMLPFunction.txx"

#endif
