#ifndef __pctSchulteMLPFunction_h
#define __pctSchulteMLPFunction_h

#include "CLHEP/Units/PhysicalConstants.h"

#include "pctMostLikelyPathFunction.h"

namespace pct
{

namespace Functor
{

namespace SchulteMLP
{

static const double aunit = 1./(CLHEP::MeV*CLHEP::MeV);
#define RIT_COEFF
#if defined(WILLIAMS_COEFF) // Were apparently .01 off
  static const double a0 =  7.507e-4 * aunit;
  static const double a1 =  3.320e-5 * aunit / (CLHEP::cm);
  static const double a2 = -4.171e-7 * aunit / (CLHEP::cm2);
  static const double a3 =  4.488e-7 * aunit / (CLHEP::cm3);
  static const double a4 = -3.739e-8 * aunit / (CLHEP::cm3 * CLHEP::cm);
  static const double a5 =  1.455e-9 * aunit / (CLHEP::cm3 * CLHEP::cm2);
#elif defined(SCHULTE_COEFF)
  static const double a0 =  7.457e-06 * aunit;
  static const double a1 =  4.548e-07 * aunit / (CLHEP::cm);
  static const double a2 = -5.777e-08 * aunit / (CLHEP::cm2);
  static const double a3 =  1.301e-08 * aunit / (CLHEP::cm3);
  static const double a4 = -9.228e-10 * aunit / (CLHEP::cm3 * CLHEP::cm);
  static const double a5 =  2.687e-11 * aunit / (CLHEP::cm3 * CLHEP::cm2);
#elif defined(RIT_COEFF)
  static const double a0 =  7.444724e-06 * aunit;
  static const double a1 =  5.463937e-07 * aunit / (CLHEP::cm);
  static const double a2 = -9.986645e-08 * aunit / (CLHEP::cm2);
  static const double a3 =  2.026409e-08 * aunit / (CLHEP::cm3);
  static const double a4 = -1.420501e-09 * aunit / (CLHEP::cm3 * CLHEP::cm);
  static const double a5 =  3.899100e-11 * aunit / (CLHEP::cm3 * CLHEP::cm2);
#elif defined(NICO_COEFF)
  static const double a0 =  3.599e-06 * aunit;
  static const double a1 =  7.303e-08 * aunit / (CLHEP::cm);
  static const double a2 =  2.419e-09 * aunit / (CLHEP::cm2);
  static const double a3 = -4.216e-11 * aunit / (CLHEP::cm3);
  static const double a4 =  1.601e-12 * aunit / (CLHEP::cm3 * CLHEP::cm);
  static const double a5 =  1.467e-13 * aunit / (CLHEP::cm3 * CLHEP::cm2);
#elif defined(KRAH_COEFF_180MEV)
  static const double a0 =  6.576282e-06 * aunit;
  static const double a1 =  5.202820e-06 * aunit / (CLHEP::cm);
  static const double a2 = -2.017266e-06 * aunit / (CLHEP::cm2);
  static const double a3 =  3.289456e-07 * aunit / (CLHEP::cm3);
  static const double a4 = -2.215321e-08 * aunit / (CLHEP::cm3 * CLHEP::cm);
  static const double a5 =  5.398933e-10 * aunit / (CLHEP::cm3 * CLHEP::cm2);
#endif

// [Schulte, Med Phys, 2008], constant part of equations 7, 8 and 9
class ConstantPartOfIntegrals
{
public:
  static double GetValue(const double ux, const double uy)
    {
    // Constant out of integral, use [Schulte, 2008] and not [Williams, 2004]
    static const double c = 13.6*CLHEP::MeV * 13.6*CLHEP::MeV / (36.1*CLHEP::cm);
    static const double invX0 = 1./(36.1*CLHEP::cm);
    double diffU = uy-ux;
    if(diffU<1e-3) diffU = 1e-3 * CLHEP::mm; // avoid divergence of log for small values
    double v = 1.+0.038*std::log((diffU)*invX0);
    return c*v*v;
    }
};

// [Schulte, Med Phys, 2008], integral for sigma in equation 8
class IntegralForSigmaSqTheta
{
public:
  static double GetValue(const double u)
    {
    // Multiplied with polynomial integral
    // Table 2, (a) in [Williams, 2004] as well as numbers in [Li, 2006]
    static const double a0Theta = a0;
    static const double a1Theta = a1 / 2.;
    static const double a2Theta = a2 / 3.;
    static const double a3Theta = a3 / 4.;
    static const double a4Theta = a4 / 5.;
    static const double a5Theta = a5 / 6.;
    return u*(a0Theta+u*(a1Theta+u*(a2Theta+u*(a3Theta+u*(a4Theta+u*a5Theta)))));
    }
};

// [Schulte, Med Phys, 2008], integral for sigma in equation 9
class IntegralForSigmaSqTTheta
{
public:
  static double GetValue(const double u)
    {
    // Multiplied with polynomial integral
    // Table 2, (a) in [Williams, 2004] as well as numbers in [Li, 2006]
    static const double a0TTheta = a0 / 2.;
    static const double a1TTheta = a1 / 3.;
    static const double a2TTheta = a2 / 4.;
    static const double a3TTheta = a3 / 5.;
    static const double a4TTheta = a4 / 6.;
    static const double a5TTheta = a5 / 7.;
    return u*u*(a0TTheta+u*(a1TTheta+u*(a2TTheta+u*(a3TTheta+u*(a4TTheta+u*a5TTheta)))));
    }
};

// [Schulte, Med Phys, 2008], integral for sigma in equation 7
class IntegralForSigmaSqT
{
public:
  static double GetValue(const double u)
    {
    // Multiplied with polynomial integral
    // Table 2, (a) in [Williams, 2004] as well as numbers in [Li, 2006]
    static const double a0T = a0 / 3.;
    static const double a1T = a1 / 4.;
    static const double a2T = a2 / 5.;
    static const double a3T = a3 / 6.;
    static const double a4T = a4 / 7.;
    static const double a5T = a5 / 8.;
    return u*u*u*(a0T+u*(a1T+u*(a2T+u*(a3T+u*(a4T+u*a5T)))));
    }
};

} // end namespace SchulteMLP

} // end namespace Functor

/** \class SchulteMLPFunction
 * \brief See [Schulte, Med Phys, 2008].
 *
 * \ingroup Functions
 */
class ITK_EXPORT SchulteMLPFunction:
    public MostLikelyPathFunction<double>
{
public:
  /** Standard class typedefs. */
  typedef SchulteMLPFunction                       Self;
  typedef MostLikelyPathFunction<double>           Superclass;
  typedef itk::SmartPointer<Self>                  Pointer;
  typedef itk::SmartPointer<const Self>            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Useful defines. */
  typedef Superclass::VectorType VectorType;

  /** Init the mlp parameters from the input and output directions and positions. */
  virtual void Init(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut) override;

  /** Init with additional parameters to consider tracker uncertainties */
  virtual void InitUncertain(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut, double dEntry, double dExit, double TrackerResolution, double TrackerPairSpacing, double MaterialBudget) override;

  /** Evaluate the coordinates (x,y) at depth u1. */
  virtual void Evaluate( const double u1, double &x, double &y, double &dx, double &dy ) override;

  /** Evaluate the error (x,y) (equation 27) at depth z. */
  void EvaluateError( const double u1, itk::Matrix<double, 2, 2> &error);

#ifdef MLP_TIMING
  /** Print timing information */
  virtual void PrintTiming(std::ostream& os) override;
#endif

protected:

  /// Implementation of 2x2 matrix inversion, faster than itk/vnl inversion
  void InverseMatrix(itk::Matrix<double, 2, 2> &mat);

  /// Constructor
  SchulteMLPFunction();

  /// Destructor
  ~SchulteMLPFunction(){}

private:
  SchulteMLPFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  // Depth position at entrance and exit, only u1 is variable
  double m_uOrigin;
  double m_u0;
  double m_u2;

  // Entrance and exit parameters (equation 1)
  itk::Vector<double, 2> m_x0;
  itk::Vector<double, 2> m_x2;
  itk::Vector<double, 2> m_y0;
  itk::Vector<double, 2> m_y2;

  // Part of the rotation matrices which is constant for the trajectory
  itk::Matrix<double, 2, 2> m_R0;
  itk::Matrix<double, 2, 2> m_R1;
  itk::Matrix<double, 2, 2> m_R0T;
  itk::Matrix<double, 2, 2> m_R1T;

  // Scattering matrices
  itk::Matrix<double, 2, 2> m_Sigma1;
  itk::Matrix<double, 2, 2> m_Sigma2;

  bool m_considerTrackerUncertainties;

  // Addtional matrices needed for MLP with tracker uncertainties
  itk::Matrix<double, 2, 2> m_SigmaIn;
  itk::Matrix<double, 2, 2> m_SigmaOut;
  itk::Matrix<double, 2, 2> m_Sin;
  itk::Matrix<double, 2, 2> m_SinT;
  itk::Matrix<double, 2, 2> m_Sout;
  itk::Matrix<double, 2, 2> m_SoutT;
  itk::Matrix<double, 2, 2> m_Sout_Inv;
  itk::Matrix<double, 2, 2> m_SoutT_Inv;
  // Part common to all positions along the trajectory
  //double m_IntForSigmaSqTheta0;  //Always 0. because m_u0=0.
  //double m_IntForSigmaSqTTheta0; //Always 0. because m_u0=0.
  //double m_IntForSigmaSqT0;      //Always 0. because m_u0=0.
  double m_IntForSigmaSqTheta2;
  double m_IntForSigmaSqTTheta2;
  double m_IntForSigmaSqT2;

#ifdef MLP_TIMING
  itk::TimeProbe m_EvaluateProbe1;
  itk::TimeProbe m_EvaluateProbe2;
#endif

};

} // end namespace pct

#include "pctSchulteMLPFunction.txx"

#endif
