#ifndef __itkSchulteMLPFunction_h
#define __itkSchulteMLPFunction_h

#include "CLHEP/Units/PhysicalConstants.h"

#include "itkMostLikelyPathFunction.h"

namespace itk
{

namespace Functor
{

namespace SchulteMLP
{

// [Schulte, Med Phys, 2008], constant part of equations 7, 8 and 9
class ConstantPartOfIntegrals
{
public:
  static double GetValue(const double ux, const double uy)
    {
    // Constant out of integral, use [Schulte, 2008] and not [Williams, 2004]
    static const double c = 13.6*CLHEP::MeV * 13.6*CLHEP::MeV / (36.1*CLHEP::cm);
    static const double invX0 = 1./(36.1*CLHEP::cm);
    double v = 1+0.038*log((uy-ux)*invX0);
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
    static const double aunit = 1./(CLHEP::MeV*CLHEP::MeV);
    static const double a0 =  7.507e-4      * aunit;
    static const double a1 =  3.320e-5 / 2. * aunit / (CLHEP::cm);
    static const double a2 = -4.171e-7 / 3. * aunit / (CLHEP::cm2);
    static const double a3 =  4.488e-7 / 4. * aunit / (CLHEP::cm3);
    static const double a4 = -3.739e-8 / 5. * aunit / (CLHEP::cm3 * CLHEP::cm);
    static const double a5 =  1.455e-9 / 6. * aunit / (CLHEP::cm3 * CLHEP::cm2);
    return u*(a0+u*(a1+u*(a2+u*(a3+u*(a4+u*a5)))));
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
    static const double aunit = 1./(CLHEP::MeV*CLHEP::MeV);
    static const double a0 =  7.507e-4 / 2. * aunit;
    static const double a1 =  3.320e-5 / 3. * aunit / (CLHEP::cm);
    static const double a2 = -4.171e-7 / 4. * aunit / (CLHEP::cm2);
    static const double a3 =  4.488e-7 / 5. * aunit / (CLHEP::cm3);
    static const double a4 = -3.739e-8 / 6. * aunit / (CLHEP::cm3 * CLHEP::cm);
    static const double a5 =  1.455e-9 / 7. * aunit / (CLHEP::cm3 * CLHEP::cm2);
    return u*u*(a0+u*(a1+u*(a2+u*(a3+u*(a4+u*a5)))));
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
    static const double aunit = 1./(CLHEP::MeV*CLHEP::MeV);
    static const double a0 =  7.507e-4 / 3. * aunit;
    static const double a1 =  3.320e-5 / 4. * aunit / (CLHEP::cm);
    static const double a2 = -4.171e-7 / 5. * aunit / (CLHEP::cm2);
    static const double a3 =  4.488e-7 / 6. * aunit / (CLHEP::cm3);
    static const double a4 = -3.739e-8 / 7. * aunit / (CLHEP::cm3 * CLHEP::cm);
    static const double a5 =  1.455e-9 / 8. * aunit / (CLHEP::cm3 * CLHEP::cm2);
    return u*u*u*(a0+u*(a1+u*(a2+u*(a3+u*(a4+u*a5)))));
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
  typedef SmartPointer<Self>                       Pointer;
  typedef SmartPointer<const Self>                 ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Useful defines. */
  typedef Superclass::VectorType VectorType;

  /** Init the mlp parameters from the input and output directions and positions. */
  void Init(VectorType posIn, VectorType posOut, VectorType dirIn, VectorType dirOut);

  /** Evaluate the coordinates (x,y) at depth z. */
  void Evaluate( const double u1, double &x, double&y );

protected:

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

  // Part common to all positions along the trajectory
  //double m_IntForSigmaSqTheta0;  //Always 0. because m_u0=0.
  //double m_IntForSigmaSqTTheta0; //Always 0. because m_u0=0.
  //double m_IntForSigmaSqT0;      //Always 0. because m_u0=0.
  double m_IntForSigmaSqTheta2;
  double m_IntForSigmaSqTTheta2;
  double m_IntForSigmaSqT2;
};

} // end namespace itk

#include "itkSchulteMLPFunction.txx"

#endif
