#ifndef __pctPolynomialMLPFunction_h
#define __pctPolynomialMLPFunction_h

#include "CLHEP/Units/PhysicalConstants.h"

#include "pctMostLikelyPathFunction.h"

namespace pct
{

namespace Functor
{

namespace PolynomialMLP
{

static const double aunit = 1./(CLHEP::MeV*CLHEP::MeV);

// *** 200 MeV ***
// fill coefficient vectors for each available polynomial degree
// // polynomial order N = 0
// static const std::vector<double> bm_0 = {9.496308e-06};
// // polynomial order N = 1
// static const std::vector<double> bm_1 = {-1.055104e-07, 8.365792e-08};
// // polynomial order N = 2
// static const std::vector<double> bm_2 = {6.218903e-06, -8.183818e-08, 7.209602e-10};
// // polynomial order N = 3
// static const std::vector<double> bm_3 = {2.522468e-06, 1.119469e-07, -1.390729e-09, 6.132850e-12};
// // polynomial order N = 4
// static const std::vector<double> bm_4 = {4.562500e-06, -6.670635e-08, 2.116152e-09, -1.764070e-11, 5.178304e-14};
// // polynomial order N = 5
// static const std::vector<double> bm_5 = {3.474283e-06, 7.665043e-08, -2.265353e-09, 3.330223e-11, -1.979538e-13, 4.351773e-16};


// *** 180 MeV ***
// polynomial order N = 0
static const std::vector<double> bm_0 = {1.125895e-05};
// polynomial order N = 1
static const std::vector<double> bm_1 = {2.221018e-07, 1.176787e-07};
// polynomial order N = 2
static const std::vector<double> bm_2 = {7.256003e-06, -1.075769e-07, 1.200877e-09};
// polynomial order N = 3
static const std::vector<double> bm_3 = {3.279817e-06, 1.475374e-07, -2.201247e-09, 1.209155e-11};
// polynomial order N = 4
static const std::vector<double> bm_4 = {5.401888e-06, -7.991486e-08, 3.262799e-09, -3.323899e-11, 1.208325e-13};
// polynomial order N = 5
static const std::vector<double> bm_5 = {4.307328e-06, 9.657939e-08, -3.338966e-09, 6.069645e-11, -4.427153e-13, 1.201749e-15};

// ADD COMMENT HERE
class FactorsABCD
{
public:
  static double GetA(const double uOut, const std::vector<double> bm)
    {
    double A = 0;
    for(std::vector<int>::size_type i = 0; i != bm.size(); i++)
    {
      A += bm[i] / (i+1) * std::pow(uOut, i+1);
    }
    return A;
    }

  static double GetB(const double uOut, const std::vector<double> bm)
    {
    double B = 0;
    for(std::vector<int>::size_type i = 0; i != bm.size(); i++)
    {
      B += bm[i] / (i+2) * std::pow(uOut, i+2);
    }
    return B;
    }

  static double GetC(const double uOut, const std::vector<double> bm)
    {
    double C = 0;
    for(std::vector<int>::size_type i = 0; i != bm.size(); i++)
    {
      C += bm[i] / (i+1) / (i+2) * std::pow(uOut, i+2);
    }
    return C;
    }

  static double GetD(const double uOut, const std::vector<double> bm)
    {
    double D = 0;
    for(std::vector<int>::size_type i = 0; i != bm.size(); i++)
    {
      D += bm[i] / (i+2) / (i+3) * std::pow(uOut, i+3);
    }
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

// class CoefficientsD
// {
// public:
//   static void GetValue(std::vector<double>& dm, const itk::Vector<double, 2> vIn, const itk::Vector<double, 2> c, std::vector<double> bm)
//   {
//     // std::vector<int>::size_type bmsize = bm.size();
//     std::cout << "bm.size() = " << bm.size() << std::endl;
//     dm.push_back(vIn[0]);
//     dm.push_back(vIn[1]);
//     dm.push_back(c[0]*bm[0]/2);
//     for(std::vector<int>::size_type i = 0; i != (bm.size()-1); i++)
//     {
//       dm.push_back((c[0]*bm[i+1] + c[1]*bm[i]) / (i+3) / (i+2));
//     }
//     std::vector<int>::size_type M = bm.size() - 1;
//     dm.push_back(c[1] * bm[M] / (M+2) / (M+3));
//     std::cout << "dm.size() = " << dm.size() << std::endl;
//
//   }
// };



} // end namespace PolynomialMLP

} // end namespace Functor

/** \class SchulteMLPFunction
 * \brief See [Schulte, Med Phys, 2008].
 *
 * \ingroup Functions
 */
class ITK_EXPORT PolynomialMLPFunction:
    public MostLikelyPathFunction<double>
{
public:
  /** Standard class typedefs. */
  typedef PolynomialMLPFunction                       Self;
  typedef MostLikelyPathFunction<double>           Superclass;
  typedef itk::SmartPointer<Self>                  Pointer;
  typedef itk::SmartPointer<const Self>            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Useful defines. */
  typedef Superclass::VectorType VectorType;

  /** Init the mlp parameters from the input and output directions and positions. */
  virtual void Init(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut) ITK_OVERRIDE;

  /** Init with additional parameters to consider tracker uncertainties */
  virtual void InitUncertain(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut, double dEntry, double dExit, double m_TrackerResolution, double m_TrackerPairSpacing, double m_MaterialBudget) ITK_OVERRIDE;

  /** Evaluate the coordinates (x,y) at depth z. */
  virtual void Evaluate( const double u1, double &x, double&y, double &dx, double&dy ) ITK_OVERRIDE;

  // vectorised version:
  virtual void Evaluate( std::vector<double> u, std::vector<double> &x, std::vector<double> &y ) ITK_OVERRIDE;

  /** Evaluate the error (x,y) (equation 27) at depth z. */
  void EvaluateError( const double u1, itk::Matrix<double, 2, 2> &error);

  void SetPolynomialDegree( const int polydeg );
  // itkSetMacro(PolynomialDegree, int)
  // itkGetMacro(PolynomialDegree, int)

#ifdef MLP_TIMING
  /** Print timing information */
  virtual void PrintTiming(std::ostream& os) override;
#endif

protected:

  /// Constructor
  PolynomialMLPFunction();
  PolynomialMLPFunction(const int polydeg);

  /// Destructor
  ~PolynomialMLPFunction(){}

private:
  PolynomialMLPFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  // parameters of the polynomial to describe 1/beta^2p^2 term
  int m_PolynomialDegree;
  int m_PolynomialDegreePlusThree;
  std::vector<double> m_bm;

  // std::vector<double> m_dm_x;
  // std::vector<double> m_dm_y;
  itk::Vector<double, 9> m_dm_x;
  itk::Vector<double, 9> m_dm_y;

  itk::Vector<double, 1> m_ScalarTest;

  // vectors holding the constants c0 and c1
  itk::Vector<double, 2> m_c_x;
  itk::Vector<double, 2> m_c_y;

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

#include "pctPolynomialMLPFunction.txx"

#endif
