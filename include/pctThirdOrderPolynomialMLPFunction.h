#ifndef __pctThirdOrderPolynomialMLPFunction_h
#define __pctThirdOrderPolynomialMLPFunction_h

#include "pctMostLikelyPathFunction.h"

namespace pct
{

/** \class ThirdOrderPolynomialMLPFunction
 * \brief Fit a third order polynomial given input and output direction and position.
 *
 * \ingroup Functions
 */
template <class TCoordRep = double>
class ITK_EXPORT ThirdOrderPolynomialMLPFunction :
    public MostLikelyPathFunction<TCoordRep>
{
public:
  /** Standard class typedefs. */
  typedef ThirdOrderPolynomialMLPFunction                       Self;
  typedef MostLikelyPathFunction<TCoordRep>                     Superclass;
  typedef itk::SmartPointer<Self>                               Pointer;
  typedef itk::SmartPointer<const Self>                         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Useful defines. */
  typedef typename Superclass::VectorType VectorType;

  /** Init the mlp parameters from the input and output directions and positions. */
  virtual void Init(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut) override;

  /** Init with additional parameters to consider tracker uncertainties */
  virtual void InitUncertain(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut, double dEntry, double dExit, double m_TrackerResolution, double m_TrackerPairSpacing, double m_MaterialBudget) override;

  /** Evaluate the coordinates (x,y) at depth z. */
  virtual void Evaluate( const TCoordRep z, TCoordRep &x, TCoordRep&y, TCoordRep &dx, TCoordRep&dy ) override;

  // vectorised version:
  virtual void Evaluate( std::vector<double> u, std::vector<double> &x, std::vector<double> &y ) override;

protected:

  /// Constructor
  ThirdOrderPolynomialMLPFunction(){}

  /// Destructor
  ~ThirdOrderPolynomialMLPFunction(){}

private:
  ThirdOrderPolynomialMLPFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  // Polynomial parameters in each direction
  TCoordRep ax, bx, cx, dx;
  TCoordRep ay, by, cy, dy;
  TCoordRep zoffset;
};

} // namespace pct

#include "pctThirdOrderPolynomialMLPFunction.txx"

#endif
