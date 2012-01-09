#ifndef __itkThirdOrderPolynomialMLPFunction_h
#define __itkThirdOrderPolynomialMLPFunction_h

#include "itkMostLikelyPathFunction.h"

namespace itk
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
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Useful defines. */
  typedef typename Superclass::VectorType VectorType;

  /** Init the mlp parameters from the input and output directions and positions. */
  void Init(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut);

  /** Evaluate the coordinates (x,y) at depth z. */
  void Evaluate( const TCoordRep z, TCoordRep &x, TCoordRep&y );

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

} // namespace itk

#include "itkThirdOrderPolynomialMLPFunction.txx"

#endif
