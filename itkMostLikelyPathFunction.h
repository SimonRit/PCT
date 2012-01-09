#ifndef __itkMostLikelyPathFunction_h
#define __itkMostLikelyPathFunction_h

#include <itkNumericTraits.h>
#include <vector>
#include <itkImageBase.h>

namespace itk
{

/** \class MostLikelyPathFunction
 * \brief Base class for computing the most likely path of a proton.
 *
 * \ingroup Functions
 */
template <class TCoordRep = double>
class ITK_EXPORT MostLikelyPathFunction :
    public LightObject
{
public:
  /** Standard class typedefs. */
  typedef MostLikelyPathFunction                                Self;
  typedef LightObject                                           Superclass;
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;

  /** Useful defines. */
  typedef Vector<TCoordRep, 3> VectorType;

  /** Init the mlp parameters from the input and output directions and positions. */
  virtual void Init(VectorType posIn, VectorType posOut, VectorType dirIn, VectorType dirOut){}

  /** Evaluate the coordinates (x,y) at depth z. */
  virtual void Evaluate( const TCoordRep z, TCoordRep &x, TCoordRep&y ){}

protected:

  /// Constructor
  MostLikelyPathFunction(){}

  /// Destructor
  ~MostLikelyPathFunction(){}

private:
  MostLikelyPathFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented
};

} // namespace itk

#endif
