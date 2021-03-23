#ifndef __pctMostLikelyPathFunction_h
#define __pctMostLikelyPathFunction_h

#include <itkNumericTraits.h>
#include <vector>
#include <itkImageBase.h>

//#define MLP_TIMING
#ifdef MLP_TIMING
#  include <itkTimeProbe.h>
#endif

namespace pct
{

/** \class MostLikelyPathFunction
 * \brief Base class for computing the most likely path of a proton.
 *
 * \ingroup Functions
 */
template <class TCoordRep = double>
class ITK_EXPORT MostLikelyPathFunction :
    public itk::LightObject
{
public:
  /** Standard class typedefs. */
  typedef MostLikelyPathFunction         Self;
  typedef itk::LightObject               Superclass;
  typedef itk::SmartPointer<Self>        Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Useful defines. */
  typedef itk::Vector<TCoordRep, 3> VectorType;

  /** Init the mlp parameters from the input and output directions and positions. */
  virtual void Init(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut)
  {
    itkGenericExceptionMacro("This version of the Init method not implemented for derived class.");
  }

  /** Init the mlp parameters from the input and output directions and positions, and energies. */
  virtual void Init(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut, double eIn, double eOut)
  {
    itkGenericExceptionMacro("This version of the Init method not implemented for derived class.");
  }

  /** Init with additional parameters to consider tracker uncertainties */
  virtual void InitUncertain(const VectorType posIn, const VectorType posOut, const VectorType dirIn, const VectorType dirOut, double dEntry, double dExit, double m_TrackerResolution, double m_TrackerPairSpacing, double m_MaterialBudget)
  {
    itkGenericExceptionMacro("Not implemented in the derived class.");
  }

  /** Evaluate the coordinates (x,y) at depth z. */
  virtual void Evaluate( const TCoordRep z, TCoordRep &x, TCoordRep&y, TCoordRep &dx, TCoordRep&dy )
  {
    itkGenericExceptionMacro("Not implemented in the derived class.");
  }

  /** Vectorised version of the above method. Implement dummy in derived class if not applicable for the type of MLP */
  // NK: maybe explicit <double> should be replaced with <TCoordRep>
  virtual void Evaluate( std::vector<double> u, std::vector<double> &x, std::vector<double> &y )
  {
    itkGenericExceptionMacro("Not implemented in the derived class.");
  }
  
  bool m_CanBeVectorised = false;

#ifdef MLP_TIMING
  /** Print timing information */
  virtual void PrintTiming(std::ostream& os){}
#endif

protected:

  /// Constructor
  MostLikelyPathFunction(){}

  /// Destructor
  ~MostLikelyPathFunction(){}

private:
  MostLikelyPathFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented
};

} // namespace pct

#endif
