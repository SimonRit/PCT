#ifndef __pctProtonPairsToDistanceDrivenProjection_h
#define __pctProtonPairsToDistanceDrivenProjection_h

#include "rtkConfiguration.h"
#include "pctBetheBlochFunctor.h"

#include <rtkRayQuadricIntersectionFunction.h>
#include <itkInPlaceImageFilter.h>

namespace pct
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT ProtonPairsToDistanceDrivenProjection :
  public itk::InPlaceImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ProtonPairsToDistanceDrivenProjection             Self;
  typedef itk::InPlaceImageFilter<TInputImage,TOutputImage> Superclass;
  typedef itk::SmartPointer<Self>                           Pointer;
  typedef itk::SmartPointer<const Self>                     ConstPointer;

  typedef itk::Vector<float, 3>                        ProtonPairsPixelType;
  typedef itk::Image<ProtonPairsPixelType,2>           ProtonPairsImageType;
  typedef ProtonPairsImageType::Pointer                ProtonPairsImagePointer;

  typedef itk::Image<float, 3>                  CountImageType;
  typedef CountImageType::Pointer                      CountImagePointer;

  typedef TOutputImage                                 OutputImageType;
  typedef typename OutputImageType::Pointer            OutputImagePointer;
  typedef typename OutputImageType::RegionType         OutputImageRegionType;

  typedef rtk::RayQuadricIntersectionFunction<double,3> RQIType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ProtonPairsToDistanceDrivenProjection, itk::InPlaceImageFilter);

  /** Get/Set image of proton pairs. */
  itkGetMacro(ProtonPairsFileName, std::string);
  itkSetMacro(ProtonPairsFileName, std::string);

  /** Get/Set the source position. */
  itkGetMacro(SourceDistance, double);
  itkSetMacro(SourceDistance, double);

  /** Get/Set the most likely path type. Can be "schulte" or "polynomial" */
  itkGetMacro(MostLikelyPathType, std::string);
  itkSetMacro(MostLikelyPathType, std::string);

  /** Get/Set the boundaries of the object. */
  itkGetMacro(QuadricIn, RQIType::Pointer);
  itkSetMacro(QuadricIn, RQIType::Pointer);
  itkGetMacro(QuadricOut, RQIType::Pointer);
  itkSetMacro(QuadricOut, RQIType::Pointer);

  /** Get/Set the count of proton pairs per pixel. */
  itkGetMacro(Count, CountImagePointer);

  /** Get/Set the ionization potential used in the Bethe-Bloch equation. */
  itkGetMacro(IonizationPotential, double);
  itkSetMacro(IonizationPotential, double);

protected:
  ProtonPairsToDistanceDrivenProjection() {}
  virtual ~ProtonPairsToDistanceDrivenProjection() {}

  virtual void BeforeThreadedGenerateData();
  virtual void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, rtk::ThreadIdType threadId );
  virtual void AfterThreadedGenerateData();

  /** The two inputs should not be in the same space so there is nothing
   * to verify. */
  virtual void VerifyInputInformation() {}

private:
  ProtonPairsToDistanceDrivenProjection(const Self&); //purposely not implemented
  void operator=(const Self&);            //purposely not implemented

  std::string m_ProtonPairsFileName;
  double m_SourceDistance;
  std::string m_MostLikelyPathType;

  /** Count event in each thread */
  CountImagePointer m_Count;
  std::vector<CountImagePointer> m_Counts;

  /** Create one output per thread */
  std::vector<OutputImagePointer> m_Outputs;

  /** The two quadric functions defining the object support. */
  RQIType::Pointer m_QuadricIn;
  RQIType::Pointer m_QuadricOut;

  /** Ionization potential used in the Bethe Bloch equation */
  double m_IonizationPotential;

  /** The functor to convert energy loss to attenuation */
  Functor::IntegratedBetheBlochProtonStoppingPowerInverse<float, double> *m_ConvFunc;

  ProtonPairsImageType::Pointer m_ProtonPairs;
};

} // end namespace pct

#ifndef ITK_MANUAL_INSTANTIATION
#include "pctProtonPairsToDistanceDrivenProjection.txx"
#endif

#endif
