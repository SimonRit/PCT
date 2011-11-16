#ifndef __itkProtonPairsToDistanceDrivenProjection_h
#define __itkProtonPairsToDistanceDrivenProjection_h

#include "rtkConfiguration.h"

#include <itkRayQuadricIntersectionFunction.h>
#include <itkInPlaceImageFilter.h>

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT ProtonPairsToDistanceDrivenProjection :
  public InPlaceImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ProtonPairsToDistanceDrivenProjection        Self;
  typedef InPlaceImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  typedef itk::Vector<float, 3>                        ProtonPairsPixelType;
  typedef itk::Image<ProtonPairsPixelType,2>           ProtonPairsImageType;
  typedef ProtonPairsImageType::Pointer                ProtonPairsImagePointer;

  typedef itk::Image<unsigned int, 3>                  CountImageType;
  typedef CountImageType::Pointer                      CountImagePointer;

  typedef TOutputImage                                 OutputImageType;
  typedef typename OutputImageType::Pointer            OutputImagePointer;
  typedef typename OutputImageType::RegionType         OutputImageRegionType;

  typedef RayQuadricIntersectionFunction<double,3>     RQIType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ProtonPairsToDistanceDrivenProjection, InPlaceImageFilter);

  /** Get/Set image of proton pairs. */
  itkGetMacro(ProtonPairsFileName, std::string);
  itkSetMacro(ProtonPairsFileName, std::string);

  /** Get/Set the source position. */
  itkGetMacro(SourceDistance, double);
  itkSetMacro(SourceDistance, double);

  /** Get/Set the boundaries of the object. */
  itkGetMacro(QuadricIn, RQIType::Pointer);
  itkSetMacro(QuadricIn, RQIType::Pointer);
  itkGetMacro(QuadricOut, RQIType::Pointer);
  itkSetMacro(QuadricOut, RQIType::Pointer);

  /** Get/Set the count of proton pairs per pixel. */
  itkGetMacro(Count, CountImagePointer);

protected:
  ProtonPairsToDistanceDrivenProjection() {
    this->SetInPlace( true );
//    this->SetNumberOfThreads(4);
  };
  virtual ~ProtonPairsToDistanceDrivenProjection() {}

  virtual void BeforeThreadedGenerateData();
  virtual void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId );
  virtual void AfterThreadedGenerateData();

  /** The two inputs should not be in the same space so there is nothing
   * to verify. */
  virtual void VerifyInputInformation() {}

private:
  ProtonPairsToDistanceDrivenProjection(const Self&); //purposely not implemented
  void operator=(const Self&);            //purposely not implemented

  std::string m_ProtonPairsFileName;
  double m_SourceDistance;

  /** Count event in each thread */
  CountImagePointer m_Count;
  std::vector<CountImagePointer> m_Counts;

  /** Create one output per thread */
  std::vector<OutputImagePointer> m_Outputs;

  /** The two quadric functions defining the object support. */
  RQIType::Pointer m_QuadricIn;
  RQIType::Pointer m_QuadricOut;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkProtonPairsToDistanceDrivenProjection.txx"
#endif

#endif
