#ifndef __pctProtonPairsToBackProjection_h
#define __pctProtonPairsToBackProjection_h

#include "rtkConfiguration.h"
#include "pctBetheBlochFunctor.h"

#include <rtkRayQuadricIntersectionFunction.h>
#include <rtkThreeDCircularProjectionGeometry.h>
#include <itkInPlaceImageFilter.h>
#include <itkBarrier.h>

namespace pct
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT ProtonPairsToBackProjection :
  public itk::InPlaceImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ProtonPairsToBackProjection                       Self;
  typedef itk::InPlaceImageFilter<TInputImage,TOutputImage> Superclass;
  typedef itk::SmartPointer<Self>                           Pointer;
  typedef itk::SmartPointer<const Self>                     ConstPointer;

  typedef  std::vector<std::string> FileNamesContainer;

  typedef itk::Vector<float, 3>                        ProtonPairsPixelType;
  typedef itk::Image<ProtonPairsPixelType,2>           ProtonPairsImageType;
  typedef ProtonPairsImageType::Pointer                ProtonPairsImagePointer;

  typedef typename itk::Image< unsigned int,
                               TInputImage::ImageDimension> CountImageType;
  typedef typename CountImageType::Pointer                  CountImagePointer;

  typedef TOutputImage                                 OutputImageType;
  typedef typename OutputImageType::Pointer            OutputImagePointer;
  typedef typename OutputImageType::RegionType         OutputImageRegionType;

  typedef rtk::RayQuadricIntersectionFunction<double,3> RQIType;

  typedef rtk::ThreeDCircularProjectionGeometry     GeometryType;
  typedef typename GeometryType::Pointer            GeometryPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ProtonPairsToBackProjection, itk::InPlaceImageFilter);

  /** Set the vector of strings that contains the file names. Files
   * are processed in sequential order. */
  void SetProtonPairsFileNames (const FileNamesContainer &name)
    {
    if ( m_ProtonPairsFileNames != name)
      {
      m_ProtonPairsFileNames = name;
      this->Modified();
      }
    }
  const FileNamesContainer & GetProtonPairsFileNames() const
    {
    return m_ProtonPairsFileNames;
    }

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
  itkSetMacro(Count, CountImagePointer);

  /** Get/Set the ionization potential used in the Bethe-Bloch equation. */
  itkGetMacro(IonizationPotential, double);
  itkSetMacro(IonizationPotential, double);

  /** Get / Set the object pointer to projection geometry */
  itkGetMacro(Geometry, GeometryPointer);
  itkSetMacro(Geometry, GeometryPointer);

protected:
  ProtonPairsToBackProjection() {}
  virtual ~ProtonPairsToBackProjection() {}

  virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;
  virtual void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, rtk::ThreadIdType threadId ) ITK_OVERRIDE;
  virtual void AfterThreadedGenerateData() ITK_OVERRIDE;

  /** The two inputs should not be in the same space so there is nothing
   * to verify. */
  virtual void VerifyInputInformation() ITK_OVERRIDE {}

private:
  ProtonPairsToBackProjection(const Self&); //purposely not implemented
  void operator=(const Self&);            //purposely not implemented

  /** A list of filenames to be processed. */
  FileNamesContainer m_ProtonPairsFileNames;
  itk::Barrier::Pointer m_Barriers[3];

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

  /** RTK geometry object */
  GeometryPointer m_Geometry;
};

} // end namespace pct

#ifndef ITK_MANUAL_INSTANTIATION
#include "pctProtonPairsToBackProjection.txx"
#endif

#endif
