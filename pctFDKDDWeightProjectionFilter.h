#ifndef __pctFDKDDWeightProjectionFilter_h
#define __pctFDKDDWeightProjectionFilter_h

#include <itkInPlaceImageFilter.h>
#include <rtkThreeDCircularProjectionGeometry.h>

/** \class FDKDDWeightProjectionFilter
 * \brief Weighting of projections to correct for the divergence in
 * filtered backprojection reconstruction algorithms.
 * The weighting comprises:
 * - the 2D weighting of the FDK algorithm [Feldkamp, 1984],
 * - the correction of the ramp factor for divergent full scan,
 * - the angular weighting for the final 3D integral of FDK.
 * Note that SourceToDetectorDistance, SourceToDetectorIsocenter
 * SouceOffsets and ProjectionOffsets are accounted for on a per
 * projection basis but InPlaneRotation and OutOfPlaneRotation are not
 * accounted for.
 * \author Simon Rit
 */
namespace pct
{

template<class TInputImage, class TOutputImage=TInputImage>
class ITK_EXPORT FDKDDWeightProjectionFilter :
  public itk::InPlaceImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef FDKDDWeightProjectionFilter Self;
  typedef itk::ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                          InputImageType;
  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(FDKDDWeightProjectionFilter, itk::ImageToImageFilter);

  /** Get/ Set geometry structure */
  itkGetMacro(Geometry, rtk::ThreeDCircularProjectionGeometry::Pointer);
  itkSetMacro(Geometry, rtk::ThreeDCircularProjectionGeometry::Pointer);

protected:
  FDKDDWeightProjectionFilter()  {}
  ~FDKDDWeightProjectionFilter() {}

  virtual void BeforeThreadedGenerateData();

  virtual void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, rtk::ThreadIdType threadId);

private:
  FDKDDWeightProjectionFilter(const Self&); //purposely not implemented
  void operator=(const Self&);            //purposely not implemented

  /** Angular weights for each projection */
  std::vector<double> m_AngularWeightsAndRampFactor;

  /** Geometrical description of the system */
  rtk::ThreeDCircularProjectionGeometry::Pointer m_Geometry;
}; // end of class

} // end namespace pct

#ifndef ITK_MANUAL_INSTANTIATION
#include "pctFDKDDWeightProjectionFilter.txx"
#endif

#endif
