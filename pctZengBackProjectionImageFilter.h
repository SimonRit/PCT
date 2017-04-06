#ifndef __pctZengBackProjectionImageFilter_h
#define __pctZengBackProjectionImageFilter_h

#include <itkImageToImageFilter.h>
#include <rtkThreeDCircularProjectionGeometry.h>
#include "rtkConfiguration.h"

/** \class ZengBackProjectionImageFilter
 *
 * From an input 4D image where the 4th dimension is the angle, computes
 * the weighted backprojection for DBP described in [Zeng, Med Phys, 2007].
 *
 * \author Simon Rit
 */
namespace pct
{

template<class TInputImage, class TOutputImage=itk::Image<typename TInputImage::PixelType, TInputImage::ImageDimension-1> >
class ITK_EXPORT ZengBackProjectionImageFilter:
  public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ZengBackProjectionImageFilter                      Self;
  typedef itk::ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef itk::SmartPointer<Self>                            Pointer;
  typedef itk::SmartPointer<const Self>                      ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                                     InputImageType;
  typedef TOutputImage                                    OutputImageType;
  typedef typename TOutputImage::RegionType               OutputImageRegionType;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(ZengBackProjectionImageFilter, itk::ImageToImageFilter);

protected:
  ZengBackProjectionImageFilter();
  ~ZengBackProjectionImageFilter(){}

  virtual void GenerateOutputInformation() ITK_OVERRIDE;
  virtual void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, rtk::ThreadIdType threadId) ITK_OVERRIDE;

  virtual itk::DataObject::Pointer MakeOutput(itk::ProcessObject::DataObjectPointerArraySizeType idx) ITK_OVERRIDE;

private:
  ZengBackProjectionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&);                //purposely not implemented
}; // end of class

} // end namespace pct

#ifndef ITK_MANUAL_INSTANTIATION
#include "pctZengBackProjectionImageFilter.txx"
#endif

#endif
