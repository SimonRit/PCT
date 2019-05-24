#ifndef __pctFDKDDBackProjectionImageFilter_h
#define __pctFDKDDBackProjectionImageFilter_h

#include <itkInPlaceImageFilter.h>
#include <rtkBackProjectionImageFilter.h>

namespace pct
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT FDKDDBackProjectionImageFilter :
  public rtk::BackProjectionImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef FDKDDBackProjectionImageFilter                           Self;
  typedef rtk::BackProjectionImageFilter<TInputImage,TOutputImage> Superclass;
  typedef itk::SmartPointer<Self>                                  Pointer;
  typedef itk::SmartPointer<const Self>                            ConstPointer;

  typedef typename Superclass::ProjectionMatrixType           ProjectionMatrixType;
  typedef typename TOutputImage::RegionType                   OutputImageRegionType;
  typedef TInputImage                                         ProjectionImageType;
  typedef typename ProjectionImageType::Pointer               ProjectionImagePointer;
  typedef typename ProjectionImageType::PixelType             ProjectionPixelType;
  typedef itk::Image<float, 4>                                ProjectionStackType;
  typedef ProjectionStackType::Pointer                        ProjectionStackPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FDKDDBackProjectionImageFilter, rtk::BackProjectionImageFilter);

  virtual ProjectionImagePointer GetDDProjection(const unsigned int iProj);

  /** Get / Set the stack of projection images */
  itkGetMacro(ProjectionStack, ProjectionStackPointer);
  itkSetMacro(ProjectionStack, ProjectionStackPointer);

  /** Creates the #iProj index to index projection matrix with current inputs
      instead of the physical point to physical point projection matrix provided by Geometry */
  ProjectionMatrixType GetIndexToIndexProjectionMatrix(const unsigned int iProj, const ProjectionImageType *proj);

protected:
  FDKDDBackProjectionImageFilter() {
    this->SetNumberOfRequiredInputs(1); this->SetInPlace( true );
  };
  virtual ~FDKDDBackProjectionImageFilter() {
  };

#if ITK_VERSION_MAJOR <= 4
  virtual void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, rtk::ThreadIdType threadId ) ITK_OVERRIDE;
#else
  virtual void DynamicThreadedGenerateData( const OutputImageRegionType& outputRegionForThread ) ITK_OVERRIDE;
#endif

  ProjectionStackPointer m_ProjectionStack;

private:
  FDKDDBackProjectionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&);               //purposely not implemented
};

} // end namespace pct

#ifndef ITK_MANUAL_INSTANTIATION
#include "pctFDKDDBackProjectionImageFilter.txx"
#endif

#endif
