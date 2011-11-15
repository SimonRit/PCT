#ifndef __itkFDKDDBackProjectionImageFilter_h
#define __itkFDKDDBackProjectionImageFilter_h

#include "itkInPlaceImageFilter.h"
#include "itkBackProjectionImageFilter.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
class ITK_EXPORT FDKDDBackProjectionImageFilter :
  public BackProjectionImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef FDKDDBackProjectionImageFilter                      Self;
  typedef BackProjectionImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  typedef typename Superclass::ProjectionMatrixType           ProjectionMatrixType;                                                                                                      typedef typename TOutputImage::RegionType                                          OutputImageRegionType;
  typedef TInputImage                                         ProjectionImageType;
  typedef typename ProjectionImageType::Pointer               ProjectionImagePointer;
  typedef typename ProjectionImageType::PixelType             ProjectionPixelType;
  typedef Image<float, 4>                                     ProjectionStackType;
  typedef ProjectionStackType::Pointer                        ProjectionStackPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FDKDDBackProjectionImageFilter, ImageToImageFilter);

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

  virtual void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId );

  ProjectionStackPointer m_ProjectionStack;

private:
  FDKDDBackProjectionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&);               //purposely not implemented
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFDKDDBackProjectionImageFilter.txx"
#endif

#endif
