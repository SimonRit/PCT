#ifndef __itkFDKDDConeBeamReconstructionFilter_h
#define __itkFDKDDConeBeamReconstructionFilter_h

#include "itkFDKDDWeightProjectionFilter.h"
#include "itkFFTRampImageFilter.h"
#include "itkFDKDDBackProjectionImageFilter.h"

#include <itkExtractImageFilter.h>
#include <itkTimeProbe.h>

/** \class FDKDDConeBeamReconstructionFilter
 * TODO
 *
 * \author Simon Rit
 */
namespace itk
{

template<class TInputImage, class TOutputImage=TInputImage, class TFFTPrecision=double>
class ITK_EXPORT FDKDDConeBeamReconstructionFilter :
  public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef FDKDDConeBeamReconstructionFilter             Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;

  typedef Image<float, 4>                                     ProjectionStackType;
  typedef ProjectionStackType::Pointer                        ProjectionStackPointer;

  /** Typedefs of each subfilter of this composite filter */
  typedef itk::ExtractImageFilter< ProjectionStackType, ProjectionStackType >                ExtractFilterType;
  typedef itk::FDKDDWeightProjectionFilter< ProjectionStackType, ProjectionStackType >       WeightFilterType;
  typedef itk::FFTRampImageFilter< ProjectionStackType, ProjectionStackType, TFFTPrecision > RampFilterType;
  typedef itk::FDKDDBackProjectionImageFilter< OutputImageType, OutputImageType >            BackProjectionFilterType;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(FDKDDConeBeamReconstructionFilter, ImageToImageFilter);

  /** Get / Set the object pointer to projection geometry */
  virtual ThreeDCircularProjectionGeometry::Pointer GetGeometry();
  virtual void SetGeometry(const ThreeDCircularProjectionGeometry::Pointer _arg);

  /** Get pointer to the ramp filter used by the feldkamp reconstruction */
  typename RampFilterType::Pointer GetRampFilter() { return m_RampFilter; }

  void PrintTiming(std::ostream& os) const;

  /** Get / Set the stack of projection images */
  itkGetMacro(ProjectionStack, ProjectionStackPointer);
  itkSetMacro(ProjectionStack, ProjectionStackPointer);

protected:
  FDKDDConeBeamReconstructionFilter();
  ~FDKDDConeBeamReconstructionFilter(){}

  virtual void GenerateInputRequestedRegion();

  void GenerateOutputInformation();

  void GenerateData();

  /** The two inputs should not be in the same space so there is nothing
   * to verify. */
  virtual void VerifyInputInformation() {}

  /** Pointers to each subfilter of this composite filter */
  typename ExtractFilterType::Pointer m_ExtractFilter;
  typename WeightFilterType::Pointer m_WeightFilter;
  typename RampFilterType::Pointer m_RampFilter;
  typename BackProjectionFilterType::Pointer m_BackProjectionFilter;

private:
  //purposely not implemented
  FDKDDConeBeamReconstructionFilter(const Self&);
  void operator=(const Self&);

  /** Probes to time reconstruction */
  itk::TimeProbe m_PreFilterProbe;
  itk::TimeProbe m_FilterProbe;
  itk::TimeProbe m_BackProjectionProbe;

  ProjectionStackPointer m_ProjectionStack;
}; // end of class

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFDKDDConeBeamReconstructionFilter.txx"
#endif

#endif
