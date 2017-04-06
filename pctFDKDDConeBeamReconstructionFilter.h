#ifndef __pctFDKDDConeBeamReconstructionFilter_h
#define __pctFDKDDConeBeamReconstructionFilter_h

#include "pctFDKDDWeightProjectionFilter.h"
#include <rtkFFTRampImageFilter.h>
#include "pctFDKDDBackProjectionImageFilter.h"

#include <itkExtractImageFilter.h>
#include <itkTimeProbe.h>

/** \class FDKDDConeBeamReconstructionFilter
 * TODO
 *
 * \author Simon Rit
 */
namespace pct
{

template<class TInputImage, class TOutputImage=TInputImage, class TFFTPrecision=double>
class ITK_EXPORT FDKDDConeBeamReconstructionFilter :
  public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef FDKDDConeBeamReconstructionFilter                  Self;
  typedef itk::ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef itk::SmartPointer<Self>                            Pointer;
  typedef itk::SmartPointer<const Self>                      ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;

  typedef itk::Image<float, 4>                               ProjectionStackType;
  typedef ProjectionStackType::Pointer                       ProjectionStackPointer;

  /** Typedefs of each subfilter of this composite filter */
  typedef itk::ExtractImageFilter< ProjectionStackType, ProjectionStackType >                ExtractFilterType;
  typedef pct::FDKDDWeightProjectionFilter< ProjectionStackType, ProjectionStackType >       WeightFilterType;
  typedef rtk::FFTRampImageFilter< ProjectionStackType, ProjectionStackType, TFFTPrecision > RampFilterType;
  typedef pct::FDKDDBackProjectionImageFilter< OutputImageType, OutputImageType >            BackProjectionFilterType;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(FDKDDConeBeamReconstructionFilter, itk::ImageToImageFilter);

  /** Get / Set the object pointer to projection geometry */
  virtual rtk::ThreeDCircularProjectionGeometry::Pointer GetGeometry();
  virtual void SetGeometry(const rtk::ThreeDCircularProjectionGeometry::Pointer _arg);

  /** Get pointer to the ramp filter used by the feldkamp reconstruction */
  typename RampFilterType::Pointer GetRampFilter() { return m_RampFilter; }

  void PrintTiming(std::ostream& os) const;

  /** Get / Set the stack of projection images */
  itkGetMacro(ProjectionStack, ProjectionStackPointer);
  itkSetMacro(ProjectionStack, ProjectionStackPointer);

protected:
  FDKDDConeBeamReconstructionFilter();
  ~FDKDDConeBeamReconstructionFilter(){}

  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

  void GenerateOutputInformation() ITK_OVERRIDE;

  void GenerateData() ITK_OVERRIDE;

  /** The two inputs should not be in the same space so there is nothing
   * to verify. */
  virtual void VerifyInputInformation() ITK_OVERRIDE {}

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

} // end namespace pct

#ifndef ITK_MANUAL_INSTANTIATION
#include "pctFDKDDConeBeamReconstructionFilter.txx"
#endif

#endif
