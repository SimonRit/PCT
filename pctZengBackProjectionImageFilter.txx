#ifndef __pctZengBackProjectionImageFilter_txx
#define __pctZengBackProjectionImageFilter_txx

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkMacro.h>
#include <itkImageRegionIteratorWithIndex.h>

namespace pct
{

template <class TInputImage, class TOutputImage>
ZengBackProjectionImageFilter<TInputImage, TOutputImage>
::ZengBackProjectionImageFilter()
{
  this->SetNumberOfRequiredOutputs(2);
  this->SetNthOutput( 0, this->MakeOutput(0) );
  this->SetNthOutput( 1, this->MakeOutput(1) );
}

template <class TInputImage, class TOutputImage>
void
ZengBackProjectionImageFilter<TInputImage, TOutputImage>
::GenerateOutputInformation()
{
  typename OutputImageType::RegionType region;
  typename OutputImageType::PointType origin;
  typename OutputImageType::SpacingType spacing;

  for(unsigned int i=0; i<TOutputImage::ImageDimension; i++)
    {
    region.SetIndex(i, this->GetInput()->GetLargestPossibleRegion().GetIndex(i));
    region.SetSize(i, this->GetInput()->GetLargestPossibleRegion().GetSize(i));
    spacing[i] = this->GetInput()->GetSpacing()[i];
    origin[i] = this->GetInput()->GetOrigin()[i];
    }

  for(unsigned int i=0; i<this->GetNumberOfRequiredOutputs(); i++)
    {
    this->GetOutput(i)->SetSpacing(spacing);
    this->GetOutput(i)->SetOrigin(origin);
    this->GetOutput(i)->SetRegions(region);
    }
}

template <class TInputImage, class TOutputImage>
void
ZengBackProjectionImageFilter<TInputImage, TOutputImage>
#if ITK_VERSION_MAJOR <= 4
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       rtk::ThreadIdType itkNotUsed(threadId) )
#else
::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread)
#endif
{
  typename TInputImage::IndexType idxIn;
  for(unsigned int i=0; i<TOutputImage::ImageDimension; i++)
    idxIn[i] = outputRegionForThread.GetIndex(i);
  idxIn[TOutputImage::ImageDimension] = 0;
  const typename TInputImage::PixelType *pIn = this->GetInput()->GetBufferPointer() + this->GetInput()->ComputeOffset(idxIn);
  typename OutputImageRegionType::SizeValueType npix = this->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels();

  // Input / ouput iterators
  itk::ImageRegionIterator<OutputImageType>     itOutC(this->GetOutput(0), outputRegionForThread);
  itk::ImageRegionIterator<OutputImageType>     itOutS(this->GetOutput(1), outputRegionForThread);

  double angspac = itk::Math::pi / this->GetInput()->GetLargestPossibleRegion().GetSize(TInputImage::ImageDimension-1);
  const double phi = 0.;
  while(!itOutC.IsAtEnd())
    {
    double ang = angspac * 0.5 + itk::Math::pi_over_2;
    double vs = 0.;
    double vc = 0.;
    const typename TInputImage::PixelType *pInCurr = pIn;
    for(unsigned int i=0; i<this->GetInput()->GetLargestPossibleRegion().GetSize(TInputImage::ImageDimension-1); i++)
      {
      while(ang>itk::Math::pi) ang -= itk::Math::pi;
      double sign = 1.;
      if(sin(ang-phi)<0.)
        sign = -1.;
      vs += sin(ang) * (*pInCurr) * sign;
      vc += cos(ang) * (*pInCurr) * sign;
      ang += angspac;
      pInCurr += npix;
      }
    itOutC.Set(vc * angspac);
    itOutS.Set(-1. * vs * angspac);
    ++itOutC;
    ++itOutS;
    ++pIn;
    }
}

template <class TInputImage, class TOutputImage>
itk::DataObject::Pointer
ZengBackProjectionImageFilter<TInputImage, TOutputImage>
::MakeOutput(itk::ProcessObject::DataObjectPointerArraySizeType idx)
{
  itk::DataObject::Pointer output;

  switch ( idx )
    {
    case 0:
      output = ( TOutputImage::New() ).GetPointer();
      break;
    case 1:
      output = ( TOutputImage::New() ).GetPointer();
      break;
    default:
      std::cerr << "No output " << idx << std::endl;
      output = NULL;
      break;
    }
  return output.GetPointer();
}

} // end namespace pct
#endif
