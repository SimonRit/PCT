#ifndef __pctFDKDDWeightProjectionFilter_txx
#define __pctFDKDDWeightProjectionFilter_txx

#include <itkImageRegionIterator.h>
#include <itkImageFileWriter.h>
#include <itkMeanProjectionImageFilter.h>

namespace pct
{
template <class TInputImage, class TOutputImage>
void
FDKDDWeightProjectionFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  // Get angular weights from geometry
  m_AngularWeightsAndRampFactor = this->GetGeometry()->GetAngularGaps();

  for(unsigned int k=0; k<m_AngularWeightsAndRampFactor.size(); k++)
    {
    // Add correction factor for ramp filter
    const double sid  = m_Geometry->GetSourceToIsocenterDistances()[k];
    const double sdd  = m_Geometry->GetSourceToDetectorDistances()[k];
    // Zoom + factor 1/2 in eq 176, page 106, Kak & Slaney
    const double rampFactor = sdd / (2 * sid);
    m_AngularWeightsAndRampFactor[k] *= rampFactor;
    }
}

template <class TInputImage, class TOutputImage>
void
FDKDDWeightProjectionFilter<TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       rtk::ThreadIdType itkNotUsed(threadId) )
{
  // Prepare point increment (TransformIndexToPhysicalPoint too slow)
  typename InputImageType::PointType pointBase, pointIncrement;
  typename InputImageType::IndexType index = outputRegionForThread.GetIndex();
  this->GetInput()->TransformIndexToPhysicalPoint( index, pointBase );
  for(int i=0; i<3; i++)
    index[i]++;
  this->GetInput()->TransformIndexToPhysicalPoint( index, pointIncrement );
  for(int i=0; i<3; i++)
    pointIncrement[i] -= pointBase[i];

  // Iterators
  typedef itk::ImageRegionConstIterator<InputImageType> InputConstIterator;
  InputConstIterator itI(this->GetInput(), outputRegionForThread);
  typedef itk::ImageRegionIterator<OutputImageType> OutputIterator;
  OutputIterator itO(this->GetOutput(), outputRegionForThread);
  itI.GoToBegin();
  itO.GoToBegin();

  // Go over output, compute weights and avoid redundant computation
  for(unsigned int l=outputRegionForThread.GetIndex(3);
                   l<outputRegionForThread.GetIndex(3)+outputRegionForThread.GetSize(3);
                   l++)
    {
    for(unsigned int k=outputRegionForThread.GetIndex(2);
                     k<outputRegionForThread.GetIndex(2)+outputRegionForThread.GetSize(2);
                     k++)
      {
      typename InputImageType::PointType point = pointBase;
      point[1] = pointBase[1]
                 + m_Geometry->GetSourceOffsetsY()[l]
                 - m_Geometry->GetProjectionOffsetsY()[l];
      const double sdd  = m_Geometry->GetSourceToDetectorDistances()[l];
      const double sdd2 = sdd * sdd;
      double weight = sdd  * m_AngularWeightsAndRampFactor[l];
      for(unsigned int j=outputRegionForThread.GetIndex(1);
                       j<outputRegionForThread.GetIndex(1)+outputRegionForThread.GetSize(1);
                       j++, point[1] += pointIncrement[1])
        {
        point[0] = pointBase[0]
                   + m_Geometry->GetSourceOffsetsX()[l]
                   - m_Geometry->GetProjectionOffsetsX()[l];
        const double sdd2y2 = sdd2 + point[1]*point[1];
        for(unsigned int i=outputRegionForThread.GetIndex(0);
                         i<outputRegionForThread.GetIndex(0)+outputRegionForThread.GetSize(0);
                         i++, ++itI, ++itO, point[0] += pointIncrement[0])
          {
          itO.Set( itI.Get() * weight / sqrt( sdd2y2 + point[0]*point[0]) );
          }
        }
      }
    }
}

} // end namespace pct
#endif
