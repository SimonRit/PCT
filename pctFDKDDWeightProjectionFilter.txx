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
  m_AngularWeightsAndRampFactor = this->GetGeometry()->GetAngularGaps( m_Geometry->GetGantryAngles() );

  for(unsigned int k=0; k<m_AngularWeightsAndRampFactor.size(); k++)
    {
    // Add correction factor for ramp filter
    const double sdd  = m_Geometry->GetSourceToDetectorDistances()[k];
    if(sdd==0.) // Parallel
      m_AngularWeightsAndRampFactor[k] *= 0.5;
    else        // Divergent
      {
      // Zoom + factor 1/2 in eq 176, page 106, Kak & Slaney
      const double sid  = m_Geometry->GetSourceToIsocenterDistances()[k];
      const double rampFactor = sdd / (2. * sid);
      m_AngularWeightsAndRampFactor[k] *= rampFactor;
      }
    }
}

template <class TInputImage, class TOutputImage>
void
FDKDDWeightProjectionFilter<TInputImage, TOutputImage>
::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread)
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
                 + m_Geometry->GetProjectionOffsetsY()[l]
                 - m_Geometry->GetSourceOffsetsY()[l];
      const double sdd  = m_Geometry->GetSourceToDetectorDistances()[l];
      const double sid  = m_Geometry->GetSourceToIsocenterDistances()[l];
      const double sdd2 = sdd * sdd;
      if(sdd != 0.) // Divergent
        {
        const double tauOverD  = m_Geometry->GetSourceOffsetsX()[l] / sid;
        const double tauOverDw = m_AngularWeightsAndRampFactor[l] * tauOverD;
        const double sddw      = m_AngularWeightsAndRampFactor[l] * sdd;
        for(unsigned int j=0;
                         j<outputRegionForThread.GetSize(1);
                         j++, point[1] += pointIncrement[1])
          {
          point[0] = pointBase[0]
                     + m_Geometry->GetProjectionOffsetsX()[l]
                     - m_Geometry->GetSourceOffsetsX()[l];
          const double sdd2y2 = sdd2 + point[1]*point[1];
          for(unsigned int i=0;
                           i<outputRegionForThread.GetSize(0);
                           i++, ++itI, ++itO, point[0] += pointIncrement[0])
            {
            // The term between parentheses comes from the publication
            // [Gullberg Crawford Tsui, TMI, 1986], equation 18
            itO.Set( itI.Get() * (sddw - tauOverDw * point[0]) / sqrt( sdd2y2 + point[0]*point[0]) );
            }
          }
        }
      else // Parallel
        {
        double weight = m_AngularWeightsAndRampFactor[l];
        for(unsigned int j=0; j<outputRegionForThread.GetSize(1); j++)
          {
          for(unsigned int i=0; i<outputRegionForThread.GetSize(0); i++, ++itI, ++itO)
            itO.Set( itI.Get() * weight);
          }
        }
      }
    }
}

} // end namespace pct
#endif
