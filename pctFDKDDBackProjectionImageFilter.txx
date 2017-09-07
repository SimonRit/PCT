#ifndef __pctFDKDDBackProjectionImageFilter_txx
#define __pctFDKDDBackProjectionImageFilter_txx

#include <itkImageRegionIteratorWithIndex.h>
#include <itkLinearInterpolateImageFunction.h>
#include <rtkThreeDCircularProjectionGeometry.h>

namespace pct
{

/**
 * GenerateData performs the accumulation
 */
template <class TInputImage, class TOutputImage>
void
FDKDDBackProjectionImageFilter<TInputImage,TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       rtk::ThreadIdType itkNotUsed(threadId) )
{
  const unsigned int Dimension = TInputImage::ImageDimension;
  const unsigned int nProj = m_ProjectionStack->GetLargestPossibleRegion().GetSize(Dimension);
  const unsigned int iFirstProj = m_ProjectionStack->GetLargestPossibleRegion().GetIndex(Dimension);
  rtk::ThreeDCircularProjectionGeometry * geometry;
  geometry = dynamic_cast<rtk::ThreeDCircularProjectionGeometry *>(this->GetGeometry().GetPointer());

  // Create interpolator, could be any interpolation
  typedef itk::LinearInterpolateImageFunction< ProjectionImageType, double > InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

  // Iterators on volume input and output
  typedef itk::ImageRegionConstIterator<TInputImage> InputRegionIterator;
  InputRegionIterator itIn(this->GetInput(), outputRegionForThread);
  typedef itk::ImageRegionIteratorWithIndex<TOutputImage> OutputRegionIterator;
  OutputRegionIterator itOut(this->GetOutput(), outputRegionForThread);

  // Initialize output region with input region in case the filter is not in
  // place
  if(this->GetInput() != this->GetOutput() )
    {
    itIn.GoToBegin();
    while(!itIn.IsAtEnd() )
      {
      itOut.Set(itIn.Get() );
      ++itIn;
      ++itOut;
      }
    }

  // Rotation center (assumed to be at 0 yet)
  typename TInputImage::PointType rotCenterPoint;
  rotCenterPoint.Fill(0.0);
  itk::ContinuousIndex<double, Dimension> rotCenterIndex;
  this->GetInput(0)->TransformPhysicalPointToContinuousIndex(rotCenterPoint, rotCenterIndex);

  // Continuous index at which we interpolate
  itk::ContinuousIndex<double, Dimension> pointProj;

  // Go over each projection
  for(unsigned int iProj=iFirstProj; iProj<iFirstProj+nProj; iProj++)
    {
    const double sid = geometry->GetSourceToIsocenterDistances()[iProj];
    const double sdd = geometry->GetSourceToDetectorDistances()[iProj];

    // Extract the current slice
    ProjectionImagePointer projection = this->GetDDProjection(iProj);
    interpolator->SetInputImage(projection);

    // Index to index matrix normalized to have a correct backprojection weight
    // (1 at the isocenter)
    ProjectionMatrixType matrix;
    if(sdd==0.)
      {
      // Create a 3D translation matrix for the projection instead of the 2D because
      // we want to keep the z coordinate
      itk::Matrix<double, Dimension+1, Dimension+1> matProjTrans;
      double sourceOffsetX = geometry->GetSourceOffsetsX()[iProj];
      double sourceOffsetY = geometry->GetSourceOffsetsY()[iProj];
      double projOffsetX = geometry->GetProjectionOffsetsX()[iProj];
      double projOffsetY = geometry->GetProjectionOffsetsY()[iProj];
      matProjTrans = rtk::ThreeDCircularProjectionGeometry::ComputeTranslationHomogeneousMatrix(sourceOffsetX-projOffsetX, sourceOffsetY-projOffsetY, 0.);

      itk::Matrix<double, Dimension+1, Dimension+1> matrixVol =
            rtk::GetIndexToPhysicalPointMatrix< TOutputImage >( this->GetOutput() );
      itk::Matrix<double, Dimension+1, Dimension+1> matrixStackProj =
            rtk::GetPhysicalPointToIndexMatrix< TOutputImage >( projection );
      matrixVol = matrixStackProj.GetVnlMatrix() *
                  matProjTrans.GetVnlMatrix() *
                  geometry->GetSourceTranslationMatrices()[iProj].GetVnlMatrix()*
                  geometry->GetRotationMatrices()[iProj].GetVnlMatrix() *
                  matrixVol.GetVnlMatrix();

      // Remove last line that is not used in the following
      for(unsigned int i=0; i<Dimension; i++)
        for(unsigned int j=0; j<Dimension+1; j++)
          matrix[i][j] = matrixVol[i][j];
      }
    else
      {
      matrix = GetIndexToIndexProjectionMatrix(iProj, projection);
      double perspFactor = matrix[Dimension-1][Dimension];
      for(unsigned int j=0; j<Dimension; j++)
        perspFactor += matrix[Dimension-1][j] * rotCenterIndex[j];
      matrix /= perspFactor;
      }


    // Go over each voxel
    itOut.GoToBegin();
    while(!itOut.IsAtEnd() )
      {
      if(sdd==0.) // Parallel
        {
        // Compute projection index
        for(unsigned int i=0; i<Dimension; i++)
          {
          pointProj[i] = matrix[i][Dimension];
          for(unsigned int j=0; j<Dimension; j++)
            pointProj[i] += matrix[i][j] * itOut.GetIndex()[j];
          }

        // Interpolate if in projection
        if( interpolator->IsInsideBuffer(pointProj) )
          {
          itOut.Set( itOut.Get() + interpolator->EvaluateAtContinuousIndex(pointProj) );
          }
        }
      else // Cone-beam
        {
        // Compute projection index
        for(unsigned int i=0; i<Dimension-1; i++)
          {
          pointProj[i] = matrix[i][Dimension];
          for(unsigned int j=0; j<Dimension; j++)
            pointProj[i] += matrix[i][j] * itOut.GetIndex()[j];
          }

        // Apply perspective
        double perspFactor = matrix[Dimension-1][Dimension];
        for(unsigned int j=0; j<Dimension; j++)
          perspFactor += matrix[Dimension-1][j] * itOut.GetIndex()[j];
        perspFactor = 1/perspFactor;
        for(unsigned int i=0; i<Dimension-1; i++)
          pointProj[i] = pointProj[i]*perspFactor;

        // Distance driven
        pointProj[Dimension-1] = (sid / perspFactor - sid - m_ProjectionStack->GetOrigin()[Dimension-1]) / m_ProjectionStack->GetSpacing()[Dimension-1];
        // Interpolate if in projection
        if( interpolator->IsInsideBuffer(pointProj) )
          {
          itOut.Set( itOut.Get() + perspFactor*perspFactor*interpolator->EvaluateAtContinuousIndex(pointProj) );
          }
        }

      ++itOut;
      }
    }
}

template <class TInputImage, class TOutputImage>
typename FDKDDBackProjectionImageFilter<TInputImage,TOutputImage>::ProjectionImagePointer
FDKDDBackProjectionImageFilter<TInputImage,TOutputImage>
::GetDDProjection(const unsigned int iProj)
{
  const int iProjBuff = m_ProjectionStack->GetBufferedRegion().GetIndex(ProjectionImageType::ImageDimension);

  ProjectionImagePointer projection = ProjectionImageType::New();
  typename ProjectionImageType::RegionType region;
  typename ProjectionImageType::SpacingType spacing;
  typename ProjectionImageType::PointType origin;

  for(unsigned int i=0; i<ProjectionImageType::ImageDimension; i++)
    {
    origin[i] = m_ProjectionStack->GetOrigin()[i];
    spacing[i] = m_ProjectionStack->GetSpacing()[i];
    region.SetSize(i, m_ProjectionStack->GetLargestPossibleRegion().GetSize()[i]);
    region.SetIndex(i, m_ProjectionStack->GetLargestPossibleRegion().GetIndex()[i]);
    }
  projection->SetSpacing(spacing);
  projection->SetOrigin(origin);
  projection->SetRegions(region);
  projection->Allocate();

  const unsigned int    npixels = projection->GetLargestPossibleRegion().GetNumberOfPixels();
  const ProjectionPixelType *pi = m_ProjectionStack->GetBufferPointer() + (iProj-iProjBuff)*npixels;
  ProjectionPixelType *      po = projection->GetBufferPointer();

  // Transpose projection for optimization
  for(unsigned int i=0; i<npixels; i++)
    *po++ = *pi++;

  return projection;
}

template <class TInputImage, class TOutputImage>
typename rtk::BackProjectionImageFilter<TInputImage,TOutputImage>::ProjectionMatrixType
FDKDDBackProjectionImageFilter<TInputImage,TOutputImage>
::GetIndexToIndexProjectionMatrix(const unsigned int iProj, const ProjectionImageType *proj)
{
  const unsigned int Dimension = TInputImage::ImageDimension;

  itk::Matrix<double, Dimension+1, Dimension+1> matrixVol;
  matrixVol = rtk::GetIndexToPhysicalPointMatrix< TOutputImage >( this->GetOutput() );
  itk::Matrix<double, Dimension+1, Dimension+1> matrixProjFull;
  matrixProjFull = rtk::GetPhysicalPointToIndexMatrix< ProjectionImageType >(proj);
  itk::Matrix<double, Dimension, Dimension> matrixProj;
  matrixProj.SetIdentity();
  for(unsigned int j=0; j<Dimension-1; j++)
    {
    for(unsigned int i=0; i<Dimension-1; i++)
      matrixProj[i][j] = matrixProjFull[i][j];
    matrixProj[j][Dimension-1] = matrixProjFull[j][Dimension];
    }

  return ProjectionMatrixType(matrixProj.GetVnlMatrix() *
                              this->GetGeometry()->GetMatrices()[iProj].GetVnlMatrix() *
                              matrixVol.GetVnlMatrix() );
}

} // end namespace pct

#endif
