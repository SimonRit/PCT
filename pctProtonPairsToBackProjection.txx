#include <itkImageFileReader.h>
#include <itkImageRegionIterator.h>

#include <rtkHomogeneousMatrix.h>

#include "pctThirdOrderPolynomialMLPFunction.h"
#include "pctSchulteMLPFunction.h"
#include "pctEnergyStragglingFunctor.h"

namespace pct
{

template <class TInputImage, class TOutputImage>
ProtonPairsToBackProjection<TInputImage, TOutputImage>
::ProtonPairsToBackProjection()
{
}

template <class TInputImage, class TOutputImage>
void
ProtonPairsToBackProjection<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  if(m_QuadricOut.GetPointer()==NULL)
    m_QuadricOut = m_QuadricIn;
  m_ConvFunc = new Functor::IntegratedBetheBlochProtonStoppingPowerInverse<float, double>(m_IonizationPotential, 600.*CLHEP::MeV, 0.1*CLHEP::keV);
}

template <class TInputImage, class TOutputImage>
void
ProtonPairsToBackProjection<TInputImage, TOutputImage>
::GenerateData()
{
  this->AllocateOutputs();
  this->BeforeThreadedGenerateData();

  // Create thread image and corresponding stack to count events
  m_Counts = CountImageType::New();
  m_Counts->SetRegions(this->GetInput()->GetLargestPossibleRegion());
  m_Counts->Allocate();
  m_Counts->FillBuffer(0);

  for(unsigned int iProj = 0; iProj<m_ProtonPairsFileNames.size(); iProj++)
    {
    std::cout << std::endl
              << "Reading "
              << m_ProtonPairsFileNames[iProj]
              << "... "
              << std::flush;
    // Read pairs
    typedef itk::ImageFileReader< ProtonPairsImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( m_ProtonPairsFileNames[iProj] );
    reader->Update();
    m_ProtonPairs = reader->GetOutput();

    std::cout << "Done !" << std::endl;

    this->GetMultiThreader()->template ParallelizeImageRegion<ProtonPairsImageType::ImageDimension>
      (
      m_ProtonPairs->GetLargestPossibleRegion(),
      [this, iProj](const ProtonPairsImageType::RegionType & outputRegionForThread)
        {
        // Create MLP depending on type
        pct::MostLikelyPathFunction<double>::Pointer mlp;
        if(m_MostLikelyPathType == "polynomial")
          mlp = pct::ThirdOrderPolynomialMLPFunction<double>::New();
        else if (m_MostLikelyPathType == "schulte")
          mlp = pct::SchulteMLPFunction::New();
        else
          {
          itkGenericExceptionMacro("MLP must either be schulte or polynomial, not [" << m_MostLikelyPathType << ']');
          }

        // Get geometry information. We need the rotation matrix alone to rotate the direction
        // and the same matrix combined with mm (physical point) to voxel conversion for the volume.
        GeometryType::ThreeDHomogeneousMatrixType rotMat, volPPToIndex;

        rotMat = m_Geometry->GetRotationMatrices()[iProj].GetInverse();

        itk::Matrix<double, 5, 5> volPPToIndex44;
        volPPToIndex44 = rtk::GetPhysicalPointToIndexMatrix( this->GetInput() );
        for(int j=0; j<3; j++)
          {
          for(int i=0; i<3; i++)
            volPPToIndex[j][i] = volPPToIndex44[j][i];
          volPPToIndex[j][3] = volPPToIndex44[j][4];
          }
        volPPToIndex[3][3] = 1.;

        GeometryType::ThreeDHomogeneousMatrixType rotAndVoxConvMat;
        rotAndVoxConvMat = volPPToIndex.GetVnlMatrix() * rotMat.GetVnlMatrix();

        // Image information constants
        const typename OutputImageType::SizeType    imgSize    = this->GetInput()->GetBufferedRegion().GetSize();
        const typename OutputImageType::SpacingType imgSpacing = this->GetInput()->GetSpacing();
        const double halfSpacingRadian = imgSpacing[3]*itk::Math::pi/360.;

        typename OutputImageType::PixelType *imgData = this->GetOutput()->GetBufferPointer();
        unsigned int *imgCountData = m_Counts->GetBufferPointer();
        itk::Vector<float, 3> imgSpacingInv;
        double minSpacing = imgSpacing[0];
        for(unsigned int i=0; i<3; i++)
          {
          imgSpacingInv[i] = 1./imgSpacing[i];
          minSpacing = std::min(imgSpacing[i], minSpacing);
          }

        // Create matrices to cancel rotation
        std::vector< GeometryType::ThreeDHomogeneousMatrixType > noRotationMatrices;
        for(unsigned int i = 0; i<imgSize[3]; i++)
          {
          GeometryType::ThreeDHomogeneousMatrixType m;
          m =  m_Geometry->ComputeRotationHomogeneousMatrix(0, -1.*i*itk::Math::pi/imgSize[3], 0);
          m =  rotAndVoxConvMat.GetVnlMatrix() * m.GetVnlMatrix();
          noRotationMatrices.push_back( m );
          }

        // Corrections
        typedef itk::Vector<double,3> VectorType;

        // Calculate corner positions and largest diagonal in axial plane
        typename TOutputImage::IndexType idxCorner1, idxCorner2;
        idxCorner1 = this->GetInput()->GetLargestPossibleRegion().GetIndex();
        idxCorner2 = idxCorner1 + this->GetInput()->GetLargestPossibleRegion().GetSize();
        for(unsigned int i=0; i<TOutputImage::ImageDimension; i++)
          idxCorner1[i] -= 1;
        typename TOutputImage::PointType corner1, corner2;
        this->GetInput()->TransformIndexToPhysicalPoint(idxCorner1, corner1);
        this->GetInput()->TransformIndexToPhysicalPoint(idxCorner2, corner2);
        const double cornerMaxX = std::max(std::fabs(corner1[0]), std::fabs(corner2[0]));
        const double cornerMaxZ = std::max(std::fabs(corner1[2]), std::fabs(corner2[2]));
        const double largestDiagonal = sqrt(cornerMaxX*cornerMaxX + cornerMaxZ*cornerMaxZ);

        // Create zmm lut (look up table)
        const double zPlaneInInMM = -1.*largestDiagonal;
        const double zPlaneOutInMM = largestDiagonal;
        if(zPlaneInInMM > zPlaneOutInMM)
          {
          itkGenericExceptionMacro("Required condition pIn[2] > pOut[2] is not met, check coordinate system.");
          }
        std::vector<double> zmm;
        zmm.push_back(zPlaneInInMM);
        while(zmm.back()+minSpacing<zPlaneOutInMM)
          zmm.push_back(zmm.back()+minSpacing);

        // Process pairs
        itk::ImageRegionIterator<ProtonPairsImageType> it(m_ProtonPairs, outputRegionForThread);
        while(!it.IsAtEnd())
          {
          if(outputRegionForThread.GetIndex(1)==0 && it.GetIndex()[1]%1000==0)
            {
            std::cout << '\r'
                      << "Pair file #"
                      << iProj+1
                      << " out of "
                      << m_ProtonPairsFileNames.size()
                      << ", "
                      << it.GetIndex()[1] << " pairs of protons processed ("
                      << 100*it.GetIndex()[1]/outputRegionForThread.GetSize(1) << "%) in thread 1"
                      << std::flush;
            }

          VectorType pIn = it.Get();
          ++it;
          VectorType pOut = it.Get();
          ++it;
          VectorType dIn = it.Get();
          ++it;
          VectorType dOut = it.Get();
          ++it;

          // Line integral
          const double eIn = it.Get()[0];
          const double eOut = it.Get()[1];
          double value = 0.;
          if(eIn==0.)
            value = eOut; // Directly read WEPL
          else
            value = m_ConvFunc->GetValue(eOut, eIn); // convert to WEPL
          ++it;

          // Move straight to entrance and exit shapes
          VectorType pSIn = pIn;
          VectorType pSOut = pOut;
          double nearDistIn, nearDistOut, farDistIn, farDistOut;
          if(m_QuadricIn.GetPointer()!=NULL)
            {
            if(m_QuadricIn->IsIntersectedByRay(pIn,dIn,nearDistIn,farDistIn) &&
               m_QuadricOut->IsIntersectedByRay(pOut,dOut,nearDistOut,farDistOut))
              {
              pSIn  = pIn  + dIn  * nearDistIn;
              if(pSIn[2]<pIn[2]  || pSIn[2]>pOut[2])
                pSIn  = pIn  + dIn  * farDistIn;
              pSOut = pOut + dOut * nearDistOut;
              if(pSOut[2]<pIn[2] || pSOut[2]>pOut[2])
                pSOut = pOut + dOut * farDistOut;
              }
            }

          // Normalize direction with respect to z
          dIn[0] /= dIn[2];
          dIn[1] /= dIn[2];
          //dIn[2] = 1.; SR: implicit in the following
          dOut[0] /= dOut[2];
          dOut[1] /= dOut[2];
          //dOut[2] = 1.; SR: implicit in the following

          VectorType dCurr = dIn;
          VectorType pCurr(0.);

          // Init MLP before mm to voxel conversion
          mlp->Init(pSIn, pSOut, dIn, dOut);

          for(unsigned int k=0; k<zmm.size(); k++)
            {
            pCurr[2] = zmm[k];
            if(pCurr[2]<=pSIn[2]) //before entrance
              {
              const double z = (pCurr[2]-pIn[2]);
              pCurr[0] = pIn[0]+z*dIn[0];
              pCurr[1] = pIn[1]+z*dIn[1];
              dCurr = dIn;
              }
            else if(pCurr[2]>=pSOut[2]) //after exit
              {
              const double z = (pCurr[2]-pSOut[2]);
              pCurr[0] = pSOut[0]+z*dOut[0];
              pCurr[1] = pSOut[1]+z*dOut[1];
              dCurr = dOut;
              }
            else //MLP
              {
              dCurr[0] = pCurr[0];
              dCurr[1] = pCurr[1];
              dCurr[2] = minSpacing;
              mlp->Evaluate(zmm[k], pCurr[0], pCurr[1]);
              dCurr[0] = pCurr[0] - dCurr[0];
              dCurr[1] = pCurr[1] - dCurr[1];
              }
            dCurr[2] *= -1.;
            pCurr[2] *= -1.;

            // estimate direction
            VectorType dCurrRot(0.);
            for(unsigned int i=0; i<3; i++)
              {
              for(unsigned int j=0; j<3; j++)
                {
                dCurrRot[i] += rotMat[i][j] * dCurr[j];
                }
              }
            double theta = std::atan( dCurrRot[0] / dCurrRot[2] ) + halfSpacingRadian;
	    theta /= itk::Math::pi;     // Convert to half turns
	    theta -= std::floor(theta); // Between 0 and 1 (i.e., 0 and pi)
            theta *= 180.;              // To degrees
            // To pixel index
            theta -= this->GetOutput()->GetOrigin()[3];
            theta /= this->GetOutput()->GetSpacing()[3];
            typename OutputImageType::IndexType idx;
            idx[3] = itk::Math::Floor<int>(theta);
            if( idx[3] < 0 || idx[3] >= (int)imgSize[3] )
              continue;

            // rotation + mm to voxel conversion
            VectorType pCurrRot;
            if(m_DisableRotation)
              {
              for(unsigned int i=0; i<3; i++)
                {
                pCurrRot[i] = noRotationMatrices[idx[3]][i][3];
                for(unsigned int j=0; j<3; j++)
                  pCurrRot[i] += noRotationMatrices[idx[3]][i][j] * pCurr[j];
                }
              }
            else
              {
              for(unsigned int i=0; i<3; i++)
                {
                pCurrRot[i] = rotAndVoxConvMat[i][3];
                for(unsigned int j=0; j<3; j++)
                  pCurrRot[i] += rotAndVoxConvMat[i][j] * pCurr[j];
                }
              }

            for(int i=0; i<3; i++)
              idx[i] = itk::Math::Round<int,double>(pCurrRot[i]);
            if(idx[0]>=0 && idx[0]<(int)imgSize[0] &&
               idx[1]>=0 && idx[1]<(int)imgSize[1] &&
               idx[2]>=0 && idx[2]<(int)imgSize[2])
              {
              typename OutputImageType::OffsetValueType offset = this->GetOutput()->ComputeOffset(idx);
              m_Mutex.lock();
              imgData[ offset ] += value;
              imgCountData[ offset ]++;
              m_Mutex.unlock();
              }
            }
          }
        },
      nullptr);
      std::cout << '\r'
                << "Pair file #"
                << iProj+1
                << " out of "
                << m_ProtonPairsFileNames.size()
                << ", "
                << m_ProtonPairs->GetLargestPossibleRegion().GetSize(1) << " pairs of protons processed (100%) in all threads."
                << std::endl;
#ifdef MLP_TIMING
      mlp->PrintTiming(std::cout);
#endif
    }
  this->AfterThreadedGenerateData();
}

template <class TInputImage, class TOutputImage>
void
ProtonPairsToBackProjection<TInputImage, TOutputImage>
::AfterThreadedGenerateData()
{
  typedef typename itk::ImageRegionIterator<TOutputImage> ImageIteratorType;
  ImageIteratorType itOut(this->GetOutput(), this->GetOutput()->GetRequestedRegion());

  typedef itk::ImageRegionIterator<CountImageType> ImageCountIteratorType;
  ImageCountIteratorType itCOut(m_Counts, this->GetOutput()->GetRequestedRegion());

  // Set count image information
  m_Counts->SetSpacing( this->GetOutput()->GetSpacing() );
  m_Counts->SetOrigin( this->GetOutput()->GetOrigin() );

  // Normalize with proton count (average)
  while(!itCOut.IsAtEnd())
    {
    if(itCOut.Get())
      itOut.Set(itOut.Get()/itCOut.Get());
    ++itOut;
    ++itCOut;
    }
}

}
