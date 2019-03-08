#include <itkImageFileReader.h>
#include <itkImageRegionIterator.h>

#include "pctThirdOrderPolynomialMLPFunction.h"
#include "pctSchulteMLPFunction.h"
#include "pctPolynomialMLPFunction.h"
#include "pctEnergyStragglingFunctor.h"

namespace pct
{

template <class TInputImage, class TOutputImage>
ProtonPairsToDistanceDrivenProjection<TInputImage, TOutputImage>
::ProtonPairsToDistanceDrivenProjection():m_Robust(false),m_ComputeScattering(false)
{
  this->DynamicMultiThreadingOff();
  this->SetNumberOfWorkUnits( itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads() );
}

template <class TInputImage, class TOutputImage>
void
ProtonPairsToDistanceDrivenProjection<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  m_Outputs.resize( this->GetNumberOfWorkUnits() );
  m_Counts.resize( this->GetNumberOfWorkUnits() );
  if(m_ComputeScattering)
    {
    m_Angles.resize( this->GetNumberOfWorkUnits() );
    m_AnglesVectors.resize( this->GetInput()->GetLargestPossibleRegion().GetNumberOfPixels() );
    m_AnglesSq.resize( this->GetNumberOfWorkUnits() );
    }

  if(m_QuadricOut.GetPointer()==NULL)
    m_QuadricOut = m_QuadricIn;
  m_ConvFunc = new Functor::IntegratedBetheBlochProtonStoppingPowerInverse<float, double>(m_IonizationPotential, 600.*CLHEP::MeV, 0.1*CLHEP::keV);

  // Read pairs
  typedef itk::ImageFileReader< ProtonPairsImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( m_ProtonPairsFileName );
  reader->Update();
  m_ProtonPairs = reader->GetOutput();
}

template <class TInputImage, class TOutputImage>
void
ProtonPairsToDistanceDrivenProjection<TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType& itkNotUsed(outputRegionForThread), rtk::ThreadIdType threadId)
{
  // Create MLP depending on type
  pct::MostLikelyPathFunction<double>::Pointer mlp;
  pct::PolynomialMLPFunction::Pointer mlp_poly; // NK: must declare outside of if clause because mlp_poly is needed later
  if(m_MostLikelyPathType == "polynomial")
    mlp = pct::ThirdOrderPolynomialMLPFunction<double>::New();
  else if(m_MostLikelyPathType == "krah")
    {
    mlp_poly = pct::PolynomialMLPFunction::New();
    mlp_poly->SetPolynomialDegree(m_MostLikelyPathPolynomialDegree);
    mlp = mlp_poly;
    }
  else if(m_MostLikelyPathType == "schulte")
    mlp = pct::SchulteMLPFunction::New();
  else
    {
    itkGenericExceptionMacro("MLP must either be schulte, polynomial, or krah, not [" << m_MostLikelyPathType << ']');
    }

  // Create thread image and corresponding stack to count events
  m_Counts[threadId] = CountImageType::New();
  m_Counts[threadId]->SetRegions(this->GetInput()->GetLargestPossibleRegion());
  m_Counts[threadId]->Allocate();
  m_Counts[threadId]->FillBuffer(0);

  if( m_ComputeScattering && (!m_Robust || threadId==0) )
    {
    m_Angles[threadId] = AngleImageType::New();
    m_Angles[threadId]->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_Angles[threadId]->Allocate();
    m_Angles[threadId]->FillBuffer(0);

    m_AnglesSq[threadId] = AngleImageType::New();
    m_AnglesSq[threadId]->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_AnglesSq[threadId]->Allocate();
    m_AnglesSq[threadId]->FillBuffer(0);
    }

  if(threadId==0)
    {
    m_Outputs[0] = this->GetOutput();
    m_Count = m_Counts[0];
    if(m_ComputeScattering)
      {
      m_Angle = m_Angles[0];
      m_AngleSq = m_AnglesSq[0];
      }
    }
  else
    {
    m_Outputs[threadId] = OutputImageType::New();
    m_Outputs[threadId]->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_Outputs[threadId]->Allocate();
    }
  m_Outputs[threadId]->FillBuffer(0.);

  size_t nprotons = m_ProtonPairs->GetLargestPossibleRegion().GetSize()[1];
  size_t nprotonsPerThread = nprotons/this->GetMultiThreader()->GetNumberOfWorkUnits();
  ProtonPairsImageType::RegionType region = m_ProtonPairs->GetLargestPossibleRegion();
  region.SetIndex(1, threadId*nprotonsPerThread);
  if(threadId == this->GetMultiThreader()->GetNumberOfWorkUnits()-1)
    region.SetSize(1, nprotons-region.GetIndex(1));
  else
    region.SetSize(1, nprotons/this->GetMultiThreader()->GetNumberOfWorkUnits());

  // Image information constants
  const typename OutputImageType::SizeType    imgSize    = this->GetInput()->GetBufferedRegion().GetSize();
  const typename OutputImageType::PointType   imgOrigin  = this->GetInput()->GetOrigin();
  const typename OutputImageType::SpacingType imgSpacing = this->GetInput()->GetSpacing();
  const unsigned long npixelsPerSlice = imgSize[0] * imgSize[1];

  typename OutputImageType::PixelType *imgData = m_Outputs[threadId]->GetBufferPointer();
  unsigned int *imgCountData = m_Counts[threadId]->GetBufferPointer();
  float *imgAngleData = NULL, *imgAngleSqData = NULL;
  if(m_ComputeScattering && !m_Robust)
    {
    imgAngleData = m_Angles[threadId]->GetBufferPointer();
    imgAngleSqData = m_AnglesSq[threadId]->GetBufferPointer();
    }

  itk::Vector<float, 3> imgSpacingInv;
  for(unsigned int i=0; i<3; i++)
    imgSpacingInv[i] = 1./imgSpacing[i];

  // Corrections
  typedef itk::Vector<double,3> VectorType;

  // Create zmm and magnitude lut (look up table)
  itk::ImageRegionIterator<ProtonPairsImageType> it(m_ProtonPairs, region);
  std::vector<double> zmm(imgSize[2]);
  std::vector<double> zmag(imgSize[2]);
  ++it;
  const double zPlaneOutInMM = it.Get()[2];
  --it;
  for(unsigned int i=0; i<imgSize[2]; i++)
    {
    zmm[i] = i*imgSpacing[2]+imgOrigin[2];
    zmag[i] = (m_SourceDistance==0.)?1:(zPlaneOutInMM-m_SourceDistance)/(zmm[i]-m_SourceDistance);
    }

  // Process pairs
  while(!it.IsAtEnd())
  {
    if(threadId==0 && it.GetIndex()[1]%10000==0)
      {
      std::cout << '\r'
                << it.GetIndex()[1] << " pairs of protons processed ("
                << 100*it.GetIndex()[1]/region.GetSize(1) << "%) in thread 1"
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

    double anglex = 0., angley = 0.;
    if( m_ComputeScattering )
      {
      typedef itk::Vector<double, 2> VectorTwoDType;

      VectorTwoDType dInX, dInY, dOutX, dOutY;
      dInX[0] = dIn[0];
      dInX[1] = dIn[2];
      dInY[0] = dIn[1];
      dInY[1] = dIn[2];
      dOutX[0] = dOut[0];
      dOutX[1] = dOut[2];
      dOutY[0] = dOut[1];
      dOutY[1] = dOut[2];

      angley = std::acos( std::min(1.,dInY*dOutY / ( dInY.GetNorm() * dOutY.GetNorm() ) ) );
      anglex = std::acos( std::min(1.,dInX*dOutX / ( dInX.GetNorm() * dOutX.GetNorm() ) ) );
      }

    if(pIn[2] > pOut[2])
      {
      itkGenericExceptionMacro("Required condition pIn[2] > pOut[2] is not met, check coordinate system.");
      }
    if(dIn[2] < 0.)
      {
      itkGenericExceptionMacro("The code assumes that protons move in positive z.");
      }

    const double eIn = it.Get()[0];
    const double eOut = it.Get()[1];
    double value = 0.;
    if(eIn==0.)
      value = eOut; // Directly read WEPL
    else
      {
      value = m_ConvFunc->GetValue(eOut, eIn); // convert to WEPL
      }
    ++it;

    VectorType nucInfo(0.);
    if(it.GetIndex()[0] != 0)
      {
      nucInfo = it.Get();
      ++it;
      }

    // Move straight to entrance and exit shapes
    VectorType pSIn  = pIn;
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

    // Init MLP before mm to voxel conversion
    mlp->Init(pSIn, pSOut, dIn, dOut);

      std::vector<double> zmmMLP;
      std::vector<unsigned int> kMLP;
      double xxArr[imgSize[2]], yyArr[imgSize[2]];

      // loop to populate a vector to be passed to Evaluate if MLP type is elgible
      for(unsigned int k=0; k<imgSize[2]; k++)
      {
        const double dk = zmm[k];
        if(dk<=pSIn[2]) //before entrance
          {
          const double z = (dk-pIn[2]);
          xxArr[k] = pIn[0]+z*dIn[0];
          yyArr[k] = pIn[1]+z*dIn[1];
          }
        else if(dk>=pSOut[2]) //after exit
          {
          const double z = (dk-pSOut[2]);
          xxArr[k] = pSOut[0]+z*dOut[0];
          yyArr[k] = pSOut[1]+z*dOut[1];
          }
        else
          {
            if(m_MostLikelyPathType == "krah") // maybe more flexible to use a m_CanBeVectorised boolean flag and set it when creating the mlp object
            {
              // stock in vector for later if vectorisable
              zmmMLP.push_back(dk);
              kMLP.push_back(k);
            }
            else
            {
              // evaluate directly if not vectorisable
              mlp->Evaluate(zmm[k], xxArr[k], yyArr[k]);
            }
          }
      }

      // call Evaluate with vector as argument and insert result into result array
      // would be prefeable to avoid the copying step and use the xxMLP and yyMLP vectors directly
      // but that requires the reste of the function further down to be restructured a bit
      if(m_MostLikelyPathType == "krah")
      {
        std::vector<double> xxMLP;
        std::vector<double> yyMLP;
        xxMLP.resize(zmmMLP.size());
        yyMLP.resize(zmmMLP.size());

        mlp->Evaluate(zmmMLP, xxMLP, yyMLP);
        for(std::vector<int>::size_type i = 0; i != kMLP.size(); i++)
          {
          xxArr[kMLP[i]] = xxMLP[i];
          yyArr[kMLP[i]] = yyMLP[i];
          }
      }

    for(unsigned int k=0; k<imgSize[2]; k++)
      {
      double xx, yy;

      xx = xxArr[k];
      yy = yyArr[k];

      // Source at (0,0,args_info.source_arg), mag then to voxel
      xx = (xx*zmag[k] - imgOrigin[0]) * imgSpacingInv[0];
      yy = (yy*zmag[k] - imgOrigin[1]) * imgSpacingInv[1];

      // Lattice conversion
      const int i = itk::Math::Round<int,double>(xx);
      const int j = itk::Math::Round<int,double>(yy);
      if(i>=0 && i<(int)imgSize[0] &&
         j>=0 && j<(int)imgSize[1])
        {
        const unsigned long idx = i+j*imgSize[0]+k*npixelsPerSlice;
        imgData[ idx ] += value;
        imgCountData[ idx ]++;
        if(m_ComputeScattering)
          {
          if(m_Robust)
            {
            m_AnglesVectorsMutex.lock();
            m_AnglesVectors[idx].push_back(anglex);
            m_AnglesVectors[idx].push_back(angley);
            m_AnglesVectorsMutex.unlock();
            }
          else
            {
            imgAngleData[ idx ] += anglex;
            imgAngleData[ idx ] += angley;
            imgAngleSqData[ idx ] += anglex*anglex;
            imgAngleSqData[ idx ] += angley*angley;
            }
          }
        }
      }
  }

  if(threadId==0)
    {
    std::cout << '\r'
              << region.GetSize(1) << " pairs of protons processed (100%) in thread 1"
              << std::endl;
#ifdef MLP_TIMING
    mlp->PrintTiming(std::cout);
#endif
    }

}

template <class TInputImage, class TOutputImage>
void
ProtonPairsToDistanceDrivenProjection<TInputImage, TOutputImage>
::AfterThreadedGenerateData()
{
  typedef typename itk::ImageRegionIterator<TOutputImage> ImageIteratorType;
  ImageIteratorType itOut(m_Outputs[0], m_Outputs[0]->GetLargestPossibleRegion());

  typedef itk::ImageRegionIterator<CountImageType> ImageCountIteratorType;
  ImageCountIteratorType itCOut(m_Counts[0], m_Outputs[0]->GetLargestPossibleRegion());

  // Merge the projection computed in each thread to the first one
  for(unsigned int i=1; i<this->GetNumberOfWorkUnits(); i++)
    {
    if(m_Outputs[i].GetPointer() == NULL)
      continue;
    ImageIteratorType itOutThread(m_Outputs[i], m_Outputs[i]->GetLargestPossibleRegion());
    ImageCountIteratorType itCOutThread(m_Counts[i], m_Outputs[i]->GetLargestPossibleRegion());

    while(!itOut.IsAtEnd())
      {
      itOut.Set(itOut.Get()+itOutThread.Get());
      ++itOutThread;
      ++itOut;

      itCOut.Set(itCOut.Get()+itCOutThread.Get());
      ++itCOutThread;
      ++itCOut;
      }

    itOut.GoToBegin();
    itCOut.GoToBegin();
    }

  // Set count image information
  m_Count->SetSpacing( this->GetOutput()->GetSpacing() );
  m_Count->SetOrigin( this->GetOutput()->GetOrigin() );

  // Normalize eloss wepl with proton count (average)
  while(!itCOut.IsAtEnd())
    {
    if(itCOut.Get())
      itOut.Set(itOut.Get()/itCOut.Get());
    ++itOut;
    ++itCOut;
    }

  if(m_ComputeScattering)
    {
    typedef itk::ImageRegionIterator<AngleImageType> ImageAngleIteratorType;
    ImageAngleIteratorType itAngleOut(m_Angles[0], m_Outputs[0]->GetLargestPossibleRegion());

    typedef itk::ImageRegionIterator<AngleImageType> ImageAngleSqIteratorType;
    ImageAngleSqIteratorType itAngleSqOut(m_AnglesSq[0], m_Outputs[0]->GetLargestPossibleRegion());

    if(!m_Robust)
      {
      for(unsigned int i=1; i<this->GetNumberOfWorkUnits(); i++)
        {
        if(m_Outputs[i].GetPointer() == NULL)
          continue;
        ImageAngleIteratorType itAngleOutThread(m_Angles[i], m_Outputs[i]->GetLargestPossibleRegion());
        ImageAngleSqIteratorType itAngleSqOutThread(m_AnglesSq[i], m_Outputs[i]->GetLargestPossibleRegion());

        while(!itAngleOut.IsAtEnd())
          {
          itAngleOut.Set(itAngleOut.Get()+itAngleOutThread.Get());
          ++itAngleOutThread;
          ++itAngleOut;

          itAngleSqOut.Set(itAngleSqOut.Get()+itAngleSqOutThread.Get());
          ++itAngleSqOutThread;
          ++itAngleSqOut;
          }
        itAngleOut.GoToBegin();
        itAngleSqOut.GoToBegin();
        }
      }

    // Set scattering wepl image information
    m_Angle->SetSpacing( this->GetOutput()->GetSpacing() );
    m_Angle->SetOrigin( this->GetOutput()->GetOrigin() );

    itCOut.GoToBegin();
    std::vector< std::vector<float> >::iterator itAnglesVectors = m_AnglesVectors.begin();
    while(!itCOut.IsAtEnd())
      {
      if(itCOut.Get())
        {
        // Calculate angular variance (sigma2) and convert to scattering wepl
        if(m_Robust)
          {
          if(itCOut.Get()==1)
            {
            itAngleOut.Set( 0. );
            }
          else
            {
            // Angle: 38.30% (0.5 sigma) with interpolation (median is 0. and we only have positive values
            double sigmaAPos = itAnglesVectors->size()*0.3830;
            unsigned int sigmaASupPos = itk::Math::Ceil<unsigned int, double>(sigmaAPos);
            std::partial_sort(itAnglesVectors->begin(),
                              itAnglesVectors->begin()+sigmaASupPos+1,
                              itAnglesVectors->end());
            double sigmaADiff = sigmaASupPos-sigmaAPos;
            double sigma = 2.*(*(itAnglesVectors->begin()+sigmaASupPos)*(1.-sigmaADiff)+
                               *(itAnglesVectors->begin()+sigmaASupPos-1)*sigmaADiff); //x2 to get 1sigma
            itAngleOut.Set( sigma * sigma );
            }
          }
        else
          {
          double sigma2 = itAngleSqOut.Get()/itCOut.Get()/2;
          itAngleOut.Set( sigma2 );
          }
        }

      ++itCOut;
      ++itAngleOut;
      ++itAngleSqOut;
      ++itAnglesVectors;
      }
    }

  // Free images created in threads
  m_Outputs.resize( 0 );
  m_Counts.resize( 0 );
  m_Angles.resize( 0 );
  m_AnglesSq.resize( 0 );
  m_AnglesVectors.resize( 0 );
}

}
