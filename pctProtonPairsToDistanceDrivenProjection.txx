#include <itkImageFileReader.h>
#include <itkImageRegionIterator.h>

#include "pctThirdOrderPolynomialMLPFunction.h"
#include "pctSchulteMLPFunction.h"
#include "pctEnergyStragglingFunctor.h"
#include "pctScatteringWEPLFunction.h"

namespace pct
{

template <class TInputImage, class TOutputImage>
void
ProtonPairsToDistanceDrivenProjection<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  m_Outputs.resize( this->GetNumberOfThreads() );
  m_Counts.resize( this->GetNumberOfThreads() );
  m_Angles.resize( this->GetNumberOfThreads() );
  m_AnglesVectors.resize( this->GetInput()->GetLargestPossibleRegion().GetNumberOfPixels() );
  m_AnglesSq.resize( this->GetNumberOfThreads() );

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
  if(m_MostLikelyPathType == "polynomial")
    mlp = pct::ThirdOrderPolynomialMLPFunction<double>::New();
  else if (m_MostLikelyPathType == "schulte")
    mlp = pct::SchulteMLPFunction::New();
  else
    {
    itkGenericExceptionMacro("MLP must either be schulte or polynomial, not [" << m_MostLikelyPathType << ']');
    }

  // Create thread image and corresponding stack to count events
  m_Counts[threadId] = CountImageType::New();
  m_Counts[threadId]->SetRegions(this->GetInput()->GetLargestPossibleRegion());
  m_Counts[threadId]->Allocate();
  m_Counts[threadId]->FillBuffer(0);

  if(!m_Robust || threadId==0)
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
    m_Angle = m_Angles[0];
    m_AngleSq = m_AnglesSq[0];
    }
  else
    {
    m_Outputs[threadId] = OutputImageType::New();
    m_Outputs[threadId]->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_Outputs[threadId]->Allocate();
    }
  m_Outputs[threadId]->FillBuffer(0.);

  size_t nprotons = m_ProtonPairs->GetLargestPossibleRegion().GetSize()[1];
  ProtonPairsImageType::RegionType region = m_ProtonPairs->GetLargestPossibleRegion();
  region.SetIndex(1, threadId*nprotons/this->GetNumberOfThreads());
  region.SetSize(1, vnl_math_min((unsigned long)nprotons/this->GetNumberOfThreads(), nprotons-region.GetIndex(1)));

  // Image information constants
  const typename OutputImageType::SizeType    imgSize    = this->GetInput()->GetBufferedRegion().GetSize();
  const typename OutputImageType::PointType   imgOrigin  = this->GetInput()->GetOrigin();
  const typename OutputImageType::SpacingType imgSpacing = this->GetInput()->GetSpacing();
  const unsigned long npixelsPerSlice = imgSize[0] * imgSize[1];

  typename OutputImageType::PixelType *imgData = m_Outputs[threadId]->GetBufferPointer();
  float *imgCountData = m_Counts[threadId]->GetBufferPointer();
  float *imgAngleData = NULL, *imgAngleSqData = NULL;
  if(!m_Robust)
    {
    imgAngleData = m_Angles[threadId]->GetBufferPointer();
    imgAngleSqData = m_AnglesSq[threadId]->GetBufferPointer();
    }

  itk::Vector<float, 3> imgSpacingInv;
  for(unsigned int i=0; i<3; i++)
    imgSpacingInv[i] = 1./imgSpacing[i];

  // Corrections
  typedef itk::Vector<double,3> VectorType;

  // Create a local copy of quadrics (surface object) for multithreading
  RQIType::Pointer quadricIn, quadricOut;
  if(m_QuadricIn.GetPointer()!=NULL)
    {
    quadricIn = RQIType::New();
    quadricIn->SetA(m_QuadricIn->GetA());
    quadricIn->SetB(m_QuadricIn->GetB());
    quadricIn->SetC(m_QuadricIn->GetC());
    quadricIn->SetD(m_QuadricIn->GetD());
    quadricIn->SetE(m_QuadricIn->GetE());
    quadricIn->SetF(m_QuadricIn->GetF());
    quadricIn->SetG(m_QuadricIn->GetG());
    quadricIn->SetH(m_QuadricIn->GetH());
    quadricIn->SetI(m_QuadricIn->GetI());
    quadricIn->SetJ(m_QuadricIn->GetJ());

    quadricOut = RQIType::New();
    quadricOut->SetA(m_QuadricOut->GetA());
    quadricOut->SetB(m_QuadricOut->GetB());
    quadricOut->SetC(m_QuadricOut->GetC());
    quadricOut->SetD(m_QuadricOut->GetD());
    quadricOut->SetE(m_QuadricOut->GetE());
    quadricOut->SetF(m_QuadricOut->GetF());
    quadricOut->SetG(m_QuadricOut->GetG());
    quadricOut->SetH(m_QuadricOut->GetH());
    quadricOut->SetI(m_QuadricOut->GetI());
    quadricOut->SetJ(m_QuadricOut->GetJ());
    }

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
    const double angle_out = it.Get()[1];

    VectorType dOut = it.Get();
    ++it;

    //typedef itk::Vector<double, 2> VectorTwoDType;
    //VectorTwoDType dInX, dOutX;
    //dInX[0] = dIn[0];
    //dInX[1] = dIn[2];
    //dOutX[0] = dOut[0];
    //dOutX[1] = dOut[2];
    //const double angle_out = vcl_acos( std::min(1.,dInX*dOutX / ( dIn.GetNorm() * dOutX.GetNorm() ) ) );

    if(pIn[2] > pOut[2])
      {
      itkGenericExceptionMacro("Required condition pIn[2] > pOut[2] is not met, check coordinate system.");
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

    VectorType nucInfo = it.Get();
    ++it;

    // Move straight to entrance and exit shapes
    VectorType pSIn  = pIn;
    VectorType pSOut = pOut;
    if(quadricIn.GetPointer()!=NULL)
      {
      quadricIn->SetRayOrigin(pIn);
      quadricOut->SetRayOrigin(pOut);
      if(quadricIn->Evaluate(dIn) && quadricOut->Evaluate(dOut))
        {
        pSIn  = pIn  + dIn  * quadricIn ->GetNearestDistance();
        if(pSIn[2]<pIn[2]  || pSIn[2]>pOut[2])
          pSIn  = pIn  + dIn  * quadricIn ->GetFarthestDistance();
        pSOut = pOut + dOut * quadricOut->GetNearestDistance();
        if(pSOut[2]<pIn[2] || pSOut[2]>pOut[2])
          pSOut = pOut + dOut * quadricOut->GetFarthestDistance();
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

    for(unsigned int k=0; k<imgSize[2]; k++)
      {
      double xx, yy;
      const double dk = zmm[k];
      if(dk<=pSIn[2]) //before entrance
        {
        const double z = (dk-pIn[2]);
        xx = pIn[0]+z*dIn[0];
        yy = pIn[1]+z*dIn[1];
        }
      else if(dk>=pSOut[2]) //after exit
        {
        const double z = (dk-pSOut[2]);
        xx = pSOut[0]+z*dOut[0];
        yy = pSOut[1]+z*dOut[1];
        }
      else //MLP
        {
        mlp->Evaluate(zmm[k], xx, yy);
        }

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
        if(m_Robust)
          {
          m_AnglesVectorsMutex.Lock();
          m_AnglesVectors[idx].push_back(vcl_abs(angle_out));
          m_AnglesVectorsMutex.Unlock();
          }
        else
          {
          imgAngleData[ idx ] += angle_out;
          imgAngleSqData[ idx ] += angle_out*angle_out;
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

  typedef itk::ImageRegionIterator<AngleImageType> ImageAngleIteratorType;
  ImageAngleIteratorType itAngleOut(m_Angles[0], m_Outputs[0]->GetLargestPossibleRegion());

  typedef itk::ImageRegionIterator<AngleImageType> ImageAngleSqIteratorType;
  ImageAngleSqIteratorType itAngleSqOut(m_AnglesSq[0], m_Outputs[0]->GetLargestPossibleRegion());

  // Merge the projection computed in each thread to the first one
  for(unsigned int i=1; i<this->GetNumberOfThreads(); i++)
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

  if(!m_Robust)
    {
    for(unsigned int i=1; i<this->GetNumberOfThreads(); i++)
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
  // Set count image information
  m_Count->SetSpacing( this->GetOutput()->GetSpacing() );
  m_Count->SetOrigin( this->GetOutput()->GetOrigin() );

  // Set scattering wepl image information
  m_Angle->SetSpacing( this->GetOutput()->GetSpacing() );
  m_Angle->SetOrigin( this->GetOutput()->GetOrigin() );

  // Initialize scattering WEPL LUT
  pct::Functor::ScatteringWEPL::ScatteringLUT::SetG4LUT(m_ProtonPairs->GetBufferPointer()[4][0]);

  // Pointer
  pct::Functor::ScatteringWEPL::ConvertToScatteringWEPL pointer;

  std::vector< std::vector<float> >::iterator itAnglesVectors = m_AnglesVectors.begin();
  while(!itCOut.IsAtEnd())
    {
    if(itCOut.Get())
      {
      // Normalize eloss wepl with proton count (average)
      itOut.Set(itOut.Get()/itCOut.Get());

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
          itAngleOut.Set( pointer.GetValue(sigma * sigma) );
          }
        }
      else
        {
        double sigma2 = itAngleSqOut.Get()/itCOut.Get() - itAngleOut.Get()*itAngleOut.Get()/itCOut.Get()/itCOut.Get() ;
        itAngleOut.Set( pointer.GetValue(sigma2) );
        }
      }

    ++itOut;
    ++itCOut;
    ++itAngleOut;
    ++itAngleSqOut;
    ++itAnglesVectors;
    }

  // Free images created in threads
  m_Outputs.resize( 0 );
  m_Counts.resize( 0 );
  m_Angles.resize( 0 );
  m_AnglesSq.resize( 0 );
  m_AnglesVectors.resize( 0 );
}

}
