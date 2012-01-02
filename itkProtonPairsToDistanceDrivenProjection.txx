#include <itkImageFileReader.h>
#include <itkImageRegionIterator.h>

#include "pctBetheBlochFunctor.h"
//#include "itkThirdOrderPolynomialMLPFunction.h"
#include "itkSchulteMLPFunction.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
void
ProtonPairsToDistanceDrivenProjection<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  m_Outputs.resize( this->GetNumberOfThreads() );
  m_Counts.resize( this->GetNumberOfThreads() );
  if(m_QuadricOut.GetPointer()==NULL)
    m_QuadricOut = m_QuadricIn;
}

template <class TInputImage, class TOutputImage>
void
ProtonPairsToDistanceDrivenProjection<TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId )
{
  // Create thread image and corresponding stack to count events
  m_Counts[threadId] = CountImageType::New();
  m_Counts[threadId]->SetRegions(this->GetInput()->GetLargestPossibleRegion());
  m_Counts[threadId]->Allocate();
  m_Counts[threadId]->FillBuffer(0);

  if(threadId==0)
    {
    m_Outputs[0] = this->GetOutput();
    m_Count = m_Counts[0];
    }
  else
    {
    m_Outputs[threadId] = OutputImageType::New();
    m_Outputs[threadId]->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_Outputs[threadId]->Allocate();
    }
  m_Outputs[threadId]->FillBuffer(0.);

  // Read pairs
  typedef itk::ImageFileReader< ProtonPairsImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( m_ProtonPairsFileName );
  reader->UpdateOutputInformation();
  size_t nprotons = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
  ProtonPairsImageType::RegionType region = reader->GetOutput()->GetLargestPossibleRegion();
  region.SetIndex(1, threadId*nprotons/this->GetNumberOfThreads());
  region.SetSize(1, vnl_math_min((unsigned long)nprotons/this->GetNumberOfThreads(), nprotons-region.GetIndex(1)));
  reader->GetOutput()->SetRequestedRegion(region);
  reader->Update();

  // Image information constants
  typename OutputImageType::SizeType    imgSize    = this->GetInput()->GetBufferedRegion().GetSize();
  typename OutputImageType::PointType   imgOrigin  = this->GetInput()->GetOrigin();
  typename OutputImageType::SpacingType imgSpacing = this->GetInput()->GetSpacing();
  unsigned long npixelsPerSlice = imgSize[0] * imgSize[1];

  typename OutputImageType::PixelType *imgData = m_Outputs[threadId]->GetBufferPointer();
  unsigned int *imgCountData = m_Counts[threadId]->GetBufferPointer();
  itk::Vector<float, 3> imgSpacingInv;
  typename OutputImageType::PointType originInVox;
  for(unsigned int i=0; i<3; i++)
    {
    imgSpacingInv[i] = 1./imgSpacing[i];
    originInVox[i] = -imgOrigin[i]*imgSpacingInv[i];
    }

  pct::Functor::IntegratedBetheBlochProtonStoppingPowerInverse<float, double> convFunc;

  // Corrections
  typedef itk::Vector<double,3> VectorType;
  VectorType source;
  source.Fill(0.);
  source[2] = m_SourceDistance;

  // Create a local copy of quadrics for multithreading
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

  // Create zmm and magnitude lut
  std::vector<double> zmm(imgSize[2]);
  std::vector<double> zmag(imgSize[2]);
  const double sourcePosInVox = (m_SourceDistance-imgOrigin[2]) * imgSpacingInv[2];
  const double zPlaneOutInMM = (*(reader->GetOutput()->GetBufferPointer()+1))[2];

  const double zPlaneOutInVox = (zPlaneOutInMM-imgOrigin[2]) * imgSpacingInv[2];
  for(unsigned int i=0; i<imgSize[2]; i++)
    {
    zmm[i] = i*imgSpacingInv[2]+imgOrigin[2];
    zmag[i] = (zPlaneOutInVox-sourcePosInVox)/(i-sourcePosInVox);
    }

  // Process pairs
  ImageRegionIterator<ProtonPairsImageType> it(reader->GetOutput(), region);
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
    const double value = convFunc.GetValue(it.Get()[1], it.Get()[0]);
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
        if(pSIn[2]<pIn[2])
          pSIn  = pIn  + dIn  * quadricIn ->GetFarthestDistance();
        pSOut = pOut + dOut * quadricOut->GetFarthestDistance();
        }
      }

    // Convert everything to voxel coordinates
    for(unsigned int i=0; i<3; i++)
      {
      pIn[i]   = (pIn[i]   - imgOrigin[i]) * imgSpacingInv[i];
      pOut[i]  = (pOut[i]  - imgOrigin[i]) * imgSpacingInv[i];
      dIn[i]   = dIn[i]  * imgSpacingInv[i];
      dOut[i]  = dOut[i] * imgSpacingInv[i];
      }
    for(unsigned int i=0; i<2; i++)
      {
      // Only first 2 coordinates, absolute depth is essential for MLP
      pSIn[i]  = (pSIn[i]  - imgOrigin[i]) * imgSpacingInv[i];
      pSOut[i] = (pSOut[i] - imgOrigin[i]) * imgSpacingInv[i];
      }

    // Normalize direction with respect to z
    dIn[0] /= dIn[2];
    dIn[1] /= dIn[2];
    //dIn[2] = 1.; SR: implicit in the following
    dOut[0] /= dOut[2];
    dOut[1] /= dOut[2];
    //dOut[2] = 1.; SR: implicit in the following

    // Init MLP before mm to voxel conversion of depth
    //itk::ThirdOrderPolynomialMLPFunction<double>::Pointer mlp;
    //mlp = itk::ThirdOrderPolynomialMLPFunction<double>::New();
    itk::SchulteMLPFunction::Pointer mlp;
    mlp = itk::SchulteMLPFunction::New();
    mlp->Init(pSIn, pSOut, dIn, dOut);

    // Finalize mm to voxel conversion
    pSIn[2]  = (pSIn[2]  - imgOrigin[2]) * imgSpacingInv[2];
    pSOut[2] = (pSOut[2] - imgOrigin[2]) * imgSpacingInv[2];

    for(unsigned int k=0; k<imgSize[2]; k+=1)
      {
      double xx, yy;
      const double dk = k; //SR make sure to convert it once only
      if(dk<=pSIn[2])
        {
        const double z = (dk-pIn[2]);
        xx = pIn[0]+z*dIn[0];
        yy = pIn[1]+z*dIn[1];
        }
      else if(dk>=pSOut[2])
        {
        const double z = (dk-pSOut[2]);
        xx = pSOut[0]+z*dOut[0];
        yy = pSOut[1]+z*dOut[1];
        }
      else
        {
        mlp->Evaluate(zmm[k], xx, yy);
        }

      xx = (xx-originInVox[0])*zmag[k]+originInVox[0];
      const int i = int(xx+0.5);
      if(i<0 || i>=(int)imgSize[0])
        continue;

      yy = (yy-originInVox[1])*zmag[k]+originInVox[1];
      const int j = int(yy+0.5);
      if(j<0 || j>=(int)imgSize[1])
        continue;

      unsigned long idx = i+j*imgSize[0]+k*npixelsPerSlice;
      imgData[ idx ] += value;
      imgCountData[ idx ]++;
      }
  }
  if(threadId==0)
    std::cout << std::endl << "Done!" << std::endl;
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
  for(int i=1; i<this->GetNumberOfThreads(); i++)
    {
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

  // Normalize with proton count (average)
  while(!itCOut.IsAtEnd())
    {
    if(itCOut.Get())
      itOut.Set(itOut.Get()/itCOut.Get());
    ++itOut;
    ++itCOut;
    }

  // Free images created in threads
  m_Outputs.resize( 0 );
  m_Counts.resize( 0 );
}

}
