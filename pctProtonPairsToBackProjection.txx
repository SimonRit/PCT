#include <itkImageFileReader.h>
#include <itkImageRegionIterator.h>

#include "pctThirdOrderPolynomialMLPFunction.h"
#include "pctSchulteMLPFunction.h"
#include "pctEnergyStragglingFunctor.h"

namespace pct
{

template <class TInputImage, class TOutputImage>
void
ProtonPairsToBackProjection<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  m_Outputs.resize( this->GetNumberOfThreads() );
  m_Counts.resize( this->GetNumberOfThreads() );
  if(m_QuadricOut.GetPointer()==NULL)
    m_QuadricOut = m_QuadricIn;
  m_ConvFunc = new Functor::IntegratedBetheBlochProtonStoppingPowerInverse<float, double>(m_IonizationPotential, 600.*CLHEP::MeV, 0.1*CLHEP::keV);

  for(unsigned int i=0; i<3; i++)
    {
    m_Barriers[i]= itk::Barrier::New();
    m_Barriers[i]->Initialize( this->GetNumberOfThreads() );
    }
}

template <class TInputImage, class TOutputImage>
void
ProtonPairsToBackProjection<TInputImage, TOutputImage>
::ThreadedGenerateData( const OutputImageRegionType& itkNotUsed(outputRegionForThread), rtk::ThreadIdType threadId)
{
  // Create thread image and corresponding stack to count events
  if(threadId==0)
    {
    m_Outputs[0] = this->GetOutput();
    m_Counts[0]  = this->GetCount();
    }
  else
    {
    m_Outputs[threadId] = OutputImageType::New();
    m_Outputs[threadId]->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_Outputs[threadId]->Allocate();
    m_Outputs[threadId]->FillBuffer(0.);

    m_Counts[threadId] = CountImageType::New();
    m_Counts[threadId]->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_Counts[threadId]->Allocate();
    m_Counts[threadId]->FillBuffer(0);
    }

  for(unsigned int iProj = 0; iProj<m_ProtonPairsFileNames.size(); iProj++)
    {
    m_Barriers[0]->Wait();
    if(threadId==0)
      {
      m_Barriers[2]= itk::Barrier::New();
      m_Barriers[2]->Initialize( this->GetNumberOfThreads() );
      }
    m_Barriers[1]->Wait();
    if(threadId==0)
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
      }
    if(threadId==0)
      {
      m_Barriers[0]= itk::Barrier::New();
      m_Barriers[0]->Initialize( this->GetNumberOfThreads() );
      }
    m_Barriers[2]->Wait();
    if(threadId==0)
      {
      m_Barriers[1]= itk::Barrier::New();
      m_Barriers[1]->Initialize( this->GetNumberOfThreads() );
      }

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

    size_t nprotons = m_ProtonPairs->GetLargestPossibleRegion().GetSize()[1];
    ProtonPairsImageType::RegionType region = m_ProtonPairs->GetLargestPossibleRegion();
    region.SetIndex(1, threadId*nprotons/this->GetNumberOfThreads());
    region.SetSize(1, vnl_math_min((unsigned long)nprotons/this->GetNumberOfThreads(), nprotons-region.GetIndex(1)));

    // Image information constants
    const typename OutputImageType::SizeType    imgSize    = this->GetInput()->GetBufferedRegion().GetSize();
    const typename OutputImageType::PointType   imgOrigin  = this->GetInput()->GetOrigin();
    const typename OutputImageType::SpacingType imgSpacing = this->GetInput()->GetSpacing();

    typename OutputImageType::PixelType *imgData = m_Outputs[threadId]->GetBufferPointer();
    unsigned int *imgCountData = m_Counts[threadId]->GetBufferPointer();
    itk::Vector<float, 3> imgSpacingInv;
    for(unsigned int i=0; i<3; i++)
      imgSpacingInv[i] = 1./imgSpacing[i];

    // Corrections
    typedef itk::Vector<double,3> VectorType;
    typedef itk::Vector<double,4> VectorHType;

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

    // Create zmm lut (look up table)
    itk::ImageRegionIterator<ProtonPairsImageType> it(m_ProtonPairs, region);
    std::vector<double> zmm[3];
    for(unsigned int d=0; d<3; d++)
      {
      zmm[d].resize(imgSize[d]);
      for(unsigned int i=0; i<imgSize[d]; i++)
        {
        zmm[d][i] = i*imgSpacing[d]+imgOrigin[d];
        }
      }

    // Process pairs
    while(!it.IsAtEnd())
      {
      if(threadId==0 && it.GetIndex()[1]%10000==0)
        {
        std::cout << '\r'
                  << "Pair file #"
                  << iProj+1
                  << " out of "
                  << m_ProtonPairsFileNames.size()
                  << ", "
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

      if(pIn[2] > pOut[2])
        {
        itkGenericExceptionMacro("Required condition pIn[2] > pOut[2] is not met, check coordinate system.");
        }

      // Apply transform
      VectorHType pInH(1.), pOutH(1.), dInH(1.), dOutH(1.);
      for(unsigned int i=0; i<3; i++)
        {
        pInH[i] = pIn[i];
        pOutH[i] = pOut[i];
        dInH[i] = dIn[i];
        dOutH[i] = dOut[i];
        }
      dInH[2]  = -1.*dInH[2];
      dOutH[2] = -1.*dOutH[2];

      pInH.SetVnlVector(m_GeometryIn->GetProjectionCoordinatesToFixedSystemMatrix(iProj).GetVnlMatrix() * pInH.GetVnlVector());
      pOutH.SetVnlVector(m_GeometryOut->GetProjectionCoordinatesToFixedSystemMatrix(iProj).GetVnlMatrix() * pOutH.GetVnlVector());
      dInH.SetVnlVector(m_GeometryIn->GetRotationMatrices()[iProj].GetInverse() * dInH.GetVnlVector());
      dOutH.SetVnlVector(m_GeometryOut->GetRotationMatrices()[iProj].GetInverse() * dOutH.GetVnlVector());

      for(unsigned int i=0; i<3; i++)
        {
        pIn[i] = pInH[i];
        pOut[i] = pOutH[i];
        dIn[i] = dInH[i];
        dOut[i] = dOutH[i];
        }

      // Select main direction
      unsigned int mainDir = 0;
      VectorType dInAbs;
      for(unsigned int i=0; i<3; i++)
        {
        dInAbs[i] = vnl_math_abs( dIn[i] );
        if(dInAbs[i]>dInAbs[mainDir])
          mainDir = i;
        }
      unsigned int dir1 = (mainDir+1)%3;
      unsigned int dir2 = (mainDir+2)%3;

      const double eIn = it.Get()[0];
      const double eOut = it.Get()[1];
      double value = 0.;
      value = m_ConvFunc->GetValue(eOut, eIn); //relative electron density
      ++it;

      // Move straight to entrance and exit shapes
      VectorType pSIn = pIn;
      VectorType pSOut = pOut;
      if(quadricIn.GetPointer()!=NULL)
        {
        quadricIn->SetRayOrigin(pSIn);
        quadricOut->SetRayOrigin(pSOut);
        if(quadricIn->Evaluate(dIn) && quadricOut->Evaluate(dOut))
          {
          pSIn  = pIn  + dIn  * quadricIn->GetNearestDistance();
          pSOut = pOut + dOut * quadricOut->GetNearestDistance();
          }
        }

      // Normalize direction with respect to z
      dIn[dir1] /= dIn[mainDir];
      dIn[dir2] /= dIn[mainDir];
      dIn[mainDir] = 1.;
      dOut[dir1] /= dOut[mainDir];
      dOut[dir2] /= dOut[mainDir];
      dOut[mainDir] = 1.;

      // Init MLP before mm to voxel conversion
      VectorType pSInO, pSOutO, dInO, dOutO;
      for(unsigned int i=0; i<3; i++)
        {
        pSInO[i] = pSIn[(mainDir+i+1)%3];
        pSOutO[i] = pSOut[(mainDir+i+1)%3];
        dInO[i] = dIn[(mainDir+i+1)%3];
        dOutO[i] = dOut[(mainDir+i+1)%3];
        }
      mlp->Init(pSInO, pSOutO, dInO, dOutO);
      VectorType dCurr = dIn;
      double xx = 0., yy = 0.;
      for(unsigned int k=0; k<imgSize[mainDir]; k++)
        {
        const double dk = zmm[mainDir][k];
        if((dInH[mainDir]>0 && dk<=pSIn[mainDir]) ||
           (dInH[mainDir]<0 && dk>=pSIn[mainDir])) //before entrance
          {
          const double z = (dk-pIn[mainDir]);
          xx = pIn[dir1]+z*dIn[dir1];
          yy = pIn[dir2]+z*dIn[dir2];
          dCurr = dIn;
          xx = (xx - imgOrigin[dir1]) * imgSpacingInv[dir1];
          yy = (yy - imgOrigin[dir2]) * imgSpacingInv[dir2];
          }
        else if((dInH[mainDir]>0 && dk>=pSOut[mainDir]) ||
                (dInH[mainDir]<0 && dk<=pSOut[mainDir])) //after exit
          {
          const double z = (dk-pSOut[mainDir]);
          xx = pSOut[dir1]+z*dOut[dir1];
          yy = pSOut[dir2]+z*dOut[dir2];
          dCurr = dOut;
          xx = (xx - imgOrigin[dir1]) * imgSpacingInv[dir1];
          yy = (yy - imgOrigin[dir2]) * imgSpacingInv[dir2];
          }
        else //MLP
          {
          dCurr[dir1] = xx;
          dCurr[dir2] = yy;
          dCurr[mainDir] = 1;
          mlp->Evaluate(zmm[mainDir][k], xx, yy);
          xx = (xx - imgOrigin[dir1]) * imgSpacingInv[dir1];
          yy = (yy - imgOrigin[dir2]) * imgSpacingInv[dir2];
          dCurr[dir1] = xx - dCurr[dir1];
          dCurr[dir2] = yy - dCurr[dir2];
          }

        // Lattice conversion
        typename OutputImageType::IndexType idx;
        idx[dir1] = itk::Math::Round<int,double>(xx);
        idx[dir2] = itk::Math::Round<int,double>(yy);
        if(idx[dir1]>=0 && idx[dir1]<(int)imgSize[dir1] &&
           idx[dir2]>=0 && idx[dir2]<(int)imgSize[dir2])
          {
          idx[mainDir] = k;
          if(dCurr[2]<0.)
            dCurr[0] *= -1.;
          double theta = acos(dCurr[0] / sqrt(dCurr[0]*dCurr[0]+dCurr[2]*dCurr[2]));
          theta *= imgSize[3] / itk::Math::pi;
          while(theta<0.)
            theta += imgSize[3];
          while(theta>=imgSize[3])
            theta -= imgSize[3];
          idx[3] = itk::Math::Floor<int, double>(theta);
          typename OutputImageType::OffsetValueType offset = m_Outputs[threadId]->ComputeOffset(idx);
          imgData[ offset ] += value;
          imgCountData[ offset ]++;
          }
        }
      }
      if(threadId==0)
        {
        std::cout << '\r'
                  << "Pair file #"
                  << iProj+1
                  << " out of "
                  << m_ProtonPairsFileNames.size()
                  << ", "
                  << region.GetSize(1) << " pairs of protons processed (100%) in thread 1"
                  << std::endl;
#ifdef MLP_TIMING
        mlp->PrintTiming(std::cout);
#endif
        }
    }
}

template <class TInputImage, class TOutputImage>
void
ProtonPairsToBackProjection<TInputImage, TOutputImage>
::AfterThreadedGenerateData()
{
  typedef typename itk::ImageRegionIterator<TOutputImage> ImageIteratorType;
  ImageIteratorType itOut(m_Outputs[0], m_Outputs[0]->GetLargestPossibleRegion());

  typedef itk::ImageRegionIterator<CountImageType> ImageCountIteratorType;
  ImageCountIteratorType itCOut(m_Counts[0], m_Outputs[0]->GetLargestPossibleRegion());

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

  // Set count image information
  m_Count->SetSpacing( this->GetOutput()->GetSpacing() );
  m_Count->SetOrigin( this->GetOutput()->GetOrigin() );

  // Normalize with proton count (average)
//  while(!itCOut.IsAtEnd())
//    {
//    if(itCOut.Get())
//      itOut.Set(itOut.Get()/itCOut.Get());
//    ++itOut;
//    ++itCOut;
//    }

  // Free images created in threads
  m_Outputs.resize( 0 );
  m_Counts.resize( 0 );
}

}
