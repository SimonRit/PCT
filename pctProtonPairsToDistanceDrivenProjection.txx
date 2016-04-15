#include <itkImageFileReader.h>
#include <itkImageRegionIterator.h>

#include "pctThirdOrderPolynomialMLPFunction.h"
#include "pctSchulteMLPFunction.h"
#include "pctEnergyStragglingFunctor.h"
#include "pctScatteringWEPLFunctor.h"


#include <iostream>
#include <TMath.h>
//#include <TROOT.h>
//#include "Math/Polynomial.h"
//#include "Math/Interpolator.h"

namespace pct
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
  itk::Vector<float, 3> imgSpacingInv;
  for(unsigned int i=0; i<3; i++)
    imgSpacingInv[i] = 1./imgSpacing[i];


  const unsigned long npixels = imgSize[0] * imgSize[1] * imgSize[2];
  std::vector< std::vector<double> > angles(npixels);       //container for exit angles
  std::vector< std::vector<double> > energies(npixels);     //container for exit energies


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

/*
  //----------------------ctq:
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

    const double anglex = vcl_acos( std::min(1.,dInX*dOutX / ( dInX.GetNorm() * dOutX.GetNorm() ) ) );
    const double angley = vcl_acos( std::min(1.,dInY*dOutY / ( dInY.GetNorm() * dOutY.GetNorm() ) ) );
    const double angle_out = sqrt(angley*angley+anglex*anglex);
  //----------------------ctq:
*/
    if(pIn[2] > pOut[2])
      {
      itkGenericExceptionMacro("Required condition pIn[2] > pOut[2] is not met, check coordinate system.");
      }

    eIn = it.Get()[0];
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
        angles[ idx ].push_back(angle_out);
        energies[ idx ].push_back(eOut);
        imgCountData[ idx ]++;
        }
      }
  }

  //----------------------ctq:  
    // Create scattering WEPL LUT
    unsigned int max_length = 500*CLHEP::mm;
    const double E0 = 13.6*CLHEP::MeV;
    const double X0 = 38.08*CLHEP::cm;
    const double proton_mass_c2 = CLHEP::proton_mass_c2;
    double integratedinvbeta2p2;
    double invbeta2p2;
    double energycurrent;
    int range;
    std::vector<double> sigma2_lut;
    std::vector<int> thickness_lut;

    for(int length = 1; length<max_length; length++)
      {    
      energycurrent = m_ConvFunc->GetEnergy(length*CLHEP::mm, eIn);
      if(energycurrent<0.001) break;        
      invbeta2p2 = (energycurrent+proton_mass_c2)*(energycurrent+proton_mass_c2) / 
                  ((energycurrent+2*proton_mass_c2)*(energycurrent+2*proton_mass_c2) * energycurrent*energycurrent); 
      integratedinvbeta2p2 += invbeta2p2;
      sigma2_lut.push_back( E0*E0/X0 * (1+0.038*vcl_log(length/X0))* (1+0.038*vcl_log(length/X0)) * integratedinvbeta2p2);

      thickness_lut.push_back(length);
      range = length;     
      }

    // Calculate angular variance per pixel
     for(unsigned int k=0; k<npixels; k++)
      {
      double meanAngle = 0.;
      double meanEnergy = 0.;
      double sigma2 = 0.;
  
      for(unsigned int i=0; i<imgCountData[k]; i++)
          {
          meanAngle += angles[k][i];
          meanEnergy += energies[k][i];
          }
       meanAngle=meanAngle/imgCountData[k];
       meanEnergy=meanEnergy/imgCountData[k];

      for(unsigned int i=0; i<imgCountData[k]; i++)
        {
        sigma2+=(angles[k][i]-meanAngle)*(angles[k][i]-meanAngle);
        }
        sigma2 = sigma2/imgCountData[k];

      double x0,y0,x1,y1;

      if(imgCountData[k] == 0.)
      imgData[ k ] = 0.;

      else
        {     
        for(int i=0; i<range; i++)
          {
          if(sigma2 < sigma2_lut[i])
          continue;
            {
            x0 = sigma2_lut[i];
            x1 = sigma2_lut[i+1];
            y0 = thickness_lut[i];
            y1 = thickness_lut[i+1];  
            }
          }
        
        // Calculate the scattering WEPL from LUT
        imgData[ k ] = y0 + (y1-y0) * (sigma2-x0) / (x1 - x0); 

/*        std::cout<<"\n"<<sigma2<<"\t"
                <<x0<<"\t"
                 <<x1<<"\t"
                <<y0<<"\t"
                <<y1<<"\t"
        <<std::endl;

        std::cout <<"\n\n"
                << imgData[ k ] << "\t"
                <<"\n"<< std::endl;
*/
        }
      }
  //----------------------ctq:


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
  while(!itCOut.IsAtEnd())
    {
    if(itCOut.Get())
//      itOut.Set(itOut.Get()/itCOut.Get());
      itOut.Set(itOut.Get());

    ++itOut;
    ++itCOut;
    }

  // Free images created in threads
  m_Outputs.resize( 0 );
  m_Counts.resize( 0 );
}

}
