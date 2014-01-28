#include "pctpaircuts_ggo.h"

#include <rtkMacro.h>
#include <rtkGgoFunctions.h>
#include <rtkConstantImageSource.h>

#include "pctProtonPairsToDistanceDrivenProjection.h"

#include <itkImageFileWriter.h>
#include <itkRegularExpressionSeriesFileNames.h>
#include <itkTimeProbe.h>

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooPlot.h>
#include <TCanvas.h>

#define PAIRS_IN_RAM 1000000

int main(int argc, char * argv[])
{
  GGO(pctpaircuts, args_info); //RTK macro parsing options from .ggo file (rtkMacro.h)

  typedef float ProjectionPixelType;
  typedef itk::Image< ProjectionPixelType, 2 > ProjectionImageType;
  typedef itk::Vector<double, 2> VectorTwoDType;

  // Create a stack of empty projection images and associated proton count
  typedef rtk::ConstantImageSource< ProjectionImageType > ConstantImageSourceType;
  ConstantImageSourceType::Pointer sumEnergy   = ConstantImageSourceType::New();
  ConstantImageSourceType::Pointer sumEnergySq = ConstantImageSourceType::New();
  ConstantImageSourceType::Pointer sumAngleSq  = ConstantImageSourceType::New();
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctpaircuts>(sumEnergy, args_info);
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctpaircuts>(sumEnergySq, args_info);
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctpaircuts>(sumAngleSq, args_info);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( sumEnergy->Update() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( sumEnergySq->Update() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( sumAngleSq->Update() );

  typedef itk::Image< unsigned int, 2 > CountImageType;
  typedef rtk::ConstantImageSource< CountImageType > CountImageSourceType;
  CountImageSourceType::Pointer counts = CountImageSourceType::New();
  rtk::SetConstantImageSourceFromGgo<CountImageSourceType, args_info_pctpaircuts>(counts, args_info);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( counts->Update() );

  // Robust case containers
  const size_t npixels = sumEnergy->GetOutput()->GetBufferedRegion().GetNumberOfPixels();
  std::vector< std::vector<double> > energies(npixels);
  std::vector< std::vector<double> > angles(npixels);

  // Read pairs
  typedef itk::Vector<float, 3> VectorType;
  typedef itk::Image<VectorType,2> PairsImageType;
  typedef itk::ImageFileReader< PairsImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( args_info.input_arg );
  reader->UpdateOutputInformation();
  size_t nprotons = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1]; //total image proton pairs number
  PairsImageType::RegionType region = reader->GetOutput()->GetLargestPossibleRegion();
  unsigned int nregions = nprotons/PAIRS_IN_RAM+1; //limit 1M proton pairs (memory)

  // Image information constants
  const ProjectionImageType::SizeType imgSize  = sumEnergy->GetOutput()->GetBufferedRegion().GetSize(); //pixels number (vector)
  const ProjectionImageType::PointType imgOrigin = sumEnergy->GetOutput()->GetOrigin(); //center position pixel (0,0) 
  const ProjectionImageType::SpacingType imgSpacing = sumEnergy->GetOutput()->GetSpacing(); //space between pixels 
  itk::Vector<float, 2> imgSpacingInv;
  for(unsigned int i=0; i<2; i++)
    imgSpacingInv[i] = 1./imgSpacing[i];

  // Data pointers
  float *pSumEnergy     = sumEnergy->GetOutput()->GetBufferPointer();
  float *pSumEnergySq   = sumEnergySq->GetOutput()->GetBufferPointer();
  float *pSumAngleSq    = sumAngleSq->GetOutput()->GetBufferPointer();
  unsigned int *pCounts = counts->GetOutput()->GetBufferPointer();

  std::cout << "Compute cuts..." << std::endl;
  for(unsigned int r=0; r<nregions; r++)
  {
    // Read r-th set of pairs
    region.SetIndex(1, r*PAIRS_IN_RAM);
    region.SetSize(1, vnl_math_min(PAIRS_IN_RAM, int(nprotons-region.GetIndex(1))));
    reader->GetOutput()->SetRequestedRegion(region); //we work on one region "r"
    TRY_AND_EXIT_ON_ITK_EXCEPTION( reader->Update() );

    // Process pairs
    itk::ImageRegionIterator<PairsImageType> it(reader->GetOutput(), region);
    for(size_t p=region.GetIndex(1); p<region.GetIndex(1)+region.GetSize(1); p++)
      {
      const VectorType pIn = it.Get();
      ++it;
      const VectorType pOut = it.Get();
      ++it;
      const VectorType dIn = it.Get();
      ++it;
      const VectorType dOut = it.Get();
      ++it;
      const VectorType data = it.Get();
      ++it;

      static double mag = (args_info.source_arg - pOut[2]) / (args_info.source_arg - pIn[2]); //Magnification (Thales theorem)

      const double xx = (pIn[0]*mag-imgOrigin[0]) * imgSpacingInv[0]; //x corrected (mag), converted in pixel units
      const int i = itk::Math::Round<int,double>(xx);
      if(i<0 || i>=(int)imgSize[0])
        continue;

      const double yy = (pIn[1]*mag-imgOrigin[1]) * imgSpacingInv[1];
      const int j = itk::Math::Round<int,double>(yy);
      if(j<0 || j>=(int)imgSize[1])
        continue;

      const unsigned long idx = i+j*imgSize[0];

      VectorTwoDType dInX, dInY, dOutX, dOutY;
      dInX[0] = dIn[0];
      dInX[1] = dIn[2];
      dInY[0] = dIn[1];
      dInY[1] = dIn[2];
      dOutX[0] = dOut[0];
      dOutX[1] = dOut[2];
      dOutY[0] = dOut[1];
      dOutY[1] = dOut[2];
      const double anglex = vcl_acos( std::min(1.,dInX*dOutX) / ( dInX.GetNorm() * dOutX.GetNorm() ) );
      const double angley = vcl_acos( std::min(1.,dInY*dOutY) / ( dInY.GetNorm() * dOutY.GetNorm() ) );
      const double energy = data[0]-data[1];

      if(args_info.robust_flag || (args_info.plotpix_given && idx==(unsigned long)args_info.plotpix_arg ) )
        {
        energies[idx].push_back(energy);
        angles  [idx].push_back(anglex);
        angles  [idx].push_back(angley);
        }
      if(!args_info.robust_flag)
        {
        pSumEnergy  [idx] += energy;
        pSumEnergySq[idx] += energy*energy;
        pSumAngleSq [idx] += anglex*anglex;
        pSumAngleSq [idx] += angley*angley;
        }
      pCounts     [idx]++;
      }
  }

  std::cout << "Finalize cuts..." << std::endl;

  // Now compute the cuts and the average
  if(args_info.robust_flag)
    {
    for(unsigned int idx=0; idx<npixels; idx++)
      {
      if(pCounts[idx]==1)
        {
        // Just one event in this pixel, keep it!
        pSumEnergy  [idx] = energies[idx][0];
        pSumEnergySq[idx] = 0.1;
        pSumAngleSq [idx] = angles[idx][0];;
        }
      else if(pCounts[idx])
        {
        // Energy: median and 30.85% (0.5 sigma) with interpolation
        double medianPos = pCounts[idx]*0.5;
        unsigned int medianSupPos = itk::Math::Ceil<unsigned int, double>(medianPos);
        std::partial_sort(energies[idx].begin(), energies[idx].begin()+medianSupPos+1, energies[idx].end());
        double medianDiff = medianSupPos-medianPos;
	//median linear interpolation
        pSumEnergy  [idx] = *(energies[idx].begin()+medianSupPos)*(1.-medianDiff)+ //tab[x][y] <=> *(tab[x]+y)
	  *(energies[idx].begin()+medianSupPos-1)*medianDiff; 

        double sigmaEPos = pCounts[idx]*0.3085; //0.5 sigma
        unsigned int sigmaESupPos = itk::Math::Ceil<unsigned int, double>(sigmaEPos);
        double sigmaEDiff = sigmaESupPos-sigmaEPos;
        pSumEnergySq[idx] = 2.*(pSumEnergy[idx]-( *(energies[idx].begin()+sigmaESupPos)*(1.-sigmaEDiff)+
                                                  *(energies[idx].begin()+sigmaESupPos-1)*sigmaEDiff) ); //x2 to get 1sigma

        // Angle: 38.30% (0.5 sigma) with interpolation (median is 0. and we only have positive values
        double sigmaAPos = angles[idx].size()*0.3830;
        unsigned int sigmaASupPos = itk::Math::Ceil<unsigned int, double>(sigmaAPos);
        std::partial_sort(angles[idx].begin(), angles[idx].begin()+sigmaASupPos+1, angles[idx].end());
        double sigmaADiff = sigmaASupPos-sigmaAPos;
        pSumAngleSq[idx] = 2.*(*(angles[idx].begin()+sigmaASupPos)*(1.-sigmaADiff)+
                               *(angles[idx].begin()+sigmaASupPos-1)*sigmaADiff); //x2 to get 1sigma
        }
      }
    }
  else
    {
    for(unsigned int idx=0; idx<npixels; idx++)
      {
      if(pCounts[idx])
        {
        pSumEnergy  [idx] /= pCounts[idx];
        pSumEnergySq[idx] /= pCounts[idx];
        pSumAngleSq [idx] /= 2*pCounts[idx];
        pSumEnergySq[idx] = sqrt( pSumEnergySq[idx]-pSumEnergy[idx]*pSumEnergy[idx] );
        pSumAngleSq [idx] = sqrt( pSumAngleSq [idx] );
        }
      }
    }

  // Weight standard deviations with parameters
  for(unsigned int idx=0; idx<npixels; idx++)
    {
    pSumEnergySq[idx] *= args_info.energycut_arg;
    pSumAngleSq[idx]  *= args_info.anglecut_arg;
    }

  std::cout << "Select pairs..." << std::endl;

  // And select the pairs
  std::vector<VectorType> pairs;
  for(unsigned int r=0; r<nregions; r++)
  {
    // Read r-th set of pairs
    region.SetIndex(1, r*PAIRS_IN_RAM);
    region.SetSize(1, vnl_math_min(PAIRS_IN_RAM, int(nprotons-region.GetIndex(1))));
    reader->GetOutput()->SetRequestedRegion(region);
    TRY_AND_EXIT_ON_ITK_EXCEPTION( reader->Update() );

    // Process pairs
    itk::ImageRegionIterator<PairsImageType> it(reader->GetOutput(), region);
    for(size_t p=region.GetIndex(1); p<region.GetIndex(1)+region.GetSize(1); p++)
      {
      const VectorType pIn = it.Get();
      ++it;
      const VectorType pOut = it.Get();
      ++it;
      const VectorType dIn = it.Get();
      ++it;
      const VectorType dOut = it.Get();
      ++it;
      const VectorType data = it.Get();
      ++it;

      static double mag = (args_info.source_arg - pOut[2]) / (args_info.source_arg - pIn[2]);

      const double xx = (pIn[0]*mag-imgOrigin[0]) * imgSpacingInv[0];
      const int i = itk::Math::Round<int,double>(xx);
      if(i<0 || i>=(int)imgSize[0])
        continue;

      const double yy = (pIn[1]*mag-imgOrigin[1]) * imgSpacingInv[1];
      const int j = itk::Math::Round<int,double>(yy);
      if(j<0 || j>=(int)imgSize[1])
        continue;

      const unsigned long idx = i+j*imgSize[0];

      VectorTwoDType dInX, dInY, dOutX, dOutY;
      dInX[0] = dIn[0];
      dInX[1] = dIn[2];
      dInY[0] = dIn[1];
      dInY[1] = dIn[2];
      dOutX[0] = dOut[0];
      dOutX[1] = dOut[2];
      dOutY[0] = dOut[1];
      dOutY[1] = dOut[2];
      const double anglex = vcl_acos( std::min(1.,dInX*dOutX) / ( dInX.GetNorm() * dOutX.GetNorm() ) );
      const double angley = vcl_acos( std::min(1.,dInY*dOutY) / ( dInY.GetNorm() * dOutY.GetNorm() ) );
      const double energy = data[0]-data[1];

      if( anglex <= pSumAngleSq [idx] &&
          angley <= pSumAngleSq [idx] &&
          vcl_abs(energy-pSumEnergy[idx]) <= pSumEnergySq[idx] )
        {
        pairs.push_back(pIn);
        pairs.push_back(pOut);
        pairs.push_back(dIn);
        pairs.push_back(dOut);
        pairs.push_back(data);
        }
    }
  }

  std::cout << "Write pairs..." << std::endl;

  itk::ImageRegion<2> pairsRegion;
  itk::ImageRegion<2>::SizeType size;
  size[0] = 5;
  size[1] = pairs.size()/5;
  pairsRegion.SetSize(size);
  PairsImageType::Pointer img = PairsImageType::New();
  img->SetRegions(pairsRegion);
  img->Allocate();

  itk::ImageRegionIterator<PairsImageType> it(img, pairsRegion);
  for(size_t i=0; i<pairs.size(); ++i, ++it)
    it.Set( pairs[i] );

  // Write
  typedef itk::ImageFileWriter<PairsImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( args_info.output_arg );
  writer->SetInput( img );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() );

  // Optional outputs
  typedef itk::ImageFileWriter<ProjectionImageType> PWriterType;
  if(args_info.sangle_given)
    {
    PWriterType::Pointer w = PWriterType::New();
    w->SetInput(sumAngleSq->GetOutput());
    w->SetFileName(args_info.sangle_arg);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(w->Update());
    }
  if(args_info.menergy_given)
    {
    PWriterType::Pointer w = PWriterType::New();
    w->SetInput(sumEnergy->GetOutput());
    w->SetFileName(args_info.menergy_arg);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(w->Update());
    }  typedef itk::ImageFileWriter<ProjectionImageType> PWriterType;
  if(args_info.senergy_given)
    {
    PWriterType::Pointer w = PWriterType::New();
    w->SetInput(sumEnergySq->GetOutput());
    w->SetFileName(args_info.senergy_arg);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(w->Update());
    }
  if(args_info.count_given)
    {
    itk::ImageFileWriter<CountImageType>::Pointer w;
    w = itk::ImageFileWriter<CountImageType>::New();
    w->SetInput(counts->GetOutput());
    w->SetFileName(args_info.count_arg);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(w->Update());
    }
  if(args_info.plotpix_given)
    {
    // Energy plot
    const unsigned int p = args_info.plotpix_arg;
    if(pSumEnergySq[p]==0.)
      std::cout << "Can not create energy.pdf, sigma is 0." << std::endl;
    else
      {
      RooRealVar rooEnergy("Energy","Energy loss",0.,500.*CLHEP::MeV, "MeV");
      RooDataSet rooEnergyData("data", "data", RooArgSet(rooEnergy));
      for(unsigned int i=0; i<pCounts[p]; i++)
        {
        rooEnergy = energies[p][i]/CLHEP::MeV;
        rooEnergyData.add(RooArgSet(rooEnergy));
        }

      TCanvas tEnergyCanvas("energyCanvas","Energy canvas",0,0,1000,500);
      RooPlot* rooEnergyPlot = rooEnergy.frame(pSumEnergy[p]-2*pSumEnergySq[p],
                                               pSumEnergy[p]+2*pSumEnergySq[p]);
      rooEnergyData.plotOn(rooEnergyPlot);

      RooRealVar rooMeanEnergy("rooMeanEnergy","mean of energy gaussian",pSumEnergy[p],0.,500.*CLHEP::MeV) ;
      RooRealVar rooSigmaEnergy("rooSigmaEnergy","width of energy gaussian",pSumEnergySq[p]/args_info.energycut_arg,0.,500.*CLHEP::MeV);
      RooGaussian rooGaussEnergy("rooGaussEnergy","energy gaussian PDF",rooEnergy,rooMeanEnergy,rooSigmaEnergy) ;

      rooGaussEnergy.plotOn(rooEnergyPlot, RooFit::LineColor(kBlue));
      //rooGaussEnergy.fitTo(rooEnergyData);
      //rooGaussEnergy.plotOn(rooEnergyPlot, LineColor(kRed));

      rooEnergyPlot->Draw();
      tEnergyCanvas.SaveAs("energy.pdf");
      }

    // Angle plot
    if(pSumAngleSq[p]==0.)
      std::cout << "Can not create angle.pdf, sigma is 0." << std::endl;
    else
      {
      RooRealVar rooAngle("Angle","Angle deviation",0.,2.*itk::Math::pi, "rad");
      RooDataSet rooAngleData("data", "data", RooArgSet(rooAngle));
      for(unsigned int i=0; i<angles[p].size(); i++)
        {
        rooAngle = angles[p][i];
        rooAngleData.add(RooArgSet(rooAngle));
        }

      TCanvas tAngleCanvas("angleCanvas","Angle canvas",0,0,1000,500);
      RooPlot* rooAnglePlot = rooAngle.frame(0., 2*pSumAngleSq[p]);
      rooAngleData.plotOn(rooAnglePlot);

      RooRealVar rooMeanAngle("rooMeanAngle","mean of angle gaussian",0.,0.,500.*CLHEP::MeV);
      RooRealVar rooSigmaAngle("rooSigmaAngle","width of angle gaussian",pSumAngleSq[p]/args_info.anglecut_arg,0.,500.*CLHEP::MeV);
      RooGaussian rooGaussAngle("rooGaussAngle","angle gaussian PDF",rooAngle,rooMeanAngle,rooSigmaAngle);

      rooGaussAngle.plotOn(rooAnglePlot, RooFit::LineColor(kBlue));
      //rooGaussAngle.fitTo(rooAngleData);
      //rooGaussAngle.plotOn(rooAnglePlot, LineColor(kRed));

      rooAnglePlot->Draw();
      tAngleCanvas.SaveAs("angle.pdf");
      }
    }

  return EXIT_SUCCESS;
}
