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
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>
#include <TLegend.h>

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
  ConstantImageSourceType::Pointer sigmaAngle  = ConstantImageSourceType::New();

  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctpaircuts>(sumEnergy, args_info);
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctpaircuts>(sumEnergySq, args_info);
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctpaircuts>(sumAngleSq, args_info);
  rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_pctpaircuts>(sigmaAngle, args_info);

  TRY_AND_EXIT_ON_ITK_EXCEPTION( sumEnergy->Update() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( sumEnergySq->Update() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( sumAngleSq->Update() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( sigmaAngle->Update() );

  typedef itk::Image< unsigned int, 2 > CountImageType;
  typedef rtk::ConstantImageSource< CountImageType > CountImageSourceType;
  CountImageSourceType::Pointer counts = CountImageSourceType::New();
  CountImageSourceType::Pointer countsCut = CountImageSourceType::New();

  rtk::SetConstantImageSourceFromGgo<CountImageSourceType, args_info_pctpaircuts>(counts, args_info);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( counts->Update() );
  rtk::SetConstantImageSourceFromGgo<CountImageSourceType, args_info_pctpaircuts>(countsCut, args_info);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( countsCut->Update() );

  // Robust case containers
  const size_t npixels = sumEnergy->GetOutput()->GetBufferedRegion().GetNumberOfPixels();
  std::vector< std::vector<double> > energies(npixels);
  std::vector< std::vector<double> > angles(npixels);
  std::vector< std::vector<double> > exitEnergies(npixels);
  std::vector< std::vector<unsigned int> > trackID(npixels);
  std::vector< std::vector<double> > anglesX(npixels);
  std::vector< std::vector<double> > anglesY(npixels);
  std::vector< std::vector<unsigned int> > creatorProcess(npixels);
  std::vector< std::vector<unsigned int> > nuclearProcess(npixels);

  std::vector< std::vector<double> > exitEnergiesCut(npixels);
  std::vector< std::vector<double> > exitAnglesCut(npixels);

  double gminEnergy[npixels];

  double gsigmaAngle[npixels];
  double gmaxAngle[npixels];

  float Ein = 0.;

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
  float *pSumEnergy       = sumEnergy->GetOutput()->GetBufferPointer();
  float *pSumEnergySq     = sumEnergySq->GetOutput()->GetBufferPointer();
  float *pSumAngleSq      = sumAngleSq->GetOutput()->GetBufferPointer();
  float *pSigmaAngle      = sigmaAngle->GetOutput()->GetBufferPointer();

  unsigned int *pCounts   = counts->GetOutput()->GetBufferPointer();
  unsigned int *pCountsCut   = countsCut->GetOutput()->GetBufferPointer();

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
      VectorType nuclearinfo;
      if(region.GetSize(0)==6)
        {
        nuclearinfo = it.Get();
        ++it;
        }

      // Magnification (Thales theorem)
      static double mag = (args_info.source_arg==0.)?1.:(args_info.source_arg - pOut[2]) / (args_info.source_arg - pIn[2]);

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

      double anglex = vcl_acos( std::min(1.,dInX*dOutX / ( dInX.GetNorm() * dOutX.GetNorm() ) ) );
      double angley = vcl_acos( std::min(1.,dInY*dOutY / ( dInY.GetNorm() * dOutY.GetNorm() ) ) );
      const double energy = data[0]-data[1];

      Ein = data[0];

      if(args_info.robust_flag
        || (args_info.plotpix_given && idx==(unsigned long)args_info.plotpix_arg )
        || args_info.gauscut_flag
        || args_info.plotSuperDistribution_flag)
        {
        energies[idx].push_back(energy);
        angles  [idx].push_back(anglex);
        angles  [idx].push_back(angley);
        exitEnergies[idx].push_back(data[1]);
        trackID[idx].push_back(data[2]);
        anglesX  [idx].push_back(anglex);
        anglesY  [idx].push_back(angley);

        creatorProcess[idx].push_back(nuclearinfo[0]);
        nuclearProcess[idx].push_back(nuclearinfo[1]);
        }

      if(!args_info.robust_flag)
        {
        pSumEnergy  [idx] += energy;
        pSumEnergySq[idx] += energy*energy;
        pSumAngleSq [idx] += anglex*anglex;
        pSumAngleSq [idx] += angley*angley;
        }
      pCounts[ idx ]++;
      }
  }

//===================================================================================================================================
// Plot raw biplots
//===================================================================================================================================
  if(args_info.plotRawDistribution_given)
    {
    gROOT->SetBatch(kTRUE);
    for(unsigned int idx=0; idx<npixels; idx++)
      {
      char biPlotX_00[200], filename_00[100], foldername_00[100];
      char biPlotX_01[200], filename_01[100], foldername_01[100];
      char biPlotX_02[200], filename_02[100], foldername_02[100];
      char biPlotX_03[200], filename_03[100], foldername_03[100];
      char biPlotX_04[200], filename_04[100], foldername_04[100];
      char biPlotX_05[200], filename_05[100], foldername_05[100];
      char biPlotX_06[200], filename_06[100], foldername_06[100];
      char biPlotX_07[200], filename_07[100], foldername_07[100];
      char biPlotX_08[200], filename_08[100], foldername_08[100];
      char biPlotX_09[200], filename_09[100], foldername_09[100];
      char biPlotX_10[200], filename_10[100], foldername_10[100];

      sprintf( biPlotX_00, "Biplot of all protons in pixel index %d", idx);
      sprintf( biPlotX_01, "Biplot of all primary protons in pixel index %d", idx);
      sprintf( biPlotX_02, "Biplot of primary protons without elastic scattering in pixel index %d", idx);
      sprintf( biPlotX_03, "Biplot of primary protons with elastic scattering in pixel index %d", idx);
      sprintf( biPlotX_04, "Biplot of all secondary protons in pixel index %d", idx);
      sprintf( biPlotX_05, "Biplot of secondary protons created from inelastic scattering in pixel index %d", idx);
      sprintf( biPlotX_06, "Biplot of secondary protons created from inelastic scattering and without elastic interaction in pixel index %d", idx);
      sprintf( biPlotX_07, "Biplot of secondary protons created from inelastic scattering and with elastic interaction in pixel index %d", idx);
      sprintf( biPlotX_08, "Biplot of secondary protons created from elastic scattering in pixel index %d", idx);
      sprintf( biPlotX_09, "Biplot of secondary protons created from elastic scattering and without elastic interaction in pixel index %d", idx);
      sprintf( biPlotX_10, "Biplot of secondary protons created from elastic scattering and with elastic interaction in pixel index %d", idx);

      sprintf( foldername_00,  "mkdir -p biplot_00");  system(foldername_00);
      sprintf( foldername_01,  "mkdir -p biplot_01");  system(foldername_01);
      sprintf( foldername_02,  "mkdir -p biplot_02");  system(foldername_02);
      sprintf( foldername_03,  "mkdir -p biplot_03");  system(foldername_03);
      sprintf( foldername_04,  "mkdir -p biplot_04");  system(foldername_04);
      sprintf( foldername_05,  "mkdir -p biplot_05");  system(foldername_05);
      sprintf( foldername_06,  "mkdir -p biplot_06");  system(foldername_06);
      sprintf( foldername_07,  "mkdir -p biplot_07");  system(foldername_07);
      sprintf( foldername_08,  "mkdir -p biplot_08");  system(foldername_08);
      sprintf( foldername_09,  "mkdir -p biplot_09");  system(foldername_09);
      sprintf( foldername_10,  "mkdir -p biplot_10");  system(foldername_10);

      sprintf( filename_00, "biplot_00/biplotX_%d.png", idx);
      sprintf( filename_01, "biplot_01/biplotX_%d.png", idx);
      sprintf( filename_02, "biplot_02/biplotX_%d.png", idx);
      sprintf( filename_03, "biplot_03/biplotX_%d.png", idx);
      sprintf( filename_04, "biplot_04/biplotX_%d.png", idx);
      sprintf( filename_05, "biplot_05/biplotX_%d.png", idx);
      sprintf( filename_06, "biplot_06/biplotX_%d.png", idx);
      sprintf( filename_07, "biplot_07/biplotX_%d.png", idx);
      sprintf( filename_08, "biplot_08/biplotX_%d.png", idx);
      sprintf( filename_09, "biplot_09/biplotX_%d.png", idx);
      sprintf( filename_10, "biplot_10/biplotX_%d.png", idx);

      TH2F *hx_00 = new TH2F("hx_00", biPlotX_00, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_01 = new TH2F("hx_01", biPlotX_01, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_02 = new TH2F("hx_02", biPlotX_02, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_03 = new TH2F("hx_03", biPlotX_03, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_04 = new TH2F("hx_04", biPlotX_04, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_05 = new TH2F("hx_05", biPlotX_05, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_06 = new TH2F("hx_06", biPlotX_06, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_07 = new TH2F("hx_07", biPlotX_07, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_08 = new TH2F("hx_08", biPlotX_08, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_09 = new TH2F("hx_09", biPlotX_09, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_10 = new TH2F("hx_10", biPlotX_10, 210, 0, 210., 90, 0, 90.);

      for(unsigned int i=0; i<pCounts[idx]; i++)
        {
        // all protons (primaries and secondaries)
        hx_00->Fill(exitEnergies[idx][i],anglesX[idx][i] * 180./TMath::Pi());

        // all primaries
        if(creatorProcess[idx][i]==0)
          {
          hx_01->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
          }

        // primaries without elastic scattering
        if(creatorProcess[idx][i]==0 && nuclearProcess[idx][i]==0)
          {
          hx_02->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
          }

        // primaries with elastic scattering
        if(creatorProcess[idx][i]==0 && nuclearProcess[idx][i]==1)
          {
          hx_03->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
          }

        // all secondaries
        if(creatorProcess[idx][i]>0)
          {
          hx_04->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
          }

        // all secondaries created from inelastic scattering
        if(creatorProcess[idx][i]==2)
          {
          hx_05->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
          }

        // secondaries created from inelastic scattering and no elastic scattering
        if(creatorProcess[idx][i]==2 && nuclearProcess[idx][i]==0)
          {
          hx_06->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
          }

        // secondaries created from inelastic scattering and with elastic scattering
        if(creatorProcess[idx][i]==2 && nuclearProcess[idx][i]==1)
          {
          hx_07->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
          }

        // all secondaries created from elastic scattering
        if(creatorProcess[idx][i]==1)
          {
          hx_08->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
          }

        // all secondaries created from elastic scattering and without elastic scattering
        if(creatorProcess[idx][i]==1 && nuclearProcess[idx][i]==0)
          {
          hx_09->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
          }

        // all secondaries created from elastic scattering and with elastic scattering
        if(creatorProcess[idx][i]==1 && nuclearProcess[idx][i]==1)
          {
          hx_10->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
          }
        }

        hx_00->Draw("colz");
        hx_00->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_00->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_00);
        delete hx_00;

        hx_01->Draw("colz");
        hx_01->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_01->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_01);
        delete hx_01;

        hx_02->Draw("colz");
        hx_02->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_02->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_02);
        delete hx_02;

        hx_03->Draw("colz");
        hx_03->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_03->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_03);
        delete hx_03;

        hx_04->Draw("colz");
        hx_04->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_04->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_04);
        delete hx_04;

        hx_05->Draw("colz");
        hx_05->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_05->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_05);
        delete hx_05;

        hx_06->Draw("colz");
        hx_06->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_06->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_06);
        delete hx_06;

        hx_07->Draw("colz");
        hx_07->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_07->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_07);
        delete hx_07;

        hx_08->Draw("colz");
        hx_08->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_08->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_08);
        delete hx_08;

        hx_09->Draw("colz");
        hx_09->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_09->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_09);
        delete hx_09;

        hx_10->Draw("colz");
        hx_10->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_10->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_10);
        delete hx_10;
        }
      }


//===================================================================================================================================
// Finalize cuts
//===================================================================================================================================

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
//===================================================================================================================================
  else if(args_info.gauscut_flag)
    {
    gROOT->SetBatch(kTRUE);
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
        TH1F *h1 = new TH1F("h1", "", 210, 0., 210.);
        TH1F *h2 = new TH1F("h2", "", 200, -1, 1);
        TF1 *gausFit = new TF1 ("gausFit", "gaus");

        for(unsigned int i=0; i<pCounts[idx]; i++)
          {
          h1->Fill(exitEnergies[idx][i]);
          h2->Fill(anglesX[idx][i]);
          }

        if (h1->GetRMS()>=1.)
          {
          h1->Fit("gausFit");
          float mean_energy = gausFit->GetParameter(1);
          float sigma_energy = gausFit->GetParameter(2);
          float fitMin_energy = mean_energy - args_info.energycut_arg*sigma_energy;
          float fitMax_energy = mean_energy + args_info.energycut_arg*sigma_energy;
          h1->Fit(gausFit, "", "", fitMin_energy, fitMax_energy);

          h2->Fit("gausFit");
          //h2->Fit(gausFit, "", "", -1., 1.);
          float mean_angle = gausFit->GetParameter(1);
          float sigma_angle = gausFit->GetParameter(2);
          float fitMin_angle = mean_angle - args_info.anglecut_arg*sigma_angle;
          float fitMax_angle = mean_angle + args_info.anglecut_arg*sigma_angle;
          h2->Fit(gausFit, "", "", fitMin_angle, fitMax_angle);

          gminEnergy[idx] = mean_energy - args_info.energycut_arg*sigma_energy;

          gsigmaAngle[idx] = gausFit->GetParameter(2);
          gmaxAngle[idx] = mean_angle + args_info.anglecut_arg*sigma_angle;

          pSigmaAngle[idx] = gsigmaAngle[idx]; //sigma_angle; //gsigmaAngle[idx];
          }

        if (h1->GetRMS()<1.)
          {
          gminEnergy[idx] = 0.;

          gmaxAngle[idx] = 2.;
          }

         if(args_info.plotproj_given)
            {
            char exitenergyfilename[100], energymakefolder[100];
            char exitanglefilename[100], anglemakefolder[100];

            sprintf( energymakefolder, "mkdir -p figs_energy_proj");
            sprintf( anglemakefolder,  "mkdir -p figs_angle_proj");

            system(energymakefolder);
            system(anglemakefolder);

            sprintf( exitenergyfilename, "figs_energy_proj/exitenergy_%d.png", idx);
            sprintf( exitanglefilename, "figs_angle_proj/exitangle_%d.png", idx);

            h1->Draw();
            h1->GetXaxis()->SetTitle("Exit energy (MeV)");
            h1->GetYaxis()->SetTitle("Number of protons");
            gPad->Print(exitenergyfilename);

            h2->Draw();
            h2->GetXaxis()->SetTitle("Exit angle (rad)");
            h2->GetYaxis()->SetTitle("Number of protons");
            gPad->Print(exitanglefilename);

            }
          delete h1;
          delete h2;
        }
      }
    }
//===================================================================================================================================
  else
    {
    for(unsigned int idx=0; idx<npixels; idx++)
      {
      //calculate mean
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

    if(pSumAngleSq[idx]==0.0)
      pSumAngleSq[idx] = 1.;        // If the sigma is below 0.01, don't perform a cut. The pixel lies on the boundaries.
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
      VectorType nuclearinfo;
      if(region.GetSize(0)==6)
        {
        nuclearinfo = it.Get();
        ++it;
        }

      static double mag = (args_info.source_arg==0.)?1.:(args_info.source_arg - pOut[2]) / (args_info.source_arg - pIn[2]);

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

      double anglex = vcl_acos( std::min(1.,dInX*dOutX / ( dInX.GetNorm() * dOutX.GetNorm() ) ) );
      double angley = vcl_acos( std::min(1.,dInY*dOutY / ( dInY.GetNorm() * dOutY.GetNorm() ) ) );
      const double energy = data[0]-data[1];

      if(args_info.gauscut_flag
         && data[1] >= gminEnergy[idx] // && data[1] <= gmaxEnergy[idx] &&
         && anglex  <= gmaxAngle[idx]  // && anglex  >= gminAngle[idx]
         && angley  <= gmaxAngle[idx]  // && angley  >= gminAngle[idx]
         )
        {
        pairs.push_back(pIn);
        pairs.push_back(pOut);
        pairs.push_back(dIn);
        pairs.push_back(dOut);
        pairs.push_back(data);
        if(region.GetSize(0)==6)
          pairs.push_back(nuclearinfo);

        exitEnergiesCut[idx].push_back(data[1]);
        exitAnglesCut[idx].push_back(anglex);
        exitAnglesCut[idx].push_back(angley);
        pCountsCut[idx]++;
        }

      else if( !args_info.gauscut_flag &&
          anglex <= pSumAngleSq [idx] && // Does not work with attenuation CT. Update: It works now, with some modifications above.
          angley <= pSumAngleSq [idx] &&
          vcl_abs(energy-pSumEnergy[idx]) <= pSumEnergySq[idx] )
        {
          if (args_info.primaries_flag && data[2]==1)
            {
              pairs.push_back(pIn);
              pairs.push_back(pOut);
              pairs.push_back(dIn);
              pairs.push_back(dOut);
              pairs.push_back(data);
              if(region.GetSize(0)==6)
                pairs.push_back(nuclearinfo);

              exitEnergiesCut[idx].push_back(data[1]);
              exitAnglesCut[idx].push_back(anglex);
              exitAnglesCut[idx].push_back(angley);
              pCountsCut[idx]++;
            }

          else if (args_info.nonuclear_flag && nuclearinfo[0]==0 && nuclearinfo[1]==0)
            {
              pairs.push_back(pIn);
              pairs.push_back(pOut);
              pairs.push_back(dIn);
              pairs.push_back(dOut);
              pairs.push_back(data);
              if(region.GetSize(0)==6)
                pairs.push_back(nuclearinfo);

              exitEnergiesCut[idx].push_back(data[1]);
              exitAnglesCut[idx].push_back(anglex);
              exitAnglesCut[idx].push_back(angley);
              pCountsCut[idx]++;
            }


          else if (!args_info.primaries_flag && !args_info.nonuclear_flag)
            {
              pairs.push_back(pIn);
              pairs.push_back(pOut);
              pairs.push_back(dIn);
              pairs.push_back(dOut);
              pairs.push_back(data);
              if(region.GetSize(0)==6)
                pairs.push_back(nuclearinfo);

              exitEnergiesCut[idx].push_back(data[1]);
              exitAnglesCut[idx].push_back(anglex);
              exitAnglesCut[idx].push_back(angley);
              pCountsCut[idx]++;
            }
        }
     }
    }


//===================================================================================================================================
// Plot superimposed raw distributions
//===================================================================================================================================
  unsigned int rawcounter00 = 0;
  if(args_info.plotSuperDistribution_given)
    {
    unsigned short int eMax[npixels];
    unsigned short int eMin[npixels];

    for(unsigned int idx=0; idx<npixels; idx++)
      {
      eMax[idx] = 0.;
      eMin[idx] = 1000.;

      for(unsigned int i=0; i<pCounts[idx]; i++)
          {
         if (eMax[idx]<exitEnergies[idx][i])
            eMax[idx] = exitEnergies[idx][i];
         if (eMin[idx]>exitEnergies[idx][i])
            eMin[idx] = exitEnergies[idx][i];
          }
       eMax[idx] +=1;
       eMin[idx] =50.;

        }

    for(unsigned int idx=0; idx<npixels; idx++)
      {
      char sPlotE[200], sfilenameE[100], sfoldernameE[100];
      sprintf( sPlotE, "Energy distributions in pixel index %d", idx);
      sprintf( sfoldernameE,  "mkdir -p splotE");  system(sfoldernameE);
      sprintf( sfilenameE, "splotE/splotE_%d.png", idx);

      char sPlotA[200], sfilenameA[100], sfoldernameA[100];
      sprintf( sPlotA, "Angular distributions in pixel index %d", idx);
      sprintf( sfoldernameA,  "mkdir -p splotA");  system(sfoldernameA);
      sprintf( sfilenameA, "splotA/splotA_%d.png", idx);

      TH1F *shistEnergy_00 = new TH1F("shistEnergy_00", sPlotE, eMax[idx]*5, eMin[idx], eMax[idx]);
      TH1F *shistEnergy_01 = new TH1F("shistEnergy_01", sPlotE, eMax[idx]*5, eMin[idx], eMax[idx]);
      TH1F *shistEnergy_02 = new TH1F("shistEnergy_02", sPlotE, eMax[idx]*5, eMin[idx], eMax[idx]);
      TH1F *shistEnergy_03 = new TH1F("shistEnergy_03", sPlotE, eMax[idx]*5, eMin[idx], eMax[idx]);
      TH1F *shistEnergy_04 = new TH1F("shistEnergy_04", sPlotE, eMax[idx]*5, eMin[idx], eMax[idx]);

      TH1F *shistAngle_00 = new TH1F("shistAngle_00", sPlotA, 300, 0., 60.);
      TH1F *shistAngle_01 = new TH1F("shistAngle_01", sPlotA, 300, 0., 60.);
      TH1F *shistAngle_02 = new TH1F("shistAngle_02", sPlotA, 300, 0., 60.);
      TH1F *shistAngle_03 = new TH1F("shistAngle_03", sPlotA, 300, 0., 60.);
      TH1F *shistAngle_04 = new TH1F("shistAngle_04", sPlotA, 300, 0., 60.);

      unsigned int rawcounter01 = 0;
      unsigned int rawcounter02 = 0;
      unsigned int rawcounter03 = 0;
      unsigned int cutcounter00 = 0;

      char legend00[30], legend01[30], legend02[60],legend03[30],legend04[30];

      for(unsigned int i=0; i<pCounts[idx]; i++)
        {
        // all protons
          shistEnergy_00->Fill(exitEnergies[idx][i]);
          shistAngle_00->Fill(anglesX[idx][i] * 180./TMath::Pi());
          rawcounter00++;

        // primaries only
        if(creatorProcess[idx][i]==0)
          {
          shistEnergy_01->Fill(exitEnergies[idx][i]);
          shistAngle_01->Fill(anglesX[idx][i] * 180./TMath::Pi());
          rawcounter01++;
          }

        // primaries without elastic scattering
        if(creatorProcess[idx][i]==0 && nuclearProcess[idx][i]==0)
          {
          shistEnergy_02->Fill(exitEnergies[idx][i]);
          shistAngle_02->Fill(anglesX[idx][i] * 180./TMath::Pi());
          rawcounter02++;
          }
        // all secondaries
        if(creatorProcess[idx][i]>0)
          {
          shistEnergy_03->Fill(exitEnergies[idx][i]);
          shistAngle_03->Fill(anglesX[idx][i] * 180./TMath::Pi());
          rawcounter03++;
          }

        if(anglesX[idx][i] <= pSumAngleSq [idx] &&
           anglesY[idx][i] <= pSumAngleSq [idx] &&
           vcl_abs((Ein-exitEnergies[idx][i])-pSumEnergy[idx]) <= pSumEnergySq[idx] )
          {
            shistEnergy_04->Fill(exitEnergies[idx][i]);
            shistAngle_04->Fill(anglesX[idx][i] * 180./TMath::Pi());
            cutcounter00++;
          }

        }

        std::cout<<rawcounter00<<"\t"<<rawcounter01<<"\t"<<rawcounter02<<"\t"<<rawcounter03<<"\t"<<std::endl;

        shistEnergy_00->Draw("");
        shistEnergy_00->SetLineColor(kBlack);
        shistEnergy_00->GetXaxis()->SetTitle("Exit energy (MeV)");
        shistEnergy_00->GetYaxis()->SetTitle("Number of protons");

        shistEnergy_01->Draw("same");
        shistEnergy_01->SetLineColor(kRed);

        shistEnergy_02->Draw("same");
        shistEnergy_02->SetLineColor(kGreen);

        shistEnergy_03->Draw("same");
        shistEnergy_03->SetLineColor(kBlue);

        shistEnergy_04->Draw("same");
        shistEnergy_04->SetLineColor(kYellow);


        sprintf( legend00, "All protons (%i)", rawcounter00);
        sprintf( legend01, "Primary protons (%.2f%%)", (float)rawcounter01/rawcounter00*100);
        sprintf( legend02, "Primary protons without nuclear interactions (%.2f%%)", (float)rawcounter02/rawcounter00*100);
        sprintf( legend03, "Secondary protons (%.2f%%)", (float)rawcounter03/rawcounter00*100);
        sprintf( legend04, "With cut (%.2f%%)", (float)cutcounter00/rawcounter00*100);

        TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
        leg->AddEntry(shistEnergy_00,legend00);
        leg->AddEntry(shistEnergy_01,legend01);
        leg->AddEntry(shistEnergy_02,legend02);
        leg->AddEntry(shistEnergy_03,legend03);
        leg->AddEntry(shistEnergy_04,legend04);
        leg->Draw();

        gPad->Print(sfilenameE);
        delete shistEnergy_00;
        delete shistEnergy_01;
        delete shistEnergy_02;
        delete shistEnergy_03;
        delete shistEnergy_04;
        delete leg;

        shistAngle_00->Draw("");
        shistAngle_00->SetLineColor(kBlack);
        shistAngle_00->GetXaxis()->SetTitle("Exit angle (deg)");
        shistAngle_00->GetYaxis()->SetTitle("Number of protons");

        shistAngle_01->Draw("same");
        shistAngle_01->SetLineColor(kRed);

        shistAngle_02->Draw("same");
        shistAngle_02->SetLineColor(kGreen);

        shistAngle_03->Draw("same");
        shistAngle_03->SetLineColor(kBlue);

        shistAngle_04->Draw("same");
        shistAngle_04->SetLineColor(kYellow);

        TLegend *legA = new TLegend(0.1,0.7,0.48,0.9);
        legA->AddEntry(shistAngle_00,legend00);
        legA->AddEntry(shistAngle_01,legend01);
        legA->AddEntry(shistAngle_02,legend02);
        legA->AddEntry(shistAngle_03,legend03);
        legA->AddEntry(shistAngle_04,legend04);
        legA->Draw();

        gPad->Print(sfilenameA);
        delete shistAngle_00;
        delete shistAngle_01;
        delete shistAngle_02;
        delete shistAngle_03;
        delete shistAngle_04;
        delete legA;
      }
    }

//===================================================================================================================================
// Plot superimposed distributions with cut
//===================================================================================================================================
  if(args_info.plotSuperDistribution_given)
    {
    for(unsigned int idx=0; idx<npixels; idx++)
      {
      char sPlotE[200], sfilenameE[100], sfoldernameE[100];
      sprintf( sPlotE, "Energy distributions with cut in pixel index %d", idx);
      sprintf( sfoldernameE,  "mkdir -p cutsplotE");  system(sfoldernameE);
      sprintf( sfilenameE, "cutsplotE/splotE_%d.png", idx);

      char sPlotA[200], sfilenameA[100], sfoldernameA[100];
      sprintf( sPlotA, "Angular distributions in pixel index %d", idx);
      sprintf( sfoldernameA,  "mkdir -p cutsplotA");  system(sfoldernameA);
      sprintf( sfilenameA, "cutsplotA/splotA_%d.png", idx);

      TH1F *shistEnergy_00 = new TH1F("shistEnergy_00", sPlotE, 420, 0., 210.);
      TH1F *shistEnergy_01 = new TH1F("shistEnergy_01", sPlotE, 420, 0., 210.);
      TH1F *shistEnergy_02 = new TH1F("shistEnergy_02", sPlotE, 420, 0., 210.);
      TH1F *shistEnergy_03 = new TH1F("shistEnergy_03", sPlotE, 420, 0., 210.);

      TH1F *shistAngle_00 = new TH1F("shistAngle_00", sPlotA, 300, 0., 60.);
      TH1F *shistAngle_01 = new TH1F("shistAngle_01", sPlotA, 300, 0., 60.);
      TH1F *shistAngle_02 = new TH1F("shistAngle_02", sPlotA, 300, 0., 60.);
      TH1F *shistAngle_03 = new TH1F("shistAngle_03", sPlotA, 300, 0., 60.);

      unsigned int cutcounter00 = 0;
      unsigned int cutcounter01 = 0;
      unsigned int cutcounter02 = 0;
      unsigned int cutcounter03 = 0;
      char legend00[30], legend01[30], legend02[60],legend03[30];

      for(unsigned int i=0; i<pCounts[idx]; i++)
        {
        if( !args_info.gauscut_flag &&
          anglesX[idx][i] <= pSumAngleSq [idx] &&
          anglesY[idx][i] <= pSumAngleSq [idx] &&
          vcl_abs((Ein-exitEnergies[idx][i])-pSumEnergy[idx]) <= pSumEnergySq[idx] )
          {

        // all protons
          shistEnergy_00->Fill(exitEnergies[idx][i]);
          shistAngle_00->Fill(anglesX[idx][i] * 180./TMath::Pi());
          cutcounter00++;

        // primaries only
        if(creatorProcess[idx][i]==0)
          {
          shistEnergy_01->Fill(exitEnergies[idx][i]);
          shistAngle_01->Fill(anglesX[idx][i] * 180./TMath::Pi());
          cutcounter01++;
          }

        // primaries without elastic scattering
        if(creatorProcess[idx][i]==0 && nuclearProcess[idx][i]==0)
          {
          shistEnergy_02->Fill(exitEnergies[idx][i]);
          shistAngle_02->Fill(anglesX[idx][i] * 180./TMath::Pi());
          cutcounter02++;
          }
        // all secondaries
        if(creatorProcess[idx][i]>0)
          {
          shistEnergy_03->Fill(exitEnergies[idx][i]);
          shistAngle_03->Fill(anglesX[idx][i] * 180./TMath::Pi());
          cutcounter03++;
          }
          }
        }

        std::cout<<cutcounter00<<"\t"<<cutcounter01<<"\t"<<cutcounter02<<"\t"<<cutcounter03<<std::endl;


        shistEnergy_00->Draw("");
        shistEnergy_00->SetLineColor(kBlack);
        shistEnergy_00->GetXaxis()->SetTitle("Exit energy (MeV)");
        shistEnergy_00->GetYaxis()->SetTitle("Number of protons");

        shistEnergy_01->Draw("same");
        shistEnergy_01->SetLineColor(kRed);

        shistEnergy_02->Draw("same");
        shistEnergy_02->SetLineColor(kGreen);

        shistEnergy_03->Draw("same");
        shistEnergy_03->SetLineColor(kBlue);

        sprintf( legend00, "All protons (%.2f%%)", (float)cutcounter00/rawcounter00*100);
        sprintf( legend01, "Primary protons (%.2f%%)", (float)cutcounter01/rawcounter00*100);
        sprintf( legend02, "Primary protons without nuclear interactions (%.2f%%)", (float)cutcounter02/rawcounter00*100);
        sprintf( legend03, "Secondary protons (%.2f%%)", (float)cutcounter03/rawcounter00*100);

        TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
        leg->AddEntry(shistEnergy_00,legend00);
        leg->AddEntry(shistEnergy_01,legend01);
        leg->AddEntry(shistEnergy_02,legend02);
        leg->AddEntry(shistEnergy_03,legend03);
        leg->Draw();

        gPad->Print(sfilenameE);
        delete shistEnergy_00;
        delete shistEnergy_01;
        delete shistEnergy_02;
        delete shistEnergy_03;
        delete leg;

        shistAngle_00->Draw("");
        shistAngle_00->SetLineColor(kBlack);
        shistAngle_00->GetXaxis()->SetTitle("Exit angle (deg)");
        shistAngle_00->GetYaxis()->SetTitle("Number of protons");

        shistAngle_01->Draw("same");
        shistAngle_01->SetLineColor(kRed);

        shistAngle_02->Draw("same");
        shistAngle_02->SetLineColor(kGreen);

        shistAngle_03->Draw("same");
        shistAngle_03->SetLineColor(kBlue);

        TLegend *legA = new TLegend(0.1,0.7,0.48,0.9);
        legA->AddEntry(shistAngle_00,legend00);
        legA->AddEntry(shistAngle_01,legend01);
        legA->AddEntry(shistAngle_02,legend02);
        legA->AddEntry(shistAngle_03,legend03);
        legA->Draw();

        gPad->Print(sfilenameA);
        delete shistAngle_00;
        delete shistAngle_01;
        delete shistAngle_02;
        delete shistAngle_03;
        delete legA;

      }
    }

//===================================================================================================================================
// biPlot cut distributions
//===================================================================================================================================

  if(args_info.plotCutDistribution_given)
    {
    gROOT->SetBatch(kTRUE);
    for(unsigned int idx=0; idx<npixels; idx++)
      {
      char cbiPlotX_00[200], filename_00[100], foldername_00[100];
      char cbiPlotX_01[200], filename_01[100], foldername_01[100];
      char cbiPlotX_02[200], filename_02[100], foldername_02[100];
      char cbiPlotX_03[200], filename_03[100], foldername_03[100];
      char cbiPlotX_04[200], filename_04[100], foldername_04[100];
      char cbiPlotX_05[200], filename_05[100], foldername_05[100];
      char cbiPlotX_06[200], filename_06[100], foldername_06[100];
      char cbiPlotX_07[200], filename_07[100], foldername_07[100];
      char cbiPlotX_08[200], filename_08[100], foldername_08[100];
      char cbiPlotX_09[200], filename_09[100], foldername_09[100];
      char cbiPlotX_10[200], filename_10[100], foldername_10[100];

      sprintf( cbiPlotX_00, "Biplot with cut of all protons in pixel index %d", idx);
      sprintf( cbiPlotX_01, "Biplot with cut of all primary protons in pixel index %d", idx);
      sprintf( cbiPlotX_02, "Biplot with cut of primary protons without elastic scattering in pixel index %d", idx);
      sprintf( cbiPlotX_03, "Biplot with cut of primary protons with elastic scattering in pixel index %d", idx);
      sprintf( cbiPlotX_04, "Biplot with cut of all secondary protons in pixel index %d", idx);
      sprintf( cbiPlotX_05, "Biplot with cut of secondary protons created from inelastic scattering in pixel index %d", idx);
      sprintf( cbiPlotX_06, "Biplot with cut of secondary protons created from inelastic scattering and without elastic interaction in pixel index %d", idx);
      sprintf( cbiPlotX_07, "Biplot with cut of secondary protons created from inelastic scattering and with elastic interaction in pixel index %d", idx);
      sprintf( cbiPlotX_08, "Biplot with cut of secondary protons created from elastic scattering in pixel index %d", idx);
      sprintf( cbiPlotX_09, "Biplot with cut of secondary protons created from elastic scattering and without elastic interaction in pixel index %d", idx);
      sprintf( cbiPlotX_10, "Biplot with cut of secondary protons created from elastic scattering and with elastic interaction in pixel index %d", idx);

      sprintf( foldername_00,  "mkdir -p cutbiplot_00");  system(foldername_00);
      sprintf( foldername_01,  "mkdir -p cutbiplot_01");  system(foldername_01);
      sprintf( foldername_02,  "mkdir -p cutbiplot_02");  system(foldername_02);
      sprintf( foldername_03,  "mkdir -p cutbiplot_03");  system(foldername_03);
      sprintf( foldername_04,  "mkdir -p cutbiplot_04");  system(foldername_04);
      sprintf( foldername_05,  "mkdir -p cutbiplot_05");  system(foldername_05);
      sprintf( foldername_06,  "mkdir -p cutbiplot_06");  system(foldername_06);
      sprintf( foldername_07,  "mkdir -p cutbiplot_07");  system(foldername_07);
      sprintf( foldername_08,  "mkdir -p cutbiplot_08");  system(foldername_08);
      sprintf( foldername_09,  "mkdir -p cutbiplot_09");  system(foldername_09);
      sprintf( foldername_10,  "mkdir -p cutbiplot_10");  system(foldername_10);

      sprintf( filename_00, "cutbiplot_00/biplotX_%d.png", idx);
      sprintf( filename_01, "cutbiplot_01/biplotX_%d.png", idx);
      sprintf( filename_02, "cutbiplot_02/biplotX_%d.png", idx);
      sprintf( filename_03, "cutbiplot_03/biplotX_%d.png", idx);
      sprintf( filename_04, "cutbiplot_04/biplotX_%d.png", idx);
      sprintf( filename_05, "cutbiplot_05/biplotX_%d.png", idx);
      sprintf( filename_06, "cutbiplot_06/biplotX_%d.png", idx);
      sprintf( filename_07, "cutbiplot_07/biplotX_%d.png", idx);
      sprintf( filename_08, "cutbiplot_08/biplotX_%d.png", idx);
      sprintf( filename_09, "cutbiplot_09/biplotX_%d.png", idx);
      sprintf( filename_10, "cutbiplot_10/biplotX_%d.png", idx);

      TH2F *hx_00 = new TH2F("hx_00", cbiPlotX_00, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_01 = new TH2F("hx_01", cbiPlotX_01, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_02 = new TH2F("hx_02", cbiPlotX_02, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_03 = new TH2F("hx_03", cbiPlotX_03, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_04 = new TH2F("hx_04", cbiPlotX_04, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_05 = new TH2F("hx_05", cbiPlotX_05, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_06 = new TH2F("hx_06", cbiPlotX_06, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_07 = new TH2F("hx_07", cbiPlotX_07, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_08 = new TH2F("hx_08", cbiPlotX_08, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_09 = new TH2F("hx_09", cbiPlotX_09, 210, 0, 210., 90, 0, 90.);
      TH2F *hx_10 = new TH2F("hx_10", cbiPlotX_10, 210, 0, 210., 90, 0, 90.);

      for(unsigned int i=0; i<pCounts[idx]; i++)
        {
        if( !args_info.gauscut_flag &&
          anglesX[idx][i] <= pSumAngleSq [idx] &&
          //anglesY[idx][i] <= pSumAngleSq [idx] &&
          vcl_abs((Ein-exitEnergies[idx][i])-pSumEnergy[idx]) <= pSumEnergySq[idx] )
          {
          // all protons (primaries and secondaries)
          hx_00->Fill(exitEnergies[idx][i],anglesX[idx][i] * 180./TMath::Pi());

          // all primaries
          if(creatorProcess[idx][i]==0)
            {
            hx_01->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
            }

          // primaries without elastic scattering
          if(creatorProcess[idx][i]==0 && nuclearProcess[idx][i]==0)
            {
            hx_02->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
            }

          // primaries with elastic scattering
          if(creatorProcess[idx][i]==0 && nuclearProcess[idx][i]==1)
            {
            hx_03->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
            }

          // all secondaries
          if(creatorProcess[idx][i]>0)
            {
            hx_04->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
            }

          // all secondaries created from inelastic scattering
          if(creatorProcess[idx][i]==2)
            {
            hx_05->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
            }

          // secondaries created from inelastic scattering and no elastic scattering
          if(creatorProcess[idx][i]==2 && nuclearProcess[idx][i]==0)
            {
            hx_06->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
            }

          // secondaries created from inelastic scattering and with elastic scattering
          if(creatorProcess[idx][i]==2 && nuclearProcess[idx][i]==1)
            {
            hx_07->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
            }

          // all secondaries created from elastic scattering
          if(creatorProcess[idx][i]==1)
            {
            hx_08->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
            }

          // all secondaries created from elastic scattering and without elastic scattering
          if(creatorProcess[idx][i]==1 && nuclearProcess[idx][i]==0)
            {
            hx_09->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
            }

          // all secondaries created from elastic scattering and with elastic scattering
          if(creatorProcess[idx][i]==1 && nuclearProcess[idx][i]==1)
            {
            hx_10->Fill(exitEnergies[idx][i],anglesX[idx][i]* 180./TMath::Pi());
            }
          }
        }

        hx_00->Draw("colz");
        hx_00->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_00->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_00);
        delete hx_00;

        hx_01->Draw("colz");
        hx_01->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_01->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_01);
        delete hx_01;

        hx_02->Draw("colz");
        hx_02->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_02->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_02);
        delete hx_02;

        hx_03->Draw("colz");
        hx_03->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_03->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_03);
        delete hx_03;

        hx_04->Draw("colz");
        hx_04->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_04->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_04);
        delete hx_04;

        hx_05->Draw("colz");
        hx_05->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_05->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_05);
        delete hx_05;

        hx_06->Draw("colz");
        hx_06->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_06->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_06);
        delete hx_06;

        hx_07->Draw("colz");
        hx_07->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_07->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_07);
        delete hx_07;

        hx_08->Draw("colz");
        hx_08->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_08->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_08);
        delete hx_08;

        hx_09->Draw("colz");
        hx_09->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_09->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_09);
        delete hx_09;

        hx_10->Draw("colz");
        hx_10->GetXaxis()->SetTitle("Exit energy (MeV)");
        hx_10->GetYaxis()->SetTitle("Exit angle (deg)");
        gPad->Print(filename_10);
        delete hx_10;
        }
      }

  std::cout << "Write pairs..." << std::endl;

  itk::ImageRegion<2> pairsRegion;
  itk::ImageRegion<2>::SizeType size;
  size[0] = region.GetSize(0);
  size[1] = pairs.size()/region.GetSize(0);

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

  if(args_info.sigmaAngle_given)
    {
    PWriterType::Pointer w = PWriterType::New();
    w->SetInput(sigmaAngle->GetOutput());
    w->SetFileName(args_info.sigmaAngle_arg);
    TRY_AND_EXIT_ON_ITK_EXCEPTION(w->Update());
    }
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
      std::cout << "Can not create energy.png, sigma is 0." << std::endl;
    else
      {
      RooRealVar rooEnergy("Energy","Energy loss",0.,500.*CLHEP::MeV, "MeV");
      RooDataSet rooEnergyData("data", "data", RooArgSet(rooEnergy));

      TH1F *h1 = new TH1F("exit energy distribution","",210, 0., 210.);
      TF1 *gausFit = new TF1 ("gausFit", "gaus");

      std::ofstream statdata_gaus, statdata_raw ;
      statdata_gaus.open ("statdata_gaus.txt",std::ofstream::app);
      statdata_raw.open ("statdata_raw.txt",std::ofstream::app);

      for(unsigned int i=0; i<pCounts[p]; i++)
        {
        rooEnergy = energies[p][i]/CLHEP::MeV;
        rooEnergyData.add(RooArgSet(rooEnergy));
        h1->Fill(exitEnergies[p][i]);
        }
      h1->Draw("");
      statdata_raw << h1->GetEntries() << "\t" << h1->GetMean() << "\t" << h1->GetRMS() << std::endl;
      h1->Fit("gausFit");

      float mean = gausFit->GetParameter(1);
      float sigma = gausFit->GetParameter(2);

      statdata_gaus << mean << "\t" << sigma << std::endl;
      gPad->SaveAs("exitenergy.png");

      TCanvas tEnergyCanvas("energyCanvas","Energy canvas",0,0,1000,500);
      RooPlot* rooEnergyPlot = rooEnergy.frame(pSumEnergy[p]-2*pSumEnergySq[p],
                                               pSumEnergy[p]+2*pSumEnergySq[p]);
      rooEnergyData.plotOn(rooEnergyPlot);

      RooRealVar rooMeanEnergy("rooMeanEnergy","mean of energy gaussian",pSumEnergy[p],0.,500.*CLHEP::MeV) ;
      RooRealVar rooSigmaEnergy("rooSigmaEnergy","width of energy gaussian",pSumEnergySq[p]/args_info.energycut_arg,0.,500.*CLHEP::MeV);
      RooGaussian rooGaussEnergy("rooGaussEnergy","energy gaussian pdf",rooEnergy,rooMeanEnergy,rooSigmaEnergy) ;

      rooGaussEnergy.plotOn(rooEnergyPlot, RooFit::LineColor(kBlue));
      //rooGaussEnergy.fitTo(rooEnergyData);
      //rooGaussEnergy.plotOn(rooEnergyPlot, LineColor(kRed));

      rooEnergyPlot->Draw();
      tEnergyCanvas.SaveAs("energy.png");
     }

    // Angle plot
    if(pSumAngleSq[p]==0.)
      std::cout << "Can not create angle.png, sigma is 0." << std::endl;
    else
      {
      TH1F *h2 = new TH1F("angular distribution","",200, -1, 1);
      RooRealVar rooAngle("Angle","Angle deviation",0.,2.*itk::Math::pi, "rad");
      RooDataSet rooAngleData("data", "data", RooArgSet(rooAngle));
      for(unsigned int i=0; i<angles[p].size(); i++)
        {
        rooAngle = angles[p][i];
        rooAngleData.add(RooArgSet(rooAngle));
        h2->Fill(angles[p][i]);
        }

      h2->Draw("");
      h2->Fit("gausFit");
      gPad->SaveAs("angulardist.png");

      TCanvas tAngleCanvas("angleCanvas","Angle canvas",0,0,1000,500);
      RooPlot* rooAnglePlot = rooAngle.frame(0., 2*pSumAngleSq[p]);
      rooAngleData.plotOn(rooAnglePlot);

      RooRealVar rooMeanAngle("rooMeanAngle","mean of angle gaussian",0.,0.,500.*CLHEP::MeV);
      RooRealVar rooSigmaAngle("rooSigmaAngle","width of angle gaussian",pSumAngleSq[p]/args_info.anglecut_arg,0.,500.*CLHEP::MeV);
      RooGaussian rooGaussAngle("rooGaussAngle","angle gaussian pdf",rooAngle,rooMeanAngle,rooSigmaAngle);

      rooGaussAngle.plotOn(rooAnglePlot, RooFit::LineColor(kBlue));
      //rooGaussAngle.fitTo(rooAngleData);
      //rooGaussAngle.plotOn(rooAnglePlot, LineColor(kRed));

      rooAnglePlot->Draw();
      tAngleCanvas.SaveAs("angle.png");
      }
    }

  return EXIT_SUCCESS;
}
