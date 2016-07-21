#include "pctpairgeometry_ggo.h"

#include <rtkMacro.h>
#include <rtkGgoFunctions.h>

#include "pctBetheBlochFunctor.h"

#define PAIRS_IN_RAM 1000000

int main(int argc, char * argv[])
{
  GGO(pctpairgeometry, args_info); //RTK macro parsing options from .ggo file (rtkMacro.h)

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

  pct::Functor::IntegratedBetheBlochProtonStoppingPowerInverse<float, double> *convFunc;
  convFunc = new pct::Functor::IntegratedBetheBlochProtonStoppingPowerInverse<float, double>(args_info.ionpot_arg * CLHEP::eV, 600.*CLHEP::MeV, 0.1*CLHEP::keV);

  double mag = 0.;
  unsigned int count = 0;
  double detDist = -1.;
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
      //const VectorType dIn = it.Get();
      ++it;
      //const VectorType dOut = it.Get();
      ++it;
      const VectorType data = it.Get();
      ++it;

      if(detDist==-1.)
        detDist = pOut[2]-pIn[2];
      if(detDist != pOut[2]-pIn[2])
        {
        std::cerr << "Error, the distance between detectors is not constant" << std::endl;
        exit(EXIT_FAILURE);
        }

      const double eIn = data[0];
      const double eOut = data[1];
      double wepl = 0.;
      if(eIn==0.)
        wepl = eOut; // Directly read WEPL
      else
        wepl = convFunc->GetValue(eOut, eIn); // convert to WEPL

      if(wepl>args_info.weplmax_arg)
        continue;
      if(abs(pOut[0])>args_info.mindist_arg && pIn[0]!=0.)
        {
        mag += double(pOut[0])/pIn[0];
        count++;
        }
      if(abs(pOut[1])>args_info.mindist_arg && pIn[1]!=0.)
        {
        mag += double(pOut[1])/pIn[1];
        count++;
        }
      }
    }
  mag /= count;
  double s = detDist / (1.-mag);
  std::cout << "Used " << count << " values to compute a " << mag << " magnification factor." << std::endl;
  std::cout << "Found source to exit detector distance: " << s+detDist << std::endl;
  std::cout << "Found source to entrance detector distance: " << s << std::endl;
  return EXIT_SUCCESS;
}
