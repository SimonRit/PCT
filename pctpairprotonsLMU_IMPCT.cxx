#include "pctpairprotonsLMU_IMPCT_ggo.h"

#include <iomanip>

#include <rtkGgoFunctions.h>
#include <itkImage.h>
#include <itkImageIterator.h>
#include <itkImageFileWriter.h>
#include <itksys/SystemTools.hxx>

// Root includes
#include <TChain.h>
#include <TROOT.h>

#define MAX_RUNS 4096

struct ParticleData
  {
  itk::Vector<float,3> position;
  itk::Vector<float,3> direction;
  float ekine;
  float wepl;
  int pdgID;
  int trackID;
  int eventID;
  int pbID;
  };
/*
struct ParticleInfo
  {
  int runID;
  int eventID;
  int trackID;
  int pdgID;
  int pbID;
  char name[256];
  };
*/

bool SetTreeBranch(TChain *tree, std::string branchName, void *add, bool mandatory=true)
{
  unsigned int found = 0;
  tree->SetBranchStatus(branchName.c_str(), 1, &found);
  if(!found)
    {
    if(mandatory)
      {
      std::cerr <<  "Could not load branch "
               << branchName << std::endl;
      exit(EXIT_FAILURE);
      }
    }
  else
    tree->SetBranchAddress(branchName.c_str(), add);
  return found;
}

void BranchParticleToPhaseSpace(/*struct ParticleInfo &pi, */struct ParticleData &pd, TChain *tree)
{
  tree->GetListOfBranches(); // force reading of chain

  SetTreeBranch(tree, "PhaseSpaceBranch", pd.position.GetDataPointer());

/*
  // WARNING: X and Z are purposely swap...
  //SetTreeBranch(tree, "X",  pd.position.GetDataPointer()+2);
  SetTreeBranch(tree, "Y",  pd.position.GetDataPointer());
  SetTreeBranch(tree, "Z",  pd.position.GetDataPointer()+1);
  SetTreeBranch(tree, "dX", pd.direction.GetDataPointer()+2);
  SetTreeBranch(tree, "dY", pd.direction.GetDataPointer());
  SetTreeBranch(tree, "dZ", pd.direction.GetDataPointer()+1);
  //SetTreeBranch(tree, "Time", &pd.time);*/
}

void WritePairs(const std::vector< std::pair<ParticleData,ParticleData> > &pairs, std::string fileName)
{
  itk::ImageRegion<2> region;
  itk::ImageRegion<2>::SizeType size;
  size[0] = 5;
  size[1] = pairs.size();
  region.SetSize(size);

  typedef itk::Vector<float,3> PixelType;
  typedef itk::Image<PixelType,2> ImageType;
  ImageType::Pointer img = ImageType::New();
  img->SetRegions(region);
  img->Allocate();

  itk::ImageRegionIterator<ImageType> it(img, region);
  PixelType eet;
  for(size_t i=0; i<pairs.size(); i++)
    {
    eet[0] = pairs[i].first.ekine;
    eet[1] = pairs[i].second.ekine;
    eet[2] = 0.;//pairs[i].second.time - pairs[i].first.time;

    it.Set( pairs[i].first.position );
    ++it;
    it.Set( pairs[i].second.position );
    ++it;
    it.Set( pairs[i].first.direction );
    ++it;
    it.Set( pairs[i].second.direction );
    ++it;
    it.Set( eet );
    ++it;
	/*std::cout << "In WritePairs: " << eet[0] << "	" << eet[1] << "	" << eet[2] << "	" 
				<< pairs[i].first.position[0] << "	" << pairs[i].first.position[1] << "	" << pairs[i].first.position[2] << "	"
				<< pairs[i].second.position[0] << "	" << pairs[i].second.position[1] << "	" << pairs[i].second.position[2] << "	"
				<< pairs[i].first.direction[0] << "	" << pairs[i].first.direction[1] << "	" << pairs[i].first.direction[2] << "	"
				<< pairs[i].second.direction[0] << "	" << pairs[i].second.direction[1] << "	" << pairs[i].second.direction[2] << "	"
				<< std::endl;*/    
    }

  // Write
  typedef itk::ImageFileWriter<  ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( fileName );
  writer->SetInput( img );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() );
}

int main(int argc, char * argv[])
{
  GGO(pctpairprotonsLMU_IMPCT, args_info); //RTK macro parsing options from .ggo file (rtkMacro.h)

  // Create root trees
  TChain *treeIn = new TChain("PhaseSpaceTree");
  TChain *treeOut = new TChain("PhaseSpaceTree");
  treeIn->AddFile(args_info.inputIn_arg);
  treeOut->AddFile(args_info.inputOut_arg);

  // Branch particles
//  struct ParticleInfo piIn, piOut;
  struct ParticleData pdIn, pdOut;
  BranchParticleToPhaseSpace(/*piIn, */pdIn, treeIn);
  BranchParticleToPhaseSpace(/*piOut, */pdOut, treeOut);

  // Init
  std::vector< std::vector< std::pair<ParticleData, ParticleData> > > pairs(MAX_RUNS);
  size_t nparticulesIn = treeIn->GetEntries();
  size_t nparticulesOut = treeOut->GetEntries();
//  std::cout << "George: Entries treeIn/treeOut " << nparticulesIn << "	"	<< nparticulesOut << std::endl;
  size_t iIn=0, iOut=0;
  int prevEventIDIn = -1;
  int prevEventIDOut = -1;
  std::cout << iIn << " particles of input phase space processed ("
            << 100*iIn/nparticulesIn << "%)" << std::flush;

  // Go over root files
  while(iIn<nparticulesIn && iOut<nparticulesOut)
    {
    if(iIn%1000000==0)
      std::cout << '\r' 
                << iIn << " particles of input phase space processed ("
                << 100*iIn/nparticulesIn << "%)"
                << std::flush;

    treeIn->GetEntry(iIn);
    treeOut->GetEntry(iOut);


/*	
	// George - test
    std::cout << "George: " << iIn << "	" << pdIn.position[0] << "	" << pdIn.position[1] << "	" << pdIn.position[2] << "	" 
				<< pdIn.direction[0] << "	" << pdIn.direction[1] << "	" << pdIn.direction[2] << "	" 
				<< pdIn.ekine << "	" << pdIn.pdgID << "	" << pdIn.trackID << "	" 
				<< pdIn.eventID << "	" << pdIn.pbID << "	" << iOut	<< "	" << pdOut.position[0] << "	" << pdOut.position[1] << "	" << pdOut.position[2]			
				<< std::endl;
*/
//    // Move to next adequate RunID
//    if(piIn.runID!=args_info.runid_arg)
//      {
//      while(piIn.runID!=args_info.runid_arg && ++iIn<nparticulesIn)
//        treeIn->GetEntry(iIn);
//      continue;
//      }
//    if(piOut.runID!=args_info.runid_arg)
//      {
//      while(piOut.runID!=args_info.runid_arg && ++iOut<nparticulesOut)
//        treeOut->GetEntry(iOut);
//      continue;
//      }
    
    
    // Manage merged root files
    if(pdIn.eventID<prevEventIDIn)
      {
      while(pdOut.eventID>=prevEventIDOut && ++iOut<nparticulesOut)
        {
        prevEventIDOut = pdOut.eventID;
        treeOut->GetEntry(iOut);
        }
      prevEventIDIn = pdIn.eventID;
      prevEventIDOut = pdOut.eventID;
      continue;
      }
    if(pdOut.eventID<prevEventIDOut)
      {
      while(pdIn.eventID>=prevEventIDIn && ++iIn<nparticulesIn)
        {
        prevEventIDIn = pdIn.eventID;
        treeIn->GetEntry(iIn);
        }
      prevEventIDIn = pdIn.eventID;
      prevEventIDOut = pdOut.eventID;
      continue;
      }
     
    prevEventIDIn = pdIn.eventID;
    prevEventIDOut = pdOut.eventID;
    

    // Condition 1: both particles must be protons
//    if( pdOut.pdgID != 2212 )
//      {
//      iOut++;
//      continue;
//      }
//    if( pdIn.pdgID != 2212 )
//      {
//      iIn++;
//      continue;
//      }

    // Condition 2: absolute time difference must be small
//    if( pIn.time-pOut.time<-100.f )
    if(pdIn.eventID < pdOut.eventID) //1 primary per event in Gate
      {
      iIn++;
      continue;
      }
//    if( pIn.time-pOut.time>100.f )
    if(pdOut.eventID < pdIn.eventID)
      {
      iOut++;
      continue;
      }
      
// In the itk coordinate system, the beam has to go from -z towards z
// In this sim it goes from -x to x
// x has to be fliped with z and y sign has to be inverted
// The phase space files for that sim are in cm!!!
	float tempIn_x = pdIn.position[0];
	float tempOut_x = pdOut.position[0];
	float tempIn_z = pdIn.position[2];
	float tempOut_z = pdOut.position[2];
	pdIn.position[0] = 10.*pdIn.position[1];
	pdOut.position[0] = 10.*pdOut.position[1];
    pdIn.position[2] = 10.*tempIn_x;
    pdOut.position[2] = 10.*tempOut_x;
    pdIn.position[1] = 10.*tempIn_z;
    pdOut.position[1] = 10.*tempOut_z;    

 	float tempIn_px = pdIn.direction[0];
	float tempOut_px = pdOut.direction[0];
	float tempIn_pz = pdIn.direction[2];
	float tempOut_pz = pdOut.direction[2];
	pdIn.direction[0] = pdIn.direction[1];
	pdOut.direction[0] = pdOut.direction[1];
    pdIn.direction[2] = tempIn_px;
    pdOut.direction[2] = tempOut_px;
    pdIn.direction[1] = tempIn_pz;
    pdOut.direction[1] = tempOut_pz; 


    // Corresponding protons found, add to vector if no nuclear interaction
    if(args_info.runID_arg>=args_info.minRun_arg && args_info.runID_arg<args_info.maxRun_arg/* &&
       (!(args_info.nonuclear_flag) || (pdIn.trackID == pdOut.trackID))*/) // Condition to remove nuclear events if flag activated
      {
			/*std::cout << "George1 iIn/iOut: " << iIn << "	" << iOut << "	" 
				<< pdIn.pbID << "	" << pdOut.pbID << "	" 
				<< pdIn.eventID << "	" << pdOut.eventID << "	"
				<< pdIn.position[0] << "	" << pdIn.position[1] << "	" << pdIn.position[2] << "	"
				<< pdOut.position[0] << "	" << pdOut.position[1] << "	" << pdOut.position[2] << "	"
				<< std::endl;    */  
      pairs[args_info.runID_arg].push_back( std::pair<ParticleData,ParticleData>(pdIn, pdOut) );
      }

    // There may be multiple protons to pair with an input proton so only
    // increment the output counter. Note that there may also be multiple input
    // protons associated to an output proton (secondary proton back scattering) but this is ignored.
    iOut++;

    }

  std::cout << "\r"
            << nparticulesIn << " particles of input phase space processed ("
            << 100 << "%)"
            << std::endl
            << "Writing..."
            << std::endl;

  for(unsigned int i=0; i<MAX_RUNS; i++)
    {
    if(pairs[i].size())
      {
      std::ostringstream os;
      os << itksys::SystemTools::GetFilenamePath(args_info.output_arg) << "/"
         << itksys::SystemTools::GetFilenameWithoutLastExtension(args_info.output_arg)
         << std::setw(4) << std::setfill ('0') << i
         << itksys::SystemTools::GetFilenameLastExtension(args_info.output_arg);
      std::cout << "Writing into file:" << os.str() << std::endl;
      WritePairs(pairs[i], os.str());
      }
    }
  return EXIT_SUCCESS;
}
