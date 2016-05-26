#include "pctpairprotonsLomaLinda_ggo.h"

#include <iomanip>

#include <rtkGgoFunctions.h>
#include <itkImage.h>
#include <itkImageIterator.h>
#include <itkImageFileWriter.h>
#include <itksys/SystemTools.hxx>

// Root includes
#include <TChain.h>
#include <TROOT.h>
#include <TVector3.h>

#define MAX_RUNS 4096

struct ParticleData
  {
  float wepl;
  itk::Vector<float,3> position0;
  itk::Vector<float,3> position1;
  itk::Vector<float,3> position2;
  itk::Vector<float,3> position3;
  };


struct ParticleDataFinal
  {
  float wepl;
  itk::Vector<float,3> position;
  itk::Vector<float,3> direction;
  float time;
  };

struct ParticleInfo
  {
  int runID;
  char name[256];
  };


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

void BranchParticleToPhaseSpace(struct ParticleInfo &piInput, struct ParticleData &pdInput, TChain *tree)
{
  std::cout << "George: In BranchParticleToPhaseSpace" << std::endl;  	
  tree->GetListOfBranches(); // force reading of chain
  if(!SetTreeBranch(tree, "ParticleName", piInput.name, false))
    strcpy(piInput.name, "proton"); // If absent, assume that particles have been filtered
//  SetTreeBranch(tree, "RunID", &piInput.runID);
//  SetTreeBranch(tree, "TrackID", &pi.trackID);
//  SetTreeBranch(tree, "EventID", &pi.eventID);
  SetTreeBranch(tree, "calculated_WEPL", &pdInput.wepl);


  // WARNING: X and Z are purposely swap...
  SetTreeBranch(tree, "v_hit0",  pdInput.position0.GetDataPointer()+1);
  SetTreeBranch(tree, "t_hit0",  pdInput.position0.GetDataPointer());
  SetTreeBranch(tree, "u_hit0",  pdInput.position0.GetDataPointer()+2);
  SetTreeBranch(tree, "v_hit1",  pdInput.position1.GetDataPointer()+1);
  SetTreeBranch(tree, "t_hit1",  pdInput.position1.GetDataPointer());
  SetTreeBranch(tree, "u_hit1",  pdInput.position1.GetDataPointer()+2);
  SetTreeBranch(tree, "v_hit2",  pdInput.position2.GetDataPointer()+1);
  SetTreeBranch(tree, "t_hit2",  pdInput.position2.GetDataPointer());
  SetTreeBranch(tree, "u_hit2",  pdInput.position2.GetDataPointer()+2);
  SetTreeBranch(tree, "v_hit3",  pdInput.position3.GetDataPointer()+1);
  SetTreeBranch(tree, "t_hit3",  pdInput.position3.GetDataPointer());
  SetTreeBranch(tree, "u_hit3",  pdInput.position3.GetDataPointer()+2);
//  SetTreeBranch(tree, "Z",  pd.position1.GetDataPointer()+1);
//  SetTreeBranch(tree, "dX", pd.direction1.GetDataPointer()+2);
//  SetTreeBranch(tree, "dY", pd.direction1.GetDataPointer());
//  SetTreeBranch(tree, "dZ", pd.direction1.GetDataPointer()+1);
//  SetTreeBranch(tree, "Time", &pd.time);

}


void WritePairs(const std::vector< std::pair<ParticleDataFinal,ParticleDataFinal> > &pairs, std::string fileName)
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
    eet[0] = pairs[i].first.wepl;
    eet[1] = pairs[i].second.wepl;
    eet[2] = pairs[i].second.time - pairs[i].first.time;

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
  GGO(pctpairprotonsLomaLinda, args_info); //RTK macro parsing options from .ggo file (rtkMacro.h)

  // Create root trees
  TChain *treeIn = new TChain("recoENTRY");
  treeIn->AddFile(args_info.input_arg);
  std::cout << "George: Reading in file:" << args_info.input_arg << std::endl;

  // Branch particles
  struct ParticleInfo pi;
  struct ParticleData pd;
  struct ParticleDataFinal pdIn;
  struct ParticleDataFinal pdOut;

  BranchParticleToPhaseSpace(pi, pd, treeIn);

  // Init
  std::vector< std::vector< std::pair<ParticleDataFinal, ParticleDataFinal> > > pairs(MAX_RUNS);
  size_t nparticulesIn = treeIn->GetEntries();
  std::cout << "Number of entries = " << nparticulesIn << std::endl;
  size_t iIn=0;
  int prevEventIDIn = -1;
  std::cout << iIn << " particles of input phase space processed ("
            << 100*iIn/nparticulesIn << "%)" << std::flush;
  // Go over root files
  while(iIn<nparticulesIn)
    {
    if(iIn%10000==0)
      std::cout << '\r'
                << iIn << " particles of input phase space processed ("
                << 100*iIn/nparticulesIn << "%)"
                << std::flush;

    treeIn->GetEntry(iIn);

#if 0
	std::cout << "George: Input Entry = " << iIn << " with data : "
	<< pd.position0[0] << "	" << pd.position0[1] << "	" << pd.position0[2] << "	"
	<< pd.position1[0] << "	" << pd.position1[1] << "	" << pd.position1[2] << "	"
	<< pd.position2[0] << "	" << pd.position2[1] << "	" << pd.position2[2] << "	" 	
	<< pd.position3[0] << "	" << pd.position3[1] << "	" << pd.position3[2] << "	"	
	<< "	" << pd.wepl << std::endl;
#endif	

// Convert info from ParticleData to info for ParticleDataFinal
// The inner tracker planes are used
	pdIn.position = pd.position1;
	pdOut.position = pd.position2;
	pdIn.direction = pd.position1 - pd.position0;
	pdIn.direction /= pdIn.direction.GetNorm();
	pdOut.direction = pd.position3 - pd.position2;
	pdOut.direction /= pdOut.direction.GetNorm();
	pdIn.wepl = 0.;
	pdOut.wepl = pd.wepl;
	pdIn.time = 0.;
	pdOut.time = 0.;

#if 0
	std::cout << "George: Output Entry = " << iIn << " with data : "
	<< pdIn.position[0] << "	" << pdIn.position[1] << "	" << pdIn.position[2] << "	"
	<< pdOut.position[0] << "	" << pdOut.position[1] << "	" << pdOut.position[2] << "	"
	<< pdIn.direction[0] << "	" << pdIn.direction[1] << "	" << pdIn.direction[2] << "	" 	
	<< pdOut.direction[0] << "	" << pdOut.direction[1] << "	" << pdOut.direction[2] << "	"	
	<< "	" << pdIn.wepl << "	" << pdOut.wepl << std::endl; 	
#endif


    // Ensuring runID is meaningful and pairing protons
    if(pi.runID>=args_info.minRun_arg && pi.runID<args_info.maxRun_arg && pdOut.wepl>=-4. && pdOut.wepl<=300.)
      {
      pairs[args_info.runID_arg].push_back( std::pair<ParticleDataFinal,ParticleDataFinal>(pdIn, pdOut) );
      }

	iIn++;
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
      os << itksys::SystemTools::GetFilenameWithoutLastExtension(args_info.output_arg)
         << std::setw(4) << std::setfill ('0') << i
         << itksys::SystemTools::GetFilenameLastExtension(args_info.output_arg);
      WritePairs(pairs[i], os.str());
      }
    }


  return EXIT_SUCCESS;
}
