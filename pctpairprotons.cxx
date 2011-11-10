#include "pctpairprotons_ggo.h"

#include <rtkGgoFunctions.h>
#include <itkImage.h>
#include <itkImageIterator.h>
#include <itkImageFileWriter.h>

// Root includes
#include <TChain.h>
#include <TROOT.h>

struct ParticleData
  {
  float ekine;
  itk::Vector<float,3> position;
  itk::Vector<float,3> direction;
  float time;
  };

struct ParticleInfo
  {
  int runID;
  int eventID;
  int trackID;
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

void BranchParticleToPhaseSpace(struct ParticleInfo &pi, struct ParticleData &pd, TChain *tree)
{
  tree->GetListOfBranches(); // force reading of chain
  if(!SetTreeBranch(tree, "ParticleName", pi.name, false))
    strcpy(pi.name, "proton"); // If absent, assume that particles have been filtered
  SetTreeBranch(tree, "RunID", &pi.runID);
  SetTreeBranch(tree, "TrackID", &pi.trackID);
  SetTreeBranch(tree, "EventID", &pi.eventID);
  SetTreeBranch(tree, "Ekine", &pd.ekine);

  // WARNING: X and Z are purposely swap...
  //SetTreeBranch(tree, "X", pd.position.GetDataPointer()+2);
  SetTreeBranch(tree, "Y", pd.position.GetDataPointer()+1);
  SetTreeBranch(tree, "Z", pd.position.GetDataPointer());
  SetTreeBranch(tree, "dX", pd.direction.GetDataPointer()+2);
  SetTreeBranch(tree, "dY", pd.direction.GetDataPointer()+1);
  SetTreeBranch(tree, "dZ", pd.direction.GetDataPointer());
  SetTreeBranch(tree, "Time", &pd.time);
}

void WritePairs(const std::vector< std::pair<ParticleData,ParticleData> > &pairs, const char *fileName)
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
  GGO(pctpairprotons, args_info);

  // Create root trees
  TChain *treeIn = new TChain("PhaseSpace");
  TChain *treeOut = new TChain("PhaseSpace");
  treeIn->AddFile(args_info.inputIn_arg);
  treeOut->AddFile(args_info.inputOut_arg);

  // Branch particles
  struct ParticleInfo piIn, piOut;
  struct ParticleData pdIn, pdOut;
  BranchParticleToPhaseSpace(piIn, pdIn, treeIn);
  BranchParticleToPhaseSpace(piOut, pdOut, treeOut);
  pdIn.position[2]  = args_info.planeIn_arg;
  pdOut.position[2] = args_info.planeOut_arg;

  // Init
  std::vector< std::pair<ParticleData, ParticleData> > pairs;
  size_t nparticulesIn = treeIn->GetEntries();
  size_t nparticulesOut = treeOut->GetEntries();
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

    // Move to next adequate RunID
    if(piIn.runID!=args_info.runid_arg)
      {
      while(piIn.runID!=args_info.runid_arg && ++iIn<nparticulesIn)
        treeIn->GetEntry(iIn);
      continue;
      }
    if(piOut.runID!=args_info.runid_arg)
      {
      while(piOut.runID!=args_info.runid_arg && ++iOut<nparticulesOut)
        treeOut->GetEntry(iOut);
      continue;
      }

    // Manage merged root files
    if(piIn.eventID<prevEventIDIn)
      {
      while(piOut.eventID>=prevEventIDOut && ++iOut<nparticulesOut)
        {
        prevEventIDOut = piOut.eventID;
        treeOut->GetEntry(iOut);
        }
      prevEventIDIn = piIn.eventID;
      prevEventIDOut = piOut.eventID;
      continue;
      }
    if(piOut.eventID<prevEventIDOut)
      {
      while(piIn.eventID>=prevEventIDIn && ++iIn<nparticulesIn)
        {
        prevEventIDIn = piIn.eventID;
        treeIn->GetEntry(iIn);
        }
      prevEventIDIn = piIn.eventID;
      prevEventIDOut = piOut.eventID;
      continue;
      }
    prevEventIDIn = piIn.eventID;
    prevEventIDOut = piOut.eventID;

    // Condition 1: both particles must be protons
    if( std::string(piOut.name) != std::string("proton") )
      {
      iOut++;
      continue;
      }
    if( std::string(piIn.name) != std::string("proton") )
      {
      iIn++;
      continue;
      }

    // Condition 2: absolute time difference must be small
//    if( pIn.time-pOut.time<-100.f )
    if(piIn.eventID < piOut.eventID)
      {
      iIn++;
      continue;
      }
//    if( pIn.time-pOut.time>100.f )
    if(piOut.eventID < piIn.eventID)
      {
      iOut++;
      continue;
      }

    // Corresponding protons found, add to vector if no nuclear interaction
    if(piIn.trackID == piOut.trackID)
      {
      // WARNING: We have swap x and z, z sign must also be changed
      pdIn.direction[2] *= -1.;
      pdOut.direction[2] *= -1.;
      pairs.push_back( std::pair<ParticleData,ParticleData>(pdIn, pdOut) );
      }
    iIn++;
    iOut++;
    }

  std::cout << "\r"
            << nparticulesIn << " particles of input phase space processed ("
            << 100 << "%)"
            << std::endl
            << "Writing..."
            << std::endl;

  WritePairs(pairs, args_info.output_arg);
  return EXIT_SUCCESS;
}
