#include "pctpairprotons_ggo.h"

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
  float ekine;
  itk::Vector<float,3> position;
  itk::Vector<float,3> direction;
  };

struct StoredParticleInfo
  {
  int trackID;
  int nuclearProcess;
  int creatorProcess;
  int order;
  };

struct ParticleInfo
  {
  int runID;
  int eventID;
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

void BranchParticleToPhaseSpace(struct ParticleInfo &pi, struct StoredParticleInfo &spi, struct ParticleData &pd, TChain *tree, args_info_pctpairprotons *args_info)
{
  tree->GetListOfBranches(); // force reading of chain
  if(!SetTreeBranch(tree, "ParticleName", pi.name, false))
    strcpy(pi.name, "proton"); // If absent, assume that particles have been filtered
  SetTreeBranch(tree, "RunID", &pi.runID);
  SetTreeBranch(tree, "EventID", &pi.eventID);
  SetTreeBranch(tree, "Ekine", &pd.ekine);

  SetTreeBranch(tree, args_info->proju_arg,  pd.position.GetDataPointer());
  SetTreeBranch(tree, args_info->projv_arg,  pd.position.GetDataPointer()+1);
  SetTreeBranch(tree, std::string("d")+std::string(args_info->proju_arg), pd.direction.GetDataPointer());
  SetTreeBranch(tree, std::string("d")+std::string(args_info->projv_arg), pd.direction.GetDataPointer()+1);
  SetTreeBranch(tree, std::string("d")+std::string(args_info->projw_arg), pd.direction.GetDataPointer()+2);
  SetTreeBranch(tree, "TrackID", &spi.trackID);

  if(!SetTreeBranch(tree, "NuclearProcess", &spi.nuclearProcess, false))
    spi.nuclearProcess = -1;
  if(!SetTreeBranch(tree, "CreatorProcess", &spi.creatorProcess, false))
    spi.creatorProcess = -1;
  if(!SetTreeBranch(tree, "Order", &spi.order, false))
    spi.order = -1;
}

void WritePairs(const std::vector< std::pair<ParticleData,ParticleData> > &pairs,
                const std::vector< StoredParticleInfo> &particlesInfo,
                std::string fileName)
{
  itk::ImageRegion<2> region;
  itk::ImageRegion<2>::SizeType size;
  if(particlesInfo.back().nuclearProcess == -1)
    size[0] = 5;
  else
    size[0] = 6;
  size[1] = pairs.size();
  region.SetSize(size);

  typedef itk::Vector<float,3> PixelType;
  typedef itk::Image<PixelType,2> ImageType;
  ImageType::Pointer img = ImageType::New();
  img->SetRegions(region);
  img->Allocate();

  itk::ImageRegionIterator<ImageType> it(img, region);
  PixelType eet;
  PixelType nuclearinfo;

  for(size_t i=0; i<pairs.size(); i++)
    {
    eet[0] = pairs[i].first.ekine;
    eet[1] = pairs[i].second.ekine;
    eet[2] = particlesInfo[i].trackID;

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
    if(size[0] == 6)
      {
      eet[0] = particlesInfo[i].creatorProcess;
      eet[1] = particlesInfo[i].nuclearProcess;
      eet[2] = particlesInfo[i].order;

      it.Set( nuclearinfo );
      ++it;
      }
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
  GGO(pctpairprotons, args_info); //RTK macro parsing options from .ggo file (rtkMacro.h)

  // Create root trees
  TChain *treeIn = new TChain(args_info.psin_arg);
  TChain *treeOut = new TChain(args_info.psout_arg);
  treeIn->AddFile(args_info.inputIn_arg);
  treeOut->AddFile(args_info.inputOut_arg);

  // Branch particles
  struct ParticleInfo piIn, piOut;
  struct StoredParticleInfo spiIn, spiOut;
  struct ParticleData pdIn, pdOut;
  BranchParticleToPhaseSpace(piIn, spiIn, pdIn, treeIn, &args_info);
  BranchParticleToPhaseSpace(piOut, spiOut, pdOut, treeOut, &args_info);
  pdIn.position[2]  = args_info.planeIn_arg;
  pdOut.position[2] = args_info.planeOut_arg;

  // Init
  std::vector< std::vector< std::pair<ParticleData, ParticleData> > > pairs(MAX_RUNS);
  std::vector< std::vector< StoredParticleInfo > > particlesInfo(MAX_RUNS);
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
    if(piIn.eventID < piOut.eventID) //1 primary per event in Gate
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
    if(piIn.runID>=args_info.minRun_arg && piIn.runID<args_info.maxRun_arg &&
       (!(args_info.nonuclear_flag) || (spiIn.trackID == spiOut.trackID))) // Condition to remove nuclear events if flag activated
      {
      // WARNING: We have swap x and z, z sign must also be changed
      pdIn.direction[2] *= args_info.wweight_arg;
      pdOut.direction[2] *=  args_info.wweight_arg;;
      pairs[piIn.runID].push_back( std::pair<ParticleData,ParticleData>(pdIn, pdOut) );
      particlesInfo[piIn.runID].push_back( spiOut );
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
      os << itksys::SystemTools::GetFilenameWithoutLastExtension(args_info.output_arg)
         << std::setw(4) << std::setfill ('0') << i
         << itksys::SystemTools::GetFilenameLastExtension(args_info.output_arg);
      WritePairs(pairs[i], particlesInfo[i], os.str());
      }
    }
  return EXIT_SUCCESS;
}
