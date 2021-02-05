#include "pctdiffxsection_ggo.h"
#include "pctGeant4.h"

#include <rtkMacro.h>

#include <G4DynamicParticle.hh>
#include <G4Proton.hh>
#include <G4Track.hh>
#include <G4ProcessTable.hh>
#include <G4HadronElasticProcess.hh>
#include <G4HadronElastic.hh>
#include <G4NistManager.hh>
#include <G4HadronicProcessStore.hh>

int main(int argc, char * argv[])
{
  GGO(pctdiffxsection, args_info);

  // Init Geant4
  pct::pctGeant4 * pctG4( pct::pctGeant4::GetInstance() );
  pctG4 = pctG4->GetInstance(); // Avoid warning

  // Create proton
  G4ThreeVector aMomentum;
  aMomentum[0] = 0;
  aMomentum[1] = 0;
  aMomentum[2] = 1;
  G4DynamicParticle *aProton = new G4DynamicParticle(G4Proton::Proton(),
                                                     aMomentum,
                                                     args_info.energy_arg);

  // Create track
  G4ThreeVector aPosition;
  aPosition[0] = 0;
  aPosition[1] = 0;
  aPosition[2] = 0;
  G4Track aTrack(aProton, 0., aPosition);
  aTrack.SetMomentumDirection(aMomentum);
  aTrack.SetKineticEnergy(args_info.energy_arg);
  G4Step * aStep = new G4Step();
  aTrack.SetStep(aStep);

  // Select material and process
  G4Material* aMaterial = G4NistManager::Instance()->FindOrBuildMaterial(args_info.material_arg);
  G4VProcess *vHadProc = G4ProcessTable::GetProcessTable()->FindProcess("hadElastic", aProton->GetParticleDefinition());
  G4HadronElasticProcess *hadProc = dynamic_cast<G4HadronElasticProcess*>(vHadProc);

  // Let's roll
  G4HadronElastic he;
  std::vector<unsigned long> hist(args_info.nbins_arg, 0);
  for(int i =0; i<args_info.npart_arg; i++){
    // Select target nucleus
    G4Nucleus aNucleus;
    hadProc->GetCrossSectionDataStore()->SampleZandA(aProton, aMaterial, aNucleus);

    // Apply process
    G4HadFinalState * hfs = he.ApplyYourself(aTrack, aNucleus);
    G4ThreeVector mc = hfs->GetMomentumChange();

    // Go to detector
    mc = mc * args_info.sdd_arg * CLHEP::mm / mc[2];

    // Measure distance from center
    unsigned long pos = floor((mc[0]*mc[0]+mc[1]*mc[1])/(args_info.binw_arg*CLHEP::mm));
    if(pos<hist.size())
      hist[pos]++;
  }

  // Account for total cross section of the process
  G4HadronicProcessStore* store = G4HadronicProcessStore::Instance();
  store->RegisterParticle(hadProc,G4Proton::Proton());
  double xs = store->GetCrossSectionPerVolume(G4Proton::Proton(),
                                              args_info.energy_arg*CLHEP::MeV,
                                              hadProc,
                                              aMaterial);

  // Output result
  std::ofstream os(args_info.output_arg);
  for(unsigned int i=0; i<hist.size(); i++)
    os << (xs*hist[i])/args_info.npart_arg << std::endl;

  return EXIT_SUCCESS;
}
