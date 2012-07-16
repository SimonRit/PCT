#include "pctGeant4.h"
#include <G4RunManager.hh>
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
//#include "RunAction.hh"

/// Unique static instance
pct::pctGeant4 *pct::pctGeant4::mSingleton = 0;
itk::FastMutexLock::Pointer pct::pctGeant4::m_Lock = itk::FastMutexLock::New();

//------------------------------------------------------------------------------
pct::pctGeant4 * pct::pctGeant4::GetInstance()
{
  m_Lock->Lock();
  if (mSingleton == 0) {
    mSingleton = new pctGeant4();
  }
  m_Lock->Unlock();
  return mSingleton;
}

//------------------------------------------------------------------------------
pct::pctGeant4::pctGeant4()
{
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;
  
  // set mandatory initialization classes
  DetectorConstruction* det;
  PrimaryGeneratorAction* prim;
  runManager->SetUserInitialization(det = new DetectorConstruction);
  runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserAction(prim = new PrimaryGeneratorAction(det));
  
  // set user action classes
  //runManager->SetUserAction(new RunAction(det,prim));
  
  //  G4String command = "/control/execute ";
  //  G4String fileName = argv[1];
  //  G4UImanager::GetUIpointer()->ApplyCommand(command+fileName);
}
  
