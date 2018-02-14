#ifndef __pctGeant4_h
#define __pctGeant4_h

#include <itkFastMutexLock.h>

class G4RunManager;
class PhysicsListMessenger;

namespace pct
{

class pctGeant4 {
public:

  /// Get or build unique instance with this method
  static pctGeant4 * GetInstance();

private:
  /// Constructor
  pctGeant4();

  /// Destructor
  ~pctGeant4();

  /// Singleton object pointer i
  static pctGeant4 * mSingleton;

  static itk::FastMutexLock::Pointer m_Lock;
  G4RunManager * m_RunManager;
  PhysicsListMessenger* m_Mess;
};

} // end namespace pct

#endif
