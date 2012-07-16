#ifndef __pctGeant4_h
#define __pctGeant4_h

#include <itkFastMutexLock.h>

namespace pct
{

class pctGeant4 {
public:

  /// Get or build unique instance with this method
  static pctGeant4 * GetInstance();

private:
  /// Constructor
  pctGeant4();

  /// Singleton object pointer i
  static pctGeant4 * mSingleton;
  
  static itk::FastMutexLock::Pointer m_Lock;
};

} // end namespace pct

#endif
