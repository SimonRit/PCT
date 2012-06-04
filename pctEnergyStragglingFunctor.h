#ifndef __pctEnergyStragglingFunctor_h
#define __pctEnergyStragglingFunctor_h

#ifdef pascal
#  undef pascal
#endif
#include "CLHEP/Units/PhysicalConstants.h"

namespace pct
{

/** \class EnergyStragglingFunctor
 * \brief Function to compute energy straggling
 * This is based on a fit of a pencil beam of 200 MeV protons in water
 */

namespace Functor
{

template< class TInput, class TOutput >
class EnergyStragglingFunctor
{
public:
  EnergyStragglingFunctor() {}
  ~EnergyStragglingFunctor() {}

  static TOutput GetValue(const TInput l)
    {
    static const double a0 =  4.901560e-02 * CLHEP::MeV;
    static const double a1 =  3.105142e-02 * CLHEP::MeV / CLHEP::mm;
    static const double a2 = -6.097551e-04 * CLHEP::MeV / CLHEP::mm2;
    static const double a3 =  6.781029e-06 * CLHEP::MeV / CLHEP::mm3;
    static const double a4 = -3.355508e-08 * CLHEP::MeV /(CLHEP::mm *CLHEP::mm3);
    static const double a5 =  6.119699e-11 * CLHEP::MeV /(CLHEP::mm2*CLHEP::mm3);
    return a0+l*(a1+l*(a2+l*(a3+l*(a4+l*a5))));
    }
};

} // end namespace Functor
} // end namespace pct

#endif
