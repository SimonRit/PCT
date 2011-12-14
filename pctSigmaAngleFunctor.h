#ifndef __pctSigmaAngleFunctor_h
#define __pctSigmaAngleFunctor_h

#ifdef pascal
#  undef pascal
#endif
#include "CLHEP/Units/PhysicalConstants.h"

namespace pct
{

namespace Functor
{

/** \class SigmaAngle
 * \brief Compute sigma angle according to equation 8 in [Schulte, Med Phys, 2008].
 */

static const double invX0 = 1./(36.1*CLHEP::cm);

#define SCHULTE
#if defined(SCHULTE)
  // Table 1 in [Schulte, 2008]
  static const double aunit = 1./(CLHEP::MeV*CLHEP::MeV);
  static const double a0 =  7.457e-6       * aunit;
  static const double a1 =  4.548e-7  / 2. * aunit /  CLHEP::cm;
  static const double a2 = -5.777e-8  / 3. * aunit / (CLHEP::cm * CLHEP::cm);
  static const double a3 =  1.301e-8  / 4. * aunit / (CLHEP::cm * CLHEP::cm * CLHEP::cm);
  static const double a4 = -9.228e-10 / 5. * aunit / (CLHEP::cm * CLHEP::cm * CLHEP::cm * CLHEP::cm);
  static const double a5 =  2.687e-11 / 6. * aunit / (CLHEP::cm * CLHEP::cm * CLHEP::cm * CLHEP::cm * CLHEP::cm);
#elif defined(WILLIAMS)
  // Table 2, (a) in [Williams, 2004]
  static const double aunit = 1./(CLHEP::MeV*CLHEP::MeV);
  static const double a0 =  7.507e-4      * aunit;
  static const double a1 =  3.320e-5 / 2. * aunit / (CLHEP::cm);
  static const double a2 = -4.171e-7 / 3. * aunit / (CLHEP::cm * CLHEP::cm);
  static const double a3 =  4.488e-7 / 4. * aunit / (CLHEP::cm * CLHEP::cm * CLHEP::cm);
  static const double a4 = -3.739e-8 / 5. * aunit / (CLHEP::cm * CLHEP::cm * CLHEP::cm * CLHEP::cm);
  static const double a5 =  1.455e-9 / 6. * aunit / (CLHEP::cm * CLHEP::cm * CLHEP::cm * CLHEP::cm * CLHEP::cm);
#endif

template< class TInput, class TOutput=TInput >
class SigmaAngle
{
public:
  SigmaAngle() {}
  ~SigmaAngle() {}

  static TOutput GetValue(const TInput l)
    {    
    // Constant out of integral
#if defined(SCHULTE)
    double c = 13.6*CLHEP::MeV*(1+0.038*log(l*invX0));
#elif defined(WILLIAMS)
    double c = 13.6*CLHEP::MeV;
#endif
    c *= c * invX0;

    // Multiplied with polynomial integral
    return sqrt(c*l*(a0+l*(a1+l*(a2+l*(a3+l*(a4+l*a5))))));
    }
};
} // end namespace Functor
} // end namespace pct

#endif
