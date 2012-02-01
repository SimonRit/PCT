#ifndef __pctBetheBlochFunctor_h
#define __pctBetheBlochFunctor_h

#ifdef pascal
#  undef pascal
#endif
#include "CLHEP/Units/PhysicalConstants.h"
#include <itkImage.h>

namespace pct
{

/** \class BetheBlochProtonStoppingPower
 * \brief Proton stopping power according to the Bethe-Bloch equation
 * The function corresponds to:
 *  - F in [Li, MIC-NSS, 2003].
 *  - S in [Schulte, Med Phys, 2005].
 *  - S in [Penfold, Med Phys, 2009].
 * Note that they have all an error. r must be squared in the equation.
 * The CLHEP system of units is used.
 */

/** Physical constants */
static const double K = 4. * CLHEP::pi *
                        CLHEP::classic_electr_radius *
                        CLHEP::classic_electr_radius *
                        CLHEP::electron_mass_c2 *
                        3.343e+23 / CLHEP::cm3;

namespace Functor
{

template< class TInput, class TOutput >
class BetheBlochProtonStoppingPower
{
public:
  BetheBlochProtonStoppingPower() {}
  ~BetheBlochProtonStoppingPower() {}

  /** Actual computation. By default, the ionization potential of ionization
   * is the one of water given in [Li, MIC-NSS, 2003].
   */
  TOutput GetValue(const TInput e, const double I = 68.9984 * CLHEP::eV) const
    {
    TOutput betasq = CLHEP::proton_mass_c2/(e + CLHEP::proton_mass_c2);
    betasq = 1.-betasq*betasq;
    return K * (log(2.*CLHEP::electron_mass_c2/I * betasq/(1-betasq))-betasq) / betasq;
    }
};

/** \class IntegratedBetheBlochProtonStoppingPowerInverse
 * \brief Numerical for integral used in proton CT.
 */
template< class TInput, class TOutput, unsigned int VMaximumkeVEnergy=500000 >
class IntegratedBetheBlochProtonStoppingPowerInverse
{
public:
  IntegratedBetheBlochProtonStoppingPowerInverse()
    {
    // Create lookup table for integer values
    m_LUT[0] = 0.;
    for(unsigned int i=1; i<VMaximumkeVEnergy; i++)
      {
      static const double dE = CLHEP::keV;
      m_LUT[i] = m_LUT[i-1] + dE / m_S.GetValue((TOutput)i*CLHEP::keV);
      }
    }
  ~IntegratedBetheBlochProtonStoppingPowerInverse() {}

  /** Get the integral from 0. to e. */
  TOutput GetValue(const TInput e) const
    {
    // SR: linear interpolation instead?
    return m_LUT[itk::Math::Round<int,TInput>(e / CLHEP::keV)];
    }

  /** Get the integral from e1 to e2. */
  TOutput GetValue(const TInput e1, const TInput e2) const
    {
    // SR: linear interpolation instead?
    return m_LUT[itk::Math::Round<int,TInput>(e2 / CLHEP::keV)] -
           m_LUT[itk::Math::Round<int,TInput>(e1 / CLHEP::keV)];
    }

private:
  BetheBlochProtonStoppingPower<TOutput, TOutput> m_S;
  TOutput m_LUT[VMaximumkeVEnergy];
};
} // end namespace Functor
} // end namespace pct

#endif
