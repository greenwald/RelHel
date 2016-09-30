/// \file
/// \brief code adapted from YAP (www.github.com/YAP/YAP)
/// \author Daniel Greenwald

#ifndef __ClebschGordan__h
#define __ClebschGordan__h

#include "RationalNumber.h"

#include <string>

namespace relhel {

/// \namespace ClebschGordan
/// \brief Clebsch-Gordan coefficient calculation
/// \author Daniel Greenwald
namespace ClebschGordan {

/// \return Clebsch-Gordan coefficient string
/// \param two_j1 2*spin of first particle
/// \param two_m1 2*spin-projection of first particle
/// \param two_j2 2*spin of second particle
/// \param two_m2 2*spin-projection of second particle
/// \param two_J  2*spin of composite system
/// \param two_M  2*spin-projection of composite system
std::string to_string(unsigned two_j1, int two_m1, unsigned two_j2, int two_m2, unsigned two_J, int two_M);

/// \return Clebsch-Gordan coefficient string
/// \param two_j1 2*spin of first particle
/// \param two_m1 2*spin-projection of first particle
/// \param two_j2 2*spin of second particle
/// \param two_m2 2*spin-projection of second particle
/// \param two_J  2*spin of composite system
inline std::string to_string(unsigned two_j1, int two_m1, unsigned two_j2, int two_m2, unsigned two_J)
{ return to_string(two_j1, two_m1, two_j2, two_m2, two_J, two_m1 + two_m2); }

/// \return consistency of spin and spin projection
/// \param two_J 2*spin
/// \param two_M 2*spin-projection
const bool consistent(unsigned two_J, int two_M);

/// \return Whether Clebsch-Gordan coefficient is nonzero
/// \param two_j1 2*spin of first particle
/// \param two_m1 2*spin-projection of first particle
/// \param two_j2 2*spin of second particle
/// \param two_m2 2*spin-projection of second particle
/// \param two_J  2*spin of composite system
/// \param two_M  2*spin-projection of composite system
const bool nonzero(unsigned two_j1, int two_m1, unsigned two_j2, int two_m2, unsigned two_J, int two_M);

/// \return Whether Clebsch-Gordan coefficient is nonzero, with M := m1 + m2.
/// \param two_j1 2*spin of first particle
/// \param two_m1 2*spin-projection of first particle
/// \param two_j2 2*spin of second particle
/// \param two_m2 2*spin-projection of second particle
/// \param two_J  2*spin of composite system
inline const bool nonzero(unsigned two_j1, int two_m1, unsigned two_j2, int two_m2, unsigned two_J)
{ return nonzero(two_j1, two_m1, two_j2, two_m2, two_J, two_m1 + two_m2); }

/// Calculates _signed_ _squared_ C-G coefficienty. For return value
/// C, the Clebsch-Gordan coefficient is sign(C)sqrt(abs(C))
/// \return signed squared Clebsch-Gordan coefficient (j1 m1 j2 m2 | J M)
/// Implemented from Eq. (16) from G. Racah, "Theory of Complex Spectra. II", Phys. Rev. 62, 438 (1942)
/// \param two_j1 2*spin of first particle
/// \param two_m1 2*spin-projection of first particle
/// \param two_j2 2*spin of second particle
/// \param two_m2 2*spin-projection of second particle
/// \param two_J  2*spin of composite system
/// \param two_M  2*spin-projection of composite system
const RationalNumber squared_coefficient(unsigned two_j1, int two_m1, unsigned two_j2, int two_m2, unsigned two_J, int two_M);

/// Calculates _signed_ _squared_ C-G coefficienty. For return value
/// C, the Clebsch-Gordan coefficient is sign(C)sqrt(abs(C))
/// \return signed squared Clebsch-Gordan coefficient (j1 m1 j2 m2 | J M)
/// Implemented from Eq. (16) from G. Racah, "Theory of Complex Spectra. II", Phys. Rev. 62, 438 (1942)
/// \param two_j1 2*spin of first particle
/// \param two_m1 2*spin-projection of first particle
/// \param two_j2 2*spin of second particle
/// \param two_m2 2*spin-projection of second particle
/// \param two_J  2*spin of composite system
inline const RationalNumber squared_coefficient(unsigned two_j1, int two_m1, unsigned two_j2, int two_m2, unsigned two_J)
{ return squared_coefficient(two_j1, two_m1, two_j2, two_m2, two_J, two_m1 + two_m2); }

}

}

#endif
