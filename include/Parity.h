/// \file
/// \author Daniel Greenwald

#ifndef relhel__Parity_h
#define relhel__Parity_h

#include <stdexcept>
#include <string>

namespace relhel {

    /// Particle parity
    enum class Parity : int
    {
        negative = -1,
        positive = +1
    };

    /// convert parity to string
    inline std::string to_string(Parity p)
    { return (p == Parity::positive) ? "+" : "-"; }

    /// Convert char (+/-) to Parity
    inline const Parity to_parity(char p)
    { if (p == '+') return Parity::positive; if (p == '-') return Parity::negative; throw std::invalid_argument("parity must be \'+\' or \'-\'"); }

}

#endif
