#ifndef relhel__SpinParity_h
#define relhel__SpinParity_h

#include "MathUtils.h"
#include "Parity.h"

#include <string>

namespace relhel {

class QuantumNumbers
{
public:
    /// constructor
    /// \param two_j twice the spin
    /// \param p the parity
    constexpr QuantumNumbers(int two_j, Parity p) : TwoJ_(two_j), Parity_(p) {}

    explicit QuantumNumbers(const std::string& s) : QuantumNumbers(std::stoi(s.substr(0, s.length() - 1)), to_parity(s.back())) {}
        
    /// \return twice the spin
    const int twoJ() const
    { return TwoJ_; }

    /// \return parity
    const Parity parity() const
    { return Parity_; }

private:
    /// twice the spin
    int TwoJ_;

    /// Parity
    Parity Parity_;
};

/// \return QuantumNumbers as string
std::string to_string(const QuantumNumbers& q)
{ return spin_to_string(q.twoJ()) + to_string(q.parity()); }

}

#endif
