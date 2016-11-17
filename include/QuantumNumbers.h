#ifndef relhel__SpinParity_h
#define relhel__SpinParity_h

#include "MathUtils.h"
#include "Parity.h"

#include <numeric>
#include <string>

namespace relhel {

class QuantumNumbers
{
public:
    /// constructor
    /// \param two_j twice the spin
    /// \param p the parity
    constexpr QuantumNumbers(unsigned two_j, Parity p) : TwoJ_(two_j), Parity_(p) {}

    /// constructor
    /// \param s string to get quantum numbers from (JP)
    explicit QuantumNumbers(const std::string& s) :
    QuantumNumbers(std::stoi(s.substr(0, s.length() - 1)), to_parity(s.back()))
    {}
        
    /// \return twice the spin
    const unsigned twoJ() const
    { return TwoJ_; }

    /// \return parity
    const Parity parity() const
    { return Parity_; }

private:
    /// twice the spin
    unsigned TwoJ_;

    /// Parity
    Parity Parity_;
};

/// \return whether spin is whole integer
constexpr bool is_boson(const QuantumNumbers& q)
{ return is_even(q.twoJ()); }

/// \return whether spin is half integer
constexpr bool is_fermion(const QuantumNumbers& q)
{ return is_odd(q.twoJ()); }

/// \return QuantumNumbers as string
inline std::string to_string(const QuantumNumbers& q)
{ return spin_to_string(q.twoJ()) + to_string(q.parity()); }

/// \return parity of a vector of QuantumNumbers
inline Parity parity(const std::vector<QuantumNumbers>& Q)
{ return std::accumulate(Q.begin(), Q.end(), Parity::positive, [](Parity& p, const QuantumNumbers& q){return p *= q.parity();}); }

/// \return sum of spins of quantum numbers
inline unsigned spin_sum(const std::vector<QuantumNumbers>& Q)
{ return std::accumulate(Q.begin(), Q.end(), 0u, [](unsigned two_j, const QuantumNumbers& q){return two_j + q.twoJ();}); }

}

#endif
