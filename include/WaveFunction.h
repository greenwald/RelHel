#ifndef relhel__WaveFunction_h
#define relhel__WaveFunction_h

/* #include "QuantumNumbers.h" */
/* #include "TJwfTensor.h" */

#include "RationalNumber.h"

#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <string>
#include <vector>

namespace relhel {

/// vector of (twice) the spin projection of unit-spin wave functions
using WaveProduct = std::vector<int>;

/// \return whether all projections are +1, 0, or -1
inline const bool valid(const WaveProduct& v)
{ return std::all_of(v.begin(), v.end(), [](int m){return std::abs(m) <= 1;}); }

/// \return projection of WaveProduct
inline const int projection(const WaveProduct& v)
{ return std::accumulate(v.begin(), v.end(), 0); }

/// \return number of projection-zero waves in WaveProduct
inline const unsigned zeroes(const WaveProduct& v)
{ return std::count(v.begin(), v.end(), 0); }

/// \return rank of WaveProduct
constexpr unsigned rank(const WaveProduct& w)
{ return w.size(); }

/// \return total spin of WaveProduct
constexpr unsigned spin(const WaveProduct& v)
{ return 2 * v.size(); }

/// \return coefficient of WaveProduct
const RationalNumber squared_coefficient(const WaveProduct& v);

/// convert to string
std::string to_string(const WaveProduct& wp);

/// terms in sum of WaveProduct's (without coefficients)
using WaveProductSum = std::vector<WaveProduct>;

/// \return rank of WaveProductSum
constexpr unsigned rank(const WaveProductSum& wps)
{ return wps.empty() ? 0 : rank(wps[0]); }

/// \return spin of WaveProductSum
constexpr unsigned spin(const WaveProductSum& wps)
{ return wps.empty() ? 0 : spin(wps[0]); }

/// convert to string
std::string to_string(const WaveProductSum& wps);

/// holds wave functions
class WaveFunction
{
public:
    /// map of spin projection to WaveProductSum
    using map_type = std::map<int, WaveProductSum>;
    
    /// Constructor
    /// \param two_j (twice) the spin of wave function
    explicit WaveFunction(unsigned two_j);

    /// \return Projections_
    const map_type& projections() const
    { return Projections_; }
    
protected:

    /// map of spin projection to WaveProductSum
    map_type Projections_;
};

/// \return spin of WaveFunction
const unsigned spin(const WaveFunction& wf);

/// \return rank of WaveFunction
const unsigned rank(const WaveFunction& wf);

/// convert to string
std::string to_string(const WaveFunction& wf);

/// holds orbital angular momentum wave function
struct OrbitalAngularMomentumWaveFunction : public WaveFunction {
    /// Constructor
    /// \param l orbital angular momentum
    OrbitalAngularMomentumWaveFunction(unsigned l) : WaveFunction(2 * l)
    {
        for (auto it = Projections_.begin(); it != Projections_.end(); )
            if (it->first != 0)
                it = Projections_.erase(it);
            else
                ++it;
    }
};

/// Coupling of two WaveFunction's
class CoupledWaveFunctions
{
public:
    
    /// Constructor
    CoupledWaveFunctions(const WaveFunction& phi1, const WaveFunction& phi2, unsigned two_S, int delta);

    /// \return array of coupled WaveFunction's
    const std::array<WaveFunction, 2>& phi() const
    { return Phi_; }

    /// \return total coupled spin
    const unsigned twoS() const
    { return TwoS_; }

    /// \return spin projection
    const int delta() const
    { return Delta_; }
    
private:

    /// First WaveFunction to couple
    std::array<WaveFunction, 2> Phi_;

    /// total coupled spin
    unsigned TwoS_;

    /// projection of total coupled spin
    int Delta_;
    
};

/// \return rank of coupled wave functions
inline unsigned rank(const CoupledWaveFunctions& cwf)
{ return std::accumulate(cwf.phi().begin(), cwf.phi().end(), 0u, [](unsigned r, const WaveFunction& w){return r + rank(w);}); }

}

#endif
