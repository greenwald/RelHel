#ifndef relhel__WaveFunction_h
#define relhel__WaveFunction_h

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
inline const unsigned rank(const WaveProduct& w)
{ return w.size(); }

/// \return total spin of WaveProduct
inline const unsigned spin(const WaveProduct& v)
{ return 2 * v.size(); }

/// \return coefficient of WaveProduct
const RationalNumber squared_coefficient(const WaveProduct& v);

/// convert to string
std::string to_string(const WaveProduct& wp);

/// terms in sum of WaveProduct's (without coefficients)
using WaveProductSum = std::vector<WaveProduct>;

/// \return rank of WaveProductSum
inline const unsigned rank(const WaveProductSum& wps)
{ return wps.empty() ? 0 : rank(wps[0]); }

/// \return spin of WaveProductSum
inline const unsigned spin(const WaveProductSum& wps)
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

    using product_type = std::vector<WaveProduct>;

    /// Constructor
    CoupledWaveFunctions(const WaveFunction& A, const WaveFunction& B, unsigned two_s, int two_m);

    /// \return Projections_
    const std::vector<product_type>& products() const
    { return Products_; }

    /// \return spin
    const int spin() const
    { return TwoS_; }

private:

    /// spin
    unsigned TwoS_;

    /// vector of products
    std::vector<product_type> Products_;
};

/// \return rank of coupled wave functions
inline unsigned rank(const CoupledWaveFunctions& cwf, size_t i)
{ return (cwf.products().empty() or cwf.products()[0].empty()) ? 0 : rank(cwf.products()[0][i]); }

/// \return rank of coupled wave functions
inline unsigned rank(const CoupledWaveFunctions& cwf)
{ return rank(cwf, 0) + rank(cwf, 1); }


}

#endif
