#ifndef relhel__Contraction_h
#define relhel__Contraction_h

#include "RationalNumber.h"
#include "WaveFunction.h"

#include <array>
#include <string>
#include <vector>

namespace relhel {

/// performs contractions over tensors
class Contraction
{
public:

    /// Contraction Term
    struct Term {

        /// index for WaveProducts
        enum index : size_t {
            /// first daughter in final state
            psi1 = 0,
            /// second daughter in final state
            psi2 = 1,
            /// orbital angular momentum
            chi = 2
        };

        /// Daughter WaveProducts
        std::vector<WaveProduct> WaveProducts;
        
        /// Square of coefficients
        RationalNumber SquaredCoefficient;

        /// exponents of gamma factors of contributing components
        std::vector<int> GammaExponents;
        
        /// Constructor
        Term(const WaveProduct& wp1, const WaveProduct& wp2, const WaveProduct& wp_chi, const RationalNumber& c2)
        : WaveProducts({wp1, wp2, wp_chi}), SquaredCoefficient(c2), GammaExponents(WaveProducts.size(), 0) {}

        // contract two components
        void contract(const std::array<index, 2>& I);
    };

    /// Constructor
    /// \param psi coupled daughter state WaveFunction's
    /// \param chi orbital angular momentum WaveFunction
    Contraction(const CoupledWaveFunctions& psi, const WaveFunction& chi);

    /// \return Terms_
    const std::vector<Term>& terms() const
    { return Terms_; }

    // prune off zero terms
    void prune()
    { Terms_.erase(std::remove_if(Terms_.begin(), Terms_.end(), [](const Term& t){return is_zero(t.SquaredCoefficient);}), Terms_.end()); }

    /// contract daughters
    void contract(const std::array<Term::index, 2>& I)
    { std::for_each(Terms_.begin(), Terms_.end(), [&I](Term& t){t.contract(I);}); prune(); }
    
private:
    /// Terms in contraction
    std::vector<Term> Terms_;
};

/// convert to string
std::string to_string(const Contraction::Term& term);

/// convert to string
std::string to_string(const std::vector<Contraction::Term>& terms);
 
}

#endif
