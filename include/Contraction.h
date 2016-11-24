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
            chi = 2,
            /// (complex conjugated) initial state
            phi = 3
        };

        /// Daughter WaveProducts
        std::vector<WaveProduct> WaveProducts;
        
        /// Square of coefficients
        RationalNumber SquaredCoefficient;

        /// exponents of gamma factors of contributing components
        std::vector<int> GammaExponents;
        
        /// Constructor
        Term(const WaveProduct& wp1, const WaveProduct& wp2, const WaveProduct& wp_chi, const WaveProduct& phi, const RationalNumber& c2)
        : WaveProducts({wp1, wp2, wp_chi, phi}), SquaredCoefficient(c2), GammaExponents(WaveProducts.size(), 0) {}

        /// contract components
        /// \param I vector of components to contract.
        /// if I is empty and only 3 components remain to be contracted, they are chosen
        void contract(std::vector<index> I = {});
    };

    /// Constructor
    /// \param psi coupled daughter state WaveFunction's
    /// \param chi orbital angular momentum WaveFunction
    Contraction(const CoupledWaveFunctions& psi, const OrbitalAngularMomentumWaveFunction& chi, const WaveFunction& phi);

    /// \return Terms_
    const std::vector<Term>& terms() const
    { return Terms_; }

    // prune off zero terms
    void prune()
    { Terms_.erase(std::remove_if(Terms_.begin(), Terms_.end(), [](const Term& t){return is_zero(t.SquaredCoefficient);}), Terms_.end()); }

    /// contract daughters
    void contract(const std::vector<Term::index>& I = {})
    { std::for_each(Terms_.begin(), Terms_.end(), [&I](Term& t){t.contract(I);}); prune(); }

    /// contract daughters
    void contract(const std::vector<Term::index>& I, unsigned n)
    { for (unsigned i = 0; i < n and !Terms_.empty(); ++i) contract(I); }

private:
    /// Terms in contraction
    std::vector<Term> Terms_;
};

namespace contractions {
    extern const std::vector<Contraction::Term::index> psi1_psi2;
    extern const std::vector<Contraction::Term::index> psi1_chi;
    extern const std::vector<Contraction::Term::index> psi2_chi;
    extern const std::vector<Contraction::Term::index> psi1_phi;
    extern const std::vector<Contraction::Term::index> psi2_phi;
    extern const std::vector<Contraction::Term::index> chi_phi;
}

/// \return vector of only the requested gamma exponents
const std::vector<int> gamma_exponents(const Contraction::Term& t, const std::vector<Contraction::Term::index>& I);

/// \typedef GammaPolynomial
 using GammaPolynomial = std::map<std::vector<int>, std::vector<std::pair<int, RationalNumber> > >;

/// add polynomials
GammaPolynomial& operator+=(GammaPolynomial& lhs, const GammaPolynomial& rhs);
 
/// add polynomials
const GammaPolynomial operator+(GammaPolynomial lhs, const GammaPolynomial& rhs)
{ return lhs += rhs; }
 
/// \return string of polynomial term
std::string to_string(const GammaPolynomial::value_type& ges_coefs);

/// \return string of gamma polynomial
std::string to_string(const GammaPolynomial& gp);
 
/// \return GammaPolynomial
const GammaPolynomial gamma_polynomial(const std::vector<Contraction::Term>& terms);
 
/// \return rank of contraction term
inline unsigned rank(const Contraction::Term& t)
{ return std::accumulate(t.WaveProducts.begin(), t.WaveProducts.end(), 0u, [](unsigned r, const WaveProduct& w){return r + rank(w);}); }

/// \return rank of contraction
inline unsigned rank(const Contraction& C)
{ return C.terms().empty() ? 0 : rank(C.terms()[0]); }

/// \return indices for triple contraction; empty if none possible
const std::vector<Contraction::Term::index> triple_contraction(const Contraction::Term& t);

/// \return indices for triple contraction; empty if none possible
inline const std::vector<Contraction::Term::index> triple_contraction(const Contraction& C)
{ if(C.terms().empty()) return std::vector<Contraction::Term::index>(); return triple_contraction(C.terms()[0]); }

/// convert to string
std::string to_string(const Contraction::Term::index i);

/// convert to string
std::string to_string(const std::vector<Contraction::Term::index>& I);
 
/// convert to string
std::string to_string(const Contraction::Term& term);

/// convert to string
std::string to_string(const std::vector<Contraction::Term>& terms);

}

#endif
