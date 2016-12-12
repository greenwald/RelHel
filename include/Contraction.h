#ifndef relhel__Contraction_h
#define relhel__Contraction_h

#include "RationalNumber.h"
#include "WaveFunction.h"

#include <algorithm>
#include <array>
#include <set>
#include <string>
#include <vector>

namespace relhel {

enum GammaIndex : size_t {
    /// first daughter in final state
    psi1 = 0,
    /// second daughter in final state
    psi2 = 1,
    /// orbital angular momentum
    chi = 2,
    /// (complex conjugated) initial state
    phi = 3,
    /// number of indices
    n_indices = 4    
};

std::string to_index_string(GammaIndex i);
std::string to_string(GammaIndex i);

/// Stores a polynomial of gamma factors
class GammaPolynomial
{
public:

    /// term in the gamma polynomial
    struct Term {
        /// coefficient
        std::vector<RationalNumber> SquaredCoefficients{std::vector<RationalNumber>(1, RationalNumber(1))};

        /// exponents of gamma terms
        std::vector<unsigned> Exponents{std::vector<unsigned>(GammaIndex::n_indices, 0)};
    };
    
    /// Constructor
    GammaPolynomial(const std::vector<GammaIndex>& I);

    /// \return Indices_
    const std::vector<GammaIndex>& indices() const
    { return Indices_; }

    const std::vector<Term>& terms() const
    { return Terms_; }
    
    GammaPolynomial& operator+=(const Term& term);
    
private:
    /// Indices to care about when building polynomial
    std::vector<GammaIndex> Indices_;

    /// Terms
    std::vector<Term> Terms_;
};

using EpsilonContraction = std::vector<GammaIndex>;
using ContractionMatrix = std::vector<std::vector<unsigned> >;

/// \return possible epsilon contractions
std::vector<EpsilonContraction> epsilon_contractions(const std::vector<unsigned>& r);

/// \return possible contraction matrices
std::set<ContractionMatrix> contraction_matrices(std::vector<unsigned> r, const std::vector<unsigned>& b);

/// \return possible contractions
std::map<EpsilonContraction, std::set<ContractionMatrix> > contractions(const CoupledWaveFunctions& psi,
                                                                        const OrbitalAngularMomentumWaveFunction& chi,
                                                                        const WaveFunction& phi);

/// convert to string
std::string to_string(const std::map<EpsilonContraction, std::set<ContractionMatrix> >& bC);

/// \return GammaPolynomial for contraction
const GammaPolynomial contract(const CoupledWaveFunctions& psi, const OrbitalAngularMomentumWaveFunction& chi,
                               const WaveFunction& phi, const std::vector<GammaIndex>& I = {GammaIndex::psi1, GammaIndex::psi2});

/// check exponents
const bool same_exponents(const GammaPolynomial::Term& A, const GammaPolynomial::Term& B, const std::vector<GammaIndex>& I)
{ return std::all_of(I.begin(), I.end(), [&](GammaIndex i){return A.Exponents[i] == B.Exponents[i];}); }

/// multiple two terms
GammaPolynomial::Term& operator*=(GammaPolynomial::Term& lhs, const GammaPolynomial::Term& rhs);

/// \return string of gamma polynomial term
std::string to_string(const GammaPolynomial::Term& t, const std::vector<GammaIndex>& I);

/// \return string of gamma polynomial
std::string to_string(const GammaPolynomial& gp);
 
}

#endif
