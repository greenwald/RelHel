#ifndef TLSCONTRIB_HH
#define TLSCONTRIB_HH

class TLSAmpl;

#include "JLS.h"
#include "TFracNum.h"

#include <algorithm>
#include <array>
#include <vector>

struct PolynomialTerm {

    TFracNum PrefactorSquared;
    std::array<int, 2> Exponents;

    /// constructor
    PolynomialTerm(const TFracNum& pf2, int exps1, int exps2)
    : PrefactorSquared(pf2), Exponents({exps1, exps2}) {}

};

/// convert to string
std::string to_string(const PolynomialTerm& p);

/// check for equality to zero
const bool is_zero(const PolynomialTerm& p)
{ return p.PrefactorSquared == TFracNum::Zero; }

/// multiplication assignment
PolynomialTerm& operator*=(PolynomialTerm& p, const TFracNum& f)
{ p.PrefactorSquared *= f; return p; }

/// copy the PolynomialTerm and swap its exponents
PolynomialTerm swap_exponents(PolynomialTerm p)
{ std::swap(p.Exponents[0], p.Exponents[1]); return p; }

/// \typedef PolynomialTerms
using PolynomialTerms = std::vector<PolynomialTerm>;

/// convert to string
std::string to_string(const PolynomialTerms& P);

/// multiplation assignment
PolynomialTerms& operator*=(PolynomialTerms& P, const TFracNum& f)
{ std::for_each(P.begin(), P.end(), [&f](PolynomialTerm& p){p *= f;}); return P; }

/// copy PolynomialTerms and swap exponents
PolynomialTerms swap_exponents(PolynomialTerms P)
{ std::for_each(P.begin(), P.end(), [](PolynomialTerm& p){p = swap_exponents(p);}); return P; }

/*!
 \class TLSContrib
 \brief Relativistic LS-coupling contributions
 \author Jan.Friedrich@ph.tum.de
 */
class TLSContrib : public JLS
{
public:

    TLSContrib(const TLSAmpl& A, int delta, const TFracNum& scfac);

    void Add(const TLSContrib& rhs, bool particle_exchange);

    const bool pureRelativistic() const
    { return PureRelativistic_; }

    const int delta() const
    { return Delta_; }
    
    const unsigned contractionNumber() const
    { return ContractionNumber_; }
    
    const TFracNum& spinCG() const
    { return SpinCG_; }
    
    const TFracNum& normFactor() const
    { return NormFactor_; }

    const PolynomialTerms& polynomialTerms() const
    { return PolynomialTerms_; }

    void Print() const;
    void PrintNR() const;
    void PrintNRG(TFracNum) const;

    friend TLSContrib exchange_particles(TLSContrib c)
    {
        for (auto& p : c.PolynomialTerms_)
            p = swap_exponents(p);
        return c;
    }


private:
    unsigned ContractionNumber_;
    int Delta_;
    TFracNum SpinCG_;

    // Square  of normalisation factor
    TFracNum  NormFactor_{TFracNum::Zero};

    PolynomialTerms PolynomialTerms_;

    bool PureRelativistic_{false};

    static bool Debug_;

};

/// convert to string
std::string to_string(const TLSContrib& C);

/// equality operator
inline const bool operator==(const TLSContrib& lhs, const TLSContrib& rhs)
{
    return static_cast<const JLS&>(lhs) ==  static_cast<const JLS&>(rhs)
        and lhs.contractionNumber() == rhs.contractionNumber();
}
    

#endif
