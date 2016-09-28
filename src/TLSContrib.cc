#include "TLSContrib.h"

#include "MathUtils.h"
#include "TLSAmpl.h"

#include <iostream>
#include <stdexcept>

using namespace std;

bool TLSContrib::Debug_ = false;

//-------------------------
TLSContrib::TLSContrib(const TLSAmpl& A, int delta, const TFracNum& scfac)
    : JLS(A),
      ContractionNumber_(A.contractionNumber()),
      Delta_(delta),
      SpinCG_(scfac)
{
    if (Debug_)
        cout << "(" << J() << ")" << L() << S() << "[" << Delta_ << "]" << endl;

    bool signReversal = Delta_ < 0 and is_odd(L() + S() + J());

    PolynomialTerms_.reserve(A.tensorSum().terms().size());
    for (const auto& t : A.tensorSum().terms()) {

        if (Debug_)
            cout << scfac.FracStringSqrt() << ","
                 << t.GetPreFac().FracStringSqrt()
                 << " L" << L() << "S" << S() << "J" << J()
                 << "Ampldelta" << A.delta()
                 << " delta" << Delta_ << ", sigrev" << signReversal;

        PolynomialTerms_.emplace_back(scfac * t.GetPreFac(), t.GetGamS(), t.GetGamSig());
        if (signReversal)
            PolynomialTerms_.back().PrefactorSquared.FlipSign();
        
        NormFactor_ = NormFactor_.SumSignedRoots(PolynomialTerms_.back().PrefactorSquared);
        if (Debug_)
            cout << " -> Normfactor: " << NormFactor_.FracStringSqrt() << endl;
    }

    PureRelativistic_ = NormFactor_ == TFracNum::Zero;
    
    if (PureRelativistic_)
        NormFactor_ = TFracNum::One;
    else
        PolynomialTerms_ *= invert(NormFactor_);
}

//-------------------------
void TLSContrib::Add(const TLSContrib& rhs, bool particleExchange)
{
    if (static_cast<const JLS&>(rhs) != *this)
        throw std::runtime_error("TLSContrib::Add : Trying to add different (J;L,S): "
                                 + to_string(static_cast<const JLS&>(rhs)) + " != "
                                 + to_string(static_cast<const JLS&>(*this)));

    //
    // Include normalisation factor and the factor (1/2)**2 in the squared
    // representation of prefactors
    //

    PolynomialTerms_ *= TFracNum::Quarter * NormFactor_;

    auto P = (particleExchange) ? swap_exponents(rhs.polynomialTerms()) : rhs.polynomialTerms();
    
    if (rhs.contractionNumber() == contractionNumber())
        for (auto it = P.begin(); it != P.end(); ) {

            // find term with matching exponents
            auto it2 = std::find(PolynomialTerms_.begin(), PolynomialTerms_.end(),
                                 [&it](const PolynomialTerm& p){return p.Exponents == it->Exponents;});

            // if match found, add term in
            if (it2 != PolynomialTerms_.end()) {
                auto bterm = TFracNum::Quarter * rhs.normFactor() * it->PrefactorSquared;
                if (is_even(J()))
                    bterm.FlipSign();
                it2->PrefactorSquared = bterm.SumSignedRoots(it2->PrefactorSquared);
                it = P.erase(it);
            } else // skip to next term
                ++it;
        }

    // loop over remaining terms not yet added in, and add them as new terms
    P *= TFracNum::Quarter * rhs.normFactor();
    for (const auto& p : P) {
        if (is_even(J()))
            if (!p.PrefactorSquared.FlipSign())
                throw std::runtime_error("Flipping sign failed.");
        PolynomialTerms_.push_back(p);
    }

    // Eliminate zero entries
    PolynomialTerms_.erase(std::find_if(PolynomialTerms_.begin(), PolynomialTerms_.end(), is_zero), PolynomialTerms_.end());
    
    // Recalculate Normalization Factor
    NormFactor_ = std::accumulate(PolynomialTerms_.begin(), PolynomialTerms_.end(), TFracNum::Zero,
                                  [](TFracNum& N, const PolynomialTerms& p)
                                  {return N = N.SumSignedRoots(p.PrefactorSquared);});

    PureRelativistic_ = NormFactor_ == TFracNum::Zero;
    
    // Apply normalization
    if (PureRelativistic_)
        NormFactor_ = TFracNum::One;
    else
        PolynomialTerms_ *= invert(NormFactor_);
}

//-------------------------
std::string exponent_string(const std::string& s, int n)
{ return (n == 0) ? "" : s + (n == 1) ? "" : "^" + std::to_string(n); }

//-------------------------
std::string to_string(const PolynomialTerms& p)
{
    return p.PrefactorSquared.FracStringSqrt()
        + exponent_string(" g_0", p.Exponents[0])
        + exponent_string(" g_1", p.Exponents[1]);
}
        
//-------------------------
std::string to_string(const PolynomialTerms& P)
{ return P.empty() ? "" : std::accumulate(P.begin(), P.end(), std::string(), to_string); }

//-------------------------
std::string to_string(const TLSContrib& C)
{
    std::string s = "";
    switch (C.contractionNumber()) {
    case 1:
        s = "g";
        break;
    case 2:
        s = "f";
        break;
    case 3:
        s = "h";
        break;
    default:
        throw;
    }
    s += "[" + ContractionNumber_ + "]" + to_string(static_cast<const JLS&>(C));
    if (!PureRelativistic_)
        s += " (" << C.normFactor().FracStringSqrt() + ")";
    s += to_string(PolynomialTerms_);
    if (PolynomialTerms_.empty())
        s += "(0)";
    return s;
}


void TLSContrib::PrintNR() const
{
    cout << _NormFactor.FracStringSqrt();
    if (_cNum == 1) {
        cout << "*g";
    }
    if (_cNum == 2) {
        cout << "*f";
    }
    if (_cNum == 3) {
        cout << "*h";
    }
}


void TLSContrib::PrintNRG(TFracNum m) const
{
    cout << (_NormFactor * m).FracStringSqrt();
    if (_cNum == 1) {
        cout << "*g";
    }
    if (_cNum == 2) {
        cout << "*f";
    }
    if (_cNum == 3) {
        cout << "*h";
    }
}
