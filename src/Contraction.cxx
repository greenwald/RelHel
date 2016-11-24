#include "Contraction.h"

#include "ClebschGordan.h"

#include <algorithm>

#include <iostream>

namespace relhel {

//-------------------------
namespace contractions {
const std::vector<Contraction::Term::index> psi1_psi2 = {Contraction::Term::index::psi1, Contraction::Term::index::psi2};
const std::vector<Contraction::Term::index> psi1_chi  = {Contraction::Term::index::psi1, Contraction::Term::index::chi};
const std::vector<Contraction::Term::index> psi2_chi  = {Contraction::Term::index::psi2, Contraction::Term::index::chi};
const std::vector<Contraction::Term::index> psi1_phi  = {Contraction::Term::index::psi1, Contraction::Term::index::phi};
const std::vector<Contraction::Term::index> psi2_phi  = {Contraction::Term::index::psi2, Contraction::Term::index::phi};
const std::vector<Contraction::Term::index> chi_phi   = {Contraction::Term::index::chi,  Contraction::Term::index::phi};
}

//-------------------------
Contraction::Contraction(const CoupledWaveFunctions& psi, const OrbitalAngularMomentumWaveFunction& chi, const WaveFunction& phi)
{
    // if chi has no projection = zero wave functions
    if (chi.projections().at(0).empty())
        return;
    
    // loop over spin projections of first daughter
    for (const auto& m_wps : psi.phi()[0].projections()) {
        if (m_wps.second.empty())
            continue;

        // find appropriate spin projection of second daughter
        auto it = psi.phi()[1].projections().find(psi.delta() - m_wps.first);
        if (it == psi.phi()[1].projections().end() or it->second.empty())
            continue;

        // loop over WaveProduct's of first daughter
        for (const auto& wp1 : m_wps.second) {
            auto two_j1 = spin(wp1);
            auto wp1_coeff = squared_coefficient(wp1);
            
            // loop over WaveProduct's of second daughter
            for (const auto& wp2 : it->second) {
                auto two_j2 = spin(wp2);
                auto psi_coeff = wp1_coeff * squared_coefficient(wp2)
                    * ClebschGordan::squared_coefficient(two_j1, m_wps.first, two_j2, it->first, psi.twoS());
                
                // loop over WaveProduct's of orbital angular momentum
                for (const auto& wp_chi : chi.projections().at(0)) {
                    auto l = spin(wp_chi) / 2;

                    auto psi_chi_coeff = psi_coeff * squared_coefficient(wp_chi) * pow(decompose(2), l);
                                
                    // loop over projections of initial state
                    for (const auto& m_wps_phi : phi.projections())
                        // loop over WaveProduct's of initial state
                        for (const auto& wp_phi : m_wps_phi.second)
                    
                            Terms_.emplace_back(wp1, wp2, wp_chi, wp_phi, squared_coefficient(wp_phi) * psi_chi_coeff);
                                 
                }
            }
        }
    }
    prune();
}

//-------------------------
const std::vector<int> gamma_exponents(const Contraction::Term& t, const std::vector<Contraction::Term::index>& I)
{
     std::vector<int> g;
     g.reserve(I.size());
     std::transform(I.begin(), I.end(), std::back_inserter(g), [&t](Contraction::Term::index i){return t.GammaExponents[static_cast<size_t>(i)]; });
     return g;
}

//-------------------------
// helper function
GammaPolynomial& add_term(GammaPolynomial& P, const GammaPolynomial::key_type& ges, const GammaPolynomial::mapped_type::value_type& coeff)
{
    // look for term with same exponents
    auto it = P.find(ges);
    // if found
    if (it != P.end()) {
        // look for coefficient term with same argument under the square root
        auto it2 = std::find_if(it->second.begin(), it->second.end(),
                                [&coeff](const GammaPolynomial::mapped_type::value_type& r)
                                {return r.second == coeff.second;});
        // if one is found, add arguments outside of square root sign
        if (it2 != it->second.end()) {
            it2->first += coeff.first;
            if (it2->first == 0)
                it->second.erase(it2);
        }
        // else add new term completely
        else it->second.push_back(coeff);
    }
    // else add new term
    else {
        P[ges].push_back(coeff);
    }
    return P;
}

//-------------------------
// helper function
GammaPolynomial& add_term(GammaPolynomial& P, const GammaPolynomial::key_type& ges, const GammaPolynomial::mapped_type& coeffs)
{
    std::for_each(coeffs.begin(), coeffs.end(), [&](const GammaPolynomial::mapped_type::value_type& c){add_term(P, ges, c);});
    return P;
}

//-------------------------
// helper function
GammaPolynomial& add_term(GammaPolynomial& P, const GammaPolynomial::value_type& ges_coeffs)
{ return add_term(P, ges_coeffs.first, ges_coeffs.second); }

//-------------------------
const GammaPolynomial gamma_polynomial(const std::vector<Contraction::Term>& terms)
{
    if (std::any_of(terms.begin(), terms.end(), [](const Contraction::Term& t){return rank(t) > 0;}))
        throw std::invalid_argument("Terms are not scalar; gamma_polynomial(...)");
    
    GammaPolynomial P;

    for (const auto& t : terms) {
        auto c = factorize_sqrt(t.SquaredCoefficient);
        add_term(P, gamma_exponents(t, contractions::psi1_psi2), std::make_pair(static_cast<double>(sqrt(c[0])), c[1]));
    }
    
    return P;
}

//-------------------------
// helper function
std::string to_exp_string(const GammaPolynomial::key_type& ges)
{
    std::string S;
    for (size_t i = 0; i < ges.size(); ++i)
        S += exponential_string(" * g_" + std::to_string(i + 1), ges[i]);
    return (S.empty()) ? S : S.erase(0, 3);
}

//-------------------------
// helper function
std::string to_coef_string(const GammaPolynomial::mapped_type& coefs)
{
    return "("
        + std::accumulate(coefs.begin(), coefs.end(), std::string(""),
                          [](std::string& S, const GammaPolynomial::mapped_type::value_type& r)
                          {return S += " + " + to_sqrt_string(RationalNumber(r.first * std::abs(r.first)) * r.second);}).erase(0, 3)
        + ")";
}

//-------------------------
GammaPolynomial& operator+=(GammaPolynomial& lhs, const GammaPolynomial& rhs)
{
    for (const auto& term : rhs)
        add_term(lhs, term);
    return lhs;
}

//-------------------------
std::string to_string(const GammaPolynomial::value_type& ges_coefs)
{
    auto s_ge = to_exp_string(ges_coefs.first);
    return s_ge + (!s_ge.empty() ? " * " : "") + to_coef_string(ges_coefs.second);
}

//-------------------------
std::string to_string(const GammaPolynomial& gp)
{
    
    return std::accumulate(gp.begin(), gp.end(), std::string(),
                           [](std::string& s, const GammaPolynomial::value_type& t){return s += " + " + to_string(t);}).erase(0, 3);
}

//-------------------------
std::string to_string(const Contraction::Term& term)
{
    return std::accumulate(term.WaveProducts.begin(), term.WaveProducts.end(), std::string(),
                           [](std::string& s, const WaveProduct& w){return s += ", " + to_string(w);}).erase(0, 2)
        + " : Gamma(" + std::accumulate(term.GammaExponents.begin(), term.GammaExponents.end(), std::string(),
                                        [](std::string& s, int g){return s += ", " + std::to_string(g);}).erase(0, 2) + ")"
        + " * " + to_sqrt_string(term.SquaredCoefficient);
}

//-------------------------
std::string to_string(const std::vector<Contraction::Term>& terms)
{
    return std::accumulate(terms.begin(), terms.end(), std::string(),
                           [](std::string& s, const Contraction::Term& t){return s += "\n" + to_string(t);}).erase(0, 1);
}

//-------------------------
std::string to_string(Contraction::Term::index i)
{
    switch(i) {
    case Contraction::Term::index::psi1:
        return "psi1";
    case Contraction::Term::index::psi2:
        return "psi2";
    case Contraction::Term::index::chi:
        return "chi";
    case Contraction::Term::index::phi:
        return "phi";
    default:
        return std::to_string(static_cast<size_t>(i));
    }
}

//-------------------------
const std::vector<Contraction::Term::index> triple_contraction(const Contraction::Term& t)
{
    std::vector<Contraction::Term::index> I;
    if (rank(t) != 3 or std::count_if(t.WaveProducts.begin(), t.WaveProducts.end(), [](const WaveProduct& w){return rank(w) == 1;}) != 3)
        return I;
    I = {Contraction::Term::index::psi1, Contraction::Term::index::psi2, Contraction::Term::index::chi, Contraction::Term::index::phi};
    I.erase(I.begin() + std::distance(t.WaveProducts.begin(), std::find_if(t.WaveProducts.begin(), t.WaveProducts.end(), [](const WaveProduct& w){return rank(w) == 0;})));
    return I;
}

//-------------------------
std::string to_string(const std::vector<Contraction::Term::index>& I)
{ return std::accumulate(I.begin(), I.end(), std::string(""), [](std::string& s, Contraction::Term::index i){return s += ", " + to_string(i);}).erase(0, 2); }

//-------------------------
void Contraction::Term::contract(std::vector<index> I)
{
    // check size
    if (I.size() < 2 or I.size() > 3)
        throw std::invalid_argument("Invalid number of contraction indices; contract(" + to_string(I) + ")");

    // check no WaveFunction is empty
    if (std::any_of(I.begin(), I.end(), [&](index i){return WaveProducts[static_cast<size_t>(i)].empty();})) {
        std::vector<index> J;
        J.reserve(I.size());
        std::copy_if(I.begin(), I.end(), std::back_inserter(J), [&](index i){return WaveProducts[static_cast<size_t>(i)].empty();});
        throw std::invalid_argument("Rank-zero Wave function(s) " + to_string(J) + "; contract(" + to_string(I) + ")");
    }
    
    // sort I (aids checking in the 3-index situation below)
    std::sort(I.begin(), I.end(), [](index a, index b){return std::less<size_t>()(static_cast<size_t>(a), static_cast<size_t>(b));});
    
    // check no index is repeated
    if (std::adjacent_find(I.begin(), I.end()) != I.end())
        throw std::invalid_argument("Can't contract WaveFunction with itself; contract(" + to_string(I) + ")");

    // get back elements and pop them off
    std::vector<int> m;
    m.reserve(I.size());
    for (auto i : I) {
        m.push_back(WaveProducts[static_cast<size_t>(i)].back());
        WaveProducts[static_cast<size_t>(i)].pop_back();
    }

    // two WaveFunction's contracted
    if (I.size() == 2) {

        // if contracting with phi
        if (I[1] == index::phi) {
            // if not the same 
            if (m[0] != m[1])
                SquaredCoefficient = RationalNumber(0);
        }
        else {
            // if opposite signs
            if (m[0] * m[1] < 0)
                SquaredCoefficient.negate();
            
            // if contracting with chi, and not opposite signed (or both zero)
            if (I[1] == index::chi and m[0] != -m[1])
                SquaredCoefficient = RationalNumber(0);
        }
    }
    // three WaveFunction's contracted
    else {

        // if not one and only one spin projection is zero
        if (std::count(m.begin(), m.end(), 0) != 1)
            SquaredCoefficient = RationalNumber(0);
        // or the first two projections are equal
        else if (m[0] == m[1])
            SquaredCoefficient = RationalNumber(0);
        else if (m[2] != 0) {
            // if contracting with initial state
            if (I[2] == index::phi) {
                if (m[0] + m[1] + m[2] == 0)
                    SquaredCoefficient = RationalNumber(0);
            } else {
                if (m[2] != 0 and m[0] + m[1] + m[2] != 0)
                    SquaredCoefficient = RationalNumber(0);
            }
        }

        // if not yet zero'ed, check negation situations
        if (!is_zero(SquaredCoefficient)) {

            // if m[2] == 0, m[0] dictates sign
            if (m[2] == 0 and m[0] < 0)
                SquaredCoefficient.negate();
            
            // if m[1] == 0, m[2] dictates sign
            if (m[1] == 0 and m[2] < 0)
                SquaredCoefficient.negate();
            
            // if m[0] == 0, m[2] dictates sign
            if (m[0] == 0 and m[2] > 0)
                SquaredCoefficient.negate();
        }
    }

    // increase relevant gamma factors
    for (size_t i = 0; i < I.size(); ++i)
        if (m[i] == 0) ++GammaExponents[static_cast<size_t>(I[i])];
}

}
