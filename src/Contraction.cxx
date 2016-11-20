#include "Contraction.h"

#include "ClebschGordan.h"

#include <algorithm>

namespace relhel {

//-------------------------
Contraction::Contraction(const CoupledWaveFunctions& psi, const WaveFunction& chi)
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

            // loop over WaveProduct's of second daughter
            for (const auto& wp2 : it->second) {
                auto two_j2 = spin(wp2);

                // loop over WaveProduct's of orbital angular momentum
                for (const auto& wp_chi : chi.projections().at(0)) {
                    auto l = spin(wp_chi) / 2;

                    auto coeff = ClebschGordan::squared_coefficient(two_j1, m_wps.first, two_j2, it->first, psi.twoS())
                        * squared_coefficient(wp1)
                        * squared_coefficient(wp2)
                        * squared_coefficient(wp_chi) * pow(decompose(2), l);
                    
                    Terms_.emplace_back(wp1, wp2, wp_chi, coeff);
                                 
                }
            }
        }
    }
    prune();
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
void Contraction::Term::contract(const std::array<index, 2>& I)
{
    if (std::any_of(I.begin(), I.end(), [&](index i){return WaveProducts[i].empty();}))
        throw std::invalid_argument("Wave function is rank 0; contract_daughters(...)");

    // get back elements and pop them off
    std::vector<int> m;
    m.reserve(2);
    for (auto i : I) {
        m.push_back(WaveProducts[i].back());
        WaveProducts[i].pop_back();
    }

    // if elts are opposite sign, flip coefficient sign
    if (m[0] * m[1] < 0)
        SquaredCoefficient.negate();
    // else if both are zero, increase gamma exponents
    else if (m[0] == 0 and m[1] == 0) {
        for (auto i : I)
            ++GammaExponents[i];
    }
    // else set coefficient zero
    else SquaredCoefficient = RationalNumber(0);
}

}
