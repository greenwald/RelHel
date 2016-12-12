#include "Contraction.h"

#include "ClebschGordan.h"

namespace relhel {

//-------------------------
GammaPolynomial::GammaPolynomial(const std::vector<GammaIndex>& I) :
    Indices_(I)
{
    // sort Indices
    std::sort(Indices_.begin(), Indices_.end());
    if (std::adjacent_find(Indices_.begin(), Indices_.end()) != Indices_.end())
        throw std::invalid_argument("indices not unique; GammaPolynomial::GammaPolynomial");
}

//-------------------------
GammaPolynomial::Term& operator*=(GammaPolynomial::Term& lhs, const GammaPolynomial::Term& rhs)
{
    std::vector<RationalNumber> C2;
    for (const auto& c2_l : lhs.SquaredCoefficients)
        for (const auto& c2_r : rhs.SquaredCoefficients)
            C2.push_back(c2_l * c2_r);
    lhs.SquaredCoefficients = C2;
    std::transform(lhs.Exponents.begin(), lhs.Exponents.end(), rhs.Exponents.begin(), lhs.Exponents.begin(), std::plus<unsigned>());
    return lhs;
}

//-------------------------
GammaPolynomial& GammaPolynomial::operator+=(const Term& term)
{
    // look for a term with the same exponents
    auto it = std::find_if(Terms_.begin(), Terms_.end(), std::bind(same_exponents, term, std::placeholders::_1, indices()));
    
    // if none found
    if (it == Terms_.end()) {
        // add new term
        Terms_.push_back(term);
        return *this;
    }

    // else add to existing term
    for (const auto& c2 : term.SquaredCoefficients) {
        
        auto sqrt_c2 = factorize_sqrt(c2);

        // look for squared coefficient with same unsqrt'able part
        for (auto& C2 : it->SquaredCoefficients) {

            auto sqrt_C2 = factorize_sqrt(C2);

            // if unsqrt'able part is the same, add sqrt's of sqrt'able part
            if (sqrt_C2[1] == sqrt_c2[1]) {
                C2 = sqrt_C2[1] * pow(sqrt(sqrt_C2[0]) + sqrt(sqrt_c2[0]) , 2);
                sqrt_c2[0] = RationalNumber(0);
            }
        }

        // if none found
        if (!is_zero(sqrt_c2[0]))
            it->SquaredCoefficients.push_back(c2);
    }

    // remove zero'ed coeffs
    it->SquaredCoefficients.erase(std::remove_if(it->SquaredCoefficients.begin(), it->SquaredCoefficients.end(), is_zero), it->SquaredCoefficients.end());

    // if term now empty, remove it
    if (it->SquaredCoefficients.empty())
        Terms_.erase(it);

    return *this;
}

//-------------------------
std::vector<EpsilonContraction> epsilon_contractions(const std::vector<unsigned>& r)
{
    if (r.size() > GammaIndex::n_indices)
        throw std::invalid_argument("rank vector too large; epsilon_contractions");

    // if even total rank
    if (is_even(std::accumulate(r.begin(), r.end(), 0u)))
        return {{}};

    // else
    std::vector<EpsilonContraction> B;

    for (unsigned i = 0; i <= r.size() - 3; ++i)
        if (r[i] > 0)
            for (unsigned j = i + 1; j <= r.size() - 2; ++j)
                if (r[j] > 0)
                    for (unsigned k = j + 1; k <= r.size() - 1; ++k)
                        if (r[k] > 0)
                            B.push_back({static_cast<GammaIndex>(i), static_cast<GammaIndex>(j), static_cast<GammaIndex>(k)});
    return B;
}

//-------------------------
// return vector of contraction matrices
std::set<ContractionMatrix> contraction_matrices(std::vector<unsigned> r, const EpsilonContraction& b)
{
    // modify ranks by epsilon contraction
    for (auto i : b)
        --r[i];

    std::set<ContractionMatrix> S;
    ContractionMatrix C(r.size() - 1, ContractionMatrix::value_type(r.size(), 0));

    // loop through contractions odometer style
    while (C.back().back() <= std::min(r[r.size() - 2], r.back())) {
        
        // check if contraction is full (over all ranks)
        std::vector<unsigned> c(r.size(), 0);
        for (size_t i = 0; i < C.size(); ++i)
            for (size_t j = i + 1; j < C[i].size(); ++j) {
                c[i] += C[i][j];
                c[j] += C[i][j];
            }
        if (c == r)
            S.insert(C);

        // increment odometer
        ++C[0][1];
        
        size_t i = 0;
        size_t j = 1;
        while(C[i][j] > std::min(r[i], r[j])) {
            C[i][j] = 0;
            
            ++j;
            
            if (j >= C[i].size()) {
                ++i;
                j = i + 1;
            }
            
            ++C[i][j];
            
            if (i == C.size() - 1 and j == C.size())
                break;
        }
    }

    return S;
}

//-------------------------
std::map<EpsilonContraction, std::set<ContractionMatrix> > contractions(const CoupledWaveFunctions& psi,
                                                                        const OrbitalAngularMomentumWaveFunction& chi,
                                                                        const WaveFunction& phi)
{
    std::map<EpsilonContraction, std::set<ContractionMatrix> > M;

    std::vector<unsigned> r = {rank(psi, 0), rank(psi, 1), rank(chi), rank(phi)};
    
    // loop over epsilon contractions
    for (const auto& b : epsilon_contractions(r))
        // loop through possible contractions
        for (const auto& C : contraction_matrices(r, b))
            M[b].insert(C);

    return M;
}

//-------------------------
std::string to_string(const std::map<EpsilonContraction, std::set<ContractionMatrix> >& bC)
{
    std::string S;
    for (const auto& b_S : bC)
        for (const auto& C : b_S.second) {
            std::string s;
            for (size_t i = 0; i < C.size(); ++i)
                for (size_t j = 0; j < C[i].size(); ++j)
                    if (C[i][j] > 0)
                        s += "  +  " + std::to_string(C[i][j]) + " * [" + to_string(static_cast<GammaIndex>(i))
                            + " x " + to_string(static_cast<GammaIndex>(j)) + "]";
            if (!b_S.first.empty())
                s += "  +  [" + std::accumulate(b_S.first.begin(), b_S.first.end(), std::string(),
                                           [](std::string& s, GammaIndex bi)
                                           {return s += " x " + to_string(bi);}).erase(0, 3) + "]";
            S += s.erase(0, 5) + "\n";
        }
    return S;
}

//-------------------------
// set to zero and return through
GammaPolynomial::Term& zero(GammaPolynomial::Term& term)
{
    for (auto& c2 : term.SquaredCoefficients)
        c2 = RationalNumber(0);
    return term;
}

//-------------------------
// negate and return through
GammaPolynomial::Term& negate(GammaPolynomial::Term& term)
{
    for (auto& c2 : term.SquaredCoefficients)
        c2.negate();
    return term;
}

//-------------------------
// perform single contraction
const GammaPolynomial::Term contract(std::vector<WaveProduct>& WP, unsigned i, unsigned j)
{
    GammaPolynomial::Term term;
    
    /////////////////////////
    // get projections to contract, popping them off their WaveProduct's
    std::vector<int> m;
    m.push_back(WP[i].back()); WP[i].pop_back();
    m.push_back(WP[j].back()); WP[j].pop_back();
    
    // increase gamma factors
    if (m[0] == 0) ++term.Exponents[i];
    if (m[1] == 0) ++term.Exponents[j];
    
    // if contracting with phi
    if (j == GammaIndex::phi) {
        // if projections not the same
        if (m[0] != m[1])
            zero(term);
        return term;
    }

    // else
    // if opposite signs
    if (m[0] * m[1] < 0)
        negate(term);
        
    // if contracting with chi, and not opposite signed (or both zero)
    if (j == GammaIndex::chi and m[0] != -m[1])
        zero(term);

    return term;
}

//-------------------------
// perform single epsilon contraction
const GammaPolynomial::Term contract(std::vector<WaveProduct>& WP, const EpsilonContraction& b)
{
    GammaPolynomial::Term term;
    
    // get projections, and pop off elements
    std::vector<int> m;
    for (auto i : b) {
        m.push_back(WP[i].back());
        WP.pop_back();
    }
    
    // increase gamma factors
    for (size_t i = 0; i < m.size(); ++i)
        if (m[i] == 0) ++term.Exponents[b[i]];
    
    // if not one and only one spin projection is zero
    if (std::count(m.begin(), m.end(), 0) != 1)
        return zero(term);

    // or the first two projections are equal
    if (m[0] == m[1])
        return zero(term);

    if (m[2] != 0) {
        // if contracting with initial state
        if (b[2] == GammaIndex::phi) {
            if (m[0] + m[1] + m[2] == 0)
                return zero(term);
        } else {
            if (m[2] != 0 and m[0] + m[1] + m[2] != 0)
                return zero(term);
        }
    }
    
    // not yet zero'ed, check negation situations

    // if m[2] == 0, m[0] dictates sign
    if (m[2] == 0 and m[0] < 0)
        negate(term);
        
    // if m[1] == 0, m[2] dictates sign
    if (m[1] == 0 and m[2] < 0)
        negate(term);
        
    // if m[0] == 0, m[2] dictates sign
    if (m[0] == 0 and m[2] > 0)
        negate(term);

    return term;
}

//-------------------------
// perform contractions
const GammaPolynomial::Term contract(std::vector<WaveProduct>& WP, const ContractionMatrix& C, const EpsilonContraction& b)
{
    GammaPolynomial::Term term;

    for (size_t i = 0; i < C.size() and !is_zero(term.SquaredCoefficients[0]); ++i)
        for (size_t j = i + 1; j < C[i].size() and !is_zero(term.SquaredCoefficients[0]); ++j)
            for (size_t c = 0; c < C[i][j] and !is_zero(term.SquaredCoefficients[0]); ++c)
                term *= contract(WP, i, j);
    
    if (is_zero(term.SquaredCoefficients[0]))
        return term;

    return term *= contract(WP, b);
}

//-------------------------
// return C-G coeff^2
const RationalNumber squared_CG(const std::vector<WaveProduct>& WP, unsigned two_S)
{
    return ClebschGordan::squared_coefficient(spin(WP[GammaIndex::psi1]), projection(WP[GammaIndex::psi1]),
                                              spin(WP[GammaIndex::psi2]), projection(WP[GammaIndex::psi2]),
                                              two_S);
}

//-------------------------
// product of squared_coefficient for each WaveProdcut * CG-Coeff * 2^L
const RationalNumber squared_coefficient(const std::vector<WaveProduct>& WP, unsigned two_S)
{
    return std::accumulate(WP.begin(), WP.end(), squared_CG(WP, two_S) * pow(decompose(2), spin(WP[GammaIndex::chi]) / 2),
                           [](RationalNumber& c2, const WaveProduct& wp){return c2 *= squared_coefficient(wp);});
}

//-------------------------
const GammaPolynomial contract(const CoupledWaveFunctions& psi, const OrbitalAngularMomentumWaveFunction& chi, const WaveFunction& phi, const std::vector<GammaIndex>& I)
{
    /////////////////////////
    // Create vector of wave product combinations to contract
    std::vector<std::vector<WaveProduct> > WPV;
    // loop over CoupledWaveProducts of psi
    for (const auto& cwp_psi : psi.products())
        // loop over WaveProducts of chi
        for (const auto& wp_chi : chi.projections().at(0))
            // loop over WaveProducts of phi
            for (const auto& wp_phi : phi.projections().at(projection(cwp_psi[0]) - projection(cwp_psi[1])))
                WPV.push_back({cwp_psi[0], cwp_psi[1], wp_chi, wp_phi});

    // store ranks
    std::vector<unsigned> r = {rank(psi, 0), rank(psi, 1), rank(chi), rank(phi)};

    GammaPolynomial P(I);
    
    // loop over epsilon contractions
    for (const auto& b : epsilon_contractions(r)) {

        /////////////////////////
        // loop through possible contractions
        for (const auto& C : contraction_matrices(r, b)) {

            // loop over WPs
            for (auto WP : WPV) {

                // calculate squared coefficient
                RationalNumber C2 = squared_coefficient(WP, psi.spin());
                
                GammaPolynomial::Term term = contract(WP, C, b);

                if (is_zero(term.SquaredCoefficients[0]))
                    continue;

                if (std::any_of(WP.begin(), WP.end(), [](const WaveProduct& wp){return !wp.empty();}))
                    throw std::runtime_error("wave functions not fully contracted; contract");

                term.SquaredCoefficients[0] *= C2;

                P += term;
            }
        }
    }
    return P;
}

//-------------------------
std::string to_index_string(GammaIndex i)
{
    switch(i) {
    case GammaIndex::psi1:
        return "1";
    case GammaIndex::psi2:
        return "2";
    case GammaIndex::chi:
        return "l";
    case GammaIndex::phi:
        return "0";
    default:
        return std::to_string(static_cast<size_t>(i));
    }
}

//-------------------------
std::string to_string(GammaIndex i)
{
    switch(i) {
    case GammaIndex::psi1:
        return "psi1";
    case GammaIndex::psi2:
        return "psi2";
    case GammaIndex::chi:
        return "chi";
    case GammaIndex::phi:
        return "phi";
    default:
        return std::to_string(static_cast<size_t>(i));
    }
}

//-------------------------
std::string to_string(const GammaPolynomial::Term& t, const std::vector<GammaIndex>& I)
{
    if (t.SquaredCoefficients.empty())
        return "(zero)";
    
    // get exponent string
    std::string exp_s;
    for (auto i : I)
        exp_s += exponential_string(" * g_" + to_index_string(i), t.Exponents[i]);
    if (!exp_s.empty())
        exp_s.erase(0, 3);

    // get coefficient string
    std::string coef_s;
    for (const auto& c2 : t.SquaredCoefficients)
        coef_s += " + " + to_sqrt_string(c2);
    coef_s.erase(0, 3);
    if (t.SquaredCoefficients.size() > 1)
        coef_s = "(" + coef_s + ")";

    if (exp_s.empty())
        return coef_s;

    if (t.SquaredCoefficients.size() == 1 and is_one(t.SquaredCoefficients[0]))
        return exp_s;
        
    return exp_s + " * " + coef_s;
}

//-------------------------
std::string to_string(const GammaPolynomial& P)
{
    return std::accumulate(P.terms().begin(), P.terms().end(), std::string(),
                           [&](std::string& s, const GammaPolynomial::Term& t)
                           {return s += " + " + to_string(t, P.indices());}).erase(0, 3);
}

}
