#include "WaveFunction.h"

#include "ClebschGordan.h"

#include <stdexcept>

namespace relhel {

//-------------------------
const RationalNumber squared_coefficient(const WaveProduct& v)
{
    unsigned j = spin(v) / 2;
    int two_m = projection(v);
    if (is_odd(two_m))
        throw std::runtime_error("projection is not integer; squared_coefficient(...)");
    int m = two_m / 2;
    return decompose_factorial(j + m) * decompose_factorial(j - m) * pow(decompose(2), zeroes(v)) / decompose_factorial(2 * j);
}

//-------------------------
std::string to_string(const WaveProduct& wp)
{
    return std::accumulate(wp.begin(), wp.end(), std::string("("), [](std::string& s, int m){return s += (m > 0) ? "+" : ((m < 0) ? "-" : "0");}) + ")";
}

//-------------------------
std::string to_string(const WaveProductSum& wps)
{
    return std::accumulate(wps.begin(), wps.end(), std::string(), [](std::string& s, const WaveProduct& w){return s += ", " + to_string(w);}).erase(0, 2);
}

//-------------------------
WaveFunction::WaveFunction(unsigned two_j)
{
    if (is_odd(two_j))
        throw std::invalid_argument("currently only supporting integer spins; WaveFunction::WaveFunction");

    WaveProduct WP(two_j / 2, -2);

    if (two_j == 0)
        Projections_[0] = WaveProductSum(1, WP);
    else {
         // loop through permutations using odometer method
        while (WP.back() <= 2) {
            // add product into sum for corresponding projection
            Projections_[projection(WP)].push_back(WP);
            // tick over the first projection
            WP[0] += 2;
            // check what projections must be reset or ticked over
            for (size_t i = 0; (i < WP.size() - 1) and (WP[i] > 2); ++i) {
                WP[i] = - 2;
                WP[i + 1] += 2;
            }
        }
    }
}

//-------------------------
const unsigned spin(const WaveFunction& wf)
{
    // look for non-empty spin projection in WaveFunction
    auto it = std::find_if(wf.projections().begin(), wf.projections().end(), [](const WaveFunction::map_type::value_type& m_wps){return !m_wps.second.empty();});
    // if none, return 0
    if (it == wf.projections().end())
        return 0;
    // return spin of WaveProductSum
    return spin(it->second);
}

//-------------------------
const unsigned rank(const WaveFunction& wf)
{
    // look for non-empty spin projection in WaveFunction
    auto it = std::find_if(wf.projections().begin(), wf.projections().end(), [](const WaveFunction::map_type::value_type& m_wps){return !m_wps.second.empty();});
    // if none, return 0
    if (it == wf.projections().end())
        return 0;
    // return tank of WaveProductSum
    return rank(it->second);
}

//-------------------------
// helper function
std::string spin_projection_string(unsigned two_j, int two_m)
{ return std::string("|") + spin_to_string(two_j) + ", " + spin_to_string(two_m) + ">"; }

//-------------------------
std::string to_string(const WaveFunction& wf)
{
    auto two_j = spin(wf);
    return std::accumulate(wf.projections().begin(), wf.projections().end(), std::string(),
                           [&](std::string& s, const WaveFunction::map_type::value_type& m_wps)
                           {return s += "\n" + spin_projection_string(two_j, m_wps.first) + " = " + to_string(m_wps.second);}).erase(0, 1);
}

//-------------------------
CoupledWaveFunctions::CoupledWaveFunctions(const WaveFunction& A, const WaveFunction& B, unsigned two_s, int two_m) :
    TwoS_(two_s)
{
    auto two_jA = ::relhel::spin(A);
    auto two_jB = ::relhel::spin(B);

    if (!triangle(two_jA, two_jB, TwoS_))
        throw std::invalid_argument("triangle(A, B, S) unfulfilled; CoupledWaveFunctions");
    
    if (is_odd(TwoS_ + two_m) or std::abs(two_m) > TwoS_)
        throw std::invalid_argument("invalid spin projection; CoupledWaveFunctions");

    // loop over A's projections
    for (const auto& m_wps : A.projections()) {
        // if empty, continue
        if (m_wps.second.empty())
            continue;
        
        // find appropriate spin projections 
        auto it = B.projections().find(two_m - m_wps.first);
        // if none, continue
        if (it == B.projections().end() or it->second.empty())
            continue;

        // loop over wave products
        for (const auto& a : m_wps.second)
            for (const auto& b : it->second)
                Products_.push_back({a, b});
    }
}

}
