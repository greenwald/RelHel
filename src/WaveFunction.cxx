#include "WaveFunction.h"

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
// helper function
std::string spin_projection_string(unsigned two_j, int two_m)
{ return std::string("|") + spin_to_string(two_j) + ", " + spin_to_string(two_m) + ">"; }

//-------------------------
std::string to_string(const WaveFunction& wf)
{
    auto two_j = spin(wf);
    if (two_j == 0)
        return spin_projection_string(0, 0);
    return std::accumulate(wf.projections().begin(), wf.projections().end(), std::string(),
                           [&](std::string& s, const WaveFunction::map_type::value_type& m_wps)
                           {return s += "\n" + spin_projection_string(two_j, m_wps.first) + " = " + to_string(m_wps.second);}).erase(0, 1);
}

//-------------------------
CoupledWaveFunctions::CoupledWaveFunctions(const WaveFunction& phi1, const WaveFunction& phi2, unsigned two_s, int delta)
    : Phi_({phi1, phi2}), TwoS_(two_s), Delta_(delta)
{
    auto two_j0 = spin(Phi_[0]);
    auto two_j1 = spin(Phi_[1]);
    
    if (!triangle(two_j0, two_j1, TwoS_))
        throw std::invalid_argument("invalid spins for coupling: triangle("
                                    + spin_to_string(two_j0) + ", "
                                    + spin_to_string(two_j1) + ", "
                                    + spin_to_string(TwoS_) + ") = false; get_spin_coupled_tensor_sum(...)");
    
    if (is_odd(two_j0 + two_j1 + Delta_) or std::abs(Delta_) > (two_j0 + two_j1))
        throw std::invalid_argument("invalid spin projection for coupled state; get_spin_coupled_tensor_sum(...)");
}

}
