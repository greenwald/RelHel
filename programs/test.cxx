#include "ClebschGordan.h"
#include "MathUtils.h"
#include "RationalNumber.h"

#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using namespace relhel;

//-------------------------
using f_map = std::map<unsigned, RationalNumber>;

//-------------------------
f_map f(unsigned two_j, int two_lambda)
{
    if (is_odd(two_j + two_lambda))
        throw std::invalid_argument("j and lambda must be either both integral or half-integral");
    
    unsigned j_p_lambda = (two_j + two_lambda) / 2;
    unsigned j_m_lambda = (two_j - two_lambda) / 2;
    
    if (abs(two_lambda) > two_j)
        throw std::runtime_error("abs(lambda) > j");

    RationalNumber a(decompose_factorial(j_p_lambda) * decompose_factorial(j_m_lambda),
                     decompose_factorial(two_j));
    
    f_map F;
    for (unsigned n = is_even(j_p_lambda) ? 0 : 1; n <= (two_j - abs(two_lambda)) / 2; n += 2)  {
        F.emplace(n, a * RationalNumber(decompose_factorial(two_j / 2) * pow(decompose(2), n),
                                        decompose_factorial(n)
                                        * decompose_factorial((j_p_lambda - n) / 2)
                                        * decompose_factorial((j_m_lambda - n) / 2)));
    }
    return F;
}

//-------------------------
using ff_map = std::map<std::array<unsigned, 2>, RationalNumber>;

//-------------------------
ff_map ff(unsigned two_j1, int two_lambda1, unsigned two_j2, int two_lambda2)
{
    auto F1 = f(two_j1, two_lambda1);
    auto F2 = f(two_j2, two_lambda2);
    ff_map FF;
    for (const auto& f1 : F1)
        for (const auto& f2 : F2)
            FF.emplace(ff_map::key_type({f1.first, f2.first}), f1.second * f2.second);
    return FF;
}

//-------------------------
const bool is_one(const ff_map& FF)
{ return FF.empty() or (FF.size() == 1 and FF.begin()->first[0] == 0 and FF.begin()->first[1] == 0); }

//-------------------------
std::string to_string(const ff_map::value_type& ff)
{
    return (!is_one(ff.second) or (ff.first[0] == 0 and ff.first[1] == 0) ? to_string(ff.second) : "")
        + (!is_one(ff.second) and (ff.first[0] != 0 or ff.first[1] != 0) ? " * " : "")
        + exponential_string("g_1", ff.first[0])
        + (ff.first[0] != 0 and ff.first[1] != 0 ? " * " : "")
        + exponential_string("g_2", ff.first[1]);
}

//-------------------------
std::string to_string(const ff_map& ff)
{
    return std::accumulate(ff.begin(), ff.end(), std::string(),
                           [](std::string& s, const ff_map::value_type& f)
                           { return s += " + " + to_string(f); }).erase(0, 3);
}

//-------------------------
std::string prefix_string(int j)
{ if (j > 0) return "+"; return (j == 0) ? " " : ""; }



//-------------------------
int main()
{

    unsigned two_J = 0;
    std::vector<unsigned> two_j = {2 , 2};
    
    for (auto two_S : triangle(two_j[0], two_j[1]))
        for (auto two_L : triangle(two_J, two_S)) {
            RationalNumber LJ2(two_L + 1, two_J + 1);
            for (auto two_lambda : projections(two_j)) {

                auto CG2 = ClebschGordan::squared_coefficient(two_L, 0, two_S, two_lambda[0] - two_lambda[1], two_J)
                    * ClebschGordan::squared_coefficient(two_j[0], two_lambda[0], two_j[1], -two_lambda[1], two_S);

                if (is_zero(CG2))
                    continue;
                
                std::cout << spin_to_string(two_J) << " -> "
                          << spin_to_string(two_j[0]) << " [" << prefix_string(two_lambda[0]) << spin_to_string(two_lambda[0]) << "]"
                          << " + "
                          << spin_to_string(two_j[1]) << " [" << prefix_string(two_lambda[1]) << spin_to_string(two_lambda[1]) << "]"
                          << ", L = " << spin_to_string(two_L)
                          << ", S = " << spin_to_string(two_S)
                          << "\t" << std::flush;

                if (is_zero(CG2)) {
                    std::cout << "0" << std::endl;
                    continue;
                }

                auto term = ff(two_j[0], two_lambda[0], two_j[1], two_lambda[1]);

                if (is_one(CG2) and is_one(LJ2)) {
                    std::cout << to_string(term) << std::endl;
                    continue;
                }

                std::cout << to_sqrt_string(LJ2 * CG2) << std::flush;

                if (!is_one(term))
                    std::cout << " * "
                              << (term.size() > 1 ? "(" : "")
                              << to_string(term)
                              << (term.size() > 1 ? ")" : "") << std::flush;

                std::cout << std::endl;
            }
        }
    
    return 0;
}
