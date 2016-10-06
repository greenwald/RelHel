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
f_map f(unsigned j, int lambda)
{
    if (abs(lambda) > j)
        throw std::runtime_error("abs(lambda) > j");

    RationalNumber a(decompose_factorial(j + lambda) * decompose_factorial(j - lambda),
                     decompose_factorial(2 * j));
    
    f_map F;
    for (unsigned n = is_even(j + lambda) ? 0 : 1; n <= j - abs(lambda); n += 2)  {
        F.emplace(n, a * RationalNumber(decompose_factorial(j) * pow(decompose(2), n),
                                        decompose_factorial(n)
                                        * decompose_factorial((j + lambda - n) / 2)
                                        * decompose_factorial((j - lambda - n) / 2)));
    }
    return F;
}

//-------------------------
using ff_map = std::map<std::array<unsigned, 2>, RationalNumber>;

//-------------------------
ff_map ff(unsigned j1, int lambda1, unsigned j2, int lambda2)
{
    auto F1 = f(j1, lambda1);
    auto F2 = f(j2, lambda2);
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

    unsigned J = 0;
    std::vector<unsigned> j = {1 , 1};
    
    for (auto S : triangle(j[0], j[1]))
        for (auto L : triangle(J, S)) {
            RationalNumber LJ2(2 * L + 1, 2 * J + 1);
            for (auto lambda : projections(j)) {

                auto CG2 = ClebschGordan::squared_coefficient(2 * L, 0, 2 * S, 2 * (lambda[0] - lambda[1]), 2 * J)
                    * ClebschGordan::squared_coefficient(2 * j[0], 2 * lambda[0], 2 * j[1], -2 * lambda[1], 2 * S);

                if (is_zero(CG2))
                    continue;
                
                std::cout << J << " -> "
                          << j[0] << " [" << prefix_string(lambda[0]) << lambda[0] << "]"
                          << " + "
                          << j[1] << " [" << prefix_string(lambda[1]) << lambda[1] << "]"
                          << ", L = " << L
                          << ", S = " << S
                          << "\t" << std::flush;

                if (is_zero(CG2)) {
                    std::cout << "0" << std::endl;
                    continue;
                }

                auto term = ff(j[0], lambda[0], j[1], lambda[1]);

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
