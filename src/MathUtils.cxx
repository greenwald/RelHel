#include "MathUtils.h"

#include <algorithm>
#include <cstdlib>
#include <functional>

namespace relhel {

//-------------------------
const std::vector<unsigned> triangle(unsigned two_s1, unsigned two_s2)
{
    unsigned two_smin = std::abs(static_cast<int>(two_s1) - static_cast<int>(two_s2));
    unsigned two_smax = two_s1 + two_s2;
    std::vector<unsigned> two_S;
    two_S.reserve((two_smax - two_smin) / 2 + 1);
    for (auto two_s = two_smin; two_s <= two_smax; two_s += 2)
        two_S.push_back(two_s);
    return two_S;
}

//-------------------------
const std::vector<int> projections(unsigned two_j)
{
    std::vector<int> two_M;
    for (int two_m = -two_j; two_m <= static_cast<int>(two_j); ++two_m)
        two_M.push_back(two_m);
    return two_M;
}
    
//-------------------------
const std::vector<std::vector<int> > projections(const std::vector<unsigned>& two_J)
{
    // initialize vector of spin projections to -two_j
    std::vector<int> two_M;
    two_M.reserve(two_J.size());
    std::transform(two_J.begin(), two_J.end(), std::back_inserter(two_M), std::negate<int>());
    
    std::vector<std::vector<int> > SPV;
    // fill SPV with "odometer"-style looping
    while (two_M.back() <= (int)two_J.back()) {
        SPV.push_back(two_M);
        two_M[0] += 2;
        for (size_t i = 0; (i < two_M.size() - 1) and (two_M[i] > (int)two_J[i]); ++i) {
            two_M[i] = -two_J[i];
            two_M[i + 1] += 2;
        }
    }
    return SPV;
}

}
