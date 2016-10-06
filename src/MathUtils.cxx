#include "MathUtils.h"

#include <algorithm>
#include <cstdlib>

namespace relhel {

//-------------------------
const std::vector<unsigned> triangle(unsigned s1, unsigned s2)
{
    unsigned smin = std::abs(static_cast<int>(s1) - static_cast<int>(s2));
    unsigned smax = s1 + s2;
    std::vector<unsigned> S;
    S.reserve(smax - smin + 1);
    for (auto s = smin; s <= smax; ++s)
        S.push_back(s);
    return S;
}

//-------------------------
const std::vector<int> projections(unsigned j)
{
    std::vector<int> M;
    for (int m = -j; m <= static_cast<int>(j); ++m)
        M.push_back(m);
    return M;
}
    
//-------------------------
const std::vector<std::vector<int> > projections(const std::vector<unsigned>& J)
{
    // initialize vector of spin projections to -j
    std::vector<int> M;
    M.reserve(J.size());
    std::transform(J.begin(), J.end(), std::back_inserter(M), std::negate<int>());
    
    std::vector<std::vector<int> > SPV;
    // fill SPV with "odometer"-style looping
    while (M.back() <= (int)J.back()) {
        SPV.push_back(M);
        ++M[0];
        for (size_t i = 0; (i < M.size() - 1) and (M[i] > (int)J[i]); ++i) {
            M[i] = -J[i];
            ++M[i + 1];
        }
    }
    return SPV;
}

}
