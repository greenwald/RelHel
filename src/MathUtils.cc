#include "MathUtils.h"

#include <cstdlib>

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
