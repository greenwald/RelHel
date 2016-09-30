#include "ClebschGordan.h"

#include "MathUtils.h"

#include <algorithm>
#include <cmath>

namespace relhel {

//-------------------------
std::string ClebschGordan::to_string(unsigned two_j1, int two_m1, unsigned two_j2, int two_m2, unsigned two_J, int two_M)
{
    return std::string("(") + spin_to_string(two_j1) + " " + spin_to_string(two_m1)
           + ", " + spin_to_string(two_j2) + " " + spin_to_string(two_m2)
           + " | " + spin_to_string(two_J) + " " + spin_to_string(two_M) + ")";
}

//-------------------------
const bool ClebschGordan::consistent(unsigned two_J, int two_M)
{
    return (std::abs(two_M) <= (int)two_J) and is_even(two_J + two_M);
}
    
//-------------------------
const bool ClebschGordan::nonzero(unsigned two_j1, int two_m1, unsigned two_j2, int two_m2, unsigned two_J, int two_M)
{
    // and that (j1+j2) and J are consistent
    if (is_odd(two_J + two_j1 + two_j2))
        return false;

    // check input spin-projection compatibilities
    if (!consistent(two_j1, two_m1) or !consistent(two_j2, two_m2) or !consistent(two_J,  two_M ))
        return false;

    // check input spin-projections
    if (two_M != two_m1 + two_m2)
        return false;

    // check whether J lies between |j1 - j2| and (j1 + j2)
    if ((int)two_J < std::abs((int)two_j1 - (int)two_j2) or two_J > (two_j1 + two_j2))
        return false;

    // when either daughter spin is zero
    if (two_j1 == 0 and two_j2 != two_J)
        return false;
    if (two_j2 == 0 and two_j1 != two_J)
        return false;

    // check for (j1 0 j2 0 | J 0), j1 + j2 + J is odd
    if (two_m1 == 0 and two_m2 == 0 and is_odd((two_j1 + two_j2 + two_J) / 2))
        return false;

    // (3/2 +-1/2, 3/2 +-1/2 | 2 +-1) == 0
    if (two_j1 == 3 and std::abs(two_m1) == 1 and two_j2 == 3 and two_m2 == two_m1 and two_J == 4)
        return false;

    // (2 +-1, 3/2 -+1/2 | 3/2 +-1/2) == 0
    if (two_j1 == 4 and std::abs(two_m1) == 2 and two_j2 == 3 and two_m2 == -two_m1 / 2 and two_J == 3)
        return false;

    // (2 +-1, 2 +-1 | 3 +-2) == 0
    if (two_j1 == 4 and std::abs(two_m1) == 2 and two_j2 == 4 and two_m2 == two_m1 and two_J == 6)
        return false;

    return true;
}
    
//-------------------------
const RationalNumber ClebschGordan::squared_coefficient(unsigned two_j1, int two_m1, unsigned two_j2, int two_m2, unsigned two_J, int two_M)
{
    if (!nonzero(two_j1, two_m1, two_j2, two_m2, two_J, two_M))
        return RationalNumber(0);

    // simple case of spin zero (must be 1 since check above already would have found 0)
    if (two_j1 == 0 or two_j2 == 0)
        return RationalNumber(1);

    // z range dictated by factorials in denominator ( 1/n! = 0 when n < 0)
    unsigned z_min = std::max<int>({0, (int)two_j2 - two_m1 - (int)two_J, (int)two_j1 + two_m2 - (int)two_J}) / 2;
    unsigned z_max = std::min<int>({(int)two_j1 + (int)two_j2 - (int)two_J, (int)two_j1 - two_m1, (int)two_j2 + two_m2}) / 2;

    RationalNumber z_sum(0);
    for (unsigned z = z_min; z <= z_max; ++z)
        // z'th term := (-)^z / z! / (j1+j2-J-z)! / (j1-m1-z)! / (j2+m2-z)! / (J-j2+m1+z)! / (J-j1-m2+z)!
        z_sum += RationalNumber(decompose(1),
                                decompose_factorial(z) *
                                decompose_factorial(((int)two_j1 + (int)two_j2 - (int)two_J) / 2 - z) *
                                decompose_factorial(((int)two_j1 - two_m1) / 2 - z) *
                                decompose_factorial(((int)two_j2 + two_m2) / 2 - z) *
                                decompose_factorial(((int)two_J - (int)two_j2 + two_m1) / 2 + z) *
                                decompose_factorial(((int)two_J - (int)two_j1 - two_m2) / 2 + z),
                                static_cast<double>(pow_negative_one(z)));
    
    // C-G coef = sqrt( (2J+1) (j1+j2-J)! (j1-j2+J)! (j2-j1+J)! / (J+j1+j2+1)! )
    //          * sqrt( (j1+m1)! (j1-m1)! (j2+m2)! (j2-m2)! (J+M)! (J-M)! )
    //          * z_sum
    return pow(z_sum, 2)
        * RationalNumber(decompose(two_J + 1)
                         * decompose_factorial(((int)two_j1 + (int)two_j2 - (int)two_J) / 2)
                         * decompose_factorial(((int)two_j1 - (int)two_j2 + (int)two_J) / 2)
                         * decompose_factorial(((int)two_j2 - (int)two_j1 + (int)two_J) / 2)
                         * decompose_factorial(((int)two_j1 + two_m1) / 2)
                         * decompose_factorial(((int)two_j1 - two_m1) / 2)
                         * decompose_factorial(((int)two_j2 + two_m2) / 2)
                         * decompose_factorial(((int)two_j2 - two_m2) / 2)
                         * decompose_factorial(((int)two_J + two_M) / 2)
                         * decompose_factorial(((int)two_J - two_M) / 2),
                         decompose_factorial(((int)two_j1 + (int)two_j2 + (int)two_J) / 2 + 1));
}
}
