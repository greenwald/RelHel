#include "ClebschGordan.h"
#include "RationalNumber.h"

#include <iostream>
#include <string>

int main()
{

    for (unsigned two_j1 = 0; two_j1 <= 2 * 2; ++two_j1)
        for (unsigned two_j2 = 0; two_j2 <= 2 * 2; ++two_j2)
            for (unsigned two_J = std::abs(two_j1 - two_j2); two_J <= two_j1 + two_j2; two_J += 2)
                for (int two_m1 = -two_j1; two_m1 <= (int)two_j2; two_m1 += 2)
                    for (int two_m2 = -two_j2; two_m2 <= (int)two_j2; two_m2 += 2)
                        std::cout << relhel::ClebschGordan::to_string(two_j1, two_m1, two_j2, two_m2, two_J) << " = " << 
                            relhel::to_sqrt_string(relhel::ClebschGordan::squared_coefficient(two_j1, two_m1, two_j2, two_m2, two_J)) << std::endl;
    
    return 0;
}
