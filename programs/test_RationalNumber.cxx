#include "RationalNumber.h"

#include <iostream>
#include <string>

int main()
{

    for (unsigned d = 0; d <= 100; ++d) {
        auto D = relhel::decompose(d);
        std::cout << d << " = " << static_cast<unsigned>(D) << " = " << to_string(D) << std::endl;
    }
    
    // relhel::RationalNumber R(456, 237);
    // std::cout << 456./237. << "\t=\t" << to_string(R) << std::endl;
    
    return 0;
}
