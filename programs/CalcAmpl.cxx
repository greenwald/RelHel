// #include "TFhh.h"
// #include "TJSS.h"

#include "QuantumNumbers.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int narg, char** args)
{

    if (narg < 4) {
        std::cout << "Usage: " << args[0] << " [Mother JP] [Daughter JP] [Daughter JP]" << std::endl
                  << "    JP takes the form [int][+/-]" << std::endl;
        return 0;
    }

    std::vector<relhel::QuantumNumbers> jp;

    for (int i = 1; i <= 3; ++i)
        jp.emplace_back(args[i]);

    std::cout << "Mother particle:   " << relhel::to_string(jp[0]) << std::endl;
    std::cout << "1. decay particle: " << relhel::to_string(jp[1]) << std::endl;
    std::cout << "2. decay particle: " << relhel::to_string(jp[2]) << std::endl;

    // TJSS jss(j[0], p[0], j[1], p[1], j[2], p[2]);
    // jss.CalcAmpl();

    return 1;
}

