// #include "TFhh.h"
// #include "TJSS.h"

#include "JSS.h"
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

    relhel::JSS jss(jp[0], jp[1], jp[2]);
    jss.CalcAmpl();

    return 0;
}

