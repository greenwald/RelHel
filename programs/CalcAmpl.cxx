// #include "TFhh.h"
// #include "TJSS.h"

#include "Parity.h"

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

    std::vector<int> j;
    std::vector<relhel::Parity> p;

    for (int i = 1; i <= 3; ++i) {
        std::string s = args[i];
        p.push_back(relhel::to_parity(s.back()));
        j.push_back(std::stoi(s.substr(0, s.length() - 1)));
    }

    std::cout << "Mother particle:   " << j[0] << relhel::to_string(p[0]) << std::endl;
    std::cout << "1. decay particle: " << j[1] << relhel::to_string(p[1]) << std::endl;
    std::cout << "2. decay particle: " << j[2] << relhel::to_string(p[2]) << std::endl;

    // TJSS jss(j[0], p[0], j[1], p[1], j[2], p[2]);
    // jss.CalcAmpl();

}

