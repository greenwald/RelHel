// #include "TFhh.h"
// #include "TJSS.h"

#include "Contraction.h"
#include "WaveFunction.h"

#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>

using namespace relhel;

const bool one_opt(unsigned a, unsigned b)
{ return abs(static_cast<int>(a) - static_cast<int>(b)) == a + b; }

int main(int narg, char** args)
{

    if (narg < 4) {
        std::cout << "Usage: " << args[0] << " (Mother 2J) (Daughter 2J) (Daughter 2J) [L] [2S] [2M]" << std::endl;
        return 1;
    }

    // get spins
    std::vector<unsigned> two_j;
    two_j.reserve(3);
    for (unsigned i = 1; i <= 3; ++i) {
        int j2 = std::atoi(args[i]);
        if (j2 < 0)
            throw std::invalid_argument("spin " + std::to_string(i) + " is negative; main");
        two_j.push_back(abs(j2));
    }

    unsigned L = 0;
    if (narg < 5) {
        // can fix S?
        if (one_opt(two_j[1], two_j[2]) and one_opt(two_j[0], two_j[1] + two_j[2]))
            L = (two_j[0] + two_j[1] + two_j[2]) / 2;
        else
            throw std::invalid_argument("L is ambiguous, you must state it explicitly; main");
    } else
        L = std::atoi(args[4]);

    unsigned two_S = 0;
    if (narg < 6) {
        if (one_opt(two_j[0], 2 * L))
            two_S = two_j[0] + 2 * L;
        else if (one_opt(two_j[1], two_j[2]))
            two_S = two_j[1] + two_j[2];
        else
            throw std::invalid_argument("S is ambiguous, you must state it explicitly; main");
    } else
        two_S = std::atoi(args[5]);

    if (narg < 7 and two_j[0] != 0)
        throw std::invalid_argument("M is ambiguous, you must state it explicitly; main");
    
    int two_M = (narg < 7) ? 0 : std::atoi(args[6]);

    std::cout << "(" + spin_to_string(two_j[0]) + ") -> (" + spin_to_string(two_j[1]) + ") + (" + spin_to_string(two_j[2]) + ")"
              << ", L = " << L
              << ", S = " << spin_to_string(two_S)
              << ", M = " << spin_to_string(two_M)
              << std::endl;
    
    if (!triangle(two_j[1], two_j[2], two_S))
        throw std::invalid_argument("j1:j2:S triangle broken; main");

    if (!triangle(two_j[0], 2 * L, two_S))
        throw std::invalid_argument("j0:L:S triangle broken; main");

    WaveFunction phi(two_j[0]);
    WaveFunction psi1(two_j[1]);
    WaveFunction psi2(two_j[2]);

    OrbitalAngularMomentumWaveFunction chi(L);

    CoupledWaveFunctions psi(psi1, psi2, two_S, two_M);

    std::cout << to_string(contractions(psi, chi, phi)) << std::endl;
    
    // auto P = contract(psi, chi, phi);

    // std::cout << to_string(P) << std::endl;
    
    return 0;
}

