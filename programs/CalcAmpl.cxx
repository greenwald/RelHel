// #include "TFhh.h"
// #include "TJSS.h"

#include "Contraction.h"
#include "JSS.h"
#include "QuantumNumbers.h"
#include "WaveFunction.h"

#include <iostream>
#include <stdexcept>
#include <vector>

using namespace relhel;

int main(int narg, char** args)
{

    if (narg < 4) {
        std::cout << "Usage: " << args[0] << " [Mother JP] [Daughter JP] [Daughter JP]" << std::endl
                  << "    JP takes the form [int][+/-]" << std::endl;
        return 0;
    }

    std::vector<QuantumNumbers> jp;

    for (int i = 1; i <= 3; ++i)
        jp.emplace_back(args[i]);

    // relhel::JSS jss(jp[0], jp[1], jp[2]);
    // jss.CalcAmpl();

    std::vector<WaveFunction> phi;
    phi.reserve(jp.size());
    std::transform(jp.begin(), jp.end(), std::back_inserter(phi), [](const QuantumNumbers& q){return WaveFunction(q.twoJ());});

    for (size_t i = 0; i < phi.size(); ++i)
        std::cout << "phi " << i << " :\n" << (spin(phi[i]) == 0 ? "|0, 0> = ()" : to_string(phi[i])) << std::endl;

    CoupledWaveFunctions psi(phi[1], phi[2], spin(phi[1]) + spin(phi[2]), 2);

    OrbitalAngularMomentumWaveFunction chi((spin(phi[0]) + psi.twoS()) / 2);
    std::cout << to_string(chi) << std::endl;
    
    Contraction C(psi, chi);

    std::cout << "Before contractions; " << C.terms().size() << " terms:" << std::endl
              << to_string(C.terms()) << std::endl;

    C.contract({Contraction::Term::psi1, Contraction::Term::psi2});

    std::cout << "After contracting psi1 and psi2 once; " << C.terms().size() << " terms:" << std::endl
              << to_string(C.terms()) << std::endl;

    C.contract({Contraction::Term::chi, Contraction::Term::psi2});

    std::cout << "After contracting chi with psi2 once; " << C.terms().size() << " terms:" << std::endl
              << to_string(C.terms()) << std::endl;

    return 0;
}

